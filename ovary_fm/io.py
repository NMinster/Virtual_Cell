"""The manifest, the readers and the obs schema.

Every dataset enters as raw integer counts with gene symbols as
``var_names`` (Ensembl ids kept in ``var['gene_ids']`` when the file has
them) and leaves ``qc`` with the obs schema below, so that ``harmonise``
can pool datasets that never shared a column name.
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from .config import data_dir, resolve_path

# The obs schema every QC'd shard carries (docs/project_plan.md 4.2).
SCHEMA = ["dataset", "species", "assay", "suspension_type", "strain", "sex", "animal",
          "author_label", "condition", "condition_family", "condition_time_h", "age_days",
          "genotype", "treatment", "estrous_stage", "is_control"]
# Added by harmonise.
CROSSWALK_COLUMNS = ["major_type", "fine_label", "cl_id", "cl_label", "crosswalk_status"]


# --- the manifest ------------------------------------------------------------

def datasets(cfg, enabled_only: bool = True, names=None) -> list[dict]:
    entries = cfg["corpus"].get("datasets") or []
    if names:
        wanted = set(names)
        missing = wanted - {e["name"] for e in entries}
        if missing:
            raise KeyError(f"not in the manifest: {sorted(missing)}")
        entries = [e for e in entries if e["name"] in wanted]
    elif enabled_only:
        entries = [e for e in entries if e.get("enabled", False)]
    return entries


def dataset_path(cfg, entry: dict) -> Path:
    """A relative ``data/...`` path follows paths.data; anything else the repo root."""
    text = str(entry["path"])
    if Path(text).is_absolute():
        return Path(text)
    if text == "data" or text.startswith("data/"):
        return Path(cfg["paths"]["data"]) / text[5:]
    return resolve_path(text)


def qc_shard_path(cfg, name: str) -> Path:
    return data_dir(cfg, "corpus", "qc") / f"{name}.h5ad"


def corpus_path(cfg) -> Path:
    return data_dir(cfg, "corpus") / "corpus.h5ad"


def heldout_path(cfg, name: str) -> Path:
    return data_dir(cfg, "corpus") / f"heldout_{name}.h5ad"


# --- readers -----------------------------------------------------------------

def _is_integer(X, n: int = 20000) -> bool:
    data = X.data if sp.issparse(X) else np.asarray(X).ravel()
    if data.size == 0:
        return True
    sample = data[:n]
    return bool(np.allclose(sample, np.rint(sample)))


def read_counts(path: Path):
    """An .h5ad, a 10x .h5, a folder with one of those or an mtx / csv matrix."""
    import anndata as ad
    import scanpy as sc

    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"{path} does not exist (run `ovary-fm fetch` or fix the manifest path)")
    if path.is_dir():
        h5ad = sorted(path.glob("*.h5ad"))
        if h5ad:
            return ad.read_h5ad(h5ad[0])
        h5 = sorted(p for p in path.glob("*.h5") if "raw" not in p.name.lower()) or sorted(path.glob("*.h5"))
        if h5:
            return sc.read_10x_h5(h5[0])
        mtx = sorted(path.rglob("matrix.mtx*"))
        if mtx:
            return _read_mtx_dir(mtx[0].parent)
        csv = sorted(p for p in list(path.glob("*.csv*")) + list(path.glob("*.tsv*"))
                     if "meta" not in p.name.lower())
        if csv:
            return _read_csv_matrix(csv[0], path)
        raise FileNotFoundError(f"no .h5ad, .h5, matrix.mtx or csv matrix under {path} "
                                f"(a Seurat .rds is converted first with R/rds_to_mtx.R)")
    if path.suffix == ".h5ad":
        return ad.read_h5ad(path)
    if path.suffix == ".h5":
        return sc.read_10x_h5(path)
    if path.suffix in (".csv", ".tsv", ".gz"):
        return _read_csv_matrix(path, path.parent)
    raise ValueError(f"cannot read {path}")


def _read_mtx_dir(folder: Path):
    import scanpy as sc

    adata = sc.read_10x_mtx(folder, var_names="gene_symbols", make_unique=True)
    meta = sorted(folder.glob("*meta*.csv*")) + sorted(folder.parent.glob("*meta*.csv*"))
    if meta:
        obs = pd.read_csv(meta[0], index_col=0)
        adata.obs = adata.obs.join(obs, how="left")
    return adata


def _read_csv_matrix(path: Path, folder: Path):
    """A genes x cells count matrix with gene names in the first column; an
    optional *meta*.csv beside it with one row per cell."""
    import anndata as ad

    sep = "\t" if ".tsv" in path.name else ","
    frame = pd.read_csv(path, sep=sep, index_col=0)
    X = sp.csr_matrix(frame.to_numpy(dtype=np.float32).T)
    adata = ad.AnnData(X=X, obs=pd.DataFrame(index=frame.columns.astype(str)),
                       var=pd.DataFrame(index=frame.index.astype(str)))
    meta = sorted(folder.glob("*meta*.csv*"))
    if meta:
        obs = pd.read_csv(meta[0], index_col=0)
        adata.obs = adata.obs.join(obs, how="left")
    return adata


def standardise(adata, gene_key: str | None = None, dataset: str | None = None):
    """Integer counts in X (csr int32), symbols in var_names, Ensembl ids in
    var['gene_ids'] when present, unique obs_names prefixed with the dataset."""
    var = adata.var
    X = adata.X
    if not _is_integer(X):
        if "counts" in adata.layers and _is_integer(adata.layers["counts"]):
            X = adata.layers["counts"]
        elif adata.raw is not None and _is_integer(adata.raw.X):
            X = adata.raw.X
            var = adata.raw.var.copy()
            adata = adata[:, :].copy()
            import anndata as ad
            adata = ad.AnnData(X=X, obs=adata.obs.copy(), var=var)
        else:
            raise ValueError(f"{dataset}: X holds normalised values and neither layers['counts'] nor "
                             f".raw carries integer counts; the corpus takes raw counts only")
    X = sp.csr_matrix(X)
    X.data = np.rint(X.data).astype(np.float32)
    X.eliminate_zeros()
    adata.X = X
    for layer in [k for k in list(adata.layers.keys()) if k is not None]:   # anndata >= 0.13 lists None (X) as a layer
        del adata.layers[layer]
    if adata.raw is not None:
        del adata.raw
    var = adata.var
    names = pd.Index(var.index.astype(str))
    looks_ensembl = names.str.match(r"^ENS[A-Z]*G\d+").mean() > 0.9
    if gene_key and gene_key in var:
        symbols = var[gene_key].astype(str)
        if looks_ensembl:
            var["gene_ids"] = names
        adata.var_names = symbols.to_numpy()
    elif looks_ensembl:
        for key in ("feature_name", "gene_symbols", "symbol", "gene_name"):
            if key in var:
                var["gene_ids"] = names
                adata.var_names = var[key].astype(str).to_numpy()
                break
        else:
            print(f"  !! {dataset}: var_names are Ensembl ids and no symbol column was found; "
                  f"set gene_key in the manifest")
    elif "gene_ids" not in var:
        for key in ("feature_id", "ensembl_id", "gene_id"):
            if key in var:
                var["gene_ids"] = var[key].astype(str)
                break
    adata.var_names = pd.Index(adata.var_names.astype(str))
    adata.var_names_make_unique()
    adata.obs_names = pd.Index([f"{dataset}:{b}" for b in adata.obs_names.astype(str)]) if dataset \
        else adata.obs_names.astype(str)
    adata.obs_names_make_unique()
    adata.var = adata.var.loc[:, [c for c in adata.var.columns if c in ("gene_ids",)]]
    return adata


# --- inspect -----------------------------------------------------------------

def inspect(cfg):
    import h5py

    rows = []
    for entry in datasets(cfg, enabled_only=False):
        path = dataset_path(cfg, entry)
        shape = None
        if path.exists() and path.suffix == ".h5ad":
            try:
                with h5py.File(path, "r") as h:
                    X = h["X"]
                    shape = tuple(X.attrs["shape"]) if "shape" in X.attrs else X.shape
            except Exception:
                shape = "?"
        rows.append({"dataset": entry["name"], "enabled": entry.get("enabled", False),
                     "tier": entry.get("tier"), "species": entry.get("species"),
                     "source": (entry.get("source") or {}).get("type", "local"),
                     "accession": (entry.get("source") or {}).get("accession"),
                     "status": entry.get("status", "ready"), "held_out": entry.get("held_out", False),
                     "path": str(path), "on_disk": path.exists(), "shape": shape,
                     "qc_shard": qc_shard_path(cfg, entry["name"]).exists()})
    frame = pd.DataFrame(rows)
    with pd.option_context("display.max_columns", 20, "display.width", 200, "display.max_colwidth", 60):
        print(frame.to_string(index=False))
    out = data_dir(cfg) / "inspect.csv"
    frame.to_csv(out, index=False)
    print(f"\n  -> {out}")
    return frame


def read_json(path: Path) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def open_text(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt", encoding="utf-8")
