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

import re

from .config import data_dir, resolve_path

# GEO "RAW" deposits: one matrix per sample, flat files named <GSM>_<title>_matrix.mtx.gz
# (+ barcodes / features or genes), or one 10x .h5 per sample
GEO_SUFFIX = {
    "matrix": re.compile(r"[_.-](matrix)\.mtx(\.gz)?$", re.I),
    "barcodes": re.compile(r"[_.-](barcodes|cells)\.tsv(\.gz)?$", re.I),
    "features": re.compile(r"[_.-](features|genes)\.tsv(\.gz)?$", re.I),
    "h5": re.compile(r"(?:[_.-](?:filtered|raw)_feature_bc_matrix)?\.h5$", re.I),
    # one count matrix per sample as a text table (GEO "TAR (of CSV)"): genes x cells or cells x genes
    "csv": re.compile(r"(?:[_.-](?:counts?|matrix|expression|umi|raw))?\.(csv|tsv|txt)(\.gz)?$", re.I),
}
_META_WORDS = ("meta", "cell_info", "annotation", "barcodes", "features", "genes", "readme", "metadata", "donor", "sample")
# a few genes every mouse or human matrix carries, to tell which axis holds the genes
_ANCHOR_GENES = {"actb", "malat1", "gapdh", "ptprc", "b2m", "mt-co1", "mt-co1", "malat1"}

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
        if len(h5ad) == 1:
            return ad.read_h5ad(h5ad[0])
        plain_mtx = sorted(p for p in path.rglob("matrix.mtx*") if p.name.lower().startswith("matrix.mtx"))
        if len(plain_mtx) == 1:
            return _read_mtx_dir(plain_mtx[0].parent)
        groups = geo_raw_groups(path)
        if groups:
            return _read_geo_raw(groups)
        h5 = sorted(p for p in path.glob("*.h5") if "raw" not in p.name.lower()) or sorted(path.glob("*.h5"))
        if len(h5) == 1:
            return sc.read_10x_h5(h5[0])
        csv = sorted(p for p in list(path.glob("*.csv*")) + list(path.glob("*.tsv*"))
                     if "meta" not in p.name.lower() and "cell_info" not in p.name.lower())
        if csv:
            return _read_csv_matrix(csv[0], path)
        raise FileNotFoundError(f"no .h5ad, .h5, matrix.mtx, GEO per-sample matrices or csv matrix under {path} "
                                f"(a Seurat .rds is converted first with R/rds_to_mtx.R)")
    if path.suffix == ".h5ad":
        return ad.read_h5ad(path)
    if path.suffix == ".h5":
        return sc.read_10x_h5(path)
    if path.suffix in (".csv", ".tsv", ".gz"):
        return _read_csv_matrix(path, path.parent)
    raise ValueError(f"cannot read {path}")


def geo_raw_groups(folder: Path) -> dict:
    """prefix -> {kind: path} for every per-sample matrix set under a folder:
    flat <GSM>_<title>_{matrix.mtx,barcodes.tsv,features.tsv}(.gz) triplets, per-sample
    10x .h5 files, or one matrix.mtx folder per sample."""
    groups = {}
    for file in sorted(folder.rglob("*")):
        if not file.is_file() or file.suffix == ".h5ad":
            continue
        name = file.name
        if name.lower().startswith("matrix.mtx"):           # a 10x folder per sample
            groups.setdefault(file.parent.name, {})["dir"] = file.parent
            continue
        for kind, rx in GEO_SUFFIX.items():
            m = rx.search(name)
            if m:
                prefix = name[: m.start()]
                if kind == "h5" and (not prefix or name.endswith(".h5ad")):
                    continue
                if kind == "csv" and (not prefix or any(w in name.lower() for w in _META_WORDS)):
                    continue
                groups.setdefault(prefix, {})[kind] = file
                break
    complete = {}
    for prefix, g in groups.items():
        if "dir" in g or "h5" in g or {"matrix", "barcodes", "features"} <= set(g):
            complete[prefix] = {k: v for k, v in g.items() if k != "csv"}
        elif "csv" in g and len(groups) > 1:
            complete[prefix] = g                     # per-sample text matrices (a lone csv is the single-matrix path)
    # a deposit that carries both the raw (every droplet) and the filtered (called cells) matrix
    # of a sample: take the filtered one
    if any("filtered" in p.lower() for p in complete) and any("raw" in p.lower() for p in complete):
        complete = {p: g for p, g in complete.items() if "raw" not in p.lower()}
    return complete


def orient_genes_by_cells(frame: pd.DataFrame, what: str = "") -> pd.DataFrame:
    """Return the table as genes x cells. The axis that carries anchor genes is the gene
    axis; with no anchor on either axis the taller axis is taken to be genes, and it is said."""
    idx = {str(v).lower() for v in frame.index[:50000]}
    cols = {str(v).lower() for v in frame.columns[:50000]}
    if idx & _ANCHOR_GENES:
        return frame
    if cols & _ANCHOR_GENES:
        return frame.T
    print(f"  !! {what}: no anchor gene on either axis; taking the taller axis ({'rows' if frame.shape[0] >= frame.shape[1] else 'columns'}) as genes")
    return frame if frame.shape[0] >= frame.shape[1] else frame.T


def _read_geo_sample(prefix: str, g: dict):
    import anndata as ad
    import scanpy as sc
    from scipy import io as sio

    if "dir" in g:
        a = sc.read_10x_mtx(g["dir"], var_names="gene_symbols", make_unique=True)
    elif "h5" in g:
        a = sc.read_10x_h5(g["h5"])
    elif "csv" in g:
        path = g["csv"]
        sep = "\t" if (".tsv" in path.name or ".txt" in path.name) else ","
        frame = orient_genes_by_cells(pd.read_csv(path, sep=sep, index_col=0), prefix)
        a = ad.AnnData(X=sp.csr_matrix(frame.to_numpy(dtype=np.float32).T),
                       obs=pd.DataFrame(index=pd.Index(frame.columns.astype(str).to_numpy())),
                       var=pd.DataFrame(index=pd.Index(frame.index.astype(str).to_numpy())))
    else:
        M = sio.mmread(str(g["matrix"]))
        X = sp.csr_matrix(M.T.astype(np.float32))                      # mtx is genes x cells
        feats = pd.read_csv(g["features"], sep="\t", header=None, dtype=str)
        symbols = feats[1] if feats.shape[1] >= 2 else feats[0]
        var = pd.DataFrame(index=pd.Index(symbols.astype(str).to_numpy()))      # a plain index: no column name
        if feats.shape[1] >= 2:
            var["gene_ids"] = feats[0].astype(str).to_numpy()
        barcodes = pd.read_csv(g["barcodes"], sep="\t", header=None, dtype=str)[0].astype(str)
        a = ad.AnnData(X=X, obs=pd.DataFrame(index=pd.Index(barcodes.to_numpy())), var=var)
    a.var_names_make_unique()
    a.obs["gsm_sample"] = prefix
    a.obs_names = pd.Index([f"{prefix}_{b}" for b in a.obs_names.astype(str)])
    return a


def _read_geo_raw(groups: dict):
    import anndata as ad

    adatas = []
    for prefix, g in sorted(groups.items()):
        a = _read_geo_sample(prefix, g)
        print(f"  sample {prefix}: {a.n_obs:,} cells x {a.n_vars:,} genes")
        adatas.append(a)
    merged = ad.concat(adatas, join="outer", merge="first", index_unique=None)
    merged.obs["gsm_sample"] = merged.obs["gsm_sample"].astype(str)
    return merged


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

    sep = "\t" if (".tsv" in path.name or ".txt" in path.name) else ","
    frame = orient_genes_by_cells(pd.read_csv(path, sep=sep, index_col=0), path.name)
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
        samples = None
        if path.is_dir():
            groups = geo_raw_groups(path)
            if groups:
                samples = sorted(groups)
        rows.append({"dataset": entry["name"], "enabled": entry.get("enabled", False),
                     "tier": entry.get("tier"), "species": entry.get("species"),
                     "source": (entry.get("source") or {}).get("type", "local"),
                     "accession": (entry.get("source") or {}).get("accession"),
                     "status": entry.get("status", "ready"), "held_out": entry.get("held_out", False),
                     "path": str(path), "on_disk": path.exists(), "shape": shape,
                     "n_samples": len(samples) if samples else None,
                     "qc_shard": qc_shard_path(cfg, entry["name"]).exists()})
        if samples:
            print(f"  {entry['name']}: {len(samples)} per-sample matrices; the gsm_sample values the manifest's "
                  f"patterns must match:")
            for name in samples[:40]:
                print(f"    {name}")
            if len(samples) > 40:
                print(f"    ... {len(samples) - 40} more")
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
