"""``ovary-fm qc``: one dataset at a time, before anything is pooled.

The Mouse_Compensatory rules, applied per dataset and per animal: the hard
filters, a mito cap that depends on the modality (a nucleus above 5 % is
cytoplasm), per-animal MAD bounds on depth and genes, an optional doublet
call, and a minimum number of cells per animal. The modality comes from
the manifest or from the tells (median mito, the Malat1 share), and the
decision is printed and written, never silent.

Every shard leaves with the obs schema (``io.SCHEMA``) and a drop-off
table, so the corpus can be audited dataset by dataset.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp

from .config import tables_dir, write_provenance
from .io import SCHEMA, dataset_path, datasets, qc_shard_path, read_counts, standardise


def cell_metrics(adata, species: str, cfg) -> pd.DataFrame:
    q = cfg["qc"]
    X = sp.csr_matrix(adata.X)
    names = pd.Index(adata.var_names.astype(str))
    total = np.asarray(X.sum(axis=1)).ravel()
    genes = np.asarray((X > 0).sum(axis=1)).ravel()
    mito = names.str.startswith(q["mito_prefix"].get(species, "mt-"))
    nuclear = names.isin(q["nuclear_genes"].get(species, []))
    ribo = np.zeros(len(names), dtype=bool)
    for prefix in q["ribo_prefix"].get(species, []):
        ribo |= names.str.startswith(prefix)
    with np.errstate(divide="ignore", invalid="ignore"):
        frame = pd.DataFrame({
            "total_counts": total,
            "n_genes": genes,
            "pct_mito": 100 * np.asarray(X[:, mito].sum(axis=1)).ravel() / np.maximum(total, 1),
            "share_nuclear": np.asarray(X[:, nuclear].sum(axis=1)).ravel() / np.maximum(total, 1),
            "share_ribo": np.asarray(X[:, ribo].sum(axis=1)).ravel() / np.maximum(total, 1),
        }, index=adata.obs_names)
    return frame


def infer_modality(metrics: pd.DataFrame, cfg) -> tuple[str, dict]:
    m = cfg["qc"]["modality"]
    tells = {"median_pct_mito": float(np.median(metrics["pct_mito"])),
             "median_share_nuclear": float(np.median(metrics["share_nuclear"])),
             "median_share_ribo": float(np.median(metrics["share_ribo"]))}
    nucleus = (tells["median_pct_mito"] < m["nucleus_if_mito_pct_below"]
               and tells["median_share_nuclear"] > m["nucleus_if_malat1_share_above"])
    return ("nucleus" if nucleus else "cell"), tells


def _mad_bounds(values: np.ndarray, n_mads: float) -> tuple[float, float]:
    med = np.median(values)
    mad = np.median(np.abs(values - med)) * 1.4826
    return med - n_mads * mad, med + n_mads * mad


def _condition_columns(obs: pd.DataFrame, entry: dict) -> pd.DataFrame:
    """condition / family / time / age from the manifest's `condition` block."""
    spec = entry.get("condition") or {}
    n = len(obs)
    if "value" in spec:
        cond = pd.Series([str(spec["value"])] * n, index=obs.index)
    elif spec.get("key") and spec["key"] in obs:
        raw = obs[spec["key"]].astype(str).str.strip()
        mapping = {str(k).strip().lower(): v for k, v in (spec.get("map") or {}).items()}
        if mapping:
            cond = raw.str.lower().map(mapping)
            unmapped = raw[cond.isna()].unique()
            if len(unmapped):
                print(f"  !! condition values not in the map, kept as given: {list(unmapped)[:10]}")
            cond = cond.fillna(raw)
        else:
            cond = raw
    else:
        if spec.get("key"):
            print(f"  !! condition key {spec['key']!r} is not an obs column; condition = 'unknown'")
        cond = pd.Series(["unknown"] * n, index=obs.index)
    time_h = {str(k): v for k, v in (spec.get("time_h") or {}).items()}
    age = {str(k): v for k, v in (spec.get("age_days") or {}).items()}
    control = spec.get("control")
    return pd.DataFrame({
        "condition": cond.astype(str).to_numpy(),
        "condition_family": spec.get("family", "none"),
        "condition_time_h": cond.map(time_h).astype(float).to_numpy(),
        "age_days": cond.map(age).astype(float).to_numpy() if age else np.full(n, entry.get("age_days", np.nan), dtype=float),
        "is_control": (cond.astype(str) == str(control)).to_numpy() if control is not None else np.zeros(n, dtype=bool),
    }, index=obs.index)


def apply_schema(adata, entry: dict, modality: str) -> None:
    obs = adata.obs
    name = entry["name"]
    new = pd.DataFrame(index=obs.index)
    new["dataset"] = name
    new["species"] = entry.get("species", "mouse")
    new["assay"] = entry.get("assay") or "unknown"
    new["suspension_type"] = modality
    new["strain"] = entry.get("strain") or "unknown"
    new["sex"] = entry.get("sex") or "female"
    key = entry.get("animal_key")
    new["animal"] = (obs[key].astype(str).to_numpy() if key and key in obs else name)
    if key and key not in obs:
        print(f"  !! animal_key {key!r} is not an obs column; the whole dataset counts as one animal")
    label = entry.get("label_key")
    new["author_label"] = (obs[label].astype(str).to_numpy() if label and label in obs else "unlabelled")
    if label and label not in obs:
        print(f"  !! label_key {label!r} is not an obs column; author_label = 'unlabelled'")
    cond = _condition_columns(obs, entry)
    for col in cond.columns:
        new[col] = cond[col].to_numpy()
    for col in ("genotype", "treatment", "estrous_stage"):
        new[col] = obs[col].astype(str).to_numpy() if col in obs else "NA"
    # the raw columns the schema was built from, for the audit
    for col in {key, label, (entry.get("condition") or {}).get("key")} - {None}:
        if col in obs:
            new[f"orig_{col}"] = obs[col].astype(str).to_numpy()
    adata.obs = new[SCHEMA + [c for c in new.columns if c.startswith("orig_")]]


def qc_one(cfg, entry: dict):
    q = cfg["qc"]
    name = entry["name"]
    species = entry.get("species", "mouse")
    adata = standardise(read_counts(dataset_path(cfg, entry)), gene_key=entry.get("gene_key"), dataset=name)
    print(f"  read {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    metrics = cell_metrics(adata, species, cfg)
    modality = entry.get("suspension_type")
    inferred, tells = infer_modality(metrics, cfg)
    if modality not in ("cell", "nucleus"):
        modality = inferred
        print(f"  modality inferred: {modality} (median mito {tells['median_pct_mito']:.2f} %, "
              f"Malat1 share {tells['median_share_nuclear']:.3f})")
    else:
        flag = "" if modality == inferred else "  <- the tells say otherwise; check the manifest"
        print(f"  modality from the manifest: {modality} (tells: {inferred}; median mito "
              f"{tells['median_pct_mito']:.2f} %, Malat1 share {tells['median_share_nuclear']:.3f}){flag}")
    for col in metrics.columns:
        adata.obs[col] = metrics[col].to_numpy()
    apply_schema(adata, entry, modality)
    for col in metrics.columns:
        adata.obs[col] = metrics[col].to_numpy()

    keep = np.ones(adata.n_obs, dtype=bool)
    steps = []

    def step(label, mask):
        nonlocal keep
        before = int(keep.sum())
        keep &= mask
        steps.append({"step": label, "cells_before": before, "cells_after": int(keep.sum()),
                      "removed": before - int(keep.sum())})

    obs = adata.obs
    step(f"min_counts >= {q['min_counts']}", obs["total_counts"].to_numpy() >= q["min_counts"])
    step(f"min_genes >= {q['min_genes']}", obs["n_genes"].to_numpy() >= q["min_genes"])
    step(f"max_genes <= {q['max_genes']}", obs["n_genes"].to_numpy() <= q["max_genes"])
    cap = q["max_mito_pct"][modality]
    step(f"mito <= {cap} % ({modality})", obs["pct_mito"].to_numpy() <= cap)
    if q["mad"].get("enabled", True):
        n_mads = float(q["mad"].get("n_mads", 3))
        mask = np.ones(adata.n_obs, dtype=bool)
        for animal in obs["animal"].unique():
            idx = np.where((obs["animal"] == animal).to_numpy() & keep)[0]
            if len(idx) < 20:
                continue
            for col in ("total_counts", "n_genes"):
                values = np.log1p(obs[col].to_numpy()[idx])
                lo, hi = _mad_bounds(values, n_mads)
                mask[idx] &= (values >= lo) & (values <= hi)
        step(f"per-animal MAD ({n_mads} MADs, log counts and genes)", mask)
    if q["doublets"].get("enabled", False):
        mask = _doublets(adata, keep, q["doublets"], cfg.get("random_state", 0))
        step("doublets (scrublet)", mask)
    counts = obs["animal"][keep].value_counts()
    small = counts[counts < q["min_cells_per_animal"]].index.tolist()
    if small:
        print(f"  animals below {q['min_cells_per_animal']} cells, dropped: {small}")
    step(f"animals with >= {q['min_cells_per_animal']} cells", ~obs["animal"].isin(small).to_numpy())

    adata = adata[keep].copy()
    import json
    adata.uns["qc"] = {"modality": modality, "tells": tells, "steps_json": json.dumps(steps)}
    return adata, pd.DataFrame(steps), tells


def _doublets(adata, keep, spec, seed) -> np.ndarray:
    import scanpy as sc

    mask = np.ones(adata.n_obs, dtype=bool)
    for animal in adata.obs["animal"].unique():
        idx = np.where((adata.obs["animal"] == animal).to_numpy() & keep)[0]
        if len(idx) < 100:
            continue
        sub = adata[idx].copy()
        expected = min(0.5, spec.get("rate_per_1000", 0.008) * len(idx) / 1000)
        try:
            sc.pp.scrublet(sub, expected_doublet_rate=max(expected, 0.01), random_state=seed, verbose=False)
            mask[idx] = ~sub.obs["predicted_doublet"].astype(bool).to_numpy()
        except Exception as exc:
            print(f"  !! scrublet failed for {animal}: {exc}")
    return mask


def run(cfg, dataset: str | None = None):
    write_provenance(cfg, "qc", extra={"dataset": dataset})
    tables = tables_dir(cfg, "qc")
    summary = []
    for entry in datasets(cfg, names=[dataset] if dataset else None):
        print(f"== {entry['name']}")
        adata, dropoff, tells = qc_one(cfg, entry)
        dropoff.to_csv(tables / f"{entry['name']}_dropoff.csv", index=False)
        out = qc_shard_path(cfg, entry["name"])
        adata.write_h5ad(out)
        obs = adata.obs
        summary.append({"dataset": entry["name"], "species": entry.get("species"),
                        "suspension_type": obs["suspension_type"].iloc[0] if len(obs) else None,
                        "cells_in": int(dropoff["cells_before"].iloc[0]) if len(dropoff) else 0,
                        "cells_out": adata.n_obs, "n_genes": adata.n_vars,
                        "n_animals": int(obs["animal"].nunique()),
                        "conditions": ";".join(sorted(obs["condition"].unique())),
                        "n_author_labels": int(obs["author_label"].nunique()),
                        "median_counts": float(obs["total_counts"].median()) if len(obs) else np.nan,
                        "median_genes": float(obs["n_genes"].median()) if len(obs) else np.nan,
                        **{f"tell_{k}": round(v, 4) for k, v in tells.items()},
                        "held_out": bool(entry.get("held_out", False)), "shard": str(out)})
        print(f"  {adata.n_obs:,} cells kept -> {out}")
    frame = pd.DataFrame(summary)
    path = tables / "qc_summary.csv"
    if path.exists() and dataset:
        old = pd.read_csv(path)
        frame = pd.concat([old[~old["dataset"].isin(frame["dataset"])], frame], ignore_index=True)
    frame.to_csv(path, index=False)
    print(f"  -> {path}")
    return frame
