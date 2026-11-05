"""``ovary-fm harmonise``: the crosswalk, the gene space, the corpus.

* Every QC'd shard gets ``major_type`` (the eight types), ``fine_label``,
  ``cl_id``, ``cl_label`` and ``crosswalk_status`` from config/crosswalk.yaml.
  A label the crosswalk does not list is kept as text with major type
  ``Unassigned`` and status ``unmapped``; it is listed in
  tables/harmonise/crosswalk_review.csv, the table that goes to the lab.
* The gene space is the union (or intersection) of symbols across the
  pooled datasets; a gene absent from a dataset is zero there and
  ``var['n_datasets_measured']`` says how many datasets measured it.
* Datasets with ``held_out: true`` are written as their own shards in the
  same gene space and never into corpus.h5ad.
"""

from __future__ import annotations

import re

import numpy as np
import pandas as pd
import scipy.sparse as sp

from .config import tables_dir, write_provenance
from .io import CROSSWALK_COLUMNS, corpus_path, datasets, heldout_path, qc_shard_path


def _norm(label) -> str:
    return str(label).strip().lower()


_STEM_SUFFIX = re.compile(r"(?:[\s_.-]*(?:cells?|cluster)|[\s_.-]*\d+[a-z]?|[\s_.-]+)+$", re.I)


def label_stem(label: str) -> str:
    """'Granulosa 3' -> 'granulosa', 'Theca_00' -> 'theca', 'Luteal Cells' -> 'luteal', 'GC_1' -> 'gc'."""
    return _STEM_SUFFIX.sub("", _norm(label)).strip()


def crosswalk_lookup(cw: dict, dataset: str, label: str) -> dict:
    """One author label -> {major, fine, cl_id, cl_label, status}: the dataset block, the
    global map, then both again with the label's trailing number / suffix removed."""
    key = _norm(label)
    per_dataset = {_norm(k): v for k, v in (cw.get("datasets", {}).get(dataset) or {}).items()}
    global_map = {_norm(k): v for k, v in (cw.get("global") or {}).items()}
    stem = label_stem(label)
    if key in per_dataset:
        entry, status = per_dataset[key], "mapped_dataset"
    elif key in global_map:
        entry, status = global_map[key], "mapped_global"
    elif stem and stem != key and stem in per_dataset:
        entry, status = per_dataset[stem], "mapped_stem"
    elif stem and stem != key and stem in global_map:
        entry, status = global_map[stem], "mapped_stem"
    else:
        return {"major": cw.get("unmapped_label", "Unassigned"), "fine": None, "cl_id": None,
                "cl_label": None, "status": "unmapped"}
    major = entry.get("major")
    fine = entry.get("fine")
    terms = cw.get("cl_terms") or {}
    if fine and fine in (terms.get("fine") or {}):
        term = terms["fine"][fine]
        cl_id, cl_label = term.get("cl_id"), fine
        if term.get("major") and major is None:
            major = term["major"]
    else:
        term = (terms.get("major") or {}).get(major, {})
        cl_id, cl_label = term.get("cl_id"), term.get("cl_label")
    if cl_id == "verify":
        cl_id, status = None, status + "_verify_cl_id"
    return {"major": major, "fine": fine, "cl_id": cl_id, "cl_label": cl_label, "status": status}


def apply_crosswalk(adata, cw: dict, dataset: str, unmapped_label: str) -> pd.DataFrame:
    majors = set(cw.get("major_types") or []) | {unmapped_label}
    labels = adata.obs["author_label"].astype(str)
    rows = []
    out = {c: [] for c in CROSSWALK_COLUMNS}
    cache = {}
    for label in labels.unique():
        hit = crosswalk_lookup(cw, dataset, label)
        if hit["major"] not in majors:
            raise ValueError(f"{dataset}: crosswalk maps {label!r} to {hit['major']!r}, not a major type")
        cache[label] = hit
        rows.append({"dataset": dataset, "author_label": label, "n_cells": int((labels == label).sum()),
                     "major_type": hit["major"], "fine_label": hit["fine"], "cl_id": hit["cl_id"],
                     "cl_label": hit["cl_label"], "crosswalk_status": hit["status"]})
    for label in labels:
        hit = cache[label]
        out["major_type"].append(hit["major"])
        out["fine_label"].append(hit["fine"] or "")
        out["cl_id"].append(hit["cl_id"] or "")
        out["cl_label"].append(hit["cl_label"] or "")
        out["crosswalk_status"].append(hit["status"])
    for col in CROSSWALK_COLUMNS:
        adata.obs[col] = pd.Categorical(out[col]) if col != "fine_label" else out[col]
    adata.obs["major_type"] = adata.obs["major_type"].astype(str)
    return pd.DataFrame(rows)


def reindex_genes(adata, genes: pd.Index):
    """Place the dataset's genes into the shared gene space (zeros elsewhere)."""
    import anndata as ad

    old = pd.Index(adata.var_names.astype(str))
    common = old.intersection(genes)
    old_idx = old.get_indexer(common)
    new_idx = genes.get_indexer(common)
    P = sp.csr_matrix((np.ones(len(common), dtype=np.float32), (old_idx, new_idx)),
                      shape=(len(old), len(genes)))
    X = sp.csr_matrix(adata.X).dot(P).tocsr()
    var = pd.DataFrame(index=genes)
    if "gene_ids" in adata.var:
        ids = pd.Series(adata.var["gene_ids"].astype(str).to_numpy(), index=old)
        var["gene_ids"] = ids.reindex(genes).fillna("").to_numpy()
    measured = np.zeros(len(genes), dtype=bool)
    measured[new_idx] = True
    var["measured_here"] = measured
    new = ad.AnnData(X=X, obs=adata.obs.copy(), var=var, uns=dict(adata.uns))
    return new


def transfer_labels(cfg, shards: dict, manifest: list) -> list:
    """Returns the names of the datasets whose labels were transferred."""
    import anndata as ad

    from .bench.t1 import _fit_predict, _subsample  # noqa: F401

    unlabelled = [n for n, (a, _) in shards.items() if (a.obs["author_label"].astype(str) == "unlabelled").all()]
    if not unlabelled:
        return []
    labelled = [a for n, (a, held) in shards.items() if n not in unlabelled and not held
                and (a.obs["crosswalk_status"].astype(str) != "unmapped").any()]
    if not labelled:
        print(f"  !! {unlabelled}: no labelled dataset to transfer from; they stay Unassigned")
        return []
    genes = None
    for a in labelled:
        genes = pd.Index(a.var_names.astype(str)) if genes is None else genes.intersection(a.var_names.astype(str))
    train = ad.concat([a[:, genes].copy() for a in labelled], join="inner", index_unique=None)
    train = train[train.obs["crosswalk_status"].astype(str) != "unmapped"].copy()
    t1cfg = dict(cfg["benchmark"]["t1"])
    min_conf = float(cfg["harmonise"].get("transfer_min_confidence", 0.5))
    for name in unlabelled:
        a, held = shards[name]
        common = genes.intersection(a.var_names.astype(str))
        test = a[:, common].copy()
        tr = train[:, common].copy()
        pred, conf = _fit_predict(tr, test, t1cfg, int(cfg.get("random_state", 0)), None, return_confidence=True)
        pred = np.where(conf >= min_conf, pred, cfg["harmonise"].get("unmapped_label", "Unassigned"))
        a.obs["major_type"] = pd.Series(pred, index=a.obs_names).astype(str).to_numpy()
        a.obs["transfer_confidence"] = conf
        a.obs["crosswalk_status"] = np.where(conf >= min_conf, "transferred", "transferred_low_confidence")
        terms = cfg["crosswalk"].get("cl_terms", {}).get("major", {})
        a.obs["cl_id"] = [str((terms.get(m) or {}).get("cl_id") or "") for m in a.obs["major_type"]]
        a.obs["cl_label"] = [str((terms.get(m) or {}).get("cl_label") or "") for m in a.obs["major_type"]]
        share = float((conf >= min_conf).mean())
        print(f"  {name}: labels transferred from {[x.obs['dataset'].iloc[0] for x in labelled]}; "
              f"{share:.1%} of cells at confidence >= {min_conf}")
        for row in manifest:
            if row["dataset"] == name:
                row["crosswalk_coverage"] = round(share, 4)
                row["labels"] = "transferred"
    return unlabelled


def run(cfg):
    write_provenance(cfg, "harmonise")
    h = cfg["harmonise"]
    cw = dict(cfg["crosswalk"])
    cw["unmapped_label"] = h.get("unmapped_label", "Unassigned")
    import anndata as ad

    tables = tables_dir(cfg, "harmonise")
    shards, review, manifest = {}, [], []
    for entry in datasets(cfg):
        path = qc_shard_path(cfg, entry["name"])
        if not path.exists():
            print(f"  !! {entry['name']}: no QC shard at {path}; run `ovary-fm qc` first (skipped)")
            continue
        adata = ad.read_h5ad(path)
        table = apply_crosswalk(adata, cw, entry["name"], cw["unmapped_label"])
        review.append(table)
        mapped = (adata.obs["crosswalk_status"] != "unmapped").mean() if adata.n_obs else 0.0
        if h.get("drop_unmapped", False):
            adata = adata[adata.obs["crosswalk_status"] != "unmapped"].copy()
        shards[entry["name"]] = (adata, bool(entry.get("held_out", False)))
        manifest.append({"dataset": entry["name"], "species": entry.get("species"), "tier": entry.get("tier"),
                         "held_out": bool(entry.get("held_out", False)), "n_cells": adata.n_obs,
                         "n_genes": adata.n_vars, "n_animals": int(adata.obs["animal"].nunique()),
                         "suspension_type": str(adata.obs["suspension_type"].iloc[0]),
                         "conditions": ";".join(sorted(adata.obs["condition"].astype(str).unique())),
                         "condition_family": str(adata.obs["condition_family"].iloc[0]),
                         "n_author_labels": int(adata.obs["author_label"].nunique()),
                         "crosswalk_coverage": round(float(mapped), 4),
                         "n_unmapped_labels": int((table["crosswalk_status"] == "unmapped").sum())})
        print(f"  {entry['name']}: {adata.n_obs:,} cells, crosswalk coverage {mapped:.1%}, "
              f"{int((table['crosswalk_status'] == 'unmapped').sum())} unmapped label(s)")
    if not shards:
        raise FileNotFoundError("no QC shards found; run `ovary-fm qc` first")

    # datasets deposited without annotations get their major type transferred from the
    # labelled, pooled shards (logistic regression on the training set's variable genes, the
    # T1 baseline); crosswalk_status = transferred, author_label stays 'unlabelled', and T1
    # never trains or tests on them
    transferred = transfer_labels(cfg, shards, manifest)
    for frame in review:
        hit = frame["dataset"].isin(transferred)
        frame.loc[hit, "crosswalk_status"] = "transferred"
        frame.loc[hit, "major_type"] = "(per cell, from the labelled datasets)"

    pooled = [a for a, held in shards.values() if not held]
    if not pooled:
        raise ValueError("every enabled dataset is held_out; nothing to pool")
    counts = pd.Series(dtype=int)
    for a in pooled:
        counts = counts.add(pd.Series(1, index=pd.Index(a.var_names.astype(str))), fill_value=0)
    if h.get("gene_space", "union") == "intersection":
        genes = counts[counts == len(pooled)].index
    else:
        genes = counts[counts >= int(h.get("min_datasets_per_gene", 1))].index
    genes = pd.Index(sorted(genes))
    print(f"  gene space ({h.get('gene_space', 'union')}): {len(genes):,} genes")

    pooled_shards, held_shards = [], []
    for name, (adata, held) in shards.items():
        new = reindex_genes(adata, genes)
        new.var[f"measured_in_{name}"] = new.var.pop("measured_here")
        (held_shards if held else pooled_shards).append((name, new))
    corpus = ad.concat([a for _, a in pooled_shards], join="inner", merge="first", uns_merge=None, index_unique=None)
    corpus.var = pd.DataFrame(index=genes)
    measured = np.zeros(len(genes), dtype=int)
    for name, a in pooled_shards:
        measured += a.var[f"measured_in_{name}"].to_numpy().astype(int)
        corpus.var[f"measured_in_{name}"] = a.var[f"measured_in_{name}"].to_numpy()
    corpus.var["n_datasets_measured"] = measured
    corpus.var["n_datasets_measured"] = counts.reindex(genes).fillna(0).astype(int).to_numpy()
    if any("gene_ids" in a.var for _, a in pooled_shards):
        ids = pd.Series("", index=genes)
        for _, a in pooled_shards:
            if "gene_ids" in a.var:
                have = a.var["gene_ids"].astype(str)
                ids[have != ""] = have[have != ""]
        corpus.var["gene_ids"] = ids.to_numpy()
    for col in corpus.obs.columns:
        if corpus.obs[col].dtype == object:
            corpus.obs[col] = corpus.obs[col].astype(str)
    corpus.uns["harmonise"] = {"gene_space": h.get("gene_space", "union"), "datasets": [n for n, _ in pooled_shards],
                               "held_out": [n for n, _ in held_shards]}
    out = corpus_path(cfg)
    corpus.write_h5ad(out)
    print(f"  corpus: {corpus.n_obs:,} cells x {corpus.n_vars:,} genes from {len(pooled_shards)} dataset(s) -> {out}")
    for name, a in held_shards:
        a.var = corpus.var.copy()
        for col in a.obs.columns:
            if a.obs[col].dtype == object:
                a.obs[col] = a.obs[col].astype(str)
        path = heldout_path(cfg, name)
        a.write_h5ad(path)
        print(f"  held out: {name}: {a.n_obs:,} cells -> {path}")

    pd.DataFrame(manifest).to_csv(tables / "manifest.csv", index=False)
    pd.concat(review, ignore_index=True).sort_values(["dataset", "n_cells"], ascending=[True, False]) \
        .to_csv(tables / "crosswalk_review.csv", index=False)
    composition = pd.crosstab(corpus.obs["dataset"], corpus.obs["major_type"])
    composition.to_csv(tables / "composition_by_dataset.csv")
    print(f"  -> {tables / 'manifest.csv'}, crosswalk_review.csv, composition_by_dataset.csv")
    return corpus
