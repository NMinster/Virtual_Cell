"""T1: leave-one-atlas-out label transfer at the major level.

Default classifier: scikit-learn logistic regression on the 2,000 most
variable genes of the training set (log-normalised, scaled), class
balanced, the training set subsampled per class. With ``--rep`` a
kNN transfer inside a given embedding (an fm or arm-C
``embedding.parquet``) replaces it, so a representation is judged on the
same split and the same metrics as the gene-space baseline.

Scores: accuracy, macro-F1 over the test types present, per-type F1, and
rare-type F1 (the mean over the benchmark's rare types present).
Unmapped cells (crosswalk_status == unmapped) are excluded from both
sides; held-out shards are scored as test sets only.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from ..config import slug, tables_dir, write_provenance
from ..io import corpus_path, datasets, heldout_path


def _features(train, test, n_top_genes: int, seed: int):
    import scanpy as sc

    tr = train.copy()
    te = test.copy()
    for a in (tr, te):
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
    n_top = min(n_top_genes, tr.n_vars - 1)
    try:
        sc.pp.highly_variable_genes(tr, n_top_genes=n_top, flavor="seurat")
        hv = tr.var["highly_variable"].to_numpy()
    except Exception:
        hv = np.ones(tr.n_vars, dtype=bool)
    Xtr = tr[:, hv].X.toarray() if hasattr(tr.X, "toarray") else tr[:, hv].X
    Xte = te[:, hv].X.toarray() if hasattr(te.X, "toarray") else te[:, hv].X
    mu = Xtr.mean(axis=0)
    sd = Xtr.std(axis=0)
    sd[sd == 0] = 1.0
    Xtr = np.clip((Xtr - mu) / sd, -10, 10)
    Xte = np.clip((Xte - mu) / sd, -10, 10)
    return Xtr, Xte


def _subsample(labels: np.ndarray, max_per_class: int, rng) -> np.ndarray:
    keep = []
    for lab in np.unique(labels):
        idx = np.where(labels == lab)[0]
        if len(idx) > max_per_class:
            idx = rng.choice(idx, max_per_class, replace=False)
        keep.append(idx)
    return np.sort(np.concatenate(keep))


def _fit_predict(train, test, t1cfg: dict, seed: int, rep: pd.DataFrame | None):
    from sklearn.linear_model import LogisticRegression
    from sklearn.neighbors import KNeighborsClassifier

    rng = np.random.default_rng(seed)
    ytr = train.obs["major_type"].astype(str).to_numpy()
    sub = _subsample(ytr, int(t1cfg.get("max_cells_per_class", 2000)), rng)
    if rep is not None:
        Xtr = rep.reindex(train.obs_names[sub]).to_numpy(dtype=float)
        Xte = rep.reindex(test.obs_names).to_numpy(dtype=float)
        if np.isnan(Xtr).any() or np.isnan(Xte).any():
            raise ValueError("the embedding does not cover every cell of the train/test split")
        clf = KNeighborsClassifier(n_neighbors=min(int(t1cfg.get("knn", 15)), len(sub)))
    else:
        Xtr, Xte = _features(train[sub], test, int(t1cfg.get("n_top_genes", 2000)), seed)
        clf = LogisticRegression(C=float(t1cfg.get("C", 0.1)), max_iter=int(t1cfg.get("max_iter", 300)),
                                 class_weight="balanced")
    clf.fit(Xtr, ytr[sub])
    return clf.predict(Xte)


def _metrics(y_true, y_pred, rare: list[str], min_cells: int) -> tuple[dict, pd.DataFrame]:
    from sklearn.metrics import f1_score

    types = [t for t, n in pd.Series(y_true).value_counts().items() if n >= min_cells]
    per_type = pd.DataFrame({"major_type": types,
                             "n_cells": [int((y_true == t).sum()) for t in types],
                             "f1": [f1_score(y_true == t, y_pred == t) for t in types],
                             "recall": [float(((y_true == t) & (y_pred == t)).sum() / max((y_true == t).sum(), 1)) for t in types]})
    present_rare = [t for t in rare if t in types]
    summary = {"accuracy": float((y_true == y_pred).mean()),
               "macro_f1": float(per_type["f1"].mean()) if len(per_type) else np.nan,
               "rare_f1": float(per_type.set_index("major_type").loc[present_rare, "f1"].mean()) if present_rare else np.nan,
               "n_types_scored": len(types), "rare_types_present": ";".join(present_rare)}
    return summary, per_type


def run(cfg, rep: str | None = None, name: str | None = None) -> pd.DataFrame:
    write_provenance(cfg, "t1", extra={"rep": rep, "name": name})
    import anndata as ad

    t1cfg = cfg["benchmark"]["t1"]
    seed = int(cfg.get("random_state", 0))
    name = name or ("logreg" if rep is None else f"rep_{Path(rep).parent.name}")
    corpus = ad.read_h5ad(corpus_path(cfg))
    corpus = corpus[corpus.obs["crosswalk_status"].astype(str) != "unmapped"].copy()
    corpus = corpus[corpus.obs["author_label"].astype(str) != "unlabelled"].copy()
    rep_frame = None
    if rep:
        rep_frame = pd.read_parquet(rep)
        rep_frame.index = rep_frame.index.astype(str)
    labelled = sorted(corpus.obs["dataset"].astype(str).unique())
    held = t1cfg.get("held_out", "all")
    held = labelled if held == "all" else [d for d in held if d in labelled]
    rare = list(t1cfg.get("rare_types", []))
    min_cells = int(t1cfg.get("min_cells_per_type", 20))
    rows, per_types = [], []
    for dset in held:
        test = corpus[corpus.obs["dataset"].astype(str) == dset]
        train = corpus[corpus.obs["dataset"].astype(str) != dset]
        if train.n_obs == 0 or test.n_obs == 0:
            continue
        pred = _fit_predict(train, test, t1cfg, seed, rep_frame)
        summary, per_type = _metrics(test.obs["major_type"].astype(str).to_numpy(), pred, rare, min_cells)
        per_type["held_out"] = dset
        rows.append({"held_out": dset, "role": "leave_one_out", "n_train": train.n_obs, "n_test": test.n_obs, **summary})
        per_types.append(per_type)
        print(f"  {dset:<24} acc {summary['accuracy']:.3f}  macro-F1 {summary['macro_f1']:.3f}  "
              f"rare-F1 {summary['rare_f1']:.3f}  ({summary['n_types_scored']} types)")
    for dset in t1cfg.get("test_only") or []:
        path = heldout_path(cfg, dset)
        if not path.exists():
            continue
        test = ad.read_h5ad(path)
        test = test[test.obs["crosswalk_status"].astype(str) != "unmapped"].copy()
        if rep_frame is not None and not test.obs_names.isin(rep_frame.index).all():
            print(f"  {dset}: the embedding does not cover the held-out shard; skipped")
            continue
        pred = _fit_predict(corpus, test, t1cfg, seed, rep_frame)
        summary, per_type = _metrics(test.obs["major_type"].astype(str).to_numpy(), pred, rare, min_cells)
        per_type["held_out"] = dset
        rows.append({"held_out": dset, "role": "test_only", "n_train": corpus.n_obs, "n_test": test.n_obs, **summary})
        per_types.append(per_type)
        print(f"  {dset:<24} acc {summary['accuracy']:.3f}  macro-F1 {summary['macro_f1']:.3f}  "
              f"rare-F1 {summary['rare_f1']:.3f}  (test only)")
    frame = pd.DataFrame(rows)
    frame["model"] = name
    tables = tables_dir(cfg, "bench")
    frame.to_csv(tables / f"t1_{slug(name)}.csv", index=False)
    if per_types:
        pd.concat(per_types, ignore_index=True).to_csv(tables / f"t1_{slug(name)}_per_type.csv", index=False)
    print(f"  -> {tables / f't1_{slug(name)}.csv'}")
    return frame
