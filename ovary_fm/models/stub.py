"""``ovary-fm predict --model b0|b1|stub``: a prediction object for every bundle.

The stub is the plumbing's model: the b1 transfer with log-fold noise, or
when no transfer source exists, a random set of genes shifted. It exists
so the scorer, the tables and the tests run end to end; it says nothing
about biology. A real model implements ``generate(ctrl, meta, genes) ->
csr integer counts`` and is registered in ``GENERATORS``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from ..config import results_dir, write_provenance
from ..bench import baselines
from ..bench.runner import load_bundle_cells, load_reference


def prediction_path(cfg, model: str) -> Path:
    return results_dir(cfg, "predictions") / f"{model}.h5ad"


def _stub_generate(cfg, ctrl, meta, genes, fc_b1, rng):
    s = cfg.get("stub", {})
    if fc_b1 is not None:
        fc = fc_b1 * np.power(2.0, rng.normal(0.0, float(s.get("noise_sd", 0.15)), size=len(genes)))
    else:
        fc = np.ones(len(genes))
        pick = rng.choice(len(genes), int(min(s.get("n_random_genes", 50), len(genes))), replace=False)
        fold = float(s.get("random_fold", 2.0))
        fc[pick] = np.where(rng.random(len(pick)) < 0.5, fold, 1.0 / fold)
    return baselines.apply_fold_change(ctrl.X, np.clip(fc, *baselines.FC_CLIP), rng)


def predict(cfg, model: str = "stub") -> Path:
    write_provenance(cfg, "predict", extra={"model": model})
    import anndata as ad

    reference = load_reference(cfg)
    genes = pd.Index(reference["genes"])
    bundles = load_bundle_cells(cfg, reference)
    fc_path = results_dir(cfg, "bench", "bundles") / "b1_fold_changes.parquet"
    fcs = pd.read_parquet(fc_path) if fc_path.exists() else pd.DataFrame(index=genes)
    rng = np.random.default_rng(int(cfg.get("random_state", 0)) + 7)
    pieces = []
    for bid, item in bundles.items():
        ctrl, meta = item["ctrl"], item["meta"]
        fc_b1 = fcs[bid].to_numpy(dtype=float) if bid in fcs else None
        if model == "b0":
            X = baselines.b0_cells(ctrl.X)
        elif model == "b1":
            if fc_b1 is None:
                continue
            X = baselines.apply_fold_change(ctrl.X, fc_b1, rng)
        elif model == "stub":
            X = _stub_generate(cfg, ctrl, meta, genes, fc_b1, rng)
        else:
            raise ValueError(f"unknown model {model!r}; b0 | b1 | stub")
        obs = pd.DataFrame({"bundle_id": bid, "context": meta["context"], "condition": meta["condition"],
                            "compartment": meta["compartment"], "model": model},
                           index=[f"{model}:{bid}:{i}" for i in range(X.shape[0])])
        pieces.append(ad.AnnData(X=sp.csr_matrix(X, dtype=np.float32), obs=obs, var=pd.DataFrame(index=genes)))
    if not pieces:
        raise ValueError(f"{model}: no bundle could be predicted (b1 needs a transfer source)")
    pred = ad.concat(pieces, join="inner", index_unique=None)
    pred.var = pd.DataFrame(index=genes)
    pred.uns["model"] = model
    out = prediction_path(cfg, model)
    pred.write_h5ad(out)
    print(f"  {model}: {pred.n_obs:,} cells over {len(pieces)} bundle(s) -> {out}")
    return out
