"""The two baselines every T3 bundle is scored against, and the cell generator
they share with the stub model.

* b0 -- the control cells copied unchanged ("no effect").
* b1 -- retrieval and transfer: the same condition family's measured
  per-gene fold change in OTHER datasets (the nearest time point when both
  sides carry one), applied to this context's control cells and resampled
  as integer counts. This is the recipe that placed in both Virtual Cell
  Challenges and the bar the preregistration sets.
"""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp

FC_CLIP = (1 / 50.0, 50.0)


def fold_change(mean_cpm_group: np.ndarray, mean_cpm_control: np.ndarray, eps: float = 1e-9) -> np.ndarray:
    fc = (mean_cpm_group + eps) / (mean_cpm_control + eps)
    return np.clip(fc, *FC_CLIP)


def combine_fold_changes(fcs: list[np.ndarray]) -> np.ndarray:
    """Geometric mean across sources."""
    stack = np.log(np.vstack(fcs))
    return np.exp(stack.mean(axis=0))


def apply_fold_change(X_control, fc: np.ndarray, rng: np.random.Generator):
    """Integer counts with the control's own sparsity and depth structure: a
    gene whose fold change is 1 keeps each cell's count exactly; a changed gene
    is Poisson-resampled at count x fold change. Resampling only what changes
    keeps the transfer from paying an expression-error penalty on the genes it
    leaves alone (the sampling step is where the integer-count member is won
    or lost)."""
    X = sp.csr_matrix(X_control, dtype=np.float64)
    out = X.multiply(fc[np.newaxis, :]).tocsr()
    changed = np.isin(out.indices, np.where(fc != 1.0)[0])
    if changed.any():
        out.data[changed] = rng.poisson(out.data[changed])
    out.data = out.data.astype(np.float32)
    out.eliminate_zeros()
    return out


def b0_cells(X_control):
    return sp.csr_matrix(X_control, dtype=np.float32).copy()


def pick_sources(bundle: dict, candidates: list[dict]) -> list[dict]:
    """Other datasets, same family and compartment; the nearest time point when known."""
    same = [c for c in candidates
            if c["dataset"] != bundle["dataset"] and c["family"] == bundle["family"]
            and c["compartment"] == bundle["compartment"]]
    if not same:
        return []
    t = bundle.get("time_h")
    timed = [c for c in same if c.get("time_h") is not None and not np.isnan(c["time_h"])]
    if t is not None and not np.isnan(t) and timed:
        best = min(abs(c["time_h"] - t) for c in timed)
        return [c for c in timed if abs(c["time_h"] - t) == best]
    return same
