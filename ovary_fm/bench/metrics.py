"""The six members of the T3 score and the rescaling, in the `lite` scorer.

Implemented from the published vcc2026 specification (cell-eval2
``docs/vcc2026_metrics``) as the research memo records it; this is NOT
Arc's code and does not reproduce their jackknife correction or their
exact DE machinery. ``ovary-fm score --scorer cell-eval2`` wraps theirs.

Conventions kept from the specification:

* the DE table is a per-gene two-sided Wilcoxon rank-sum on CPM-normalised
  counts, genes tested only when their mean in the REAL control cells is
  above ``min_cpm_control`` (so every submission is scored on the same
  genes), log2 fold change on arithmetic means with a small epsilon,
  Benjamini-Hochberg within each bundle; the measured control is the
  comparison group on both sides;
* the pseudobulk for the expression and discrimination members is the
  group-summed counts scaled to ``target_sum`` and log1p'd;
* excluded genes (the condition-defining list) leave every member;
* ``pds_cosine`` ranks the true bundle among every bundle of the pool by
  cosine distance of pseudobulk deltas; a zero delta scores 0.5;
* ``de_direction_fidelity`` = k / max(n_pred, n_real); ``de_direction_reach``
  = the deepest prefix of the real DE list, ordered by the submission's
  confidence, with >= 90 % sign agreement, over n_real;
  ``de_sig_jaccard`` on the significant sets (empty union = 1);
  ``de_lfc_nmae`` over the real DE genes when there are >= 10 of them;
* the replicate is measured at HALF depth and this implementation carries no
  jackknife correction for it, so on a small effect the half-vs-half expression
  error can exceed the control copy's (b0 < r); Arc's expression member is
  noise-corrected, which is one reason to run cell-eval2 beside this scorer;
* every member is rescaled s = (u - b) / (r - b); the expression member is
  clamped to [0, 1] and the NMAE member floored at -6; nothing else is
  clipped.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats

MEMBERS = ["pds_cosine", "expr_mse", "de_direction_fidelity", "de_direction_reach", "de_sig_jaccard",
           "de_lfc_nmae"]


# --- pseudobulk and DE ------------------------------------------------------------

def pseudobulk_log(X, target_sum: float = 5e4) -> np.ndarray:
    s = np.asarray(sp.csr_matrix(X).sum(axis=0)).ravel().astype(np.float64)
    total = s.sum()
    return np.log1p(s / total * target_sum) if total > 0 else s


def cpm(X, target: float = 1e6):
    X = sp.csr_matrix(X, dtype=np.float64)
    rows = np.asarray(X.sum(axis=1)).ravel()
    rows[rows == 0] = 1.0
    return sp.diags(target / rows).dot(X).tocsr()


def bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    n = p.size
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.minimum(ranked, 1.0)
    return out


def de_table(X_group, X_control, genes, de_cfg: dict, tested_mask: np.ndarray | None = None,
             chunk: int = 2000) -> pd.DataFrame:
    """Wilcoxon rank-sum per gene, group vs control, on CPM; BH within the table."""
    target = float(de_cfg.get("normalise_to", 1e6))
    eps = float(de_cfg.get("lfc_eps", 1e-9))
    fdr = float(de_cfg.get("fdr", 0.05))
    G = cpm(X_group, target)
    C = cpm(X_control, target)
    mean_g = np.asarray(G.mean(axis=0)).ravel()
    mean_c = np.asarray(C.mean(axis=0)).ravel()
    if tested_mask is None:
        tested_mask = mean_c > float(de_cfg.get("min_cpm_control", 5.0))
    idx = np.where(tested_mask)[0]
    pvals = np.ones(len(genes))
    for start in range(0, len(idx), chunk):
        cols = idx[start:start + chunk]
        a = G[:, cols].toarray()
        b = C[:, cols].toarray()
        if a.shape[0] < 2 or b.shape[0] < 2:
            continue
        res = stats.mannwhitneyu(a, b, axis=0, alternative="two-sided", method="asymptotic")
        p = np.asarray(res.pvalue, dtype=float)
        # a gene constant in both groups has no test: p = 1
        constant = (a.max(axis=0) == a.min(axis=0)) & (b.max(axis=0) == b.min(axis=0)) & (a.max(axis=0) == b.max(axis=0))
        p[constant] = 1.0
        pvals[cols] = p
    padj = np.ones(len(genes))
    padj[idx] = bh(pvals[idx])
    lfc = np.log2((mean_g + eps) / (mean_c + eps))
    frame = pd.DataFrame({"gene": np.asarray(genes, dtype=str), "mean_cpm_group": mean_g, "mean_cpm_control": mean_c,
                          "lfc": lfc, "pval": pvals, "padj": padj, "tested": tested_mask})
    frame["sig"] = frame["tested"] & (frame["padj"] < fdr)
    return frame.set_index("gene")


# --- members ------------------------------------------------------------------------

def _restrict(de: pd.DataFrame, keep: np.ndarray) -> pd.DataFrame:
    return de[keep]


def de_members(pred_de: pd.DataFrame, real_de: pd.DataFrame, keep: np.ndarray, purity: float = 0.9,
               nmae_min_genes: int = 10) -> dict:
    """fidelity, reach, jaccard, nmae from two DE tables over the kept genes."""
    p = pred_de[keep]
    r = real_de[keep]
    P = set(p.index[p["sig"].to_numpy()])
    R = set(r.index[r["sig"].to_numpy()])
    inter = P & R
    union = P | R
    same_sign = sum(1 for g in inter if np.sign(p.at[g, "lfc"]) == np.sign(r.at[g, "lfc"]) and p.at[g, "lfc"] != 0)
    denom = max(len(P), len(R))
    fidelity = same_sign / denom if denom else np.nan
    jaccard = len(inter) / len(union) if union else 1.0
    n_real = len(R)
    if n_real:
        real_genes = sorted(R)
        conf = p.loc[real_genes, "lfc"].abs().to_numpy()
        order = np.lexsort((np.asarray(real_genes), -conf))
        ordered = np.asarray(real_genes)[order]
        agree = ((np.sign(p.loc[ordered, "lfc"].to_numpy()) == np.sign(r.loc[ordered, "lfc"].to_numpy()))
                 & (p.loc[ordered, "lfc"].to_numpy() != 0)).astype(float)
        cummean = np.cumsum(agree) / (np.arange(n_real) + 1)
        ok = np.where(cummean >= purity)[0]
        k_star = int(ok.max() + 1) if ok.size else 0
        reach = k_star / n_real
        if n_real >= nmae_min_genes:
            diff = np.abs(p.loc[real_genes, "lfc"].to_numpy() - r.loc[real_genes, "lfc"].to_numpy()).sum()
            scale = np.abs(r.loc[real_genes, "lfc"].to_numpy()).sum()
            nmae = diff / scale if scale > 0 else np.nan
        else:
            nmae = np.nan
    else:
        reach, nmae = np.nan, np.nan
    return {"de_direction_fidelity": fidelity, "de_direction_reach": reach, "de_sig_jaccard": jaccard,
            "de_lfc_nmae": nmae, "n_pred_de": len(P), "n_real_de": n_real, "n_overlap_de": len(inter)}


def expr_mse(pb_pred: np.ndarray, pb_real: np.ndarray, keep: np.ndarray) -> float:
    d = pb_pred[keep] - pb_real[keep]
    return float(np.mean(d * d)) if d.size else np.nan


def cosine_distance(a: np.ndarray, b: np.ndarray) -> float:
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    if na == 0 or nb == 0:
        return np.nan
    return float(1.0 - np.dot(a, b) / (na * nb))


def pds_cosine(delta_pred: np.ndarray, true_id: str, real_deltas: dict, keep: np.ndarray) -> float:
    """Mid-rank of the true bundle among the pool by cosine distance; zero delta -> 0.5."""
    ids = list(real_deltas)
    if len(ids) < 2 or true_id not in real_deltas:
        return np.nan
    if not np.any(delta_pred[keep]):
        return 0.5
    dists = {i: cosine_distance(delta_pred[keep], real_deltas[i][keep]) for i in ids}
    if any(np.isnan(v) for v in dists.values()):
        return 0.5
    d_true = dists[true_id]
    closer = sum(1 for i in ids if i != true_id and dists[i] < d_true)
    tied = sum(1 for i in ids if i != true_id and dists[i] == d_true)
    k = closer + 0.5 * tied
    return 1.0 - k / (len(ids) - 1)


# --- rescaling ------------------------------------------------------------------------

def rescale(member: str, u: float, b: float, r: float, clamp: dict | None = None) -> float:
    if any(v is None or (isinstance(v, float) and np.isnan(v)) for v in (u, b, r)):
        return np.nan
    if r == b:
        return np.nan
    s = (u - b) / (r - b)
    clamp = clamp or {}
    if member == "expr_mse" and "expr_mse" in clamp:
        lo, hi = clamp["expr_mse"]
        s = float(min(max(s, lo), hi))
    if member == "de_lfc_nmae" and "de_lfc_nmae_floor" in clamp:
        s = float(max(s, clamp["de_lfc_nmae_floor"]))
    return float(s)
