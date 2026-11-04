"""``ovary-fm bench`` builds the T3 reference bundles; ``ovary-fm score`` scores
a prediction object against them.

A bundle is one (context, perturbed condition, compartment): the context's
control cells of that compartment and the real perturbed cells. For each
bundle the builder stores the cells and their split halves (cross-animal
when both groups have >= 2 animals, else by cells, and it says which),
the real DE table, the two baselines b0 and b1 as scored submissions, and
the replicate r (half A scored against half B, five seeds). The ceilings
table lists all of it, including the bundles that fall below the gates.

A prediction object is an .h5ad of integer counts with obs ``bundle_id``;
``score`` computes the six members for every bundle it covers, rescales
each against the bundle's baseline and replicate, and aggregates: bundle
score = mean of the available members, context score = mean over its
bundles, overall = mean over contexts.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from ..config import results_dir, slug, tables_dir, write_provenance
from ..io import corpus_path, datasets, heldout_path
from . import baselines
from .metrics import MEMBERS, de_members, de_table, expr_mse, pds_cosine, pseudobulk_log, rescale


# --- loading the corpus by dataset -------------------------------------------------------

def load_shards(cfg) -> dict:
    """dataset name -> AnnData (the corpus split by dataset, plus the held-out shards)."""
    import anndata as ad

    shards = {}
    path = corpus_path(cfg)
    if path.exists():
        corpus = ad.read_h5ad(path)
        for name in corpus.obs["dataset"].unique():
            shards[str(name)] = corpus[corpus.obs["dataset"] == name].copy()
    for entry in datasets(cfg):
        if entry.get("held_out") and heldout_path(cfg, entry["name"]).exists():
            shards[entry["name"]] = ad.read_h5ad(heldout_path(cfg, entry["name"]))
    if not shards:
        raise FileNotFoundError(f"no corpus at {path}; run `ovary-fm harmonise` first")
    return shards


def excluded_genes(t3: dict, family: str, context: dict) -> set:
    ex = t3.get("exclude_genes") or {}
    out = set(ex.get("global") or []) | set(ex.get(family) or []) | set(context.get("exclude") or [])
    return out


# --- splits ----------------------------------------------------------------------------------

def split_halves(animals: np.ndarray, rng: np.random.Generator, min_animals: int) -> tuple[np.ndarray, str]:
    """A/B per cell: by animal when there are enough animals, else by cell."""
    unique = np.unique(animals)
    if len(unique) >= min_animals:
        order = rng.permutation(unique)
        a_set = set(order[: len(order) // 2]) if len(order) % 2 == 0 else set(order[: (len(order) + 1) // 2])
        # alternate so the halves are balanced in animal count; with an odd
        # count the extra animal goes to A on even seeds and B on odd ones
        half = np.where(np.isin(animals, list(a_set)), "A", "B")
        return half, "animals"
    half = np.where(rng.random(len(animals)) < 0.5, "A", "B")
    return half, "cells"


# --- build ------------------------------------------------------------------------------------

def build(cfg):
    write_provenance(cfg, "bench")
    import anndata as ad

    bench = cfg["benchmark"]
    t3 = bench["t3"]
    de_cfg = t3["de"]
    shards = load_shards(cfg)
    genes = None
    out_dir = results_dir(cfg, "bench", "bundles")
    de_dir = tables_dir(cfg, "bench", "de")
    tables = tables_dir(cfg, "bench")
    key = t3.get("compartment_key", "major_type")
    min_hard = int(t3.get("min_cells_hard", 50))
    min_animals = int(t3.get("min_animals_for_cross_animal_halves", 2))
    n_splits = int(t3.get("n_splits", 5))
    base_seed = int(t3.get("base_seed", 0))
    purity = float(t3.get("reach_purity", 0.9))
    nmae_min = int(de_cfg.get("nmae_min_genes", 10))
    pds_target = float(t3.get("pds", {}).get("target_sum", 5e4))

    bundles, cells_by_context, skipped = [], {}, []
    for context in t3["contexts"]:
        name, dset = context["name"], context["dataset"]
        if dset not in shards:
            skipped.append({"context": name, "reason": f"dataset {dset} not in the corpus or held-out shards"})
            print(f"  - {name}: dataset {dset} is not available; skipped")
            continue
        adata = shards[dset]
        if genes is None:
            genes = pd.Index(adata.var_names.astype(str))
        elif not genes.equals(pd.Index(adata.var_names.astype(str))):
            raise ValueError(f"{dset} is not in the corpus gene space; rerun `ovary-fm harmonise`")
        obs = adata.obs
        exclude = excluded_genes(t3, context["family"], context)
        keep_mask = ~genes.isin(exclude)
        pieces = []
        for compartment in t3["compartments"]:
            ctrl_idx = np.where((obs["condition"].astype(str) == str(context["control"])).to_numpy()
                                & (obs[key].astype(str) == compartment).to_numpy())[0]
            if len(ctrl_idx) == 0:
                continue
            ctrl = adata[ctrl_idx].copy()
            ctrl.obs["role"] = "control"
            ctrl.obs["compartment"] = compartment
            ctrl.obs["bundle_id"] = f"{name}|control|{compartment}"
            for seed in range(n_splits):
                rng = np.random.default_rng(base_seed + seed)
                ctrl.obs[f"split_{seed}"], ctrl_level = split_halves(ctrl.obs["animal"].astype(str).to_numpy(), rng, min_animals)
            pieces.append(ctrl)
            for p in context["perturbed"]:
                pert_idx = np.where((obs["condition"].astype(str) == str(p)).to_numpy()
                                    & (obs[key].astype(str) == compartment).to_numpy())[0]
                bid = f"{name}|{p}|{compartment}"
                record = {"bundle_id": bid, "context": name, "dataset": dset, "family": context["family"],
                          "condition": p, "control": context["control"], "compartment": compartment,
                          "n_control": int(len(ctrl_idx)), "n_pert": int(len(pert_idx)),
                          "n_animals_control": int(ctrl.obs["animal"].nunique()),
                          "n_animals_pert": int(obs["animal"].iloc[pert_idx].nunique()) if len(pert_idx) else 0,
                          "held_out": bool(context.get("held_out", False)),
                          "n_excluded_genes": int((~keep_mask).sum())}
                time_h = obs["condition_time_h"].iloc[pert_idx] if len(pert_idx) else pd.Series(dtype=float)
                record["time_h"] = float(time_h.iloc[0]) if len(time_h) and pd.notna(time_h.iloc[0]) else np.nan
                if len(pert_idx) < min_hard or len(ctrl_idx) < min_hard:
                    record["scorable"] = False
                    record["reason"] = f"fewer than {min_hard} cells (control {len(ctrl_idx)}, perturbed {len(pert_idx)})"
                    bundles.append(record)
                    continue
                pert = adata[pert_idx].copy()
                pert.obs["role"] = "perturbed"
                pert.obs["compartment"] = compartment
                pert.obs["bundle_id"] = bid
                levels = []
                for seed in range(n_splits):
                    rng = np.random.default_rng(base_seed + seed + 1000)
                    pert.obs[f"split_{seed}"], level = split_halves(pert.obs["animal"].astype(str).to_numpy(), rng, min_animals)
                    levels.append(level)
                record["replicate_level"] = "animals" if all(l == "animals" for l in levels) and ctrl_level == "animals" else "cells"
                record["scorable"] = True
                record["reason"] = ""
                record["_ctrl"] = ctrl
                record["_pert"] = pert
                record["_keep"] = keep_mask
                bundles.append(record)
                pieces.append(pert)
        if pieces:
            cells_by_context[name] = ad.concat(pieces, join="inner", merge="first", index_unique=None)
    if genes is None:
        raise ValueError("no context could be built; check benchmark.contexts against the manifest")

    # real DE, pseudobulks, b0, r --------------------------------------------------------------
    scorable = [b for b in bundles if b.get("scorable")]
    print(f"  {len(scorable)} scorable bundle(s) of {len(bundles)}; {len(skipped)} context(s) skipped")
    for b in scorable:
        ctrl, pert, keep = b["_ctrl"], b["_pert"], b["_keep"]
        real = de_table(pert.X, ctrl.X, genes, de_cfg)
        b["_real_de"] = real
        b["_pb_ctrl"] = pseudobulk_log(ctrl.X, pds_target)
        b["_pb_real"] = pseudobulk_log(pert.X, pds_target)
        # the fold change a source hands on: its SIGNIFICANT genes only, 1.0 elsewhere (the
        # retrieval-and-transfer recipe pools tested summaries, not every gene's noisy ratio)
        fc_all = baselines.fold_change(real["mean_cpm_group"].to_numpy(), real["mean_cpm_control"].to_numpy())
        b["_fc_real"] = np.where(real["sig"].to_numpy(), fc_all, 1.0)
        b["n_de_real"] = int(real["sig"].sum())
        b["n_de_real_kept"] = int((real["sig"].to_numpy() & keep).sum())
        b["gate_nmae"] = bool(b["n_de_real_kept"] >= nmae_min)
        real.to_csv(de_dir / f"{slug(b['bundle_id'])}.csv")
    # b1 needs every bundle's real fold change first
    for b in scorable:
        sources = baselines.pick_sources(b, scorable)
        b["b1_sources"] = [s["bundle_id"] for s in sources]
        b["n_b1_sources"] = len(sources)
        b["_fc_b1"] = baselines.combine_fold_changes([s["_fc_real"] for s in sources]) if sources else None

    # the pools for PDS: every bundle of a compartment
    pools = {}
    for b in scorable:
        pools.setdefault(b["compartment"], []).append(b)

    rng = np.random.default_rng(base_seed)
    for b in scorable:
        ctrl, pert, keep, real = b["_ctrl"], b["_pert"], b["_keep"], b["_real_de"]
        pool = {p["bundle_id"]: p["_pb_real"] - p["_pb_ctrl"] for p in pools[b["compartment"]]}
        pool_keep = np.logical_and.reduce([p["_keep"] for p in pools[b["compartment"]]])
        # b0: the control copied
        X0 = baselines.b0_cells(ctrl.X)
        u0 = _members_for(X0, ctrl.X, genes, de_cfg, real, b, pool, pool_keep, purity, nmae_min, pds_target)
        # b1: transfer
        if b["_fc_b1"] is not None:
            X1 = baselines.apply_fold_change(ctrl.X, b["_fc_b1"], rng)
            u1 = _members_for(X1, ctrl.X, genes, de_cfg, real, b, pool, pool_keep, purity, nmae_min, pds_target)
        else:
            u1 = {m: np.nan for m in MEMBERS}
        # r: half A as the submission, half B as the reference, per seed
        rs = {m: [] for m in MEMBERS}
        for seed in range(n_splits):
            col = f"split_{seed}"
            cA = ctrl[ctrl.obs[col] == "A"]
            cB = ctrl[ctrl.obs[col] == "B"]
            pA = pert[pert.obs[col] == "A"]
            pB = pert[pert.obs[col] == "B"]
            if min(cA.n_obs, cB.n_obs, pA.n_obs, pB.n_obs) < 5:
                continue
            tested = real["tested"].to_numpy()
            deA = de_table(pA.X, cA.X, genes, de_cfg, tested_mask=tested)
            deB = de_table(pB.X, cB.X, genes, de_cfg, tested_mask=tested)
            pbA, pbB = pseudobulk_log(pA.X, pds_target), pseudobulk_log(pB.X, pds_target)
            pb_cA, pb_cB = pseudobulk_log(cA.X, pds_target), pseudobulk_log(cB.X, pds_target)
            vals = de_members(deA, deB, keep, purity, nmae_min)
            vals["expr_mse"] = expr_mse(pbA, pbB, keep)
            # the pool for the replicate's PDS is the B halves of every bundle in the compartment
            poolB = {}
            for p in pools[b["compartment"]]:
                pc = p["_ctrl"][p["_ctrl"].obs[col] == "B"]
                pp = p["_pert"][p["_pert"].obs[col] == "B"]
                if pc.n_obs and pp.n_obs:
                    poolB[p["bundle_id"]] = pseudobulk_log(pp.X, pds_target) - pseudobulk_log(pc.X, pds_target)
            vals["pds_cosine"] = pds_cosine(pbA - pb_cA, b["bundle_id"], poolB, pool_keep)
            for m in MEMBERS:
                rs[m].append(vals[m])
        for m in MEMBERS:
            b[f"b0_{m}"] = u0[m]
            b[f"b1_{m}"] = u1[m]
            arr = np.asarray([v for v in rs[m] if v is not None and not np.isnan(v)], dtype=float)
            b[f"r_{m}"] = float(arr.mean()) if arr.size else np.nan
            b[f"r_{m}_sd"] = float(arr.std(ddof=0)) if arr.size else np.nan
            b[f"r_{m}_n"] = int(arr.size)
        print(f"  {b['bundle_id']}: control {b['n_control']} / perturbed {b['n_pert']} cells, "
              f"{b['n_animals_control']}+{b['n_animals_pert']} animals ({b['replicate_level']}), "
              f"{b['n_de_real']} DE genes; b1 from {b['n_b1_sources']} source(s)")

    # persist ---------------------------------------------------------------------------------
    for name, cells in cells_by_context.items():
        for col in cells.obs.columns:
            if cells.obs[col].dtype == object:
                cells.obs[col] = cells.obs[col].astype(str)
        cells.write_h5ad(out_dir / f"{slug(name)}.h5ad")
    fc_frame = pd.DataFrame({b["bundle_id"]: b["_fc_b1"] for b in scorable if b["_fc_b1"] is not None},
                            index=genes)
    fc_frame.to_parquet(out_dir / "b1_fold_changes.parquet")
    public = [{k: v for k, v in b.items() if not k.startswith("_")} for b in bundles]
    reference = {"genes": list(genes), "members": MEMBERS, "n_splits": n_splits, "base_seed": base_seed,
                 "bundles": public, "skipped_contexts": skipped,
                 "pools": {c: [p["bundle_id"] for p in ps] for c, ps in pools.items()},
                 "pool_keep": {c: [bool(v) for v in np.logical_and.reduce([p["_keep"] for p in ps])] for c, ps in pools.items()},
                 "keep": {b["bundle_id"]: [bool(v) for v in b["_keep"]] for b in scorable}}
    (results_dir(cfg, "bench") / "reference.json").write_text(json.dumps(reference, indent=1, default=_json_default),
                                                              encoding="utf-8")
    ceilings = pd.DataFrame(public)
    front = ["bundle_id", "context", "dataset", "family", "condition", "control", "compartment", "scorable", "reason",
             "n_control", "n_pert", "n_animals_control", "n_animals_pert", "replicate_level", "n_de_real",
             "n_de_real_kept", "gate_nmae", "n_b1_sources", "b1_sources", "held_out", "time_h", "n_excluded_genes"]
    cols = [c for c in front if c in ceilings] + [c for c in ceilings.columns if c not in front]
    ceilings = ceilings[cols]
    ceilings.to_csv(tables / "ceilings.csv", index=False)
    print(f"  -> {tables / 'ceilings.csv'} ({len(ceilings)} rows; {len(scorable)} scorable)")
    return ceilings


def _json_default(o):
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (pd.Series, pd.Index)):
        return o.tolist()
    return str(o)


def _members_for(X_pred, X_ctrl, genes, de_cfg, real, b, pool, pool_keep, purity, nmae_min, pds_target) -> dict:
    keep = b["_keep"]
    pred_de = de_table(X_pred, X_ctrl, genes, de_cfg, tested_mask=real["tested"].to_numpy())
    vals = de_members(pred_de, real, keep, purity, nmae_min)
    pb_pred = pseudobulk_log(X_pred, pds_target)
    vals["expr_mse"] = expr_mse(pb_pred, b["_pb_real"], keep)
    vals["pds_cosine"] = pds_cosine(pb_pred - b["_pb_ctrl"], b["bundle_id"], pool, pool_keep)
    return vals


# --- load the bundles back ----------------------------------------------------------------

def load_reference(cfg) -> dict:
    path = results_dir(cfg, "bench") / "reference.json"
    if not path.exists():
        raise FileNotFoundError(f"{path} is missing; run `ovary-fm bench` first")
    return json.loads(path.read_text(encoding="utf-8"))


def load_bundle_cells(cfg, reference: dict) -> dict:
    """bundle_id -> {'ctrl': AnnData, 'pert': AnnData} for the scorable bundles."""
    import anndata as ad

    folder = results_dir(cfg, "bench", "bundles")
    by_context = {}
    out = {}
    for b in reference["bundles"]:
        if not b.get("scorable"):
            continue
        ctx = b["context"]
        if ctx not in by_context:
            by_context[ctx] = ad.read_h5ad(folder / f"{slug(ctx)}.h5ad")
        cells = by_context[ctx]
        ctrl = cells[cells.obs["bundle_id"] == f"{ctx}|control|{b['compartment']}"]
        pert = cells[cells.obs["bundle_id"] == b["bundle_id"]]
        out[b["bundle_id"]] = {"ctrl": ctrl, "pert": pert, "meta": b}
    return out


# --- score -------------------------------------------------------------------------------------

def score(cfg, pred_path: Path, name: str, scorer: str = "lite") -> pd.DataFrame:
    write_provenance(cfg, "score", extra={"pred": str(pred_path), "name": name, "scorer": scorer})
    if scorer == "cell-eval2":
        return _score_cell_eval2(cfg, pred_path, name)
    import anndata as ad

    bench = cfg["benchmark"]
    t3 = bench["t3"]
    de_cfg = t3["de"]
    purity = float(t3.get("reach_purity", 0.9))
    nmae_min = int(de_cfg.get("nmae_min_genes", 10))
    pds_target = float(t3.get("pds", {}).get("target_sum", 5e4))
    clamp = t3.get("clamp") or {}
    policy = t3.get("baseline", "b1_else_b0")

    reference = load_reference(cfg)
    genes = pd.Index(reference["genes"])
    bundles = load_bundle_cells(cfg, reference)
    pred = ad.read_h5ad(pred_path)
    if "bundle_id" not in pred.obs:
        raise ValueError("the prediction object needs obs['bundle_id']")
    # align genes: missing ones are zero, extra ones dropped
    pred_genes = pd.Index(pred.var_names.astype(str))
    common = pred_genes.intersection(genes)
    P = sp.csr_matrix((np.ones(len(common), dtype=np.float32), (pred_genes.get_indexer(common), genes.get_indexer(common))),
                      shape=(len(pred_genes), len(genes)))
    Xp = sp.csr_matrix(pred.X).dot(P).tocsr()
    if len(common) < len(genes):
        print(f"  {len(genes) - len(common)} corpus genes are absent from the prediction (zero)")

    # real deltas per pool, from the bundle cells
    pool_delta, pool_keep = {}, {}
    pb_cache = {}
    for bid, item in bundles.items():
        pb_c = pseudobulk_log(item["ctrl"].X, pds_target)
        pb_r = pseudobulk_log(item["pert"].X, pds_target)
        pb_cache[bid] = (pb_c, pb_r)
        pool_delta.setdefault(item["meta"]["compartment"], {})[bid] = pb_r - pb_c
    for comp, ids in reference["pools"].items():
        pool_keep[comp] = np.asarray(reference["pool_keep"][comp], dtype=bool)

    rows = []
    ids = pred.obs["bundle_id"].astype(str).to_numpy()
    for bid, item in bundles.items():
        mask = ids == bid
        meta = item["meta"]
        if not mask.any():
            continue
        X_pred = Xp[mask]
        keep = np.asarray(reference["keep"][bid], dtype=bool)
        real = de_table(item["pert"].X, item["ctrl"].X, genes, de_cfg)
        pred_de = de_table(X_pred, item["ctrl"].X, genes, de_cfg, tested_mask=real["tested"].to_numpy())
        vals = de_members(pred_de, real, keep, purity, nmae_min)
        pb_c, pb_r = pb_cache[bid]
        pb_p = pseudobulk_log(X_pred, pds_target)
        vals["expr_mse"] = expr_mse(pb_p, pb_r, keep)
        comp = meta["compartment"]
        vals["pds_cosine"] = pds_cosine(pb_p - pb_c, bid, pool_delta[comp], pool_keep[comp])
        row = {"bundle_id": bid, "context": meta["context"], "dataset": meta["dataset"], "family": meta["family"],
               "condition": meta["condition"], "compartment": comp, "n_pred_cells": int(mask.sum()),
               "n_pert": meta["n_pert"], "n_control": meta["n_control"], "replicate_level": meta.get("replicate_level"),
               "gate_nmae": meta.get("gate_nmae"), "n_pred_de": vals["n_pred_de"], "n_real_de": vals["n_real_de"]}
        has_b1 = meta.get("n_b1_sources", 0) > 0 and not _isnan(meta.get("b1_expr_mse"))
        use = "b1" if (policy == "b1_else_b0" and has_b1) else "b0"
        row["baseline_used"] = use
        ss = []
        for m in MEMBERS:
            u = vals[m]
            row[f"u_{m}"] = u
            s_b0 = rescale(m, u, meta.get(f"b0_{m}"), meta.get(f"r_{m}"), clamp)
            s_b1 = rescale(m, u, meta.get(f"b1_{m}"), meta.get(f"r_{m}"), clamp) if has_b1 else np.nan
            row[f"s_b0_{m}"] = s_b0
            row[f"s_b1_{m}"] = s_b1
            s = s_b1 if use == "b1" else s_b0
            row[f"s_{m}"] = s
            if not _isnan(s):
                ss.append(s)
        row["n_members"] = len(ss)
        row["members_positive"] = int(sum(1 for s in ss if s > 0))
        row["bundle_score"] = float(np.mean(ss)) if ss else np.nan
        rows.append(row)
    if not rows:
        raise ValueError("the prediction covers none of the scorable bundles")
    frame = pd.DataFrame(rows)
    tables = tables_dir(cfg, "bench")
    frame.to_csv(tables / f"scores_{slug(name)}.csv", index=False)
    summary = frame.groupby("context", as_index=False).agg(n_bundles=("bundle_id", "size"),
                                                           context_score=("bundle_score", "mean"),
                                                           members_positive=("members_positive", "mean"))
    summary["model"] = name
    overall = float(summary["context_score"].mean())
    summary.loc[len(summary)] = {"context": "OVERALL", "n_bundles": int(frame.shape[0]), "context_score": overall,
                                 "members_positive": float(frame["members_positive"].mean()), "model": name}
    summary.to_csv(tables / f"scores_{slug(name)}_summary.csv", index=False)
    print(f"  {name}: overall {overall:.3f} over {summary.shape[0] - 1} context(s), {frame.shape[0]} bundle(s)")
    for _, r in summary.iterrows():
        print(f"    {r['context']:<28} {r['context_score']:.3f}  (members > 0: {r['members_positive']:.1f} of 6)")
    print(f"  -> {tables / f'scores_{slug(name)}.csv'}")
    return frame


def _isnan(v) -> bool:
    return v is None or (isinstance(v, float) and np.isnan(v))


def _score_cell_eval2(cfg, pred_path: Path, name: str):
    """Arc's scorer. UNTESTED: written from the cell-eval2 README's command
    shape; check the flags against `cell-eval2 run --help` in the
    ovary_fm_arc env before trusting the table it writes."""
    import shutil
    import subprocess

    import anndata as ad

    reference = load_reference(cfg)
    bundles = load_bundle_cells(cfg, reference)
    pred = ad.read_h5ad(pred_path)
    work = results_dir(cfg, "bench", "cell_eval2", slug(name))
    env = "ovary_fm_arc"
    base = ["conda", "run", "--no-capture-output", "-n", env] if shutil.which("conda") else []
    rows = []
    for ctx in sorted({b["meta"]["context"] for b in bundles.values()}):
        items = {bid: it for bid, it in bundles.items() if it["meta"]["context"] == ctx}
        real = ad.concat([it["pert"] for it in items.values()] + [next(iter(items.values()))["ctrl"]],
                         join="inner", index_unique=None)
        real.obs["target_gene"] = np.where(real.obs["role"] == "control", "non-targeting",
                                           real.obs["bundle_id"].astype(str))
        sub = pred[pred.obs["bundle_id"].isin(list(items))].copy()
        sub.obs["target_gene"] = sub.obs["bundle_id"].astype(str)
        rp, sp_ = work / f"{slug(ctx)}_real.h5ad", work / f"{slug(ctx)}_pred.h5ad"
        real.write_h5ad(rp)
        sub.write_h5ad(sp_)
        cmd = base + ["cell-eval2", "run", "--preset", "vcc2026", "--pred", str(sp_), "--real", str(rp),
                      "--control-pert", "non-targeting", "--pert-col", "target_gene", "--outdir", str(work / slug(ctx))]
        print(f"  $ {' '.join(cmd)}")
        res = subprocess.run(cmd)
        rows.append({"context": ctx, "exit": res.returncode, "outdir": str(work / slug(ctx))})
    frame = pd.DataFrame(rows)
    frame.to_csv(tables_dir(cfg, "bench") / f"scores_{slug(name)}_cell_eval2_runs.csv", index=False)
    print("  read each outdir's tables; this adapter does not parse them yet")
    return frame
