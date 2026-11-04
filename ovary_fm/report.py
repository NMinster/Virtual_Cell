"""``ovary-fm report``: the first-results page, built from every table that exists."""

from __future__ import annotations

from datetime import datetime, timezone

import numpy as np
import pandas as pd

from .config import git_revision, results_dir, tables_dir, write_provenance


def md_table(frame: pd.DataFrame, max_rows: int = 60, floatfmt: str = "{:.3f}") -> str:
    if frame is None or frame.empty:
        return "_(empty)_\n"
    frame = frame.head(max_rows).copy()
    for col in frame.columns:
        if frame[col].dtype.kind == "f":
            frame[col] = frame[col].map(lambda v: "" if pd.isna(v) else floatfmt.format(v))
    cols = [str(c) for c in frame.columns]
    lines = ["| " + " | ".join(cols) + " |", "| " + " | ".join("---" for _ in cols) + " |"]
    for _, row in frame.iterrows():
        lines.append("| " + " | ".join(str(v) for v in row.tolist()) + " |")
    return "\n".join(lines) + "\n"


def _read(path):
    return pd.read_csv(path) if path.exists() else None


def run(cfg):
    write_provenance(cfg, "report")
    tables = tables_dir(cfg)
    out = results_dir(cfg, "report") / "first_results.md"
    parts = [f"# First results\n\n{datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')} · "
             f"git {git_revision() or '?'} · results root `{cfg['paths']['results']}`\n"]

    qc = _read(tables / "qc" / "qc_summary.csv")
    if qc is not None:
        parts.append("## QC per dataset\n\n" + md_table(qc[[c for c in ("dataset", "species", "suspension_type", "cells_in",
                                                                     "cells_out", "n_animals", "conditions", "n_author_labels",
                                                                     "median_counts", "median_genes") if c in qc]]))
    manifest = _read(tables / "harmonise" / "manifest.csv")
    if manifest is not None:
        parts.append("## The corpus\n\n" + md_table(manifest))
    review = _read(tables / "harmonise" / "crosswalk_review.csv")
    if review is not None:
        unmapped = review[review["crosswalk_status"] == "unmapped"]
        parts.append(f"## Crosswalk\n\n{len(review)} author labels across the datasets; "
                     f"{len(unmapped)} unmapped (kept as text, out of the scoring).\n\n" + md_table(unmapped))
    ceilings = _read(tables / "bench" / "ceilings.csv")
    if ceilings is not None:
        show = [c for c in ("bundle_id", "scorable", "reason", "n_control", "n_pert", "n_animals_control", "n_animals_pert",
                            "replicate_level", "n_de_real", "n_de_real_kept", "gate_nmae", "n_b1_sources") if c in ceilings]
        parts.append("## T3 bundles and ceilings\n\n" + md_table(ceilings[show]))
        members = [c for c in ceilings.columns if c.startswith("r_") and not c.endswith(("_sd", "_n"))]
        if members:
            rc = ceilings[ceilings.get("scorable", True) == True]
            cols = ["bundle_id"]
            for m in members:
                cols += [c for c in (m.replace("r_", "b0_"), m.replace("r_", "b1_"), m) if c in rc]
            parts.append("### Baselines and replicate per member (raw u)\n\n" + md_table(rc[cols]))
    scores = sorted((tables / "bench").glob("scores_*_summary.csv")) if (tables / "bench").exists() else []
    if scores:
        frames = [pd.read_csv(p) for p in scores]
        allscores = pd.concat(frames, ignore_index=True)
        wide = allscores.pivot_table(index="context", columns="model", values="context_score")
        parts.append("## T3 scores (s against b1 where a source exists, else b0; 0 = baseline, 1 = replicate)\n\n"
                     + md_table(wide.reset_index()))
    t1s = sorted((tables / "bench").glob("t1_*.csv")) if (tables / "bench").exists() else []
    t1s = [p for p in t1s if not p.name.endswith("_per_type.csv")]
    if t1s:
        frames = [pd.read_csv(p) for p in t1s]
        allt1 = pd.concat(frames, ignore_index=True)
        wide = allt1.pivot_table(index="held_out", columns="model", values=["macro_f1", "rare_f1"])
        wide.columns = [f"{a}:{b}" for a, b in wide.columns]
        parts.append("## T1 label transfer (leave one atlas out)\n\n" + md_table(wide.reset_index()))
    fmc = tables / "bench" / "fm_compare.csv"
    if fmc.exists():
        parts.append(f"## Embeddings compared\n\nSee `{fmc}`.\n")
    parts.append("## How to read this\n\n"
                 "* A T3 score of 0 is the baseline the bundle was scored against (b1 when another dataset carried the same "
                 "condition family, else the control copied); 1 is the cross-animal split-half replicate. Values above 1 and "
                 "below 0 are reported, not clipped.\n"
                 "* `replicate_level = cells` means a group had a single animal, so r is a technical, not biological, ceiling.\n"
                 "* `gate_nmae = False` means fewer than the required real DE genes after exclusions; the fold-change member "
                 "is absent for that bundle.\n"
                 "* The stub model is plumbing. Its score means the pipeline runs, nothing more.\n")
    out.write_text("\n".join(parts), encoding="utf-8")
    print(f"  -> {out}")
    return out
