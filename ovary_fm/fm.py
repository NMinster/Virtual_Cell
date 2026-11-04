"""The zero-shot arm: one generalist's embedding of the corpus, then T1 on it.

``ovary-fm fm --tool <model>`` runs ``scripts/run_fm_embed.py`` in that
model's own conda env (the Mouse_Compensatory pattern: the main env never
imports torch) and gets back ``results/fm/<model>/embedding.parquet`` plus
a provenance json. ``ovary-fm fm-compare`` then runs T1 inside every
embedding found under results/fm/ and results/train/ and collects the
macro-F1 and rare-type F1 per held-out atlas beside the gene-space
logistic regression.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

from .config import REPO_ROOT, resolve_path, results_dir, tables_dir, write_provenance
from .io import corpus_path

STUB = "stub"


def embed(cfg, tool: str | None, adata_path: str | None = None):
    write_provenance(cfg, "fm", extra={"tool": tool})
    f = cfg["fm"]
    models = f.get("models") or {}
    if tool is None or tool not in models:
        raise ValueError(f"--tool is required: one of {list(models)}")
    spec = models[tool] or {}
    source = Path(adata_path) if adata_path else corpus_path(cfg)
    if not source.exists():
        raise FileNotFoundError(f"{source} is missing; run `ovary-fm harmonise` first")
    out = results_dir(cfg, "fm", tool)
    script = REPO_ROOT / "scripts" / "run_fm_embed.py"
    cmd = ["python", str(script), "--model", tool, "--adata", str(source), "--out", str(out),
           "--batch-size", str(f.get("batch_size", 8)), "--device", str(f.get("device", "cuda")),
           "--species", str(f.get("species", "mouse")), "--n-comps", str(f.get("n_comps", 50))]
    if spec.get("checkpoint"):
        cmd += ["--checkpoint", str(resolve_path(spec["checkpoint"]))]
    if spec.get("uce_dir"):
        cmd += ["--uce-dir", str(resolve_path(spec["uce_dir"]))]
    env_name = spec.get("env")
    if env_name and shutil.which("conda"):
        cmd = ["conda", "run", "--no-capture-output", "-n", env_name] + cmd
    elif env_name:
        print(f"  conda is not on PATH; running {tool} in the current env")
    if tool == STUB or not env_name:
        cmd[cmd.index("python")] = sys.executable
    print(f"  $ {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=str(REPO_ROOT))
    if result.returncode != 0:
        raise RuntimeError(f"{tool} exited with status {result.returncode}")
    prov = out / "fm_provenance.json"
    if prov.exists():
        record = json.loads(prov.read_text(encoding="utf-8"))
        print(f"  {tool}: {record.get('n_cells')} cells, {record.get('dims')} dims, {record.get('wall_seconds')} s")
    print(f"  now `ovary-fm t1 --rep {out / 'embedding.parquet'}` or `ovary-fm fm-compare`")


def compare(cfg) -> pd.DataFrame:
    from .bench import t1

    write_provenance(cfg, "fm-compare")
    found = {}
    for root in (results_dir(cfg, "fm"), results_dir(cfg, "train")):
        for p in sorted(root.glob("*/embedding.parquet")):
            found[p.parent.name] = p
    if not found:
        raise FileNotFoundError("no results/fm/<model>/embedding.parquet or results/train/<arm>/embedding.parquet")
    frames = []
    base = tables_dir(cfg, "bench") / "t1_logreg.csv"
    if base.exists():
        frames.append(pd.read_csv(base))
    else:
        frames.append(t1.run(cfg, rep=None, name="logreg"))
    for model, path in found.items():
        print(f"== {model}")
        try:
            frames.append(t1.run(cfg, rep=str(path), name=f"rep_{model}"))
        except Exception as exc:
            print(f"  !! {model}: {exc}")
    table = pd.concat(frames, ignore_index=True)
    wide = table.pivot_table(index="held_out", columns="model", values=["macro_f1", "rare_f1"])
    out = tables_dir(cfg, "bench") / "fm_compare.csv"
    wide.to_csv(out)
    print(wide.round(3).to_string())
    print(f"  -> {out}")
    return wide
