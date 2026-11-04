#!/usr/bin/env python
"""One foundation model's cell embedding, in that model's own conda env.

    python scripts/run_fm_embed.py --model stub \
        --adata data/corpus/corpus.h5ad --out results/fm/stub
    conda run -n ovary_fm_scprint2 python scripts/run_fm_embed.py \
        --model scprint2 --adata data/corpus/corpus.h5ad --out results/fm/scprint2 \
        --checkpoint small-v2
    conda run -n ovary_fm_uce python scripts/run_fm_embed.py \
        --model uce --adata ... --out results/fm/uce \
        --uce-dir data/fm/UCE --checkpoint data/fm/uce/33l_8ep_1024t_1280.torch
    conda run -n ovary_fm_transcriptformer python scripts/run_fm_embed.py \
        --model transcriptformer --adata ... --out results/fm/transcriptformer \
        --checkpoint data/fm/transcriptformer/tf_exemplar

`ovary-fm fm --tool <model>` builds these commands. Whatever the model,
the script writes the same two files, which `ovary-fm fm-compare` reads:

    <out>/embedding.parquet     obs_names x dim_0..dim_{d-1}
    <out>/fm_provenance.json    model, checkpoint (+ sha256), genes matched,
                                batch size, device, wall time, versions

The object is the corpus (or any h5ad given to --adata): integer counts in X (or layers['counts']),
var_names = mouse gene symbols, var['gene_ids'] = Ensembl ids (kept by the
10x reader). Each model gets the identifier it wants:

* stub -- scanpy normalise / log1p / 2000 HVGs / PCA. No model; this is
  what the test suite runs to exercise the plumbing and the kill rule.
* scprint2 -- `scprint2 embed` (cantinilab/scPRINT-2) with
  obs['organism_ontology_term_id'] = NCBITaxon:10090 and Ensembl ids as
  var_names; the checkpoint name resolves through the Hugging Face hub
  (jkobject/scPRINT) or a local .ckpt path.
* uce -- shells out to `eval_single_anndata.py` in a clone of
  snap-stanford/UCE (`--uce-dir`), symbols as var_names, `--species mouse`;
  reads obsm['X_uce'] from the h5ad it writes.
* transcriptformer -- `transcriptformer inference` with var['ensembl_id'],
  `--gene-col-name ensembl_id`, `--emb-type cell`; reads obsm['embeddings'].

UNTESTED against the real models: the three model branches are written
from the tools' READMEs (September 2026) and the stub is the only branch
the tests run. Read the log the first time; the flag names above are the
ones to check when a model's CLI has moved.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

MOUSE_TAXON = "NCBITaxon:10090"
HUMAN_TAXON = "NCBITaxon:9606"


def sha256_of(path) -> str | None:
    path = Path(path)
    if not path.is_file():
        return None
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def versions(*names) -> dict:
    import importlib

    out = {}
    for name in names:
        try:
            out[name] = getattr(importlib.import_module(name), "__version__", "?")
        except Exception:
            out[name] = None
    return out


def counts_first(adata):
    """Counts into X; symbols in var_names; Ensembl ids in var['gene_ids'] when present."""
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
    for layer in [k for k in list(adata.layers.keys()) if k is not None]:   # anndata >= 0.13 lists None (X) as a layer
        del adata.layers[layer]
    adata.obsm.clear()
    adata.obsp.clear()
    adata.var_names_make_unique()
    return adata


def write_outputs(out: Path, frame, record: dict):
    out.mkdir(parents=True, exist_ok=True)
    frame.columns = [f"dim_{i}" for i in range(frame.shape[1])]
    frame.index = frame.index.astype(str)
    frame.to_parquet(out / "embedding.parquet")
    (out / "fm_provenance.json").write_text(json.dumps(record, indent=2, default=str),
                                            encoding="utf-8")
    print(f"[fm] {frame.shape[0]:,} cells x {frame.shape[1]} dims -> {out / 'embedding.parquet'}")


# --- models ---------------------------------------------------------------------------------

def embed_stub(adata, args):
    import pandas as pd
    import scanpy as sc

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    n_top = min(2000, adata.n_vars - 1)
    try:
        sc.pp.highly_variable_genes(adata, n_top_genes=n_top, flavor="seurat")
        sub = adata[:, adata.var["highly_variable"]].copy()
    except Exception:
        sub = adata.copy()
    sc.pp.scale(sub, max_value=10)
    n = max(2, min(int(args.n_comps), sub.n_obs - 1, sub.n_vars - 1))
    sc.tl.pca(sub, n_comps=n, random_state=0)
    frame = pd.DataFrame(sub.obsm["X_pca"], index=adata.obs_names)
    return frame, {"n_genes_matched": int(sub.n_vars)}


def embed_scprint2(adata, args, out: Path):
    import anndata as ad
    import pandas as pd

    if shutil.which("scprint2") is None:
        try:
            import scprint2  # noqa: F401
        except ImportError as exc:
            raise SystemExit("[fm] scprint2 is not installed in this env: build it with "
                             "`bash biowulf/setup_env.sh --fm scprint2` "
                             "(environment-fm-scprint2.yml)") from exc
    taxon = MOUSE_TAXON if args.species == "mouse" else HUMAN_TAXON
    adata.obs["organism_ontology_term_id"] = taxon
    print("[fm] scprint2 reads its ontologies from a laminDB instance: once, `lamin init --storage "
          "/data/$USER/ovary_fm/lamin --name ovary_fm --modules bionty && scprint2 easy_setup`")
    if "gene_ids" in adata.var:
        adata.var["symbol"] = adata.var_names
        adata.var_names = adata.var["gene_ids"].astype(str).to_numpy()
        adata.var_names_make_unique()
    work = out / "scprint2_work"
    work.mkdir(parents=True, exist_ok=True)
    src = work / "input.h5ad"
    adata.write_h5ad(src)
    checkpoint = args.checkpoint or "small-v2"
    ckpt = checkpoint if checkpoint.endswith(".ckpt") else f"{checkpoint}.ckpt"
    if not Path(ckpt).exists():
        # the checkpoint name resolves through the Hugging Face hub (jkobject/scPRINT)
        pull = ["hf", "download", "jkobject/scPRINT", ckpt, "--local-dir", str(work)]
        if shutil.which("hf") is None:
            pull = ["huggingface-cli", "download", "jkobject/scPRINT", ckpt, "--local-dir", str(work)]
        print(f"[fm] $ {' '.join(pull)}")
        subprocess.run(pull, check=True)
        ckpt = str(work / ckpt)
    output = work / "embedded.h5ad"
    cmd = ["scprint2", "embed", "--adata", str(src), "--ckpt_path", ckpt, "--species", taxon,
           "--output_filename", str(output), "--batch_size", str(args.batch_size)]
    if args.device == "cpu":
        cmd += ["--precision", "32"]
    print(f"[fm] $ {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    emb = ad.read_h5ad(output)
    key = next((k for k in emb.obsm if "scprint" in k.lower() or "emb" in k.lower()), None)
    if key is None:
        raise SystemExit(f"[fm] no embedding in {output}: obsm keys {list(emb.obsm)}")
    frame = pd.DataFrame(emb.obsm[key], index=emb.obs_names)
    labels = [c for c in emb.obs.columns if c.startswith("conv_pred_") or c.endswith("_pred")]
    if labels:
        emb.obs[labels].to_csv(out / "predicted_labels.csv")
        print(f"[fm] predicted labels ({labels}) -> {out / 'predicted_labels.csv'}")
    return frame, {"n_genes_matched": None, "obsm_key": key, "checkpoint": ckpt}


def embed_uce(adata, args, out: Path):
    import anndata as ad
    import pandas as pd

    if not args.uce_dir:
        raise SystemExit("[fm] --uce-dir (a clone of github.com/snap-stanford/UCE) is required")
    uce_dir = Path(args.uce_dir)
    script = uce_dir / "eval_single_anndata.py"
    if not script.exists():
        raise SystemExit(f"[fm] {script} not found; clone snap-stanford/UCE there")
    work = out / "uce_work"
    work.mkdir(parents=True, exist_ok=True)
    src = work / "input.h5ad"
    adata.write_h5ad(src)      # symbols as var_names, counts in X, as UCE wants
    cmd = [sys.executable, str(script), "--adata_path", str(src), "--dir", str(work) + "/",
           "--species", args.species, "--batch_size", str(args.batch_size)]
    if args.checkpoint:
        cmd += ["--model_loc", str(args.checkpoint)]
        if "33l" in Path(args.checkpoint).name:
            cmd += ["--nlayers", "33"]
    if args.device == "cpu":
        cmd += ["--cpu"]
    print(f"[fm] $ {' '.join(cmd)}")
    subprocess.run(cmd, check=True, cwd=str(uce_dir))
    output = work / "input.h5ad"
    written = sorted(p for p in work.glob("*.h5ad") if p != src) or [output]
    emb = ad.read_h5ad(written[-1])
    if "X_uce" not in emb.obsm:
        raise SystemExit(f"[fm] obsm['X_uce'] missing in {written[-1]}: {list(emb.obsm)}")
    frame = pd.DataFrame(emb.obsm["X_uce"], index=emb.obs_names)
    return frame, {"n_genes_matched": None, "obsm_key": "X_uce"}


def embed_transcriptformer(adata, args, out: Path):
    import anndata as ad
    import pandas as pd

    if shutil.which("transcriptformer") is None:
        raise SystemExit("[fm] the transcriptformer CLI is not on PATH: build the env with "
                         "`bash biowulf/setup_env.sh --fm transcriptformer`")
    if "gene_ids" not in adata.var:
        raise SystemExit("[fm] var['gene_ids'] (Ensembl ids) is needed for TranscriptFormer")
    adata.var["ensembl_id"] = adata.var["gene_ids"].astype(str).to_numpy()
    work = out / "transcriptformer_work"
    work.mkdir(parents=True, exist_ok=True)
    src = work / "input.h5ad"
    adata.write_h5ad(src)
    checkpoint = args.checkpoint or "checkpoints/tf_exemplar"
    cmd = ["transcriptformer", "inference", "--checkpoint-path", str(checkpoint),
           "--data-file", str(src), "--output-path", str(work), "--output-filename",
           "embeddings.h5ad", "--batch-size", str(args.batch_size), "--gene-col-name",
           "ensembl_id", "--emb-type", "cell"]
    emb_file = Path(checkpoint).parent / "all_embeddings" / f"{'mus_musculus' if args.species == 'mouse' else 'homo_sapiens'}_gene.h5"
    if emb_file.exists():
        cmd += ["--pretrained-embedding", str(emb_file)]
    print(f"[fm] $ {' '.join(cmd)}")
    subprocess.run(cmd, check=True)
    emb = ad.read_h5ad(work / "embeddings.h5ad")
    if "embeddings" not in emb.obsm:
        raise SystemExit(f"[fm] obsm['embeddings'] missing: {list(emb.obsm)}")
    frame = pd.DataFrame(emb.obsm["embeddings"], index=emb.obs_names)
    return frame, {"n_genes_matched": None, "obsm_key": "embeddings"}


# --- main -----------------------------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--model", required=True, choices=["stub", "scprint2", "uce", "transcriptformer"])
    parser.add_argument("--adata", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--checkpoint", default=None, help="checkpoint path or hub name")
    parser.add_argument("--uce-dir", default=None, help="uce: a clone of snap-stanford/UCE")
    parser.add_argument("--batch-size", type=int, default=8)
    parser.add_argument("--device", default="cuda", choices=["cuda", "cpu"])
    parser.add_argument("--species", default="mouse")
    parser.add_argument("--n-comps", type=int, default=50, help="stub: PCA components")
    args = parser.parse_args(argv)

    import anndata as ad

    start = time.time()
    out = Path(args.out)
    adata = counts_first(ad.read_h5ad(args.adata))
    print(f"[fm] {args.model}: {adata.n_obs:,} cells x {adata.n_vars:,} genes from {args.adata}")
    extra = {}
    if args.model == "stub":
        frame, extra = embed_stub(adata, args)
    elif args.model == "scprint2":
        frame, extra = embed_scprint2(adata, args, out)
    elif args.model == "uce":
        frame, extra = embed_uce(adata, args, out)
    else:
        frame, extra = embed_transcriptformer(adata, args, out)
    record = {"model": args.model, "checkpoint": extra.get("checkpoint", args.checkpoint),
              "checkpoint_sha256": sha256_of(extra.get("checkpoint", args.checkpoint) or ""),
              "n_cells": int(frame.shape[0]), "n_genes_in": int(adata.n_vars),
              "n_genes_matched": extra.get("n_genes_matched"), "dims": int(frame.shape[1]),
              "batch_size": args.batch_size, "device": args.device, "species": args.species,
              "wall_seconds": round(time.time() - start, 1),
              "package_versions": versions("scanpy", "anndata", "torch", "scprint2", "transcriptformer"),
              "command": sys.argv, **{k: v for k, v in extra.items() if k not in ("checkpoint", "n_genes_matched")}}
    write_outputs(out, frame, record)


if __name__ == "__main__":
    main()
