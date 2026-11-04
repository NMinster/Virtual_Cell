"""``ovary-fm train --arm c``: scVI, then scANVI on the crosswalk labels.

The primary baseline of the preregistration. Writes
results/train/arm_c/embedding.parquet in the fm format, so
``ovary-fm t1 --rep results/train/arm_c/embedding.parquet`` scores it on
the same splits as the gene-space baseline and the generalists. Needs
scvi-tools (environment-scvi.yml); the call is written from the
scvi-tools 1.x API and was not run on the machine this was written on.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import pandas as pd

from ..config import results_dir, write_provenance
from ..io import corpus_path


def train(cfg, adata_path: str | None = None) -> Path:
    write_provenance(cfg, "train", extra={"arm": "c"})
    try:
        import scvi
    except ImportError:
        raise SystemExit("scvi-tools is not installed in this env: `bash biowulf/setup_env.sh --train scvi` "
                         "and run the stage through `conda run -n ovary_fm_scvi`, or pip install -e .[scvi]") from None
    import anndata as ad
    import scanpy as sc

    spec = cfg["models"]["arm_c"]
    seed = int(cfg.get("random_state", 0))
    scvi.settings.seed = seed
    start = time.time()
    adata = ad.read_h5ad(adata_path or corpus_path(cfg))
    adata.layers["counts"] = adata.X.copy()
    sc.pp.highly_variable_genes(adata, n_top_genes=int(spec.get("n_top_genes", 2000)), flavor="seurat_v3",
                                layer="counts", batch_key=spec.get("batch_key", "dataset"), subset=True)
    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key=spec.get("batch_key", "dataset"))
    model = scvi.model.SCVI(adata, n_latent=int(spec.get("n_latent", 30)), n_layers=int(spec.get("n_layers", 2)),
                            n_hidden=int(spec.get("n_hidden", 128)))
    model.train(max_epochs=int(spec.get("max_epochs", 50)))
    out = results_dir(cfg, "train", "arm_c")
    model.save(str(out / "scvi"), overwrite=True)
    latent = model.get_latent_representation()
    labels_key = spec.get("labels_key", "major_type")
    if labels_key in adata.obs:
        lvae = scvi.model.SCANVI.from_scvi_model(model, labels_key=labels_key,
                                                 unlabeled_category=spec.get("unlabeled_category", "Unassigned"))
        lvae.train(max_epochs=int(spec.get("scanvi_epochs", 20)), n_samples_per_label=100)
        lvae.save(str(out / "scanvi"), overwrite=True)
        latent = lvae.get_latent_representation()
        adata.obs["scanvi_prediction"] = lvae.predict()
        adata.obs[["dataset", labels_key, "scanvi_prediction"]].to_csv(out / "scanvi_predictions.csv")
    frame = pd.DataFrame(latent, index=adata.obs_names.astype(str))
    frame.columns = [f"dim_{i}" for i in range(frame.shape[1])]
    frame.to_parquet(out / "embedding.parquet")
    record = {"model": "scvi_scanvi", "n_cells": int(frame.shape[0]), "dims": int(frame.shape[1]),
              "n_genes": int(adata.n_vars), "wall_seconds": round(time.time() - start, 1), "spec": spec,
              "scvi_version": getattr(scvi, "__version__", "?")}
    (out / "fm_provenance.json").write_text(json.dumps(record, indent=2, default=str), encoding="utf-8")
    print(f"  arm C: {frame.shape[0]:,} cells x {frame.shape[1]} dims -> {out / 'embedding.parquet'}")
    return out / "embedding.parquet"
