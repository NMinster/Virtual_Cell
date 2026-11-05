"""Unit + end-to-end tests.

    pytest -q

The end-to-end test builds the synthetic corpus from
``scripts/make_synthetic.py`` and runs ``ovary-fm all`` on it (qc,
harmonise, bench, the three baseline predictions scored, t1, report),
then the fm stub and t1 inside it. It asserts on the CONTENT of the
outputs: the modality is inferred right, the crosswalk maps what it
should and flags what it should not, the planted responses come out
differential, b0 scores 0 by construction, the transfer baseline beats
the control copy where a source exists, the replicate ceiling is
cross-animal where the design allows, and label transfer recovers the
types.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from ovary_fm.bench.metrics import bh, de_members, pds_cosine, rescale  # noqa: E402
from ovary_fm.config import apply_override, load_config  # noqa: E402
from ovary_fm.harmonise import crosswalk_lookup  # noqa: E402


# --- config and crosswalk -------------------------------------------------------------

def test_config_carries_the_plan():
    cfg = load_config()
    assert cfg["qc"]["max_mito_pct"] == {"cell": 20, "nucleus": 5}
    assert cfg["qc"]["mad"]["n_mads"] == 3
    assert cfg["harmonise"]["gene_space"] == "union"
    assert cfg["benchmark"]["t3"]["n_splits"] == 5 and cfg["benchmark"]["t3"]["base_seed"] == 0
    assert cfg["benchmark"]["t3"]["min_cells"] == 400
    assert cfg["benchmark"]["t3"]["de"]["min_cpm_control"] == 5.0
    assert cfg["benchmark"]["t3"]["de"]["nmae_min_genes"] == 10
    assert cfg["benchmark"]["t3"]["members"] == ["pds_cosine", "expr_mse", "de_direction_fidelity",
                                                 "de_direction_reach", "de_sig_jaccard", "de_lfc_nmae"]
    assert cfg["benchmark"]["t3"]["baseline"] == "b0"                        # the scale is anchored on the control copy
    assert cfg["benchmark"]["t3"]["min_de_genes_to_score"] == 5             # the signal gate
    assert cfg["benchmark"]["t3"]["pds"]["min_pool"] == 5
    assert cfg["benchmark"]["t3"]["replicate"]["expr_mse_depth_correction"] == 4
    assert "Egr1" in cfg["benchmark"]["t3"]["exclude_genes"]["ovulation"]
    assert cfg["benchmark"]["margins"]["t3_context_score"] == 0.05
    assert cfg["crosswalk"]["major_types"] == ["Oocytes", "Granulosa", "Luteal", "Theca", "Stromal", "Endothelial",
                                               "Epithelial", "Immune"]
    names = [d["name"] for d in cfg["corpus"]["datasets"]]
    assert {"isola2024", "morris2022", "huang2025", "mouse_comp_30d"} <= set(names)
    own = next(d for d in cfg["corpus"]["datasets"] if d["name"] == "mouse_comp_30d")
    assert own["held_out"] is True and own["enabled"] is False


def test_override_refuses_unknown_keys():
    cfg = load_config()
    apply_override(cfg, "qc.min_counts=300")
    assert cfg["qc"]["min_counts"] == 300
    with pytest.raises(KeyError):
        apply_override(cfg, "qc.no_such_key=1")


def test_crosswalk_lookup_levels():
    cw = load_config()["crosswalk"]
    hit = crosswalk_lookup(cw, "isola2024", "Epithelium B")
    assert hit["major"] == "Endothelial" and hit["cl_id"] == "CL:0002138"
    hit = crosswalk_lookup(cw, "morris2022", "GC_Mural")
    assert hit["major"] == "Granulosa" and hit["cl_id"] == "CL:4033067"
    hit = crosswalk_lookup(cw, "morris2022", "M_immature theca")      # a `verify` id: kept as None, flagged
    assert hit["major"] == "Theca" and hit["cl_id"] is None and hit["status"].endswith("verify_cl_id")
    hit = crosswalk_lookup(cw, "nobody", "granulosa cell")              # the global map
    assert hit["major"] == "Granulosa" and hit["status"] == "mapped_global"
    hit = crosswalk_lookup(cw, "nobody", "mystery cells")
    assert hit["major"] == "Unassigned" and hit["status"] == "unmapped"
    # the stem rule: Huang 2025's numbered clusters and time-split names
    for label, major in (("Granulosa 3", "Granulosa"), ("Stroma 1", "Stromal"), ("Theca_00", "Theca"),
                         ("Luteal Cells", "Luteal"), ("Endothelial 2", "Endothelial"), ("Myeloid", "Immune")):
        hit = crosswalk_lookup(cw, "huang2025", label)
        assert hit["major"] == major, (label, hit)
        assert hit["status"] in ("mapped_dataset", "mapped_global", "mapped_stem")


# --- the members -----------------------------------------------------------------------

def _de(genes, sig, lfc):
    return pd.DataFrame({"sig": sig, "lfc": lfc, "tested": True}, index=pd.Index(genes, name="gene"))


def test_de_members_on_toy_tables():
    genes = [f"g{i}" for i in range(12)]
    real = _de(genes, [True] * 10 + [False] * 2, [1.0] * 5 + [-1.0] * 5 + [0.0, 0.0])
    perfect = de_members(real, real, np.ones(12, dtype=bool))
    assert perfect["de_direction_fidelity"] == 1.0 and perfect["de_sig_jaccard"] == 1.0
    assert perfect["de_direction_reach"] == 1.0 and perfect["de_lfc_nmae"] == 0.0
    nothing = _de(genes, [False] * 12, [0.0] * 12)
    none = de_members(nothing, real, np.ones(12, dtype=bool))
    assert none["de_direction_fidelity"] == 0.0 and none["de_sig_jaccard"] == 0.0
    assert none["de_direction_reach"] == 0.0 and none["de_lfc_nmae"] == 1.0
    flipped = _de(genes, [True] * 10 + [False] * 2, [-1.0] * 5 + [1.0] * 5 + [0.0, 0.0])
    bad = de_members(flipped, real, np.ones(12, dtype=bool))
    assert bad["de_direction_fidelity"] == 0.0 and bad["de_sig_jaccard"] == 1.0 and bad["de_lfc_nmae"] == 2.0
    # the exclusion mask leaves the excluded genes out of every member
    keep = np.array([True] * 6 + [False] * 6)
    half = de_members(real, real, keep)
    assert half["n_real_de"] == 6 and half["de_lfc_nmae"] is not None


def test_pds_and_rescale():
    keep = np.ones(3, dtype=bool)
    pool = {"a": np.array([1.0, 0.0, 0.0]), "b": np.array([0.0, 1.0, 0.0]), "c": np.array([0.0, 0.0, 1.0])}
    assert pds_cosine(np.array([0.9, 0.1, 0.0]), "a", pool, keep) == 1.0
    assert pds_cosine(np.array([0.0, 0.0, 0.0]), "a", pool, keep) == 0.5
    assert pds_cosine(np.array([0.0, 0.0, 1.0]), "a", pool, keep) == 0.25     # c closer, b tied: mid-rank 1.5 of 2
    assert pds_cosine(np.array([0.0, 0.2, 1.0]), "a", pool, keep) == 0.0      # both others closer
    assert rescale("x", 0.5, 0.5, 0.9) == 0.0 and rescale("x", 0.9, 0.5, 0.9) == 1.0
    assert rescale("expr_mse", 0.02, 0.10, 0.01) == pytest.approx(0.888, abs=1e-2)
    assert rescale("expr_mse", 0.5, 0.10, 0.01, {"expr_mse": [0, 1]}) == 0.0
    assert rescale("de_lfc_nmae", 10.0, 1.0, 0.4, {"de_lfc_nmae_floor": -6}) == -6.0
    assert np.isnan(rescale("x", 0.5, 0.5, 0.5))
    p = bh(np.array([0.01, 0.04, 0.03, 0.5]))
    assert p.max() <= 1.0 and p[0] <= p[3]


# --- end to end ------------------------------------------------------------------------

def _run(args, cwd):
    cmd = [sys.executable, "-m", "ovary_fm.cli"] + args
    res = subprocess.run(cmd, cwd=str(cwd), capture_output=True, text=True)
    if res.returncode != 0:
        raise AssertionError(f"{' '.join(cmd)}\n--- stdout ---\n{res.stdout[-6000:]}\n--- stderr ---\n{res.stderr[-6000:]}")
    return res.stdout


@pytest.fixture(scope="module")
def synthetic(tmp_path_factory):
    out = tmp_path_factory.mktemp("fake")
    res = subprocess.run([sys.executable, str(REPO_ROOT / "scripts" / "make_synthetic.py"), "--out", str(out),
                          "--n-per-animal", "260", "--n-extra-genes", "500"], capture_output=True, text=True)
    assert res.returncode == 0, res.stderr[-4000:]
    sets = [f"--set=paths.corpus_manifest={out / 'corpus.yaml'}", f"--set=paths.crosswalk={out / 'crosswalk.yaml'}",
            f"--set=paths.benchmark={out / 'benchmark.yaml'}", "--data", str(out / "data"), "--results", str(out / "results")]
    log = _run(["all"] + sets, REPO_ROOT)
    return {"out": out, "sets": sets, "log": log, "results": out / "results", "data": out / "data"}


def test_qc_infers_the_modality_and_keeps_the_design(synthetic):
    qc = pd.read_csv(synthetic["results"] / "tables" / "qc" / "qc_summary.csv").set_index("dataset")
    assert qc.loc["synth_isola", "suspension_type"] == "cell"
    assert qc.loc["synth_benayoun", "suspension_type"] == "nucleus"
    assert qc.loc["synth_mouse_comp", "suspension_type"] == "nucleus"
    assert qc.loc["synth_isola", "n_animals"] == 8 and qc.loc["synth_huang", "n_animals"] == 9
    assert qc.loc["synth_benayoun", "n_animals"] == 6            # read from the GEO RAW layout, one matrix per animal
    assert qc.loc["synth_benayoun", "n_author_labels"] == 9      # the labels came from the cell_info table
    assert qc.loc["synth_unlabelled", "n_author_labels"] == 1    # 'unlabelled'
    assert qc.loc["synth_unlabelled", "n_animals"] == 6          # twelve csv libraries, animal_regex folds them to six
    assert set(qc.loc["synth_unlabelled", "conditions"].split(";")) == {"young", "old"}   # regex patterns on the sample names
    assert set(qc.loc["synth_huang", "conditions"].split(";")) == {"hcg_0h", "hcg_4h", "hcg_12h"}
    for name in qc.index:
        assert qc.loc[name, "cells_out"] >= 0.8 * qc.loc[name, "cells_in"]
    drop = pd.read_csv(synthetic["results"] / "tables" / "qc" / "synth_benayoun_dropoff.csv")
    assert drop["step"].str.contains("nucleus").any()


def test_harmonise_crosswalk_and_corpus(synthetic):
    import anndata as ad

    manifest = pd.read_csv(synthetic["results"] / "tables" / "harmonise" / "manifest.csv").set_index("dataset")
    assert manifest.loc["synth_isola", "crosswalk_coverage"] == 1.0
    assert manifest.loc["synth_huang", "crosswalk_coverage"] == 1.0
    assert manifest.loc["synth_benayoun", "n_unmapped_labels"] == 1          # "mystery cells"
    assert bool(manifest.loc["synth_mouse_comp", "held_out"]) is True
    review = pd.read_csv(synthetic["results"] / "tables" / "harmonise" / "crosswalk_review.csv")
    unmapped = review[review["crosswalk_status"] == "unmapped"]
    assert list(unmapped["author_label"]) == ["mystery cells"]
    assert (review.loc[review["author_label"] == "Epithelium B", "major_type"] == "Endothelial").all()
    corpus = ad.read_h5ad(synthetic["data"] / "corpus" / "corpus.h5ad")
    assert set(corpus.obs["dataset"]) == {"synth_isola", "synth_benayoun", "synth_unlabelled", "synth_huang"}
    # the unlabelled dataset got its types transferred, and they are right
    un = corpus[corpus.obs["dataset"] == "synth_unlabelled"]
    assert set(un.obs["crosswalk_status"].astype(str)) <= {"transferred", "transferred_low_confidence"}
    assert manifest.loc["synth_unlabelled", "labels"] == "transferred"
    cell_info = pd.read_csv(synthetic["data"] / "raw" / "synth_unlabelled" / "synth_unlabelled_cell_info.txt.gz", sep="\t").set_index("cell")
    bare = [i.split(":", 1)[1] for i in un.obs_names]
    assert (un.obs["major_type"].astype(str).to_numpy() == cell_info["true_type"].reindex(bare).to_numpy()).mean() > 0.95
    assert set(un.obs["assay"].astype(str)) == {"10x 3'", "Smart-seq2"}     # assay_patterns per library
    assert (synthetic["data"] / "corpus" / "heldout_synth_mouse_comp.h5ad").exists()
    assert "suspension_type" in corpus.obs and set(corpus.obs["suspension_type"]) == {"cell", "nucleus"}
    truth = corpus.obs["major_type"].astype(str)
    assert (truth == "Unassigned").mean() < 0.15
    # the crosswalk recovers the planted type wherever it maps
    raw = ad.read_h5ad(synthetic["data"] / "raw" / "synth_isola" / "synth_isola.h5ad")
    joined = corpus.obs.loc[corpus.obs["dataset"] == "synth_isola", "major_type"].astype(str)
    planted = raw.obs["true_type"].reindex([i.split(":", 1)[1] for i in joined.index]).to_numpy()
    assert (joined.to_numpy() == planted).mean() > 0.99


def test_bench_ceilings_and_baselines(synthetic):
    ceilings = pd.read_csv(synthetic["results"] / "tables" / "bench" / "ceilings.csv")
    scorable = ceilings[ceilings["scorable"]]
    assert len(scorable) >= 8
    # the held-out surgery context is built (its own replicate ceiling), the rare compartments are listed not scored
    assert (ceilings["context"] == "synth_mouse_comp_surgery").any()
    oocytes = ceilings[ceilings["compartment"] == "Oocytes"]
    assert len(oocytes) and not oocytes["scorable"].any()
    # cross-animal halves where the design allows (4 + 4, 3 + 3, 3 per time point, 2 + 2)
    assert (scorable["replicate_level"] == "animals").all()
    # the planted ageing response comes out differential in the stroma, and the b1 transfer exists across the two ageing sets
    iso = scorable[(scorable["context"] == "synth_isola_ageing") & (scorable["compartment"] == "Stromal")].iloc[0]
    assert iso["n_de_real"] >= 4 and iso["n_b1_sources"] >= 1          # the other ageing datasets (labelled or transferred)
    de = pd.read_csv(synthetic["results"] / "tables" / "bench" / "de" / "synth_isola_ageing_age_9mo_Stromal.csv").set_index("gene")
    assert de.loc["Cdkn1a", "sig"] and de.loc["Cdkn1a", "lfc"] > 1 and de.loc["Col1a1", "lfc"] < -0.5
    # ovulation has no other dataset to transfer from
    hu = scorable[scorable["context"] == "synth_huang_ovulation"]
    assert (hu["n_b1_sources"] == 0).all()
    # the excluded immediate-early genes are out of the kept DE count
    g4 = hu[(hu["condition"] == "hcg_4h") & (hu["compartment"] == "Granulosa")].iloc[0]
    assert g4["n_de_real_kept"] < g4["n_de_real"]
    # b0 is a worse fit than the replicate on the DE members wherever there are DE genes; on the
    # expression member only where the effect is large: the lite replicate is measured at half
    # depth and carries no jackknife correction, so a small effect can sit below the noise
    with_de = scorable[scorable["n_de_real_kept"] >= 10]
    assert (with_de["b0_de_direction_fidelity"] <= with_de["r_de_direction_fidelity"]).all()
    big = scorable[scorable["n_de_real_kept"] >= 50]
    assert len(big) >= 2 and (big["b0_expr_mse"] >= big["r_expr_mse"]).all()


def test_scores_are_anchored(synthetic):
    tables = synthetic["results"] / "tables" / "bench"
    b0 = pd.read_csv(tables / "scores_b0.csv")
    b1 = pd.read_csv(tables / "scores_b1.csv")
    stub = pd.read_csv(tables / "scores_stub.csv")
    # b0 scored against itself is exactly 0 on every deterministic member (the scale is anchored on it)
    for m in ("expr_mse", "de_direction_fidelity", "de_sig_jaccard"):
        assert np.allclose(b0[f"s_{m}"].dropna(), 0.0, atol=1e-9)
    assert np.allclose(b0["bundle_score"].dropna(), 0.0, atol=1e-9)
    # b1 is scored as a model: positive where a source carries the same response
    assert (b1["bundle_score"] > 0).mean() >= 0.75
    # where a transfer source exists the transfer beats the control copy on the DE members
    merged = b1.merge(b0, on="bundle_id", suffixes=("_b1", "_b0"))
    merged = merged[merged["n_real_de_b1"] >= 10]
    assert len(merged) >= 2
    assert (merged["u_de_direction_fidelity_b1"] > merged["u_de_direction_fidelity_b0"]).mean() >= 0.75
    assert (merged["u_de_sig_jaccard_b1"] > merged["u_de_sig_jaccard_b0"]).mean() >= 0.75
    # the transfer must not blow up the expression error (a dozen planted genes barely move the pseudobulk)
    assert (merged["u_expr_mse_b1"] <= merged["u_expr_mse_b0"] * 1.5).all()
    # the stub covers every scorable bundle, and the summary has one row per context plus OVERALL
    ceilings = pd.read_csv(tables / "ceilings.csv")
    assert set(stub["bundle_id"]) == set(ceilings.loc[ceilings["scored"] == True, "bundle_id"])
    assert (ceilings["scored"] == True).sum() < ceilings["scorable"].sum()      # the gate removes the null compartments
    summary = pd.read_csv(tables / "scores_stub_summary.csv")
    assert summary["context"].iloc[-1] == "OVERALL" and summary.shape[0] == 6


def test_t1_recovers_the_types(synthetic):
    t1 = pd.read_csv(synthetic["results"] / "tables" / "bench" / "t1_logreg.csv").set_index("held_out")
    assert set(t1.index) == {"synth_isola", "synth_benayoun", "synth_huang", "synth_mouse_comp"}   # never the unlabelled one
    assert (t1["macro_f1"] > 0.85).all()
    assert t1.loc["synth_mouse_comp", "role"] == "test_only"
    per_type = pd.read_csv(synthetic["results"] / "tables" / "bench" / "t1_logreg_per_type.csv")
    assert per_type.loc[per_type["major_type"] == "Oocytes", "f1"].mean() > 0.7


def test_fm_stub_and_report(synthetic):
    log = _run(["fm", "--tool", "stub"] + synthetic["sets"], REPO_ROOT)
    emb = synthetic["results"] / "fm" / "stub" / "embedding.parquet"
    assert emb.exists()
    _run(["t1", "--rep", str(emb)] + synthetic["sets"], REPO_ROOT)
    rep = pd.read_csv(synthetic["results"] / "tables" / "bench" / "t1_rep_stub.csv").set_index("held_out")
    assert (rep["macro_f1"] > 0.6).all()
    _run(["fm-compare"] + synthetic["sets"], REPO_ROOT)
    assert (synthetic["results"] / "tables" / "bench" / "fm_compare.csv").exists()
    _run(["train", "--arm", "a"] + synthetic["sets"], REPO_ROOT)
    assert (synthetic["results"] / "train" / "arm_a" / "train.sbatch").exists()
    _run(["report"] + synthetic["sets"], REPO_ROOT)
    text = (synthetic["results"] / "report" / "first_results.md").read_text()
    assert "T3 bundles and ceilings" in text and "T1 label transfer" in text and "mystery cells" in text
    prov = list((synthetic["results"] / "provenance").glob("bench_*.json"))
    assert prov
