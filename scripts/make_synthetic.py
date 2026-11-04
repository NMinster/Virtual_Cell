#!/usr/bin/env python
"""A synthetic four-dataset ovary corpus shaped like the real manifest.

    python scripts/make_synthetic.py --out /tmp/fake
    ovary-fm all --set paths.corpus_manifest=/tmp/fake/corpus.yaml \\
        --set paths.crosswalk=/tmp/fake/crosswalk.yaml \\
        --set paths.benchmark=/tmp/fake/benchmark.yaml \\
        --data /tmp/fake/data --results /tmp/fake/results

Four datasets, each with its own label vocabulary, modality and design:

* synth_isola      cells, 4 + 4 animals, 3 vs 9 months (ageing), Isola-like labels
* synth_benayoun   nuclei, 3 + 3 animals, young vs old (ageing), free-text labels
                   with one label the crosswalk does not know
* synth_huang      cells, 3 animals x hCG 0 / 4 / 12 h (ovulation), Morris-like labels;
                   the immediate-early genes the benchmark excludes are planted too
* synth_mouse_comp nuclei, 2 + 2 animals, sham vs ooph (surgery), the eight types;
                   held out, like the real one

Each type's panel genes are genuinely enriched and the planted responses
are real in the pseudobulk, so the output LOOKS right: a clean run here
says the code works end to end and says nothing about the real data.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent

TYPES = ["Oocytes", "Granulosa", "Luteal", "Theca", "Stromal", "Endothelial", "Epithelial", "Immune"]
SHARE = {"Oocytes": 0.03, "Granulosa": 0.27, "Luteal": 0.08, "Theca": 0.14, "Stromal": 0.18,
         "Endothelial": 0.11, "Epithelial": 0.06, "Immune": 0.13}
PANELS = {
    "Oocytes": ["Zp3", "Zp2", "Gdf9", "Bmp15", "Ooep", "Nlrp5", "Padi6", "Dppa3", "Zar1", "Nobox"],
    "Granulosa": ["Fshr", "Esr2", "Nppc", "Cyp19a1", "Hcn1", "Fst", "Hsd17b1", "Inhbb", "Bmpr1b", "Amh"],
    "Luteal": ["Ptgfr", "Sfrp4", "Akr1c18", "Ccn5", "Gm2a", "Ephx2", "Runx2", "Parm1", "Trpm6", "Nrcam"],
    "Theca": ["Cyp17a1", "Nckap5", "Srd5a1", "Mgarp", "Hhip", "Hao2", "Serpina5", "Aldh1a1", "Enpep", "Ptch1"],
    "Stromal": ["Dcn", "Lum", "Col1a1", "Col1a2", "Col3a1", "Pdgfra", "Tcf21", "Lama2", "Mfap4", "Rgs5"],
    "Endothelial": ["Pecam1", "Cdh5", "Kdr", "Flt1", "Emcn", "Cldn5", "Ptprb", "Egfl7", "Cd34", "Mecom"],
    "Epithelial": ["Epcam", "Lgr5", "Gpm6a", "Upk1b", "Upk3b", "Lgals7", "Msln", "Krt7", "Krt8", "Krt19"],
    "Immune": ["Ptprc", "Dock2", "Arhgap15", "Skap1", "Cd68", "C1qa", "C1qb", "Adgre1", "Cd3e", "Cd79a"],
}
MITO = ["mt-Co1", "mt-Co2", "mt-Co3", "mt-Nd1", "mt-Nd2", "mt-Atp6", "mt-Cytb"]
NUCLEAR = ["Malat1", "Neat1", "Meg3"]
RIBO = ["Rpl13a", "Rpl3", "Rps6", "Rps18", "Rplp0"]
MODALITY = {"cell": {"mito": 0.05, "nuclear": 0.005, "ribo": 0.08, "depth": 4000},
            "nucleus": {"mito": 0.005, "nuclear": 0.04, "ribo": 0.005, "depth": 2500}}

# planted responses (fold change on the rate), per condition family
AGEING = {"Stromal": {"Col1a1": 0.5, "Cdkn1a": 3.0, "Cdkn2a": 3.0, "C3": 2.5, "Lgals3": 2.0, "Xist": 2.0},
          "Immune": {"Lgals3": 3.0, "Spp1": 3.0, "Cd74": 2.0, "Xist": 2.0},
          "Granulosa": {"Amh": 0.5, "Xist": 2.0}}
OVULATION = {
    4: {"Granulosa": {"Egr1": 8, "Fos": 6, "Nr4a1": 5, "Ptgs2": 8, "Areg": 6, "Tnfaip6": 4, "Has2": 3, "Runx2": 3,
                      "Sult1e1": 3, "Adamts1": 2.5, "Cyp19a1": 0.4, "Lhcgr": 0.5, "Fshr": 0.5},
        "Theca": {"Star": 2.5, "Nr4a1": 3}},
    12: {"Granulosa": {"Egr1": 1.5, "Ptgs2": 3, "Tnfaip6": 6, "Has2": 4, "Adamts1": 5, "Runx2": 4, "Star": 4,
                       "Cyp11a1": 3, "Cyp19a1": 0.2, "Lhcgr": 0.3, "Pgr": 3},
         "Theca": {"Star": 3, "Cyp11a1": 2}},
}
SURGERY = {"Granulosa": {"Cyp19a1": 3.0, "Inhba": 2.5, "Star": 2.0, "Fshr": 0.4, "Amh": 0.35},
           "*": {"Cdkn1a": 2.5}}

VOCAB = {
    "isola": {"Granulosa": ["Granulosa A", "Granulosa B"], "Luteal": ["Luteal"], "Theca": ["Theca"],
              "Stromal": ["Stroma A", "Stroma B", "Stroma C"], "Endothelial": ["Endothelium"],
              "Epithelial": ["Epithelium A"], "Immune": ["Phagocytes A", "T-Lymphocytes"], "Oocytes": ["Oocytes"]},
    "morris": {"Granulosa": ["GC_Preantral", "GC_Antral", "GC_Mural"], "Luteal": ["GC_CL_Active"],
               "Theca": ["M_Steroidogenic theca"], "Stromal": ["M_Cortical stroma", "M_Pericyte"],
               "Endothelial": ["EN_Blood"], "Epithelial": ["Epithelium"], "Immune": ["Macrophage", "T_cell"],
               "Oocytes": ["Oocyte"]},
    "free": {"Granulosa": ["granulosa cell"], "Luteal": ["luteal cell"], "Theca": ["theca cell"],
             "Stromal": ["fibroblast", "mystery cells"], "Endothelial": ["endothelial cell"],
             "Epithelial": ["ovarian surface epithelial cell"], "Immune": ["macrophage"], "Oocytes": ["oocyte"]},
    "own": {t: [t] for t in TYPES},
}

DATASETS = [
    {"name": "synth_isola", "modality": "cell", "vocab": "isola", "label_col": "cluster.names",
     "cond_col": "development_stage", "family": "ageing",
     "conditions": [("3-month-old stage", "age_3mo", 4, None), ("9-month-old stage", "age_9mo", 4, AGEING)],
     "map": {"3-month-old stage": "age_3mo", "9-month-old stage": "age_9mo"}, "control": "age_3mo",
     "age_days": {"age_3mo": 90, "age_9mo": 270}},
    {"name": "synth_benayoun", "modality": "nucleus", "vocab": "free", "label_col": "cell_type",
     "cond_col": "age_group", "family": "ageing",
     "conditions": [("young", "young", 3, None), ("old", "old", 3, "AGEING_0.8")],
     "map": {"young": "young", "old": "old"}, "control": "young", "age_days": {"young": 120, "old": 540}},
    {"name": "synth_huang", "modality": "cell", "vocab": "morris", "label_col": "Level1",
     "cond_col": "timepoint", "family": "ovulation",
     "conditions": [("0h", "hcg_0h", 3, None), ("4h", "hcg_4h", 3, OVULATION[4]), ("12h", "hcg_12h", 3, OVULATION[12])],
     "map": {"0h": "hcg_0h", "4h": "hcg_4h", "12h": "hcg_12h"}, "control": "hcg_0h",
     "time_h": {"hcg_0h": 0, "hcg_4h": 4, "hcg_12h": 12}},
    {"name": "synth_mouse_comp", "modality": "nucleus", "vocab": "own", "label_col": "major_type",
     "cond_col": "condition", "family": "surgery",
     "conditions": [("sham", "sham", 2, None), ("ooph", "ooph_30d", 2, SURGERY)],
     "map": {"sham": "sham", "ooph": "ooph_30d"}, "control": "sham", "held_out": True},
]


def gene_universe(n_extra: int) -> list[str]:
    genes = [g for panel in PANELS.values() for g in panel]
    for fam in (AGEING, SURGERY, *OVULATION.values()):
        for changes in fam.values():
            genes += list(changes)
    genes += ["Serpina3n"]          # the dataset-specific ageing gene of synth_benayoun
    genes += MITO + NUCLEAR + RIBO
    genes += [f"Gene{i:04d}" for i in range(n_extra)]
    return list(dict.fromkeys(genes))


def shared_biology(rng, genes):
    n = len(genes)
    base = rng.gamma(0.6, 1.0, size=n) * 0.4
    responsive = {g for fam in (AGEING, SURGERY, *OVULATION.values()) for ch in fam.values() for g in ch} | {"Serpina3n"}
    for g in responsive:
        base[genes.index(g)] = max(base[genes.index(g)], 0.5)
    programme = {t: rng.choice(n, 60, replace=False) for t in TYPES}
    return base, programme


def resolve_response(spec):
    if spec == "AGEING_0.8":
        out = {t: {g: float(f) ** 0.8 for g, f in ch.items()} for t, ch in AGEING.items()}
        out["Stromal"]["Serpina3n"] = 3.0
        return out
    return spec or {}


def make_dataset(rng, spec: dict, genes: list[str], base, programme, n_per_animal: int, batch_genes):
    import anndata as ad

    n_genes = len(genes)
    gi = {g: i for i, g in enumerate(genes)}
    mod = MODALITY[spec["modality"]]
    vocab = VOCAB[spec["vocab"]]
    obs_rows, blocks = [], []
    for raw_value, cond, n_animals, response in spec["conditions"]:
        response = resolve_response(response)
        for a in range(n_animals):
            animal = f"{spec['name']}_{cond}_{a + 1}"
            n_cells = int(rng.integers(int(n_per_animal * 0.85), int(n_per_animal * 1.15)))
            labels = rng.choice(TYPES, size=n_cells, p=[SHARE[t] for t in TYPES])
            rates = np.tile(base, (n_cells, 1))
            for t in TYPES:
                m = labels == t
                for g in PANELS[t]:
                    rates[m, gi[g]] *= 12.0
                rates[np.ix_(m, programme[t])] *= 2.5
                for g, fold in response.get(t, {}).items():
                    rates[m, gi[g]] *= fold
            for g, fold in response.get("*", {}).items():
                rates[:, gi[g]] *= fold
            rates[:, batch_genes] *= 1.5
            animal_genes = rng.choice(n_genes, 40, replace=False)
            rates[:, animal_genes] *= rng.uniform(0.8, 1.25, size=40)
            # modality tells: mito, nuclear lncRNA and ribosomal shares
            special = [gi[g] for g in MITO + NUCLEAR + RIBO]
            rates[:, special] = 0.0
            other = rates.sum(axis=1, keepdims=True)
            rest = 1.0 - mod["mito"] - mod["nuclear"] - mod["ribo"]
            for group, share in ((MITO, mod["mito"]), (NUCLEAR, mod["nuclear"]), (RIBO, mod["ribo"])):
                idx = [gi[g] for g in group]
                per_cell = other * share / rest
                rates[:, idx] = per_cell / len(idx) * rng.uniform(0.6, 1.4, size=(n_cells, len(idx)))
            depth = rng.lognormal(np.log(mod["depth"]), 0.35, size=n_cells) * rng.uniform(0.85, 1.15)
            probs = rates / rates.sum(axis=1, keepdims=True)
            counts = rng.poisson(probs * depth[:, None]).astype(np.float32)
            blocks.append(sp.csr_matrix(counts))
            author = [rng.choice(vocab[t]) for t in labels]
            for i in range(n_cells):
                obs_rows.append({"animal": animal, spec["label_col"]: author[i], spec["cond_col"]: raw_value,
                                 "true_type": labels[i]})
    X = sp.vstack(blocks).tocsr()
    obs = pd.DataFrame(obs_rows)
    obs.index = [f"{spec['name']}_c{i:06d}" for i in range(len(obs))]
    var = pd.DataFrame({"gene_ids": [f"ENSMUSG{i:011d}" for i in range(n_genes)]}, index=genes)
    return ad.AnnData(X=X, obs=obs, var=var)


def write_side_files(out: Path, specs: list[dict]):
    repo_cw = yaml.safe_load((REPO_ROOT / "config" / "crosswalk.yaml").read_text(encoding="utf-8"))
    repo_bench = yaml.safe_load((REPO_ROOT / "config" / "benchmark.yaml").read_text(encoding="utf-8"))
    datasets = []
    for spec in specs:
        entry = {"name": spec["name"], "enabled": True, "tier": "synthetic", "species": "mouse",
                 "suspension_type": None,          # inferred from the tells, on purpose
                 "assay": "synthetic", "strain": "synthetic", "sex": "female",
                 "source": {"type": "local"}, "path": str(out / "data" / "raw" / spec["name"] / f"{spec['name']}.h5ad"),
                 "label_key": spec["label_col"], "gene_key": None, "animal_key": "animal",
                 "condition": {"key": spec["cond_col"], "map": spec["map"], "family": spec["family"],
                               "control": spec["control"]},
                 "status": "ready", "notes": "synthetic"}
        if spec.get("age_days"):
            entry["condition"]["age_days"] = spec["age_days"]
        if spec.get("time_h"):
            entry["condition"]["time_h"] = spec["time_h"]
        if spec.get("held_out"):
            entry["held_out"] = True
        datasets.append(entry)
    (out / "corpus.yaml").write_text(yaml.safe_dump({"datasets": datasets}, sort_keys=False), encoding="utf-8")

    cw = dict(repo_cw)
    cw["datasets"] = {"synth_isola": repo_cw["datasets"]["isola2024"],
                      "synth_huang": repo_cw["datasets"]["morris2022"],
                      "synth_benayoun": {},                       # the global map, with one unmapped label
                      "synth_mouse_comp": repo_cw["datasets"]["mouse_comp_30d"]}
    (out / "crosswalk.yaml").write_text(yaml.safe_dump(cw, sort_keys=False, allow_unicode=True), encoding="utf-8")

    bench = dict(repo_bench)
    bench["t1"] = dict(repo_bench["t1"], min_cells_per_type=10, max_cells_per_class=800,
                       test_only=["synth_mouse_comp"])
    bench["t3"] = dict(repo_bench["t3"], min_cells_hard=40, min_de_genes_to_score=4)
    bench["t3"]["contexts"] = [
        {"name": "synth_isola_ageing", "dataset": "synth_isola", "family": "ageing", "control": "age_3mo",
         "perturbed": ["age_9mo"]},
        {"name": "synth_benayoun_ageing", "dataset": "synth_benayoun", "family": "ageing", "control": "young",
         "perturbed": ["old"]},
        {"name": "synth_huang_ovulation", "dataset": "synth_huang", "family": "ovulation", "control": "hcg_0h",
         "perturbed": ["hcg_4h", "hcg_12h"]},
        {"name": "synth_mouse_comp_surgery", "dataset": "synth_mouse_comp", "family": "surgery", "control": "sham",
         "perturbed": ["ooph_30d"], "held_out": True},
    ]
    (out / "benchmark.yaml").write_text(yaml.safe_dump(bench, sort_keys=False), encoding="utf-8")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--out", required=True)
    parser.add_argument("--n-per-animal", type=int, default=350)
    parser.add_argument("--n-extra-genes", type=int, default=800)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args(argv)

    rng = np.random.default_rng(args.seed)
    out = Path(args.out)
    genes = gene_universe(args.n_extra_genes)
    base, programme = shared_biology(rng, genes)
    for spec in DATASETS:
        batch_genes = rng.choice(len(genes), 60, replace=False)
        adata = make_dataset(rng, spec, genes, base, programme, args.n_per_animal, batch_genes)
        folder = out / "data" / "raw" / spec["name"]
        folder.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(folder / f"{spec['name']}.h5ad")
        print(f"  {spec['name']}: {adata.n_obs:,} {spec['modality']}s x {adata.n_vars:,} genes, "
              f"{adata.obs['animal'].nunique()} animals -> {folder}")
    write_side_files(out, DATASETS)
    (out / "README.md").write_text(__doc__, encoding="utf-8")
    print(f"  side files: {out / 'corpus.yaml'}, crosswalk.yaml, benchmark.yaml")


if __name__ == "__main__":
    main()
