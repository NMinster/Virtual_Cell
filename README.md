# Virtual_Cell

An ovary-specific foundation model, with the Arc Institute **Virtual Cell
Challenge 2026** as the motivation and the evaluation template rather
than a competition to enter. The experimental anchor is the
[Mouse_Compensatory](https://github.com/NMinster/Mouse_Compensatory)
study (oophorectomy vs sham mouse ovary: single-nucleus RNA-seq plus
Visium HD), which this project treats as a held-out test context.

This repo is laid out like Mouse_Compensatory: one YAML of parameters,
side files for the things that grow, one stage per step, a provenance
json before every output, the same three Biowulf scripts, and a synthetic
dataset the tests run end to end.

## Documents

| Document | What it is |
| --- | --- |
| [`docs/ovary_foundation_model_research.md`](docs/ovary_foundation_model_research.md) | The research memo (2026-11-03): the state of single-cell foundation models, the measured ovary gap (Census, scBaseCount, Tahoe-100M, the model corpora), the public ovary data inventory, the Challenge's 2026 scoring contract, cross-context perturbation methods, licences, compute and scaling, the ovary tools gap. Measured tables in [`docs/data/`](docs/data/). |
| [`docs/project_plan.md`](docs/project_plan.md) | The plan: corpus, held-out-context benchmark, model arms, gates and kill rules, the spatial arm, how the own data plug in, a 12-month timeline. |
| [`docs/preregistration.md`](docs/preregistration.md) | The rules, frozen before any model is trained: hypotheses, tasks, members, baselines, replicate, exclusions, margins, deviations log. **Draft until `frozen_at` is filled.** |
| [`docs/first_results_checklist.md`](docs/first_results_checklist.md) | Step by step from a clean clone to the three first tables (ceilings, baseline scores, label transfer), with the commands and what good looks like. |

## The gap, in three numbers

In the newest CELLxGENE Census (LTS 2025-11-08; no 2026 build exists)
human ovary is 373,467 unique cells, 0.375 % of the corpus and rank 31
of 71 tissues, half of it tumour; mouse ovary is 17,408 cells, 0.083 %,
rank 27 of 36, from two studies, with no nuclei, no spatial data and
nothing older than 9 months; Tahoe-100M's 50 cell lines include no
ovarian line. Well over a million mouse ovary cells, nuclei and spots and
more than half a million human ovary nuclei exist in GEO and SRA outside
the Census, which is why the plan starts with curation.

## The stages

| Stage | Reads | Writes |
| --- | --- | --- |
| `inspect` | `config/corpus.yaml` | what is on disk, QC'd, or still `fill_in` |
| `fetch` | the manifest's sources | `data/raw/<dataset>/` (+ `download.json`); Census, plain URLs, scBaseCount, local |
| `qc` | each raw dataset | `data/corpus/qc/<dataset>.h5ad` with the obs schema; `tables/qc/qc_summary.csv`, `<dataset>_dropoff.csv`; the modality inferred from the tells when the manifest says null |
| `harmonise` | the QC shards + `config/crosswalk.yaml` | `data/corpus/corpus.h5ad` (pooled), `heldout_<dataset>.h5ad`; `tables/harmonise/manifest.csv`, `crosswalk_review.csv` (the table for the lab), `composition_by_dataset.csv` |
| `bench` | the corpus + `config/benchmark.yaml` | the T3 bundles (`results/bench/bundles/`, `reference.json`), the real DE tables, **`tables/bench/ceilings.csv`**: b0, b1 and the cross-animal replicate r per member, per (context, condition, compartment) |
| `predict --model b0\|b1\|stub` | the bundles | `results/predictions/<model>.h5ad` |
| `score --pred PATH --name NAME` | a prediction + the bundles | `tables/bench/scores_<name>.csv`, `_summary.csv`: the six members rescaled 0 = baseline, 1 = replicate (`--scorer cell-eval2` wraps Arc's scorer, untested) |
| `t1 [--rep embedding.parquet]` | the corpus | `tables/bench/t1_<name>.csv`, `_per_type.csv`: leave-one-atlas-out macro-F1 and rare-type F1 |
| `fm --tool stub\|scprint2\|uce\|transcriptformer` | the corpus | `results/fm/<model>/embedding.parquet` + provenance (the zero-shot arm; the real-model branches are ported untested from Mouse_Compensatory) |
| `fm-compare` | every embedding found | `tables/bench/fm_compare.csv` |
| `train --arm a\|b` | the corpus | the recipe only: `results/train/arm_<x>/{config.yml,train.sbatch,README.md}` |
| `train --arm c` | the corpus | scVI then scANVI (needs scvi-tools): `results/train/arm_c/embedding.parquet` |
| `report` | every table above | `results/report/first_results.md` |
| `all` | — | `qc harmonise bench predict score t1 report` with the stub model |

Every stage writes `results/provenance/<stage>_<utc>.json` first.

## Quick start

```bash
git clone https://github.com/NMinster/Virtual_Cell.git && cd Virtual_Cell
bash biowulf/setup_env.sh              # conda env `ovary_fm` + editable install
bash biowulf/run.sh pytest -q          # the synthetic corpus through the whole chain (~2 min)

export OVARY_FM_DATA=/data/$USER/ovary_fm/data
export OVARY_FM_RESULTS=/data/$USER/ovary_fm/results
bash biowulf/run.sh inspect
bash biowulf/submit.sh fetch --dataset isola2024
bash biowulf/submit.sh first           # qc harmonise bench, then predict score t1 report, as two chained jobs
```

The checklist walks through filling the manifest for the first three
datasets, the crosswalk review, the ceilings table and the baselines.

## Layout

```
config/          params.yaml, corpus.yaml (manifest), crosswalk.yaml, benchmark.yaml, train/ (recipe templates)
ovary_fm/        config, io, fetch, qc, harmonise, fm, report, cli; bench/ (metrics, baselines, runner, t1); models/ (stub, recipes, scvi_arm)
scripts/         make_synthetic.py, run_fm_embed.py, scbasecount_metadata.py
R/               rds_to_mtx.R (a Seurat / SCE object -> an mtx folder qc reads)
biowulf/         setup_env.sh, run.sh, submit.sh
environment*.yml one conda env per stack: ovary_fm, ovary_fm_arc (cell-eval2), ovary_fm_scprint2, ovary_fm_tahoe, ovary_fm_scvi, ovary_fm_uce, ovary_fm_transcriptformer
tests/           test_smoke.py: units + the synthetic end-to-end run
docs/            the memos and the tracked data tables
```

`data/`, `results/`, `logs/`, `research_notes/` and `reports/` are not
tracked (see `.gitignore`); the memos under `docs/` are.

## The Challenge, as the template

The 2026 Challenge asks participants to predict how six unseen cell lines
respond to gene knockdowns given only each line's non-targeting control
cells, and scores six metrics each rescaled so that 0 is the
context-mean baseline and 1 is a split-half biological replicate. The
benchmark here copies that contract onto ovarian contexts (a compartment
in a dataset's baseline state) and conditions (hours after hCG, age,
chemical menopause, genotype, drug, surgery), with whole datasets held
out and the condition-defining genes excluded the way Arc excludes the
knockdown targets.
