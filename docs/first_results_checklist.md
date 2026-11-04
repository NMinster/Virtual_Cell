# Checklist: from a clean clone to the first results

Nov 4, 2026. Step by step, with the command, what it writes and what
good looks like. The commands are for Biowulf; on a laptop drop the
`biowulf/` wrappers and run `ovary-fm ...` in the activated env. Times
are rough. Tick the box in a copy of this file, or in the issue tracker.

The "first results" this leads to are three tables, none of which needs a
model:

1. **the ceilings table** (`tables/bench/ceilings.csv`): per compartment
   and condition, how reproducible an ovarian response is between animals,
   which bundles pass the gates, and the two baselines' raw values;
2. **the baseline scores** (`tables/bench/scores_b0*.csv`,
   `scores_b1*.csv`): how far retrieval-and-transfer gets on each held-out
   context on the 0-to-1 scale;
3. **the label-transfer table** (`tables/bench/t1_logreg.csv`,
   `t1_rep_*.csv`): macro-F1 and rare-type F1 per held-out atlas for the
   gene-space baseline and, once they run, the generalists and scANVI.

Everything a model later produces is read against those three.

## 0. Before anything (10 min)

- [ ] You have a Biowulf account, `/data/$USER` has room (the corpus will be
      tens of GB; weights and the Census cache another 20-50 GB), and you
      can reach GitHub from a login node.
- [ ] Clone both repos side by side (the manifest points at
      `../Mouse_Compensatory/results/objects/subtypes.h5ad` for the own data):

```bash
cd /data/$USER
git clone https://github.com/NMinster/Virtual_Cell.git
git clone https://github.com/NMinster/Mouse_Compensatory.git   # if not already there
cd Virtual_Cell
export OVARY_FM_DATA=/data/$USER/ovary_fm/data
export OVARY_FM_RESULTS=/data/$USER/ovary_fm/results
```

Put the two exports in your login shell's profile; `submit.sh` bakes them
into every job script.

Clone on Biowulf itself, not on a Windows machine and then copy. A
checkout made with `core.autocrlf=true` arrives with `\r` line endings
and every shell script dies at its first line with
`set: pipefail\r: invalid option name` (the terminal garbles the message
because of the carriage return). The repo's `.gitattributes` now forces
LF on checkout and the scripts re-run themselves cleaned when they
detect `\r`, but if you ever see that message:

```bash
git config --global core.autocrlf        # should be empty or `input` on Biowulf
sed -i 's/\r$//' biowulf/*.sh R/*.R      # strip the carriage returns in place
git pull && git checkout -- .            # or re-checkout with the LF attributes
```

## 1. Build the environment and prove the plumbing (30-60 min)

- [ ] `bash biowulf/setup_env.sh` (in an `sinteractive --mem=16g
      --cpus-per-task=4` if you are not already in a job). It creates the
      `ovary_fm` conda env under `/data/$USER/conda` and installs the
      package editable.
- [ ] `bash biowulf/run.sh pytest -q` → every test green, about two minutes.
      The suite builds the synthetic four-dataset corpus and runs the whole
      chain on it.
- [ ] Look at what the chain produces before touching real data:

```bash
bash biowulf/run.sh synth --out /data/$USER/ovary_fm/fake
bash biowulf/run.sh all \
    --set paths.corpus_manifest=/data/$USER/ovary_fm/fake/corpus.yaml \
    --set paths.crosswalk=/data/$USER/ovary_fm/fake/crosswalk.yaml \
    --set paths.benchmark=/data/$USER/ovary_fm/fake/benchmark.yaml \
    --data /data/$USER/ovary_fm/fake/data --results /data/$USER/ovary_fm/fake/results
less /data/$USER/ovary_fm/fake/results/report/first_results.md
```

What good looks like: the modality is inferred (cell / nucleus) from the
tells and printed; one label (`mystery cells`) is reported as unmapped;
the ceilings table has cross-animal replicate levels; b0 scores 0 by
construction; b1 beats b0 on the DE members where a source exists; T1 is
near-perfect (the synthetic types are easy on purpose).

## 2. Verify the numbers the plan leans on (done for five records on 2026-11-04)

The GEO records were read on 2026-11-04 and the manifest carries what they
say. Settled:

- Huang 2025 (GSE294534): six scRNA libraries, 0 / 4 / 12 h post-hCG,
  two pools (ABC, DEF) of three mice per time point; the deposit is one
  5.5 GB Seurat object, `GSE294534_10x_update.rds`. PLOS Biology
  23(6):e3003193.
- Wang/Suh (GSE267729): 18 per-sample matrices (young regular, peri-
  estropausal regular vs irregular, post-estropausal, SPiDER-βGal-high
  sorted cells), public since 2026-07-01, PNAS 123(30):e2600323123. No
  annotations deposited: the major types are transferred in `harmonise`.
- The ex vivo ovulation series is GSE255274 (Thomas, Marx, Penir, Schuh;
  Nature Cell Biology 2024): 28 per-sample matrices, in vivo and in vitro,
  0 / 3 / 6 / 9 / 12 h post-hCG, replicates. Ovulating follicles, not
  whole ovaries. The second ovulation dataset: b1 for Huang exists.
- GSE329417: 439,011 nuclei, 64 samples, 28 postmenopausal donors
  (Byrne, Martin, Soygur et al. 2026, Melov lab); a final annotated h5ad
  (5.4 GB) is deposited. GSE267315: 57,970 nuclei, ages 12-54, plus
  Stereo-seq (Zhang, Guo et al. 2026, Aging Cell). Both are Tier C and
  stay disabled until the species ablation.
- Gaylord 2025 mouse data are GSE290742 (Science 390:eadx0659): three
  young and three aged mice, a 10x run plus a Smart-seq2 plate of large
  oocytes each, twelve csv matrices (215 MB), no annotations. The human
  half is GSE290743. Enabled, with its own ageing context.
- **GSE192722 is not Morris 2022.** It is Lustgarten et al., human theca
  cells (Weill Cornell). The number came from a docstring in
  Mouse_Compensatory and is wrong there too.
- Morris 2022 is not on GEO: the deposit is OSF project 924fz (David
  Pépin, "Single cell sequencing of the mouse ovary in diverse
  reproductive states"), one 301.8 MB inDrop count matrix of > 34,000
  cells across randomly cycling, proestrus, estrus, metestrus, diestrus
  and 10 d post-partum lactating and non-lactating ovaries. The manifest
  entry points at it; the OSF files pane asks for a sign-in (step 4).
- scPRINT's model card says Apache-2.0 (read 2026-11-04): a derived
  checkpoint can be released under Apache-2.0 with attribution. The code
  stays GPL-3. Nothing to ask anyone.

- Biowulf (`batchlim`, `freen`, 2026-11-04): the gpu partition has 320
  A100s on 80 nodes of 4 GPUs, 64 CPUs, 243 GB and 3.2 TB SSD each; 56
  GPUs per user; walltime 2 h by default, 10 days at most. Whole nodes
  are rarely free (0 of 80 that day, 13 single GPUs free), so 4-GPU jobs
  queue. `quick` has 12 of the same nodes with a 4 h cap, right for the
  single-GPU tests; `gpuh200` offers 8 GPUs per user for a day. `freen`
  does not print memory per card; `nvidia-smi` in an `sinteractive
  --gres=gpu:a100:1` says whether they are 40 or 80 GB. The arm A and B
  sbatch scripts ask for 230 GB and 32 CPUs on one node, the most the
  node allows.

Still to do, no compute: nothing. The remaining items are the searches below.

Run the searches the memo never got to, and record "found / not found"
with the query in the memo's section 8: PCOS, POI and chemotherapy mouse
scRNA-seq; KGN / COV434 / HGrC1 pooled screens; FSH, androgen and AMH
single-cell treatments; reproductive-tissue foundation-model abstracts
from 2025-2026 meetings.

## 3. The manifest (filled on 2026-11-04; three things left)

`config/corpus.yaml` is filled for Isola 2024 (Census), Huang 2025,
Wang/Suh, the Schuh series and Gaylord 2025 (all `enabled: true`), and for
the two human atlases (`enabled: false`). What the download will tell you:

- [ ] `bash biowulf/run.sh inspect` → the five mouse datasets are enabled;
      the `fill_in` statuses are what step 5 settles.
- [ ] Huang 2025: `label_key` and `animal_key` come from the Seurat
      object's metadata columns, printed by `qc` as `obs columns` after
      the conversion in step 4. The paper's 16 named clusters are the
      labels; if there is a per-mouse (hashtag) column use it as
      `animal_key`, otherwise the library (`orig.ident`) is the replicate.
- [ ] Wang/Suh, Schuh and Gaylord: the `condition.patterns` (and Gaylord's
      `animal_regex`, which folds the 10x and plate libraries of one mouse
      into one animal) are written against the sample titles GEO shows;
      the real file names may differ. `inspect` prints the `gsm_sample`
      values once the tars are unpacked, and `qc` prints any value no
      pattern matches. Adjust the regexes, rerun.

## 4. Fetch (one batch job per dataset; minutes to an hour each; ~11 GB)

```bash
bash biowulf/submit.sh fetch --dataset isola2024              # Census: needs the extra below, once
bash biowulf/submit.sh fetch --dataset huang2025              # 5.5 GB .rds
bash biowulf/submit.sh fetch --dataset wang_suh_estropause    # 1.2 GB tar, unpacked by fetch
bash biowulf/submit.sh fetch --dataset schuh2024_ovulation    # 4.4 GB tar, unpacked by fetch
bash biowulf/submit.sh fetch --dataset gaylord2025_mouse      # 215 MB tar of csv matrices
bash biowulf/submit.sh fetch --dataset morris2022             # OSF: the public API, or the message to download by hand
bash biowulf/submit.sh python scripts/scbasecount_metadata.py --species mouse
```

Morris 2022 sits on OSF (project 924fz) and its files pane asks for a
sign-in. If the fetch job prints "did not answer the public API", sign in
at https://osf.io/924fz/files with ORCID or an institutional login,
download the 301.8 MB file, put it under `$OVARY_FM_DATA/raw/morris2022/`
(a `.rds` is converted as for Huang; a csv or h5ad is read as is), then
set `enabled: true` on `morris2022` and fill its `label_key`, `animal_key`
and condition column from what `qc` prints.

Yes, download these five now; they are the whole mouse Tier A that is
reachable today. Leave the human atlases (another 9 GB) until the species
ablation. The Census route needs the proxy and should never run on the
login node (the job script carries both).

- [ ] `data/raw/<name>/download.json` exists for each, with sizes and
      sha256; `GSE267729_RAW/` and `GSE255274_RAW/` are unpacked beside
      their tars.
- [ ] Huang 2025 is a Seurat object; convert it once (the R module's
      Seurat; if `library(Seurat)` fails, install it into `R_LIBS_USER`):

```bash
module load R
Rscript R/rds_to_mtx.R $OVARY_FM_DATA/raw/huang2025/GSE294534_10x_update.rds $OVARY_FM_DATA/raw/huang2025/
head -3 $OVARY_FM_DATA/raw/huang2025/meta.csv      # the columns for label_key / animal_key
```

- [ ] `bash biowulf/run.sh inspect` → `on_disk: True` for all five, the
      per-sample counts (18, 28 and 12) and the `gsm_sample` values for
      the three tar datasets.
- [ ] `data/scbasecount/mouse_ovary_candidates.csv` lists the ovary
      experiments of scBaseCount (132 experiments, 1.1 M cells, 86,346
      nuclei in six experiments on 2026-11-04), including a 2025 in vivo
      hCG series, SRX28395245-50, to resolve later as a third ovulation
      dataset.

What good looks like: four raw folders, the candidates table, and
`inspect` showing the sample names you can match the patterns against.

## 5. QC, one dataset at a time (minutes each)

```bash
bash biowulf/submit.sh qc --dataset isola2024
```

- [ ] Read the log: `obs columns` tells you the real `label_key`,
      `animal_key` and condition column; `condition values no pattern
      matches` tells you which regexes to fix; fix them in the manifest and
      rerun until the summary row is right.
- [ ] `tables/qc/qc_summary.csv`: modality agrees with the manifest (the
      tells are printed either way); `n_animals` is the design's (Isola 8;
      Huang 6 libraries unless a mouse column exists; Wang/Suh 18; Schuh
      28; Gaylord 6 animals from 12 libraries); `conditions` are the
      mapped names (`age_3mo;age_9mo`, `young;aged`,
      `hcg_0h;hcg_4h;hcg_12h`, `young;peri_regular;peri_irregular;
      post_estropausal;senescent_sorted`, `invivo_hcg_0h;...`).
- [ ] `tables/qc/<name>_dropoff.csv`: no step removes more than you can
      explain; nuclei get the 5 % mito cap, cells 20 %.
- [ ] Repeat for `huang2025`, `wang_suh_estropause`, `schuh2024_ovulation`,
      `gaylord2025_mouse`. Gaylord's plate libraries are Smart-seq2 read
      counts, not UMIs; they carry `assay: Smart-seq2` in obs so a model
      can token them apart, and they are mostly oocytes.

## 6. Harmonise and review the crosswalk (an hour, plus the lab's review)

```bash
bash biowulf/submit.sh harmonise
```

- [ ] `tables/harmonise/manifest.csv`: `crosswalk_coverage` per dataset.
      Isola is 1.0 (its block is written); Huang is low until its clusters
      are mapped in step 6; Wang/Suh and Schuh show `labels: transferred`
      and their coverage is the share of cells above the transfer
      confidence (0.5), with `crosswalk_status` `transferred` per cell.
      Datasets with transferred labels are never T1 train or test sets.
- [ ] `tables/harmonise/crosswalk_review.csv`: every author label with its
      cell count and mapping. Map Huang's clusters in
      `config/crosswalk.yaml` under `datasets.huang2025` from the paper's
      marker table; rerun.
- [ ] Every `cl_id: verify` term: look the label up at
      https://www.ebi.ac.uk/ols4 and replace `verify` with the id (or
      `null` with a note if the term does not exist).
- [ ] Send `crosswalk_review.csv` and the marker evidence to Raghu's lab;
      the reviewed table is the one that gets committed.

What good looks like: `data/corpus/corpus.h5ad` with three datasets,
coverage ≥ 0.9 everywhere, `composition_by_dataset.csv` with granulosa the
largest compartment in each.

## 7. The first result: the ceilings table (a 128 GB job; minutes to an hour)

```bash
bash biowulf/submit.sh bench
```

- [ ] `tables/bench/ceilings.csv`: one row per (context, condition,
      compartment). For Isola (4 + 4 animals) and Huang (3 per time
      point) `replicate_level` must read `animals`; `n_de_real` per
      compartment is the first measurement of how much each compartment
      responds; `gate_nmae` says where the fold-change member is scorable.
- [ ] `tables/bench/de/*.csv`: the real DE tables; open the granulosa
      hcG 4 h one and check that the excluded immediate-early genes carry
      large fold changes (they are the leak the rule removes) and that
      known LH targets outside the list are there.
- [ ] Note which contexts have `n_b1_sources = 0`: with only Isola in the
      ageing family and only Huang in ovulation, b1 exists nowhere yet.
      That is expected; step 9 fixes it.

Read it as: r per member is the ceiling a perfect model could reach at
half depth; b0 is where "no effect" sits; the distance between them is
the room the benchmark has in that compartment. `scored` says whether
the compartment responded enough (≥ 5 real DE genes after exclusions)
to be aggregated at all; the rest are listed with their ceilings and
left out of the scores and the PDS pools, and the paper says so. Expect
most compartments outside granulosa, theca and stroma to sit below the
gate in a hormone time course, and the immune and stromal compartments
to carry ageing.

## 8. The second result: the baselines scored (minutes)

```bash
bash biowulf/submit.sh chain predict score     # predict --model stub by default; see below
bash biowulf/run.sh predict --model b0
bash biowulf/run.sh score --pred $OVARY_FM_RESULTS/predictions/b0.h5ad --name b0
```

- [ ] `scores_b0_summary.csv`: 0.000 on every context (the scale is
      anchored on the control copy; this is the scorer's consistency check).
- [ ] Once b1 exists (step 9): `predict --model b1`, `score ... --name b1`;
      b1's row is the bar every model has to clear. Expect it well above 0
      on the DE members in every compartment with ≥ 10 real DE genes; where
      it is near 0 or negative the two datasets disagree on the response,
      which is itself a finding to write down.

## 9. b1 already exists for both families; a third dataset each is the next step

With Isola, Wang/Suh and Gaylord in the ageing family and Huang and Schuh
in the ovulation family, every scored bundle of those contexts has a
transfer source and the b1 row is the bar. To make the bar harder:

- [ ] Ageing: `galligos2026` or `benayoun_vcd` (fill the source, fetch, qc,
      harmonise, bench again).
- [ ] Ovulation: the 2025 in vivo hCG series the scBaseCount table surfaced
      (SRX28395245-50, "0hr post-hCG", "hCG treatment", 14-18 K cells each):
      resolve its GEO series and time points, add it as a third ovulation
      dataset.
- [ ] Cycle: Morris 2022 once the OSF file is in place (step 4); it is the
      only estrous-staged whole-ovary atlas and the `morris_cycle` context
      is already defined (diestrus against the other stages and the two
      post-partum states).

## 10. The third result: label transfer (minutes; an hour with the generalists)

```bash
bash biowulf/submit.sh t1
```

- [ ] `tables/bench/t1_logreg.csv`: macro-F1 and rare-type F1 per
      held-out atlas for the gene-space baseline. Expect high macro-F1
      and a lower rare-type F1 (oocytes, luteal, theca) where the
      vocabularies disagree; `t1_logreg_per_type.csv` says which type.
- [ ] The zero-shot arm, one model (scPRINT-2 small is the cheapest):

```bash
bash biowulf/setup_env.sh --fm scprint2
# once, inside the env: the ontologies scPRINT-2 reads and the public checkpoint
bash biowulf/run.sh python -c "print()" >/dev/null   # (just to source conda); then:
conda activate ovary_fm_scprint2
lamin init --storage /data/$USER/ovary_fm/lamin --name ovary_fm --modules bionty
scprint2 easy_setup
hf download jkobject/scPRINT small-v2.ckpt --local-dir $HF_HOME/scprint
conda deactivate
bash biowulf/submit.sh fm --tool scprint2 --set fm.models.scprint2.checkpoint=$HF_HOME/scprint/small-v2.ckpt
bash biowulf/submit.sh fm-compare                  # t1 inside every embedding found
```

      `lamin init` needs the proxy (login node, or `run.sh python` to source
      it); the instance lives under `/data/$USER`, not the home quota.

      `tables/bench/fm_compare.csv` puts the generalist beside the
      logistic regression per held-out atlas. The September memo predicts
      it will not win; if it does, say so loudly.
- [ ] The primary baseline, scANVI:

```bash
bash biowulf/setup_env.sh --train scvi
bash biowulf/submit.sh train --arm c
bash biowulf/run.sh t1 --rep $OVARY_FM_RESULTS/train/arm_c/embedding.parquet
```

## 11. Two single-GPU checks before any training budget is asked for (half a day)

- [ ] `bash biowulf/setup_env.sh --train scprint2` adds the dev clones
      (scPRINT-2, scDataLoader, GRnnData, benGRN under
      `/data/$USER/ovary_fm/src`): pretraining lives in the repo, not in
      the PyPI package. Register the QC shards as a laminDB collection
      named `ovary_fm_corpus` (and the Isola shard alone as
      `ovary_fm_isola_only`) with scDataLoader's Preprocessor, as their
      `docs/pretrain.md` describes.
- [ ] `bash biowulf/run.sh train --arm a` writes
      `results/train/arm_a/{config.yml,train.sbatch,README.md}`. The config
      is passed as the last `--config` to `scprint2 fit` after their
      `base_v2.yml` and `pretrain_medium.yml`; check its keys against the
      clone's `config/` folder, point their periodic `test()` at the test
      files on hf.co/jkobject/scPRINT-2 or disable it (their issue #12),
      then run a few hundred steps on the Isola collection on one GPU in
      the `quick` partition (`sbatch --partition=quick --gres=gpu:a100:1
      --mem=60g --cpus-per-task=16 --time=04:00:00`, or `sinteractive`
      with the same resources). Record the step time and the memory; that
      is the number the 200-500 A100-hour estimate is checked against,
      and 60 GB per GPU is what a shared node leaves.
- [ ] `bash biowulf/setup_env.sh --arc`, then in the `ovary_fm_arc` env
      `cell-eval2 run --help`. Align the flags in
      `ovary_fm/bench/runner.py::_score_cell_eval2` with what you see,
      score one context with `--scorer cell-eval2` and with the `lite`
      scorer, and write the differences per member into the deviations
      log of `docs/preregistration.md`. From then on, report Arc's numbers
      where both exist.

## 12. Write it down (an hour)

```bash
bash biowulf/run.sh report
```

- [ ] `results/report/first_results.md` carries all three tables. Copy it
      to `docs/first_results_<date>.md` and commit it with the config
      files as they were; the provenance json under `results/provenance/`
      names the git revision each table came from.
- [ ] Freeze `docs/preregistration.md` (fill `frozen_at:` with the commit
      hash) before step 13 or any model arm.

## 13. The own data (half a day)

- [ ] Enable `mouse_comp_30d` in the manifest (path to the final
      Mouse_Compensatory object), `qc --dataset mouse_comp_30d`,
      `harmonise`, `bench`: the surgery context appears with its own
      2 + 2-mouse ceiling (`replicate_level: animals`) and, almost
      certainly, several compartments below the gates. That row is the
      honest statement of what the 30 d data can score, and the reason
      the 120 d cohort is designed with ≥ 400 nuclei per compartment per
      arm and three mice per arm if the budget allows.
- [ ] `t1` scores it as a test-only set: how well the public atlases'
      labels transfer onto nuclei the lab annotated by hand.

## Where things land

| Path | What |
| --- | --- |
| `$OVARY_FM_DATA/raw/<dataset>/` | downloads, `download.json` |
| `$OVARY_FM_DATA/corpus/qc/<dataset>.h5ad` | QC'd shards with the obs schema |
| `$OVARY_FM_DATA/corpus/corpus.h5ad`, `heldout_<dataset>.h5ad` | the corpus and the held-out shards |
| `$OVARY_FM_RESULTS/tables/qc/`, `harmonise/`, `bench/` | the tables named above |
| `$OVARY_FM_RESULTS/bench/bundles/`, `reference.json` | the T3 bundles, splits and reference values |
| `$OVARY_FM_RESULTS/predictions/<model>.h5ad` | prediction objects |
| `$OVARY_FM_RESULTS/fm/<model>/`, `train/<arm>/` | embeddings and recipes |
| `$OVARY_FM_RESULTS/report/first_results.md` | the page to read |
| `$OVARY_FM_RESULTS/provenance/*.json` | git revision and live config per stage |
| `logs/<stage>_<jobid>.out` | the job logs |
