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

## 2. Verify the numbers the plan leans on (half a day, no compute)

The research memo could not open the journals, so these are tagged as
search-snippet claims. Open the paper or the GEO record and write the
confirmed number into `config/corpus.yaml` `notes:` (and fix the memo if
it is wrong):

- [ ] Morris 2022: cell count (the paper says > 34,000; the reused object
      has 28,014), GEO GSE192722 file list and format.
- [ ] Huang 2025 (PLOS Biology): 25,714 cells, 3 animals per time point,
      GSE294534 file list.
- [ ] Wang/Suh estropausal transition: accession GSE267729 and the design.
- [ ] Gaylord 2025: which accession holds the mouse object.
- [ ] GSE329417 (439,011 nuclei), GSE267315, the Slide-seq ageing atlas
      (610,620 spots): sizes and whether raw counts are deposited.
- [ ] Nature Cell Biology 2024 ex vivo ovulation series: the GEO series
      behind scBaseCount SRX23566711-38 and the time-point mapping.
- [ ] scPRINT-2's checkpoints are public on the hub (`hf download
      jkobject/scPRINT small-v2.ckpt`), so nothing is needed to use them.
      The open question is only what licence a *derived* checkpoint can
      be released under at month 12: open https://huggingface.co/jkobject/scPRINT
      in a browser and read the model card's licence field (unreachable
      from the research sandbox). If it names a licence, follow it; if it
      is empty, a GitHub issue on cantinilab/scPRINT-2 asking which
      licence covers the weights is the quick public route. Not a blocker
      for anything before the release.
- [ ] Biowulf: `freen | grep a100` and `batchlim` for the current A100
      memory mix and walltime limits.

Run the searches the memo never got to, and record "found / not found"
with the query in the memo's section 8: PCOS, POI and chemotherapy mouse
scRNA-seq; KGN / COV434 / HGrC1 pooled screens; FSH, androgen and AMH
single-cell treatments; reproductive-tissue foundation-model abstracts
from 2025-2026 meetings.

## 3. Fill the manifest for the first three datasets (1-2 h)

- [ ] `bash biowulf/run.sh inspect` → every enabled dataset is `status:
      ready` or says what is `fill_in`.
- [ ] Morris 2022: list
      `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE192nnn/GSE192722/suppl/`,
      put the file URLs under `source.urls`. If the deposit is a Seurat
      `.rds`, it is converted after download (step 4).
- [ ] Huang 2025: same for GSE294534.
- [ ] Isola 2024 needs nothing: it comes through the Census
      (`pip install -e .[census]` inside the env: `bash biowulf/run.sh pip
      install -e .[census]`).

## 4. Fetch (one batch job per dataset; minutes to an hour each)

```bash
bash biowulf/submit.sh fetch --dataset isola2024
bash biowulf/submit.sh fetch --dataset morris2022
bash biowulf/submit.sh fetch --dataset huang2025
bash biowulf/submit.sh python scripts/scbasecount_metadata.py --species mouse
```

- [ ] `data/raw/<name>/download.json` exists for each, with sizes and sha256.
- [ ] A Seurat object: `Rscript R/rds_to_mtx.R data/raw/morris2022/<file>.rds
      data/raw/morris2022/` (module load R first) → `matrix.mtx.gz`,
      `features.tsv.gz`, `barcodes.tsv.gz`, `meta.csv`.
- [ ] `data/scbasecount/mouse_ovary_candidates.csv` lists the ovary
      experiments of scBaseCount (the memo counted 53 by ontology id,
      488,541 cells, 86,346 nuclei in six experiments); the script prints
      one h5ad URL to test a download with `curl -O`.

What good looks like: three raw folders, the candidates table, and
`inspect` showing `on_disk: True` for all three.

## 5. QC, one dataset at a time (minutes each)

```bash
bash biowulf/submit.sh qc --dataset isola2024
```

- [ ] Read the log: the printed obs columns tell you the real `label_key`,
      `animal_key` and condition column; fix them in the manifest and
      rerun until the summary row is right.
- [ ] `tables/qc/qc_summary.csv`: modality agrees with the manifest (the
      tells are printed either way), `n_animals` is the paper's (Isola:
      8 ovaries → check whether `donor_id` is the animal or the library),
      `conditions` are the mapped names (`age_3mo;age_9mo`).
- [ ] `tables/qc/<name>_dropoff.csv`: no step removes more than you can
      explain; nuclei get the 5 % mito cap, cells 20 %.
- [ ] Repeat for `morris2022` and `huang2025`.

## 6. Harmonise and review the crosswalk (an hour, plus the lab's review)

```bash
bash biowulf/submit.sh harmonise
```

- [ ] `tables/harmonise/manifest.csv`: `crosswalk_coverage` per dataset.
      Isola and Morris should be 1.0 (their blocks are written); Huang
      will be low until its clusters are mapped.
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

## 9. Make b1 exist: a second dataset per family (hours)

- [ ] Ageing: enable `wang_suh_estropause` (or `benayoun_vcd` /
      `galligos2026`), fill its source, fetch, qc, harmonise, bench again.
- [ ] Ovulation: enable `ovulation_exvivo_2024` from the scBaseCount
      candidates (fill `source.srx` from the candidates table), same. The
      candidates table also lists a 2025 in vivo series, SRX28395245-50
      ("0hr post-hCG", "hCG treatment", 14-18 K cells each): resolve its
      GEO series and time points and add it as a second ovulation dataset.
- [ ] Rerun step 7 and 8. The ceilings table now has `n_b1_sources ≥ 1`
      for those families, and the b1 row is the bar every model must beat.

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
      then in `sinteractive --gres=gpu:a100:1 --mem=90g --cpus-per-task=16`
      run a few hundred steps on the Isola collection. Record the step time
      and the memory; that is the number the 200-500 A100-hour estimate is
      checked against.
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
