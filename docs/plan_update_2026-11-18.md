# Plan update: from "ovary foundation model" to "virtual ovary v0", intern-ready by May 2027

Nov 18, 2026 · supersedes the timeline in [`project_plan.md`](project_plan.md)
§11 and amends its §3, §5-6 and §9. The research memo
([`ovary_foundation_model_research.md`](ovary_foundation_model_research.md))
and the preregistration draft stand except where this file says otherwise;
every change to the preregistration is logged in its §10.

## 0. Summary

**Where things are.** The repo is a complete, tested scaffold and nothing
more: 11/11 tests pass on the synthetic corpus (re-run 2026-11-18), but no
real dataset has been fetched, no QC run, no ceiling computed, no model
trained, the preregistration is unfrozen, and the phase-0 verification
checklist has not been started. Two weeks after the plan was written the
project is at month 0 of its 12-month timeline.

**What the target is now.** The lab's stated goal is a *virtual ovary*: a
model that, given the cells of any ovarian compartment in any state, emits
that compartment's cells under a named intervention (hormone, age, drug,
chemotherapy, genotype, oophorectomy). Translated into the project's own
evaluation contract that is one function,

    f(control cells of compartment c in context k, condition prompt p) -> integer-count cells,

scored on held-out contexts with the six `vcc2026` members against b1 and
r. The plan already had this function as T3; what changes is that it moves
from "the third layer that a kill rule may remove" to **the product**, with
the corpus and benchmark as its substrate, and the prompt p gains a second
form (section 2.2) that is the only route to "any perturbation".

**What changed in the evidence since Nov 3** (section 2): Biowulf has 400
A100-80 GB and 112 H200-141 GB GPUs, not the 40/80 GB mix the memo could not
read; Arc's Stack (Jan 2026) made *in-context* condition transfer (prompt
cells → counterfactual query cells) a shipped, trainable paradigm; three
2026 cross-context benchmarks confirm that scVI-class representations and
retrieval baselines remain the bar; and the ovary condition data the memo
never searched turn out to exist, including the first **chemotherapy**
snRNA-seq time course and a three-model **menopause** atlas (VCD at two
induction ages × two time points, Foxl2 haploinsufficiency, 20-month
ageing), which together more than double the T3 context list and bring the
"drug" and "chemo" families the lab asked for.

**What is realistic in six months** (section 4): a mouse corpus passing G0;
the ceilings table and baseline scores for ≥ 8 contexts across ≥ 6
condition families; arm C (scVI/scANVI) and a trained in-context condition
model (arm B′) scored against b1; a `ovary-fm simulate` command and a
notebook an intern can run on day one; the ooph 30 d context scored. Not
realistic: validated predictions for arbitrary drugs or for human; the
spatial arm; arm A (scPRINT-2 continued pretraining). Those are the
intern's and the second half-year's work, and the plan says which.

## 1. Status of the repo, verified 2026-11-18

| Item | State |
| --- | --- |
| Package, CLI, 12 stages, Biowulf scripts, synthetic corpus | present; `pytest -q` → 11 passed in 74 s |
| Research memo, project plan, preregistration (draft), first-results checklist | present (Nov 3-4) |
| Real data on disk | none; `config/corpus.yaml` has 3 datasets `enabled`, 9 `fill_in` |
| Phase-0 checklist (`first_results_checklist.md` §0-§2) | not started |
| Ceilings table, baseline scores, T1 table | not produced |
| Models | stub only; scPRINT-2 / UCE / TranscriptFormer branches ported untested; arm C untested on real data |
| Preregistration | draft, `frozen_at` empty |
| Commits since scaffold | none (last: 2026-11-04) |

The scaffold is good; the project is un-started. The six-month plan below
is built around that fact and around the intern deadline.

## 2. Literature and infrastructure update (what the memo did not have)

### 2.1 Compute: Biowulf is much larger than assumed

The NIH HPC hardware page now lists, in the gpu partition: **100 nodes ×
4 A100-80 GB (400 GPUs; NVLink; 32 cores, 256 GB RAM per node)**, **14
nodes × 8 H200-141 GB (112 GPUs; 192 cores, 2.26 TB RAM per node)**, 20
nodes × 8 L40-48 GB (160 GPUs), plus 60 × 4 V100-32 GB and the P100s
([hpc.nih.gov/systems/hardware.html](https://hpc.nih.gov/systems/hardware.html)).
No 40 GB A100s remain. Consequences:

- The memo's extrapolation of 200-500 A100-hours for a 20-90 M model on
  1-5 M cells stands, but on an H200 node (8 × 141 GB) a run that needed
  four A100-days needs roughly one H200-day; the whole 2,500-GPU-hour
  budget of `project_plan.md` §6.2 is a few node-weeks, not four months.
- The binding constraint is therefore **per-user caps and queue time, not
  GPU count**. `batchlim` and `freen` on Biowulf give the current limits;
  they could not be read from outside and are the first thing to record
  (section 5, week 1). `biowulf/submit.sh` requests `gpu:a100:1`; add an
  `SUBMIT_GRES=gpu:h200:N` override and a `--constraint=gpuh200` path.
- A 256 GB A100 node cannot hold a 2 M × 30 K dense corpus; the H200
  nodes (2.26 TB) can, which matters for the `harmonise` and `bench`
  stages as much as for training.

### 2.2 Modelling: in-context condition transfer is now the paradigm

**Stack** (Arc, bioRxiv 2026-01-09; model and code released) is an
encoder-decoder trained on ~150 M human cells with tabular attention over
*sets* of cells; after a post-training alignment it does in-context
learning: given *prompt* cells that exhibit a condition and *query* cells
of interest, it predicts the query cells' state under the prompt's
condition, "to cell types never experimentally perturbed, across
perturbation categories", with no fine-tuning ([preprint](https://www.biorxiv.org/content/10.64898/2026.01.09.698608v1);
[repo](https://github.com/arcinstitute/stack); [HF Stack-Large](https://huggingface.co/arcinstitute/Stack-Large)).
Arc used it to generate "Perturb Sapiens", 28 tissues × 201 perturbations
of simulated cells. Facts that matter here, read from the repo:

- Human only (`human:` dataset prefix; no mouse). Mouse ovary would enter
  through one-to-one orthologs (77-82 % of genes, memo §6) as a query, and
  human prompts (Tahoe drugs, cytokines) could be applied to it; this is a
  zero-shot arm, not a base to train on.
- Code CC-BY-NC-SA 4.0; weights under the Arc Stack Model Non-Commercial
  License + AUP. NIH intramural research is non-commercial, so inference
  and fine-tuning are allowed; the same industry-sponsorship clause as
  STATE applies to any later collaboration, and derived weights inherit
  the licence. This rules Stack out as a *base* for a released ovary model
  and in as a *comparator*.
- Training (`stack-train`) and a "frozen-teacher" fine-tuning
  (`stack-finetune`) are shipped; tested on one H100-80 GB. The
  architecture is the set-based transformer the plan called "the second
  design" (§6.1).
- Arc has said a STATE 2 is in development, building on Stack.

**What this changes.** The plan's T3 head was a decoder conditioned on
(control-set embedding, condition *token*). A token can only name
conditions in the training corpus. The virtual-ovary goal ("any number of
perturbations, drugs, chemo") is reachable only if the condition can also
be given as **cells**: granulosa cells from a chemotherapy study, Tahoe
cells under doxorubicin, hCG-treated cells from another atlas. So:

- The condition head becomes **in-context by design**: the model receives
  (query = this context's control cells of compartment c; prompt = a set
  of cells carrying condition p from any source) and emits perturbed
  query cells. A learned condition token is the degenerate prompt for
  in-corpus families and stays as the first thing trained.
- Training pairs come for free from the corpus: every (dataset, family,
  compartment) with a control and a perturbed group is a prompt/query
  pair; leave-one-dataset-out gives the held-out contexts. Cross-dataset
  pairs (prompt from Huang hCG 4 h, query from Morris 0 h) are the
  training signal for transfer, exactly the signal the retrieval
  baseline b1 uses by hand.
- Stack zero-shot (mouse query via orthologs, human or mouse prompts)
  joins arm Z, so the comparison "our ovary-trained set model vs Arc's
  150 M-cell human generalist with in-context learning" is in every T3
  table. If Stack wins, that is the paper's result and the ovary corpus
  is still the benchmark.

### 2.3 Benchmarks: the 2026 verdict has not moved

- PertDiffBench (bioRxiv 2026-06): across diffusion models for
  perturbation response, "pretrained encoders did not consistently improve
  performance, with the classical scVI representation slightly exceeding
  STATE in seen-condition and unseen-cell-type settings"
  ([preprint](https://www.biorxiv.org/content/10.64898/2026.06.13.732013v1)).
- "Benchmarking virtual cell models for in-the-wild perturbation response"
  (arXiv 2604.27646): under unseen contexts and cross-dataset splits
  "model performance drops markedly".
- A DE-classification benchmark (bioRxiv 2026-07-20) argues evaluation
  should be on DE calls, not mean expression, under both unseen-perturbation
  and unseen-context regimes — which is what four of the six `vcc2026`
  members already do.
- VCC 2026 closed 2026-11-05; winners and the final leaderboard are due
  late November and were not out on 2026-11-18. The plan's risk row
  stands: re-read before the condition head is frozen (section 4, month 3).
- The Tahoe–Arc–Biohub 120 M-cell, 225 K-perturbation dataset (announced
  2026-01-12, to be released open) will be the drug prompt library when it
  lands; like Tahoe-100M it is cancer lines, so it supplies prompts, not
  ovary contexts.

Reading: the benchmark design (whole datasets held out, b1 retrieval as
the bar, scANVI as the embedding bar, DE-based members) is the right one
and does not change. What the literature adds is that the *learned* model
should be a set model with in-context prompts, because that is the one
design with a mechanism for transfer.

### 2.4 Ovary condition data: the searches the memo did not run

| Dataset | Design | Size | Access | New T3 contexts | Family |
| --- | --- | --- | --- | --- | --- |
| **Benayoun lab menopause models** (Kim et al., bioRxiv 2025; PMC12259092) | natural ageing 4 vs 20 mo; **VCD** induced at 3 or 10 mo, profiled 30 and 90 d post-injection, with vehicle; **Foxl2+/−** vs +/+ at ~4 and ~9 mo | 9,559 + 47,386 + 21,301 cells | SRA PRJNA863443 (already in the manifest as `benayoun_vcd`); annotated Seurat objects on Figshare 10.6084/m9.figshare.28934126; code github.com/BenayounLaboratory/Mouse_Menopause_Models | VCD × 2 ages × 2 times (4); Foxl2 × 2 ages (2); ageing 20 mo (1) | chemical_menopause, genotype, ageing |
| **Chemotherapy-induced POI snRNA-seq** (Stem Cell Res Ther 2026, s13287-026-05064-8) | busulfan 30 mg/kg + cyclophosphamide 120 mg/kg i.p. once; sham, 3 d, 21 d post-chemo; n = 3/group, both ovaries pooled per mouse; 8 major types incl. granulosa sub-states and oocyte stages | 70,520 nuclei (+ a 43,420-cell transplantation arm) | **OMIX015442** (CNCB, not GEO): download route to confirm | chemo 3 d, chemo 21 d (2) | **chemotherapy** (new family) |
| DHEA-induced PCOS scRNA-seq + spatial (Int J Womens Health 2026; PMC13695363) | DHEA vs vehicle, mouse ovary, spatial as annotation reference | not read | GEO accession to confirm in the paper | DHEA (1) | drug (androgen excess) |
| IL-11 / matrix stiffness ageing (Nat Aging 2026, s43587-026-01159-2) | snRNA-seq mouse ovarian ageing; Il11ra1 KO | not read | to confirm | ageing nuclei; Il11ra1 KO (1-2) | ageing, genotype |
| FOXO3 POI spatial atlas (2026; PMC13355553) | spatial transcriptomics of POI models | — | to confirm | spatial arm only | — |
| Oocyte development atlas (bioRxiv 2026-03-11) | Stereo-cell oocytes + DNBelab somatic cells, growth to maturation | not read | to confirm | oocyte-stage contexts if counts allow | development |
| 10x Visium HD mouse ovary demo (published 2026-07-14, Space Ranger 4.1) | 8-week C57BL/6, fresh frozen | 1 section | 10x datasets page | spatial arm comparator (same platform as ours) | — |
| McKey lab Visium HD developing/adult ovary (bioRxiv 2025-09-16) | E12-P90 | — | GSE307688 (already in plan §4) | spatial arm | — |
| Cisplatin-POI mechanism paper (Cell Signal 2026) | mouse + KGN; bulk | — | — | none (bulk) | — |

With the four contexts already in `benchmark.yaml`, the list is now **≥ 14
T3 contexts in 8 families** (ovulation, cycle, ageing, chemical menopause,
genotype, chemotherapy, drug, surgery). Gate G1 (≥ 5 contexts with a
computable r) is no longer in doubt on paper; what remains uncertain is
how many have ≥ 2 animals per group *and* ≥ 50 cells per compartment, which
the ceilings table decides. The chemotherapy dataset is nuclei, like the
own data, which makes it the natural modality twin of the ooph context.
Cross-check the Benayoun counts against the figure legends and the
Figshare objects on download; the OMIX route for the chemo data is the one
access risk and is tried in week 1.

### 2.5 Tooling

`cell-eval2` is current on PyPI with the `vcc2026` preset
(`EvalConfig.from_preset("vcc2026")`), GPU extras and out-of-core
`cellstream` scoring; `vcc-cli` 0.2.2 is MIT. The `--scorer cell-eval2`
path in `ovary-fm score` is still untested and is a week-2 item.

## 3. The virtual ovary, translated into deliverables

| Phrase in the goal | What it means in this project | By May 2027 | Later |
| --- | --- | --- | --- |
| "simulate any cell in the ovary" | emit integer-count cells for a compartment given its controls | the 7 somatic compartments where a context has ≥ 400 cells; oocytes and luteal cells scored where the gate allows, else descriptive | fine-level types once the crosswalk covers them |
| "predict any number of perturbations" | the condition prompt p | **token prompts** for the 8 in-corpus families; **cell prompts** evaluated for cross-dataset transfer within a family (e.g. Huang hCG → Morris) | cross-tissue and cross-species cell prompts (Tahoe drugs → ovary) once the Tahoe–Arc–Biohub release exists |
| "drugs" | a chemical the ovary has not been profiled under | VCD and DHEA token prompts; doxorubicin ± AMH, metformin, LPS from scBaseCount as contexts if their SRA metadata holds up | the general drug prompt library, unvalidated for ovary until a wet-lab test exists |
| "chemo" | chemotherapy | busulfan + cyclophosphamide 3 d and 21 d as scored contexts (the first in existence) | cisplatin, doxorubicin as they appear |
| "oophorectomy" | the own data | ooph 30 d scored as a held-out context with its (low, 2 vs 2) ceiling | 120 d when it lands, designed to be scorable (plan §9) |
| "virtual" | a thing someone else can run | `ovary-fm simulate` + a notebook + a frozen checkpoint under `/data/$USER` | a Hugging Face release past KR2 |

Honesty clause for the write-up: a prediction for a prompt whose family
has no held-out ovary context is a *hypothesis list*, labelled as such in
the output (`prediction_status: unvalidated_family`), never a result.

## 4. The six-month plan (2026-11-18 → 2027-05-18)

One computational lead; the lab for the crosswalk review and gene lists;
Biowulf gpu partition from month 3. Months are four-week blocks from today.

| Month | Weeks | Work | Gate / output |
| --- | --- | --- | --- |
| **1** (Nov 18 - Dec 15) | 1 | `batchlim`, `freen`, record GPU caps; checklist §0-§2 (verify the snippet numbers that drive the design); fetch isola, morris, huang; try the OMIX download; download the Benayoun Figshare objects; `cell-eval2` on the synthetic object | caps recorded; 3 datasets raw on disk |
| | 2-4 | `qc` and `harmonise` for the first 3; crosswalk v0 to the lab; manifest entries for the 2.4 datasets filled from their downloads; `fm --tool scprint2` single-GPU check | **first ceilings table** (3 contexts); crosswalk_review.csv out for review |
| **2** (Dec 15 - Jan 12) | 5-8 | fetch + QC: Benayoun (3 models), chemo-POI, Wang/Suh, ex-vivo ovulation series, Gaylord mouse, Garcia-Alonso; scBaseCount mouse ovary ingest (the 53 UBERON experiments) with hand-curated condition metadata; b1 sources exist for every family with ≥ 2 datasets | **G0** checked (≥ 0.75 M mouse obs; crosswalk ≥ 90 %); ceilings for ≥ 8 contexts; **G1** decided |
| **3** (Jan 12 - Feb 9) | 9-12 | baselines scored (b0, b1, ridge, scVI counterfactual); arm C scVI/scANVI on the corpus (GPU, < 20 h); T1 table with scANVI, logreg, CellTypist, and the three generalists; **freeze the preregistration** (with the 2.4 contexts and the in-context head written in; re-read the VCC 2026 results first) | baseline and T1 tables; `frozen_at` filled |
| **4** (Feb 9 - Mar 9) | 13-16 | **arm B′**: the in-context set model (≤ 30 M params; encoder from arm B's scPRINT-2-small / Tahoe-x1-70M configuration or a from-scratch count encoder; set-attention transition; NB count decoder with control-matched dispersion), token prompts first, trained leave-one-dataset-out on A100/H200; Stack zero-shot via orthologs as arm Z | first T3 scores against b1, 3 seeds |
| **5** (Mar 9 - Apr 6) | 17-20 | cell prompts (cross-dataset within family); data-size ablation 50 K / 200 K / all (arm B′, 3 seeds); species ablation as arm C only (cheap); ooph 30 d scored | **KR1 / KR2** read; the paper-1 tables (corpus, crosswalk, ceilings, baselines, scaling) complete |
| **6** (Apr 6 - May 18) | 21-26 | `ovary-fm simulate` (section 4.1); frozen checkpoint + model card; the intern notebook and the intern project menu; write-up of paper 1 drafted; buffer | **intern-ready**: a clean clone + one script reproduces every table; `simulate` runs on a laptop against the Biowulf checkpoint |

Dropped from the first six months, with the reason:

- **Arm A** (scPRINT-2 medium continued pretraining): GPL-3 on derived
  code, weights licence unstated, and the scaling literature says a 90 M
  model buys little at 1-2 M cells. It becomes an optional intern project
  if the licence question is answered.
- **The spatial arm** (phase 4): no spatial model beats naive baselines on
  local structure (memo), and the intern deadline leaves no month for it.
  The 10x Visium HD demo section and GSE307688 are catalogued for it.
- **Human as a second species in pretraining**: kept only as the cheap arm
  C ablation (KR4); the human corpus (Tier C) is not curated before May.
- **T2 and T4** stay out of the scaffold.

Compute, revised: arm C < 20 GPU-h; arm B′ ≈ 100-300 A100-h per run × (1
main + 3 seeds + 3 data-size points × 3 seeds + leave-one-out folds) ≈
1,500-2,500 A100-h ≈ 400-700 H200-h. At a cap of 8 GPUs per user that is
2-4 node-weeks of wall time in months 4-5; at 4 GPUs it still fits. The
memory-bound stages (`harmonise`, `bench`) go to an H200 node or a largemem
node, not the 256 GB A100 nodes.

### 4.1 What the intern receives

- `ovary-fm simulate --context <dataset>:<compartment> --condition
  <name> [--prompt cells.h5ad] --n 400 --out pred.h5ad`, which loads the
  frozen checkpoint, builds the query from the context's control cells and
  the prompt from the condition token or the supplied cells, samples
  integer counts, and writes `prediction_status` and provenance.
- `ovary-fm score --pred pred.h5ad --name <run>` against the bundles, so a
  new idea is one command from a number on the same scale as every
  baseline.
- A notebook: one context end to end (controls → simulate → score → the
  DE table read against the lab's gene sets).
- The project menu, each item scorable with the existing runner: (a) add a
  context (the 120 d ooph cohort when it lands; a scBaseCount SRA
  experiment); (b) a prompt-library experiment (Tahoe doxorubicin cells as
  a prompt vs the measured doxorubicin ± AMH context); (c) the modality
  ablation (nuclei-only vs mixed) that is still pre-registered and unrun;
  (d) Stack fine-tuning with the frozen-teacher recipe on the ovary corpus
  as a licence-permitting comparison; (e) arm A if GPL is acceptable.

## 5. Changes made to the repo with this update

- This file; `README.md` points to it.
- `config/corpus.yaml`: `benayoun_vcd` replaced by `benayoun_menopause`
  (the three models, the Figshare objects, the SRA project); new entries
  `chemo_poi_2026` (OMIX015442), `dhea_pcos_2026`, `il11_ageing_2026`; all
  `enabled: false`, `status: fill_in`, counts marked † until downloaded.
- `config/benchmark.yaml`: families `chemotherapy` and `drug` populated;
  contexts `benayoun_vcd_3mo_30d`, `benayoun_vcd_3mo_90d`,
  `benayoun_vcd_10mo_30d`, `benayoun_vcd_10mo_90d`, `benayoun_foxl2_4mo`,
  `benayoun_foxl2_9mo`, `benayoun_ageing_20mo`, `chemo_poi_3d_21d`,
  `dhea_pcos`; the `Foxl2` exclusion for the genotype contexts.
- `docs/preregistration.md` §10: one deviations row (the file is a draft,
  so this is a recorded edit, not a deviation from a frozen rule): the
  contexts above, the in-context head as the first T3 design, Stack as an
  arm Z member, arm A and the spatial arm deferred.
- `biowulf/submit.sh`: `SUBMIT_GRES` override so a stage can ask for
  `gpu:h200:N`.

Nothing here touches the metrics, margins or kill rules.

## 6. Risks specific to the six-month plan

| Risk | Mitigation |
| --- | --- |
| OMIX (CNCB) download fails or is slow from Biowulf | try week 1; fall back to the authors; the chemotherapy family is then descriptive until a GEO set appears |
| Benayoun groups have < 2 animals per cell × condition after QC | the ceilings table says so in month 2; contexts without r are listed, not scored |
| GPU caps per user are low (e.g. 4) | arm B′ is ≤ 30 M params and fits one GPU; serialise seeds; drop data-size points before seeds |
| In-context head underperforms b1 everywhere | KR2 says it is a negative result; the intern deliverable is unchanged (simulate then emits b1's prediction with `prediction_status: baseline_b1`) |
| VCC 2026 results (late Nov) show something that changes the head design | re-read before the month-3 freeze; the benchmark does not change |
| Nick's time is split across projects | the month-1 and month-2 work is pipeline running on Biowulf, not modelling; each month has one output that stands alone |

## 7. Sources consulted for this update (beyond the memo's register)

- NIH HPC hardware: https://hpc.nih.gov/systems/hardware.html ; user guide https://hpc.nih.gov/docs/userguide.html
- Stack: https://www.biorxiv.org/content/10.64898/2026.01.09.698608v1 ; https://github.com/arcinstitute/stack ; https://arcinstitute.org/news/foundation-model-stack
- PertDiffBench: https://www.biorxiv.org/content/10.64898/2026.06.13.732013v1
- In-the-wild VC benchmark: https://arxiv.org/abs/2604.27646
- DE-classification benchmark: https://www.biorxiv.org/content/10.64898/2026.07.20.739620v1
- VCC 2026 status: https://virtualcellchallenge.org/ ; https://arcinstitute.org/news/virtual-cell-challenge-2026
- Tahoe–Arc–Biohub dataset: https://biohub.org/news/arc-tahoe-perturbation-dataset/
- Benayoun menopause models: https://pmc.ncbi.nlm.nih.gov/articles/PMC12259092/
- Chemotherapy POI snRNA-seq: https://link.springer.com/article/10.1186/s13287-026-05064-8
- DHEA PCOS scRNA-seq: https://pmc.ncbi.nlm.nih.gov/articles/PMC13695363
- IL-11 ovarian ageing: https://www.nature.com/articles/s43587-026-01159-2
- 10x Visium HD mouse ovary: https://www.10xgenomics.com/datasets/visium-hd-three-prime-mouse-ovary-fresh-frozen
- cell-eval2: https://pypi.org/project/cell-eval2/
