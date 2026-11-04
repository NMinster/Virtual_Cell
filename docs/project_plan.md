# Project plan: an ovary-specific foundation model

Nov 3, 2026 · project plan, drafted from
[`ovary_foundation_model_research.md`](ovary_foundation_model_research.md)
(the research memo) and the `raghu-update` branch of
[Mouse_Compensatory](https://github.com/NMinster/Mouse_Compensatory).

The Virtual Cell Challenge 2026 is the motivation and the template, not a
target: nothing here is submitted to it. What it gives this project is a
task shape (predict a cellular state in a context the model never saw,
given only that context's controls) and an evaluation contract (six
metrics, each rescaled so that 0 is a trivial baseline and 1 is a
biological replicate) that a small lab can run on its own data and that
reviewers sceptical of foundation models will accept.

## Summary

**The decision.** Build a narrow, well-evaluated ovary model in three
layers that are each worth publishing on their own, and let a
pre-registered kill rule decide whether the third layer survives:

1. **The corpus**: the >1 million mouse ovary cells, nuclei and spots and
   the >0.5 million human ovary nuclei that exist in GEO and SRA but not
   in CELLxGENE, re-processed from raw counts into one object with a
   Cell Ontology-aligned label crosswalk and condition metadata (estrous
   stage, age, hormone time point, genotype, surgery).
2. **The benchmark**: held-out-context tasks scored with Arc's
   `cell-eval2` contract (`vcc2026` preset: 0 = context-mean baseline,
   1 = split-half replicate, raw counts, condition-defining genes
   excluded, no clipping), plus leave-one-atlas-out label transfer with
   rare-type F1, each with its measured data ceiling.
3. **The model**: a mouse-native model of at most about 100 M parameters,
   either scPRINT-2 medium continued-pretrained on the corpus or a
   ≤30 M-parameter model trained from scratch in the scPRINT-2 or
   Tahoe-x1 configuration, judged against scVI/scANVI, linear and
   retrieval-and-transfer baselines, with the existing generalist pilot
   (scPRINT-2, UCE, TranscriptFormer) as the zero-shot arm.

**Why this shape.** The research memo measured the gap rather than
asserting it: human ovary is 0.375 % of the Census (rank 31 of 71) and
half of it is tumour; mouse ovary is 0.083 % (rank 27 of 36), two
studies, zero nuclei, zero spatial, nothing older than 9 months, 237
oocytes; Tahoe-100M has no ovarian line; no reproductive-tissue
foundation model exists in any form. The same memo found that the field's
2026 benchmarks put scVI, PCA and linear baselines at or above generalist
models on embedding, integration and trajectory tasks, that pretraining
scale saturates early, and that the only places foundation models win
are perturbation prediction with in-context data, low-supervision
annotation and gene-level tasks. So the honest differentiator of an
ovary model is not "an ovary embedding" (the September memo's verdict on
that stands) but **condition and perturbation coverage the generalists
never saw**: the cycling, ovulating, ageing, knocked-out and surgically
altered ovary, evaluated on held-out contexts.

**What is not in scope.** Entering the Challenge; training a large model
from scratch (the only organ-specific precedent, a kidney model, had
39.5 M cells, twenty times what ovary has); claiming the model improves
the ooph vs sham analysis itself (that dataset is a held-out test
context, not training data); any in vivo perturbation experiment (an
option for the lab, raised in section 9, not a dependency).

## 1. Goals, ranked, and what each is worth alone

| # | Goal | Deliverable | Worth publishing if the model fails? |
| --- | --- | --- | --- |
| G1 | A curated, count-level, ontology-aligned ovary corpus | the object, its manifest, the label crosswalk, a QC report | yes: a resource note; nothing like it exists |
| G2 | A held-out-context benchmark with measured ceilings | the task definitions, the reference bundles (baseline b and replicate r per context), the baseline scores, the ceilings | yes: no ovary benchmark exists, and the data ceilings alone are informative |
| G3 | A model that beats the baselines on G2 by a pre-registered margin | weights, code, the comparison tables | only as a negative result, which the field is currently publishing |
| G4 | Two ablations nobody has run for any organ | the organ-specific data-size curve; mouse-only vs mouse + human at matched compute | yes, either way |
| G5 | The Mouse_Compensatory use | ooph vs sham scored as a held-out context; the 120 d cohort designed so its replicate ceiling is computable | part of the ooph paper |

## 2. The template: what the Challenge contributes and how it translates

The 2026 Challenge gives participants, per context, 18,400 non-targeting
control cells and a list of 300 target genes, and asks for 400
integer-count cells per perturbation; nothing is trained on the
challenge contexts. Scoring is six metrics, each rescaled as
s = (u − b) / (r − b): b is the context's mean-perturbation-response
baseline, r is the real data scored against its own split half (five
seeds), a context's score is the equal-weight mean of the six, the
overall score the mean over contexts, and target genes are excluded from
every distance. The 2025 profile (three metrics clipped at zero) is not
reused because Arc's own post-mortem found it gameable.

| Challenge element | Ovary translation |
| --- | --- |
| context = an anonymised cancer cell line | context = one ovarian compartment in one dataset's baseline state (granulosa, theca, stroma, luteal, endothelial, immune, epithelial; oocytes when numerous enough) |
| perturbation = a CRISPRi knockdown | condition with a measured response: hours after hCG, aged vs young, VCD vs vehicle, knockout vs wild type, drug vs vehicle, ooph vs sham |
| control cells supplied as input | the context's baseline cells (0 h, young, wild type, vehicle, sham) |
| target gene list; targets excluded from scoring | the condition's defining genes are excluded: immediate-early genes for hCG time points, Xist and strain markers, the knocked-out gene and its neighbours |
| 400 cells per perturbation, 18,400 controls | ≥ 400 cells per (context, condition) wherever possible; the ≥ 10-significant-gene gate on the fold-change metric decides which rare compartments are scorable |
| b = context mean response | the scale is anchored on b0 = the control copied (no effect); b1 = the same condition family's measured response transferred from other datasets (the "retrieval and transfer" recipe that placed in both Challenges) is scored like a model and its row is the bar to beat. Bundles whose compartment did not respond (fewer than 5 real DE genes after exclusions) are listed with their ceilings but not scored |
| r = split-half replicate of the real perturbed cells | the same, with halves drawn **across animals** whenever a condition has ≥ 2 animals, so r reflects biological not technical replication |
| held out: whole cell lines | held out: whole datasets (leave-one-atlas-out) and, for the own data, the whole experiment |
| raw integer counts in, pseudobulk and DE scored | the same; a count-likelihood model with a sampling step is required even for a pseudobulk-style predictor |
| perturbation discrimination among 300 targets | discrimination among all conditions pooled within a compartment across datasets (a single dataset has too few conditions for the rank to mean much) |

Four of the six members are differential-expression based and their
baselines sit near zero, so recovering the identity and sign of the
responding genes dominates the score; fold-change magnitude and
integer-count error are the hard terms. That is the right emphasis for
ovary biology, where the questions are which granulosa genes move at 4 h
after hCG, or which stromal genes move with age, not the exact UMI count.

## 3. Phase 0: decisions and pre-registration (weeks 1-3)

Everything below is written down before any model is trained, in
`docs/preregistration.md`, and the file is frozen by commit hash.

**Scope decisions, made here unless the lab overrides them.**

| Decision | Choice | Reason |
| --- | --- | --- |
| Species | mouse first; human as the second species in a pre-registered ablation, never silently mixed in | the own data and most condition-structured public data are mouse; cross-species benefit is unmeasured for any organ |
| Modality | single-cell and single-nucleus together, with a suspension-type token; spatial bins kept out of core pretraining (section 7) | the field encodes modality as a token rather than filtering genes; no model includes Visium-type spots in pretraining except scGPT-spatial, with a dedicated decoder |
| Gene space | Ensembl ids; protein-embedding gene tokens (ESM2-style) rather than a one-to-one ortholog table | one-to-one orthologs keep only 77 % of mouse genes (82 % with rescue); embedding tokens let ovary-specific genes in without retraining a table |
| Objective | count likelihood on raw counts | avoids binning or ranking across the depth differences of nuclei, cells and bins; it is also what the scorer consumes |
| Base model, arm A | scPRINT-2 medium (90 M), continued pretraining | mouse-native, count-based, ESM-style gene embeddings, ships the pretraining recipe and SLURM configs, already saw the public ovary data; GPL-3 on derived code, weights licence unstated (ask the authors) |
| Base model, arm B | ≤ 30 M parameters from scratch in the scPRINT-2 small or Tahoe-x1 70 M configuration | the scaling literature says capacity buys little at 1-2 M cells; Apache-2 (Tahoe-x1) or GPL-3 (scPRINT-2) |
| Primary baseline, arm C | scVI/scANVI on the same corpus, plus a Richter-style masked-autoencoder MLP | the 2026 benchmarks say this is the model to beat |
| Excluded bases | Arc STATE weights (non-commercial licence that explicitly covers fine-tunings and excludes industry-sponsored research), scFoundation weights (same class), scGPT and Geneformer checkpoints (human-only vocabularies), UCE and TranscriptFormer (inference-only code) | licence or vocabulary; all three generalists stay as the zero-shot arm through the existing `fm` stage |
| Compute | NIH Biowulf gpu partition, one 4 × A100 node per run | the maintainers' scPRINT-2 recipe is 4 × A40-40 GB with 16 CPU and 90 GB RAM per GPU |

**Pre-registered metrics, baselines, margins and kill rules.** Section 8
has the table. The rule of the September memo (`fm.compare` in
Mouse_Compensatory: keep an embedding only if it mixes samples at least
as well as Harmony, reaches ARI ≥ 0.8 and transfers labels no worse)
carries over as the embedding-level rule.

**Verification checklist.** The research memo tags every number it could
not open at the source as a search-snippet claim; before any of them
goes into a grant or a manuscript they are re-read in the paper. The
ones the plan actually leans on:

- cell and nucleus counts for Morris 2022 (the repo's reused object has
  28,014 cells; the paper reports >34,000; GEO GSE192722), Wagner 2020,
  Huang 2025 (GSE294534), Wang/Suh (GSE267729), Gaylord 2025 (GSE290743),
  GSE329417 (439,011 nuclei), the Slide-seq ageing atlas (610,620 spots);
- the STATE margins (63 % / 47 % / 5 %) and the cancer
  continued-pretraining delta (0.86 vs 0.81);
- Geneformer's licence (model-card metadata only) and scPRINT-2's weights
  licence (none stated);
- Biowulf's current A100 40 vs 80 GB mix and walltime limits.

Searches the memo could not run before its budget ran out, to be run in
week 1: PCOS, POI and chemotherapy mouse scRNA-seq; KGN / COV434 / HGrC1
pooled screens; FSH, androgen and AMH single-cell treatments; a
PubMed and bioRxiv query for reproductive-tissue foundation-model
abstracts from 2025-2026 conferences; the organ lists of Nicheformer and
Tahoe-x1; the scPerturBench and PertEval numbers.

Two things to try on a single GPU in week 2, before anything is
committed to: `scprint2` pretraining for a few hundred steps on the
Isola 2024 object (does the recipe run on Biowulf at all, and at what
memory), and `cell-eval2 run --preset vcc2026` on a synthetic
two-condition object (does the reference-bundle machinery accept
cross-animal split halves).

## 4. Phase 1: the corpus (months 1-4)

### 4.1 Sources, tiered by provenance

Counts in brackets are from the research memo; those marked † are
snippet-derived and get confirmed on download. Everything enters as raw
counts; nothing enters from an integrated or normalised object.

**Tier A, mouse, curated atlases and condition studies**

| Dataset | What it contributes | Accession |
| --- | --- | --- |
| Isola/Ocañas 2024 (Stout lab) | 3 vs 9 months, 14,349 cells, the authors' 16 labels; already the repo's reference atlas | GSE232309 (also Census) |
| Morris/Meinsohn 2022 (Pépin lab) | the only estrous-staged whole-ovary atlas; PMSG/hCG 8 h arm; lactation; the Level 0/1 vocabulary the Cell Ontology terms were anchored to | GSE192722 (28,014-cell reused object; >34,000† in the paper) |
| Huang 2025 (Shalek lab) | hCG 0 / 4 / 12 h, n = 3 per time point, 25,714† cells, paired MERFISH | GSE294534 |
| ex vivo ovulation series (Nature Cell Biology 2024) | follicles at 0 / 3 / 6 / 9 / 12 h after hCG, with in vitro arms | scBaseCount SRX23566711-38 (confirm the GEO series) |
| Wang/Suh estropausal transition (PNAS 2026) | young, peri- and post-estropausal, estrous-staged, senescence-sorted | GSE267729† |
| Gaylord 2025 (Laird lab, Science) | mouse and human across age, ~100,000† cells, 11 types incl. glia | GSE290743 (human object linked; confirm the mouse part) |
| Benayoun cross-species ageing inputs | natural ageing cohorts and the VCD chemical-menopause model | PRJNA863443, E-MTAB-12889, GSE290742 |
| Galligos/Gerton 2026 inflammaging | ageing with immune expansion; spatial twin | GSE317144 (+ GSE317297 spatial) |
| lifespan and puberty sets | P1-12 months stroma; 2-6 weeks weekly; ovary vs oviduct vs uterus | GSE127106, GSE268466, GSE241318 |
| developmental | pregranulosa 52,542† cells (E11.5-P5); Garcia-Alonso mouse gonad 63,929 | Niu & Spradling 2020 (accession to confirm); Census "Mouse Female" |
| genetic models | ERα knockout; Cacna2d3 knockout (16,705† cells); Trim28 knockout | PLOS ONE 2024 (accession to confirm); PMC10626498; SRA |
| **own data** | ooph vs sham 30 d, 4 mice, ~50,000 nuclei, Visium HD twins; 120 d cohort to come | Mouse_Compensatory delivery (held out, never trained on) |

**Tier B, mouse, SRA-wide via scBaseCount** (release 2026-01-12): 53
UBERON-ovary experiments with 488,541 cells, 870,640 cells by text match,
86,346 nuclei in 6 experiments; includes superovulation, metformin
(SRX14708250), doxorubicin ± AMH, melatonin + DEHP, LPS, high-fat diet,
Trim28 knockout, rete ovarii and capsule, and ID8 lines (excluded). These
counts come from one uniform STARsolo pipeline with intronic reads
included, which is the right setting for mixing with nuclei. Each
experiment's condition metadata is curated by hand from its SRA record.

**Tier C, human** (second species, ablation only): Census normal ovary
175,002 cells (Fan 2019, Lengyel 2022, Jones 2024, the 42,568-nucleus Jin
atlas GSE202601, Tabula Sapiens), the 439,011-nucleus† postmenopausal
atlas GSE329417, GSE267315, Wu/Tang 2024 (92,965† cells), Wagner 2020
(>24,000†), the Garcia-Alonso fetal gonad (213,898), and scBaseCount's
239,418 non-cancer cells. Ovarian cancer (half of Census ovary; 86 % of
scBaseCount human ovary) is excluded from the core and kept as an
out-of-distribution test set.

**Tier D, other species** (optional, cross-species ablation): macaque
71,721, pig 34,758, sheep 26,055 via scBaseCount; the goat ageing set
PRJNA1010653.

**Spatial** (the spatial arm, section 7, not core pretraining): the McKey
lab Visium HD atlas GSE307688 (same platform as ours, E12 to P90),
Mantri 2024 Slide-seq GSE240271 (eight ovulation stages), the Shalek/Chung
Slide-seq cycling and ageing atlas (610,620† spots), the 10x Visium HD
demo section, our own Visium HD captures.

Target size: 0.9-1.2 M mouse observations, 0.5-1.0 M human; 1-2 M in all.

### 4.2 Processing

One pipeline, parameterised in `config/corpus.yaml`, with the rules
Mouse_Compensatory already has: modality detection (mito, Malat1,
ribosomal shares, intronic fraction), per-sample MAD bounds, a 5 % mito
cap for nuclei, kit-rate doublet calls, ambient-gene flags. Every dataset
gets a provenance record (accession, file, md5, pipeline version,
`git rev-parse HEAD`).

Metadata schema per cell: `dataset`, `species`, `assay`, `suspension_type`,
`strain`, `sex`, `age_days`, `estrous_stage` (when staged), `condition`,
`condition_time_h`, `genotype`, `treatment`, `animal`, `author_label`,
`cl_id`, `cl_label`, `crosswalk_status`.

### 4.3 The label crosswalk

Nearly every ovary atlas predates the fourteen-plus ovary terms the Cell
Ontology gained between April 2024 and June 2026 (pre-granulosa
CL:4033066, mural granulosa CL:4033067, cumulus, squamous and cuboidal
granulosa, atretic granulosa and theca, pre-theca, theca interna and
externa, steroidogenic stromal cell, early luteal CL:4052046 and late
luteal CL:4052047), and no crosswalk from atlas labels to them has been
published. Building it is a prerequisite for the benchmark and a
contribution in itself.

- Seed: `config/reference_labels.yaml` from Mouse_Compensatory already
  maps Isola 2024 and Morris 2022 Level 0/1 onto the repo's eight types;
  extend it to CL ids at the finest level each atlas supports.
- Rule: a label maps to the most specific CL term its authors' markers
  justify; a label with no term is kept as author free text with
  `crosswalk_status: unmapped` and listed, never folded silently (the
  same rule the repo's agreement tables use).
- Review: the crosswalk table goes to Raghu's lab with the marker
  evidence per row; the reviewed table is the one released.
- Check: train the repo's CellTypist/logistic reference on each atlas
  under the crosswalk and confirm the cross-atlas confusion is
  concentrated where the biology says it should be (granulosa stages,
  theca vs steroidogenic stroma, early vs late luteal).

### 4.4 Deliverables and gate G0

Deliverables: the corpus object (h5ad shards or a TileDB-SOMA store under
`/data/$USER`), `config/corpus.yaml` (the manifest), the crosswalk
(`config/crosswalk.yaml` and the reviewed CSV), a QC report, and a short
resource note.

Gate G0: at least 0.75 M curated mouse observations with raw counts and
condition metadata, and crosswalk coverage of at least 90 % of cells at
the eight-type level. Below that, publish the inventory and stop.

## 5. Phase 2: the benchmark (months 2-5)

### 5.1 Tasks

| Task | Question | Held out | Score | Baselines |
| --- | --- | --- | --- | --- |
| T1 label transfer | does the representation carry ovary identity across datasets? | one whole atlas at a time (Isola, Morris, Huang, Wang/Suh); one species direction | macro-F1 and rare-type F1 (oocyte, luteal, theca subtypes, atretic granulosa) under the crosswalk | scANVI, scPoli, the repo's CellTypist model, kNN on PCA/Harmony, the three generalists zero-shot |
| T2 integration, with caveats | does the representation mix datasets without erasing condition? | — | scIB-style metrics stratified by suspension type; the repo's observed/expected same-sample mixing table | Harmony, scVI; silhouette-based metrics reported but not trusted alone |
| T3 condition effect in a held-out context | given a compartment's baseline cells only, predict its perturbed state | whole datasets | the six `vcc2026` members, s = (u − b1) / (r − b1), per context and averaged | b0 control copy; b1 retrieval-and-transfer; ridge/linear on the control pseudobulk (Ahlmann-Eltze style); scVI counterfactual decoding |
| T4 biological priors | do gene embeddings and programs recover ovary biology? | — | recovery of the lab's gene sets, the follicle-stage hallmarks and the Morris trajectories in gene-embedding neighbourhoods; a leave-genes-out imputation test | cNMF programs; co-expression on the same corpus |

T3 contexts, from the condition-structured datasets in section 4:

| Context family | Control | Conditions | Animals per condition |
| --- | --- | --- | --- |
| ovulation (granulosa, theca, cumulus) | hCG 0 h | 3, 4, 6, 9, 12 h (Huang 2025; the ex vivo series; Morris 8 h) | 3 (Huang); confirm others |
| ageing (every compartment) | 3 months | 9 months (Isola); peri- and post-estropausal (Wang/Suh); Benayoun cohorts | 4 (Isola) |
| chemical menopause | vehicle | VCD (PRJNA863443 subset) | confirm |
| genotype | wild type | ERα KO, Cacna2d3 KO, Trim28 KO | confirm |
| drug and diet | vehicle | metformin, doxorubicin ± AMH, LPS, high-fat diet (SRA) | confirm |
| **surgery (own)** | sham | ooph 30 d; ooph 120 d when it lands | 2 (30 d) |

Condition-defining genes excluded from scoring per family: the
immediate-early and known LH-target genes for ovulation time points
(Egr1, Fos, Jun, Ptgs2, Areg, Ereg, Pgr, Runx1; the list is finalised
from the literature before scoring and written into the preregistration);
Xist and strain-polymorphic genes everywhere; the knocked-out gene and
its immediate neighbours for genotype contexts.

### 5.2 Ceilings and gates

For every context: compute r with five cross-animal split halves, report
the seed spread, list the compartments that fail the ≥ 10-significant-gene
gate (oocytes and luteal cells will, often) and report them descriptively
rather than scored. The ceilings table is itself a result: it says how
much of each ovarian response is reproducible between animals at all,
which no one has measured.

Gate G1: at least five T3 contexts with a computable r and at least two
compartments each. Fewer, and T3 shrinks to a descriptive comparison
while T1, T2 and T4 carry the benchmark.

### 5.3 Deliverables

The task definitions and splits (`config/benchmark.yaml`), the reference
bundles (b0, b1, r per context), baseline scores, the ceilings table, and
a runner (`ovary-fm bench`) that scores any model emitting the standard
prediction object. The runner wraps `cell-eval2` rather than
reimplementing it.

## 6. Phase 3: the model (months 4-9)

### 6.1 Arms

| Arm | What | Data | Compute (A100-hours, extrapolated) |
| --- | --- | --- | --- |
| A | scPRINT-2 medium, continued pretraining, a few epochs | mouse corpus; then mouse + human | 200-500 per run |
| B | ≤ 30 M from scratch, scPRINT-2 small or Tahoe-x1 70 M config | same | 100-300 per run |
| C | scVI/scANVI (+ masked-autoencoder MLP) | same | < 20 |
| Z | zero-shot generalists through the existing `fm` stage (scPRINT-2, UCE, TranscriptFormer) | — | inference only |

For T3 the model needs a conditional, count-level generator: predict the
perturbed cells of a context from its control cells and a condition
token. The simplest honest version is a decoder head on the pretrained
encoder conditioned on (control-set embedding, condition embedding), with
control-matched dispersion in the sampling step; a set-attention
transition module in the STATE style is the second version if the first
beats b1. Either way the sampling step is required by the scorer and is
where the integer-count member is won or lost.

### 6.2 Ablations (pre-registered, three seeds each)

- **Data size**: 50 K / 200 K / 1 M / all cells, arm B. The organ-specific
  scaling curve nobody has published; if it is flat past 200 K, the
  corpus stops growing and the paper says so.
- **Species**: mouse-only vs mouse + human at matched compute, arms A and
  B, scored on mouse tasks only. The only direct evidence that human
  helps mouse is one kidney number; this settles it for ovary.
- **Modality**: nuclei-only vs cells + nuclei with the suspension token,
  scored on nuclei tasks (the own data are nuclei).
- **Base**: arm A vs arm B, same data, same compute budget.

Budget: roughly 2,500 A100-hours in all (12 arm-B runs for the size
curve, four arm-A runs, the species and modality runs, inference and the
baselines). On one 4 × A100 node that is about four months of wall time
if runs are serialised, or weeks if Biowulf allows two nodes at once.

### 6.3 Engineering

- One conda stack per model, the Mouse_Compensatory "fourth stack"
  pattern (`environment-fm-*.yml`): the training env never shares torch
  pins with the analysis env. scPRINT-2 pins torch 2.8; Tahoe-x1 needs
  CUDA 12.1+ and FlashAttention.
- Data loading: scDataLoader (lamindb-backed, what scPRINT-2 trains
  from) for arm A; the Tahoe-x1 Composer pipeline for arm B; one
  exporter from the corpus object to each.
- SLURM: `biowulf/submit.sh` extended with a `train` stage on the gpu
  partition; checkpoints and `HF_HOME` under `/data/$USER`, not the home
  quota; proxy exported for weight downloads.
- Provenance: the repo's `write_provenance` pattern (git revision, live
  config, checkpoint hash, genes matched, wall time) on every run.
- Tests: the stub model exercises the plumbing on the synthetic
  two-condition object; real-model tests are skipped unless the
  checkpoint is present (as in `test_smoke.py`).

### 6.4 Gate G2 (the kill rule, in full in section 8)

The model survives only if it beats scANVI on T1 rare-type F1 and beats
b1 on T3 on at least two held-out contexts by the pre-specified margins.
Otherwise the corpus, crosswalk, benchmark, ceilings and ablations are
the paper, and the model is reported as the negative result it is.

## 7. Phase 4: the spatial arm (months 7-10)

Visium HD 8 µm bins are partial cells or two-cell mixtures; every spatial
model except scGPT-spatial keeps spot-level data out of pretraining, and
the memo found no spatial foundation model beating naive baselines on
local structure. So:

- Bins do not enter core pretraining.
- Two cheap designs are compared: a separate spatial head with a
  neighbourhood loss over the transcript model's embeddings (the
  scGPT-spatial pattern), and a graph model over those embeddings (the
  Novae pattern).
- Evaluation: deconvolution agreement with cell2location and RCTD on our
  captures, using the segmented-cell arbiter Mouse_Compensatory already
  has (`cells`, `deconvolution_arbiter.csv`); domain recovery on the
  McKey GSE307688 atlas (same platform, eight time points), with its
  published domains as labels.
- Kill rule: dropped unless it improves the arbiter agreement or the
  McKey domain recovery over cell2location/RCTD plus BANKSY by the
  pre-specified margin.

## 8. Gates, kill rules and margins (the pre-registration table)

| Gate | When | Rule | If failed |
| --- | --- | --- | --- |
| G0 corpus | end of phase 1 | ≥ 0.75 M curated mouse observations with raw counts and condition metadata; crosswalk ≥ 90 % at the eight-type level | publish the inventory; stop |
| G1 benchmark | end of phase 2 | ≥ 5 T3 contexts with computable r and ≥ 2 compartments each | T3 becomes descriptive; T1/T2/T4 carry the benchmark |
| KR1 embedding | phase 3, T1 and T2 | the model beats scANVI on macro-F1 by ≥ 0.02 **and** on rare-type F1 by ≥ 0.05 on ≥ 2 of 3 held-out mouse atlases; mixing no worse than Harmony on the repo's observed/expected table (the September rule) | the "ovary embedding" claim is dropped; the model is evaluated on T3 only |
| KR2 condition effect | phase 3, T3 | the model's context score (s against b1) exceeds the best baseline's by ≥ 0.05 and is > 0 on ≥ 4 of 6 members, on ≥ 2 held-out contexts | the model is reported as a negative result; G1-G4 deliverables stand |
| KR3 scale | data-size ablation | if the primary metrics change by < 0.01 between 200 K and 1 M cells, corpus growth stops | reported as the organ-specific scaling result |
| KR4 species | species ablation | mouse + human is kept only if it beats mouse-only on mouse tasks by ≥ 0.02 at matched compute | mouse-only is the released model |
| KR5 spatial | phase 4 | beats cell2location/RCTD + BANKSY on the arbiter or on McKey domains by ≥ 0.05 | spatial arm dropped |

Margins are deliberately modest and in the units of each metric; they
are the smallest differences the lab would act on. They are frozen in
`docs/preregistration.md` before phase 3 starts.

## 9. Phase 5: the Mouse_Compensatory use (months 8-11)

- **ooph vs sham as a held-out context.** The model (trained with the own
  data excluded) receives the sham nuclei of each compartment and emits
  the ooph state; the score is against the real ooph nuclei with r from
  the within-arm split half across the two mice. With two mice per arm
  the ceiling will be low and some compartments unscorable; that is
  reported, not hidden. The prediction itself is a hypothesis list: the
  genes the model expects to move in a compensating ovary, read against
  the lab's monotonic and pathway results.
- **The 120 d cohort.** Designed so the benchmark can score it: at least
  three mice per arm if the budget allows, ≥ 400 nuclei per major
  compartment per arm, vaginal cytology at collection (the repo's
  `estrous_stage` field exists for this), and the same library protocol
  as 30 d. It becomes the only two-time-point surgical context in
  existence.
- **The generalist pilot.** The `fm` and `fm-compare` stages run as
  designed and their table is the zero-shot arm's row in every
  comparison; nothing new is built for it.
- **The option the lab owns.** Every 2026 result says in-context
  perturbation data is what moves the score, and no ovary has any. An in
  vivo pooled screen in the mouse ovary (the liver and brain AAV
  precedents are the feasibility anchors) or a denser hormonal time
  course with ≥ 3 animals per point would make this group the only holder
  of in-context ovarian perturbation signal. It is raised here as a
  decision for the lab, with no part of this plan depending on it.

## 10. Phase 6: write-up and release (months 10-12)

- Paper 1, whatever happens to the model: the corpus, the crosswalk, the
  benchmark, the ceilings, the scaling and species ablations.
- Paper 2, only past KR1 or KR2: the model and its held-out-context
  results, with the negative-result alternative drafted in parallel.
- Release: corpus manifest and crosswalk under CC-BY; code under GPL-3 if
  derived from scPRINT-2 (arm A) or Apache-2 otherwise; weights on
  Zenodo and Hugging Face with a model card that states the ovary share
  of pretraining and every excluded dataset.

## 11. Timeline

| Month | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Phase 0 decisions, verification, pre-registration | ■ | | | | | | | | | | | |
| Phase 1 corpus and crosswalk | ■ | ■ | ■ | ■ | | | | | | | | |
| Phase 2 benchmark, ceilings | | ■ | ■ | ■ | ■ | | | | | | | |
| Phase 3 model arms and ablations | | | | ■ | ■ | ■ | ■ | ■ | ■ | | | |
| Phase 4 spatial arm | | | | | | | ■ | ■ | ■ | ■ | | |
| Phase 5 ooph vs sham; 120 d design | | | | | | | | ■ | ■ | ■ | ■ | |
| Phase 6 write-up, release | | | | | | | | | | ■ | ■ | ■ |

Effort: one computational lead throughout; the lab for the crosswalk
review (phase 1), the condition-defining gene lists (phase 2) and the
120 d design (phase 5); Biowulf for roughly 2,500 A100-hours over months
4-10.

## 12. Risks

| Risk | Why it matters | Mitigation |
| --- | --- | --- |
| Snippet-derived numbers in the memo are wrong | counts, margins and licences drive the design | the phase 0 checklist; nothing quoted outside the repo until re-read at the source |
| Mixing nuclei and cells hurts more than a token fixes | the own data are nuclei; Census mouse ovary has none | the modality ablation; intronic-inclusive counts from scBaseCount; the repo's nuclei QC rules |
| Label noise across atlases | T1 scores reflect the crosswalk, not the model | lab-reviewed crosswalk; confusion reported at the finest level; unmapped labels excluded from scoring |
| Too few animals per condition in public data | r is uncomputable, T3 shrinks | the ceilings table decides early (gate G1); the 120 d cohort is designed to be scorable |
| Compute or queue limits on Biowulf | the ablation grid is the expensive part | arm B first (cheap); serialise; drop seeds before dropping ablations |
| GPL-3 on arm A | derived code must be GPL; weights licence unstated | ask the scPRINT-2 authors in phase 0; arm B is the permissive fallback |
| Over-claiming on 2 vs 2 mice | the own data cannot carry a model claim | it is a held-out context with a reported ceiling; the claims rest on the public contexts |
| A 2026 Challenge result that changes the picture | if learned models beat retrieval on unseen contexts in late November, the T3 head design should follow them | re-read the results before phase 3; the benchmark does not change |

## 13. Relation to the Mouse_Compensatory repo

| Reused from Mouse_Compensatory | Where it goes here |
| --- | --- |
| nuclei QC rules, MAD bounds, doublet calls, ambient flags (`qc`) | the corpus pipeline |
| `config/markers.yaml` panels, negatives, hallmarks; `config/reference_labels.yaml` | the crosswalk seed and the T4 gene sets |
| `reference` stage (CellTypist / logistic on an atlas) | the T1 baseline and the crosswalk check |
| `fm`, `fm-compare` stages and `environment-fm-*.yml` | the zero-shot arm; the per-model conda pattern |
| `cells`, `deconvolution_arbiter.csv`, `niche` | the spatial arm's evaluation |
| `condition` stage, pseudobulk `de`, the monotonic and lab gene sets | the reading of the ooph prediction |
| `biowulf/{setup_env,run,submit}.sh`, `write_provenance` | the training jobs |
| the synthetic data generator and `test_smoke.py` pattern | the plumbing tests |

Proposed layout of this repo:

```
Virtual_Cell/
  README.md
  docs/            ovary_foundation_model_research.md, project_plan.md, preregistration.md, data/
  config/          corpus.yaml (manifest), crosswalk.yaml, benchmark.yaml, params.yaml
  ovary_fm/        corpus/ (download, qc, harmonise), bench/ (tasks, bundles, runner), models/ (arms A-C, heads), cli.py
  scripts/         one-off importers (scBaseCount metadata, GEO fetches)
  biowulf/         setup_env.sh, run.sh, submit.sh (train and bench stages)
  environment*.yml one per stack
  tests/           synthetic two-condition object, stub model, runner
```

The package, the config side files, the Biowulf scripts, the synthetic
corpus and the tests were scaffolded on 2026-11-04 (the stage table is in
the README); the preregistration draft is
[`preregistration.md`](preregistration.md) and the step-by-step route to
the first three tables is
[`first_results_checklist.md`](first_results_checklist.md). What the
scaffold does not do yet: run a real base model (the scPRINT-2, UCE and
TranscriptFormer branches are ported untested), parse Arc's scorer
output (`--scorer cell-eval2` writes the files and calls the CLI with
flags to be checked), the T2 integration and T4 gene-program tasks, the
condition head for arms A and B, and the spatial arm.
