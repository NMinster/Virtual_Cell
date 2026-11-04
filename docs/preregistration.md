# Preregistration: the ovary model, its benchmark and its kill rules

Nov 4, 2026 · **DRAFT, not yet frozen.** This file is frozen by one commit
before any model arm is trained (phase 3 of
[`project_plan.md`](project_plan.md)); the freezing commit's hash goes
here: `frozen_at: <hash>`. After that, every change is a row in the
deviations log (section 10) with a reason, and the tables are reported
under both the frozen and the changed rule.

The machine-readable copy of sections 3-6 is
[`config/benchmark.yaml`](../config/benchmark.yaml); where the two
disagree, this file wins and the yaml is fixed.

## 1. Hypotheses

| # | Hypothesis | Predicted direction | What falsifies it |
| --- | --- | --- | --- |
| H1 | A mouse-native model pretrained on the ovary corpus carries ovary cell identity across datasets better than scANVI trained on the same corpus | higher macro-F1 and rare-type F1 on leave-one-atlas-out transfer | KR1 margins not met on ≥ 2 of 3 held-out mouse atlases |
| H2 | Given a compartment's control cells only, the model predicts the perturbed state of a held-out context better than retrieval-and-transfer of the same condition family from other datasets | higher `vcc2026`-style context score against b1 | KR2 margins not met on ≥ 2 held-out contexts |
| H3 | At ovary scale (0.2-2 M cells) model performance saturates early | a flat data-size curve past 200 K cells | a gain ≥ 0.01 per primary metric between 200 K and 1 M cells |
| H4 | Adding human ovary data to pretraining does not help mouse tasks at matched compute | mouse-only ≥ mouse + human | mouse + human beats mouse-only by ≥ 0.02 |
| H5 | A spatial head over the transcript model does not beat cell2location/RCTD + BANKSY on Visium HD | no gain on the segmented-cell arbiter or on the McKey domains | a gain ≥ 0.05 |

H1 is expected to fail (the September 2026 memo and the 2026 benchmarks say
so); it is kept because its failure is the result that licenses dropping
the "ovary embedding" claim. H2 is the hypothesis the project exists for.

## 2. Data

- **Corpus**: the datasets of `config/corpus.yaml` with `enabled: true`
  and `held_out: false`, re-processed from raw counts through `ovary-fm
  qc` with the parameters in `config/params.yaml` as committed at the
  freeze. Cancer and cell-line samples are excluded from the core corpus.
  Spatial bins never enter core pretraining.
- **Held-out datasets** (never in any training run, including scVI and
  the gene-space baseline when they are the test set): every dataset in
  turn for T1 (leave-one-atlas-out); every context listed in section 4
  for T3; the own oophorectomy data always (`held_out: true`).
- **Hyperparameter selection** uses a validation split inside the
  training datasets only (one dataset held out of the training pool,
  rotated), never a benchmark held-out set.
- **Unmapped labels** (crosswalk status `unmapped`) stay in the corpus for
  pretraining and are excluded from T1 scoring on both sides.
- **Gene space**: the union of symbols over the pooled datasets; a gene a
  dataset did not measure is zero there and is flagged in `var`.

## 3. T1: label transfer

- Level: the eight major types (`Oocytes, Granulosa, Luteal, Theca,
  Stromal, Endothelial, Epithelial, Immune`); the fine level is reported
  descriptively once the crosswalk covers ≥ 50 % of cells at that level.
- Split: leave one labelled atlas out; the held-out own data is scored as
  a test set only.
- Classifier for an embedding: kNN (k = 15) inside the embedding. Gene-space
  baseline: logistic regression (C = 0.1, class-balanced) on the 2,000
  most variable genes of the training set, log-normalised and scaled.
  Primary baseline: scANVI (arm C). Secondary: CellTypist, scPoli.
- Metrics: macro-F1 over the test types with ≥ 20 cells; rare-type F1 =
  mean F1 over `Oocytes, Luteal, Theca, Epithelial` where present;
  accuracy reported, not used in a rule.

## 4. T3: condition effect in a held-out context

**Unit.** A bundle is one (context, perturbed condition, compartment).
The model receives the context's control cells of that compartment and
must emit integer-count cells for the perturbed condition.

**Contexts** (frozen list; one row per entry of `config/benchmark.yaml`
`t3.contexts`; a context whose dataset is not in the corpus at the
freeze is dropped, not replaced):

| Context | Dataset | Family | Control | Perturbed |
| --- | --- | --- | --- | --- |
| isola_ageing | isola2024 | ageing | age_3mo | age_9mo |
| wang_suh_ageing | wang_suh_estropause (GSE267729) | ageing | young | peri_regular, peri_irregular, post_estropausal |
| gaylord_ageing | gaylord2025_mouse (GSE290742) | ageing | young | aged |
| huang_ovulation | huang2025 (GSE294534) | ovulation | hcg_0h | hcg_4h, hcg_12h |
| schuh_ovulation_invivo | schuh2024_ovulation (GSE255274) | ovulation | invivo_hcg_0h | invivo 3, 6, 9, 12 h |
| schuh_ovulation_invitro | schuh2024_ovulation (GSE255274) | ovulation | invitro_hcg_0h | invitro 3, 6, 9, 12 h |
| morris_cycle | morris2022 (OSF 924fz, inDrop) | cycle | diestrus | proestrus, estrus, metestrus, pp10_lactating, pp10_nonlactating |
| mouse_comp_surgery (own, held out) | mouse_comp_30d | surgery | sham | ooph_30d |

Contexts added before the freeze (section 10 lists them) follow the same
template; the planned additions are the VCD subset, the Galligos 2026
cohort and the knockout datasets. Datasets deposited without annotations
(Wang/Suh, Schuh) carry major types transferred from the labelled
datasets; the transfer is the T1 logistic regression at confidence
≥ 0.5, and those datasets are never T1 train or test sets.

**Compartments**: the eight major types. A bundle is built only when
both groups have ≥ 50 cells (`min_cells_hard`); the target is Arc's 400.
A built bundle is **scored** only when ≥ 5 real DE genes remain after
the exclusions (`min_de_genes_to_score`): a compartment that does not
respond cannot tell models apart, and the empty-set conventions (Jaccard
of two empty sets = 1, a zero delta's PDS = 0.5) would otherwise make the
control copy look perfect there. Bundles below either gate are listed
in the ceilings table with their reason; only scored bundles enter the
context score and the PDS pools, and PDS itself is reported only when
its pool holds ≥ 5 scored bundles (`pds.min_pool`).

**Members** (six; the `lite` implementation in
`ovary_fm/bench/metrics.py` follows the published `vcc2026`
specification and is replaced by `cell-eval2` wherever the two can be
run side by side):

| Member | Definition | Direction |
| --- | --- | --- |
| `pds_cosine` | cosine distance between the predicted and every real pseudobulk delta (group-summed counts → 5 × 10⁴ → log1p, control subtracted) of the compartment's pool; 1 − mid-rank/(n − 1); a zero delta scores 0.5 | higher better |
| `expr_mse` | mean squared error between predicted and real log pseudobulk | lower better |
| `de_direction_fidelity` | k / max(n_pred, n_real), k = genes significant in both with the same sign | higher better |
| `de_direction_reach` | the deepest prefix of the real DE list, ordered by the prediction's \|log2FC\|, with ≥ 90 % sign agreement, over n_real | higher better |
| `de_sig_jaccard` | Jaccard of the significant sets; empty union = 1 | higher better |
| `de_lfc_nmae` | Σ\|lfc_pred − lfc_real\| / Σ\|lfc_real\| over the real DE genes; only when ≥ 10 of them remain after exclusions | lower better |

DE: two-sided Wilcoxon rank-sum per gene on CPM (10⁶), genes tested only
when their mean in the real control cells exceeds 5 CPM, log2 fold change
on arithmetic means (ε = 10⁻⁹), Benjamini-Hochberg within the bundle at
0.05; the real control is the comparison group for both the prediction
and the real cells.

**Rescaling.** Every member: s = (u − b0) / (r − b0). b0 is the control
cells copied (no effect), r the replicate. The expression member is
clamped to [0, 1], the NMAE member floored at −6, nothing else is
clipped. Bundle score = mean of the available members; context score =
mean of its scored bundles; overall = mean over contexts. The scale is
anchored on b0, not on the transfer baseline, because wherever the
transfer already sits near the ceiling the denominator r − b1 collapses
and the scores explode (seen on the first synthetic run); Arc anchors on
a near-null baseline for the same reason.

**Baselines.** b0 = the control cells copied: 0 on every member by
construction. b1 = retrieval and transfer: the per-gene fold change of
the same condition family in other datasets (the nearest time point when
both carry one; the source's significant genes only, 1.0 elsewhere;
geometric mean over sources), applied to this context's control cells
and Poisson-resampled. b1 is scored exactly like a model and its row is
the bar: where b1 exists a model must exceed b1's score, where it does
not a model must exceed 0.

**Replicate.** The real perturbed cells are split into two halves five
times (seeds 0-4); the halves follow animals whenever a group has ≥ 2
animals, otherwise cells, and the ceilings table says which. Half A is
scored as a submission against half B with each half's own control cells;
r is the mean over seeds and its spread is reported. For the expression
member the half-vs-half error is divided by 4 (it is the sum of two
half-depth variances, each twice the full group's), so r is the
noiseless full-depth floor; the DE members keep the half-depth value,
which is conservative.

**Excluded genes** (removed from every member; frozen lists):

- global: `Xist, Tsix`;
- ovulation: `Egr1, Fos, Fosb, Jun, Junb, Nr4a1, Nr4a2, Nr4a3, Ier2, Ier3,
  Atf3, Btg2, Ptgs2, Areg, Ereg, Pgr` (the LH/hCG immediate-early
  transcription factors and the four hallmark LH targets: the genes a
  model told "4 h after hCG" could write down from a textbook);
- genotype contexts: the knocked-out gene and its immediate pathway
  neighbours, listed per context when the context is added;
- cycle, ageing, chemical menopause, drug, surgery: no family-specific
  exclusions.

A list is extended only before the family's first scoring and the
extension is logged.

## 5. Models and arms

| Arm | What | Compute cap |
| --- | --- | --- |
| A | scPRINT-2 medium, continued pretraining (a few epochs) on the corpus | 500 A100-hours per run |
| B | ≤ 30 M parameters from scratch, scPRINT-2 small or Tahoe-x1 70 M configuration | 300 A100-hours per run |
| C | scVI then scANVI on the corpus (the primary baseline), plus a masked-autoencoder MLP | 20 A100-hours |
| Z | zero-shot scPRINT-2, UCE, TranscriptFormer through `ovary-fm fm` | inference |

For T3 the condition head is a conditional count decoder on (control-set
embedding, condition token) with control-matched dispersion in the
sampling step; a set-attention transition module is the second design
and is tried only if the first beats b1. Seeds: 0, 1, 2 for every arm;
tables report the mean and the range.

## 6. Decision rules

| Rule | Reads | Passes when | If it fails |
| --- | --- | --- | --- |
| G0 corpus | `tables/harmonise/manifest.csv` | ≥ 0.75 M curated mouse cells with raw counts and condition metadata; crosswalk coverage ≥ 90 % at the major level | publish the inventory; stop |
| G1 benchmark | `tables/bench/ceilings.csv` | ≥ 5 T3 contexts with a computable r and ≥ 2 scorable compartments each | T3 becomes descriptive; T1 carries the benchmark |
| KR1 embedding | `tables/bench/t1_*.csv` | the model beats scANVI on macro-F1 by ≥ 0.02 **and** on rare-type F1 by ≥ 0.05 on ≥ 2 of 3 held-out mouse atlases, and mixes datasets no worse than Harmony on the observed/expected table | the embedding claim is dropped |
| KR2 condition effect | `tables/bench/scores_*_summary.csv` | the model's context score exceeds b1's (or 0 where no b1 exists) by ≥ 0.05 and the model exceeds that baseline on ≥ 4 of 6 members, on ≥ 2 held-out contexts | the model is a negative result; corpus, crosswalk, benchmark and ceilings are the paper |
| KR3 scale | the data-size ablation | the primary metrics change by ≥ 0.01 between 200 K and 1 M cells | corpus growth stops; the flat curve is reported |
| KR4 species | the species ablation | mouse + human beats mouse-only on mouse tasks by ≥ 0.02 at matched compute | mouse-only is the released model |
| KR5 spatial | the spatial arm | ≥ 0.05 gain on the arbiter or the McKey domains over cell2location/RCTD + BANKSY | spatial arm dropped |

"Beats" means the mean over three seeds exceeds the comparator's by the
margin; no significance test is attached to 2 vs 2 animals or to three
seeds, and none is claimed.

## 7. Ablations (all pre-registered, three seeds each)

1. Data size: 50 K / 200 K / 1 M / all cells, arm B.
2. Species: mouse-only vs mouse + human at matched compute, arms A and B,
   scored on mouse tasks only.
3. Modality: nuclei-only vs cells + nuclei with the suspension token,
   scored on nuclei tasks.
4. Base: arm A vs arm B at the same data and compute.

## 8. What is not allowed

- Looking at a held-out context's perturbed cells before the model that
  is scored on it is frozen.
- Changing an exclusion list, a gate or a margin after a score has been
  seen, except through section 10.
- Choosing which baseline to report after the fact: both b0 and b1 are
  always in the table.
- Reporting a clipped 2025-style profile instead of the unclipped one.
- Treating the stub model's score as evidence of anything.

## 9. What is reported regardless of outcome

The ceilings table (what is reproducible between animals per compartment
and condition), the baseline scores, the T1 table with the gene-space and
scANVI rows, the data-size and species curves, and every deviation.

## 10. Deviations log

| Date | Section | Change | Reason | Scores seen before the change? |
| --- | --- | --- | --- | --- |
| 2026-11-04 | all | file created as a draft | — | no |
| 2026-11-04 | §4 | morris_cycle: the deposit is OSF project 924fz (inDrop; cycling, four estrous stages, two post-partum states), not GEO; the hcg_8h arm the memo listed is not in the deposit's description and leaves the context; the two post-partum states enter | the OSF record was read | no |
| 2026-11-04 | §4 | contexts added from the GEO records read that day: wang_suh_ageing (GSE267729), gaylord_ageing (GSE290742), schuh_ovulation_invivo and _invitro (GSE255274); morris_cycle's accession marked as unconfirmed (GSE192722 is a different study) | the records became readable; both families now have a second dataset, so b1 exists | no |
| 2026-11-04 | §4, §6 | the scale is anchored on b0; b1 is scored as a model and is the bar; the expression replicate is corrected to full depth (÷ 4) | anchoring on b1 made the denominator collapse where the transfer sits near the ceiling (b0 at −34 on one member of the synthetic run); the half-depth expression replicate was noisier than the control copy on small effects | synthetic scores only |
| 2026-11-04 | §4 | signal gate added (≥ 5 real DE genes after exclusions to score a bundle); PDS requires a pool of ≥ 5 scored bundles | the first synthetic run on Biowulf showed null compartments dominating the aggregate through the empty-set conventions and two-bundle PDS pools | synthetic scores only |
