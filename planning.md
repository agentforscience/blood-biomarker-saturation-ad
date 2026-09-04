# Planning

> **Phase 0/1 note (experiment_runner, 2026-09-04):** Sections 0.x below were added at the start of
> the experiment_runner phase. Sections 1–5 are the resource_finder direction budget and
> pre-registered commitments, retained verbatim.

---

## 0. Motivation & Novelty Assessment

### 0.1 Why this research matters

Anti-amyloid therapies (lecanemab, donanemab) work best before symptoms appear, which turns
**preclinical amyloid detection in cognitively unimpaired (CU) people into a population-screening
problem**. Amyloid PET (~$5,000, limited scanners) and CSF (lumbar puncture) cannot scale to the
tens of millions of at-risk adults, so the field has moved to blood. Two incompatible strategies
are being pursued simultaneously: cheap targeted immunoassays measuring 1–5 analytes
(p-tau217, Aβ42/40, GFAP, NfL), and broad discovery proteomics measuring 120–7,000 proteins
(SomaScan, Olink, NULISAseq). These differ by **two orders of magnitude in cost per sample**
(~$50–150 vs ~$500–2,500) and by an equally large factor in analytical and regulatory complexity.
If breadth buys little, hundreds of millions of screening dollars and years of assay-development
effort are being misallocated. The beneficiaries of an answer are health systems designing
screening pathways, trial sponsors doing enrichment, and payers deciding what to reimburse.

### 0.2 Gap in existing work (from `literature_review.md`)

Three specific gaps, all confirmed across the 548-record corpus:

1. **No published AUC(k) saturation curve exists for amyloid positivity.** Every paper in the
   corpus compares a handful of hand-picked model nestings (e.g. "base vs base+p-tau217 vs
   base+p-tau217+GFAP"). Nobody plots out-of-sample AUC as an explicit function of panel size k
   over a real high-dimensional plasma proteome. The shape of that curve — where it rises, where
   it saturates, where it *declines* — is exactly what the hypothesis is about, and it has never
   been measured.
2. **The mandatory baseline is routinely omitted.** Trelle et al. (2026) show age+sex+APOE-ε4
   alone reaches AUC 0.748–0.773 in CU cohorts. Most panel papers never report it, so a reported
   AUC of 0.80 for a 295-protein panel (Bio-Hermes) is presented as a success when it is in fact
   barely distinguishable from three variables obtainable from a questionnaire and a cheek swab.
3. **"90% of the performance" is never defined.** The claim the hypothesis makes — and that the
   field makes informally — is stated as a ratio of AUCs. Nobody has noted that this ratio is
   **degenerate**: AUC has a floor of 0.5, so a coin flip already achieves 53% of AUC 0.94, and
   *any* model with AUC ≥ 0.846 automatically clears "90% of 0.94". The claim as literally stated
   is nearly unfalsifiable.

### 0.3 Our novel contribution

Four things, in increasing order of generality:

- **C1 — The first empirical AUC(k) saturation curve** for preclinical amyloid positivity on a
  real 1,305-protein plasma proteome (GSE275392), with feature selection strictly inside CV folds,
  repeated nested CV, bootstrap CIs, and a label-permutation null. We report the whole curve
  including its decline.
- **C2 — A metric correction that makes the hypothesis falsifiable.** We evaluate the "≥90%"
  claim under three ratios: raw `AUC(k)/AUC(k_max)`; **excess/Gini ratio**
  `(AUC(k)−0.5)/(AUC(k_max)−0.5)`, which removes the chance floor; and **covariate-anchored ratio**
  `(AUC(k)−AUC_cov)/(AUC(k_max)−AUC_cov)`, which asks how much of the *biomarker-attributable*
  signal a compact panel recovers. We show the three ratios can disagree by >40 percentage points
  on the same data, and argue the second and third are the ones that carry decision-relevant
  meaning. To our knowledge this reframing has not been made in this literature.
- **C3 — A ground-truth boundary analysis.** A literature-calibrated generative model, anchored to
  published ADNI/SAMS/AMASS effect sizes, that identifies the conditions (number and strength of
  additional informative proteins, their correlation with the core markers, training-set size,
  prevalence) under which the ≥90% claim holds or breaks. This is the only route to the
  hypothesis's *named* analytes, because no public individual-level dataset contains them.
- **C4 — A quantitative synthesis** of the (panel size → AUC) relationship across the corpus,
  stratified by cohort stage, endpoint, and platform, with heterogeneity reported rather than
  pooled away.

### 0.4 Experiment justification

| Exp | What it is | Why it is *necessary* |
|---|---|---|
| **D1** | AUC(k) saturation curve on GSE275392 (n=36 APOE33 primary; n=53 secondary), k∈{1,2,3,5,10,20,50,100,300,1305}, L2-logistic + elastic net + random forest, selection inside folds, 20×5 repeated CV, bootstrap + permutation | The only direction using **real individual-level measurements**. Tests the structural claim (signal concentrates in few analytes) on a genuine plasma proteome in the target population. Without it, the project is entirely simulation and citation. |
| **D1b** | Baselines: age+sex(+APOE), best single protein, a priori AD-protein panel, full panel | Establishes the reference every panel must beat. Gap #2 above. Without it, an AUC of 0.76 looks like a result rather than a tie with demographics. |
| **D2** | Literature-calibrated simulation with known ground truth; sweeps over n_informative_extras, extra effect size, core–extra correlation, n, prevalence; 5 metric variants | Covers the **named analytes** (p-tau217, Aβ42/40, GFAP, NfL, MTBR-tau243), which D1's platform physically cannot measure, and converts "what happened in one small cohort" into "under what conditions is the claim true". Also supplies the power analysis D1's n=36 cannot. |
| **D3** | Corpus extraction of (panel size, AUC, population, endpoint, platform) from 548 abstracts + 44 PDFs, plus a hand-verified core table; stratified synthesis | Supplies the **large-n external evidence** and checks whether D1/D2 conclusions match what the field has actually observed. Failure modes of D3 (publication bias, cross-study incomparability) are uncorrelated with those of D1 (small n) and D2 (calibration assumptions), so convergence is informative. |
| **D0** | Metric-behaviour analysis (analytic + applied to D1/D2/D3 outputs) | Necessary because without it the hypothesis cannot be falsified (gap #3). This experiment defines what the other three are even measuring. |

**Pre-registered predictions** (recorded before running the full experiments; the feasibility probe
in `tools/validate_dataset.py` had already shown a peak AUC ≈0.76 at k≈10–20 on the APOE33 stratum):
- P1: D1's AUC(k) rises to k≈5–20 then flattens/declines; full 1,305-protein panel does not win.
- P2: D1's entire curve sits below the published single-analyte p-tau217 benchmark (~0.90).
- P3: Under the raw AUC ratio the ≥90% claim is trivially satisfied nearly everywhere.
- P4: Under the covariate-anchored ratio it is satisfied for p-tau217-containing panels and
  **fails** for GFAP/NfL-only panels.
- P5: D3 shows no positive relationship between panel size and AUC for amyloid positivity in CU.

---

**Phase:** resource_finder (Phase 1)
**Date:** 2026-09-04
**Hypothesis under test:** A panel of 3–5 blood proteins (p-tau217, GFAP, NfL, Aβ42/40, possibly
MTBR-tau243) achieves ≥90% of the predictive performance of 20+ protein panels for detecting
preclinical amyloid positivity.

---

## 1. The binding constraint discovered in Phase 1

The single most important Phase-1 finding is a **data-access constraint** that shapes every
downstream decision:

> There is **no public, individual-level dataset** that contains *both* (a) the hypothesis's core
> analytes (p-tau217, GFAP, NfL, Aβ42/40) *and* (b) a broad 20+ protein panel *and* (c) amyloid
> status in a preclinical (cognitively unimpaired) cohort.

Every cohort that has all three — ADNI, BioFINDER, A4/LEARN, Knight ADRC, Bio-Hermes, UK Biobank,
SAMS/AMASS, Stanford — is behind a **data-use application** (weeks-to-months turnaround, signed
DUA, institutional affiliation). This was verified empirically, not assumed (see
`resources.md` §Challenges for the full search log: GEO, Dryad, Zenodo, figshare, OSF,
HuggingFace, GitHub, Kaggle, PRIDE, ArrayExpress).

What *is* publicly available at individual level:
- **GSE275392** — SomaScan 1,305 plasma proteins × 53 non-demented elderly, with `amyloid_status`.
  The only public plasma-proteome + amyloid-status dataset found. Downloaded and assembled.

What is publicly available at **aggregate** level (rich, and directly usable):
- Published AUCs, sensitivities/specificities, correlation matrices, and incremental-value deltas
  from ~30 downloaded papers, including ADNI benchmarks (p-tau217 AUC 0.904, Aβ42/40 AUC 0.831,
  n=1317) and A4-trial incremental-value estimates.

The three directions below are chosen to be **jointly sufficient** under this constraint: each
attacks the hypothesis from a different angle, and their failure modes are uncorrelated.

---

## 2. Directions retained (top 3)

### D1 — Empirical panel-size saturation curve on real proteomic data (GSE275392)
**What:** Measure out-of-sample AUC for amyloid positivity as a function of panel size k
(k = 1, 2, 3, 5, 10, 20, 50, 100, 300, 1305) on the SomaScan 1,305-protein plasma proteome,
using nested cross-validation with feature selection performed *strictly inside* the training
folds. Report AUC(k) / AUC(k_max) — the exact quantity the hypothesis is about.

**Why retained:** Only direction grounded in real individual-level measurements. Directly tests
the *structural* claim (performance saturates at small k) on a genuine high-dimensional plasma
proteome in a non-demented cohort — the hypothesis's target population.

**Known limitations (must be handled, not hidden):**
- **n = 53** (35 amyloid+, 18 amyloid−). Small. Requires repeated nested CV + bootstrap CIs;
  point estimates will be wide and must be reported as such.
- **APOE is almost perfectly confounded with amyloid status**: all 18 amyloid-negatives are
  APOE33, and all 17 APOE44 subjects are amyloid-positive (crosstab verified). A naive classifier
  can score well by learning APOE4, not amyloid. **Mitigation:** the pre-registered primary
  analysis is the APOE33-only stratum (18 amyloid+ vs 18 amyloid−, perfectly balanced, n=36);
  the full-cohort analysis is secondary and explicitly labelled confounded.
- SomaScan **cannot measure p-tau217** (no phospho-epitope aptamers), and NEFL is absent from
  this 1,305-plex. GFAP, MAPT, APP, APOE are present. So D1 tests panel-size *scaling*, not the
  specific named core-4 panel. That gap is exactly what D2 covers.

### D2 — Literature-calibrated simulation of the core-4 vs. large-panel comparison
**What:** Build a generative model of (p-tau217, GFAP, NfL, Aβ42/40, MTBR-tau243) plus a large
set of weakly-informative proteins, with marginal effect sizes and inter-marker correlations
calibrated to published values (ADNI S7 table, A4 incremental-value deltas, head-to-head AUC
papers). Then compute the panel-size curve under *known ground truth* and sweep the parameters
that matter: number of truly-informative extra proteins, their effect sizes, their correlation
with the core-4, sample size, and outcome prevalence.

**Why retained:** The only direction that can address the hypothesis's *named analytes*, because
no public individual-level data contains them. Converts the question from "what happened in one
cohort" into "under what conditions is the 90% claim true or false" — which is the more useful
and more generalizable answer. Also supplies the power analysis that D1's n=53 cannot.

**Limitation:** Simulation results are conditional on calibration assumptions. Must be presented
as a sensitivity/boundary analysis, never as empirical evidence. Calibration targets must be
traceable to specific downloaded papers.

### D3 — Quantitative meta-analytic synthesis of published panel comparisons
**What:** Extract, from the downloaded corpus, every reported (panel composition → AUC) pair for
detecting amyloid positivity, and compute the achieved ratio AUC_small / AUC_large. Stratify by
cohort stage (cognitively unimpaired vs. mixed/impaired), reference standard (amyloid PET vs.
CSF), and platform (Simoa/Lumipulse/NULISAseq/SomaScan/MS).

**Why retained:** Highest evidence-per-unit-effort, and it is the direction whose data actually
exists at scale. Early extraction already shows the hypothesis is broadly supported and, in one
case, *understated*:
- A4 trial (n=866+343, cognitively unimpaired): adding Aβ42/40 + GFAP + NfL on top of p-tau217
  gave only **1–3% AUC gain**.
- p-tau231 + GFAP study (n=155, dementia-free): "adding β-amyloid42/40 and NfL did **not** produce
  a better fitting model" — 2 markers sufficed.
- Bio-Hermes (n=988, **295 proteins**): discovery-proteomics RF/GB reached AUC **0.79–0.81** —
  *below* single-analyte plasma p-tau217 (~0.90 in ADNI). A 295-protein panel underperforming one
  targeted assay is direct, strong evidence for the hypothesis.
- Knight ADRC (n=3,232, 120-plex NULISAseq): only **8 of 120** proteins associated with amyloid PET.

**Limitation:** Cross-study AUCs are not directly comparable (different cohorts, prevalence,
reference standards, assay platforms). Must report heterogeneity and avoid naive pooling.

---

## 3. Directions considered and PRUNED

| # | Direction | Reason for pruning |
|---|-----------|--------------------|
| P1 | Apply for ADNI / A4 / BioFINDER / Bio-Hermes / UK Biobank individual-level data | **Infeasible within pipeline horizon.** All require signed DUA + institutional review; turnaround is weeks to months. Verified: no anonymous download path exists for any of them. |
| P2 | Use UK Biobank Olink proteomics (~3,000 proteins, large n) | Same gating as P1, plus access fee. Only GWAS/pQTL *summary* statistics are public — no individual-level protein × amyloid data. Amyloid PET is also unavailable for nearly all UKB participants. |
| P3 | Re-analyse CSF (rather than plasma) proteomics, which is more publicly available | Changes the research question. The hypothesis is explicitly about **blood-based** screening; CSF requires lumbar puncture and so cannot address the practical screening claim that motivates the work. |
| P4 | Use blood **transcriptomic** (RNA-seq/microarray) AD datasets, which are abundant in GEO | Wrong analyte class. Blood transcriptomics has consistently weaker amyloid discrimination than targeted protein assays, and the hypothesis is specifically about *protein* panels. Would not be an honest proxy. |
| P5 | Predict clinical AD diagnosis instead of amyloid positivity (far more public data) | Changes the endpoint. Hypothesis targets **preclinical amyloid positivity** in cognitively unimpaired people; diagnosis-based endpoints are contaminated by symptom-driven signal and would inflate all panels' apparent performance. |
| P6 | Build a new plasma proteomics dataset | Out of scope — no wet-lab capability. |
| P7 | Deep-learning / representation-learning architectures over the proteome | Not informative for the hypothesis and actively harmful at n=53 (guaranteed overfitting). Panel-size saturation is a question about *signal concentration*, not model capacity; regularised linear models + tree ensembles are the appropriate and standard comparators. |
| P8 | Dryad/Zenodo "Data from:" plasma-amyloid deposits (Schindler, AIBL Aβ42/40, ARIC) | **Inspected and rejected on contents.** All four candidate Dryad records contain only supplementary PDFs/DOCX/TIFFs (summary tables and figures), not individual-level data. File listings verified via Dryad API. |

**Standing rule:** the search space is not to be expanded later unless new evidence invalidates
this ranking. If that happens, the ranking must be updated here and the change explained in
STATE.md.

---

## 4. How the three directions combine into one answer

- **D3** establishes what the field has already observed (breadth, external validity).
- **D1** tests the saturation claim on real high-dimensional plasma proteomic data (internal
  validity on genuine measurements, weak power).
- **D2** determines the conditions under which the ≥90% claim holds or breaks for the specific
  named analytes (mechanistic generality, no empirical grounding).

Convergence across all three would be a substantially stronger result than any one alone.
Divergence is itself the interesting finding and must be reported rather than reconciled away.

## 5. Pre-specified analysis commitments (to prevent post-hoc drift)

1. Primary D1 stratum is **APOE33-only (n=36, 18/18)**; full-cohort result is secondary and
   labelled confounded.
2. Feature selection occurs **inside** training folds only. Any selection on the full dataset is
   reported separately and explicitly labelled as optimistically biased.
3. Primary metric: **AUC ratio** AUC(k) / AUC(k_max), with bootstrap CIs — this is the quantity
   the "≥90%" hypothesis names.
4. Report the **whole curve**, not just the k values that favour the hypothesis.
5. Because n=53 yields wide CIs, D1 conclusions are stated as compatible/incompatible with the
   hypothesis, never as confirmation.
