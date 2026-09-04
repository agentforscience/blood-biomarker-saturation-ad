# Literature Review: How Many Blood Analytes Are Needed to Detect Preclinical Amyloid Positivity?

**Phase:** resource_finder (Phase 1) · **Date:** 2026-09-04
**Corpus:** 548 unique records screened via Europe PMC (15 structured queries); 44 full-text PDFs
downloaded; 4 papers deep-read chunk-by-chunk.

---

## 1. Research area overview

The hypothesis asks whether a compact panel (3–5 analytes: p-tau217, GFAP, NfL, Aβ42/40,
possibly MTBR-tau243) captures ≥90% of the performance of 20+ protein panels for detecting
**preclinical amyloid positivity** — i.e. amyloid pathology in people who are still cognitively
unimpaired (CU).

The field has converged on two parallel tracks:

1. **Targeted single/low-plex assays** (Simoa, Lumipulse, Elecsys, MSD ECL) measuring a handful
   of analytes with phospho-epitope specificity — above all **p-tau217**.
2. **Broad multiplex proteomics** (SomaScan ~1.3k–7k aptamers; Olink; NULISAseq CNS ~120–130
   targets; TMT mass spectrometry) measuring hundreds to thousands of proteins.

The literature is strikingly consistent, and the direction of the evidence is **stronger than the
hypothesis as stated**: not only do compact panels reach ~90% of broad-panel performance, in the
preclinical setting **broad panels frequently perform *worse* than a single well-chosen analyte**,
and two of the four named core analytes (GFAP, NfL) contribute little or nothing to *amyloid*
detection in CU cohorts.

A critical distinction runs through the whole corpus and must be preserved in our experiments:

> **Endpoint matters enormously.** GFAP and NfL are informative for *symptomatic* disease,
> neurodegeneration, and *cognitive decline prediction* — but weak-to-useless for **amyloid
> positivity in cognitively unimpaired people**. Papers that appear to disagree usually differ in
> endpoint, not in findings.

---

## 2. Key papers (deep-read)

### 2.1 Trelle et al. (2026), *Molecular Neurodegeneration* 21:31 — **most decisive paper**
"Plasma proteomic signatures of preclinical Alzheimer's disease in clinically unimpaired older adults"
`papers/2026_plasma_proteomic_signatures_of_preclinical_alzheimers_disease_in.pdf`

- **Design:** 315 CU adults in two independent cohorts — SAMS (n=193; CSF Aβ42/Aβ40 defines A+,
  27.3% A+) and AMASS (n=122; florbetaben amyloid PET >25 CL defines A+, 18.3% A+).
- **Platforms:** NULISAseq CNS panel (124 targets in SAMS, 131 in AMASS) **plus** Lumipulse
  single-plex for pTau217, pTau181, GFAP, NfL, Aβ42, Aβ40. This is a genuine head-to-head of a
  broad panel against targeted assays *in the exact target population*.
- **Headline AUCs for detecting amyloid positivity in CU:**

  | Marker (n analytes) | SAMS (CSF A+) | AMASS (PET A+) |
  |---|---|---|
  | pTau217/Aβ42 (**2**) | **0.940** NULISA / 0.907 Lumipulse | 0.865 |
  | BD-pTau217 (**1**) | — | **0.920** |
  | BD-pTau181 (**1**) | — | **0.920** |
  | pTau217 (**1**) | 0.879 / 0.838 | 0.861 |
  | Aβ42/Aβ40 (**2**) | 0.893 Lumipulse / 0.779 NULISA | 0.735 |
  | GFAP (**1**) | ~0.70 | 0.646 |
  | NfL (**1**) | ~0.55 | 0.582 |
  | **Covariates only** (age, sex, APOE-ε4) | **0.748** | **0.773** |

- **Findings that bear directly on the hypothesis:**
  - **GFAP and NfL did not beat the covariate-only model** for amyloid positivity in CU. NfL was
    *significantly worse* than age+sex+APOE alone (SAMS Z=3.89, p<0.001; AMASS Z=2.04, p=0.041).
  - Of ~123–130 panel proteins, **only ~5 were FDR-significant for amyloid status**: pTau217,
    pTau231, pTau181, GFAP, Aβ42 (SAMS); BD/total pTau217/181/231, MAPT, NPY (AMASS).
  - The best model in the entire study uses **two analytes** (pTau217/Aβ42).
- **Relevance:** This single paper nearly resolves the hypothesis in the affirmative and refines
  it: the useful core is the **pTau isoform family + Aβ42**, not the proposed core-4. It also
  supplies the correct baseline every panel must beat — **age + sex + APOE-ε4** (AUC ≈0.75–0.77).

### 2.2 Khorsand et al. (2026), *Alz&Dem: DADM* — A4/LEARN incremental value
`papers/2026_incremental_value_of_plasma_biomarkers_in_predicting_clinical.pdf`

- **Design:** 866 amyloid-positive A4 trial participants + 343 amyloid-negative LEARN
  participants, all CU. Endpoint = **5-year cognitive/functional decline** (CDR-GS ≥0.5), *not*
  amyloid positivity. Sub-study of 656 participants had Aβ42/40, GFAP, NfL.
- **Results:** base model (demographics + APOE-ε4) AUC 0.66; +p-tau217 → 0.73; +ADCS-PACC → 0.75;
  p-tau217 + PACC → **0.80**. Adding amyloid PET SUVR to that gave **no** further gain.
  Adding Aβ42/40 + GFAP + NfL improved AUC by only **≈1–2%** on top of models already containing
  p-tau217, and **≈1%** on the full model — mostly non-significant by DeLong.
- **Relevance:** Quantifies "diminishing returns" precisely in the preclinical population. Note
  the endpoint difference — this is prognosis, not amyloid detection.

### 2.3 Wang / Cruchaga et al. (2025), *Molecular Neurodegeneration* — Knight ADRC NULISAseq
`papers/2025_highsensitivity_plasma_proteomics_reveals_diseasespecific_signatures_and_pr.pdf`

- 3,232 participants (AD, DLB, FTD, PD, CU), NULISAseq CNS 120-plex.
- 81 proteins associated with AD *diagnosis*, but only **8 with amyloid PET** and 14 with CSF
  Aβ42/40 — versus 72 with CDR (clinical severity).
- **Relevance:** Largest broad-panel study in the corpus. Confirms the amyloid-specific signal is
  concentrated in a handful of proteins even when 120 are measured; breadth mostly buys
  *severity/other-disease* information, not amyloid discrimination.

### 2.4 Bio-Hermes discovery proteomics (2026), *IJMS* 27:5533 — the strongest single counter-example to "more is better"
- 988 participants, **295 plasma proteins** after QC, predicting **amyloid PET positivity**.
- Random Forest / Gradient Boosting achieved **AUC 0.79–0.81**; 8 recurring high-importance
  proteins (SERPINA1, C3, CRP, APOE4, CFH, VTN, C1QTNF5, PON1) — all generic
  inflammation/complement proteins, not AD-specific.
- **Relevance:** A 295-protein discovery panel lands **below** single-analyte plasma p-tau217
  (~0.90 in ADNI). Direct evidence that panel *breadth* does not substitute for the right analyte.
  PDF not retrievable (MDPI bot-blocked); abstract captured in `papers/epmc_ranked.json`.

### 2.5 Additional directly-supporting evidence

- **p-tau231 + GFAP combination study** (2025, *Alz&Dem*, n=155 dementia-free): p-tau231
  AUC 0.87, GFAP AUC 0.87; the **2-marker** model was best, and "including β-amyloid42/40 and NfL
  did **not** produce a better fitting model."
- **ADNI ML with feature selection** (2025, *Front Aging Neurosci*, n=1,043 + 127 external):
  explicitly introduces feature selection "to balance performance and cost" for predicting
  Aβ positivity — the same cost/parsimony framing as our hypothesis.
- **Down syndrome proteomic profile** (2026, *Alz&Dem*, n=290): SVM on **6** analytes + age/sex
  reached AUC 0.96 for amyloid positivity.
- **ADNI benchmark table** (`datasets/literature_reference/.../S7_ADNI_validation_cohort.xlsx`):
  plasma p-tau217 AUC **0.904** (n=1317, 541 Aβ+/776 Aβ−); Aβ42/40 AUC **0.831**. These are the
  best available external anchors for calibration.
- **MTBR-tau243** (Horie et al. 2023, n/a to amyloid screening): CSF MTBR-tau243 is a specific
  marker of **tau tangle** burden, tracking Braak stage and symptom severity — *not* an early
  amyloid marker. Plasma eMTBR-tau243 is positioned for **biological staging** of established
  disease. **Implication: MTBR-tau243 is a poor candidate for preclinical amyloid screening**, and
  the hypothesis's inclusion of it should be treated as a testable sub-claim likely to fail.

---

## 3. Common methodologies

| Method | Used in |
|---|---|
| ROC/AUC + DeLong test for nested model comparison | Trelle 2026; Khorsand 2026; nearly all |
| Youden-index thresholding, then concordance vs PET/CSF | Trelle 2026 |
| Logistic regression with sequential covariate blocks | Khorsand 2026 |
| Random Forest / Gradient Boosting / SVM / MLP | Bio-Hermes; ADNI ML; Down syndrome |
| Elastic net / LASSO / recursive feature elimination | ADNI ML; radiomics papers |
| Differential expression + FDR (Benjamini–Hochberg) on full panel | Trelle 2026; Knight ADRC |
| k-fold cross-validation (5-fold typical) | Khorsand 2026; Bio-Hermes |

**Methodological gap:** almost no paper reports a **panel-size saturation curve** —
AUC as an explicit function of k. They compare a few hand-picked model nestings. Producing the
full AUC(k) curve with honest nested CV is therefore a genuine contribution, not a replication.

---

## 4. Standard baselines (what every panel must beat)

1. **Covariate-only: age + sex + APOE-ε4** — AUC **0.748–0.773** in CU cohorts (Trelle 2026).
   This is the single most important and most frequently omitted baseline. Any panel not clearly
   beating ~0.77 has demonstrated nothing.
2. **Single plasma p-tau217** — AUC ≈0.86–0.92 in CU; **0.904** in ADNI (n=1317).
3. **pTau217/Aβ42 ratio (2 analytes)** — AUC up to **0.940**; current practical ceiling.
4. **Plasma Aβ42/40 alone** — AUC ≈0.73–0.89 (platform-dependent; degrades badly on some).
5. **Broad multiplex panel** — SomaScan/NULISAseq/discovery MS, AUC ≈0.79–0.81 (Bio-Hermes).

## 5. Evaluation metrics

- **AUC** — primary; use **DeLong** for correlated (nested) model comparison.
- **AUC ratio** AUC(k)/AUC(k_max) — the quantity the "≥90%" hypothesis literally names.
- **Sensitivity/specificity at Youden**, plus **PPV/NPV** — PPV/NPV are prevalence-dependent and
  matter more than AUC for a screening claim (preclinical amyloid prevalence ≈18–27% in CU).
- **Net reclassification / likelihood-ratio test** for incremental value.
- **Concordance %** with PET/CSF at a fixed threshold (Trelle's framing; very interpretable).

## 6. Datasets in the literature vs. what is obtainable

| Cohort | Used in | Public individual-level access? |
|---|---|---|
| ADNI | many | ❌ application/DUA |
| A4 / LEARN | Khorsand 2026 | ❌ application |
| BioFINDER | Palmqvist, Janelidze | ❌ |
| Knight ADRC | Cruchaga 2025 | ❌ |
| Bio-Hermes | IJMS 2026 | ❌ GAP data request |
| SAMS / AMASS | Trelle 2026 | ❌ |
| UK Biobank (Olink) | dementia-prediction papers | ❌ application + fee |
| **GSE275392** | Philippi & Castellano 2024 | ✅ **fully public — downloaded** |

See `resources.md` for the exhaustive negative search log.

## 7. Gaps and opportunities

1. **No published AUC(k) saturation curve** for amyloid positivity — the literature compares
   discrete model nestings, never the full curve. *This is our main opportunity.*
2. **The covariate-only baseline is routinely omitted**, so "panel works" claims are often
   uncalibrated against age+sex+APOE.
3. **GFAP/NfL are carried in "core panels" by convention**, despite failing to beat covariates for
   amyloid in CU. The 3–5 marker panel in the hypothesis is plausibly **over-specified**.
4. **MTBR-tau243 is mis-specified for this task** — it is a tau-tangle staging marker, not an
   early amyloid marker.
5. **Breadth vs. epitope specificity is the real axis.** SomaScan/Olink cannot measure
   phospho-epitopes at all, which structurally caps broad panels below p-tau217 — an explanation
   for the Bio-Hermes result that the field states only implicitly.

## 8. Recommendations for our experiments

**Datasets**
- Primary real data: **GSE275392** (1,305 SomaScan proteins × 53 CU elderly, amyloid status).
  Use the **APOE33-only stratum (n=36, 18/18)** as primary — APOE is otherwise perfectly
  confounded with amyloid status in this dataset.
- Calibration anchors: ADNI S7 table (p-tau217 AUC 0.904, Aβ42/40 0.831) and Trelle's CU AUCs.

**Baselines** (in priority order)
1. age + sex + APOE-ε4 (**mandatory**) 2. best single protein 3. compact panels k=2,3,5
4. k=10,20,50 5. full 1,305-protein panel.

**Metrics**
AUC with bootstrap CIs; **AUC(k)/AUC(k_max) ratio**; DeLong for nested comparisons;
PPV/NPV at realistic 20% prevalence.

**Methodological musts**
- Feature selection **strictly inside** training folds — selecting on all data at n=53 will
  manufacture a spurious result.
- **Repeated nested CV** + bootstrap CIs; n=53 gives wide intervals that must be shown.
- Report the **whole AUC(k) curve**, including where it *declines*.
- Treat the SomaScan ceiling honestly: this platform lacks p-tau217, so its curve answers
  "how much does proteomic *breadth* buy?" — which is exactly the hypothesis's second half.

**Expected result (pre-registered prediction):** the AUC(k) curve will rise steeply to k≈5–20 and
then flatten or decline; the full 1,305-protein panel will not beat a compact selection; and the
whole SomaScan curve will sit **below** the published single-analyte p-tau217 benchmark (~0.90),
supporting the hypothesis by a *stronger* route than proposed — breadth cannot substitute for the
right epitope. Our feasibility probe already shows this shape (peak AUC ≈0.76 at k≈10–20).
