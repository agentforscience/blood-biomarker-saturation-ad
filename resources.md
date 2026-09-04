# Resources Catalog

**Phase:** resource_finder (Phase 1) · **Date:** 2026-09-04
**Topic:** Benchmarking Blood-Based Protein Biomarker Panels for Preclinical Alzheimer's Detection

## Summary

| Resource | Count | Location |
|---|---|---|
| Papers screened (Europe PMC, with abstracts) | **548** | `papers/epmc_ranked.json` |
| Papers downloaded as full-text PDF | **44** | `papers/` (159 MB) |
| Papers deep-read chunk-by-chunk | **4** | see `literature_review.md` §2 |
| Individual-level datasets downloaded | **1** | `datasets/GSE275392/` |
| Reference/aggregate datasets downloaded | **10 tables** | `datasets/literature_reference/` |
| Code repositories cloned | **0** (justified) | `code/README.md` |
| Analysis scripts written | **6** | `tools/` |

---

## Papers

44 unique open-access PDFs (45 records; the Knight ADRC study appears as both preprint and
journal version, which resolve to the same file). Top 20 by relevance score; full list in `papers/README.md`, complete
metadata for all 548 screened records (incl. abstracts) in `papers/epmc_ranked.json`.

| Title | First author | Year | File | Key info |
|---|---|---|---|---|
| Optimizing timing and cost-effective use of plasma biomarkers in Alzhe | Chang HI et al. | 2025 | `2025_optimizing_timing_and_costeffective_use_of_plas` | cit=2 |
| Incremental value of plasma biomarkers in predicting clinical decline  | Khorsand B et al. | 2026 | `2026_incremental_value_of_plasma_biomarkers_in_predi` | cit=2 |
| High-sensitivity plasma proteomics reveals disease-specific signatures | Gong K et al. | 2025 | `2025_highsensitivity_plasma_proteomics_reveals_disea` | cit=13 |
| Plasma proteomic signatures of preclinical Alzheimer's disease in clin | Trelle AN et al. | 2026 | `2026_plasma_proteomic_signatures_of_preclinical_alzh` | cit=2 |
| Prospective study on clinical utility of plasma p-Tau217 and other bio | Ishiguro T et al. | 2026 | `2026_prospective_study_on_clinical_utility_of_plasma` | cit=1 |
| Real-world diagnostic performance of blood-based biomarkers for Alzhei | Vigneswaran S et al. | 2026 | `2026_realworld_diagnostic_performance_of_bloodbased_` | cit=1 |
| High-sensitivity plasma proteomics reveals disease-specific signatures | Gong K et al. | 2025 | `2025_highsensitivity_plasma_proteomics_reveals_disea` | cit=0 |
| Relative importance of blood-based biomarkers for Alzheimer's disease- | Kim KY et al. | 2026 | `2026_relative_importance_of_bloodbased_biomarkers_fo` | cit=0 |
| Breakpoints in Alzheimer's disease biomarkers and cognition across the | Hu M et al. | 2026 | `2026_breakpoints_in_alzheimers_disease_biomarkers_an` | cit=1 |
| Independent study demonstrates amyloid probability score accurately in | Fogelman I et al. | 2023 | `2023_independent_study_demonstrates_amyloid_probabil` | cit=30 |
| Development of thresholds and a visualization tool for use of a blood  | Verberk IMW et al. | 2024 | `2024_development_of_thresholds_and_a_visualization_t` | cit=18 |
| Predicting continuous amyloid PET levels with CSF and plasma brain-der | Trudel L et al. | 2026 | `2026_predicting_continuous_amyloid_pet_levels_with_c` | cit=0 |
| Sex differences in Alzheimer's disease plasma biomarker levels and cli | Milà-Alomà M et al. | 2026 | `2026_sex_differences_in_alzheimers_disease_plasma_bi` | cit=1 |
| Individualized prediction of clinical progression to dementia using pl | Honey MIJ et al. | 2025 | `2025_individualized_prediction_of_clinical_progressi` | cit=1 |
| Plasma p-tau217, p-tau181, and Aβ42 predict amyloid PET positivity in  | Bao R et al. | 2026 | `2026_plasma_ptau217_ptau181_and_a42_predict_amyloid_` | cit=0 |
| Blood biomarkers to improve dementia diagnostic accuracy: a cross-sect | Kwon J et al. | 2026 | `2026_blood_biomarkers_to_improve_dementia_diagnostic` | cit=0 |
| Longitudinal plasma p-tau217 as a marker for tracking progression and  | Ahn J et al. | 2026 | `2026_longitudinal_plasma_ptau217_as_a_marker_for_tra` | cit=1 |
| Body mass index and blood volume influence plasma biomarkers and posit | Jacobs T et al. | 2025 | `2025_body_mass_index_and_blood_volume_influence_plas` | cit=13 |
| Differential associations of NFL and GFAP with neuropsychiatric sympto | Wu J et al. | 2026 | `2026_differential_associations_of_nfl_and_gfap_with_` | cit=0 |
| Plasma phosphorylated tau 217 detects amyloid-β in neuronal synuclein  | Smith AM et al. | 2026 | `2026_plasma_phosphorylated_tau_217_detects_amyloid_i` | cit=3 |
**Nine highly-ranked papers had no retrievable OA PDF** (paywalled, preprint-server-only, or
publisher bot-blocked — notably the Bio-Hermes 295-protein study on MDPI). Their abstracts are in
`papers/epmc_ranked.json` and were used in the review; they are listed at the end of
`papers/README.md`.

---

## Datasets

| Name | Source | Size | Task | Location | Notes |
|---|---|---|---|---|---|
| **GSE275392** | NCBI GEO | 53 samples × 1,305 proteins | Binary: amyloid positivity | `datasets/GSE275392/` | **Only public individual-level plasma-proteome + amyloid-status dataset found.** SomaScan RFU. ⚠️ APOE perfectly confounded with amyloid; no p-tau217/NfL. |
| Age-dependent multiplex supp. tables | figshare (Alz Res Ther 2026) | 10 files, ~180 KB | Reference/calibration | `datasets/literature_reference/age_dependent_multiplex/` | **Aggregate only.** Contains ADNI benchmark: p-tau217 AUC 0.904, Aβ42/40 AUC 0.831 (n=1317). |

figshare download ids: S1=66944695, S2=66944698, S3=66944701, S4=66944704, S5=66944707,
S6=66944713, S7=66944716, S8=66944719, S9=66944722, S15=66944725
(`https://ndownloader.figshare.com/files/<id>`).

Full details, limitations, loading code, and validation results: `datasets/README.md`.
Data files are git-ignored (`datasets/.gitignore`); `download.sh` reproduces them.

---

## Code Repositories

**None cloned.** No `code_references` were specified in the research topic, and a five-query
GitHub search surfaced no repository that would be used. Full justification and a description of
the six analysis scripts written in this phase: `code/README.md`.

---

## Resource Gathering Notes

### Search strategy
1. **Paper-finder service was unavailable** (`localhost:8000` not running) — fell back to manual
   search as instructed.
2. Chose **Europe PMC** as primary source rather than arXiv/Papers-with-Code: this is a
   biomedical topic with essentially no arXiv presence, and Europe PMC provides abstracts,
   citation counts, OA flags and PMC full-text links in one API.
3. Ran **15 structured queries** (`tools/queries.json`) spanning: core analytes; preclinical/CU
   populations; broad proteomics platforms (SomaScan/Olink/NULISA); ML and feature selection;
   panel-size and parsimony language; screening/trial-enrichment; head-to-head comparisons;
   named cohorts. Cursor-paged, 2 pages × 25 per query → **548 unique records**.
4. Scored and ranked by weighted keyword match + recency + log-citations + cross-query agreement
   (`tools/rank.py`), then downloaded the top 45 OA PDFs.
5. Deep-read the four most decisive papers with the PDF chunker.

### Selection criteria
Prioritised papers that (a) study **cognitively unimpaired / preclinical** populations,
(b) use **amyloid positivity** as the endpoint, (c) report **AUCs for panels of differing size**,
and (d) compare targeted assays against broad proteomic panels. Papers on animal models,
in vitro work, or other diseases were down-weighted.

### Challenges encountered

1. **Paper-finder service down** → manual Europe PMC pipeline built instead.
2. **Silent API failure:** Europe PMC returns HTTP 200 with *zero* results when passed
   `sort=CITED_BY_COUNT desc` (invalid syntax). The first search run returned 0 hits for all 15
   queries and looked like a legitimate empty result. Diagnosed by parameter bisection; removed
   the sort and ranked client-side.
3. **The central obstacle — data access.** No public individual-level dataset exists containing
   the core analytes + a broad panel + amyloid status in a preclinical cohort. Verified
   exhaustively, not assumed:

   | Source probed | Result |
   |---|---|
   | GEO (6 query formulations) | 1 usable dataset: **GSE275392**. Others were transcriptomics or wrong disease. |
   | Dryad (7 queries) | 4 candidates; **all inspected via API → supplementary PDFs/DOCX/TIFFs only**, no individual data |
   | Zenodo (3 queries) | Mirrors of the same Dryad records |
   | figshare (3 queries) | Aggregate supplementary tables only (downloaded as reference) |
   | OSF, HuggingFace, Kaggle | Nothing relevant (HF returned only MRI *image* datasets) |
   | GitHub (5 queries) | 4 queries returned **zero** repos; 1 returned an unrelated GWAS repo |
   | ADNI / A4 / BioFINDER / Knight ADRC / Bio-Hermes / SAMS-AMASS / UK Biobank / OASIS-3 / NACC | All require signed DUA + institutional review; **no anonymous download path** |
   | MDPI (Bio-Hermes paper), GAP Bio-Hermes portal | Bot-blocked / DNS-unresolvable in this environment |

4. **Encoding trap in GSE275392:** one of 53 per-sample files is UTF-16 while the other 52 are
   UTF-8. A naive `pd.read_csv` loop crashes with `UnicodeDecodeError`. Handled in
   `tools/build_gse275392.py`.
5. **Confounding in GSE275392:** discovered on inspection that all 18 amyloid-negatives are
   APOE33 and all 17 APOE44 subjects are amyloid-positive (by design — the study sampled
   homozygotes). This would have silently inflated any full-cohort result. Mitigation
   (APOE33-only primary stratum) is pre-registered in `planning.md`.
6. **`pyproject.toml` build failure:** the prescribed hatchling build backend fails without a
   package directory. Removed the `[build-system]` block; `uv add` then worked and still keeps
   dependency resolution local to this workspace (the stated purpose).

### Gaps and workarounds

| Gap | Workaround |
|---|---|
| No individual-level data with p-tau217/GFAP/NfL/Aβ42:40 | Direction **D2**: literature-calibrated simulation using ADNI/published AUCs and correlation matrices as anchors (`planning.md`) |
| GSE275392 has only n=53 | Repeated nested CV + bootstrap CIs; conclusions framed as compatible/incompatible, never confirmatory. Direction **D3** (meta-analytic synthesis) supplies the large-n evidence |
| GSE275392 confounded by APOE | Pre-registered APOE33-only primary stratum (n=36, 18/18 balanced) |
| SomaScan cannot measure p-tau217 | Reframed as a *feature*: it isolates the question "how much does proteomic **breadth** buy, absent the key epitope?" — directly the hypothesis's second half |

---

## Recommendations for Experiment Design

### 1. Primary datasets
- **GSE275392**, **APOE33-only stratum (n=36)** as primary; full cohort (n=53) secondary and
  labelled confounded.
- **ADNI/published AUCs** (`S7`, `S1`, `S2` tables) as calibration anchors and external
  benchmarks.

### 2. Baseline methods (priority order)
1. **age + sex + APOE-ε4 covariate-only model** — *mandatory*; published CU AUC 0.748–0.773.
   Routinely omitted in the literature; any panel not beating it has shown nothing.
2. Best single protein (k=1)
3. Compact panels k=2, 3, 5
4. Mid panels k=10, 20, 50
5. Full 1,305-protein panel
Models: L2/elastic-net logistic regression (primary) + random forest (secondary). **Not** deep
learning — at n=36 it only overfits (`planning.md` P7).

### 3. Evaluation metrics
- **AUC(k)/AUC(k_max) ratio with bootstrap CIs** — the quantity the "≥90%" claim names.
- AUC with CIs; DeLong for nested model comparison.
- PPV/NPV at a realistic ~20% preclinical amyloid prevalence — decisive for a *screening* claim,
  and not implied by AUC.

### 4. Code to adapt/reuse
`tools/build_gse275392.py` (data assembly, handles the encoding trap) and
`tools/validate_dataset.py` (correct CV structure with feature selection **inside** folds —
extend this into the full experiment rather than starting fresh).

### 5. Non-negotiable methodological constraints
- Feature selection **strictly inside** training folds. At n=36, selecting on all data will
  manufacture a result.
- Report the **entire AUC(k) curve**, including the region where it declines.
- Keep **endpoint discipline**: amyloid positivity ≠ cognitive decline. GFAP/NfL look useful for
  the latter and are near-useless for the former; conflating them will produce a wrong conclusion.

### 6. Expected outcome
The hypothesis is likely to be **supported, and by a stronger route than proposed**: the evidence
indicates 1–2 analytes (p-tau217, or pTau217/Aβ42) already reach the practical ceiling
(AUC 0.90–0.94), while 120–295-protein panels reach only 0.79–0.81. Two sub-claims are likely to
**fail** and should be reported as such: **NfL** does not beat age+sex+APOE for amyloid in CU
(it performs *worse*), and **MTBR-tau243** is a tau-tangle staging marker, not an early amyloid
marker.
