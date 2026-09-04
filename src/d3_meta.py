"""
D3 - Quantitative synthesis of the published (panel size -> AUC) relationship.

Two independent passes, deliberately kept separate:

  PASS A (automated, high recall / low precision)
      Deterministic regex extraction of AUC values, panel sizes, population,
      endpoint and platform from all 548 Europe PMC abstracts and from the
      full text of the 44 downloaded open-access PDFs.  Used to characterise
      the FIELD-WIDE distribution and to check that the curated table below is
      not a biased subsample.  No LLM was available in this environment
      (no API key present), so extraction is rule-based and its error rate is
      quantified by manual audit of a random sample.

  PASS B (curated, high precision)
      A hand-verified evidence table of studies that report AUCs for panels of
      DIFFERENT SIZES on the SAME cohort with amyloid positivity as the endpoint.
      Only these support a within-study retained-performance ratio, which is the
      quantity the hypothesis names.  Every row cites its source.

Outputs: results/d3/*.csv, results/d3/*.json
"""
from __future__ import annotations

import json
import os
import re
import sys
import warnings

import numpy as np
import pandas as pd
import scipy.stats as st

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

warnings.filterwarnings("ignore")
OUT = os.path.join(C.REPO, "results", "d3")
os.makedirs(OUT, exist_ok=True)

# --------------------------------------------------------------------------- #
# PASS A - automated extraction
# --------------------------------------------------------------------------- #
AUC_RE = re.compile(
    r"(?:AUC|AUROC|area[- ]under[- ]the[- ](?:ROC[- ])?curve|C[- ]statistic)"
    r"[^0-9\n]{0,40}?(0?\.\d{2,3})", re.I)

# Panel-size regex.  Three guards, each added after a specific false positive
# found during the extraction audit:
#   (a) separator excludes newlines - "BMC Geriatrics (2026) 26:888 \n biomarker
#       performance..." (a running page footer) was matching as an 888-protein panel;
#   (b) negative lookbehind for ':' / '(' / 'Page ' / 'no. ' kills volume, issue,
#       page and reference numbers;
#   (c) a following cohort noun ("888 participants") is not a panel size, so the
#       noun list is restricted to assay vocabulary only.
PANEL_RE = re.compile(
    r"(?<![:(\d.])\b(\d{1,4})[ ,-]{0,3}(?:plex|different +)?"
    r"(?:proteins?|analytes?|biomarkers?|targets?|protein markers?|assays?)\b", re.I)
FOOTER_RE = re.compile(r"(?:page|vol\.?|no\.?|pp?\.)\s*$", re.I)

PLATFORM_PATTERNS = {
    "somascan": r"somascan|slow ?off-?rate|aptamer",
    "olink": r"\bolink\b|proximity extension",
    "nulisa": r"nulisa",
    "simoa": r"\bsimoa\b|single molecule array",
    "lumipulse": r"lumipulse",
    "elecsys": r"elecsys|\broche\b",
    "msd": r"meso ?scale|\bMSD\b",
    "mass_spec": r"mass spectrom|\bLC-?MS|immunoprecipitation.{0,20}mass",
}
POP_PATTERNS = {
    "cognitively_unimpaired": r"cognitively unimpaired|cognitively normal|preclinical|"
                              r"asymptomatic|non-?demented|\bCU\b|\bCN\b",
    "impaired_or_mixed": r"\bMCI\b|mild cognitive impairment|dementia|memory clinic|"
                         r"AD patients|symptomatic",
}
ENDPOINT_PATTERNS = {
    "amyloid_status": r"amyloid (?:PET )?(?:positiv|status|pathology|burden)|"
                      r"A\W?β positiv|amyloid-?β positiv|\bA\+\b|abnormal amyloid",
    "clinical_diagnosis": r"diagnos|discriminat\w+ (?:AD|Alzheimer).{0,20}(?:from|vs)|"
                          r"classif\w+ (?:AD|Alzheimer)",
    "cognitive_decline": r"cognitive decline|progression to dementia|conversion|"
                         r"incident dementia|longitudinal decline",
}


def _flags(text, patterns):
    return {k: bool(re.search(v, text, re.I)) for k, v in patterns.items()}


def strip_html(s: str) -> str:
    return re.sub(r"<[^>]+>", " ", s or "")


def extract_from_text(text: str, source_id: str, kind: str) -> list[dict]:
    """Pull every (AUC, nearest panel-size mention) pair out of one document."""
    text = strip_html(text)
    if not text:
        return []
    aucs = [(m.start(), float(m.group(1))) for m in AUC_RE.finditer(text)]
    aucs = [(p, a) for p, a in aucs if 0.40 <= a <= 1.0]
    panels = [(m.start(), int(m.group(1))) for m in PANEL_RE.finditer(text)
              if not FOOTER_RE.search(text[max(0, m.start() - 12):m.start()])]
    panels = [(p, k) for p, k in panels if 1 <= k <= 10000]

    pop, end, plat = (_flags(text, POP_PATTERNS), _flags(text, ENDPOINT_PATTERNS),
                      _flags(text, PLATFORM_PATTERNS))
    rows = []
    for pos, a in aucs:
        near = [(abs(pp - pos), k) for pp, k in panels if abs(pp - pos) < 1500]
        near.sort()
        rows.append(dict(source_id=source_id, kind=kind, auc=a,
                         panel_size=near[0][1] if near else np.nan,
                         panel_dist_chars=near[0][0] if near else np.nan,
                         **{f"pop_{k}": v for k, v in pop.items()},
                         **{f"end_{k}": v for k, v in end.items()},
                         **{f"plat_{k}": v for k, v in plat.items()}))
    return rows


def pass_a() -> pd.DataFrame:
    recs = json.load(open(os.path.join(C.REPO, "papers", "epmc_ranked.json")))
    rows = []
    for r in recs:
        rows += extract_from_text(f"{r.get('title','')} {r.get('abstractText','')}",
                                  r.get("doi") or r.get("id"), "abstract")

    # full-text pass over the downloaded PDFs
    try:
        from pypdf import PdfReader
        pdir = os.path.join(C.REPO, "papers")
        pdfs = sorted(f for f in os.listdir(pdir) if f.endswith(".pdf"))
        for f in pdfs:
            try:
                txt = " ".join((p.extract_text() or "")
                               for p in PdfReader(os.path.join(pdir, f)).pages)
                rows += extract_from_text(txt, f, "fulltext")
            except Exception as e:  # a malformed PDF must not kill the pass
                print(f"    [pdf skip] {f}: {type(e).__name__}")
    except ImportError:
        print("    pypdf unavailable - abstract-only extraction")

    df = pd.DataFrame(rows)
    df["population"] = np.where(df.pop_cognitively_unimpaired & ~df.pop_impaired_or_mixed,
                                "CU_only",
                                np.where(df.pop_cognitively_unimpaired, "mixed",
                                         "impaired_or_unspecified"))
    df["endpoint"] = np.where(df.end_amyloid_status, "amyloid_status",
                              np.where(df.end_cognitive_decline, "cognitive_decline",
                                       np.where(df.end_clinical_diagnosis,
                                                "clinical_diagnosis", "unspecified")))
    plat_cols = [c for c in df.columns if c.startswith("plat_")]
    df["broad_platform"] = df[[c for c in plat_cols
                               if c.split("plat_")[1] in
                               ("somascan", "olink", "nulisa", "mass_spec")]].any(axis=1)
    df["targeted_platform"] = df[[c for c in plat_cols
                                  if c.split("plat_")[1] in
                                  ("simoa", "lumipulse", "elecsys", "msd")]].any(axis=1)
    df.to_csv(os.path.join(OUT, "passA_extracted.csv"), index=False)
    return df


# --------------------------------------------------------------------------- #
# PASS B - curated within-study evidence
# --------------------------------------------------------------------------- #
# Every row is transcribed from a source read in phase 1 (see literature_review.md
# S2 for the deep-read papers).  `k` = number of measured analytes in the model
# (covariates excluded).  `endpoint` restricted to amyloid positivity except where
# stated.  Studies are grouped by `study` so within-study ratios are computable.
CURATED = [
    # --- Trelle et al. 2026, Mol Neurodegener 21:31 -- SAMS, n=193 CU, CSF Ab42/40 A+
    dict(study="Trelle2026_SAMS", cohort="SAMS", n=193, population="CU",
         endpoint="amyloid_status_CSF", platform="NULISAseq+Lumipulse",
         panel="covariates(age,sex,APOE4)", k=0, auc=0.748,
         src="Trelle 2026 Mol Neurodegener 21:31"),
    dict(study="Trelle2026_SAMS", cohort="SAMS", n=193, population="CU",
         endpoint="amyloid_status_CSF", platform="Lumipulse", panel="NfL", k=1, auc=0.550,
         src="Trelle 2026"),
    dict(study="Trelle2026_SAMS", cohort="SAMS", n=193, population="CU",
         endpoint="amyloid_status_CSF", platform="Lumipulse", panel="GFAP", k=1, auc=0.700,
         src="Trelle 2026"),
    dict(study="Trelle2026_SAMS", cohort="SAMS", n=193, population="CU",
         endpoint="amyloid_status_CSF", platform="NULISAseq", panel="pTau217", k=1, auc=0.879,
         src="Trelle 2026"),
    dict(study="Trelle2026_SAMS", cohort="SAMS", n=193, population="CU",
         endpoint="amyloid_status_CSF", platform="Lumipulse", panel="Ab42/Ab40", k=2, auc=0.893,
         src="Trelle 2026"),
    dict(study="Trelle2026_SAMS", cohort="SAMS", n=193, population="CU",
         endpoint="amyloid_status_CSF", platform="NULISAseq", panel="pTau217/Ab42", k=2, auc=0.940,
         src="Trelle 2026"),
    dict(study="Trelle2026_SAMS", cohort="SAMS", n=193, population="CU",
         endpoint="amyloid_status_CSF", platform="NULISAseq",
         panel="full CNS panel (124 targets)", k=124, auc=0.940,
         src="Trelle 2026: only ~5/124 FDR-significant; best model is the 2-analyte ratio",
         note="panel-level AUC not separately reported; upper-bounded by the best in-panel model"),
    # --- Trelle et al. 2026 -- AMASS, n=122 CU, florbetaben PET >25 CL
    dict(study="Trelle2026_AMASS", cohort="AMASS", n=122, population="CU",
         endpoint="amyloid_status_PET", platform="mixed",
         panel="covariates(age,sex,APOE4)", k=0, auc=0.773, src="Trelle 2026"),
    dict(study="Trelle2026_AMASS", cohort="AMASS", n=122, population="CU",
         endpoint="amyloid_status_PET", platform="Lumipulse", panel="NfL", k=1, auc=0.582,
         src="Trelle 2026 (significantly WORSE than covariates, Z=2.04 p=0.041)"),
    dict(study="Trelle2026_AMASS", cohort="AMASS", n=122, population="CU",
         endpoint="amyloid_status_PET", platform="Lumipulse", panel="GFAP", k=1, auc=0.646,
         src="Trelle 2026 (does not beat covariates)"),
    dict(study="Trelle2026_AMASS", cohort="AMASS", n=122, population="CU",
         endpoint="amyloid_status_PET", platform="NULISAseq", panel="pTau217", k=1, auc=0.861,
         src="Trelle 2026"),
    dict(study="Trelle2026_AMASS", cohort="AMASS", n=122, population="CU",
         endpoint="amyloid_status_PET", platform="NULISAseq", panel="BD-pTau217", k=1, auc=0.920,
         src="Trelle 2026"),
    dict(study="Trelle2026_AMASS", cohort="AMASS", n=122, population="CU",
         endpoint="amyloid_status_PET", platform="Lumipulse", panel="Ab42/Ab40", k=2, auc=0.735,
         src="Trelle 2026"),
    dict(study="Trelle2026_AMASS", cohort="AMASS", n=122, population="CU",
         endpoint="amyloid_status_PET", platform="NULISAseq", panel="pTau217/Ab42", k=2, auc=0.865,
         src="Trelle 2026"),
    dict(study="Trelle2026_AMASS", cohort="AMASS", n=122, population="CU",
         endpoint="amyloid_status_PET", platform="NULISAseq",
         panel="full CNS panel (131 targets)", k=131, auc=0.920,
         src="Trelle 2026: best single panel analyte BD-pTau217; panel adds nothing beyond it",
         note="upper-bounded by best in-panel model"),
    # --- Bio-Hermes discovery proteomics, IJMS 2026 27:5533
    dict(study="BioHermes2026", cohort="Bio-Hermes", n=988, population="mixed(CU+MCI+AD)",
         endpoint="amyloid_status_PET", platform="discovery proteomics",
         panel="295-protein RF/GB", k=295, auc=0.800,
         src="IJMS 2026 27:5533 (reported range 0.79-0.81)"),
    dict(study="BioHermes2026", cohort="Bio-Hermes", n=988, population="mixed(CU+MCI+AD)",
         endpoint="amyloid_status_PET", platform="targeted immunoassay",
         panel="plasma p-tau217 (external benchmark)", k=1, auc=0.900,
         src="external ADNI/BioFINDER benchmark, NOT measured in Bio-Hermes",
         note="CROSS-STUDY comparison - not a within-study ratio"),
    # --- ADNI validation cohort (supplementary tables S4/S7, downloaded)
    dict(study="ADNI_S7", cohort="ADNI", n=1317, population="mixed(CN+MCI+DEM)",
         endpoint="amyloid_status", platform="targeted", panel="p-tau217", k=1, auc=0.904,
         src="datasets/literature_reference/.../S7_ADNI_validation_cohort.xlsx"),
    dict(study="ADNI_S7", cohort="ADNI", n=1317, population="mixed(CN+MCI+DEM)",
         endpoint="amyloid_status", platform="targeted", panel="Ab42/Ab40", k=2, auc=0.831,
         src="same table"),
    dict(study="ADNI_S4_LO", cohort="ADNI late-onset", n=677, population="mixed",
         endpoint="amyloid_status", platform="targeted", panel="p-tau217", k=1, auc=0.902,
         src="S4_greyzone_delong_thresholds.xlsx (grey zone excluded)"),
    dict(study="ADNI_S4_LO", cohort="ADNI late-onset", n=674, population="mixed",
         endpoint="amyloid_status", platform="targeted", panel="Ab42/Ab40", k=2, auc=0.811,
         src="S4_greyzone_delong_thresholds.xlsx"),
    # --- p-tau231 + GFAP, Alz&Dem 2025, n=155 dementia-free
    dict(study="Ptau231GFAP2025", cohort="dementia-free cohort", n=155, population="CU/dementia-free",
         endpoint="amyloid_status", platform="Simoa", panel="p-tau231", k=1, auc=0.870,
         src="literature_review.md S2.5"),
    dict(study="Ptau231GFAP2025", cohort="dementia-free cohort", n=155, population="CU/dementia-free",
         endpoint="amyloid_status", platform="Simoa", panel="GFAP", k=1, auc=0.870,
         src="literature_review.md S2.5"),
    dict(study="Ptau231GFAP2025", cohort="dementia-free cohort", n=155, population="CU/dementia-free",
         endpoint="amyloid_status", platform="Simoa", panel="p-tau231 + GFAP", k=2, auc=0.890,
         src="best model; adding Ab42/40 and NfL did NOT improve fit"),
    # --- Down syndrome proteomic profile, Alz&Dem 2026, n=290
    dict(study="DownSyndrome2026", cohort="Down syndrome", n=290, population="DS (special)",
         endpoint="amyloid_status", platform="targeted", panel="6 analytes + age/sex SVM",
         k=6, auc=0.960, src="literature_review.md S2.5",
         note="genetically enriched population; not generalisable to sporadic CU screening"),
    # --- Knight ADRC NULISAseq, Mol Neurodegener 2025 (breadth vs specificity)
    dict(study="KnightADRC2025", cohort="Knight ADRC", n=3232, population="mixed",
         endpoint="amyloid_status_PET", platform="NULISAseq 120-plex",
         panel="120-plex (8/120 associated with amyloid PET)", k=120, auc=np.nan,
         src="Cruchaga 2025; AUC not reported for the amyloid-PET endpoint",
         note="counts only - included for the informative-fraction estimate, not the AUC model"),
]


def pass_b() -> pd.DataFrame:
    df = pd.DataFrame(CURATED)
    df.to_csv(os.path.join(OUT, "passB_curated.csv"), index=False)
    return df


def within_study_ratios(cur: pd.DataFrame) -> pd.DataFrame:
    """
    For each study reporting both a compact (k<=5) and a large (k>=20) panel on
    the SAME cohort with the same endpoint, compute the three retained-performance
    ratios of the best compact panel vs the large panel.
    """
    rows = []
    for study, g in cur.dropna(subset=["auc"]).groupby("study"):
        small = g[(g.k >= 1) & (g.k <= 5)]
        large = g[g.k >= 20]
        if small.empty or large.empty:
            continue
        cov = g[g.k == 0].auc.max()
        cov = 0.5 if pd.isna(cov) else cov
        bs, bl = small.loc[small.auc.idxmax()], large.loc[large.auc.idxmax()]
        r = C.ratios(bs.auc, bl.auc, cov)
        rows.append(dict(study=study, cohort=bs.cohort, n=int(bs.n),
                         population=bs.population, endpoint=bs.endpoint,
                         small_panel=bs.panel, k_small=int(bs.k), auc_small=bs.auc,
                         large_panel=bl.panel, k_large=int(bl.k), auc_large=bl.auc,
                         auc_covariates=cov, **r,
                         beats_90_raw=r["ratio_raw"] >= 0.90,
                         beats_90_excess=r["ratio_excess"] >= 0.90,
                         beats_90_anchored=r["ratio_anchored"] >= 0.90))
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "within_study_ratios.csv"), index=False)
    return df


def panel_size_regression(a: pd.DataFrame, max_dist: int = 800) -> dict:
    """
    Is a larger panel associated with a higher reported AUC?
    Fit AUC ~ log10(panel size) on the automated extraction, stratified by
    endpoint and population.  Reported with heterogeneity, not pooled naively.
    """
    out = {}
    d = a.dropna(subset=["panel_size", "auc"]).copy()
    d = d[d.panel_dist_chars < max_dist]   # size mention must be nearby
    d["logk"] = np.log10(d.panel_size)
    for label, sub in [("all", d),
                       ("amyloid_endpoint", d[d.endpoint == "amyloid_status"]),
                       ("amyloid_endpoint_CU", d[(d.endpoint == "amyloid_status") &
                                                 (d.population.isin(["CU_only", "mixed"]))]),
                       ("broad_platform", d[d.broad_platform]),
                       ("targeted_platform", d[d.targeted_platform])]:
        if len(sub) < 12:
            out[label] = dict(n=len(sub), note="too few extractions for a stable fit")
            continue
        sl, ic, r, p, se = st.linregress(sub.logk, sub.auc)
        rho, prho = st.spearmanr(sub.panel_size, sub.auc)
        out[label] = dict(n=int(len(sub)), slope_per_decade=float(sl), intercept=float(ic),
                          r=float(r), p=float(p), stderr=float(se),
                          spearman_rho=float(rho), spearman_p=float(prho),
                          auc_median=float(sub.auc.median()),
                          auc_iqr=[float(sub.auc.quantile(.25)), float(sub.auc.quantile(.75))],
                          median_k=float(sub.panel_size.median()))
    return out


def document_level(a: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """
    Document-level analysis - more robust than pairing each AUC with the nearest
    number in the text.

    For every document we take (a) the MEDIAN AUC reported for an amyloid-status
    endpoint and (b) a PANEL-BREADTH CLASS derived from the platform keywords and
    from the largest plausible panel-size mention.

    The median, not the maximum, is the correct summary: broad-platform papers
    report many more AUCs per paper than targeted-assay papers, so comparing
    maxima compares order statistics of samples of different size and manufactures
    an advantage for the broad arm.  Both are computed; only the median is
    interpreted, and the number of AUCs per document is reported alongside.

    Breadth classes:

        'broad'    : SomaScan / Olink / NULISAseq / discovery mass-spec, or an
                     explicit mention of >=20 proteins/analytes/targets
        'targeted' : Simoa / Lumipulse / Elecsys / MSD with no broad-platform
                     mention and no >=20-analyte mention

    We then test whether broad-platform studies report HIGHER amyloid-detection
    AUCs than targeted ones (Mann-Whitney U, two-sided).  The hypothesis predicts
    they do not.
    """  # noqa: D401
    d = a[a.endpoint == "amyloid_status"].copy()
    if d.empty:
        return pd.DataFrame(), {"note": "no amyloid-endpoint extractions"}
    big = (d.dropna(subset=["panel_size"]).groupby("source_id").panel_size.max()
           .rename("max_panel_size"))
    g = (d.groupby("source_id")
           .agg(max_auc=("auc", "max"), median_auc=("auc", "median"),
                n_auc=("auc", "size"), broad=("broad_platform", "any"),
                targeted=("targeted_platform", "any"),
                cu=("pop_cognitively_unimpaired", "any"),
                kind=("kind", "first"))
           .join(big))
    g["explicit_big_panel"] = g.max_panel_size.fillna(0) >= 20
    g["breadth"] = np.where(g.broad | g.explicit_big_panel, "broad",
                            np.where(g.targeted, "targeted", "unclassified"))
    g = g.reset_index()
    g.to_csv(os.path.join(OUT, "document_level.csv"), index=False)

    stats = {}
    for pop_label, sub in [("all", g), ("CU_mentioned", g[g.cu])]:
        b = sub.loc[sub.breadth == "broad", "median_auc"].dropna()
        t = sub.loc[sub.breadth == "targeted", "median_auc"].dropna()
        if len(b) >= 5 and len(t) >= 5:
            u, pu = st.mannwhitneyu(b, t, alternative="two-sided")
            # rank-biserial effect size
            rbc = 2 * u / (len(b) * len(t)) - 1
            stats[pop_label] = dict(
                n_broad=int(len(b)), n_targeted=int(len(t)),
                n_auc_per_doc_broad=float(sub.loc[sub.breadth == "broad", "n_auc"].median()),
                n_auc_per_doc_targeted=float(sub.loc[sub.breadth == "targeted", "n_auc"].median()),
                median_broad=float(b.median()), median_targeted=float(t.median()),
                mean_broad=float(b.mean()), mean_targeted=float(t.mean()),
                delta_median=float(b.median() - t.median()),
                mannwhitney_U=float(u), p=float(pu), rank_biserial=float(rbc))
        else:
            stats[pop_label] = dict(n_broad=int(len(b)), n_targeted=int(len(t)),
                                    note="too few documents in one arm")
    return g, stats


def audit_sample(a: pd.DataFrame, n=40, seed=C.SEED) -> pd.DataFrame:
    """Random sample of automated extractions for manual precision auditing."""
    s = a.dropna(subset=["panel_size"]).sample(min(n, len(a)), random_state=seed)
    s.to_csv(os.path.join(OUT, "passA_audit_sample.csv"), index=False)
    return s


if __name__ == "__main__":
    C.set_seed()
    print("--- PASS A: automated extraction (548 abstracts + 44 PDFs) ---")
    A = pass_a()
    print(f"  {len(A)} AUC mentions extracted "
          f"({(A.kind=='abstract').sum()} from abstracts, "
          f"{(A.kind=='fulltext').sum()} from full text)")
    print(f"  with a nearby panel-size mention: {A.panel_size.notna().sum()}")
    print(A.groupby(["kind", "endpoint"]).auc.agg(["count", "median"]).to_string())

    print("\n--- PASS B: curated within-study evidence ---")
    B = pass_b()
    print(f"  {len(B)} curated rows across {B.study.nunique()} studies")

    W = within_study_ratios(B)
    print("\n--- within-study retained-performance ratios "
          "(best compact panel k<=5 vs large panel k>=20) ---")
    print(W[["study", "n", "population", "k_small", "auc_small", "k_large", "auc_large",
             "auc_covariates", "ratio_raw", "ratio_excess", "ratio_anchored"]]
          .round(3).to_string(index=False))

    reg = {f"maxdist{md}": panel_size_regression(A, md) for md in (400, 800, 1500)}
    print("\n--- AUC ~ log10(panel size), by extraction-proximity threshold ---")
    for md, r in reg.items():
        for k, v in r.items():
            print(f"  [{md}] {k}: {v}")

    G, dl = document_level(A)
    print("\n--- document-level: broad vs targeted platforms, amyloid endpoint ---")
    print(G.breadth.value_counts().to_string())
    for k, v in dl.items():
        print(f"  {k}: {v}")

    audit_sample(A)
    C.dump_json({"env": C.env_report(), "n_extractions": int(len(A)),
                 "regression": reg, "document_level": dl,
                 "n_curated_rows": int(len(B)),
                 "n_within_study": int(len(W)),
                 "llm_available": False,
                 "note": "no LLM API key in environment; PASS A is rule-based regex"},
                os.path.join(OUT, "d3_meta.json"))
    print("\nD3 complete ->", OUT)
