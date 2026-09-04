"""
D2 - Literature-calibrated boundary analysis: when does a compact core panel fail
     to retain >=90% of a large panel's performance?

No public individual-level dataset contains the hypothesis' NAMED analytes
(p-tau217, Ab42/40, GFAP, NfL, MTBR-tau243) together with a broad panel and
amyloid status in a cognitively unimpaired cohort (planning.md S1).  D2 therefore
builds a generative model whose marginal effect sizes are pinned to published
AUCs, and asks the question under KNOWN ground truth:

    under what conditions - how many extra informative proteins, of what effect
    size, at what correlation with the core markers, at what training-set size -
    does the "compact panel retains >=90%" claim hold or break?

This is a sensitivity / boundary analysis, NOT empirical evidence.  Every
calibration constant is traceable to a source listed in CALIBRATION below.

Calibration anchors
-------------------
  p-tau217      AUC 0.904  ADNI n=1317   datasets/literature_reference/.../S7 (also S4: 0.902)
  Ab42/40       AUC 0.831  ADNI n=1317   same table (S4 LO stratum: 0.811)
  GFAP          AUC 0.646  AMASS n=122   Trelle et al. 2026 Mol Neurodegener 21:31
  NfL           AUC 0.582  AMASS n=122   Trelle et al. 2026 (worse than covariates)
  MTBR-tau243   AUC 0.60   assumed       Horie 2023: tau-TANGLE staging marker, not an
                                         early amyloid marker -> weak for this endpoint
  age+sex+APOE4 AUC 0.773  AMASS n=122   Trelle et al. 2026 (the mandatory baseline)
  n informative / panel size: 8/120       Knight ADRC NULISAseq, Cruchaga 2025
                             ~5/130       Trelle 2026
  extras' AUC range 0.60-0.70            Bio-Hermes top proteins (SERPINA1/C3/CRP/...)

Effect sizes are converted from AUC by the equal-variance binormal identity
    d = sqrt(2) * Phi^-1(AUC).

Outputs: results/d2/*.csv, results/d2/*.json
"""
from __future__ import annotations

import os
import sys
import time
import warnings

import numpy as np
import pandas as pd
import scipy.stats as st
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

warnings.filterwarnings("ignore")
OUT = os.path.join(C.REPO, "results", "d2")
os.makedirs(OUT, exist_ok=True)

# --------------------------------------------------------------------------- #
# Calibration
# --------------------------------------------------------------------------- #
def auc_to_d(a: float) -> float:
    """Binormal equal-variance conversion AUC -> Cohen's d."""
    return float(np.sqrt(2) * st.norm.ppf(a))


def d_to_auc(d: float) -> float:
    return float(st.norm.cdf(d / np.sqrt(2)))


CALIBRATION = {
    "ptau217":    dict(auc=0.904, src="ADNI n=1317 (S7/S4 supplementary tables)"),
    "ab42_40":    dict(auc=0.831, src="ADNI n=1317 (S7)"),
    "gfap":       dict(auc=0.646, src="Trelle 2026 AMASS n=122"),
    "nfl":        dict(auc=0.582, src="Trelle 2026 AMASS n=122"),
    "mtbr243":    dict(auc=0.600, src="assumed; Horie 2023 = tau-tangle staging marker"),
    "covariates": dict(auc=0.773, src="Trelle 2026 AMASS, age+sex+APOE4"),
}
CORE = ["ptau217", "ab42_40", "gfap", "nfl", "mtbr243"]
CORE_D = np.array([auc_to_d(CALIBRATION[m]["auc"]) for m in CORE])

# Assumed BIOLOGICAL correlation structure among the core markers (all AD-biology
# driven, hence positively correlated once signs are aligned).  Observed
# correlations are attenuated by assay noise (PSI_ASSAY below).  Sweep S4 varies
# the core<->extra shared-factor loading.
CORE_R = np.array([
    # ptau217 ab42/40  gfap   nfl   mtbr
    [1.00,   0.45,   0.35,  0.25,  0.55],
    [0.45,   1.00,   0.25,  0.20,  0.30],
    [0.35,   0.25,   1.00,  0.45,  0.25],
    [0.25,   0.20,   0.45,  1.00,  0.20],
    [0.55,   0.30,   0.25,  0.20,  1.00],
])

# Fraction of each core marker's variance that is ASSAY noise - independent,
# and therefore not removable by measuring any number of other proteins.
# Immunoassay CVs of 5-15% against typical between-subject biological spread put
# this at roughly 0.20-0.35; 0.25 is used and swept in the sensitivity run.
PSI_ASSAY = 0.25

N_EXTRA_FACTORS = 8   # shared non-AD factors among the extras (age, BMI, renal, plate, ...)

PANELS = {
    "covariates_only": [],
    "core1_ptau217":   ["ptau217"],
    "core2_ptau217_ab": ["ptau217", "ab42_40"],
    "core3":           ["ptau217", "ab42_40", "gfap"],
    "core4_hypothesis": ["ptau217", "ab42_40", "gfap", "nfl"],
    "core5_hypothesis": ["ptau217", "ab42_40", "gfap", "nfl", "mtbr243"],
    "gfap_nfl_only":   ["gfap", "nfl"],
}


# --------------------------------------------------------------------------- #
# Generative model
# --------------------------------------------------------------------------- #
def simulate(n, prevalence, n_extra, n_informative, d_extra, rho_core_extra,
             rho_extra, rng, core_d=CORE_D, core_r=CORE_R, psi=PSI_ASSAY):
    """
    Draw a cohort under a factor model.

    Structure
    ---------
    y ~ Bernoulli(prevalence).

    Core markers   m_i = d_i * y + sqrt(1-psi) * s_i + sqrt(psi) * eta_i
        s ~ N(0, core_r)  : shared AD biology, partially measurable by other proteins
        eta ~ N(0, I)     : ASSAY noise, independent, NOT removable by any panel

    Extra proteins e_j = d_e * y * 1[j informative]
                       + a * (u_j . G) + b * H_{g(j)} + c * eps_j
        G           : whitened 5-dim basis of the core biology s
        u_j         : a RANDOM unit direction, different for every extra
        H           : N_EXTRA_FACTORS shared non-AD nuisance factors
        a=rho_core_extra, b=rho_extra, c=sqrt(1-a^2-b^2)

    Why the random direction u_j matters
    ------------------------------------
    An exchangeable correlation (every extra correlating identically with every
    core marker) is degenerate: the average of 120 such extras estimates the
    shared factor almost noiselessly, so a large panel can subtract the core
    markers' biological noise and drive AUC to 1.0.  That is a modelling
    artifact, not biology.  Random per-protein directions plus an irreducible
    assay-noise floor reproduce the empirically observed ceiling instead (see
    the realism check in __main__: broad panel without p-tau217 -> AUC ~0.80,
    matching Bio-Hermes' 295-protein result).

    Every marker is marginally unit-variance, so each core marker's marginal AUC
    equals Phi(d_i/sqrt2), i.e. its calibration target, by construction.
    """
    y = (rng.random(n) < prevalence).astype(int)

    # --- core markers ---
    Lc = np.linalg.cholesky(core_r)                 # rows are unit-norm loadings
    G = rng.standard_normal((n, 5))                 # whitened core-biology factors
    s = G @ Lc.T
    eta = rng.standard_normal((n, 5))
    M_core = np.outer(y, core_d) + np.sqrt(1 - psi) * s + np.sqrt(psi) * eta

    # --- extra proteins ---
    if n_extra > 0:
        a = float(np.clip(rho_core_extra, 0, 0.9))
        b = float(np.clip(rho_extra, 0, np.sqrt(max(1e-9, 1 - a ** 2))))
        c = float(np.sqrt(max(1e-9, 1 - a ** 2 - b ** 2)))
        U = rng.standard_normal((n_extra, 5))
        U /= np.linalg.norm(U, axis=1, keepdims=True)      # random unit directions
        H = rng.standard_normal((n, N_EXTRA_FACTORS))
        grp = rng.integers(0, N_EXTRA_FACTORS, n_extra)
        M_ex = (a * (G @ U.T) + b * H[:, grp] + c * rng.standard_normal((n, n_extra)))
        if n_informative > 0:
            k = min(n_informative, n_extra)
            M_ex[:, :k] += np.outer(y, np.full(k, d_extra))
    else:
        M_ex = np.zeros((n, 0))

    cols = CORE + [f"extra{i:04d}" for i in range(n_extra)]
    X = pd.DataFrame(np.hstack([M_core, M_ex]), columns=cols)

    # --- covariates ---
    p_e4 = _APOE_P0 * np.exp(_APOE_B * y) / (1 - _APOE_P0 + _APOE_P0 * np.exp(_APOE_B * y))
    apoe4 = rng.binomial(2, np.clip(p_e4, 0, 1))
    age = 70 + 7 * rng.standard_normal(n) + _AGE_D * 7 * y
    sex = rng.integers(0, 2, n)
    COV = pd.DataFrame({"apoe4_dose": apoe4, "age": age, "sex_M": sex})
    return X, COV, y


# Solve the APOE4 effect so covariate-only AUC == calibrated 0.773 -------------
_APOE_P0 = 0.14      # allele freq of e4 in amyloid-negative CU (population value)
_AGE_D = 0.35        # amyloid+ are modestly older
_APOE_B = 1.4        # placeholder; replaced by the calibration below


def _calibrate_apoe(target=CALIBRATION["covariates"]["auc"], n=200_000, seed=7):
    """Numerically solve the APOE4 log-odds giving the published covariate AUC."""
    global _APOE_B
    lo, hi = 0.0, 6.0
    for _ in range(40):
        mid = (lo + hi) / 2
        _APOE_B = mid
        rng = np.random.default_rng(seed)
        _, COV, y = simulate(n, 0.25, 0, 0, 0.0, 0.0, 0.0, rng)
        m = LogisticRegression(max_iter=2000).fit(
            StandardScaler().fit_transform(COV.values), y)
        a = C.auc(y, m.decision_function(StandardScaler().fit_transform(COV.values)))
        if a < target:
            lo = mid
        else:
            hi = mid
    _APOE_B = (lo + hi) / 2
    return _APOE_B, a


# --------------------------------------------------------------------------- #
# Evaluation
# --------------------------------------------------------------------------- #
def fit_eval(Xtr, Ctr, ytr, Xte, Cte, yte, cols, l1_ratio=0.0, Cparam=1.0):
    """Fit an L2 (or elastic-net) logistic model on `cols` + covariates; test AUC."""
    A = np.hstack([Xtr[cols].values, Ctr.values]) if cols else Ctr.values
    B = np.hstack([Xte[cols].values, Cte.values]) if cols else Cte.values
    sc = StandardScaler().fit(A)
    solver = "lbfgs" if l1_ratio == 0.0 else "saga"
    m = LogisticRegression(l1_ratio=l1_ratio, C=Cparam, solver=solver,
                           max_iter=8000, tol=1e-3).fit(sc.transform(A), ytr)
    return C.auc(yte, m.decision_function(sc.transform(B)))


def oracle_auc(d_vec, R):
    """Bayes-optimal AUC for the binormal model: Phi( sqrt(d' R^-1 d) / sqrt2 )."""
    m = float(d_vec @ np.linalg.solve(R, d_vec))
    return float(st.norm.cdf(np.sqrt(max(m, 0)) / np.sqrt(2)))


def one_condition(rep, n_train, prevalence, n_extra, n_informative, d_extra,
                  rho_core_extra, rho_extra, n_test=8000, seed0=1000):
    """One simulation replicate: returns a dict of AUCs for every panel."""
    rng = np.random.default_rng(seed0 + rep)
    Xtr, Ctr, ytr = simulate(n_train, prevalence, n_extra, n_informative,
                             d_extra, rho_core_extra, rho_extra, rng)
    Xte, Cte, yte = simulate(n_test, prevalence, n_extra, n_informative,
                             d_extra, rho_core_extra, rho_extra, rng)
    if ytr.sum() < 5 or (1 - ytr).sum() < 5:
        return None

    res = {}
    for name, cols in PANELS.items():
        res[name] = fit_eval(Xtr, Ctr, ytr, Xte, Cte, yte, cols)
    allc = list(Xtr.columns)
    res["full_panel"] = fit_eval(Xtr, Ctr, ytr, Xte, Cte, yte, allc)
    res["full_panel_enet"] = fit_eval(Xtr, Ctr, ytr, Xte, Cte, yte, allc, l1_ratio=0.5)
    # the SomaScan / Bio-Hermes situation: broad panels that CANNOT measure
    # phospho-epitopes.  `discovery_only` is the realistic aptamer/PEA case
    # (extras + GFAP + NfL, no p-tau217, no Ab42/40 ratio); `extras_only` is the
    # pure discovery-proteomics case used as the realism-calibration target
    # (Bio-Hermes: 295 proteins -> AUC 0.79-0.81).
    res["broad_no_ptau"] = fit_eval(Xtr, Ctr, ytr, Xte, Cte, yte,
                                    [c for c in allc if c not in ("ptau217", "mtbr243")])
    res["discovery_only"] = fit_eval(Xtr, Ctr, ytr, Xte, Cte, yte,
                                     [c for c in allc if c.startswith("extra")]
                                     + ["gfap", "nfl"])
    res["extras_only"] = fit_eval(Xtr, Ctr, ytr, Xte, Cte, yte,
                                  [c for c in allc if c.startswith("extra")])
    # ground-truth reference: the Bayes-optimal AUC obtainable from the 5 core
    # markers alone, under their OBSERVED covariance (biology attenuated by assay noise)
    Rc = (1 - PSI_ASSAY) * CORE_R + PSI_ASSAY * np.eye(5)
    res["oracle_core5"] = oracle_auc(CORE_D, Rc)
    res.update(rep=rep, n_train=n_train, prevalence=prevalence, n_extra=n_extra,
               n_informative=n_informative, d_extra=d_extra,
               rho_core_extra=rho_core_extra, rho_extra=rho_extra)
    return res


def run_sweep(name, conditions, n_reps, tag=""):
    rows = []
    t0 = time.time()
    for ci, cond in enumerate(conditions):
        for r in range(n_reps):
            out = one_condition(r, **cond)
            if out:
                out["condition"] = ci
                rows.append(out)
        print(f"  [{name}] cond {ci+1}/{len(conditions)} {cond} "
              f"({time.time()-t0:.0f}s)", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, f"raw_{name}.csv"), index=False)
    return df


def summarise(df: pd.DataFrame, group_cols) -> pd.DataFrame:
    """Aggregate replicate AUCs and compute the three retained-performance ratios."""
    panel_cols = [c for c in df.columns
                  if c in list(PANELS) + ["full_panel", "full_panel_enet",
                                          "broad_no_ptau", "discovery_only",
                                          "extras_only", "oracle_core5"]]
    g = df.groupby(group_cols)[panel_cols].agg(["mean", "std"])
    g.columns = [f"{a}_{b}" for a, b in g.columns]
    g = g.reset_index()
    ref, cov = g["full_panel_mean"], g["covariates_only_mean"]
    for p in ["core1_ptau217", "core2_ptau217_ab", "core4_hypothesis",
              "core5_hypothesis", "gfap_nfl_only", "broad_no_ptau",
              "discovery_only", "extras_only"]:
        s = g[f"{p}_mean"]
        g[f"{p}__ratio_raw"] = s / ref
        g[f"{p}__ratio_excess"] = (s - 0.5) / (ref - 0.5)
        g[f"{p}__ratio_anchored"] = (s - cov) / (ref - cov)
    return g


# --------------------------------------------------------------------------- #
# Analytic boundary (closed form, no simulation)
# --------------------------------------------------------------------------- #
def analytic_boundary():
    """
    Under independence, discriminant information adds: D^2_total = sum d_i^2 and
    AUC = Phi(D/sqrt2).  For a core panel with D_core and m extras each of effect
    d_e, the excess ratio is
        (Phi(D_core/sqrt2) - 0.5) / (Phi(sqrt(D_core^2 + m d_e^2)/sqrt2) - 0.5).
    We solve for the smallest m at which each ratio definition drops below 0.90.
    This gives the *most favourable possible* case for large panels (fully
    independent extras); correlation can only reduce their contribution.
    """
    rows = []
    # observed core covariance = biological correlation attenuated by assay noise
    Rc = (1 - PSI_ASSAY) * CORE_R + PSI_ASSAY * np.eye(5)
    d_core_full = float(np.sqrt(CORE_D @ np.linalg.solve(Rc, CORE_D)))
    variants = {"core5_hypothesis": d_core_full,
                "core1_ptau217": auc_to_d(CALIBRATION["ptau217"]["auc"]),
                "core2_ptau217_ab": float(np.sqrt(
                    CORE_D[:2] @ np.linalg.solve(Rc[:2, :2], CORE_D[:2])))}
    for vname, Dc in variants.items():
        for d_e in [0.2, 0.3, 0.4, 0.5, 0.75, 1.0]:
            row = dict(core=vname, d_core=Dc, auc_core=d_to_auc(Dc),
                       d_extra=d_e, auc_extra=d_to_auc(d_e))
            for rname in ["raw", "excess"]:
                m_break = np.nan
                for m in range(0, 5001):
                    Df = np.sqrt(Dc ** 2 + m * d_e ** 2)
                    a_c, a_f = d_to_auc(Dc), d_to_auc(Df)
                    r = (a_c / a_f if rname == "raw"
                         else (a_c - 0.5) / (a_f - 0.5))
                    if r < 0.90:
                        m_break = m
                        break
                row[f"m_extras_to_break_{rname}"] = m_break
            rows.append(row)
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    C.set_seed()
    t_all = time.time()
    b, a_ach = _calibrate_apoe()
    print(f"APOE4 log-odds calibrated to {b:.3f} -> covariate-only AUC {a_ach:.3f} "
          f"(target {CALIBRATION['covariates']['auc']})")
    print("Core effect sizes (Cohen's d):",
          {m: round(float(dd), 3) for m, dd in zip(CORE, CORE_D)})

    # BASE condition, chosen so the simulation reproduces THREE independent
    # published quantities simultaneously (verified by the realism check below):
    #   (i)   each core marker's marginal AUC == its calibration target
    #   (ii)  covariate-only AUC == 0.773            (Trelle 2026 AMASS)
    #   (iii) a large discovery panel without phospho-epitopes == ~0.80
    #         (Bio-Hermes, 295 proteins, IJMS 2026)
    # n_informative=8/120 matches the Knight ADRC rate (8/120 associated with
    # amyloid PET); d_extra=0.35 corresponds to a per-protein AUC of 0.598.
    BASE = dict(prevalence=0.20, n_extra=120, n_informative=8, d_extra=0.35,
                rho_core_extra=0.15, rho_extra=0.15)
    REPS, REPS_SWEEP = 200, 100

    # ---- realism check: does the calibrated model reproduce the published
    #      broad-panel ceiling it was NOT fitted to marker-by-marker? --------- #
    chk = [one_condition(r, n_train=1000, **BASE, n_test=8000) for r in range(30)]
    realism = {k: float(np.mean([c[k] for c in chk])) for k in
               ["covariates_only", "core1_ptau217", "core2_ptau217_ab",
                "core5_hypothesis", "extras_only", "discovery_only", "full_panel"]}
    print("\n--- realism check at the BASE condition (30 reps, n_train=1000) ---")
    for k, v in realism.items():
        print(f"    {k:22s} {v:.3f}")
    print("    targets: covariates_only~0.773 (Trelle), core1_ptau217~0.904 (ADNI, "
          "here jointly with covariates so slightly higher), extras_only~0.80 (Bio-Hermes)")

    sweeps = {}
    # S1 - training-set size
    sweeps["S1_n"] = ([{**BASE, "n_train": n} for n in [100, 250, 500, 1000, 2500, 5000]],
                      REPS, ["n_train"])
    # S2 - how many extra proteins are truly informative
    sweeps["S2_ninf"] = ([{**BASE, "n_train": 1000, "n_informative": k}
                          for k in [0, 2, 5, 8, 15, 30, 60, 120]], REPS_SWEEP, ["n_informative"])
    # S3 - how strong those extras are
    sweeps["S3_deffect"] = ([{**BASE, "n_train": 1000, "d_extra": d}
                             for d in [0.0, 0.2, 0.35, 0.45, 0.6, 0.8, 1.0, 1.5]],
                            REPS_SWEEP, ["d_extra"])
    # S4 - correlation between extras and the core markers (redundancy)
    sweeps["S4_rho"] = ([{**BASE, "n_train": 1000, "rho_core_extra": r}
                         for r in [0.0, 0.15, 0.30, 0.50, 0.70]], REPS_SWEEP,
                        ["rho_core_extra"])
    # S5 - panel size at fixed informative fraction (~6.7%, the Knight ADRC rate)
    sweeps["S5_panelsize"] = ([{**BASE, "n_train": 1000, "n_extra": p,
                                "n_informative": max(1, int(round(0.067 * p)))}
                               for p in [20, 60, 120, 300, 600, 1200]], REPS_SWEEP,
                              ["n_extra"])
    # S6 - prevalence (drives PPV/NPV, not AUC)
    sweeps["S6_prev"] = ([{**BASE, "n_train": 1000, "prevalence": p}
                          for p in [0.10, 0.20, 0.30, 0.50]], REPS_SWEEP, ["prevalence"])
    # S7 - adversarial: MANY strong independent extras (best case for big panels)
    sweeps["S7_adversarial"] = ([{**BASE, "n_train": 2500, "n_informative": k,
                                  "d_extra": 0.8, "rho_core_extra": 0.0, "rho_extra": 0.0}
                                 for k in [0, 5, 10, 20, 40, 80, 120]], REPS_SWEEP,
                                ["n_informative"])

    summaries = {}
    for name, (conds, reps, gcols) in sweeps.items():
        print(f"\n--- sweep {name} ({len(conds)} conditions x {reps} reps) ---")
        df = run_sweep(name, conds, reps)
        s = summarise(df, gcols)
        s.to_csv(os.path.join(OUT, f"summary_{name}.csv"), index=False)
        summaries[name] = s
        show = gcols + ["covariates_only_mean", "core1_ptau217_mean",
                        "core4_hypothesis_mean", "core5_hypothesis_mean",
                        "full_panel_mean", "broad_no_ptau_mean",
                        "core5_hypothesis__ratio_raw", "core5_hypothesis__ratio_excess",
                        "core5_hypothesis__ratio_anchored"]
        print(s[show].round(3).to_string(index=False))

    ab = analytic_boundary()
    ab.to_csv(os.path.join(OUT, "analytic_boundary.csv"), index=False)
    print("\n--- analytic boundary: number of INDEPENDENT extra proteins of effect "
          "d_extra needed to push the retained-performance ratio below 0.90 ---")
    print(ab.round(3).to_string(index=False))

    C.dump_json({"env": C.env_report(), "calibration": CALIBRATION,
                 "core_d": {m: float(d) for m, d in zip(CORE, CORE_D)},
                 "apoe_logodds": float(b), "covariate_auc_achieved": float(a_ach),
                 "base_condition": BASE, "reps": REPS, "reps_sweep": REPS_SWEEP,
                 "core_corr": CORE_R.tolist(), "psi_assay": PSI_ASSAY,
                 "n_extra_factors": N_EXTRA_FACTORS, "realism_check": realism,
                 "runtime_s": time.time() - t_all},
                os.path.join(OUT, "d2_meta.json"))
    print(f"\nD2 complete in {time.time()-t_all:.0f}s ->", OUT)
