"""
D1 - Empirical panel-size saturation curve on real plasma proteomic data.

Dataset : GSE275392 (Philippi & Castellano 2024) - SomaScan 1,305 plasma proteins
          x 53 non-demented older adults with amyloid status.
Primary : APOE e3/e3 stratum, n=36 (18 amyloid+ / 18 amyloid-).  Pre-registered in
          planning.md S5 because APOE genotype is otherwise perfectly confounded
          with amyloid status in this cohort.
Secondary: full cohort n=53, explicitly labelled confounded.

Question: how does out-of-sample AUC for amyloid positivity behave as a function of
          panel size k, when features are selected STRICTLY inside training folds?

Design (all pre-registered):
  * repeated stratified K-fold CV (N_REPEATS x N_FOLDS), out-of-fold (OOF) scores
    pooled within each repeat, AUC computed per repeat -> distribution over repeats
  * feature ranking (absolute two-sample t statistic) fitted on the TRAINING fold only
  * standardisation fitted on the training fold only
  * models: L2 logistic regression, elastic net logistic, random forest
  * k in {1,2,3,5,10,20,50,100,300,1305}
  * baselines: covariates (age, sex[, APOE4 dose]); a-priori AD protein panel;
    the Bio-Hermes top-8 inflammation panel (external panel replication)
  * label-permutation null (N_PERM repeats) to establish what AUC is achievable
    by chance at this sample size with this selection procedure
  * stratified bootstrap CIs on pooled OOF scores; paired DeLong vs the full panel

Outputs: results/d1/*.csv, results/d1/*.json
"""
from __future__ import annotations

import os
import sys
import time
import warnings

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

OUT = os.path.join(C.REPO, "results", "d1")
os.makedirs(OUT, exist_ok=True)

K_GRID = [1, 2, 3, 5, 10, 20, 50, 100, 300, 1305]
N_REPEATS = 20      # repeated CV
N_FOLDS = 5
N_PERM = 200        # label permutations for the null
N_BOOT = 2000

# A-priori compact panels, defined BEFORE looking at any performance number.
# Only proteins actually present on this 1,305-plex are listed.
APRIORI_AD = ["GFAP", "MAPT", "APP", "APOE", "CLU", "CST3", "CHI3L1", "SMOC1"]
# Bio-Hermes (IJMS 2026) reported these as its top recurring importances from a
# 295-protein discovery panel; C1QTNF5 is absent from this assay.
BIOHERMES8 = ["SERPINA1", "C3", "CRP", "APOE", "CFH", "VTN", "PON1"]
# The subset of the hypothesis' named core-4 that this platform can measure at all.
CORE_AVAILABLE = ["GFAP", "MAPT", "APP"]   # no p-tau217 / no NfL / no Ab42:40 on SomaScan


# --------------------------------------------------------------------------- #
def make_model(kind: str, seed: int):
    """Model zoo. Hyperparameters are fixed a priori (no tuning at n=36)."""
    # sklearn >=1.8 deprecates `penalty=`; regularisation type is set via l1_ratio.
    if kind == "l2":
        return LogisticRegression(l1_ratio=0.0, C=1.0, max_iter=5000, solver="lbfgs")
    if kind == "enet":
        return LogisticRegression(l1_ratio=0.5, C=1.0, solver="saga",
                                  max_iter=20000, tol=1e-3)
    if kind == "rf":
        # n_jobs=1: at n=36 the joblib dispatch overhead exceeds the fit cost
        return RandomForestClassifier(n_estimators=300, min_samples_leaf=2,
                                      max_features="sqrt", random_state=seed, n_jobs=1)
    raise ValueError(kind)


def rank_features(Xtr: np.ndarray, ytr: np.ndarray) -> np.ndarray:
    """
    Rank feature columns by |two-sample t| computed on TRAINING data only.

    Welch t is used (unequal variances); columns with zero variance in either
    class get t=0 and sort last.  Returns indices in descending |t| order.
    """
    a, b = Xtr[ytr == 1], Xtr[ytr == 0]
    ma, mb = a.mean(0), b.mean(0)
    va, vb = a.var(0, ddof=1), b.var(0, ddof=1)
    denom = np.sqrt(va / len(a) + vb / len(b))
    t = np.divide(ma - mb, denom, out=np.zeros_like(ma), where=denom > 1e-12)
    return np.argsort(-np.abs(t), kind="stable")


def oof_scores(X: np.ndarray, y: np.ndarray, k: int, kind: str, seed: int,
               n_folds: int = N_FOLDS, fixed_cols: np.ndarray | None = None,
               extra: np.ndarray | None = None) -> np.ndarray:
    """
    One repeat of stratified K-fold CV; returns pooled out-of-fold decision scores.

    Feature selection (top-k by |t|) and standardisation are fitted inside each
    training fold - never on the held-out fold.  `fixed_cols` bypasses selection
    (used for a-priori panels).  `extra` supplies always-included covariate
    columns (already numeric) appended after the selected proteins.
    """
    skf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=seed)
    s = np.full(len(y), np.nan)
    for tr, te in skf.split(X, y):
        if fixed_cols is not None:
            cols = fixed_cols
        elif X.shape[1] == 0:
            cols = np.array([], dtype=int)
        else:
            cols = rank_features(X[tr], y[tr])[:k]

        Xtr, Xte = X[np.ix_(tr, cols)], X[np.ix_(te, cols)]
        if extra is not None:
            Xtr = np.hstack([Xtr, extra[tr]])
            Xte = np.hstack([Xte, extra[te]])
        if Xtr.shape[1] == 0:
            s[te] = 0.0
            continue

        sc = StandardScaler().fit(Xtr)
        Xtr, Xte = sc.transform(Xtr), sc.transform(Xte)
        m = make_model(kind, seed).fit(Xtr, y[tr])
        s[te] = (m.predict_proba(Xte)[:, 1] if kind == "rf"
                 else m.decision_function(Xte))
    assert not np.isnan(s).any(), "some samples never appeared in a test fold"
    return s


def repeated_cv(X, y, k, kind, n_repeats=N_REPEATS, fixed_cols=None, extra=None,
                base_seed=C.SEED):
    """Run `n_repeats` independent repeated-CV passes. Returns (aucs, mean_oof_scores)."""
    aucs, S = [], []
    for r in range(n_repeats):
        s = oof_scores(X, y, k, kind, seed=base_seed + r,
                       fixed_cols=fixed_cols, extra=extra)
        # rank-normalise before averaging: decision-function scales differ per repeat
        S.append((np.argsort(np.argsort(s)) + 0.5) / len(s))
        aucs.append(C.auc(y, s))
    return np.array(aucs), np.mean(S, axis=0)


# --------------------------------------------------------------------------- #
def run_stratum(stratum: str) -> dict:
    t0 = time.time()
    coh = C.load_gse275392(stratum)
    X = coh.X.values.astype(float)
    y = coh.y
    names = np.array(coh.X.columns)
    print(f"\n=== D1 stratum={stratum} :: {coh} :: {coh.note}")

    # covariates: APOE4 dose is only meaningful (and non-degenerate) in the full cohort
    cov_cols = ["age", "sex_M"] + (["apoe4_dose"] if stratum == "full" else [])
    COV = coh.cov[cov_cols].values.astype(float)

    rows, oof_store = [], {}

    # ---- baseline 0: covariates only -------------------------------------- #
    a_cov, s_cov = repeated_cv(np.zeros((len(y), 0)), y, 0, "l2", extra=COV)
    oof_store["covariates"] = s_cov
    rows.append(dict(model="l2", panel="covariates_only", k=len(cov_cols),
                     auc_mean=a_cov.mean(), auc_sd=a_cov.std(ddof=1),
                     auc_pooled=C.auc(y, s_cov)))
    print(f"  covariates({'+'.join(cov_cols)}): AUC {a_cov.mean():.3f} +- {a_cov.std(ddof=1):.3f}")

    # ---- a-priori fixed panels -------------------------------------------- #
    for pname, plist in [("apriori_AD8", APRIORI_AD),
                         ("biohermes7", BIOHERMES8),
                         ("core_available3", CORE_AVAILABLE)]:
        cols = np.array([i for i, n in enumerate(names) if n in plist])
        assert len(cols) == len(plist), f"{pname}: missing proteins on this platform"
        for kind in ["l2", "rf"]:
            a, s = repeated_cv(X, y, len(cols), kind, fixed_cols=cols)
            if kind == "l2":
                oof_store[pname] = s
            rows.append(dict(model=kind, panel=pname, k=len(cols),
                             auc_mean=a.mean(), auc_sd=a.std(ddof=1),
                             auc_pooled=C.auc(y, s)))
        print(f"  {pname} (k={len(cols)}): l2 AUC "
              f"{rows[-2]['auc_mean']:.3f}, rf AUC {rows[-1]['auc_mean']:.3f}")

    # ---- the saturation curve --------------------------------------------- #
    for kind in ["l2", "enet", "rf"]:
        for k in K_GRID:
            a, s = repeated_cv(X, y, k, kind)
            rows.append(dict(model=kind, panel=f"topk", k=k,
                             auc_mean=a.mean(), auc_sd=a.std(ddof=1),
                             auc_pooled=C.auc(y, s),
                             auc_q025=np.percentile(a, 2.5),
                             auc_q975=np.percentile(a, 97.5)))
            oof_store[f"topk_{kind}_{k}"] = s
            print(f"  {kind:5s} k={k:5d}: AUC {a.mean():.3f} +- {a.std(ddof=1):.3f} "
                  f"(pooled {C.auc(y, s):.3f})  [{time.time()-t0:.0f}s]")

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, f"curve_{stratum}.csv"), index=False)

    # ---- permutation null (primary model only, to bound compute) ---------- #
    print("  permutation null (l2)...")
    rng = np.random.default_rng(C.SEED)
    perm = {k: [] for k in K_GRID}
    for p in range(N_PERM):
        yp = rng.permutation(y)
        for k in K_GRID:
            perm[k].append(C.auc(yp, oof_scores(X, yp, k, "l2", seed=C.SEED + 1000 + p)))
    pn = pd.DataFrame({"k": list(perm), **{
        "perm_mean": [np.mean(perm[k]) for k in perm],
        "perm_sd": [np.std(perm[k], ddof=1) for k in perm],
        "perm_q95": [np.percentile(perm[k], 95) for k in perm],
        "perm_q975": [np.percentile(perm[k], 97.5) for k in perm]}})
    obs = df[(df.model == "l2") & (df.panel == "topk")].set_index("k")["auc_mean"]
    pn["observed_auc"] = pn.k.map(obs)
    pn["p_perm"] = [(1 + np.sum(np.array(perm[k]) >= obs[k])) / (1 + N_PERM) for k in pn.k]
    pn.to_csv(os.path.join(OUT, f"permutation_{stratum}.csv"), index=False)
    print(pn.to_string(index=False))

    # ---- bootstrap CIs + DeLong vs full panel + retained-performance ratios #
    ref_key = f"topk_l2_{K_GRID[-1]}"          # the "20+ protein" comparator: full panel
    best_key = max([f"topk_l2_{k}" for k in K_GRID], key=lambda kk: C.auc(y, oof_store[kk]))
    auc_cov = C.auc(y, oof_store["covariates"])
    comp = []
    for key, s in oof_store.items():
        lo, hi, _ = C.boot_auc_ci(y, s, n_boot=N_BOOT)
        a1, a2, z, p = C.delong_test(y, s, oof_store[ref_key])
        r = C.ratios(a1, C.auc(y, oof_store[ref_key]), auc_cov)
        pr, plo, phi, _ = C.boot_paired_ratio(y, s, oof_store[ref_key], "excess",
                                              n_boot=N_BOOT)
        sens, spec, _ = C.youden_operating_point(y, s)
        ppv20, npv20 = C.ppv_npv_at_prevalence(sens, spec, 0.20)
        comp.append(dict(panel=key, auc=a1, ci_lo=lo, ci_hi=hi,
                         delong_z_vs_full=z, delong_p_vs_full=p, **r,
                         ratio_excess_lo=plo, ratio_excess_hi=phi,
                         sens_youden=sens, spec_youden=spec,
                         ppv_at_20pct=ppv20, npv_at_20pct=npv20))
    cdf = pd.DataFrame(comp).sort_values("auc", ascending=False)
    cdf.to_csv(os.path.join(OUT, f"panels_{stratum}.csv"), index=False)
    print(cdf.head(12).to_string(index=False))

    # ---- selection stability: which proteins get picked, how often? -------- #
    counts = {}
    for r in range(N_REPEATS):
        skf = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=C.SEED + r)
        for tr, _ in skf.split(X, y):
            for i in rank_features(X[tr], y[tr])[:10]:
                counts[names[i]] = counts.get(names[i], 0) + 1
    total = N_REPEATS * N_FOLDS
    stab = (pd.DataFrame({"protein": list(counts), "n_selected_top10": list(counts.values())})
            .assign(selection_frequency=lambda d: d.n_selected_top10 / total)
            .sort_values("selection_frequency", ascending=False))
    stab.to_csv(os.path.join(OUT, f"stability_{stratum}.csv"), index=False)
    print("\n  top-10 selection stability (fraction of the "
          f"{total} training folds in which a protein was in the top 10):")
    print(stab.head(12).to_string(index=False))

    np.savez_compressed(os.path.join(OUT, f"oof_{stratum}.npz"), y=y,
                        **{k: v for k, v in oof_store.items()})
    return dict(stratum=stratum, n=int(len(y)), n_pos=int(y.sum()),
                auc_cov=auc_cov, ref_key=ref_key, best_key=best_key,
                runtime_s=time.time() - t0)


if __name__ == "__main__":
    C.set_seed()
    meta = {"env": C.env_report(), "config": dict(
        K_GRID=K_GRID, N_REPEATS=N_REPEATS, N_FOLDS=N_FOLDS,
        N_PERM=N_PERM, N_BOOT=N_BOOT, APRIORI_AD=APRIORI_AD,
        BIOHERMES8=BIOHERMES8, CORE_AVAILABLE=CORE_AVAILABLE)}
    meta["runs"] = [run_stratum(s) for s in ["apoe33", "full"]]
    C.dump_json(meta, os.path.join(OUT, "d1_meta.json"))
    print("\nD1 complete ->", OUT)
