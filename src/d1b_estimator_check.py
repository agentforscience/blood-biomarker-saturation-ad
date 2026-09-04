"""
D1b - Why the headline AUC depends on how you average the folds.

The phase-1 feasibility probe (`tools/validate_dataset.py`) reported AUC ~0.75 at
k=10 on the APOE33 stratum.  The full D1 experiment reports ~0.64 for the same
data, model class and selection procedure.  This module identifies why, testing
two candidate explanations.  The answer turned out to be the second one.

  CANDIDATE 1 - AUC estimator choice.  Measured below; it accounts for only
  +0.005 to +0.015 AUC.  Real but far too small to explain the gap.

  CANDIDATE 2 - single-split sampling noise.  THIS IS THE EXPLANATION.  Rerunning
  the probe's exact configuration over 200 random CV seeds gives mean 0.641
  (SD 0.065, range 0.487-0.804) at k=10.  The probe used `random_state=0`, which
  happens to land at the 94.5th percentile of that distribution.  Only 5.5% of
  seeds reach >=0.75.  At n=36 a single 5-fold CV estimate has a standard
  deviation of ~0.065 AUC, so any single-split number is uninterpretable.
  This is the concrete reason D1 uses 20x5 repeated CV plus a permutation null.

The estimator comparison itself remains worth reporting:

  per-fold averaging   (sklearn's cross_val_score(scoring='roc_auc'))
      computes an AUC inside each held-out fold and averages the five values.
      At n=36 each fold has ~7 subjects (3-4 per class), so each fold AUC takes
      one of a handful of discrete values and is enormously variable.  Averaging
      such AUCs is a biased estimator of the pooled discrimination, because AUC
      is a non-linear functional of the score distribution and small-sample AUCs
      are not centred on the population value.

  pooled out-of-fold   (used throughout D1)
      concatenates every held-out prediction and computes ONE AUC over all n
      subjects.  This is the standard recommendation for small samples
      (Airola et al. 2011; Forman & Scholz 2010).

This module quantifies the gap on the real data and, decisively, under a NULL
where the true AUC is known to be 0.5 by construction (permuted labels).  Any
estimator that does not return ~0.5 under the null is biased.

Outputs: results/d1/estimator_comparison.csv, results/d1/single_split_variability.csv,
         figures/fig7_estimator_bias.png
"""
from __future__ import annotations

import os
import sys
import warnings

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C          # noqa: E402
import d1_saturation as D1  # noqa: E402

warnings.filterwarnings("ignore")
OUT = os.path.join(C.REPO, "results", "d1")
FIG = os.path.join(C.REPO, "figures")

K_GRID = [1, 5, 10, 20, 100, 1305]
N_REPEATS = 50
N_NULL = 300


def per_fold_auc(X, y, k, kind, seed, n_folds=D1.N_FOLDS):
    """The sklearn-style estimator: AUC computed inside each fold, then averaged."""
    from sklearn.preprocessing import StandardScaler

    skf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=seed)
    vals = []
    for tr, te in skf.split(X, y):
        cols = D1.rank_features(X[tr], y[tr])[:k]
        sc = StandardScaler().fit(X[np.ix_(tr, cols)])
        m = D1.make_model(kind, seed).fit(sc.transform(X[np.ix_(tr, cols)]), y[tr])
        vals.append(C.auc(y[te], m.decision_function(sc.transform(X[np.ix_(te, cols)]))))
    return float(np.nanmean(vals)), vals


def single_split_variability(n_seeds=200, ks=(10, 20)):
    """
    Rerun the phase-1 probe's EXACT configuration (SelectKBest f_classif,
    StandardScaler, LogisticRegression(C=0.1), 5-fold, scoring='roc_auc') over
    `n_seeds` random CV splits, to show how much a single-split estimate moves.
    """
    from sklearn.feature_selection import SelectKBest, f_classif
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import cross_val_score
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler

    coh = C.load_gse275392("apoe33")
    X, y = coh.X.values.astype(float), coh.y
    rows = []
    for k in ks:
        vals = []
        for seed in range(n_seeds):
            pipe = Pipeline([("sel", SelectKBest(f_classif, k=k)),
                             ("sc", StandardScaler()),
                             ("lr", LogisticRegression(C=0.1, max_iter=5000))])
            vals.append(cross_val_score(
                pipe, X, y, cv=StratifiedKFold(5, shuffle=True, random_state=seed),
                scoring="roc_auc").mean())
        v = np.array(vals)
        rows.append(dict(k=k, n_seeds=n_seeds, mean=v.mean(), sd=v.std(ddof=1),
                         lo=v.min(), hi=v.max(), seed0_as_used_in_phase1=vals[0],
                         pct_seeds_ge_0p75=float((v >= 0.75).mean()),
                         percentile_of_seed0=float((v <= vals[0]).mean())))
        print(f"  probe config, k={k}: mean {v.mean():.3f} (sd {v.std(ddof=1):.3f}), "
              f"range [{v.min():.3f}, {v.max():.3f}]; seed 0 = {vals[0]:.3f} "
              f"({(v <= vals[0]).mean():.1%} percentile)", flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "single_split_variability.csv"), index=False)
    return df


def main():
    C.set_seed()
    print("--- single-split variability of the phase-1 probe configuration ---")
    single_split_variability()
    print("\n--- estimator comparison ---")
    coh = C.load_gse275392("apoe33")
    X, y = coh.X.values.astype(float), coh.y

    rows = []
    for k in K_GRID:
        pf = [per_fold_auc(X, y, k, "l2", C.SEED + r)[0] for r in range(N_REPEATS)]
        po = [C.auc(y, D1.oof_scores(X, y, k, "l2", seed=C.SEED + r))
              for r in range(N_REPEATS)]
        rows.append(dict(k=k, labels="observed",
                         per_fold_mean=np.mean(pf), per_fold_sd=np.std(pf, ddof=1),
                         pooled_mean=np.mean(po), pooled_sd=np.std(po, ddof=1),
                         gap=np.mean(pf) - np.mean(po)))
        print(f"  k={k:5d} observed : per-fold {np.mean(pf):.3f}  "
              f"pooled {np.mean(po):.3f}  gap {np.mean(pf)-np.mean(po):+.3f}", flush=True)

    # --- the decisive test: permuted labels, where the truth is AUC = 0.5 ------
    rng = np.random.default_rng(C.SEED)
    for k in K_GRID:
        pf, po = [], []
        for _ in range(N_NULL):
            yp = rng.permutation(y)
            pf.append(per_fold_auc(X, yp, k, "l2", C.SEED)[0])
            po.append(C.auc(yp, D1.oof_scores(X, yp, k, "l2", seed=C.SEED)))
        rows.append(dict(k=k, labels="permuted_null_true_auc_0.5",
                         per_fold_mean=np.mean(pf), per_fold_sd=np.std(pf, ddof=1),
                         pooled_mean=np.mean(po), pooled_sd=np.std(po, ddof=1),
                         gap=np.mean(pf) - np.mean(po)))
        print(f"  k={k:5d} NULL     : per-fold {np.mean(pf):.3f} (sd {np.std(pf,ddof=1):.3f})  "
              f"pooled {np.mean(po):.3f} (sd {np.std(po,ddof=1):.3f})  "
              f"-- truth is 0.500", flush=True)

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "estimator_comparison.csv"), index=False)

    plt = C.use_style()
    fig, ax = plt.subplots(1, 2, figsize=(10.5, 4.1))
    for a, lab, title in [(ax[0], "observed", "A. Observed labels"),
                          (ax[1], "permuted_null_true_auc_0.5",
                           "B. Permuted labels (true AUC = 0.5)")]:
        d = df[df.labels == lab].sort_values("k")
        a.plot(d.k, d.per_fold_mean, "o-", color=C.PAL["red"], lw=1.8, ms=4,
               label="per-fold averaged (sklearn default)")
        a.fill_between(d.k, d.per_fold_mean - d.per_fold_sd,
                       d.per_fold_mean + d.per_fold_sd, color=C.PAL["red"], alpha=0.12, lw=0)
        a.plot(d.k, d.pooled_mean, "s-", color=C.PAL["primary"], lw=1.8, ms=4,
               label="pooled out-of-fold (used in D1)")
        a.fill_between(d.k, d.pooled_mean - d.pooled_sd, d.pooled_mean + d.pooled_sd,
                       color=C.PAL["primary"], alpha=0.12, lw=0)
        a.set_xscale("log")
        a.set_xlabel("panel size $k$")
        a.set_title(title, fontsize=10)
        if lab != "observed":
            a.axhline(0.5, color="k", ls="--", lw=1.2)
            a.text(1.1, 0.505, "truth", fontsize=8)
    ax[0].set_ylabel("estimated AUC")
    ax[0].legend(fontsize=8, loc="upper left")
    fig.suptitle("At n=36, averaging per-fold AUCs inflates the estimate; "
                 "pooling out-of-fold predictions does not", fontsize=10.5, y=1.02)
    fig.savefig(os.path.join(FIG, "fig7_estimator_bias.png"))
    plt.close(fig)
    print("\nfigure -> figures/fig7_estimator_bias.png")
    return df


if __name__ == "__main__":
    main()
