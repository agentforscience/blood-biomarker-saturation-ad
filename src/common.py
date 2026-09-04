"""
Shared utilities for the panel-size benchmarking study.

Contains:
  * reproducible seeding and environment capture
  * GSE275392 loading / preprocessing
  * AUC machinery: fast Mann-Whitney AUC, DeLong variance + paired DeLong test,
    stratified bootstrap CIs
  * the three "fraction of performance retained" ratio definitions that the
    >=90% hypothesis can be read under (see planning.md S0.3 C2)
  * matplotlib defaults

Everything here is deterministic given a seed.
"""
from __future__ import annotations

import json
import os
import platform
import random
import sys
from dataclasses import dataclass

import numpy as np
import pandas as pd
import scipy.stats as st

# --------------------------------------------------------------------------- #
# Reproducibility
# --------------------------------------------------------------------------- #
SEED = 42
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def set_seed(seed: int = SEED) -> None:
    """Seed every stochastic source we use."""
    random.seed(seed)
    np.random.seed(seed)
    os.environ["PYTHONHASHSEED"] = str(seed)


def env_report() -> dict:
    """Capture the software/hardware environment for the reproducibility section."""
    import sklearn
    import scipy
    import matplotlib

    return {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scipy": scipy.__version__,
        "scikit_learn": sklearn.__version__,
        "matplotlib": matplotlib.__version__,
        "seed": SEED,
    }


# --------------------------------------------------------------------------- #
# Data loading
# --------------------------------------------------------------------------- #
@dataclass
class Cohort:
    """A design matrix + label vector + covariate frame for one analysis stratum."""

    name: str
    X: pd.DataFrame          # samples x proteins, log2 RFU, z-scored later inside folds
    y: np.ndarray            # 1 = amyloid positive
    cov: pd.DataFrame        # covariates (age, sex, apoe4 dose)
    note: str = ""

    def __repr__(self) -> str:  # pragma: no cover - display only
        return (f"<Cohort {self.name}: n={len(self.y)} "
                f"({int(self.y.sum())}+/{int((1 - self.y).sum())}-), "
                f"p={self.X.shape[1]}>")


def load_gse275392(stratum: str = "apoe33") -> Cohort:
    """
    Load the assembled GSE275392 plasma SomaScan matrix.

    Parameters
    ----------
    stratum : {'apoe33', 'full'}
        'apoe33' (PRE-REGISTERED PRIMARY) keeps only APOE e3/e3 subjects.  In the
        full cohort APOE genotype is *perfectly confounded* with amyloid status
        (all 18 amyloid-negatives are APOE33; all 17 APOE44 subjects are
        amyloid-positive), so a classifier can score well by learning APOE4
        rather than amyloid.  The APOE33 stratum is exactly balanced 18/18.
        'full' is secondary and must be labelled confounded wherever reported.

    Notes
    -----
    RFUs are log2-transformed (SomaScan RFUs are strongly right-skewed).
    Standardisation is deliberately NOT done here - it happens inside CV folds.
    """
    base = os.path.join(REPO, "datasets", "GSE275392", "processed")
    X = pd.read_csv(os.path.join(base, "proteins_rfu.csv"), index_col=0)
    ph = pd.read_csv(os.path.join(base, "phenotypes.csv"), index_col=0)
    assert (X.index == ph.index).all(), "sample index mismatch between matrix and phenotypes"

    if stratum == "apoe33":
        keep = ph["apoe_genotype"].astype(str) == "33"
        note = "PRIMARY: APOE e3/e3 only, deconfounded, 18+/18-"
    elif stratum == "full":
        keep = pd.Series(True, index=ph.index)
        note = "SECONDARY: full cohort - APOE genotype CONFOUNDED with amyloid status"
    else:
        raise ValueError(f"unknown stratum {stratum!r}")

    X, ph = X[keep.values], ph[keep.values]

    # log2 transform (RFUs are positive, right-skewed); guard against non-positives
    assert (X.values > 0).all(), "non-positive RFU encountered; log2 would fail"
    Xl = np.log2(X)

    y = (ph["amyloid_status"] == "positive").astype(int).values
    apoe = ph["apoe_genotype"].astype(str)
    cov = pd.DataFrame(
        {
            "age": ph["age"].values,
            "sex_M": (ph["Sex"] == "M").astype(int).values,
            # e4 allele dose: 33 -> 0, 34 -> 1, 44 -> 2
            "apoe4_dose": apoe.map(lambda g: sum(c == "4" for c in g)).values,
        },
        index=ph.index,
    )
    return Cohort(name=stratum, X=Xl, y=y, cov=cov, note=note)


# --------------------------------------------------------------------------- #
# AUC machinery
# --------------------------------------------------------------------------- #
def auc(y: np.ndarray, s: np.ndarray) -> float:
    """Mann-Whitney U AUC (ties handled via midranks). Returns nan if one class."""
    y = np.asarray(y).astype(int)
    s = np.asarray(s, dtype=float)
    n1, n0 = int(y.sum()), int((1 - y).sum())
    if n1 == 0 or n0 == 0:
        return float("nan")
    r = st.rankdata(s)
    return (r[y == 1].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0)


def _midrank(x: np.ndarray) -> np.ndarray:
    return st.rankdata(x)


def delong_cov(y: np.ndarray, scores: np.ndarray):
    """
    DeLong (1988) / Sun & Xu (2014) fast implementation.

    Parameters
    ----------
    y : (n,) binary labels
    scores : (m, n) array of m correlated score vectors on the SAME samples

    Returns
    -------
    aucs : (m,) AUC per score vector
    S    : (m, m) estimated covariance matrix of the AUCs
    """
    y = np.asarray(y).astype(int)
    scores = np.atleast_2d(np.asarray(scores, dtype=float))
    pos, neg = scores[:, y == 1], scores[:, y == 0]
    m, n = pos.shape[1], neg.shape[1]
    k = scores.shape[0]

    tx = np.vstack([_midrank(pos[i]) for i in range(k)])
    ty = np.vstack([_midrank(neg[i]) for i in range(k)])
    tz = np.vstack([_midrank(scores[i]) for i in range(k)])

    aucs = np.array([auc(y, scores[i]) for i in range(k)])

    v01 = (tz[:, :m] - tx) / n          # placement values for positives
    v10 = 1.0 - (tz[:, m:] - ty) / m    # placement values for negatives
    sx = np.cov(v01, ddof=1) if k > 1 else np.array([[np.var(v01[0], ddof=1)]])
    sy = np.cov(v10, ddof=1) if k > 1 else np.array([[np.var(v10[0], ddof=1)]])
    S = np.atleast_2d(sx) / m + np.atleast_2d(sy) / n
    return aucs, S


def delong_test(y: np.ndarray, s1: np.ndarray, s2: np.ndarray):
    """
    Paired DeLong test for two correlated ROC curves on the same samples.

    Returns (auc1, auc2, z, two_sided_p).  Returns z=nan, p=1.0 when the
    variance of the difference is numerically zero (identical score vectors).
    """
    # order positives first, as delong_cov expects when slicing tz
    y = np.asarray(y).astype(int)
    order = np.argsort(-y, kind="stable")
    yo = y[order]
    sc = np.vstack([np.asarray(s1, float)[order], np.asarray(s2, float)[order]])
    aucs, S = delong_cov(yo, sc)
    var = S[0, 0] + S[1, 1] - 2 * S[0, 1]
    if not np.isfinite(var) or var <= 0:
        return float(aucs[0]), float(aucs[1]), float("nan"), 1.0
    z = (aucs[0] - aucs[1]) / np.sqrt(var)
    p = 2 * st.norm.sf(abs(z))
    return float(aucs[0]), float(aucs[1]), float(z), float(p)


def boot_auc_ci(y, s, n_boot=2000, alpha=0.05, seed=SEED):
    """Stratified (class-preserving) bootstrap percentile CI for a single AUC."""
    rng = np.random.default_rng(seed)
    y = np.asarray(y).astype(int)
    s = np.asarray(s, float)
    ip, ineg = np.where(y == 1)[0], np.where(y == 0)[0]
    out = np.empty(n_boot)
    for b in range(n_boot):
        idx = np.concatenate([rng.choice(ip, ip.size, replace=True),
                              rng.choice(ineg, ineg.size, replace=True)])
        out[b] = auc(y[idx], s[idx])
    lo, hi = np.nanpercentile(out, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return float(lo), float(hi), out


def boot_paired_ratio(y, s_small, s_large, ratio="excess", anchor=0.5,
                      n_boot=2000, alpha=0.05, seed=SEED):
    """
    Bootstrap CI for the *ratio of performance retained* by a small panel
    relative to a large one, resampling subjects jointly (paired) so the
    correlation between the two score vectors is preserved.

    ratio : 'raw'     -> AUC_s / AUC_l                       (chance floor 0.5 kept)
            'excess'  -> (AUC_s-0.5)/(AUC_l-0.5)             (Gini / chance-corrected)
            'anchored'-> (AUC_s-a)/(AUC_l-a) with a=`anchor` (e.g. covariate-only AUC)
    """
    rng = np.random.default_rng(seed)
    y = np.asarray(y).astype(int)
    ss, sl = np.asarray(s_small, float), np.asarray(s_large, float)
    ip, ineg = np.where(y == 1)[0], np.where(y == 0)[0]
    a = {"raw": 0.0, "excess": 0.5, "anchored": anchor}[ratio]
    out = np.empty(n_boot)
    for b in range(n_boot):
        idx = np.concatenate([rng.choice(ip, ip.size, replace=True),
                              rng.choice(ineg, ineg.size, replace=True)])
        num, den = auc(y[idx], ss[idx]) - a, auc(y[idx], sl[idx]) - a
        out[b] = num / den if den > 1e-9 else np.nan
    lo, hi = np.nanpercentile(out, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    point = (auc(y, ss) - a) / max(auc(y, sl) - a, 1e-9)
    return float(point), float(lo), float(hi), out


def ratios(auc_small: float, auc_large: float, auc_anchor: float = 0.5) -> dict:
    """All three readings of 'fraction of performance retained'."""
    return {
        "ratio_raw": auc_small / auc_large if auc_large else np.nan,
        "ratio_excess": (auc_small - 0.5) / (auc_large - 0.5) if auc_large > 0.5 else np.nan,
        "ratio_anchored": ((auc_small - auc_anchor) / (auc_large - auc_anchor)
                           if auc_large > auc_anchor else np.nan),
    }


# --------------------------------------------------------------------------- #
# Screening operating characteristics
# --------------------------------------------------------------------------- #
def ppv_npv_at_prevalence(sens: float, spec: float, prev: float) -> tuple[float, float]:
    """Bayes-transport PPV/NPV to a target prevalence (AUC does not imply these)."""
    tp, fn = sens * prev, (1 - sens) * prev
    tn, fp = spec * (1 - prev), (1 - spec) * (1 - prev)
    ppv = tp / (tp + fp) if (tp + fp) > 0 else np.nan
    npv = tn / (tn + fn) if (tn + fn) > 0 else np.nan
    return float(ppv), float(npv)


def youden_operating_point(y, s):
    """Sensitivity/specificity at the Youden-J-maximising threshold."""
    from sklearn.metrics import roc_curve

    fpr, tpr, thr = roc_curve(y, s)
    j = int(np.argmax(tpr - fpr))
    return float(tpr[j]), float(1 - fpr[j]), float(thr[j])


# --------------------------------------------------------------------------- #
# Plot style
# --------------------------------------------------------------------------- #
def use_style():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "figure.dpi": 130, "savefig.dpi": 200, "savefig.bbox": "tight",
        "font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.6,
        "legend.frameon": False, "figure.facecolor": "white",
    })
    return plt


# Brand-neutral categorical palette (colour-blind safe, consistent across figures)
PAL = {
    "primary": "#2B6CB0", "accent": "#D97706", "green": "#059669",
    "red": "#DC2626", "purple": "#7C3AED", "grey": "#6B7280",
    "light": "#93C5FD", "teal": "#0D9488",
}


def dump_json(obj, path):
    """Write JSON, creating parent dirs, with numpy types coerced."""
    def conv(o):
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        raise TypeError(f"not serialisable: {type(o)}")

    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=2, default=conv)
    return path
