"""
D0 - What does "90% of the predictive performance" actually mean?

The hypothesis states that a compact panel achieves ">=90% of the predictive
performance" of a large panel.  In the biomarker literature this is read as a
ratio of AUCs.  That reading is degenerate, because AUC has a chance floor of
0.5: a coin flip already achieves 0.5/0.94 = 53% of an AUC of 0.94, and a model
that is *significantly worse than age+sex+APOE* can still clear the 90% bar.

This module characterises three readings analytically and shows where they
disagree.  Nothing here is fitted to data; it is a property of the metric.

    raw       R_raw(s,l)      = s / l
    excess    R_exc(s,l)      = (s - 0.5) / (l - 0.5)          [Gini ratio]
    anchored  R_anc(s,l;c)    = (s - c) / (l - c)              [c = covariate-only AUC]

Outputs: results/d0/*.csv, figures/fig1_metric_behaviour.png
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

OUT = os.path.join(C.REPO, "results", "d0")
FIG = os.path.join(C.REPO, "figures")
os.makedirs(OUT, exist_ok=True)
os.makedirs(FIG, exist_ok=True)

COV_AUC = 0.773   # age+sex+APOE4 in cognitively unimpaired adults (Trelle 2026, AMASS)


def required_small_auc(auc_large: float, thresh: float = 0.90,
                       cov: float = COV_AUC) -> dict:
    """Minimum compact-panel AUC needed to clear `thresh` under each ratio."""
    return {
        "auc_large": auc_large,
        "req_raw": thresh * auc_large,
        "req_excess": 0.5 + thresh * (auc_large - 0.5),
        "req_anchored": (cov + thresh * (auc_large - cov)
                         if auc_large > cov else np.nan),
    }


def reference_models(auc_large: float, cov: float = COV_AUC) -> dict:
    """What do trivial / baseline models score under each ratio?"""
    def r(s):
        return C.ratios(s, auc_large, cov)
    return {
        "auc_large": auc_large,
        **{f"coinflip_{k}": v for k, v in r(0.50).items()},
        **{f"covariates_{k}": v for k, v in r(cov).items()},
        **{f"nfl_alone_{k}": v for k, v in r(0.582).items()},   # Trelle AMASS
        **{f"gfap_alone_{k}": v for k, v in r(0.646).items()},  # Trelle AMASS
        **{f"ptau217_{k}": v for k, v in r(0.904).items()},     # ADNI
    }


def main():
    C.set_seed()
    grid = np.round(np.arange(0.60, 1.00, 0.02), 3)
    req = pd.DataFrame([required_small_auc(a) for a in grid])
    req.to_csv(os.path.join(OUT, "required_small_auc.csv"), index=False)

    ref = pd.DataFrame([reference_models(a) for a in [0.80, 0.90, 0.94]])
    ref.to_csv(os.path.join(OUT, "reference_models.csv"), index=False)

    print("--- Minimum compact-panel AUC required to clear the '>=90%' bar ---")
    print(req[req.auc_large.isin([0.80, 0.86, 0.90, 0.94, 0.98])]
          .round(3).to_string(index=False))

    print("\n--- What trivial and baseline models score under each ratio ---")
    show = ["auc_large"] + [c for c in ref.columns if c != "auc_large"]
    print(ref[show].round(3).to_string(index=False))

    # The headline number: a coin flip's raw ratio against a strong large panel
    print("\nA COIN FLIP (AUC 0.500) scores "
          f"{0.5/0.94:.1%} of an AUC-0.94 panel under the raw ratio, "
          f"{0.0:.1%} under the excess ratio, and "
          f"{(0.5-COV_AUC)/(0.94-COV_AUC):+.1%} under the covariate-anchored ratio.")
    print("Plasma NfL alone (AUC 0.582, significantly WORSE than age+sex+APOE for "
          "amyloid in CU adults) scores "
          f"{0.582/0.94:.1%} / {(0.582-0.5)/(0.94-0.5):.1%} / "
          f"{(0.582-COV_AUC)/(0.94-COV_AUC):+.1%} under raw / excess / anchored.")

    # ---- figure -------------------------------------------------------- #
    plt = C.use_style()
    fig, ax = plt.subplots(1, 2, figsize=(10.5, 4.1))

    a = ax[0]
    a.plot(req.auc_large, req.req_raw, color=C.PAL["grey"], lw=2, label="raw  $s/l$")
    a.plot(req.auc_large, req.req_excess, color=C.PAL["primary"], lw=2,
           label="excess  $(s-0.5)/(l-0.5)$")
    a.plot(req.auc_large, req.req_anchored, color=C.PAL["accent"], lw=2,
           label=f"anchored  $(s-c)/(l-c)$, $c$={COV_AUC}")
    a.plot(req.auc_large, req.auc_large, ls=":", color="k", lw=1, label="parity")
    a.axhline(0.5, color=C.PAL["red"], ls="--", lw=1)
    a.text(0.605, 0.512, "chance", color=C.PAL["red"], fontsize=8)
    a.axhline(COV_AUC, color=C.PAL["green"], ls="--", lw=1)
    a.text(0.605, COV_AUC + 0.012, "age+sex+APOE baseline", color=C.PAL["green"], fontsize=8)
    a.set_xlabel("AUC of the large (20+ protein) panel")
    a.set_ylabel("Compact-panel AUC needed to clear '90%'")
    a.set_title("A. The '90%' bar under three readings")
    a.legend(loc="upper left", fontsize=8)
    a.set_ylim(0.45, 1.0)

    b = ax[1]
    models = [("coin flip", 0.500), ("NfL alone", 0.582), ("GFAP alone", 0.646),
              ("age+sex+APOE", COV_AUC), ("p-tau217", 0.904)]
    x = np.arange(len(models))
    w = 0.26
    L = 0.94   # comparison large panel: NULISAseq 124-plex best model (Trelle SAMS)
    for i, (key, col, lab) in enumerate([
            ("ratio_raw", C.PAL["grey"], "raw"),
            ("ratio_excess", C.PAL["primary"], "excess"),
            ("ratio_anchored", C.PAL["accent"], "anchored")]):
        vals = [C.ratios(v, L, COV_AUC)[key] for _, v in models]
        bars = b.bar(x + (i - 1) * w, vals, w, color=col, label=lab)
        for xi, v in zip(x + (i - 1) * w, vals):
            b.text(xi, max(v, 0) + 0.02, f"{v:.2f}", ha="center", fontsize=7,
                   color=col)
    b.axhline(0.90, color=C.PAL["red"], ls="--", lw=1.2)
    b.text(len(models) - 0.55, 0.92, "90% bar", color=C.PAL["red"], fontsize=8)
    b.axhline(0, color="k", lw=0.8)
    b.set_xticks(x)
    b.set_xticklabels([m for m, _ in models], rotation=18, ha="right", fontsize=8)
    b.set_ylabel(f"fraction of a large panel's performance\n(large panel AUC = {L})")
    b.set_title("B. Same models, three verdicts")
    b.legend(fontsize=8, loc="lower right")
    b.set_ylim(-0.9, 1.25)

    fig.suptitle("The '≥90% of performance' claim is only falsifiable once the "
                 "chance floor is removed", fontsize=11, y=1.02)
    p = os.path.join(FIG, "fig1_metric_behaviour.png")
    fig.savefig(p)
    print("\nfigure ->", p)
    C.dump_json({"env": C.env_report(), "cov_auc": COV_AUC},
                os.path.join(OUT, "d0_meta.json"))


if __name__ == "__main__":
    main()
