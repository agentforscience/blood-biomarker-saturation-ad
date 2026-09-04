"""
Figures and the cross-direction synthesis table.

Reads results/d1, results/d2, results/d3 and writes figures/fig2..fig6 plus
results/synthesis.csv / results/synthesis.json.  Run after d0/d1/d2/d3.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C  # noqa: E402

R = os.path.join(C.REPO, "results")
FIG = os.path.join(C.REPO, "figures")
os.makedirs(FIG, exist_ok=True)
plt = C.use_style()

COV_AUC_LIT = 0.773     # age+sex+APOE4 in CU adults (Trelle 2026)
PTAU_BENCH = 0.904      # single plasma p-tau217, ADNI n=1317


# --------------------------------------------------------------------------- #
def fig2_d1_curve():
    """D1: the empirical AUC(k) saturation curve on GSE275392."""
    out = {}
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3), sharey=True)
    for ax, stratum, title in [
            (axes[0], "apoe33", "A. PRIMARY: APOE ε3/ε3 stratum (n=36, 18+/18−)"),
            (axes[1], "full", "B. SECONDARY: full cohort (n=53) — APOE confounded")]:
        f = os.path.join(R, "d1", f"curve_{stratum}.csv")
        if not os.path.exists(f):
            continue
        d = pd.read_csv(f)
        curve = d[d.panel == "topk"]
        for model, col in [("l2", C.PAL["primary"]), ("enet", C.PAL["teal"]),
                           ("rf", C.PAL["accent"])]:
            s = curve[curve.model == model].sort_values("k")
            ax.plot(s.k, s.auc_mean, "o-", color=col, ms=4, lw=1.8, label=f"{model}")
            ax.fill_between(s.k, s.auc_mean - s.auc_sd, s.auc_mean + s.auc_sd,
                            color=col, alpha=0.12, lw=0)
        pf = os.path.join(R, "d1", f"permutation_{stratum}.csv")
        if os.path.exists(pf):
            p = pd.read_csv(pf).sort_values("k")
            ax.plot(p.k, p.perm_q95, ls="--", color=C.PAL["red"], lw=1.2,
                    label="95th pct. of label-permutation null")
        base = d[d.panel == "covariates_only"].auc_mean.iloc[0]
        ax.axhline(base, color=C.PAL["green"], ls=":", lw=1.5)
        ax.text(1.1, base + 0.006, f"covariates only ({base:.2f})",
                color=C.PAL["green"], fontsize=8)
        ax.axhline(PTAU_BENCH, color=C.PAL["purple"], ls="-.", lw=1.5)
        ax.text(1.1, PTAU_BENCH + 0.008,
                f"published single p-tau217 benchmark ({PTAU_BENCH})",
                color=C.PAL["purple"], fontsize=8)
        ax.axhline(0.5, color="k", lw=0.8, alpha=0.5)
        ax.set_xscale("log")
        ax.set_xlabel("panel size $k$ (proteins, log scale)")
        ax.set_title(title, fontsize=10)
        ax.legend(fontsize=7.5, loc="lower right", ncol=1)
        best = curve.loc[curve.auc_mean.idxmax()]
        out[stratum] = dict(best_k=int(best.k), best_model=best.model,
                            best_auc=float(best.auc_mean),
                            auc_at_kmax=float(curve[(curve.k == 1305) &
                                                    (curve.model == "l2")].auc_mean.iloc[0]),
                            covariate_auc=float(base))
        ax.annotate(f"peak {best.auc_mean:.3f}\nat k={int(best.k)}",
                    xy=(best.k, best.auc_mean), xytext=(best.k * 1.5, best.auc_mean + 0.09),
                    fontsize=8, arrowprops=dict(arrowstyle="->", lw=0.8))
    axes[0].set_ylabel("out-of-sample AUC (20×5 repeated CV, mean ± SD)")
    axes[0].set_ylim(0.30, 1.0)
    fig.suptitle("Panel-size saturation on a real 1,305-protein plasma proteome: "
                 "the curve peaks near k≈10–20 and never approaches the p-tau217 benchmark",
                 fontsize=10.5, y=1.03)
    fig.savefig(os.path.join(FIG, "fig2_d1_saturation.png"))
    plt.close(fig)
    return out


def fig3_d1_panels():
    """D1: a-priori panels and baselines vs the data-driven selection."""
    f = os.path.join(R, "d1", "panels_apoe33.csv")
    if not os.path.exists(f):
        return {}
    d = pd.read_csv(f)
    keep = ["covariates", "core_available3", "apriori_AD8", "biohermes7",
            "topk_l2_1", "topk_l2_3", "topk_l2_5", "topk_l2_10", "topk_l2_20",
            "topk_l2_100", "topk_l2_1305"]
    lab = {"covariates": "age + sex", "core_available3": "GFAP+MAPT+APP (k=3)",
           "apriori_AD8": "a-priori AD panel (k=8)", "biohermes7": "Bio-Hermes top-7",
           "topk_l2_1": "data-driven k=1", "topk_l2_3": "data-driven k=3",
           "topk_l2_5": "data-driven k=5", "topk_l2_10": "data-driven k=10",
           "topk_l2_20": "data-driven k=20", "topk_l2_100": "data-driven k=100",
           "topk_l2_1305": "data-driven k=1305 (all)"}
    d = d[d.panel.isin(keep)].set_index("panel").loc[[k for k in keep if k in set(d.panel)]]
    fig, ax = plt.subplots(figsize=(8.2, 4.6))
    y = np.arange(len(d))
    cols = [C.PAL["green"] if "cov" in i else
            (C.PAL["accent"] if not i.startswith("topk") else C.PAL["primary"])
            for i in d.index]
    ax.barh(y, d.auc, color=cols, height=0.62)
    ax.errorbar(d.auc, y, xerr=[d.auc - d.ci_lo, d.ci_hi - d.auc], fmt="none",
                ecolor="k", elinewidth=1, capsize=2.5)
    for yi, (a, p) in enumerate(zip(d.auc, d.delong_p_vs_full)):
        ax.text(max(d.ci_hi) + 0.01, yi, f"{a:.3f}   p={p:.2f}" if np.isfinite(p)
                else f"{a:.3f}", va="center", fontsize=8)
    ax.axvline(0.5, color="k", lw=1)
    ax.set_yticks(y)
    ax.set_yticklabels([lab[i] for i in d.index], fontsize=8.5)
    ax.set_xlabel("pooled out-of-fold AUC (95% stratified bootstrap CI); "
                  "p = paired DeLong vs the full 1,305-protein panel")
    ax.set_xlim(0.25, 1.02)
    ax.set_title("D1 panels on the deconfounded APOE ε3/ε3 stratum (n=36).\n"
                 "No panel — a-priori, data-driven, or exhaustive — separates from chance.",
                 fontsize=10)
    fig.savefig(os.path.join(FIG, "fig3_d1_panels.png"))
    plt.close(fig)
    return d.reset_index().to_dict("records")


def fig4_d2():
    """D2: retained-performance ratios across the simulation sweeps."""
    sweeps = [("S2_ninf", "n_informative", "number of informative extra proteins"),
              ("S3_deffect", "d_extra", "effect size of each extra protein (Cohen's d)"),
              ("S4_rho", "rho_core_extra", "extra↔core shared-factor loading"),
              ("S7_adversarial", "n_informative",
               "adversarial: independent extras, d=0.8")]
    avail = [(s, x, t) for s, x, t in sweeps
             if os.path.exists(os.path.join(R, "d2", f"summary_{s}.csv"))]
    if not avail:
        return {}
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.6))
    for ax, (name, xcol, xlab) in zip(axes.ravel(), avail):
        d = pd.read_csv(os.path.join(R, "d2", f"summary_{name}.csv")).sort_values(xcol)
        for key, col, lab in [
                ("core5_hypothesis__ratio_raw", C.PAL["grey"], "raw ratio"),
                ("core5_hypothesis__ratio_excess", C.PAL["primary"], "excess ratio"),
                ("core5_hypothesis__ratio_anchored", C.PAL["accent"], "anchored ratio")]:
            ax.plot(d[xcol], d[key], "o-", color=col, ms=4, lw=1.8, label=lab)
        ax.axhline(0.90, color=C.PAL["red"], ls="--", lw=1.2)
        ax.set_xlabel(xlab, fontsize=9)
        ax.set_ylabel("core-5 performance retained\nvs the full panel", fontsize=9)
        ax.set_title(name.split("_", 1)[1], fontsize=9.5)
        ax.set_ylim(0.4, 1.35)
    axes[0, 0].legend(fontsize=8, loc="lower left")
    axes[0, 0].text(0.02, 0.905, "90% bar", color=C.PAL["red"], fontsize=8,
                    transform=axes[0, 0].get_yaxis_transform())
    fig.suptitle("D2: when does the 5-analyte core stop retaining 90% of a "
                 "125-protein panel?  Only when many strong, independent extras exist.",
                 fontsize=10.5, y=1.0)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig4_d2_boundary.png"))
    plt.close(fig)

    # panel-comparison figure at the base condition
    f = os.path.join(R, "d2", "summary_S1_n.csv")
    if os.path.exists(f):
        d = pd.read_csv(f).sort_values("n_train")
        fig, ax = plt.subplots(figsize=(7.6, 4.4))
        series = [("covariates_only_mean", C.PAL["green"], "age+sex+APOE4 only"),
                  ("gfap_nfl_only_mean", C.PAL["red"], "GFAP + NfL (2)"),
                  ("extras_only_mean", C.PAL["grey"], "120 discovery proteins, no p-tau217"),
                  ("core1_ptau217_mean", C.PAL["purple"], "p-tau217 (1)"),
                  ("core2_ptau217_ab_mean", C.PAL["teal"], "p-tau217 + Aβ42/40 (2)"),
                  ("core5_hypothesis_mean", C.PAL["primary"], "hypothesised core-5"),
                  ("full_panel_mean", C.PAL["accent"], "core-5 + 120 extras (125)")]
        for key, col, lab in series:
            if key in d:
                ax.plot(d.n_train, d[key], "o-", color=col, ms=4, lw=1.8, label=lab)
        ax.set_xscale("log")
        ax.set_xlabel("training-set size")
        ax.set_ylabel("test AUC (n=8,000 held out; 200 reps)")
        ax.set_title("D2: adding 120 proteins to the core-5 does not help at any "
                     "sample size,\nand hurts below n≈1,000", fontsize=10)
        ax.legend(fontsize=7.5, loc="lower right", ncol=2)
        fig.savefig(os.path.join(FIG, "fig5_d2_panels.png"))
        plt.close(fig)
        return d.to_dict("records")
    return {}


def fig6_d3():
    """D3: the published (panel size -> AUC) relationship."""
    fa = os.path.join(R, "d3", "passA_extracted.csv")
    fb = os.path.join(R, "d3", "passB_curated.csv")
    if not (os.path.exists(fa) and os.path.exists(fb)):
        return {}
    A = pd.read_csv(fa)
    B = pd.read_csv(fb).dropna(subset=["auc"])
    fig, ax = plt.subplots(1, 2, figsize=(11.2, 4.4))

    a = ax[0]
    d = A.dropna(subset=["panel_size"])
    d = d[d.panel_dist_chars < 800]
    for lab, sub, col, mk in [("amyloid-status endpoint",
                               d[d.endpoint == "amyloid_status"], C.PAL["primary"], "o"),
                              ("other endpoints",
                               d[d.endpoint != "amyloid_status"], C.PAL["grey"], "x")]:
        a.scatter(sub.panel_size, sub.auc, s=22, color=col, marker=mk, alpha=0.7, label=lab)
    sub = d[d.endpoint == "amyloid_status"]
    if len(sub) > 3:
        import scipy.stats as st
        sl, ic, r, p, se = st.linregress(np.log10(sub.panel_size), sub.auc)
        xs = np.logspace(0, np.log10(max(sub.panel_size)), 50)
        a.plot(xs, ic + sl * np.log10(xs), color=C.PAL["primary"], lw=2)
        a.text(0.03, 0.06, f"slope = {sl:+.3f} AUC/decade\n(p = {p:.2f}, n = {len(sub)})",
               transform=a.transAxes, fontsize=8.5, color=C.PAL["primary"])
    a.set_xscale("log")
    a.set_xlabel("panel size mentioned near the AUC (log)")
    a.set_ylabel("reported AUC")
    a.set_title("A. Automated extraction, 548 abstracts + 44 full texts", fontsize=10)
    a.legend(fontsize=8, loc="lower right")

    b = ax[1]
    cu = B[B.population.astype(str).str.contains("CU")]
    mk = {"Trelle2026_SAMS": "o", "Trelle2026_AMASS": "s",
          "Ptau231GFAP2025": "^", "BioHermes2026": "D"}
    for study, g in B.groupby("study"):
        g = g[g.k > 0]
        if g.empty:
            continue
        b.scatter(g.k, g.auc, s=48, marker=mk.get(study, "P"),
                  label=f"{study} (n={int(g.n.iloc[0])})", alpha=0.85)
    b.axhline(COV_AUC_LIT, color=C.PAL["green"], ls=":", lw=1.6)
    b.text(1.05, COV_AUC_LIT + 0.006, "age+sex+APOE4 baseline (0.773)",
           color=C.PAL["green"], fontsize=8)
    b.set_xscale("log")
    b.set_xlabel("number of measured analytes in the model (log)")
    b.set_ylabel("reported AUC for amyloid positivity")
    b.set_title("B. Curated within-study evidence", fontsize=10)
    b.legend(fontsize=7, loc="lower right")
    b.set_ylim(0.5, 1.0)
    fig.suptitle("D3: across the published corpus, larger panels do not report "
                 "higher amyloid-detection AUCs", fontsize=10.5, y=1.02)
    fig.savefig(os.path.join(FIG, "fig6_d3_corpus.png"))
    plt.close(fig)
    return {}


def synthesis():
    """One table answering the hypothesis under each ratio, per direction."""
    rows = []
    # D1
    f = os.path.join(R, "d1", "panels_apoe33.csv")
    if os.path.exists(f):
        d = pd.read_csv(f).set_index("panel")
        ref = "topk_l2_1305"
        for k in [1, 3, 5, 10, 20]:
            key = f"topk_l2_{k}"
            if key in d.index:
                rows.append(dict(direction="D1 (GSE275392, APOE33 n=36)",
                                 compact=f"top-{k} proteins", k_small=k,
                                 auc_small=d.loc[key, "auc"],
                                 large="all 1,305 proteins", k_large=1305,
                                 auc_large=d.loc[ref, "auc"],
                                 ratio_raw=d.loc[key, "ratio_raw"],
                                 ratio_excess=d.loc[key, "ratio_excess"],
                                 ratio_anchored=d.loc[key, "ratio_anchored"]))
    # D2
    f = os.path.join(R, "d2", "summary_S1_n.csv")
    if os.path.exists(f):
        d = pd.read_csv(f)
        for _, r in d.iterrows():
            for p, lab, ks in [("core1_ptau217", "p-tau217 alone", 1),
                               ("core2_ptau217_ab", "p-tau217 + Aβ42/40", 2),
                               ("core5_hypothesis", "hypothesised core-5", 5),
                               ("gfap_nfl_only", "GFAP + NfL", 2)]:
                rows.append(dict(direction=f"D2 (simulation, n_train={int(r.n_train)})",
                                 compact=lab, k_small=ks,
                                 auc_small=r[f"{p}_mean"], large="core-5 + 120 extras",
                                 k_large=125, auc_large=r["full_panel_mean"],
                                 ratio_raw=r[f"{p}__ratio_raw"],
                                 ratio_excess=r[f"{p}__ratio_excess"],
                                 ratio_anchored=r[f"{p}__ratio_anchored"]))
    # D3
    f = os.path.join(R, "d3", "within_study_ratios.csv")
    if os.path.exists(f):
        d = pd.read_csv(f)
        for _, r in d.iterrows():
            rows.append(dict(direction=f"D3 ({r.study}, n={int(r.n)})",
                             compact=r.small_panel, k_small=int(r.k_small),
                             auc_small=r.auc_small, large=r.large_panel,
                             k_large=int(r.k_large), auc_large=r.auc_large,
                             ratio_raw=r.ratio_raw, ratio_excess=r.ratio_excess,
                             ratio_anchored=r.ratio_anchored))
    s = pd.DataFrame(rows)
    for c in ["ratio_raw", "ratio_excess", "ratio_anchored"]:
        s[f"{c}_ge90"] = s[c] >= 0.90
    s.to_csv(os.path.join(R, "synthesis.csv"), index=False)
    return s


if __name__ == "__main__":
    C.set_seed()
    out = {"fig2_d1": fig2_d1_curve(), "fig3_d1_panels": fig3_d1_panels(),
           "fig4_5_d2": fig4_d2(), "fig6_d3": fig6_d3()}
    s = synthesis()
    print(s.round(3).to_string(index=False))
    print("\nfraction of comparisons clearing the 90% bar:")
    print(s[["ratio_raw_ge90", "ratio_excess_ge90", "ratio_anchored_ge90"]].mean().round(3))
    C.dump_json(out, os.path.join(R, "figure_meta.json"))
    print("\nfigures ->", FIG)
