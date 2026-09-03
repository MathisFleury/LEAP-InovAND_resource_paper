#!/usr/bin/env python3
# =============================================================================
# 07 - IQ figures: (1) IQ x PGS-intelligence scatter frame, (2) IQ ~ LOEUF betas
# =============================================================================
# Two Section-2 figures, both keyed on measured IQ (the `IQ` column: most
# complete, ~1.7k non-null; NOT total/full-scale IQ) and gnomAD v4 genetics:
#
#   Figure 1 (frame only; plotted by 07b_plot_iq_pgs_cluster_loeuf.R):
#     Cluster-coloured, LOEUF-sized IQ (y) vs PGS-intelligence (x) scatter --
#     the analogue of 10_clinical_analysis's figure_IQ_SRS_cluster_loeuf_16cm
#     with the SRS-2 x-axis swapped for the intelligence PGS (int_savage2018).
#     This step reuses the tested cluster/LOEUF/gene-set frame logic from
#     10_clinical_analysis/04_plot_iq_srs_scatter.py and joins the PGS on
#     barcode; R adds the OLS regression line + R^2 / p annotation.
#     -> tables_v4/iq_pgs_plot_frame.csv
#
#   Figure 2 (computed AND plotted here):
#     Per-gene-list beta of (IQ ~ -log10 LOEUF best-score), +/-1.96 SE, with
#     permutation p -- the analogue of 4_anatomical_analysis's
#     beta_coefficients_*_loeuf.pdf with cortical thickness swapped for IQ.
#     Same OLS + gene-list-specific permutation null as
#     08_loeuf_mri_regression_permutation_v2.py.
#     -> tables_v4/iq_loeuf_regperm.csv, figures_v4/beta_coefficients_iq_loeuf.pdf
#
# Run:  python3.11 07_iq_pgs_loeuf_figures.py
# =============================================================================

import os
import sys

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats as sps
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _config import (
    CARRIER_V4, CLUSTERS_CURATED_FILE, PGS_EUR_V4,
    TABLES_DIR_V4, FIGURES_DIR_V4, N_PERM, SEED,
)
PGS_FILE = PGS_EUR_V4  # European-ancestry v4 PGS (int_savage2018 on sample_id)

PGS_INT_COL = "int_savage2018"  # SBayesRC intelligence PGS (Savage 2018)

# --- LOEUF binning (copied from 10_clinical_analysis/04, band 5 = most constrained)
LOEUF_EDGES = [0.10, 0.20, 0.36, 0.45]
POP_RELABEL = {"TD": "NT", "ID": "IDD",
               "Autism with ID": "Autism with IDD",
               "Autism without ID": "Autism without IDD"}

# Fig-2 gene lists: (label, lof best-score col, del best-score col). dellof = fmin.
GENE_SETS = [
    ("Protein coding", "lof_proteincoding_genes_best_score",
     "del_proteincoding_genes_best_score"),
    ("SynGO", "lof_syngo_constraint_genes_best_score",
     "del_syngo_constraint_genes_best_score"),
    ("CHROM", "lof_chromepitf_constraint_genes_best_score",
     "del_chromepitf_constraint_genes_best_score"),
]


def neg_log10(series):
    """-log10 of the LOEUF best-score; drops NaN and non-positive (paper method)."""
    s = series.dropna()
    return -np.log10(s[s > 0])


def load_cohort():
    """Curated clusters (IQ, cluster, population) inner-joined to the v4 carrier
    matrix (LOEUF best-scores + gene-set flags), then PGS on barcode."""
    df = pd.read_csv(CLUSTERS_CURATED_FILE, low_memory=False)
    df["ID"] = df["ID"].astype(str)
    df = df.drop_duplicates("ID")
    df["Population1"] = df["Population1"].replace(POP_RELABEL)
    df.loc[df["Population1"] == "NT", "Cluster"] = "NT"    # NT gets its own colour
    df.loc[df["Population1"] == "IDD", "Cluster"] = "IDD"  # IDD too (like NT)

    gen = pd.read_table(CARRIER_V4, low_memory=False)
    gen["ID"] = gen["ID"].astype(str)
    gene_cols = ["dellof_constraint_genes_best_gene", "dellof_constraint_genes_best_score",
                 "dellof_syngo_constraint_genes_best_gene", "dellof_syngo_constraint_genes_best_score",
                 "dellof_chromepitf_constraint_genes_best_gene", "dellof_chromepitf_constraint_genes_best_score",
                 "dellof_syngo_chromepitf_constraint_genes_best_gene",
                 "dellof_syngo_chromepitf_constraint_genes_best_score"]
    carrier_cols = ["dellof_syngo_contraint_carrier", "dellof_chromepitf_contraint_carrier"]
    score_cols = sorted({c for _, lof, dele in GENE_SETS for c in (lof, dele)})
    keep = ["ID", "barcode"] + gene_cols + carrier_cols + score_cols
    df = df.merge(gen[keep], on="ID", how="inner")

    pgs = pd.read_table(PGS_FILE, low_memory=False)
    pgs = pgs[["sample_id", PGS_INT_COL]].rename(columns={"sample_id": "barcode"})
    df = df.merge(pgs, how="left", on="barcode")
    print(f"  cohort: {len(df)} · PGS-int non-null: {df[PGS_INT_COL].notna().sum()}")
    return df


# ---------------------------------------------------------------------------
# Figure 1 frame
# ---------------------------------------------------------------------------
def build_pgs_frame(df):
    """Reproduce the 04 IQ×SRS frame's gene/LOEUF/gene-set columns, but keyed on
    PGS-intelligence instead of SRS. Written for the R (ggrepel) plotter."""
    d = df.dropna(subset=["IQ", PGS_INT_COL, "Cluster"]).copy()

    sy = d["dellof_syngo_contraint_carrier"] == 1
    ch = d["dellof_chromepitf_contraint_carrier"] == 1
    d["gs_syngo"] = (sy & ~ch).astype(int)
    d["gs_chrom"] = (ch & ~sy).astype(int)
    d["gs_both"] = (sy & ch).astype(int)

    d["gene"] = d["dellof_constraint_genes_best_gene"]
    d["LOEUF_raw"] = d["dellof_constraint_genes_best_score"].replace([np.inf, -np.inf], np.nan)
    for mask, gcol, scol in [
        (d["gs_syngo"] == 1, "dellof_syngo_constraint_genes_best_gene",
         "dellof_syngo_constraint_genes_best_score"),
        (d["gs_chrom"] == 1, "dellof_chromepitf_constraint_genes_best_gene",
         "dellof_chromepitf_constraint_genes_best_score"),
        (d["gs_both"] == 1, "dellof_syngo_chromepitf_constraint_genes_best_gene",
         "dellof_syngo_chromepitf_constraint_genes_best_score"),
    ]:
        d.loc[mask, "gene"] = d.loc[mask, gcol]
        d.loc[mask, "LOEUF_raw"] = d.loc[mask, scol].replace([np.inf, -np.inf], np.nan)

    d["was_nan_LOEUF"] = d["LOEUF_raw"].isna()
    lo = d["LOEUF_raw"].fillna(2.0)  # missing -> smallest band
    d["LOEUF"] = np.select(
        [lo < LOEUF_EDGES[0], lo < LOEUF_EDGES[1], lo < LOEUF_EDGES[2], lo < LOEUF_EDGES[3]],
        [5, 4, 3, 2], default=1)

    cols = ["ID", PGS_INT_COL, "IQ", "Cluster", "Population1",
            "LOEUF", "was_nan_LOEUF", "LOEUF_raw", "gene",
            "gs_syngo", "gs_chrom", "gs_both"]
    out = os.path.join(TABLES_DIR_V4, "iq_pgs_plot_frame.csv")
    os.makedirs(TABLES_DIR_V4, exist_ok=True)
    d[cols].rename(columns={PGS_INT_COL: "PGS_int"}).to_csv(out, index=False)
    print(f"  Saved: {out}  ({len(d)} rows)")


# ---------------------------------------------------------------------------
# Figure 2: IQ ~ -log10 LOEUF beta coefficients + permutation test
# ---------------------------------------------------------------------------
XLABEL = {"Protein coding": "All protein\ncoding genes", "CHROM": "ChromEpiTF", "SynGO": "SYNGO"}
COLORS = {"Protein coding": "gray", "CHROM": "#8E44AD", "SynGO": "#4A90E2"}
ORDER = ["Protein coding", "CHROM", "SynGO"]


def _residuals_not_normal(resid):
    """Four complementary normality tests (paper method); flag if any p<0.05."""
    flags = []
    if len(resid) < 5000:
        flags.append(sps.shapiro(resid).pvalue < 0.05)
    flags.append(sps.normaltest(resid).pvalue < 0.05)          # D'Agostino
    flags.append(sps.jarque_bera(resid).pvalue < 0.05)
    ad = sps.anderson(resid, dist="norm")                      # 5% crit = index 2
    flags.append(ad.statistic > ad.critical_values[2])
    return any(flags)


def _bootstrap_median_beta(X, y, n=1000):
    """Non-parametric robust beta: median of 1000 resample-with-replacement OLS
    slopes (seed 42, paper method)."""
    np.random.seed(42)
    betas = []
    for _ in range(n):
        idx = np.random.choice(len(X), size=len(X), replace=True)
        try:
            betas.append(sm.OLS(y[idx], sm.add_constant(X[idx])).fit().params[1])
        except Exception:
            continue
    return float(np.median(betas)) if betas else None


def _perm_slopes(Xp, y):
    xc = Xp - Xp.mean(axis=1, keepdims=True)
    yc = y - y.mean()
    denom = (xc ** 2).sum(axis=1)
    denom[denom == 0] = np.nan
    return (xc * yc).sum(axis=1) / denom


def iq_loeuf_regression(df):
    """OLS of IQ on -log10(LOEUF best-score) per gene list, across ALL carriers
    (whole cohort, no population filter), with the same gene-list-specific
    permutation null as the MRI regperm (script 08)."""
    aut = df[df["IQ"].notna()].copy()
    print(f"  cohort for IQ~LOEUF (all carriers): {len(aut)}")

    for label, lof, dele in GENE_SETS:
        aut[label] = np.fmin(aut[lof].values, aut[dele].values)

    rng = np.random.default_rng(SEED)
    background = neg_log10(aut["Protein coding"]).dropna().values
    iq = aut["IQ"]
    rows = []
    for label, *_ in GENE_SETS:
        x_full = neg_log10(aut[label])
        common = x_full.index.intersection(iq.dropna().index)
        if len(common) < 10:
            print(f"  skip {label}: n<10"); continue
        X = x_full.loc[common].values.astype(float)
        y = iq.loc[common].values.astype(float)
        model = sm.OLS(y, sm.add_constant(X)).fit()
        beta, se, p_param = float(model.params[1]), float(model.bse[1]), float(model.pvalues[1])

        # Permutation p (gene-list-specific null); observed beta on OLS scale.
        if label == "Protein coding":
            Xp = np.array([rng.permutation(X) for _ in range(N_PERM)])
        else:
            Xp = rng.choice(background, size=(N_PERM, len(X)), replace=True)
        perm = _perm_slopes(Xp, y)
        perm = perm[np.isfinite(perm)]
        n_ext = int(np.sum(perm >= beta) if beta >= 0 else np.sum(perm <= beta))
        p_perm = min((n_ext + 1) / (N_PERM + 1) * 2, 1.0) if perm.size else np.nan

        # Robust beta: bootstrap median when OLS residuals fail normality.
        not_normal = _residuals_not_normal(model.resid)
        boot_beta = _bootstrap_median_beta(X, y) if not_normal else None
        final_beta = boot_beta if boot_beta is not None else beta
        rows.append({"gene_list_label": label, "n": len(common),
                     "beta": final_beta, "beta_ols": beta, "se": se,
                     "p_value": p_param, "p_permutation": p_perm,
                     "residuals_normal": not not_normal})
        print(f"    {label:14s} n={len(common):4d} beta={final_beta:+.3f} "
              f"se(1x)={se:.3f} p_perm={p_perm:.4f} "
              f"resid_normal={not not_normal}")
    return pd.DataFrame(rows).set_index("gene_list_label")


def plot_beta(stats):
    present = [g for g in ORDER if g in stats.index]
    fig, ax = plt.subplots(figsize=(4.2, 4.6))
    ax.axhline(0, ls="--", color="grey", lw=1)
    for i, g in enumerate(present):
        r = stats.loc[g]
        beta, se, p = r["beta"], r["se"], r["p_permutation"]
        sig = p < 0.05
        ax.errorbar(i, beta, yerr=se, fmt="o", color=COLORS[g],  # +/-1 SE (paper)
                    capsize=5, capthick=2, markersize=12 if sig else 10,
                    markeredgecolor="black" if sig else COLORS[g],
                    markeredgewidth=1.5 if sig else 0.5,
                    elinewidth=2.5 if sig else 2, alpha=0.85, zorder=3)
        if sig:
            ax.text(i, beta + se + 0.03 * abs(beta), f"p = {p:.3f}",
                    ha="center", va="bottom", fontsize=10, fontweight="bold")

    ymin, ymax = ax.get_ylim()
    ax.set_ylim(ymin, ymax + 0.10 * (ymax - ymin))
    ax.set_xticks(range(len(present)))
    ax.set_xticklabels([XLABEL[g] for g in present])
    ax.set_xlim(-0.6, len(present) - 0.4)
    ax.set_xlabel("Gene List")
    ax.set_ylabel("Beta Coefficient\n(Measured IQ vs -log LOEUF)", fontsize=10)
    ax.set_title("Beta Coefficients\nMeasured IQ vs -log LOEUF",
                 fontsize=12, fontweight="bold")
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    handles = [Line2D([0], [0], marker="o", color="w", label=XLABEL[g],
                      markerfacecolor=COLORS[g], markersize=10,
                      markeredgecolor="black" if stats.loc[g, "p_permutation"] < 0.05 else COLORS[g])
               for g in present]
    ax.legend(handles=handles, title="Gene List", loc="center left",
              bbox_to_anchor=(1.02, 0.5), frameon=True, fontsize=9)
    fig.tight_layout()
    os.makedirs(FIGURES_DIR_V4, exist_ok=True)
    out = os.path.join(FIGURES_DIR_V4, "beta_coefficients_iq_loeuf.pdf")
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out}")


def main():
    print("Loading cohort (curated clusters x v4 carriers x PGS)...")
    df = load_cohort()

    print("\n=== Figure 1: IQ x PGS-intelligence frame ===")
    build_pgs_frame(df)

    print("\n=== Figure 2: IQ ~ -log10 LOEUF beta coefficients ===")
    stats = iq_loeuf_regression(df)
    if stats.empty:
        print("  No stats produced."); return
    stats.to_csv(os.path.join(TABLES_DIR_V4, "iq_loeuf_regperm.csv"))
    plot_beta(stats)


if __name__ == "__main__":
    main()
