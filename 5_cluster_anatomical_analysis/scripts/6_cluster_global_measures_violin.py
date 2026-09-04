#!/usr/bin/env python3
"""
06 - Global structural-MRI measures by cluster (violin) — REVISION

Violin plots of three global FreeSurfer measures — total mean cortical
thickness, total intracranial volume (eTIV), and total white surface area —
z-scored, comparing autism clusters C1/C2/C3 and the IDD group with NT.

Each non-NT group (C1/C2/C3 and IDD) is compared to NT with a Welch t-test or
Mann-Whitney U (chosen per comparison by a Shapiro-Wilk normality check on both
groups), FDR-BH corrected across all measure x group comparisons. PFDR < 0.05
is annotated.

Source: the SAME regressed z-score FreeSurfer table used by the rest of this
section (1_generate_cluster_mri_inputs.py) — QC + ComBat + age/sex/eTIV
regression, full-sample z-scored. The three globals are read from that file:
    total mean cortical thickness = mean(lh_MeanThickness, rh_MeanThickness)
    total intracranial volume     = eTIV
    total white surface area       = mean(lh_WhiteSurfArea, rh_WhiteSurfArea)
Values are already z-scored in that table, so no further z-scoring is applied.
Data prep (wave selection, join keys, NT pool) is reused verbatim from 01.

Cluster labels: curated k-means (individuals_metrics_with_clusters_curated.csv).
Aesthetic matches 10_clinical_analysis/.../01_plot_psychomotor_milestones.py.

Output: outputs/figures/cluster_global_measures_violin.pdf
        outputs/tables/cluster_global_measures_stats.csv
"""

import os
import sys
import importlib.util
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Patch
from scipy.stats import ttest_ind, mannwhitneyu, shapiro
from statsmodels.stats.multitest import multipletests

_SCRIPT_DIR = Path(__file__).parent
_SECTION_DIR = _SCRIPT_DIR.parent
_PROJECT_DIR = _SECTION_DIR.parent

# --- palette from project config ---------------------------------------------
sys.path.insert(0, str(_PROJECT_DIR / "2_genetic_analysis" / "scripts"))
from _config import PALETTE_CLUSTERS  # noqa: E402

# --- default to curated k-means labels unless the caller overrode CLUSTER_FILE-
os.environ.setdefault("CLUSTER_FILE", str(
    _PROJECT_DIR / "1_clustering" / "outputs" / "tables" /
    "individuals_metrics_with_clusters_curated.csv"))

# --- reuse 01's loader verbatim (same MRI table, QC, wave-selection, NT pool) -
_spec = importlib.util.spec_from_file_location(
    "_gen01", _SCRIPT_DIR / "1_generate_cluster_mri_inputs.py")
_gen01 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_gen01)

_OUT_BASE = _SECTION_DIR / os.environ.get("CLUSTER_OUT_DIR", "outputs")
OUT_DIR = _OUT_BASE / "figures"
TABLES_DIR = _OUT_BASE / "tables"

# measure key -> (source columns to average, axis label)
MEASURES = {
    "totalMeanCorticalThickness": (["lh_MeanThickness", "rh_MeanThickness"],
                                   "Total mean\ncortical thickness"),
    "EstimatedTotalIntraCranialVol": (["eTIV"],
                                      "Total intracranial\nvolume (eTIV)"),
    "totalWhiteSurfArea": (["lh_WhiteSurfArea", "rh_WhiteSurfArea"],
                           "Total white\nsurface area"),
}
# Non-NT groups compared against NT. IDD is its own category (PopulationS1=="IDD"),
# separate from the autism clusters — same convention as the milestone figure.
COMPARE = ["C1", "C2", "C3", "IDD"]
GROUPS = ["NT"] + COMPARE

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["font.family"] = "sans-serif"


def load_data():
    df = _gen01.load_and_prepare_data()
    # Build the three globals from the (already z-scored) source columns.
    for key, (cols, _) in MEASURES.items():
        have = [c for c in cols if c in df.columns]
        if not have:
            raise KeyError(f"missing source columns for {key}: {cols}")
        df[key] = df[have].apply(pd.to_numeric, errors="coerce").mean(axis=1)

    # group: NT / IDD from phenotype, C1/C2/C3 for autism subjects with a cluster.
    df["group"] = np.select(
        [df["PopulationS1"] == "NT", df["PopulationS1"] == "IDD",
         df["PopulationS1"] == "Autism"],
        ["NT", "IDD", df["Cluster"]], default=None)
    df = df[df["group"].isin(GROUPS)].copy()
    print("Group Ns:", df["group"].value_counts().reindex(GROUPS).to_dict())
    return df


def compute_stats(df):
    """cluster-vs-NT per measure; auto t-test/MWU by Shapiro; FDR-BH over all."""
    nt = df[df["group"] == "NT"]
    rows = []
    for key, (_, label) in MEASURES.items():
        b = nt[key].dropna().values
        for cl in COMPARE:
            a = df.loc[df["group"] == cl, key].dropna().values
            if len(a) < 3 or len(b) < 3:
                continue
            normal = shapiro(a)[1] > 0.05 and shapiro(b)[1] > 0.05
            if normal:
                stat, p = ttest_ind(a, b, equal_var=False)
                test = "t-test"
            else:
                stat, p = mannwhitneyu(a, b, alternative="two-sided")
                test = "MWU"
            rows.append({"measure": key, "label": label.replace("\n", " "),
                         "cluster": cl, "test": test, "statistic": stat,
                         "p_value": p, "n_cluster": len(a), "n_nt": len(b),
                         "median_cluster": np.median(a), "median_nt": np.median(b)})
    res = pd.DataFrame(rows)
    res["p_fdr"] = multipletests(res["p_value"], method="fdr_bh")[1]
    return res


def plot(df, res, out_pdf):
    """Violin plot styled to match plot_combined_milestones()."""
    n_g = len(GROUPS)
    spacing = n_g + 2  # gap between measure groups (milestone convention)
    fig, ax = plt.subplots(figsize=(max(14, len(MEASURES) * 3.5), 6))
    colors = [PALETTE_CLUSTERS[g] for g in GROUPS]

    # --- violins + median dot / mean±95%CI line + n= labels ------------------
    for i, key in enumerate(MEASURES):
        for j, g in enumerate(GROUPS):
            vals = df.loc[df["group"] == g, key].dropna().values
            if len(vals) < 2:
                continue
            x = j + i * spacing
            vp = ax.violinplot([vals], positions=[x], widths=0.8,
                               showmeans=False, showmedians=False)
            vp["bodies"][0].set_facecolor(colors[j])
            vp["bodies"][0].set_alpha(0.7)
            for k in ("cmaxes", "cmins", "cbars"):
                if k in vp:
                    vp[k].set_visible(False)

            # mean + 95% CI (dot + vertical line — same marker style as milestone)
            m = np.mean(vals)
            ci = 1.959964 * np.std(vals, ddof=1) / np.sqrt(len(vals))
            ax.plot([x, x], [m - ci, m + ci], "k-", linewidth=1)
            ax.plot(x, m, "ko", markersize=5)

            # n= label
            y_top = np.percentile(vals, 98)
            ax.text(x, y_top + 0.15, f"n={len(vals)}", ha="center", va="bottom",
                    fontweight="bold", fontsize=8)

    # --- significance brackets (cluster vs NT, p-value text) -----------------
    y_min_all, y_max_all = ax.get_ylim()
    span = y_max_all - y_min_all
    for i, key in enumerate(MEASURES):
        sig = res[(res["measure"] == key) & (res["p_fdr"] < 0.05)]
        if sig.empty:
            continue
        local_max = max(np.percentile(df.loc[df["group"] == g, key].dropna().values, 98)
                        for g in GROUPS
                        if (df["group"] == g).sum() >= 2)
        y_base = local_max + span * 0.06
        step = span * 0.06
        x_nt = 0 + i * spacing
        for k, (_, row) in enumerate(sig.iterrows()):
            x_cl = GROUPS.index(row["cluster"]) + i * spacing
            y = y_base + k * step
            ax.plot([x_nt, x_cl], [y, y], "k-", linewidth=1.5)
            ax.text((x_nt + x_cl) / 2, y + step * 0.12, f"{row['p_fdr']:.2e}",
                    ha="center", va="bottom", fontsize=8, fontweight="bold")

    # --- axes ----------------------------------------------------------------
    centers = [(n_g - 1) / 2 + i * spacing for i in range(len(MEASURES))]
    ax.set_xticks(centers)
    ax.set_xticklabels([lbl for _, lbl in MEASURES.values()],
                       fontweight="bold", fontsize=10)
    ax.set_ylabel("z-score", fontsize=12)
    ax.set_title("Global structural MRI measures by cluster",
                 fontsize=16, fontweight="bold")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1] * 1.15)
    ax.legend(handles=[Patch(facecolor=PALETTE_CLUSTERS[g], alpha=0.7, label=g)
                       for g in GROUPS], loc="upper right")

    plt.tight_layout()
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_pdf, dpi=500, format="pdf", transparent=True, bbox_inches="tight")
    plt.close()
    print(f"  Saved: {out_pdf}")


def main():
    print("=" * 60)
    print("Global structural-MRI measures by cluster (violin)")
    print("=" * 60)
    df = load_data()
    res = compute_stats(df)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    csv = TABLES_DIR / "cluster_global_measures_stats.csv"
    res.to_csv(csv, index=False)
    print(f"  Saved: {csv}")
    print(res[["label", "cluster", "test", "p_value", "p_fdr"]].to_string(index=False))
    plot(df, res, OUT_DIR / "cluster_global_measures_violin.pdf")
    print("\nDone.")


if __name__ == "__main__":
    main()
