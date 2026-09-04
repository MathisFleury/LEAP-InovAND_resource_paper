#!/usr/bin/env python3
"""
Euler number by population — REVISION

Reports the FreeSurfer Euler number (``mean_euler`` column in the new
QC+ComBat+regression TSV) **per population type**, using the same wave
selection as 1_anatomical_mri_autism_nt_v2.py (one row per subject;
LEAP_W1 / INOVAND_T1 priority).

Output is analogous to the manuscript's existing Supplementary Fig. 19a
(mean-FD-by-group plot), but for the FreeSurfer Euler number — a single
per-subject scalar that summarises cortical-surface reconstruction quality.
Lower (more negative) values reflect more topological defects in the
reconstructed surfaces.

Populations covered (from df_clusters_complete_kmeans.csv):
  Autism, NT, ID, Relatives (when available).

Outputs (under ../outputs/):
- tables/euler_by_population.csv          per-subject Euler values + label
- tables/euler_by_population_summary.csv  per-group n / mean / SD / median /
                                           p (Welch's t vs NT) / Cohen's d
- figures/euler_by_population.pdf         violin + stripplot + group mean ±
                                           CI (manuscript palette, vector PDF)
"""

from pathlib import Path
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# Wave selection helpers + paths are inherited from the v2 generator.
_THIS_DIR = Path(__file__).parent
sys.path.insert(0, str(_THIS_DIR))
from importlib import import_module
_gen = import_module("1_anatomical_mri_autism_nt_v2")  # noqa: E402

MRI_FILE = _gen.MRI_FILE
DF_CLUSTERS_FILE = _gen.DF_CLUSTERS_FILE
WAVE_PRIORITY = _gen.WAVE_PRIORITY

_SECTION_DIR = _THIS_DIR.parent
FIG_DIR = _SECTION_DIR / "outputs" / "figures"
TAB_DIR = _SECTION_DIR / "outputs" / "tables"
FIG_DIR.mkdir(parents=True, exist_ok=True)
TAB_DIR.mkdir(parents=True, exist_ok=True)

sns.set_theme(style="white", context="paper")
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["font.family"] = "sans-serif"

# Palette matches the manuscript EEG figures
# (eeg_mri-pipeline/.../01_eeg_alpha_peak_autism_td.py): pink for IDD.
GROUP_COLORS = {
    "Autism":    "#5CAEE1",
    "NT":        "#C1C2BC",
    "IDD":       "#D8A4CB",
    "Relatives": "#9AD5D3",
}
ACCENT = "#8A0201"

POP_ORDER = ["NT", "Autism", "IDD", "Relatives"]
# Rename TD→NT for display, and the dataset's "ID" label to the manuscript's
# "IDD" (intellectual / developmental disability) for figure / table use.
POP_RENAME = {"TD": "NT", "ID": "IDD"}


def _cohens_d(a, b):
    a, b = np.asarray(a), np.asarray(b)
    if len(a) < 2 or len(b) < 2:
        return np.nan
    sd = np.sqrt(((len(a) - 1) * np.var(a, ddof=1)
                  + (len(b) - 1) * np.var(b, ddof=1))
                 / (len(a) + len(b) - 2))
    if sd == 0:
        return np.nan
    return (np.mean(a) - np.mean(b)) / sd


def _despine(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_facecolor("white")


def load_euler_with_population() -> pd.DataFrame:
    """Load the new TSV, apply wave selection, join with df_clusters_complete
    on the canonical key, and return one row per subject with the columns
    needed for the Euler analysis."""
    print(f"Loading {MRI_FILE}")
    df = pd.read_csv(MRI_FILE, sep="\t", low_memory=False,
                     usecols=["ID", "MRI_ID", "cohort",
                              "session_id", "QC_seg", "mean_euler"])
    df["ID"] = df["ID"].astype(str)
    print(f"  raw rows: {len(df)}  cohorts: "
          f"{df['cohort'].value_counts().to_dict()}")

    df = _gen.select_one_row_per_subject(df)
    print(f"  after wave selection: {len(df)} rows")

    df["_join_key"] = df.apply(_gen._canonical_join_key_mri, axis=1)
    df = df.dropna(subset=["_join_key", "mean_euler"])
    print(f"  with join key + Euler: {len(df)} rows")

    clu = pd.read_csv(DF_CLUSTERS_FILE, low_memory=False)
    clu["_join_key"] = clu.apply(_gen._canonical_join_key_clusters, axis=1)
    keep = ["ID", "_join_key", "cohort", "PopulationS1", "Cluster"]
    keep = [c for c in keep if c in clu.columns]
    clu = clu[keep].dropna(subset=["_join_key"])
    print(f"  clusters file: {len(clu)} rows with join key")

    df = df.drop(columns=[c for c in ["ID"] if c in df.columns])
    merged = clu.merge(df, on="_join_key", how="inner", suffixes=("", "_mri"))
    print(f"  merged: {len(merged)} rows")

    merged["PopulationS1"] = merged["PopulationS1"].replace(POP_RENAME)
    merged = merged.dropna(subset=["mean_euler"])
    print(f"\nFinal: {len(merged)} subjects with Euler number")
    print(f"  per-population n: "
          f"{merged['PopulationS1'].value_counts(dropna=False).to_dict()}")
    return merged


def summarise(df: pd.DataFrame) -> pd.DataFrame:
    """Per-population descriptive stats + Welch t vs NT + Cohen's d vs NT."""
    if "NT" not in df["PopulationS1"].unique():
        print("[warn] No NT subjects — vs-NT contrasts skipped.")
        nt_vals = np.array([])
    else:
        nt_vals = df.loc[df["PopulationS1"] == "NT", "mean_euler"].values
    rows = []
    for pop in sorted(df["PopulationS1"].dropna().unique()):
        vals = df.loc[df["PopulationS1"] == pop, "mean_euler"].values
        if len(vals) == 0:
            continue
        row = {
            "population": pop,
            "n": int(len(vals)),
            "mean": float(np.mean(vals)),
            "sd": float(np.std(vals, ddof=1)) if len(vals) > 1 else np.nan,
            "median": float(np.median(vals)),
            "min": float(np.min(vals)),
            "max": float(np.max(vals)),
        }
        if pop != "NT" and len(nt_vals) > 1 and len(vals) > 1:
            t, p = stats.ttest_ind(vals, nt_vals, equal_var=False)
            row["t_vs_NT"] = float(t)
            row["p_vs_NT"] = float(p)
            row["cohens_d_vs_NT"] = _cohens_d(vals, nt_vals)
        else:
            row["t_vs_NT"] = np.nan
            row["p_vs_NT"] = np.nan
            row["cohens_d_vs_NT"] = np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def plot_distribution(df: pd.DataFrame, summary: pd.DataFrame) -> None:
    pops_present = [p for p in POP_ORDER if p in df["PopulationS1"].unique()]
    sub = df[df["PopulationS1"].isin(pops_present)].copy()
    sub["PopulationS1"] = pd.Categorical(sub["PopulationS1"],
                                          categories=pops_present, ordered=True)
    palette = {p: GROUP_COLORS.get(p, "#888888") for p in pops_present}

    fig, ax = plt.subplots(figsize=(7, 4.8))
    sns.violinplot(data=sub, x="PopulationS1", y="mean_euler",
                   order=pops_present, palette=palette,
                   inner=None, linewidth=0, ax=ax, cut=0)
    sns.stripplot(data=sub, x="PopulationS1", y="mean_euler",
                  order=pops_present, color="black",
                  size=1.6, alpha=0.4, jitter=0.25, ax=ax)
    sns.pointplot(data=sub, x="PopulationS1", y="mean_euler",
                  order=pops_present, color=ACCENT,
                  errorbar=("ci", 95), markers="D",
                  linestyles="", ax=ax)
    # Annotate group counts and p (vs NT)
    summ = summary.set_index("population")
    labels = []
    for p in pops_present:
        n = int(summ.loc[p, "n"])
        if p == "NT":
            labels.append(f"NT\n(n={n})")
        else:
            p_val = summ.loc[p].get("p_vs_NT", np.nan)
            d = summ.loc[p].get("cohens_d_vs_NT", np.nan)
            tag = f"\np={p_val:.2g}  d={d:+.2f}" if np.isfinite(p_val) else ""
            labels.append(f"{p}\n(n={n}){tag}")
    ax.set_xticks(range(len(pops_present)))
    ax.set_xticklabels(labels, fontsize=10)
    ax.set_xlabel("")
    ax.set_ylabel("FreeSurfer Euler number  (per subject)", fontsize=12)
    ax.set_title("Cortical surface reconstruction quality by population",
                 fontsize=13)
    ax.tick_params(axis="y", labelsize=10)
    _despine(ax)
    plt.tight_layout()
    out = FIG_DIR / "euler_by_population.pdf"
    plt.savefig(out, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"[ok] {out}")


def plot_by_cohort(df: pd.DataFrame) -> None:
    """Companion figure: Euler × population, split by cohort (LEAP vs INOVAND)."""
    cohorts = [c for c in ("LEAP", "INOVAND") if c in df["cohort"].unique()]
    if not cohorts:
        return
    pops_present = [p for p in POP_ORDER if p in df["PopulationS1"].unique()]
    palette = {p: GROUP_COLORS.get(p, "#888888") for p in pops_present}

    fig, axes = plt.subplots(1, len(cohorts), figsize=(5 * len(cohorts), 4.2),
                             sharey=True)
    if len(cohorts) == 1:
        axes = [axes]
    for ax, coh in zip(axes, cohorts):
        sub = df[(df["cohort"] == coh)
                 & df["PopulationS1"].isin(pops_present)].copy()
        if sub.empty:
            ax.set_title(f"{coh} (no data)")
            _despine(ax)
            continue
        order_here = [p for p in pops_present if p in sub["PopulationS1"].values]
        sns.violinplot(data=sub, x="PopulationS1", y="mean_euler",
                       order=order_here, palette=palette, ax=ax,
                       inner=None, linewidth=0, cut=0)
        sns.stripplot(data=sub, x="PopulationS1", y="mean_euler",
                      order=order_here, color="black",
                      size=1.4, alpha=0.4, jitter=0.25, ax=ax)
        sns.pointplot(data=sub, x="PopulationS1", y="mean_euler",
                      order=order_here, color=ACCENT,
                      errorbar=("ci", 95), markers="D",
                      linestyles="", ax=ax)
        counts = sub["PopulationS1"].value_counts()
        ax.set_xticks(range(len(order_here)))
        ax.set_xticklabels([f"{p}\n(n={int(counts.get(p, 0))})"
                            for p in order_here], fontsize=9)
        ax.set_title(coh, fontsize=12)
        ax.set_xlabel("")
        ax.set_ylabel("Euler number" if coh == cohorts[0] else "", fontsize=11)
        _despine(ax)
    plt.tight_layout()
    out = FIG_DIR / "euler_by_population_by_cohort.pdf"
    plt.savefig(out, dpi=600, bbox_inches="tight")
    plt.close(fig)
    print(f"[ok] {out}")


def main():
    print("=" * 60)
    print("Euler number by population (REVISION)")
    print("=" * 60)
    df = load_euler_with_population()
    df[["ID", "cohort", "PopulationS1", "Cluster", "mean_euler"]].to_csv(
        TAB_DIR / "euler_by_population.csv", index=False)
    print(f"[ok] {TAB_DIR / 'euler_by_population.csv'}")

    summary = summarise(df)
    print("\n=== Summary ===")
    pd.set_option("display.max_colwidth", 30)
    print(summary.to_string(index=False))
    summary.to_csv(TAB_DIR / "euler_by_population_summary.csv", index=False)
    print(f"\n[ok] {TAB_DIR / 'euler_by_population_summary.csv'}")

    plot_distribution(df, summary)
    plot_by_cohort(df)


if __name__ == "__main__":
    main()
