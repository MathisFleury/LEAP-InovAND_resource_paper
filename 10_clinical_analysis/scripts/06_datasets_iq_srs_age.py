#!/usr/bin/env python3
# =============================================================================
# 06 - Cross-dataset IQ / SRS / Age by sex (LEAP-InovAND vs SPARK vs ABIDE)
# =============================================================================
# Port of eeg_mri-pipeline/.../fig_1_clinics_datasets.py (plot_iq_srs_v3):
# split-violin (Male/Female) of full-scale IQ, SRS-2 t-score and age for the
# three autism cohorts.
#
# LEAP-InovAND is taken from the CURATED cohort (autism only). SPARK and ABIDE
# live on external mounts (ghfc_wgs / Imaging5-ABIDE); if a mount is absent the
# cohort is skipped with a warning and the figure is drawn from what is present.
#
# Output: figures/fig_1_iq_srs_age_by_sex.pdf (+ .png)
# =============================================================================

import os
import sys
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _config import CLUSTERS_CURATED_FILE, FIGURES_DIR

# SPARK + ABIDE now live under the mounted IMG5 volume.
IMG5 = "/Volumes/Imaging5/EEG_MRI-MF"
_SPARK_RELEASES = sorted(glob.glob(f"{IMG5}/SPARK/_clinical_data/raw/SPARKDataRelease_*"))
SPARK_DIR = _SPARK_RELEASES[-1] if _SPARK_RELEASES else None   # latest release
ABIDE1 = f"{IMG5}/ABIDE/_clinical_data/raw/Phenotypic_V1_0b.csv"
ABIDE2 = f"{IMG5}/ABIDE2/_clinical_data/raw/ABIDEII_Composite_Phenotypic.csv"

COHORT_ORDER = ["LEAP-InovAND", "SPARK", "ABIDE"]
SEX_PALETTE = {"Male": "#927FB7", "Female": "#EDDED9"}
COLS = ["ID", "Cohort", "total_IQ", "SRS_tscore", "age_yrs", "Sex"]


def load_leap_inovand():
    """Curated autism cohort — ALL autistic individuals (incl. 'Autism to
    exclude', which still carry valid IQ/SRS/age), using the project IQ
    (total_IQ with performance_IQ fallback) for full coverage. Sex is coded
    1=Male / 0=Female here (curated). ponytail: flip if the codebook disagrees."""
    df = pd.read_csv(CLUSTERS_CURATED_FILE, low_memory=False).drop_duplicates("ID")
    df = df[df["PopulationS1"].isin(["Autism", "Autism to exclude"])].copy()
    df["total_IQ"] = df["IQ"]                      # project IQ (fallback) — max coverage
    df["Sex"] = df["Sex"].map({1: "Male", 0: "Female"})
    df["Cohort"] = "LEAP-InovAND"
    for c in ("total_IQ", "SRS_tscore"):
        df[c] = df[c].replace(999, np.nan)
    return df.reindex(columns=COLS)


def load_spark():
    if SPARK_DIR is None or not os.path.isdir(SPARK_DIR):
        print("SKIP SPARK — no release dir found under IMG5")
        return None
    date = os.path.basename(SPARK_DIR).split("_")[-1]   # e.g. 2026-03-23
    iq = pd.read_csv(f"{SPARK_DIR}/iq-{date}.csv")
    srs = pd.read_csv(f"{SPARK_DIR}/srs-2_adult_self_report-{date}.csv")
    pheno = pd.read_csv(f"{SPARK_DIR}/individuals_registration-{date}.csv")
    # keep both scores per subject (outer merge), then attach phenotype
    s = (iq[["subject_sp_id", "fsiq_score"]]
         .merge(srs[["subject_sp_id", "total_t_score"]], on="subject_sp_id", how="outer")
         .merge(pheno, on="subject_sp_id", how="inner"))
    s = s[s["asd"] == True]
    s = s.rename(columns={"subject_sp_id": "ID", "fsiq_score": "total_IQ",
                          "total_t_score": "SRS_tscore",
                          "age_at_registration_years": "age_yrs"})
    s["Sex"] = s["sex"]  # already 'Male'/'Female'
    s["Cohort"] = "SPARK"
    print(f"SPARK release {date}: {len(s)} ASD subjects")
    return s.reindex(columns=COLS)


def load_abide():
    if not (os.path.isfile(ABIDE1) and os.path.isfile(ABIDE2)):
        print(f"SKIP ABIDE — not mounted ({ABIDE1})")
        return None
    a = pd.concat([pd.read_csv(ABIDE1), pd.read_csv(ABIDE2, sep=";")])
    a = a[a["DX_GROUP"] == 1].rename(columns={
        "SUB_ID": "ID", "FIQ": "total_IQ", "SRS_TOTAL_T": "SRS_tscore",
        "SEX": "Sex", "AGE_AT_SCAN": "age_yrs"})
    a["SRS_tscore"] = a["SRS_tscore"].replace(-9999, np.nan)
    a["total_IQ"] = a["total_IQ"].replace(-9999, np.nan)
    a["Sex"] = a["Sex"].map({1: "Male", 2: "Female"})   # ABIDE coding
    a["Cohort"] = "ABIDE"
    return a.reindex(columns=COLS)


def plot(df):
    order = [c for c in COHORT_ORDER if (df["Cohort"] == c).any()]

    def panel(ax, y, title):
        d = df.dropna(subset=[y, "Sex"])
        sns.violinplot(data=d, x="Cohort", y=y, hue="Sex", order=order,
                       palette=SEX_PALETTE, split=True, inner="quartile",
                       linewidth=1.2, cut=0, ax=ax)
        ax.set_title(title, fontsize=15, weight="bold")
        ax.set_xlabel(""); ax.set_ylabel("")
        n = d.groupby("Cohort", observed=True)["ID"].count().reindex(order).fillna(0).astype(int)
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels([f"{c}\n(n={n[c]})" for c in order], fontsize=11)
        ax.tick_params(axis="y", labelsize=11)
        ax.grid(axis="y", linestyle="--", alpha=0.5)
        if ax.legend_:
            ax.legend_.remove()

    fig, axes = plt.subplots(1, 3, figsize=(15, 7))
    axes[0].axhline(70, color="black", ls="--", lw=0.5)
    axes[0].axhline(130, color="black", ls="--", lw=0.5)
    axes[1].axhline(60, color="black", ls="--", lw=0.5)
    axes[1].axhline(75, color="black", ls="--", lw=0.5)
    panel(axes[0], "total_IQ", "Full-scale IQ")
    panel(axes[1], "SRS_tscore", "SRS-2 t-score")
    panel(axes[2], "age_yrs", "Age (years)")

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False, fontsize=12)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    sns.despine(trim=True)

    os.makedirs(FIGURES_DIR, exist_ok=True)
    out = os.path.join(FIGURES_DIR, "fig_1_iq_srs_age_by_sex.pdf")
    fig.savefig(out, dpi=300, bbox_inches="tight")
    fig.savefig(out.replace(".pdf", ".png"), dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out}  cohorts={order}")


def main():
    parts = [load_leap_inovand(), load_spark(), load_abide()]
    df = pd.concat([p for p in parts if p is not None], ignore_index=True).drop_duplicates("ID")
    print("N per cohort:", df["Cohort"].value_counts().to_dict())
    plot(df)


if __name__ == "__main__":
    main()
