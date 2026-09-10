#!/usr/bin/env python3
"""
Per-site effect size (Autism vs NT) for left superior temporal thickness.

For the ROI lh_superiortemporal_thickness, computes Cohen's d (Autism - NT)
within each acquisition site, with analytic 95% CI (Hedges-Olkin) and a Welch
t-test, to show whether the cortical-thickness result is consistent across
sites or driven by one. Site comes from the curated clinical `site_id`
(LEAP: STOCKHOLM/KCL/UTRECHT/NIJMEGEN/MANNHEIM/ROME; INOVAND: its hospitals).

Reuses the MRI loader / wave-selection / join keys from 01 and the curated
clinical loader from curated_clinical (so site_id is available).

Outputs (under <base>/):
  tables/site_effect_lh_superiortemporal_thickness.csv
  figures/site_effect_lh_superiortemporal_thickness.pdf   (forest plot)
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
_m01 = import_module("1_anatomical_mri_autism_nt")
import _config
import curated_clinical

REGION = "lh_superiortemporal_thickness"   # left superior temporal gyrus thickness
MIN_PER_GROUP = 3                          # need >= this many Autism AND NT per site
OUT_BASE = Path(_config.OUTPUT_BASE)
TABLES_DIR = OUT_BASE / "tables"
FIGS_DIR = OUT_BASE / "figures"


def load() -> pd.DataFrame:
    df = pd.read_csv(_config.MRI_FILE, sep="\t", low_memory=False)
    df["ID"] = df["ID"].astype(str)
    df = _m01._rename_freesurfer_cols(df)
    df = _m01.select_one_row_per_subject(df)
    df["_join_key"] = df.apply(_m01._canonical_join_key_mri, axis=1)
    df = df.dropna(subset=["_join_key"])

    clu = curated_clinical.load_all_cohorts(require_features=False)
    clu["_join_key"] = clu.apply(_m01._canonical_join_key_clusters, axis=1)
    pheno = "population_group" if "population_group" in clu.columns else "PopulationS1"
    keep = ["_join_key", "cohort", "site_id", pheno]
    keep = [c for c in keep if c in clu.columns]
    clu = clu[keep].rename(columns={pheno: "PopulationS1"}).dropna(subset=["_join_key"])

    df = df.drop(columns=[c for c in ["ID"] if c in df.columns])
    m = clu.merge(df, on="_join_key", how="inner", suffixes=("", "_mri"))
    m["PopulationS1"] = m["PopulationS1"].replace("TD", "NT")
    m = m[m["PopulationS1"].isin(["Autism", "NT"])]
    m = m.dropna(subset=[REGION])
    m["site_id"] = m["site_id"].fillna("Unknown").astype(str)
    print(f"Merged Autism/NT with {REGION}: {len(m)}  "
          f"(Autism {int((m.PopulationS1=='Autism').sum())}, "
          f"NT {int((m.PopulationS1=='NT').sum())})")
    return m


def effect(autism, nt):
    a = autism.dropna().values
    b = nt.dropna().values
    if len(a) < 2 or len(b) < 2:
        return None
    d = _m01._cohens_d(a, b)
    lo, hi = _m01._cohens_d_ci(d, len(a), len(b))
    t, p = ttest_ind(a, b, equal_var=False)
    return dict(cohens_d=d, ci_lo=lo, ci_hi=hi, t=float(t), p=float(p),
                n_autism=len(a), n_nt=len(b))


def main():
    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    FIGS_DIR.mkdir(parents=True, exist_ok=True)
    m = load()

    rows = []
    # Overall (all sites pooled) for reference.
    ov = effect(m.loc[m.PopulationS1 == "Autism", REGION],
                m.loc[m.PopulationS1 == "NT", REGION])
    if ov:
        rows.append({"site": "ALL (pooled)", "cohort": "all", **ov})

    for site, sub in m.groupby("site_id"):
        e = effect(sub.loc[sub.PopulationS1 == "Autism", REGION],
                   sub.loc[sub.PopulationS1 == "NT", REGION])
        if e is None or e["n_autism"] < MIN_PER_GROUP or e["n_nt"] < MIN_PER_GROUP:
            continue
        cohort = sub["cohort"].mode().iat[0] if "cohort" in sub else ""
        rows.append({"site": site, "cohort": cohort, **e})

    res = pd.DataFrame(rows)
    out_csv = TABLES_DIR / "site_effect_lh_superiortemporal_thickness.csv"
    res.to_csv(out_csv, index=False)
    print(f"\nWrote {out_csv}\n")
    print(res.to_string(index=False,
          columns=["site", "cohort", "n_autism", "n_nt", "cohens_d",
                   "ci_lo", "ci_hi", "p"], float_format=lambda x: f"{x:.3f}"))

    # --- Forest plot ---------------------------------------------------------
    plot_df = res[res.site != "ALL (pooled)"].sort_values("cohens_d")
    fig, ax = plt.subplots(figsize=(7, 0.5 * len(plot_df) + 1.5))
    y = np.arange(len(plot_df))
    ax.errorbar(plot_df["cohens_d"], y,
                xerr=[plot_df["cohens_d"] - plot_df["ci_lo"],
                      plot_df["ci_hi"] - plot_df["cohens_d"]],
                fmt="o", color="#324095", capsize=3, lw=1)
    ax.axvline(0, color="black", lw=0.6)
    if ov:
        ax.axvline(ov["cohens_d"], color="grey", ls="--", lw=0.8,
                   label=f"pooled d = {ov['cohens_d']:.2f}")
        ax.legend(loc="best", fontsize=8)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r.site} (n={r.n_autism}/{r.n_nt})"
                        for r in plot_df.itertuples()], fontsize=8)
    ax.set_xlabel("Cohen's d (Autism - NT), 95% CI")
    ax.set_title("Left superior temporal thickness — effect size by site")
    plt.tight_layout()
    out_pdf = FIGS_DIR / "site_effect_lh_superiortemporal_thickness.pdf"
    plt.savefig(out_pdf, bbox_inches="tight")
    print(f"\nWrote {out_pdf}")

    # --- Lollipop plot (point estimate only, no CI whiskers) -----------------
    fig2, ax2 = plt.subplots(figsize=(7, 0.5 * len(plot_df) + 1.5))
    colors = ["#324095" if p < 0.05 else "#9AA0C0" for p in plot_df["p"]]
    ax2.hlines(y, 0, plot_df["cohens_d"], color=colors, lw=2, zorder=1)
    ax2.scatter(plot_df["cohens_d"], y, s=80, color=colors, zorder=2,
               edgecolors="black", linewidths=0.5)
    ax2.axvline(0, color="black", lw=0.6)
    if ov:
        ax2.axvline(ov["cohens_d"], color="grey", ls="--", lw=0.8,
                   label=f"pooled d = {ov['cohens_d']:.2f}")
        ax2.legend(loc="best", fontsize=8)
    ax2.set_yticks(y)
    ax2.set_yticklabels([f"{r.site} (n={r.n_autism}/{r.n_nt})"
                        for r in plot_df.itertuples()], fontsize=8)
    ax2.set_xlabel("Cohen's d (Autism - NT)")
    ax2.set_title("Left superior temporal thickness — effect size by site (lollipop)")
    for s in ("top", "right"): ax2.spines[s].set_visible(False)
    plt.tight_layout()
    out_pdf2 = FIGS_DIR / "site_effect_lh_superiortemporal_thickness_lollipop.pdf"
    plt.savefig(out_pdf2, bbox_inches="tight")
    print(f"Wrote {out_pdf2}")


def main_by_scanner():
    """Same effect-size breakdown, but by scanner (the `batch` column -- the
    same scanner-level ComBat batch label used throughout this repo, e.g.
    22_anat_sequence_combat_svm.py) rather than by site_id. This splits
    INOVAND (one site_id, several scanner generations: Intera/Ingenia/Ingenia
    3T) and shows LEAP at the same granularity it already had (its batch is
    one scanner per site)."""
    df = _m01.load_data()
    df = df.dropna(subset=[REGION, "batch"])

    rows = []
    ov = effect(df.loc[df.PopulationS1 == "Autism", REGION],
                df.loc[df.PopulationS1 == "NT", REGION])
    if ov:
        rows.append({"scanner": "ALL (pooled)", **ov})
    for scanner, sub in df.groupby("batch"):
        e = effect(sub.loc[sub.PopulationS1 == "Autism", REGION],
                   sub.loc[sub.PopulationS1 == "NT", REGION])
        if e is None or e["n_autism"] < MIN_PER_GROUP or e["n_nt"] < MIN_PER_GROUP:
            continue
        rows.append({"scanner": scanner, **e})

    res = pd.DataFrame(rows)
    out_csv = TABLES_DIR / "scanner_effect_lh_superiortemporal_thickness.csv"
    res.to_csv(out_csv, index=False)
    print(f"\nWrote {out_csv}\n")
    print(res.to_string(index=False,
          columns=["scanner", "n_autism", "n_nt", "cohens_d", "ci_lo", "ci_hi", "p"],
          float_format=lambda x: f"{x:.3f}"))

    plot_df = res[res.scanner != "ALL (pooled)"].sort_values("cohens_d")
    y = np.arange(len(plot_df))

    # --- Forest plot (matches main()'s style) --------------------------------
    fig0, ax0 = plt.subplots(figsize=(7, 0.5 * len(plot_df) + 1.5))
    ax0.errorbar(plot_df["cohens_d"], y,
                xerr=[plot_df["cohens_d"] - plot_df["ci_lo"],
                      plot_df["ci_hi"] - plot_df["cohens_d"]],
                fmt="o", color="#324095", capsize=3, lw=1)
    ax0.axvline(0, color="black", lw=0.6)
    if ov:
        ax0.axvline(ov["cohens_d"], color="grey", ls="--", lw=0.8,
                   label=f"pooled d = {ov['cohens_d']:.2f}")
        ax0.legend(loc="best", fontsize=8)
    ax0.set_yticks(y)
    ax0.set_yticklabels([f"{r.scanner} (n={r.n_autism}/{r.n_nt})"
                        for r in plot_df.itertuples()], fontsize=8)
    ax0.set_xlabel("Cohen's d (Autism - NT), 95% CI")
    ax0.set_title("Left superior temporal thickness — effect size by scanner")
    plt.tight_layout()
    out_pdf0 = FIGS_DIR / "scanner_effect_lh_superiortemporal_thickness.pdf"
    plt.savefig(out_pdf0, bbox_inches="tight")
    print(f"\nWrote {out_pdf0}")

    # --- Lollipop plot (kept alongside, in case still useful) ---------------
    fig, ax = plt.subplots(figsize=(7, 0.5 * len(plot_df) + 1.5))
    colors = ["#7a0010" if p < 0.05 else "#D9A0A8" for p in plot_df["p"]]
    ax.hlines(y, 0, plot_df["cohens_d"], color=colors, lw=2, zorder=1)
    ax.scatter(plot_df["cohens_d"], y, s=80, color=colors, zorder=2,
              edgecolors="black", linewidths=0.5)
    ax.axvline(0, color="black", lw=0.6)
    if ov:
        ax.axvline(ov["cohens_d"], color="grey", ls="--", lw=0.8,
                  label=f"pooled d = {ov['cohens_d']:.2f}")
        ax.legend(loc="best", fontsize=8)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r.scanner} (n={r.n_autism}/{r.n_nt})"
                       for r in plot_df.itertuples()], fontsize=8)
    ax.set_xlabel("Cohen's d (Autism - NT)")
    ax.set_title("Left superior temporal thickness — effect size by scanner (lollipop)")
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    plt.tight_layout()
    out_pdf = FIGS_DIR / "scanner_effect_lh_superiortemporal_thickness_lollipop.pdf"
    plt.savefig(out_pdf, bbox_inches="tight")
    print(f"\nWrote {out_pdf}")


if __name__ == "__main__":
    main()
    main_by_scanner()
