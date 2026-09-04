#!/usr/bin/env python3
"""
Anatomical MRI Autism vs NT — REVISION

v2 of ../scripts/01_anatomical_mri_autism_td_analysis.py, re-run on the
QC + ComBat + age/sex/eTIV-regression z-scored FreeSurfer table:

  /Volumes/Imaging5/EEG_MRI-MF/ALL/results/tabular/anat/
      z_scoring_qc+combat+regression_ndd/output/freesurfer_zscore_qc1_combat_regress.tsv

Conventions match the 5_cluster_anatomical_analysis pipeline:
- Wave selection: one row per subject. LEAP prefers W1 → W2 → W3; INOVAND
  prefers T1 → T2; other cohorts deduped by ID.
- Phenotype labels come from
  eeg_mri-pipeline/results/dataset_paper/dataframes/df_clusters_complete_kmeans.csv
  (canonical Autism / TD label set; TD renamed to NT downstream).
- Canonical join key:
    LEAP    → int(MRI_ID)        (clusters side: int(ID))
    INOVAND → int(MRI_ID)        (clusters side: int(MRI_ID.strip('sub-')))
- Column rename so the existing R plotter naming holds:
    lh_<region>_ThickAvg  → lh_<region>_thickness
    lh_<region>_SurfArea  → lh_<region>_area
    <Region>_Volume_mm3   → <Region>                (aseg)
- Per-region Welch's t-test (Autism vs NT) + Cohen's d + analytic 95 % CI
  (Hedges-Olkin 1985) + FDR-BH.

Outputs (under ../outputs/):
- tables/r_input_files/t_stat_anat_dk_thickness_mri_autism_vs_control.csv
- tables/r_input_files/t_stat_anat_dk_area_mri_autism_vs_control.csv
- tables/r_input_files/t_stat_anat_aseg_volume_mri_autism_vs_control.csv
- tables/autism_vs_nt_summary.csv  — per-metric headline counts
"""

import os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

# =============================================================================
# PATHS
# =============================================================================
import _config  # shared paths + hyperparameters (see _config.py)

MRI_FILE = _config.MRI_FILE
# Phenotype-label file (CLUSTER_FILE env overridable, handled in _config). Only
# PopulationS1 + the join keys are used here (autism is pooled across clusters).
DF_CLUSTERS_FILE = _config.DF_CLUSTERS_FILE

_SCRIPT_DIR = Path(__file__).parent
_SECTION_DIR = _SCRIPT_DIR.parent  # 4_anatomical_analysis/
R_INPUT_DIR = _SECTION_DIR / "outputs" / "tables" / "r_input_files"
TABLES_DIR = _SECTION_DIR / "outputs" / "tables"

DESIKAN_REGIONS = [
    "bankssts", "caudalanteriorcingulate", "caudalmiddlefrontal", "cuneus",
    "entorhinal", "fusiform", "inferiorparietal", "inferiortemporal",
    "isthmuscingulate", "lateraloccipital", "lateralorbitofrontal", "lingual",
    "medialorbitofrontal", "middletemporal", "parahippocampal", "paracentral",
    "parsopercularis", "parsorbitalis", "parstriangularis", "pericalcarine",
    "postcentral", "posteriorcingulate", "precentral", "precuneus",
    "rostralanteriorcingulate", "rostralmiddlefrontal", "superiorfrontal",
    "superiorparietal", "superiortemporal", "supramarginal", "frontalpole",
    "temporalpole", "transversetemporal", "insula",
]

ASEG_REGIONS = [
    "Left-Lateral-Ventricle", "Left-Thalamus", "Left-Caudate", "Left-Putamen",
    "Left-Pallidum", "3rd-Ventricle", "4th-Ventricle", "Brain-Stem",
    "Left-Hippocampus", "Left-Amygdala", "Left-VentralDC",
    "Right-Lateral-Ventricle", "Right-Thalamus", "Right-Caudate", "Right-Putamen",
    "Right-Pallidum", "Right-Hippocampus", "Right-Amygdala", "Right-VentralDC",
    "CC_Posterior", "CC_Mid_Posterior", "CC_Central", "CC_Mid_Anterior",
    "CC_Anterior", "Right-Cerebellum-Cortex", "Right-Cerebellum-White-Matter",
    "Left-Cerebellum-Cortex", "Left-Cerebellum-White-Matter",
]

WAVE_PRIORITY = {
    "LEAP":    ["LEAP_W1", "LEAP_W2", "LEAP_W3"],
    "INOVAND": ["INOVAND_T1", "INOVAND_T2"],
}


# =============================================================================
# LOADERS — identical to 5_cluster_anatomical_analysis/scripts/...
# =============================================================================
def _rename_freesurfer_cols(df: pd.DataFrame) -> pd.DataFrame:
    rename = {}
    for r in DESIKAN_REGIONS:
        for hemi in ("lh", "rh"):
            rename[f"{hemi}_{r}_ThickAvg"] = f"{hemi}_{r}_thickness"
            rename[f"{hemi}_{r}_SurfArea"] = f"{hemi}_{r}_area"
            rename[f"{hemi}_{r}_GrayVol"]  = f"{hemi}_{r}_grayvol"
    for r in ASEG_REGIONS:
        rename[f"{r}_Volume_mm3"] = r
    return df.rename(columns=rename)


def _select_wave_for_cohort(df: pd.DataFrame, prefix: str, priority: list) -> pd.DataFrame:
    sub = df[df["cohort"].astype(str).str.startswith(prefix)].copy()
    if sub.empty:
        return sub
    rank = {w: i for i, w in enumerate(priority)}
    sub["_wave_rank"] = sub["cohort"].map(rank).fillna(len(rank))
    sub = (sub.sort_values(["ID", "_wave_rank"])
              .drop_duplicates(subset=["ID"], keep="first")
              .drop(columns="_wave_rank"))
    return sub


def select_one_row_per_subject(df: pd.DataFrame) -> pd.DataFrame:
    chunks = []
    handled = pd.Series(False, index=df.index)
    for prefix, priority in WAVE_PRIORITY.items():
        chunk = _select_wave_for_cohort(df, prefix, priority)
        if not chunk.empty:
            chunks.append(chunk)
            handled |= df["cohort"].astype(str).str.startswith(prefix)
    other = df.loc[~handled].copy()
    if not other.empty:
        other = other.drop_duplicates(subset=["ID"], keep="first")
        chunks.append(other)
    return pd.concat(chunks, ignore_index=True)


def _canonical_join_key_mri(row) -> str | None:
    raw = row.get("MRI_ID", row.get("ID"))
    try:
        f = float(raw)
        if not np.isfinite(f):
            return None
        return str(int(f))
    except (TypeError, ValueError):
        return None


def _canonical_join_key_clusters(row) -> str | None:
    cohort = str(row.get("cohort", ""))
    if cohort == "LEAP":
        try:
            f = float(row["ID"])
            return str(int(f)) if np.isfinite(f) else None
        except (TypeError, ValueError):
            return None
    if cohort == "INOVAND":
        mri = str(row.get("MRI_ID", ""))
        if mri.startswith("sub-"):
            mri = mri[4:]
        # Accept both "894" (paper file) and "894.0" (curated float-inferred).
        try:
            return str(int(float(mri)))
        except (TypeError, ValueError):
            return None
    return None


def load_data() -> pd.DataFrame:
    print(f"Loading {MRI_FILE}")
    df = pd.read_csv(MRI_FILE, sep="\t", low_memory=False)
    df["ID"] = df["ID"].astype(str)
    df = _rename_freesurfer_cols(df)
    print(f"  raw shape: {df.shape}  cohorts: "
          f"{df['cohort'].value_counts().to_dict()}")

    df = select_one_row_per_subject(df)
    print(f"  after wave selection: {len(df)} rows")

    df["_join_key"] = df.apply(_canonical_join_key_mri, axis=1)
    df = df.dropna(subset=["_join_key"])
    print(f"  with canonical join key: {len(df)} rows")

    # Phenotype source: curated per-cohort clinical TSVs (mirrors 1_clustering)
    # when CURATED_CLINICAL is set, else the DF_CLUSTERS_FILE roster.
    if _config.CURATED_CLINICAL:
        import curated_clinical
        print("  phenotype source: curated clinical TSVs (LEAP+INOVAND+INFOR)")
        # No IQ/SRS requirement for the MRI case-control: keep every scanned +
        # diagnosed subject (769 Autism / 349 NT). Set require_features=True to
        # restrict to the clustering-eligible IQ+SRS-complete cohort (747 / 319).
        clu = curated_clinical.load_all_cohorts(require_features=False)
    else:
        clu = pd.read_csv(DF_CLUSTERS_FILE, low_memory=False)
    clu["_join_key"] = clu.apply(_canonical_join_key_clusters, axis=1)
    # Phenotype source: the curated file's `population_group` is the clean
    # grouping (it folds "Autism to exclude" / "Autism with|without IDD" into
    # "Autism"); the frozen paper file only has `PopulationS1`. Prefer the
    # former when present so curated runs keep those cases.
    pheno_col = "population_group" if "population_group" in clu.columns else "PopulationS1"
    print(f"  phenotype column: {pheno_col}")
    keep = ["ID", "_join_key", "cohort", pheno_col]
    keep = [c for c in keep if c in clu.columns]
    clu = clu[keep].rename(columns={pheno_col: "PopulationS1"}).dropna(subset=["_join_key"])
    print(f"  clusters file with join key: {len(clu)} rows  "
          f"per cohort: {clu['cohort'].value_counts().to_dict()}")

    df = df.drop(columns=[c for c in ["ID"] if c in df.columns])
    merged = clu.merge(df, on="_join_key", how="inner", suffixes=("", "_mri"))
    print(f"  merged shape: {merged.shape}")

    merged["PopulationS1"] = merged["PopulationS1"].replace("TD", "NT")
    keep = merged[merged["PopulationS1"].isin(["Autism", "NT"])].copy()
    print(f"\nFinal Autism vs NT sample: {len(keep)} participants")
    print(f"  Autism: {(keep['PopulationS1']=='Autism').sum()}")
    print(f"  NT    : {(keep['PopulationS1']=='NT').sum()}")
    print(f"  per cohort: {keep['cohort'].value_counts().to_dict()}")
    return keep


# =============================================================================
# STATS
# =============================================================================
def _cohens_d(a, b):
    n_a, n_b = len(a), len(b)
    if n_a < 2 or n_b < 2:
        return np.nan
    sd = np.sqrt(((n_a - 1) * np.var(a, ddof=1) + (n_b - 1) * np.var(b, ddof=1))
                 / (n_a + n_b - 2))
    if sd == 0:
        return np.nan
    return (np.mean(a) - np.mean(b)) / sd


def _cohens_d_ci(d, n_a, n_b):
    if not np.isfinite(d) or n_a < 2 or n_b < 2:
        return np.nan, np.nan
    se = np.sqrt((n_a + n_b) / (n_a * n_b) + d ** 2 / (2.0 * (n_a + n_b)))
    z = 1.959963984540054
    return float(d - z * se), float(d + z * se)


def _run_metric(df: pd.DataFrame, metric: str, atlas_type: str) -> pd.DataFrame:
    """metric in {'thickness','area','grayvol','volume'}; atlas in {'dk','aseg'}."""
    autism = df[df["PopulationS1"] == "Autism"]
    nt = df[df["PopulationS1"] == "NT"]
    rows = []

    if atlas_type == "dk":
        cols = [c for c in df.columns
                if c.endswith(f"_{metric}")
                and (c.startswith("lh_") or c.startswith("rh_"))]
        for c in cols:
            a = autism[c].dropna().values
            b = nt[c].dropna().values
            if len(a) < 2 or len(b) < 2:
                continue
            t_stat, p_val = ttest_ind(a, b, equal_var=False)
            d = _cohens_d(a, b)
            d_lo, d_hi = _cohens_d_ci(d, len(a), len(b))
            region = (c.replace("_thickness", "")
                       .replace("_area", "")
                       .replace("_grayvol", ""))
            # 5_cluster_anatomical_analysis uses lh_<region> as the label;
            # ../scripts/03_plot_anatomical_mri_brain_visualizations.R also
            # converts _left/_right back to lh_/rh_ — pick the lh_/rh_ form
            # so the R script's branch for ggseg works unchanged.
            label = region
            rows.append({
                "label": label, "t_stat": t_stat, "p_val": p_val,
                "cohens_d": d, "cohens_d_lo": d_lo, "cohens_d_hi": d_hi,
                "n_a": int(len(a)), "n_b": int(len(b)),
            })
    else:  # aseg
        for c in ASEG_REGIONS:
            if c not in df.columns:
                continue
            a = autism[c].dropna().values
            b = nt[c].dropna().values
            if len(a) < 2 or len(b) < 2:
                continue
            t_stat, p_val = ttest_ind(a, b, equal_var=False)
            d = _cohens_d(a, b)
            d_lo, d_hi = _cohens_d_ci(d, len(a), len(b))
            rows.append({
                "label": c, "t_stat": t_stat, "p_val": p_val,
                "cohens_d": d, "cohens_d_lo": d_lo, "cohens_d_hi": d_hi,
                "n_a": int(len(a)), "n_b": int(len(b)),
            })

    if not rows:
        return pd.DataFrame()
    res = pd.DataFrame(rows)
    res["p_fdr"] = multipletests(res["p_val"].values, method="fdr_bh")[1]
    res["p_bonf"] = multipletests(res["p_val"].values, method="bonferroni")[1]
    return res


# =============================================================================
# MAIN
# =============================================================================
def main() -> None:
    R_INPUT_DIR.mkdir(parents=True, exist_ok=True)
    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    print("=" * 60)
    print("Anatomical MRI · Autism vs NT (REVISION)")
    print(f"Output : {R_INPUT_DIR}")
    print("=" * 60)

    df = load_data()

    summary = []
    for metric, atlas in [
        ("thickness", "dk"),
        ("area", "dk"),
        ("grayvol", "dk"),    # per-Desikan cortical gray-matter volume
        ("volume", "aseg"),
    ]:
        if metric == "volume":
            res = _run_metric(df, "volume", "aseg")
            fname = "t_stat_anat_aseg_volume_mri_autism_vs_control.csv"
        else:
            res = _run_metric(df, metric, "dk")
            fname = f"t_stat_anat_dk_{metric}_mri_autism_vs_control.csv"

        if res.empty:
            print(f"\n[skip] {metric}: no testable region")
            continue
        out_path = R_INPUT_DIR / fname
        res.to_csv(out_path, index=False)
        n_fdr = int((res["p_fdr"] < 0.05).sum())
        n_bonf = int((res["p_bonf"] < 0.05).sum())
        top = res.iloc[res["cohens_d"].abs().argmax()]
        print(f"\n{metric:10s} ({atlas}): {len(res):3d} regions · "
              f"FDR<0.05 {n_fdr} · Bonferroni<0.05 {n_bonf} · "
              f"top |d| {top['label']} = {top['cohens_d']:+.3f} "
              f"(p_fdr={top['p_fdr']:.4g}, p_bonf={top['p_bonf']:.4g})")
        print(f"  → {out_path.name}")
        summary.append({
            "metric": metric, "atlas": atlas,
            "n_regions": len(res),
            "n_fdr_sig": n_fdr,
            "n_bonf_sig": n_bonf,
            "n_autism": int(res["n_a"].iat[0]),
            "n_nt": int(res["n_b"].iat[0]),
            "top_region": top["label"],
            "top_cohens_d": float(top["cohens_d"]),
            "top_p_fdr": float(top["p_fdr"]),
            "top_p_bonf": float(top["p_bonf"]),
        })

    pd.DataFrame(summary).to_csv(TABLES_DIR / "autism_vs_nt_summary.csv",
                                 index=False)
    print(f"\nWrote {TABLES_DIR / 'autism_vs_nt_summary.csv'}")


if __name__ == "__main__":
    main()
