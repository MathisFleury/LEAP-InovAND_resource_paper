#!/usr/bin/env python3
"""
Cluster MRI Input Generation — REVISION (Reviewer 3)

Re-runs the per-cluster anatomical-MRI contrasts (CT, SA, subcortical
volume) against pooled NT controls and the pairwise cluster comparisons,
using the **new QC + ComBat + age/sex/eTIV-regressed z-scored** FreeSurfer
table:

  /Volumes/Imaging5/EEG_MRI-MF/ALL/results/tabular/anat/
      z_scoring_qc+combat+regression_ndd/output/freesurfer_zscore_qc1_combat_regress.tsv

Differences vs. ../scripts/01_generate_cluster_mri_inputs.py:

- New input file (better QC, longitudinal-aware, ComBat-harmonised across
  cohorts).
- Per-cohort wave selection: one row per subject, preferring the earliest
  wave and falling back when missing.
    LEAP:    LEAP_W1   → LEAP_W2   → LEAP_W3
    InovAND: INOVAND_T1 → INOVAND_T2
  All other cohorts (NT-source datasets like ABIDE, BRAINMAP, SEARCHLITE,
  CANDY_PIP) are de-duplicated by ID, keeping the first observed row.
- Column naming rewrite so the downstream R brain plotter
  (../scripts/02_plot_cluster_mri_brain.R) can be reused without changes:
    lh_<region>_ThickAvg  → lh_<region>_thickness
    lh_<region>_SurfArea  → lh_<region>_area
    <Region>_Volume_mm3   → <Region>            (aseg)
- Cohen's d + analytic 95 % CI (Hedges-Olkin 1985) — same convention as
  the existing pipeline.

Outputs are written to:
  5_cluster_anatomical_analysis/revision/outputs/figures/r_input_files/
"""

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

_parser = argparse.ArgumentParser(description=__doc__)
_parser.add_argument(
    "--nt-c1", action="store_true",
    help="Restrict the NT control pool to NT subjects whose own Cluster == 'C1' "
         "(the 'NT-like' cluster). Mirrors the _nt_c1 variant used by the "
         "cluster-genetic scripts in 3_cluster_genetic_analysis.",
)
_args, _ = _parser.parse_known_args()
NT_C1_ONLY = _args.nt_c1

# =============================================================================
# PATHS
# =============================================================================
DF_CLUSTERS_FILE = (
    "/Users/mfleury/POSTDOC/LIBRAIRY/LEAP-InovAND_resource/"
    "1_clustering/outputs/tables/individuals_metrics_with_clusters.csv"
)
MRI_FILE = (
    "/Volumes/Imaging5/EEG_MRI-MF/ALL/results/tabular/anat/"
    "z_scoring_qc+combat+regression_ndd/output/freesurfer_zscore_qc1_combat_regress.tsv"
)
_SCRIPT_DIR = Path(__file__).parent
_SECTION_DIR = _SCRIPT_DIR.parent  # = 5_cluster_anatomical_analysis/revision/
OUTPUT_DIR = _SECTION_DIR / "outputs" / "figures" / "r_input_files"

# =============================================================================
# REGION LABEL MAPS
# (Re-used by the R brain plotter; we strip the metric-specific suffix from
# the original column names so the labels match the existing data_dk/data_aseg
# dictionaries in 02_plot_cluster_mri_brain.R.)
# =============================================================================
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

# Wave preference per cohort family.  The keys are *prefixes*: any cohort
# string starting with the prefix is grouped together for wave-selection.
WAVE_PRIORITY = {
    "LEAP":    ["LEAP_W1", "LEAP_W2", "LEAP_W3"],
    "INOVAND": ["INOVAND_T1", "INOVAND_T2"],
}


# =============================================================================
# LOADER
# =============================================================================
def load_freesurfer_table() -> pd.DataFrame:
    """Load the new TSV and rename ROI columns to match the canonical schema."""
    print(f"Loading: {MRI_FILE}")
    df = pd.read_csv(MRI_FILE, sep="\t", low_memory=False)
    df["ID"] = df["ID"].astype(str)
    print(f"  raw shape: {df.shape}   cohorts: "
          f"{df['cohort'].value_counts().to_dict()}")

    rename = {}
    for r in DESIKAN_REGIONS:
        for hemi in ("lh", "rh"):
            rename[f"{hemi}_{r}_ThickAvg"] = f"{hemi}_{r}_thickness"
            rename[f"{hemi}_{r}_SurfArea"] = f"{hemi}_{r}_area"
            rename[f"{hemi}_{r}_GrayVol"]  = f"{hemi}_{r}_grayvol"
    for r in ASEG_REGIONS:
        rename[f"{r}_Volume_mm3"] = r
    df = df.rename(columns=rename)
    n_th = sum(1 for c in df.columns if c.endswith("_thickness"))
    n_ar = sum(1 for c in df.columns if c.endswith("_area"))
    n_as = sum(1 for c in df.columns if c in ASEG_REGIONS)
    print(f"  after rename: {n_th} thickness · {n_ar} area · {n_as} aseg cols")
    return df


def _select_wave_for_cohort(df: pd.DataFrame, prefix: str,
                            priority: list[str]) -> pd.DataFrame:
    """One row per ID for cohort strings starting with ``prefix``.

    Picks the earliest wave from ``priority`` that the subject has data for.
    Returns the filtered subset (only cohort==prefix... rows).
    """
    sub = df[df["cohort"].astype(str).str.startswith(prefix)].copy()
    if sub.empty:
        return sub
    n_before = len(sub)
    n_unique = sub["ID"].nunique()
    print(f"\n{prefix} rows: {n_before}  ({n_unique} unique IDs)")
    rank = {w: i for i, w in enumerate(priority)}
    sub["_wave_rank"] = sub["cohort"].map(rank).fillna(len(rank))
    sub = (sub.sort_values(["ID", "_wave_rank"])
              .drop_duplicates(subset=["ID"], keep="first")
              .drop(columns="_wave_rank"))
    chosen = sub["cohort"].value_counts().to_dict()
    print(f"  after wave selection: {len(sub)} rows · per-wave kept: {chosen}")
    return sub


def select_one_row_per_subject(df: pd.DataFrame) -> pd.DataFrame:
    """Apply wave-selection per cohort; de-duplicate everything else by ID."""
    df = df.copy()
    handled_mask = pd.Series(False, index=df.index)
    chunks: list[pd.DataFrame] = []
    for prefix, priority in WAVE_PRIORITY.items():
        chunk = _select_wave_for_cohort(df, prefix, priority)
        if not chunk.empty:
            chunks.append(chunk)
            # Mark every row whose cohort starts with this prefix as handled,
            # even those de-duplicated away (so they are not duplicated later).
            handled_mask |= df["cohort"].astype(str).str.startswith(prefix)
    other = df.loc[~handled_mask].copy()
    if not other.empty:
        print(f"\nOther cohorts (deduped by ID): "
              f"{other['cohort'].value_counts().to_dict()}")
        other = other.drop_duplicates(subset=["ID"], keep="first")
        print(f"  after dedup: {len(other)} rows")
        chunks.append(other)
    return pd.concat(chunks, ignore_index=True)


def _canonical_join_key_mri(row) -> str | None:
    """Build a canonical join key from a *new-TSV* row.

    LEAP    → int(MRI_ID), e.g. '100693509718'
    INOVAND → int(MRI_ID), e.g. '894'    (clusters side has 'sub-0894')
    others  → None (not used for cluster matching)
    """
    raw = row.get("MRI_ID", row.get("ID"))
    try:
        f = float(raw)
        if not np.isfinite(f):
            return None
        return str(int(f))
    except (TypeError, ValueError):
        return None


def _canonical_join_key_clusters(row) -> str | None:
    """Build a canonical join key from a *df_clusters_complete* row.

    LEAP    → int(ID), e.g. '700148775290'
    INOVAND → int(strip 'sub-' from MRI_ID), e.g. '894'
    """
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
        # Accept both "894" (clean int string) and "894.0" (pandas float-inferred).
        try:
            return str(int(float(mri)))
        except (TypeError, ValueError):
            return None
    return None


def load_and_prepare_data() -> pd.DataFrame:
    df_mri = load_freesurfer_table()
    df_mri = select_one_row_per_subject(df_mri)
    df_mri["_join_key"] = df_mri.apply(_canonical_join_key_mri, axis=1)
    df_mri = df_mri.dropna(subset=["_join_key"])
    print(f"\nMRI rows with a canonical join key: {len(df_mri)}")

    df_clu = pd.read_csv(DF_CLUSTERS_FILE, low_memory=False)
    df_clu["_join_key"] = df_clu.apply(_canonical_join_key_clusters, axis=1)
    keep = [c for c in ["ID", "_join_key", "cohort", "Cluster",
                        "PopulationS1", "Population1"]
            if c in df_clu.columns]
    df_clu = df_clu[keep].dropna(subset=["_join_key"])
    print(f"Cluster rows with a canonical join key: {len(df_clu)}")
    print(f"  per cohort: {df_clu['cohort'].value_counts(dropna=False).to_dict()}")

    # The MRI file's `cohort` clashes with the clusters' `cohort` after merge;
    # we keep the clusters version as authoritative.
    df_mri = df_mri.drop(columns=[c for c in ["ID"] if c in df_mri.columns])

    print(f"\nMerging on canonical join key (inner)...")
    df = df_clu.merge(df_mri, on="_join_key", how="inner",
                      suffixes=("", "_mri"))
    print(f"  merged shape: {df.shape}")
    print(f"  per source cohort: {df['cohort'].value_counts().to_dict()}")

    # QC is handled upstream (during FreeSurfer ComBat/regression z-scoring) —
    # the per-script LEAP-only blacklist filter has been removed.

    # Keep clustered subjects AND the NT pool (in autism-only clustering mode,
    # NT subjects have NaN Cluster but are still needed for cluster-vs-NT stats).
    df = df[df["Cluster"].notna() | (df["PopulationS1"] == "NT")
            | (df["PopulationS1"] == "TD")]
    df["PopulationS1"] = df["PopulationS1"].replace("TD", "NT")
    print(f"  after cluster filter + label rename: {df.shape}")
    print(f"  PopulationS1 vc: {df['PopulationS1'].value_counts().to_dict()}")
    print(f"  Cluster vc (Autism): "
          f"{df.loc[df['PopulationS1']=='Autism','Cluster'].value_counts().to_dict()}")
    return df


# =============================================================================
# STATS HELPERS
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


def _run_tests(group_a, group_b, columns, label_map):
    rows = []
    for col in columns:
        if col not in group_a.columns or col not in group_b.columns:
            continue
        a = group_a[col].dropna().values
        b = group_b[col].dropna().values
        if len(a) < 2 or len(b) < 2:
            continue
        t_stat, p_val = ttest_ind(a, b, equal_var=False)
        d = _cohens_d(a, b)
        d_lo, d_hi = _cohens_d_ci(d, len(a), len(b))
        rows.append({
            "label": label_map.get(col, col),
            "t_stat": t_stat,
            "p_val": p_val,
            "cohens_d": d,
            "cohens_d_lo": d_lo,
            "cohens_d_hi": d_hi,
            "n_a": int(len(a)),
            "n_b": int(len(b)),
        })
    if not rows:
        return pd.DataFrame()
    out = pd.DataFrame(rows)
    out["p_fdr"] = multipletests(out["p_val"].values, method="fdr_bh")[1]
    return out


# =============================================================================
# COMPARISONS
# =============================================================================
def cluster_vs_nt(df: pd.DataFrame, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    nt = df[df["PopulationS1"] == "NT"]
    if NT_C1_ONLY:
        before = len(nt)
        nt = nt[nt["Cluster"] == "C1"]
        print(f"\n[--nt-c1] NT pool restricted to NT subjects with Cluster == 'C1' "
              f"({before} → {len(nt)})")
    print(f"\nNT pool: {len(nt)} participants")

    thickness_cols = [c for c in df.columns
                      if c.endswith("_thickness")
                      and (c.startswith("lh_") or c.startswith("rh_"))]
    area_cols = [c for c in df.columns
                 if c.endswith("_area")
                 and (c.startswith("lh_") or c.startswith("rh_"))]
    grayvol_cols = [c for c in df.columns
                    if c.endswith("_grayvol")
                    and (c.startswith("lh_") or c.startswith("rh_"))]
    volume_cols = [c for c in df.columns if c in ASEG_REGIONS]

    metrics = [
        ("thickness", "dk", thickness_cols,
         {c: c.replace("_thickness", "") for c in thickness_cols}),
        ("area", "dk", area_cols,
         {c: c.replace("_area", "") for c in area_cols}),
        ("grayvol", "dk", grayvol_cols,
         {c: c.replace("_grayvol", "") for c in grayvol_cols}),
        ("volume", "aseg", volume_cols, {c: c for c in volume_cols}),
    ]

    counts = []
    for cluster in ["C1", "C2", "C3"]:
        sub = df[df["Cluster"] == cluster]
        autism = sub[sub["PopulationS1"] == "Autism"]
        counts.append({"Cluster": cluster,
                       "N_autism": int(len(autism)),
                       "N_total": int(len(sub))})
    pd.DataFrame(counts).to_csv(output_dir / "cluster_counts_mri.csv", index=False)
    print(f"  Wrote cluster_counts_mri.csv  — counts: {counts}")

    for cluster in ["C1", "C2", "C3"]:
        autism = df[(df["Cluster"] == cluster) & (df["PopulationS1"] == "Autism")]
        print(f"\nCluster {cluster}: n autism = {len(autism)}")
        for metric_short, atlas, cols, lmap in metrics:
            if not cols:
                continue
            res = _run_tests(autism, nt, cols, lmap)
            if res.empty:
                continue
            fname = (f"t_stat_cluster_{cluster}_{metric_short}_{atlas}"
                     f"_mri_cluster_vs_td.csv")
            res.to_csv(output_dir / fname, index=False)
            sig = int((res["p_fdr"] < 0.05).sum())
            print(f"  {metric_short:>10s}: {len(res):3d} regions · "
                  f"FDR<0.05: {sig}")


def pairwise(df: pd.DataFrame, output_dir: Path) -> None:
    pairs = [("C1", "C2"), ("C1", "C3"), ("C2", "C3")]
    thickness_cols = [c for c in df.columns
                      if c.endswith("_thickness")
                      and (c.startswith("lh_") or c.startswith("rh_"))]
    area_cols = [c for c in df.columns
                 if c.endswith("_area")
                 and (c.startswith("lh_") or c.startswith("rh_"))]
    grayvol_cols = [c for c in df.columns
                    if c.endswith("_grayvol")
                    and (c.startswith("lh_") or c.startswith("rh_"))]
    volume_cols = [c for c in df.columns if c in ASEG_REGIONS]
    metrics = [
        ("thickness", "dk", thickness_cols,
         {c: c.replace("_thickness", "") for c in thickness_cols}),
        ("area", "dk", area_cols,
         {c: c.replace("_area", "") for c in area_cols}),
        ("grayvol", "dk", grayvol_cols,
         {c: c.replace("_grayvol", "") for c in grayvol_cols}),
        ("volume", "aseg", volume_cols, {c: c for c in volume_cols}),
    ]
    for c1, c2 in pairs:
        a = df[(df["Cluster"] == c1) & (df["PopulationS1"] == "Autism")]
        b = df[(df["Cluster"] == c2) & (df["PopulationS1"] == "Autism")]
        print(f"\nPairwise {c1} vs {c2}: n={len(a)} vs n={len(b)}")
        for metric_short, atlas, cols, lmap in metrics:
            if not cols:
                continue
            res = _run_tests(a, b, cols, lmap)
            if res.empty:
                continue
            fname = f"pairwise_{c1}_vs_{c2}_{metric_short}_{atlas}.csv"
            res.to_csv(output_dir / fname, index=False)
            sig = int((res["p_fdr"] < 0.05).sum())
            print(f"  {metric_short:>10s}: {len(res):3d} regions · "
                  f"FDR<0.05: {sig}")


def main() -> None:
    print("=" * 60)
    print("Cluster MRI Input Generation — REVISION (Reviewer 3)")
    print(f"Output : {OUTPUT_DIR}")
    print("=" * 60)
    df = load_and_prepare_data()
    cluster_vs_nt(df, OUTPUT_DIR)
    pairwise(df, OUTPUT_DIR)
    print("\nDone.")


if __name__ == "__main__":
    main()
