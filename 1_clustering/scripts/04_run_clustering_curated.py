#!/usr/bin/env python3
# =============================================================================
# 04 - Run Clustering (k=3) — UP-TO-DATE CURATED CLINICAL DATA
# =============================================================================
# Mirror of 04_run_clustering.py but reads from the curated per-cohort
# clinical TSVs on Imaging5 (see eeg_mri-pipeline/CLAUDE.md):
#
#   LEAP:    /Volumes/Imaging5/EEG_MRI-MF/LEAP/_clinical_data/curated/LEAP_clinical_curated_t1.tsv
#   INOVAND: /Volumes/Imaging5/EEG_MRI-MF/INOVAND/_clinical_data/curated/INOVAND_clinical_curated.tsv
#   INFOR:   /Volumes/Imaging5/EEG_MRI-MF/INFOR/_clinical_data/curated/INFOR_clinical_curated.tsv
#
# INFOR is collapsed into the INOVAND cohort. Features are standardised
# IQ × SRS-2 (IQ = total_IQ fillna performance_IQ), restricted to
# `relation_to_proposant == "participant"`, deduplicated by ID, with no
# EXCLUDED_IDS filter (since the original 3 IDs only mattered for matching
# the paper's frozen df_clusters_complete.csv; this script is a fresh
# clustering on the new data).
#
# Same methodology as 04_run_clustering.py:
#   - Primary: K-means (n_init=25, random_state=42).
#   - Sensitivity: Ward (op.heatmap convention) and GMM (covariance_type="full").
#   - Labels relabeled by descending cluster size.
#
# Outputs go to 1_clustering/outputs/curated/ to avoid overwriting the
# paper-aligned outputs.
# =============================================================================

import argparse
import os
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist, squareform
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings(
    "ignore",
    message=".*looks suspiciously like an uncondensed distance matrix.*",
)

# CLI
_parser = argparse.ArgumentParser(description=__doc__)
_parser.add_argument(
    "--autism-only", action="store_true",
    help="Cluster only on PopulationS1 == 'Autism' subjects. Non-autistic "
         "subjects get NaN cluster labels. Outputs go to outputs/curated_autism/ "
         "and sibling files become df_clusters_complete_curated_autism{,_ward,_gmm}.csv.",
)
_args = _parser.parse_args()
AUTISM_ONLY = _args.autism_only
# LEAP wave backfill is now the default: for any LEAP participant lacking IQ
# or SRS-2 at T1, the script falls back to T2 then T3 (earliest-wins). This
# adds ~12 LEAP participants vs T1-only loading and the new sample is more
# representative of the full LEAP cohort.
_RUN_TAG = "curated_autism" if AUTISM_ONLY else "curated"
print(f"=== Mode: {'AUTISM-ONLY' if AUTISM_ONLY else 'ALL POPULATIONS'} "
      f"(LEAP T1→T2→T3 backfill: ON) ===")


# -----------------------------------------------------------------------------
# Clustering helpers (identical to 04_run_clustering.py)
# -----------------------------------------------------------------------------
def _label_by_size_desc(raw_labels, k):
    """Relabel cluster IDs by descending cluster size: largest → 0 (= C1)."""
    sizes = pd.Series(raw_labels).value_counts()
    order = list(sizes.index[:k])
    remap = {old: new for new, old in enumerate(order)}
    return np.array([remap[l] for l in raw_labels])


def _kmeans_assignments(X_scaled, k, random_state=42):
    km = KMeans(n_clusters=k, n_init=25, random_state=random_state).fit(X_scaled)
    return _label_by_size_desc(km.labels_, k)


def _ward_assignments(X_scaled, k):
    """omniplot.plot.heatmap convention: Ward on the n×n distance matrix."""
    D = squareform(pdist(X_scaled, metric="euclidean"))
    Z = linkage(D, method="ward")
    raw = fcluster(Z, t=k, criterion="maxclust") - 1
    return _label_by_size_desc(raw, k)


def _gmm_assignments(X_scaled, k, random_state=42):
    gmm = GaussianMixture(
        n_components=k, covariance_type="full", n_init=10, random_state=random_state,
    ).fit(X_scaled)
    return _label_by_size_desc(gmm.predict(X_scaled), k)


# -----------------------------------------------------------------------------
# Curated clinical data loading
# -----------------------------------------------------------------------------
DATA_BASE_PATH = "/Volumes/Imaging5/EEG_MRI-MF"

CURATED_FILES = {
    # cohort_label : (path_relative_to_DATA_BASE_PATH, source_tag)
    "LEAP":    ("LEAP/_clinical_data/curated/LEAP_clinical_curated_t1.tsv",  "LEAP"),
    "INOVAND": ("INOVAND/_clinical_data/curated/INOVAND_clinical_curated_merged.tsv", "INOVAND"),
    "INFOR":   ("INFOR/_clinical_data/curated/INFOR_clinical_curated.tsv",   "INOVAND"),
    # INFOR is collapsed into INOVAND per the eeg_mri-pipeline convention.
}


def _load_curated(path, source_tag):
    """Load one curated clinical TSV, filter to participants, compute IQ."""
    if not os.path.exists(path):
        print(f"  WARNING: {path} not found — skipping.")
        return pd.DataFrame()
    df = pd.read_csv(path, sep="\t", low_memory=False)
    df = df.drop_duplicates(subset=["ID"])
    if "relation_to_proposant" in df.columns:
        df = df[df["relation_to_proposant"] == "participant"].copy()
    # Full-scale IQ with performance_IQ fallback (matches the paper pipeline)
    df["IQ"] = df["total_IQ"].fillna(df.get("performance_IQ", np.nan))
    # Carry forward the population/diagnosis columns under their canonical names.
    if "population" in df.columns and "Population1" not in df.columns:
        df["Population1"] = df["population"]
    if "Population1" in df.columns and "PopulationS1" not in df.columns:
        df["PopulationS1"] = df["Population1"]
    # Normalise PopulationS1 to the convention expected by downstream cluster-
    # aware scripts: {"Autism", "NT", "IDD"}.  The curated TSVs split autism
    # into "Autism with IDD" / "Autism without IDD" — collapse them.
    if "PopulationS1" in df.columns:
        df["PopulationS1_raw"] = df["PopulationS1"]
        df["PopulationS1"] = (
            df["PopulationS1"].astype(str)
              .str.replace("Autism with IDD", "Autism", regex=False)
              .str.replace("Autism without IDD", "Autism", regex=False)
              .str.replace("TD", "NT", regex=False)
              .str.replace(r"^ID$", "IDD", regex=True)
        )

    # Normalise MRI_ID for INOVAND/INFOR (curated TSV stores it as a float
    # string like "894.0"; v2 join expects an int-parseable token, optionally
    # prefixed with "sub-").  Convert numeric MRI_IDs to clean ints; leave
    # "sub-XXXX" strings (INFOR) untouched.
    if "MRI_ID" in df.columns:
        def _norm_mri(x):
            s = str(x).strip()
            if s.lower() in ("", "nan"):
                return None
            if s.startswith("sub-"):
                return s
            try:
                return str(int(float(s)))
            except (TypeError, ValueError):
                return s
        df["MRI_ID"] = df["MRI_ID"].apply(_norm_mri)
    # Two cohort columns: keep the original (e.g. "LEAP_T1") under
    # `cohort_raw` and overwrite the canonical `cohort` to match the convention
    # used by downstream cluster-aware scripts (`LEAP`, `INOVAND`).
    if "cohort" in df.columns:
        df["cohort_raw"] = df["cohort"]
    df["cohort"] = source_tag
    df["Cohort"] = source_tag
    df["ID"] = df["ID"].astype(str)
    print(f"  Loaded {len(df):>5} participant rows from {os.path.basename(path)}  → tagged '{source_tag}'")
    return df


def _load_leap_with_wave_backfill():
    """Load LEAP T1/T2/T3; for each ID keep the EARLIEST wave with complete
    IQ + SRS_tscore. Tagged with cohort = 'LEAP' (wave-of-origin in
    `leap_wave_used` column). Returns a DataFrame with one row per ID."""
    waves = [("t1", "LEAP/_clinical_data/curated/LEAP_clinical_curated_t1.tsv"),
             ("t2", "LEAP/_clinical_data/curated/LEAP_clinical_curated_t2.tsv"),
             ("t3", "LEAP/_clinical_data/curated/LEAP_clinical_curated_t3.tsv")]
    frames = []
    print("[--leap-backfill-waves] Loading T1, T2, T3 with priority backfill...")
    for wave_id, rel in waves:
        p = os.path.join(DATA_BASE_PATH, rel)
        df = _load_curated(p, source_tag="LEAP")
        if df.empty:
            continue
        # Keep only rows where IQ AND SRS-2 are both available
        has_iq = pd.to_numeric(df["IQ"], errors="coerce").notna()
        has_srs = pd.to_numeric(df.get("SRS_tscore", np.nan), errors="coerce").notna()
        df = df[has_iq & has_srs].copy()
        df["_wave_order"] = {"t1": 0, "t2": 1, "t3": 2}[wave_id]
        df["leap_wave_used"] = wave_id.upper()
        print(f"    {wave_id.upper()}: {len(df)} participants with complete IQ+SRS")
        frames.append(df)
    if not frames:
        return pd.DataFrame()
    all_df = pd.concat(frames, ignore_index=True)
    # Priority backfill: T1 wins, else T2, else T3.
    out = (all_df.sort_values(["ID", "_wave_order"])
                 .drop_duplicates(subset=["ID"], keep="first")
                 .drop(columns="_wave_order"))
    breakdown = out["leap_wave_used"].value_counts().to_dict()
    print(f"  → Backfilled LEAP: {len(out)} unique participants "
          f"(wave-of-origin: {breakdown})")
    return out


def load_all_cohorts():
    """Load LEAP (T1→T2→T3 backfill), INOVAND, INFOR; collapse INFOR into
    INOVAND; dedup by ID."""
    frames = []
    for label, (rel_path, tag) in CURATED_FILES.items():
        if label == "LEAP":
            df = _load_leap_with_wave_backfill()
        else:
            full_path = os.path.join(DATA_BASE_PATH, rel_path)
            df = _load_curated(full_path, source_tag=tag)
        if not df.empty:
            frames.append(df)
    if not frames:
        raise FileNotFoundError(
            f"No curated clinical TSVs found under {DATA_BASE_PATH}. "
            "Ensure the Imaging5 volume is mounted."
        )
    combined = pd.concat(frames, ignore_index=True).drop_duplicates(subset=["ID"])
    print(f"\n  Combined unique IDs: {len(combined):,}   (Cohort tags: "
          f"{combined['Cohort'].value_counts().to_dict()})")
    return combined


# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
_script_dir = os.path.dirname(os.path.abspath(__file__))
OUTPUT_BASE = os.path.normpath(os.path.join(_script_dir, "..", "outputs", _RUN_TAG))
FIGURES_DIR = os.path.join(OUTPUT_BASE, "figures")
TABLES_DIR  = os.path.join(OUTPUT_BASE, "tables")
os.makedirs(FIGURES_DIR, exist_ok=True)
os.makedirs(TABLES_DIR, exist_ok=True)

# Sibling-repo propagation (parallel file alongside df_clusters_complete_kmeans.csv).
SIBLING_CLUSTERS_DIR  = os.path.normpath(os.path.join(
    _script_dir, "..", "..", "..", "eeg_mri-pipeline",
    "results", "dataset_paper", "dataframes",
))
SIBLING_REF_FILE         = os.path.join(SIBLING_CLUSTERS_DIR, "df_clusters_complete.csv")
_SUFFIX = "_curated_autism" if AUTISM_ONLY else "_curated"
SIBLING_CURATED_FILE     = os.path.join(SIBLING_CLUSTERS_DIR, f"df_clusters_complete{_SUFFIX}.csv")
SIBLING_CURATED_WARD     = os.path.join(SIBLING_CLUSTERS_DIR, f"df_clusters_complete{_SUFFIX}_ward.csv")
SIBLING_CURATED_GMM      = os.path.join(SIBLING_CLUSTERS_DIR, f"df_clusters_complete{_SUFFIX}_gmm.csv")

N_CLUSTERS   = 3
RANDOM_STATE = 42


# -----------------------------------------------------------------------------
# Data preparation
# -----------------------------------------------------------------------------
print("=== Loading curated clinical data (LEAP + INOVAND + INFOR) ===")
df = load_all_cohorts()

clinical_features = ["IQ", "SRS_tscore"]
# Make sure both clinical features exist (and are numeric)
for c in clinical_features:
    if c not in df.columns:
        df[c] = np.nan
    df[c] = pd.to_numeric(df[c], errors="coerce")

# Carry forward any diagnosis columns we have available downstream
meta_cols = [c for c in ["ID", "Cohort", "Population1", "PopulationS1",
                          "Sex", "age_yrs", "pheno_ASD", "pheno_ID", "pheno_ADHD"]
             if c in df.columns]
df_clust = df[clinical_features + meta_cols].dropna(subset=clinical_features).copy()

if AUTISM_ONLY:
    before = len(df_clust)
    if "PopulationS1" not in df_clust.columns:
        raise RuntimeError(
            "--autism-only requested but PopulationS1 column is missing from the loaded data."
        )
    df_clust = df_clust[df_clust["PopulationS1"] == "Autism"].copy()
    print(f"\n=== --autism-only: filtered to PopulationS1 == 'Autism' "
          f"({before} → {len(df_clust)} participants) ===")

X = df_clust[clinical_features].values
X_scaled = StandardScaler().fit_transform(X)
print(f"\n=== Clustering on {len(df_clust)} participants (IQ × SRS-2 complete"
      f"{', Autism only' if AUTISM_ONLY else ''}) ===")
print(f"  Cohort breakdown of clustering set: "
      f"{df_clust['Cohort'].value_counts().to_dict() if 'Cohort' in df_clust else '<no Cohort col>'}")


# -----------------------------------------------------------------------------
# Clustering — primary k-means + sensitivity Ward + sensitivity GMM
# -----------------------------------------------------------------------------
labels_kmeans = _kmeans_assignments(X_scaled, N_CLUSTERS, random_state=RANDOM_STATE)
labels_ward   = _ward_assignments(X_scaled, N_CLUSTERS)
labels_gmm    = _gmm_assignments(X_scaled, N_CLUSTERS, random_state=RANDOM_STATE)

df_clust["Cluster"]      = [f"C{c + 1}" for c in labels_kmeans]
df_clust["Cluster_ward"] = [f"C{c + 1}" for c in labels_ward]
df_clust["Cluster_gmm"]  = [f"C{c + 1}" for c in labels_gmm]


def _emit(method_name, label_col, scatter_pdf, assign_csv, full_csv):
    """Write scatter + per-method tables for one labelling method."""
    # Same palette as 2_genetic_analysis/scripts/_config.py::PALETTE_CLUSTERS
    palette = {"C1": "#7A8B47", "C2": "#ff9fa0", "C3": "#e7ba52"}
    fig, ax = plt.subplots(figsize=(10, 8))
    sns.scatterplot(
        data=df_clust, x="SRS_tscore", y="IQ",
        hue=label_col, palette=palette, alpha=0.7, s=50, ax=ax,
    )
    ax.set_xlabel("SRS-2 t-score")
    ax.set_ylabel("Full-scale IQ")
    ax.set_title(f"Clinical Clustering — curated data (k=3, {method_name})")
    ax.legend(title="Cluster")
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, scatter_pdf), dpi=300, bbox_inches="tight")
    plt.close()

    keep_cols = [c for c in ["ID", "IQ", "SRS_tscore", label_col, "Cohort",
                              "Population1", "PopulationS1"] if c in df_clust.columns]
    df_clust[keep_cols].rename(columns={label_col: "Cluster"}).to_csv(
        os.path.join(TABLES_DIR, assign_csv), index=False,
    )
    full = df.merge(df_clust[["ID", label_col]].rename(columns={label_col: "Cluster"}),
                    on="ID", how="left")
    full.to_csv(os.path.join(TABLES_DIR, full_csv), index=False)

    print(f"  [{method_name}] sizes:", df_clust[label_col].value_counts().sort_index().to_dict())
    print(f"  [{method_name}] saved: {scatter_pdf}, {assign_csv}, {full_csv}")


print("\n--- Writing per-method outputs ---")
_fp = _SUFFIX[1:]  # "curated" or "curated_autism"
_emit("k-means", "Cluster",
      f"clustering_scatter_{_fp}.pdf",
      f"cluster_assignments_{_fp}.csv",
      f"individuals_metrics_with_clusters_{_fp}.csv")
_emit("Ward (sensitivity)", "Cluster_ward",
      f"clustering_scatter_{_fp}_ward.pdf",
      f"cluster_assignments_{_fp}_ward.csv",
      f"individuals_metrics_with_clusters_{_fp}_ward.csv")
_emit("GMM (sensitivity)", "Cluster_gmm",
      f"clustering_scatter_{_fp}_gmm.pdf",
      f"cluster_assignments_{_fp}_gmm.csv",
      f"individuals_metrics_with_clusters_{_fp}_gmm.csv")


# -----------------------------------------------------------------------------
# Cohort × cluster cross-tab (sanity check / table for the appendix)
# -----------------------------------------------------------------------------
print("\n--- Cluster × cohort cross-tab (k-means primary) ---")
if "Cohort" in df_clust.columns:
    ct = pd.crosstab(df_clust["Cluster"], df_clust["Cohort"], margins=True)
    print(ct)
    ct.to_csv(os.path.join(TABLES_DIR, f"cluster_by_cohort_{_fp}.csv"))

# Sibling-repo propagation: write df_clusters_complete_curated{,_ward,_gmm}.csv
# Each file has the same row structure as df_clusters_complete.csv but with
# the Cluster column replaced by the labels of the corresponding method.
#
# We also patch MRI_ID for INFOR-source subjects, because the paper's
# df_clusters_complete.csv encodes INFOR with "sub-2XXX" MRI_IDs while the
# freesurfer ALL/.../qc12 file uses "1, 2, … 35".  The clinical
# INFOR_clinical_curated.tsv uses "sub-000X" which strips to the same int —
# so we override MRI_ID with the clinical-TSV value.  Without this, INFOR
# subjects get silently dropped from the v2 anatomical merge.
print("\n--- Writing sibling-repo df_clusters_complete_curated{,_ward,_gmm}.csv ---")
if os.path.exists(SIBLING_REF_FILE):
    ref_df = pd.read_csv(SIBLING_REF_FILE, low_memory=False)
    ref_df["ID"] = ref_df["ID"].astype(str)
    ref_df = ref_df.drop(columns=["Cluster"], errors="ignore")

    # Build an INFOR ID→MRI_ID map from the clinical TSV (the freesurfer-
    # compatible encoding: "sub-0001" … "sub-0035").
    infor_path = os.path.join(DATA_BASE_PATH, CURATED_FILES["INFOR"][0])
    if os.path.exists(infor_path):
        infor_df = pd.read_csv(infor_path, sep="\t", low_memory=False)
        infor_df["ID"] = infor_df["ID"].astype(str)
        infor_map = dict(zip(infor_df["ID"], infor_df["MRI_ID"].astype(str)))
        # Override MRI_ID for any INFOR-source row in the sibling table
        infor_ids = set(infor_map)
        mask = ref_df["ID"].isin(infor_ids)
        if mask.any():
            n_before = ref_df.loc[mask, "MRI_ID"].notna().sum()
            ref_df.loc[mask, "MRI_ID"] = ref_df.loc[mask, "ID"].map(infor_map)
            n_after = ref_df.loc[mask, "MRI_ID"].notna().sum()
            print(f"  INFOR MRI_ID alignment: patched {mask.sum()} rows "
                  f"(non-null MRI_ID: {n_before} → {n_after})")

    # Some labelled rows in the paper sibling file have a missing `cohort`
    # tag (~91 rows) — set it to "INOVAND" for any ID that the curated
    # INOVAND or INFOR clinical TSV recognises so the v2 join key logic
    # picks them up.
    inovand_path = os.path.join(DATA_BASE_PATH, CURATED_FILES["INOVAND"][0])
    known_inovand_ids = set()
    if os.path.exists(inovand_path):
        known_inovand_ids |= set(pd.read_csv(inovand_path, sep="\t",
                                              usecols=["ID"], low_memory=False)["ID"].astype(str))
    if os.path.exists(infor_path):
        known_inovand_ids |= set(infor_df["ID"])
    if "cohort" in ref_df.columns:
        nan_mask = ref_df["cohort"].isna() & ref_df["ID"].isin(known_inovand_ids)
        n_fix = int(nan_mask.sum())
        if n_fix:
            ref_df.loc[nan_mask, "cohort"] = "INOVAND"
            print(f"  cohort=NaN rows tagged as 'INOVAND' (known IDs): {n_fix}")

    # Outer-merge so the ~25 newly clustered IDs from the merged INOVAND TSV
    # that aren't in the frozen paper sibling file are not silently dropped.
    # Backfill key clinical / cohort metadata for those new rows from the
    # combined curated `df` so downstream joins (on cohort + MRI_ID) resolve.
    ref_ids = set(ref_df["ID"])
    # Mapping: curated-df column -> sibling-file column.
    _BACKFILL = {
        "Cohort": "cohort", "MRI_ID": "MRI_ID", "Sex": "Sex",
        "age_yrs": "age_yrs", "pheno_ASD": "pheno_ASD",
        "pheno_ID": "pheno_ID", "pheno_ADHD": "pheno_ADHD",
        "Population1": "Population1", "PopulationS1": "PopulationS1",
        "total_IQ": "total_IQ", "verbal_IQ": "verbal_IQ",
        "performance_IQ": "performance_IQ", "IQ": "IQ",
        "SRS_tscore": "SRS_tscore", "SRS_rawscore": "SRS_rawscore",
        "FID": "FID", "relation_to_proposant": "Relation_to_proposant",
    }
    df_lookup = df.drop_duplicates(subset=["ID"]).set_index("ID")

    for label_col, out_path in [
        ("Cluster",      SIBLING_CURATED_FILE),
        ("Cluster_ward", SIBLING_CURATED_WARD),
        ("Cluster_gmm",  SIBLING_CURATED_GMM),
    ]:
        new_labels = df_clust[["ID", label_col]].rename(columns={label_col: "Cluster"}).copy()
        new_labels["ID"] = new_labels["ID"].astype(str)
        out = ref_df.merge(new_labels, on="ID", how="outer")

        new_mask = ~out["ID"].isin(ref_ids)
        n_new = int(new_mask.sum())
        if n_new:
            for src, dst in _BACKFILL.items():
                if src not in df_lookup.columns:
                    continue
                if dst not in out.columns:
                    out[dst] = np.nan
                out.loc[new_mask, dst] = out.loc[new_mask, "ID"].map(df_lookup[src])

        out.to_csv(out_path, index=False)
        n_lab = int(out["Cluster"].notna().sum())
        n_paper_unlab = int((~new_mask & out["Cluster"].isna()).sum())
        print(f"  wrote: {os.path.basename(out_path)}  "
              f"({len(out):,} rows, {n_lab:,} labelled, "
              f"{n_new} new-from-curated, {n_paper_unlab} paper-only-unlabelled)")
else:
    print(f"  WARNING: {SIBLING_REF_FILE} not found — skipping sibling-repo propagation.")

print("\nDone. Outputs in:", OUTPUT_BASE)
