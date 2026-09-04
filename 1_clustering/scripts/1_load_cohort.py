#!/usr/bin/env python3
# =============================================================================
# 1 - Load Cohort — CURATED (current, priority regime)
# =============================================================================
# Loads and cleans the curated per-cohort clinical TSVs (LEAP, INOVAND, INFOR)
# from the data directory pointed to by LEAP_INOVAND_CURATED_DATA, and writes
# one combined, participant-only, deduplicated cohort table. This is the
# shared input for both 2_pca_features.R (feature justification) and
# 4_run_clustering.py (the clustering itself) -- kept as its own step so PCA
# can run before clustering without duplicating the loading logic.
# INFOR is collapsed into the INOVAND cohort. LEAP uses T1->T2->T3 backfill:
# for any participant lacking IQ or SRS-2 at T1, falls back to T2 then T3
# (earliest-wins).
# =============================================================================

import os

import numpy as np
import pandas as pd

# -----------------------------------------------------------------------------
# Curated clinical data loading
# -----------------------------------------------------------------------------
DATA_BASE_PATH = os.environ.get("LEAP_INOVAND_CURATED_DATA", "/Volumes/Imaging5/EEG_MRI-MF")

CURATED_FILES = {
    # cohort_label : (path_relative_to_DATA_BASE_PATH, source_tag)
    "LEAP":    ("LEAP/_clinical_data/curated/LEAP_clinical_curated_t1.tsv",  "LEAP"),
    "INOVAND": ("INOVAND/_clinical_data/curated/INOVAND_clinical_curated.tsv", "INOVAND"),
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
            "Set LEAP_INOVAND_CURATED_DATA to your data directory."
        )
    combined = pd.concat(frames, ignore_index=True).drop_duplicates(subset=["ID"])
    print(f"\n  Combined unique IDs: {len(combined):,}   (Cohort tags: "
          f"{combined['Cohort'].value_counts().to_dict()})")
    return combined


if __name__ == "__main__":
    _script_dir = os.path.dirname(os.path.abspath(__file__))
    TABLES_DIR = os.path.join(_script_dir, "..", "outputs", "tables")
    os.makedirs(TABLES_DIR, exist_ok=True)
    COHORT_OUTPUT = os.path.join(TABLES_DIR, "cohort_curated.csv")

    print("=== Loading curated clinical data (LEAP + INOVAND + INFOR) ===")
    df = load_all_cohorts()
    for c in ["IQ", "SRS_tscore"]:
        if c not in df.columns:
            df[c] = np.nan
        df[c] = pd.to_numeric(df[c], errors="coerce")

    df.to_csv(COHORT_OUTPUT, index=False)
    print(f"\nDone. Wrote {len(df):,} rows to {COHORT_OUTPUT}")
