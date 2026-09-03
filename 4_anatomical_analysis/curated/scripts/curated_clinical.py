#!/usr/bin/env python3
"""
Curated clinical loader for the curated-pipeline anatomical scripts.

Mirrors the curated-clinical loading in
1_clustering/scripts/04_run_clustering_curated.py (load_all_cohorts): merges the
per-cohort curated clinical TSVs — LEAP (T1->T2->T3 backfill), INOVAND, INFOR —
collapsing INFOR into INOVAND, normalising PopulationS1 and MRI_ID, deduping by
ID. Returns one row per ID carrying ID / MRI_ID / cohort / PopulationS1 /
population_group, suitable as the phenotype source for the Autism-vs-NT script.

Difference vs. 1_clustering: the LEAP wave-backfill there keeps only rows with
complete IQ + SRS-2 (a *clustering* requirement). The MRI case-control needs the
diagnosis, not the clustering features, so by default we backfill on diagnosis
presence (require_features=False) to avoid dropping scanned subjects who lack
IQ/SRS. Pass require_features=True to reproduce the clustering sample exactly.
"""

import os
import numpy as np
import pandas as pd

DATA_BASE_PATH = "/Volumes/Imaging5/EEG_MRI-MF"

CURATED_FILES = {
    # cohort_label : (path_relative_to_DATA_BASE_PATH, source_tag)
    "LEAP":    ("LEAP/_clinical_data/curated/LEAP_clinical_curated_t1.tsv",   "LEAP"),
    "INOVAND": ("INOVAND/_clinical_data/curated/INOVAND_clinical_curated.tsv", "INOVAND"),
    "INFOR":   ("INFOR/_clinical_data/curated/INFOR_clinical_curated.tsv",    "INOVAND"),
}


def _load_curated(path, source_tag):
    """Load one curated clinical TSV, filter to participants, normalise labels."""
    if not os.path.exists(path):
        print(f"  WARNING: {path} not found — skipping.")
        return pd.DataFrame()
    df = pd.read_csv(path, sep="\t", low_memory=False)
    df = df.drop_duplicates(subset=["ID"])
    if "relation_to_proposant" in df.columns:
        df = df[df["relation_to_proposant"] == "participant"].copy()
    df["IQ"] = df["total_IQ"].fillna(df.get("performance_IQ", np.nan)) \
        if "total_IQ" in df.columns else np.nan
    if "population" in df.columns and "Population1" not in df.columns:
        df["Population1"] = df["population"]
    if "Population1" in df.columns and "PopulationS1" not in df.columns:
        df["PopulationS1"] = df["Population1"]
    # Normalise PopulationS1 to {"Autism", "NT", "IDD"} (collapse autism subtypes).
    if "PopulationS1" in df.columns:
        df["PopulationS1_raw"] = df["PopulationS1"]
        df["PopulationS1"] = (
            df["PopulationS1"].astype(str)
              .str.replace("Autism with IDD", "Autism", regex=False)
              .str.replace("Autism without IDD", "Autism", regex=False)
              .str.replace("TD", "NT", regex=False)
              .str.replace(r"^ID$", "IDD", regex=True)
        )
    # Normalise MRI_ID (curated TSV stores "894.0"; keep "sub-XXXX" as-is).
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
    if "cohort" in df.columns:
        df["cohort_raw"] = df["cohort"]
    df["cohort"] = source_tag
    df["Cohort"] = source_tag
    df["ID"] = df["ID"].astype(str)
    print(f"  Loaded {len(df):>5} participant rows from {os.path.basename(path)}"
          f"  -> tagged '{source_tag}'")
    return df


def _load_leap_with_wave_backfill(require_features: bool = False):
    """LEAP T1/T2/T3; keep the earliest wave that satisfies the keep condition.

    require_features=True  -> earliest wave with complete IQ + SRS-2 (clustering).
    require_features=False -> earliest wave with a non-null diagnosis (MRI case
                              control: keep scanned subjects regardless of IQ/SRS).
    """
    waves = [("t1", "LEAP/_clinical_data/curated/LEAP_clinical_curated_t1.tsv"),
             ("t2", "LEAP/_clinical_data/curated/LEAP_clinical_curated_t2.tsv"),
             ("t3", "LEAP/_clinical_data/curated/LEAP_clinical_curated_t3.tsv")]
    frames = []
    for wave_id, rel in waves:
        df = _load_curated(os.path.join(DATA_BASE_PATH, rel), source_tag="LEAP")
        if df.empty:
            continue
        if require_features:
            ok = (pd.to_numeric(df["IQ"], errors="coerce").notna()
                  & pd.to_numeric(df.get("SRS_tscore", np.nan), errors="coerce").notna())
        else:
            ok = df["PopulationS1"].notna() if "PopulationS1" in df.columns else True
        df = df[ok].copy()
        df["_wave_order"] = {"t1": 0, "t2": 1, "t3": 2}[wave_id]
        df["leap_wave_used"] = wave_id.upper()
        frames.append(df)
    if not frames:
        return pd.DataFrame()
    out = (pd.concat(frames, ignore_index=True)
             .sort_values(["ID", "_wave_order"])
             .drop_duplicates(subset=["ID"], keep="first")
             .drop(columns="_wave_order"))
    print(f"  -> Backfilled LEAP: {len(out)} unique participants "
          f"(waves: {out['leap_wave_used'].value_counts().to_dict()})")
    return out


def load_all_cohorts(require_features: bool = False) -> pd.DataFrame:
    """Merge LEAP (wave backfill), INOVAND, INFOR (collapsed into INOVAND)."""
    frames = []
    for label, (rel_path, tag) in CURATED_FILES.items():
        if label == "LEAP":
            df = _load_leap_with_wave_backfill(require_features=require_features)
        else:
            df = _load_curated(os.path.join(DATA_BASE_PATH, rel_path), source_tag=tag)
        if not df.empty:
            frames.append(df)
    if not frames:
        raise FileNotFoundError(
            f"No curated clinical TSVs found under {DATA_BASE_PATH}. "
            "Ensure the Imaging5 volume is mounted."
        )
    combined = pd.concat(frames, ignore_index=True).drop_duplicates(subset=["ID"])
    print(f"  Combined curated-clinical unique IDs: {len(combined):,}  "
          f"(cohorts: {combined['cohort'].value_counts().to_dict()})")
    return combined


if __name__ == "__main__":
    d = load_all_cohorts()
    print(d[["ID", "MRI_ID", "cohort", "PopulationS1"]].head())
