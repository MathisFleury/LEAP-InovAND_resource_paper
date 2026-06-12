"""
Centralised MRI data loader for 4-1_anatomical_analysis.

Source: /Volumes/Imaging5/EEG_MRI-MF/ALL/results/tabular/anat/QC_demo_only/
        freesurfer_zscore_demo_only.tsv

Key differences vs 4_anatomical_analysis (COMBAT file):
  - ID column is `subject_id` (renamed to `ID` on load)
  - Subcortical columns have `_Volume_mm3` suffix (stripped on load)
  - Values are already z-scored (demo-corrected) — no further normalisation needed
  - No PopulationS1 / diagnosis column — merged from df_clusters_complete_kmeans.csv
  - No QC filtering (file already QC-passed)
"""

import os
import pandas as pd

_lib_dir = os.path.normpath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..")
)

MRI_ZSCORE_FILE = (
    "/Volumes/Imaging5/EEG_MRI-MF/ALL/results/"
    "tabular/anat/all_demo_only/freesurfer_zscore_demo_only.tsv"
)

CLUSTERS_COMPLETE_FILE = os.path.join(
    _lib_dir, "eeg_mri-pipeline", "results", "dataset_paper", "dataframes",
    "df_clusters_complete_kmeans.csv",
)


def load_mri(merge_clusters_cols=None):
    """
    Load and normalise the z-scored MRI file.

    Parameters
    ----------
    merge_clusters_cols : list[str] or None
        Extra columns to bring in from df_clusters_complete_kmeans.csv (e.g.
        ['PopulationS1', 'total_IQ', 'SRS_tscore']). 'ID' is always used as
        the merge key and need not be listed.

    Returns
    -------
    pd.DataFrame
        MRI data with:
          - `ID` column (renamed from `subject_id`)
          - subcortical columns without the `_Volume_mm3` suffix
          - any requested columns from df_clusters_complete_kmeans.csv
    """
    df = pd.read_table(MRI_ZSCORE_FILE, low_memory=False)

    # Rename ID column
    df = df.rename(columns={"subject_id": "ID"})
    df["ID"] = df["ID"].astype(str)

    # Strip _Volume_mm3 suffix from subcortical columns
    df.columns = [
        c.replace("_Volume_mm3", "") if c.endswith("_Volume_mm3") else c
        for c in df.columns
    ]

    # Standardise cohort column name
    df = df.rename(columns={"Cohort": "cohort"})

    if merge_clusters_cols:
        df_c = pd.read_csv(CLUSTERS_COMPLETE_FILE, low_memory=False)
        df_c["ID"] = df_c["ID"].astype(str)
        # Also try MRI_ID_clean as fallback merge key
        keep = list({*merge_clusters_cols, "ID"})
        df_c = df_c[[c for c in keep if c in df_c.columns]]
        df = df.merge(df_c, on="ID", how="left")

    return df
