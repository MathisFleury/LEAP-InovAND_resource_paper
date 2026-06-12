"""
Configuration for the LEAP + INOVAND fMRI connectivity preprocessing pipeline.

Modular port of `eeg_mri-pipeline/analysis/fmri/data_preprocessing/preprocessing_cohort.py`
modelled on the structural pipeline at
`eeg_mri-pipeline/analysis/mri/zscore_autism_td_analysis/pipeline_ndd`.

    load (LEAP + INOVAND) → QC filter (mean_fd) → assign machine batch
        → load connectivity tsv per subject (Fisher z, NaN ROIs handled)
        → ComBat (neuroHarmonize, batch = machine, biological covars preserved)
        → per-feature linear regression (age, age², sex, cohort)
        → z-score → save df_conn_cohort_norm.csv

Defaults assume the user's data layout under /Volumes/Imaging5; paths can be
overridden by setting LEAP_INOVAND_FMRI_* environment variables.
"""
from __future__ import annotations

import os
from pathlib import Path


def _env_path(key: str, default: str) -> Path:
    return Path(os.environ.get(key, default))


# --------------------------------------------------------------------- paths
XCP_D_SUBJECTS_LEAP = _env_path(
    "LEAP_INOVAND_FMRI_XCP_LEAP",
    "/Users/mfleury/POSTDOC/LIBRAIRY/eeg_mri-pipeline/results/dataframes/fmri/df_xcp_subjects_LEAP.csv",
)
XCP_D_SUBJECTS_INOVAND = _env_path(
    "LEAP_INOVAND_FMRI_XCP_INOVAND",
    "/Users/mfleury/POSTDOC/LIBRAIRY/eeg_mri-pipeline/results/dataframes/fmri/df_xcp_subjects_INOVAND.csv",
)
QC_XCP_LEAP = _env_path(
    "LEAP_INOVAND_FMRI_QC_LEAP",
    "/Users/mfleury/POSTDOC/LIBRAIRY/eeg_mri-pipeline/results/func_qc/dataframes/qc_xcp_d_LEAP.csv",
)
QC_XCP_INOVAND = _env_path(
    "LEAP_INOVAND_FMRI_QC_INOVAND",
    "/Users/mfleury/POSTDOC/LIBRAIRY/eeg_mri-pipeline/results/func_qc/dataframes/qc_xcp_d_INOVAND.csv",
)

OUT_DIR = _env_path(
    "LEAP_INOVAND_FMRI_OUT",
    str(Path(__file__).resolve().parent / "outputs"),
)
OUT_DIR.mkdir(parents=True, exist_ok=True)

# --------------------------------------------------------------- atlas / parcellation
ATLAS = "4S156Parcels"
# The connectivity .tsv files referenced by `path` in df_xcp_subjects_* are
# already per-subject correlation matrices in this parcellation. The atlas
# TSV is only used by downstream scripts (network labels).

# --------------------------------------------------------------- cohort filter
# Diagnostic codes treated as autism / ASD-related.
ASD_STATUSES = (
    "autistic_index",
    "autistic_relatives",
    "other_NDD",
    "relatives_famASD",
)
# Plus rows where `control_status` is NaN — preserved in the original pipeline
# as ASD-bucket. Keep flag here so it's discoverable.
KEEP_NAN_STATUS_AS_ASD = True

CONTROL_STATUSES = (
    "control_index",
    "control_relatives",
    "control_relatives_famASD",
)

# --------------------------------------------------------------- QC policy
# Drop scans with mean framewise displacement above this threshold (in mm).
# 1.0 mm matches the user's current operating choice; lower the value for a
# stricter QC run (e.g. 0.3 to recover the FD<0.3 sensitivity branch).
FD_THRESHOLD_MM = 0.5
QC_TAG = f"fd{int(FD_THRESHOLD_MM * 10):02d}"  # e.g. fd10 / fd03 / fd05

# When a subject has multiple sessions/runs surviving QC, keep one row.
# `lowest_fd` keeps the run with the smallest mean_fd; `first` keeps the
# first row (alphabetical session/run). The original pipeline used a
# `drop_duplicates(subset=['ID'])` which is non-deterministic — we make it
# explicit.
DEDUP_STRATEGY = "lowest_fd"

# --------------------------------------------------------------- batch column
# `machine` is populated for INOVAND directly; for LEAP we fall back on
# (LEAP_t1_site + 5) — the same offset trick as preprocessing_cohort.py,
# which guarantees LEAP/INOVAND batches do not collide.
BATCH_COL = "machine_batch"
LEAP_SITE_OFFSET = 5

# --------------------------------------------------------------- covariates
AGE_COL = "age_yrs"
AGE2_COL = "age_yrs_sq"
SEX_COL = "Sex"
COHORT_COL = "cohort"
ID_COL = "ID"

# Biological + nuisance covariates preserved by ComBat (so the batch model
# does not absorb their signal).
COMBAT_COVARS = [AGE_COL, AGE2_COL, SEX_COL, COHORT_COL]

# Same covariates are regressed out of the harmonized features in step 05.
REGRESS_COVARS = [AGE_COL, AGE2_COL, SEX_COL, COHORT_COL]

# --------------------------------------------------------------- normalisation
# Apply median-centred / IQR-scaled within-subject rescale before Fisher z.
# Replicates the original preprocessing_cohort.py behaviour. Disable to
# study raw correlation distributions.
APPLY_IQR_RESCALE = True
APPLY_FISHER_Z = True

# Final per-feature z-score across subjects after ComBat + regression.
APPLY_ZSCORE = True

# --------------------------------------------------------------- output names
OUT_NORM_CSV = OUT_DIR / "df_conn_cohort_norm.csv"
OUT_BAD_SUBJECTS = OUT_DIR / "bad_subjects.txt"
OUT_NAN_MASK = OUT_DIR / "nan_mask_combined.npy"
OUT_AUDIT = OUT_DIR / f"n_at_each_stage_{QC_TAG}.tsv"
OUT_FEATURES = OUT_DIR / f"features_used_{QC_TAG}.txt"
OUT_COMBAT_MODEL = OUT_DIR / f"combat_model_{QC_TAG}.pkl"
OUT_REGRESS_BETAS = OUT_DIR / f"regression_betas_{QC_TAG}.npy"
