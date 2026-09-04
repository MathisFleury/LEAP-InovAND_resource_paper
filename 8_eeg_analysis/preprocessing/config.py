"""
EEG preprocessing / correction config.

Retake of eeg_mri-pipeline/.../correct_alpha_peak_LEAP_INOVAND.py, but curated
throughout: clinical demographics + clusters come ONLY from the curated k-means
cluster file, and alpha peak is rebuilt from the raw $IMG5 feature files.

Correction policy = regress age/age²/sex out of each feature, then z-score,
BOTH fit on the FULL sample (DX_COL=None) — matches the anat/func modalities.
To revert to an NT-controls reference, set DX_COL="PopulationS1", CONTROL_VALUE="NT".
"""
import importlib.util
import os
from pathlib import Path

# Root config.py = single source of truth. Loaded by path (not `import config`)
# because THIS file is also named config.py — a bare import would self-import.
_root = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "config.py"))
_spec = importlib.util.spec_from_file_location("_project_config", _root)
_pc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_pc)

# --- raw EEG feature dirs ($IMG5, must be mounted) ---
IMG5 = Path("/Volumes/Imaging5/EEG_MRI-MF")
LEAP_FEATURES = IMG5 / "LEAP/_converted_BIDS/derivatives/eeg_features/features"
INOVAND_FEATURES = IMG5 / "INOVAND/EEG_TSA/_converted_BIDS/derivatives/eeg_features/features"

# --- curated clinical: demographics come from the BROAD per-cohort clinical TSVs
#     (age_yrs, Sex, population_group->PopulationS1, population->Population1);
#     the curated k-means cluster table adds only the Cluster label. Using the
#     clinical TSVs (not the clustering table) avoids dropping every EEG subject
#     that lacked IQ+SRS to be clustered. ---
LEAP_CLINICAL = Path(_pc.LEAP_CLINICAL_CURATED)
INOVAND_CLINICAL = Path(_pc.INOVAND_CLINICAL_CURATED)
CURATED_CLUSTERS = Path(_pc.CLUSTERS_CURATED)

# --- INOVAND sub-XXXX -> clinical ID crosswalk (ID LOOKUP ONLY; every clinical
#     value still comes from CURATED_CLUSTERS). LEAP raw IDs already equal the
#     clinical ID, so only INOVAND needs a crosswalk, and it lives only here. ---
INOVAND_BIDS_MAP = Path(
    "/Users/mfleury/POSTDOC/LIBRAIRY/eeg_mri-pipeline/results/dataset_paper/"
    "dataframes/df_clusters_complete_kmeans.csv")

# --- frozen combined power spectrum: LEAP has no raw power on $IMG5, so power is
#     NOT rebuilt from raw — this long-format corrected file stays the source. ---
POWER_COMBINED = Path(
    "/Users/mfleury/POSTDOC/LIBRAIRY/eeg_mri-pipeline/results/dataset_paper/"
    "dataframes/Power_spectrum_combined_corrected.csv")

# --- outputs: corrected features consumed by BOTH 8_ and 9_ analysis ---
OUT_DIR = Path(__file__).resolve().parent / "outputs"
ALPHA_PEAK_CORRECTED = OUT_DIR / "Alpha_peak_combined_corrected.csv"

# --- demographic columns pulled from the curated file ---
DEMO_COLS = ["ID", "Cluster", "PopulationS1", "Population1", "age_yrs", "Sex", "cohort"]

# --- correction policy (single-sourced from root config.py) ---
REGRESS_COVARS = _pc.EEG_REGRESS_COVARS      # ["age_yrs", "age_yrs_sq", "Sex"]
DX_COL = _pc.DX_COL                          # None -> full-sample regression + z-score
CONTROL_VALUE = _pc.CONTROL_VALUE
