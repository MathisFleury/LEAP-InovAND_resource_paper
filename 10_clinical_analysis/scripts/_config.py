# =============================================================================
# Shared Configuration for Clinical Analysis Scripts
# =============================================================================
# Paths and constants for 10_clinical_analysis.
# Inherits shared palettes and cluster paths from 2_genetic_analysis.
# =============================================================================

import os
import importlib.util
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

_script_dir = os.path.dirname(os.path.abspath(__file__))
_section_dir = os.path.dirname(_script_dir)
_project_dir = os.path.dirname(_section_dir)
_lib_dir = os.path.dirname(_project_dir)

# --- Load shared config from 2_genetic_analysis ---
_shared_spec = importlib.util.spec_from_file_location(
    "_shared_config",
    os.path.join(_project_dir, "2_genetic_analysis", "scripts", "_config.py"),
)
_shared = importlib.util.module_from_spec(_shared_spec)
_shared_spec.loader.exec_module(_shared)

# Re-export shared constants
PALETTE_CLUSTERS = _shared.PALETTE_CLUSTERS
ORDER_CLUSTERS = _shared.ORDER_CLUSTERS
CLUSTERS_FILE = _shared.CLUSTERS_FILE

# --- Matplotlib defaults (same as project) ---
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Helvetica", "Arial", "DejaVu Sans", "Liberation Sans"]
mpl.rcParams["pdf.fonttype"] = 42
sns.set_context("paper")

# --- Output directories ---
OUTPUT_BASE = os.path.join(_section_dir, "outputs")
FIGURES_DIR = os.path.join(OUTPUT_BASE, "figures")
TABLES_DIR = os.path.join(OUTPUT_BASE, "tables")

# --- Input data files ---
FONDA_CSV = os.path.join(
    "/Volumes", "Imaging5", "EEG_MRI-MF", "INOVAND",
    "_clinical_data", "raw", "concat_fonda_FIRST_ACQUISITIONS_251112.csv",
)
VERBAL_TSV = os.path.join(
    "/Volumes", "Imaging5", "EEG_MRI-MF", "INOVAND",
    "_clinical_data", "raw", "df_verbalornot.tsv",
)
LEAP_ADIR_CSV = os.path.join(
    _lib_dir, "eeg_mri-pipeline", "results", "dataset_paper", "dataframes",
    "T1-ADIR-earlydevelconcernskills-forthomas.csv",
)

# --- LEAP → FONDA column mapping (same measures, different naming) ---
LEAP_TO_FONDA = {
    "5: First walked unaided [Other]": "AGEPAS",
    "9: Age of first single words [Other]": "AGEMOTS",
    "10: Age of first phrases (if ever used) [Other]": "AGEPHRAS",
}

# --- Psychomotor milestone columns and labels ---
MILESTONE_COLS = ["AGEASSIS", "AGEPAS", "AGEMOTS", "AGEPHRAS"]

MILESTONE_LABELS = {
    "AGEASSIS": "Age Sitting\n(months)",
    "AGEPAS": "Age Walking\n(months)",
    "AGEMOTS": "Age First Words\n(months)",
    "AGEPHRAS": "Age First Phrases\n(months)",
}

# --- Verbal status labels ---
VERBAL_LABELS = {
    "Verbal": "Verbal",
    "Speech impairment/non verbal": "Non-verbal /\nSpeech impairment",
}

PALETTE_VERBAL = {
    "Verbal": "#5CAEE1",
    "Speech impairment/non verbal": "#324095",
}
