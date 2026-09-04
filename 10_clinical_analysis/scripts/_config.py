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
INDIVIDUALS_METRICS = _shared.INDIVIDUALS_METRICS
CARRIER_V4 = _shared.CARRIER_V4
PALETTE_POPULATION1 = _shared.PALETTE_POPULATION1

# Curated k-means cluster solution produced WITHIN this project by
# 1_clustering/scripts/4_run_clustering.py. Carries the curated
# clinical columns (IQ, SRS_tscore, vabsabcabc_standard, Population1) with the
# k-means Cluster label. (In-project output — not the sibling-repo copy.)
CLUSTERS_CURATED_FILE = os.path.join(
    _project_dir, "1_clustering", "outputs", "tables",
    "individuals_metrics_with_clusters_curated.csv",
)

# --- Clinical scores (panel C): column -> axis label ---
# Full-scale IQ uses total_IQ with performance_IQ fallback (as in 1_clustering).
CLINICAL_SCORES = {
    "IQ": "Full-scale IQ",
    "SRS_tscore": "SRS-2 total t-score",
    "ssp_total": "SSP total score",
    "vabsabcabc_standard": "VABS-II ABC score",
}
# Reference median dashed lines (whole-cohort)
COLOR_AUTISM_MEDIAN = "#8991FA"  # autistic (with or without IDD) — canonical "Autism" group color
COLOR_NT_MEDIAN = "#C1C2BC"      # neurotypical — canonical "NT" group color

COLORS = {
    "Autism": "#8991FA",
    "Relatives": "#9AD5D3", 
    "IDD": "#D8A4CB",
    "TD": "#C1C2BC",
    "Autism with IDD": "#324095",
    "Autism without IDD": "#5CAEE1",
    "NDD": "#8991FA"  # Use autism color for NDD
}

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
# Broad curated clinical rosters — NOT gated on IQ+SRS completeness, so used for
# NT/IDD group membership (the clustering table only supplies C1/C2/C3 labels).
INOVAND_CLINICAL_CURATED = os.path.join(
    "/Volumes", "Imaging5", "EEG_MRI-MF", "INOVAND",
    "_clinical_data", "curated", "INOVAND_clinical_curated_merged.tsv",
)
LEAP_CLINICAL_CURATED = os.path.join(
    "/Volumes", "Imaging5", "EEG_MRI-MF", "LEAP",
    "_clinical_data", "curated", "LEAP_clinical_curated_t1.tsv",
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
