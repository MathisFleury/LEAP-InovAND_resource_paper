# =============================================================================
# Shared Configuration for Cluster Genetic Analysis Scripts
# =============================================================================
# Paths and constants for 3_cluster_genetic_analysis.
# Inherits shared palettes, column lists, and labels from 2_genetic_analysis.
# =============================================================================

import os
import importlib.util

_script_dir = os.path.dirname(os.path.abspath(__file__))
_lib_dir = os.path.normpath(os.path.join(_script_dir, "..", "..", ".."))

# --- Load shared config from 2_genetic_analysis ---
_shared_spec = importlib.util.spec_from_file_location(
    "_shared_config",
    os.path.join(
        _lib_dir, "LEAP-InovAND_resource", "2_genetic_analysis", "scripts", "_config.py"
    ),
)
_shared = importlib.util.module_from_spec(_shared_spec)
_shared_spec.loader.exec_module(_shared)

# Re-export all shared constants (including those needed by 02_carrier_freq_or.py)
INDIVIDUALS_METRICS = _shared.INDIVIDUALS_METRICS
DIAG_FILE = _shared.DIAG_FILE
PGS_FILE = _shared.PGS_FILE
PALETTE_CLUSTERS = _shared.PALETTE_CLUSTERS
ORDER_CLUSTERS = _shared.ORDER_CLUSTERS
PALETTE_FREQ = _shared.PALETTE_FREQ
PALETTE_POPULATION1 = _shared.PALETTE_POPULATION1
ORDER_FREQ = _shared.ORDER_FREQ
ORDER_OR = _shared.ORDER_OR
ORDER_POPULATION1 = _shared.ORDER_POPULATION1
COLS_DELLOF_CONSTRAINED = _shared.COLS_DELLOF_CONSTRAINED
COLS_DELLOF_ALL = _shared.COLS_DELLOF_ALL
COLS_DELLOFMISS_CONSTRAINED = _shared.COLS_DELLOFMISS_CONSTRAINED
COLS_DELLOFMISS_ALL = _shared.COLS_DELLOFMISS_ALL
LABELS_CONSTRAINED = _shared.LABELS_CONSTRAINED
LABELS_ALL_GENES = _shared.LABELS_ALL_GENES
LABELS_DELLOFMISS_CONSTRAINED = _shared.LABELS_DELLOFMISS_CONSTRAINED
LABELS_DELLOFMISS_ALL = _shared.LABELS_DELLOFMISS_ALL
LABELS_PGS = _shared.LABELS_PGS
PGS_TRAITS = _shared.PGS_TRAITS

# FIGURES_DIR / TABLES_DIR: point to 2_genetic_analysis outputs as fallback
# (not used directly — cluster scripts always pass explicit dirs)
FIGURES_DIR = _shared.FIGURES_DIR
TABLES_DIR = _shared.TABLES_DIR

# --- Carrier annotations from 2_genetic_analysis outputs ---
CARRIER_ANNOTATIONS = os.path.join(
    _lib_dir, "LEAP-InovAND_resource", "2_genetic_analysis", "outputs", "tables",
    "carrier_annotations.tsv",
)
CARRIER_ANNOTATIONS_HG38 = os.path.join(
    _lib_dir, "LEAP-InovAND_resource", "2_genetic_analysis", "outputs", "tables_hg38",
    "carrier_annotations_hg38.tsv",
)

# --- Cluster assignment files ---
# Legacy: from 1_clustering pipeline
CLUSTERS_LEGACY_FILE = os.path.join(
    _lib_dir, "LEAP-InovAND_resource", "1_clustering", "outputs", "tables",
    "cluster_assignments.csv",
)

# New GMM method
CLUSTERS_GMM_FILE = os.path.join(
    _lib_dir, "eeg_mri-pipeline", "results", "clinical",
    "large_clustering_datasets", "results", "df_multi_dataset_with_clusters.csv",
)

# --- Output directories ---
OUTPUT_BASE = os.path.normpath(os.path.join(_script_dir, "..", "outputs"))

FIGURES_LEGACY_DIR = os.path.join(OUTPUT_BASE, "figures_legacy")
TABLES_LEGACY_DIR = os.path.join(OUTPUT_BASE, "tables_legacy")
FIGURES_GMM_DIR = os.path.join(OUTPUT_BASE, "figures_gmm")
TABLES_GMM_DIR = os.path.join(OUTPUT_BASE, "tables_gmm")

# NT-restricted-to-C1 variants
FIGURES_LEGACY_NTC1_DIR = os.path.join(OUTPUT_BASE, "figures_legacy_nt_c1")
TABLES_LEGACY_NTC1_DIR = os.path.join(OUTPUT_BASE, "tables_legacy_nt_c1")
FIGURES_GMM_NTC1_DIR = os.path.join(OUTPUT_BASE, "figures_gmm_nt_c1")
TABLES_GMM_NTC1_DIR = os.path.join(OUTPUT_BASE, "tables_gmm_nt_c1")

# hg38 variants
FIGURES_LEGACY_HG38_DIR = os.path.join(OUTPUT_BASE, "figures_legacy_hg38")
TABLES_LEGACY_HG38_DIR = os.path.join(OUTPUT_BASE, "tables_legacy_hg38")
FIGURES_GMM_HG38_DIR = os.path.join(OUTPUT_BASE, "figures_gmm_hg38")
TABLES_GMM_HG38_DIR = os.path.join(OUTPUT_BASE, "tables_gmm_hg38")
FIGURES_LEGACY_HG38_NTC1_DIR = os.path.join(OUTPUT_BASE, "figures_legacy_hg38_nt_c1")
TABLES_LEGACY_HG38_NTC1_DIR = os.path.join(OUTPUT_BASE, "tables_legacy_hg38_nt_c1")
FIGURES_GMM_HG38_NTC1_DIR = os.path.join(OUTPUT_BASE, "figures_gmm_hg38_nt_c1")
TABLES_GMM_HG38_NTC1_DIR = os.path.join(OUTPUT_BASE, "tables_gmm_hg38_nt_c1")

# --- SPARK carrier file ---
SPARK_LOF_CARRIER = os.path.join(
    _lib_dir, "eeg_mri-pipeline", "analysis", "genetic", "SPARK_clinics",
    "tables", "SPARK_lof_carrier_annotation.tsv",
)

SPARK_PGS_FILE = os.path.join(
    "/Volumes", "Imaging5", "EEG_MRI-MF", "SPARK", "results", "tabular", "genetics",
    "SPARK-iWES-v3.pgs.20250411.tsv",
)

# SPARK output directories
FIGURES_SPARK_DIR = os.path.join(OUTPUT_BASE, "figures_spark")
TABLES_SPARK_DIR = os.path.join(OUTPUT_BASE, "tables_spark")
FIGURES_SPARK_MANUAL_DIR = os.path.join(OUTPUT_BASE, "figures_spark_manual")
TABLES_SPARK_MANUAL_DIR = os.path.join(OUTPUT_BASE, "tables_spark_manual")

# SPARK vs LEAP-InovAND NT supplementary analysis
FIGURES_SPARK_LEAPNT_DIR = os.path.join(OUTPUT_BASE, "figures_spark_leapnt")
TABLES_SPARK_LEAPNT_DIR = os.path.join(OUTPUT_BASE, "tables_spark_leapnt")
FIGURES_SPARK_LEAPNT_MANUAL_DIR = os.path.join(OUTPUT_BASE, "figures_spark_leapnt_manual")
TABLES_SPARK_LEAPNT_MANUAL_DIR = os.path.join(OUTPUT_BASE, "tables_spark_leapnt_manual")

# Manual cluster output directories (LEAP-InovAND, hg19 + hg38)
FIGURES_MANUAL_DIR = os.path.join(OUTPUT_BASE, "figures_manual")
TABLES_MANUAL_DIR = os.path.join(OUTPUT_BASE, "tables_manual")
FIGURES_MANUAL_HG38_DIR = os.path.join(OUTPUT_BASE, "figures_manual_hg38")
TABLES_MANUAL_HG38_DIR = os.path.join(OUTPUT_BASE, "tables_manual_hg38")

# --- SPARK LoF carrier columns ---
COLS_SPARK_LOF_CONSTRAINED = [
    "lof_contraint_carrier",
    "lof_hcndddomv6_contraint_carrier",
    "lof_hcnddxlinked_contraint_carrier",
    "lof_sparksfari1_contraint_carrier",
    "lof_eagle_contraint_carrier",
    "lof_syngo_contraint_carrier",
    "lof_chromepitf_contraint_carrier",
    "lof_syngo_chromepitf_contraint_carrier",
]

COLS_SPARK_LOF_ALL = [
    "lof_hcndddomv6_carrier",
    "lof_hcnddxlinked_carrier",
    "lof_sparksfari1_carrier",
    "lof_eagle_carrier",
    "lof_syngo_carrier",
    "lof_chromepitf_carrier",
    "lof_syngo_chromepitf_carrier",
]

LABELS_SPARK_LOF_CONSTRAINED = {
    "lof_contraint_carrier": "Constrained",
    "lof_hcndddomv6_contraint_carrier": "Constrained\nHCNDD DOM",
    "lof_hcnddxlinked_contraint_carrier": "Constrained\nHCNDD X-linked",
    "lof_sparksfari1_contraint_carrier": "Constrained\nSPARK SFARI1",
    "lof_eagle_contraint_carrier": "Constrained\nEAGLE d&s",
    "lof_syngo_contraint_carrier": "Constrained\nSynGO",
    "lof_chromepitf_contraint_carrier": "Constrained\nChromEpiTF",
    "lof_syngo_chromepitf_contraint_carrier": "Constrained\nSynGO+ChromEpiTF",
}

LABELS_SPARK_LOF_ALL = {
    "lof_hcndddomv6_carrier": "HCNDD DOM",
    "lof_hcnddxlinked_carrier": "HCNDD X-linked",
    "lof_sparksfari1_carrier": "SPARK SFARI1",
    "lof_eagle_carrier": "EAGLE d&s",
    "lof_syngo_carrier": "SynGO",
    "lof_chromepitf_carrier": "ChromEpiTF",
    "lof_syngo_chromepitf_carrier": "SynGO+ChromEpiTF",
}
