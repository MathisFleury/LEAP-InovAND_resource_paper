# =============================================================================
# 4_anatomical_analysis/curated config — thin layer over the project config.py
# =============================================================================
# Shared constants (MRI table, v4 carrier file, AUTISM_ONLY / N_PERM / SEED)
# come from the root config.py. Only the env-driven, section-local bits
# (output base, cluster roster, curated-clinical toggle) stay here.
#
# Used by: 01_anatomical_mri_autism_nt_v2.py, 05_loeuf_mri_correlations_hg38_v2.py,
#          08_loeuf_mri_regression_permutation_v2.py.
# =============================================================================

import os
import sys

_script_dir  = os.path.dirname(os.path.abspath(__file__))
_section_dir = os.path.normpath(os.path.join(_script_dir, ".."))              # .../curated
_project_dir = os.path.normpath(os.path.join(_script_dir, "..", "..", ".."))  # LEAP-InovAND_resource
if _project_dir not in sys.path:
    sys.path.insert(0, _project_dir)

# Shared: MRI table (qc1<->qc12 switch lives in config.py), v4 carrier, toggles.
from config import MRI_FILE, CARRIER_V4, AUTISM_ONLY, N_PERM, SEED, CLUSTERS_KMEANS_FROZEN  # noqa: E402,F401

# --- Output base (env-overridable) ------------------------------------------
OUTPUT_BASE_NAME = os.environ.get("LOEUF_OUT_DIR", "outputs")
OUTPUT_BASE = os.path.join(_section_dir, OUTPUT_BASE_NAME)

# --- Phenotype / cluster rosters --------------------------------------------
# Autism-vs-NT roster; CLUSTER_FILE env overrides it (curated runs).
DF_CLUSTERS_FILE = os.environ.get("CLUSTER_FILE", CLUSTERS_KMEANS_FROZEN)
# Genetics<->MRI bridge: always the paper genetic roster (barcodes), independent
# of CLUSTER_FILE (see 05 notes).
GENETIC_ROSTER = CLUSTERS_KMEANS_FROZEN

# When True, 01 loads curated per-cohort clinical TSVs instead of DF_CLUSTERS_FILE.
# Curated is the priority regime (see CLAUDE.md), so default ON; set
# CURATED_CLINICAL=0 to reproduce the frozen paper roster.
CURATED_CLINICAL = os.environ.get("CURATED_CLINICAL", "1") == "1"
