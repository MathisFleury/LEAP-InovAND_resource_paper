# =============================================================================
# 11_age_sex_analysis config — thin layer over the project config.py
# =============================================================================
# Shared constants (MRI table, v4 carrier file, AUTISM_ONLY / N_PERM / SEED)
# come from the root config.py. Only the env-driven, section-local bits
# (output base, curated-clinical toggle) stay here.
# =============================================================================

import os
import sys

_script_dir  = os.path.dirname(os.path.abspath(__file__))
_section_dir = os.path.normpath(os.path.join(_script_dir, ".."))
_project_dir = os.path.normpath(os.path.join(_script_dir, "..", ".."))  # LEAP-InovAND_resource
if _project_dir not in sys.path:
    sys.path.insert(0, _project_dir)

from config import MRI_FILE, CARRIER_V4, AUTISM_ONLY, N_PERM, SEED  # noqa: E402,F401

# --- Output base (env-overridable) ------------------------------------------
OUTPUT_BASE_NAME = os.environ.get("LOEUF_OUT_DIR", "outputs")
OUTPUT_BASE = os.path.join(_section_dir, OUTPUT_BASE_NAME)

# When True, 01 loads curated per-cohort clinical TSVs instead of the roster.
CURATED_CLINICAL = os.environ.get("CURATED_CLINICAL", "0") == "1"
