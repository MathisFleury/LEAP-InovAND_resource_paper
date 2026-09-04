# =============================================================================
# 2_genetic_analysis config — thin re-export of the project-wide config.py
# =============================================================================
# All shared constants (paths, palettes, gene-column lists, labels) now live in
# the root config.py (single source of truth). This file only adds the
# section-relative output directories, which must stay section-local.
# =============================================================================

import os
import sys

_script_dir = os.path.dirname(os.path.abspath(__file__))
_project_dir = os.path.normpath(os.path.join(_script_dir, "..", ".."))
if _project_dir not in sys.path:
    sys.path.insert(0, _project_dir)

from config import *  # noqa: F401,F403,E402  (shared constants + mpl defaults)

# --- Section-relative output directories ------------------------------------
OUTPUT_BASE = os.path.normpath(os.path.join(_script_dir, "..", "outputs"))
# Legacy (hg19 / hg38-v2) output base defaults to the same OUTPUT_BASE, but
# legacy/2_genetic_analysis/run_all.sh overrides it to a self-contained
# location so legacy outputs don't sit inside the main tree's outputs/
# alongside the current v4 ones.
_LEGACY_OUTPUT_BASE = os.environ.get("GENETIC_LEGACY_OUTPUT_BASE", OUTPUT_BASE)
FIGURES_DIR = os.path.join(_LEGACY_OUTPUT_BASE, "figures")
TABLES_DIR = os.path.join(_LEGACY_OUTPUT_BASE, "tables")
FIGURES_DIR_HG38 = os.path.join(_LEGACY_OUTPUT_BASE, "figures_hg38")
TABLES_DIR_HG38 = os.path.join(_LEGACY_OUTPUT_BASE, "tables_hg38")
# Variable names keep the _V4 suffix (still distinct from the legacy
# FIGURES_DIR/TABLES_DIR above, which this same file also exports) even
# though gnomAD v4 is now the only/current regime -- the folder names
# themselves don't carry "v4" since there's nothing left to disambiguate
# from in the main tree (legacy lives fully self-contained under legacy/).
FIGURES_DIR_V4 = os.path.join(OUTPUT_BASE, "figures")
TABLES_DIR_V4 = os.path.join(OUTPUT_BASE, "tables")
