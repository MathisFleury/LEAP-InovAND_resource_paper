"""
Config for the XCP-D v0.11 post-processing (6_functional_analysis/concat).

XCP-D 0.11 no longer writes connectivity matrices — only parcel mean timeseries.
So connectivity is computed here (Pearson -> Fisher z) from the concatenated
4S156 `stat-mean_timeseries.tsv`, then (downstream) ComBat + regression + z-score.

Two pipelines, selected by the GSR env var / VARIANT:
  - nogsr : XCP-D params 24P  -> idp/xcp_d-0.11.0       (exists)
  - gsr   : XCP-D params 36P  -> idp/xcp_d_gsr-0.11.0   (after the GSR nipoppy run)

Reads the nipoppy IDP outputs on the mounted volume; override roots via env.
"""
from __future__ import annotations

import os
from pathlib import Path

# ---- which denoising variant to post-process ----
# nogsr = 24P/8mm (original) | gsr = 36P/6mm | nogsr6 = 24P/6mm
VARIANT = os.environ.get("XCPD_VARIANT", "nogsr")           # nogsr | gsr | nogsr6
VARIANT_DIR = {
    "nogsr":  "xcp_d-0.11.0",        # 24P, 8mm
    "gsr":    "xcp_d_gsr-0.11.0",    # 36P, 6mm
    "nogsr6": "xcp_d_sm6-0.11.0",    # 24P, 6mm
}[VARIANT]

# ---- nipoppy IDP roots (per cohort) ----
# INFOR is collapsed into INOVAND clinically but is a separate acquisition site
# (Créteil / Henri-Mondor); it's globbed as its own source so its scanners form
# distinct ComBat batches. Missing dirs are skipped by the builder, so INFOR is
# harmless until its XCP-D v0.11 IDP is populated.
_IMG5 = Path(os.environ.get("IMG5_ROOT", "/Volumes/Imaging5/EEG_MRI-MF"))
COHORT_IDP = {
    c: _IMG5 / c / "_nipoppy" / "derivatives" / "fmriprep" / "23.0.0" / "idp" / VARIANT_DIR
    for c in ("INOVAND", "LEAP", "INFOR")
}

ATLAS = "4S156Parcels"
ATLAS_FILE = Path(
    "/Users/mfleury/POSTDOC/LIBRAIRY/eeg_mri-pipeline/ressources/atlases/"
    "SCHAEFER/atlas-4S156Parcels/atlas-4S156Parcels_dseg.tsv"
)
ATLAS_NAME = "4S156"   # value in the dseg 'atlas_name' column

# timeseries filename pattern — matches BOTH per-run (..._run-01_space-...) and
# the no-run combined file; the builder keeps per-run and skips the combined.
TS_GLOB = f"*task-rest*seg-{ATLAS}_stat-mean_timeseries.tsv"

# ---- outputs (per variant) ----
OUT_BASE = Path(__file__).resolve().parent / "outputs"
OUT_DIR = OUT_BASE / VARIANT
OUT_DIR.mkdir(parents=True, exist_ok=True)
CONN_RAW_CSV = OUT_DIR / f"df_conn_raw_{VARIANT}.csv"

# ---- harmonization covariates (used by step_03_combat / step_04_regress) ----
# Single-sourced from the root config.py so the covariate policy stays identical
# to the rest of the project. Loaded by path because this file is also named
# config.py (a bare `import config` would self-import). FUNC_* = the functional
# covariate set (cohort + mean_fd), distinct from the EEG set.
import importlib.util as _ilu
_root = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "config.py"))
_spec = _ilu.spec_from_file_location("_project_config", _root)
_pc = _ilu.module_from_spec(_spec)
_spec.loader.exec_module(_pc)

BATCH_COL = "machine_batch"          # ComBat batch = scanner (step_02_build_cohort assigns it)
AGE_COL = "age_yrs"
AGE2_COL = "age_yrs_sq"
COHORT_COL = "cohort"
COMBAT_COVARS = _pc.FUNC_COMBAT_COVARS   # biological covariates ComBat PRESERVES
REGRESS_COVARS = _pc.FUNC_REGRESS_COVARS  # nuisance covariates step_04 regresses OUT
