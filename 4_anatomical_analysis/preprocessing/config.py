"""
Configuration for the NDD-focused anatomical z-scoring pipeline:
INFOR + INOVAND + LEAP only, with scanner-level ComBat batches.

Ported in-repo from eeg_mri-pipeline/analysis/mri/zscore_autism_td_analysis/
pipeline_ndd (2026-09), so the raw-data -> z-scored-table pipeline lives next
to the analysis scripts that consume it, matching the 6_functional_analysis
and 8_eeg_analysis preprocessing/ pattern.

Differences vs. the all-cohort `pipeline/` (frozen, not ported):
  • COHORT_KEEP filters to INFOR / INOVAND_* / LEAP_W* only
  • MIN_BATCH_N = 2 (minimum needed for ComBat variance) — keep every site

Output stays at the existing $IMG5 canonical location (unchanged) — root
config.py's MRI_CURATED already points here and is consumed by sections
4/5/11, so moving the output in-repo (unlike 6_/8_'s local preprocessing/
outputs/) would ripple far beyond this port. Only the CODE moved in-repo.
"""
import os
from pathlib import Path

# Root config.py = single source of truth. Loaded by path (not `import config`)
# because THIS file is also named config.py — a bare import would self-import.
import importlib.util as _ilu
_root = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "config.py"))
_spec = _ilu.spec_from_file_location("_project_config", _root)
_pc = _ilu.module_from_spec(_spec)
_spec.loader.exec_module(_pc)

# --------------------------------------------------------------------- paths
# $IMG5 must be mounted. Override root via env for a different mount point.
IMG5_ROOT = Path(os.environ.get("IMG5_ROOT", "/Volumes/Imaging5/EEG_MRI-MF"))

INPUT_TSV = (
    IMG5_ROOT / "ALL" / "results" / "tabular" / "anat" / "freesurfer" / "8.1.0" /
    "output" / "ALL_QC_FS_clinical_ML_corrected.tsv"
)
# Per-cohort curated clinical files (preferred over the merged ALL file).
# LEAP-t1 reuses the same curated file root config.py already names (verified
# identical path/content). INOVAND deliberately uses the PLAIN (non-"_merged")
# file here, NOT root config.py's INOVAND_CLINICAL_CURATED (= "_merged") — this
# anatomical pipeline was built and its canonical qc1 output was last frozen
# against the plain file; switching to "_merged" pulls in ~370 additional
# INOVAND_T1 matches and changes N by +26% (1424 -> 1795 subjects). That is a
# real sample-size/results change, not a path cleanup, so it is left exactly as
# the original pipeline had it. Revisit deliberately (not as a silent side
# effect of this port) if the project wants the anat pipeline on "_merged" too.
CLINICAL_TSVS = [
    IMG5_ROOT / "INFOR" / "_clinical_data" / "curated" / "INFOR_clinical_curated.tsv",
    IMG5_ROOT / "INOVAND" / "_clinical_data" / "curated" / "INOVAND_clinical_curated.tsv",
    Path(_pc.LEAP_CLINICAL_CURATED),
    IMG5_ROOT / "LEAP" / "_clinical_data" / "curated" / "LEAP_clinical_curated_t2.tsv",
    IMG5_ROOT / "LEAP" / "_clinical_data" / "curated" / "LEAP_clinical_curated_t3.tsv",
]


def load_clinical(**kwargs):
    """Read + concat the per-cohort curated clinical TSVs into one frame."""
    import pandas as pd
    return pd.concat(
        [pd.read_csv(f, sep="\t", **kwargs) for f in CLINICAL_TSVS],
        ignore_index=True,
    )


OUT_DIR = (
    IMG5_ROOT / "ALL" / "results" / "tabular" / "anat" /
    "z_scoring_qc+combat+regression_ndd" / "output"
)
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ----------------------------------------------------------- cohort filter
# Keep only NDD cohorts of interest. Wave/timepoint variants included.
COHORT_KEEP = {
    "INFOR",
    "INOVAND_T1", "INOVAND_T2",
    "LEAP_W1", "LEAP_W2", "LEAP_W3",
}

# ------------------------------------------------------------------- columns
ID_COLS      = ["ID", "MRI_ID", "cohort", "session_id"]
BATCH_COL    = "batch"           # built in step_02b from machine info
QC_COL       = "QC_seg"
SEX_COL      = "Sex"
AGE_COL      = "age_yrs"
AGE2_COL     = "age2"
ETIV_COL     = "eTIV"
LH_HOLES_COL = "lhSurfaceHoles"
RH_HOLES_COL = "rhSurfaceHoles"

# ---------------------------------------------------------------- QC policy
# Override per-run with env var, e.g. QC_KEEP=1,2 -> keep grades 1 and 2.
QC_KEEP        = [float(x) for x in os.environ.get("QC_KEEP", "1,2,3,4").split(",")]
MIN_BATCH_N    = 2               # ComBat-minimum — keep every scanner

# Multi-run dedup: when a subject (MRI_ID) has several runs/waves, keep the
# single BEST run — lowest QC_seg grade, tie-broken by highest mean_euler,
# then earliest wave. Set False to keep every run.
KEEP_BEST_RUN = True

QC_TAG         = "qc" + "".join(str(int(v)) for v in sorted(QC_KEEP))

# ----------------------------------------------- regression fitting policy
# Which subset the nuisance-regression beta (and the z-score reference) are FIT
# on. Defaults to the project-wide policy in root config.py (DX_COL/CONTROL_VALUE
# = None/None -> fit on the FULL sample, no group privileged as reference; this
# is normative-modelling-neutral, so group means stay symmetric about 0).
# To fit on a subgroup instead (e.g. the old NT-normative behaviour), set both
# env vars: FIT_REFERENCE_COL=population_group  FIT_REFERENCE_VALUE=NT.
# NB: fitting on a subgroup is NOT scale-invariant — if the reference group
# differs from the test groups in a covariate (eTIV/age), it inflates the
# volumetric contrasts (see 5_cluster_anatomical_analysis report Sec 2.5).
FIT_REFERENCE_COL   = os.environ.get("FIT_REFERENCE_COL")   or _pc.DX_COL
FIT_REFERENCE_VALUE = os.environ.get("FIT_REFERENCE_VALUE") or _pc.CONTROL_VALUE

# --------------------------------------------------------------- covariates
# When USE_EULER_COVARIATE = True the pipeline reads `mean_euler` (mean of
# lh + rh orig.nofix Euler numbers) from EULER_TSV — produced by
# build_euler_table.py walking the recon-all logs — and adds it as a quality
# covariate to ComBat and the regression. Subjects without a parseable
# recon-all log are dropped. Env-overridable (USE_EULER=0 for the no-Euler
# sensitivity variant, which is saved to a `_noeuler`-tagged filename).
USE_EULER_COVARIATE = os.environ.get("USE_EULER", "1") == "1"
EULER_TSV           = OUT_DIR / "euler_numbers_recon_all.tsv"
EULER_COL           = "mean_euler"
# Filename tag so the no-Euler variant never clobbers the main output.
EULER_TAG           = "" if USE_EULER_COVARIATE else "_noeuler"

# When REGRESS_AGE = False, age/age2 are still passed to ComBat (protecting
# real age signal from being absorbed into the batch/site correction) but are
# dropped from the per-feature nuisance regression, so age-related variance
# stays in the output. Needed for age-stratified (age-bin) robustness checks:
# testing whether a group effect holds across age bins is close to meaningless
# on a table where the linear age effect has already been regressed out
# globally. Env-overridable (REGRESS_AGE=0), saved to a `_noage`-tagged file.
REGRESS_AGE         = os.environ.get("REGRESS_AGE", "1") == "1"
AGE_TAG             = "" if REGRESS_AGE else "_noage"
VARIANT_TAG         = EULER_TAG + AGE_TAG

_BASE_COVARS      = [AGE_COL, AGE2_COL, SEX_COL, ETIV_COL]
_EXTRA            = [EULER_COL] if USE_EULER_COVARIATE else []
COMBAT_COVARS     = _BASE_COVARS + _EXTRA
_REGRESS_BASE     = _BASE_COVARS if REGRESS_AGE else [SEX_COL, ETIV_COL]
REGRESS_COVARS    = _REGRESS_BASE + _EXTRA

# ---------------------------------------------------------------- features
FEATURE_PREFIX = ""

FS_SUFFIXES = (
    "_NVoxels", "_Volume_mm3", "_normMean", "_normStdDev", "_normMin",
    "_normMax", "_normRange", "_NumVert", "_SurfArea", "_GrayVol",
    "_ThickAvg", "_ThickStd", "_MeanCurv", "_GausCurv", "_FoldInd",
    "_CurvInd",
)
FS_GLOBALS = {
    "BrainSegVol", "BrainSegVolNotVent", "VentricleChoroidVol",
    "lhCortexVol", "rhCortexVol", "CortexVol",
    "lhCerebralWhiteMatterVol", "rhCerebralWhiteMatterVol",
    "CerebralWhiteMatterVol", "SubCortGrayVol", "TotalGrayVol",
    "SupraTentorialVol", "SupraTentorialVolNotVent", "MaskVol",
    "BrainSegVol-to-eTIV", "MaskVol-to-eTIV", "eTIV", "sTIV",
    "lh_NumVert", "lh_WhiteSurfArea", "lh_MeanThickness",
    "rh_NumVert", "rh_WhiteSurfArea", "rh_MeanThickness",
}
FEATURE_EXCLUDE = {
    # eTIV is kept as a feature — step_05 regresses it without self-reference.
    "lhSurfaceHoles", "rhSurfaceHoles", "SurfaceHoles",
    "lh_eTIV", "rh_eTIV",
    "lh_BrainSegVol", "rh_BrainSegVol",
    "lh_BrainSegVolNotVent", "rh_BrainSegVolNotVent",
    "lh_BrainSegVolNotVentSurf", "rh_BrainSegVolNotVentSurf",
    "lh_CortexVol", "rh_CortexVol",
    "lh_SupraTentorialVol", "rh_SupraTentorialVol",
    "lh_SupraTentorialVolNotVent", "rh_SupraTentorialVolNotVent",
}
