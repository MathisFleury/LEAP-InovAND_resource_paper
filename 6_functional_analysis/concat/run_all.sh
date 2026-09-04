#!/bin/bash
# =============================================================================
# 6-2 functional analysis: preprocessing + core analysis suite + Figure 6b/3
# brain-map chain + QC, all on the manuscript's primary regime (nogsr_concat:
# run-concatenated, no-GSR, >6 min -- verified n=393 autism/327 NT).
#
# This is the primary entry point; run_gsr_all_modes.sh separately sweeps the
# other 8 sensitivity variants (gsr, mixed, regfirst, per-cohort, anatQCpass).
# =============================================================================
set -e
SEC="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PRE="$SEC/preprocessing"
ANA="$SEC/scripts"
PY=/usr/local/bin/python3.11
OUT="$SEC/outputs/nogsr_concat"
FMRI_CSV="$PRE/outputs/nogsr_concat/df_conn_cohort_norm_nogsr_concat.csv"

echo "=== Preprocessing: build + harmonise nogsr_concat connectivity ==="
cd "$PRE"
[ -f "outputs/nogsr_concat/df_conn_raw_nogsr_concat.csv" ] || $PY step_01_build_connectivity_concat.py
CONCAT=1 $PY run_pipeline.py

echo ""
echo "=== Core analysis suite (01-07): edge/network-block/NBS tests + figures ==="
cd "$ANA"
AUTISM_TD_FMRI_CSV="$FMRI_CSV" AUTISM_TD_OUTPUT_DIR="$OUT" $PY run_functional_analysis.py

echo ""
echo "=== Figure 6b/3 brain-map chain (08-11) ==="
NBS_SUBDIR=nbs_wholegroup_nogsr_concat FMRI_CONN_FILE="$FMRI_CSV" $PY 08_nbs_brainmap_edges.py
# 10 reads/writes outputs/<v>/nbs_double/ directly (04's NBS-component edges),
# per its own docstring and 11's input contract -- NOT AUTISM_TD_OUTPUT_DIR="$OUT"
# (that silently fed it the FDR-significant edges from 01 instead, a near-empty
# subcortical inset: 1 hyper/0 hypo region vs. the intended 21 hyper/3 hypo).
AUTISM_TD_OUTPUT_DIR="$OUT/nbs_double" $PY 10_nbs_brainmap_subcortical.py
NBS_SUBDIR=nbs_wholegroup_nogsr_concat VARIANTS=nogsr_concat Rscript 09_nbs_brainmap_cortical.R
VARIANTS=nogsr_concat Rscript 11_nbs_combined_figure.R

echo ""
echo "=== QC (12-15): motion confound, sequence/ComBat SVM, FD x connectivity, network proportions ==="
$PY 12_motion_confound_check.py
$PY 13_sequence_combat_svm.py
$PY 14_fd_connectivity_all_edges.py
$PY 15_nbs_network_proportions.py

echo ""
echo "=== Done. Outputs in 6_functional_analysis/concat/outputs/nogsr_concat/ ==="
