#!/bin/bash
# =============================================================================
# 7. Cluster functional MRI analysis: per-cluster Schaefer / hyper-hypo /
# network-matrix / yabplot chain (01-05, + 01b/01c sensitivity) and the
# Figure 6b per-cluster NBS chain (06-09), on the manuscript's primary regime
# (6_functional_analysis/concat's nogsr_concat -- same connectivity table
# section 6's run_all.sh uses).
# =============================================================================
set -e
SEC="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ANA="$SEC/scripts"
PY=/usr/local/bin/python3.11
FMRI_CSV="$SEC/../6_functional_analysis/concat/preprocessing/outputs/nogsr_concat/df_conn_cohort_norm_nogsr_concat.csv"

cd "$ANA"

echo "=== Core cluster analysis (01-05): Schaefer / hyper-hypo / network-matrix / yabplot ==="
FMRI_CONN_FILE="$FMRI_CSV" $PY 01_generate_cluster_fmri_inputs.py
Rscript 02_plot_cluster_fmri_schaefer.R
$PY 03_cluster_hyper_hypo_analysis.py
Rscript 04_plot_cluster_network_matrix.R
$PY 05_plot_cluster_subcortical_yabplot.py

echo ""
echo "=== Sensitivity checks (01b/01c): 6-min cutoff, covariate-adjusted ==="
FMRI_CONN_FILE="$FMRI_CSV" $PY 01b_generate_cluster_fmri_inputs_6min.py
FMRI_CONN_FILE="$FMRI_CSV" $PY 01c_generate_cluster_fmri_inputs_covariate.py

echo ""
echo "=== Figure 6b chain (06-09): per-cluster NBS ==="
NBS_SUBDIR=nbs_nogsr_concat FMRI_CONN_FILE="$FMRI_CSV" $PY 06_cluster_nbs.py
NBS_SUBDIR=nbs_nogsr_concat $PY 08_figure6b_nbs_subcortical.py
NBS_SUBDIR=nbs_nogsr_concat Rscript 07_figure6b_nbs_cortical.R
NBS_SUBDIR=nbs_nogsr_concat $PY 09_nbs_prominent_networks.py

echo ""
echo "=== Done. Outputs in 7_cluster_functional_analysis/outputs/ ==="
