#!/bin/bash
# Run the gnomAD v4 population-level genetic analysis pipeline -- the data
# used in the paper. Cluster-level analysis lives in
# 3_cluster_genetic_analysis/ (matches the anatomical/functional section
# convention). No hg19/hg38-v2 reference anywhere in this tree; that legacy
# (superseded) pipeline lives in ../legacy/2_genetic_analysis/ with its own
# run_all.sh.
set -e
cd "$(dirname "$0")/scripts"

echo "=== 1. gnomAD v4 carrier frequencies & odds ratios (population-level) ==="
python3.11 1_carrier_freq_or.py

echo "=== 2. IQ x PGS-intelligence x LOEUF (population-level frame + betas) ==="
python3.11 2_iq_pgs_loeuf_figures.py

echo "=== Done. Outputs in 2_genetic_analysis/outputs/figures/ + tables/ ==="
echo "Cluster-coloured IQ x PGS x LOEUF panel: 3_cluster_genetic_analysis/scripts/4_plot_iq_pgs_cluster_loeuf.R"
echo "Cluster-level genetic analysis: cd ../3_cluster_genetic_analysis && ./run_all.sh"
echo "For the legacy hg19/hg38-v2 pipeline: cd ../legacy/2_genetic_analysis && ./run_all.sh"
