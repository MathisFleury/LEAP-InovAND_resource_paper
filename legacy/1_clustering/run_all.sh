#!/bin/bash
# =============================================================================
# 1_clustering, frozen-only leaf scripts with zero downstream dependents.
# =============================================================================
# These produce paper-reproduction supplementary figures/tables but nothing
# else in the repo reads their output, unlike ../../1_clustering/scripts/
# {01_pca_features.R, 04_run_clustering.py, 05_cluster_stability.py,
# 06_method_comparison.py}, which stayed there because curated mode sources
# or subprocess-calls those same files (they're shared, not separately
# frozen). Run ../../1_clustering/run_all.sh's frozen block FIRST -- these
# read its output (04_run_clustering.py's cluster_assignments/
# individuals_metrics_with_clusters tables in 1_clustering/outputs/).
# =============================================================================
set -e
cd "$(dirname "$0")"

echo "--- 2. Feature selection rationale (frozen) ---"
Rscript 02_feature_selection_rationale.R

echo "--- 3. Cluster validation, NbClust/gap statistic (frozen) ---"
Rscript 03_cluster_validation.R

echo "--- 8. Cluster heatmap, Figure 4a style (frozen) ---"
python3.11 08_cluster_heatmap.py

echo "--- 9. Cluster concept map (frozen) ---"
python3.11 09_cluster_jointplot.py

echo "--- 12 (reval). k-selection stability-based validation ---"
echo "    NOTE: reads the DEPRECATED frozen individuals_metrics.tsv directly"
echo "    (see root CLAUDE.md) -- kept as a historical robustness check, not"
echo "    migrated to curated data."
python3.11 12_reval_kselection.py

echo ""
echo "=== Done. Outputs in legacy/1_clustering/outputs/ ==="
