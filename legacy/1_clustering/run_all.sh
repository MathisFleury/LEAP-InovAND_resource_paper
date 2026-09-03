#!/bin/bash
# =============================================================================
# 1_clustering, frozen-only scripts kept out of the main tree.
# =============================================================================
# Two kinds of file here:
#   - 01_pca_features_frozen.R / 05_cluster_stability_frozen.py: thin
#     invokers for the shared PCA/stability engines in
#     ../../1_clustering/scripts/ -- those engines are curated-only now, so
#     these wrappers are what keep the deprecated individuals_metrics.tsv
#     reference (see root CLAUDE.md) out of the main tree.
#   - 02/03/08/09/12: genuinely frozen-only leaf scripts with zero downstream
#     dependents anywhere in the repo.
# Run ../../1_clustering/run_all.sh's frozen block FIRST -- 08/09 read its
# output (04_run_clustering.py's cluster_assignments/
# individuals_metrics_with_clusters tables in 1_clustering/outputs/).
# =============================================================================
set -e
cd "$(dirname "$0")"

echo "--- 1. PCA (frozen) ---"
Rscript 01_pca_features_frozen.R

echo "--- 2. Feature selection rationale (frozen) ---"
Rscript 02_feature_selection_rationale.R

echo "--- 3. Cluster validation, NbClust/gap statistic (frozen) ---"
Rscript 03_cluster_validation.R

echo "--- 5. Cluster stability, frozen cohort (k-means/Ward/GMM) ---"
python3.11 05_cluster_stability_frozen.py

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
echo "=== Done. Outputs in 1_clustering/outputs/ (01/05/04/06) + legacy/1_clustering/outputs/ (02/03/08/09/12) ==="
