#!/bin/bash
# Run all clustering analysis scripts.
# Curated (current, priority regime per project convention) runs by default;
# the original frozen/paper-reproduction pipeline is a separate, optional
# block below it -- kept only to reproduce the published figures, not the
# path new work should build on.
set -e
cd "$(dirname "$0")/scripts"

echo "=== CURATED regime (default) ==="

echo "--- 4. Run clustering (k=3), curated cohort ---"
python3.11 04_run_clustering_curated.py

echo "--- 1. PCA on clinical features, curated cohort ---"
Rscript 01_pca_features.R curated
Rscript 01b_sample_size_table.R curated

echo "--- 5. Cluster stability, curated cohort (k-means/Ward/GMM) ---"
python3.11 05_cluster_stability.py curated
python3.11 05b_stability_figure_merged.py

echo "--- 6. Method comparison (k-means vs Ward vs GMM), curated ---"
python3.11 06_method_comparison.py curated
Rscript 06b_method_comparison_table.R curated

echo "--- 7. Autism-only clustering sensitivity (Reviewer #4.6) ---"
python3.11 07_autism_only_clustering.py
Rscript 07b_autism_only_table.R

echo "--- 10. Cluster scatter panels (IQ x SRS, 4 colourings) ---"
python3.11 10_cluster_scatter_panels.py

echo ""
echo "=== FROZEN / paper-reproduction regime (optional -- reproduces the"
echo "    original submission's figures; not the default for new work) ==="
echo "    Only the genuinely shared scripts run here (curated mode sources/"
echo "    subprocess-calls these same files). Pure frozen-only leaf scripts"
echo "    (feature selection rationale, cluster validation, heatmap, concept"
echo "    map, reval k-selection) moved to ../legacy/1_clustering/ -- run"
echo "    that folder's own run_all.sh separately if you need them."

echo "--- 1. PCA (frozen) ---"
Rscript 01_pca_features.R

echo "--- 4. Run clustering (k=3), frozen cohort ---"
python3.11 04_run_clustering.py

echo "--- 5. Cluster stability, frozen cohort ---"
python3.11 05_cluster_stability.py

echo "--- 6. Method comparison, frozen cohort ---"
python3.11 06_method_comparison.py

echo ""
echo "=== Done. Outputs in 1_clustering/outputs/ (curated/ subdir + root) ==="
