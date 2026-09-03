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
Rscript 01_pca_features_curated.R
Rscript 01b_sample_size_table.R curated

echo "--- 5. Cluster stability, curated cohort (k-means/Ward/GMM) ---"
python3.11 05_cluster_stability_curated.py
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

echo "--- 1-3. PCA, feature selection, cluster validation (frozen) ---"
Rscript 01_pca_features.R
Rscript 02_feature_selection_rationale.R
Rscript 03_cluster_validation.R

echo "--- 4. Run clustering (k=3), frozen cohort ---"
python3.11 04_run_clustering.py

echo "--- 5. Cluster stability, frozen cohort ---"
python3.11 05_cluster_stability.py

echo "--- 6. Method comparison, frozen cohort ---"
python3.11 06_method_comparison.py

echo "--- 8-9. Cluster heatmap (Figure 4a style) + concept map ---"
python3.11 08_cluster_heatmap.py
python3.11 09_cluster_jointplot.py

echo "--- 7 (reval). k-selection stability-based validation ---"
echo "    NOTE: reads the DEPRECATED frozen individuals_metrics.tsv directly"
echo "    (see root CLAUDE.md) -- kept as a historical robustness check, not"
echo "    migrated to curated data. Flag for removal if no longer needed."
python3.11 07_reval_kselection.py

echo ""
echo "=== Done. Outputs in 1_clustering/outputs/ (curated/ subdir + root) ==="
