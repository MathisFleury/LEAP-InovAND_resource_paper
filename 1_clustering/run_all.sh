#!/bin/bash
# Run all clustering analysis scripts -- curated (current, priority regime)
# only. No frozen/deprecated-data reference anywhere in this tree; the
# entire frozen/paper-reproduction pipeline lives in ../legacy/1_clustering/
# with its own self-contained run_all.sh.
set -e
cd "$(dirname "$0")/scripts"

echo "--- 4. Run clustering (k=3), curated cohort ---"
python3.11 04_run_clustering.py

echo "--- 1. PCA on clinical features, curated cohort ---"
Rscript 01_pca_features.R
Rscript 01b_sample_size_table.R curated

echo "--- 5. Cluster stability, curated cohort (k-means/Ward/GMM) ---"
python3.11 05_cluster_stability.py
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
echo "=== Done. Outputs in 1_clustering/outputs/curated/ ==="
echo "For the frozen/paper-reproduction pipeline: cd ../legacy/1_clustering && ./run_all.sh"
