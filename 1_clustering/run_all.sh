#!/bin/bash
# Run all clustering analysis scripts -- curated (current, priority regime)
# only. No frozen/deprecated-data reference anywhere in this tree; the
# entire frozen/paper-reproduction pipeline lives in ../legacy/1_clustering/
# with its own self-contained run_all.sh.
set -e
cd "$(dirname "$0")/scripts"

echo "--- 1. Load and clean the curated cohort (LEAP+INOVAND+INFOR) ---"
python3.11 1_load_cohort.py

echo "--- 2. PCA on clinical features, justifies the IQ/SRS feature choice ---"
Rscript 2_pca_features.R
Rscript 3_sample_size_table.R curated

echo "--- 4. Run clustering (k=3), curated cohort ---"
python3.11 4_run_clustering.py

echo "--- 5. Cluster stability, curated cohort (k-means/Ward/GMM) ---"
python3.11 5_cluster_stability.py
python3.11 6_stability_figure_merged.py

echo "--- 7. Method comparison (k-means vs Ward vs GMM), curated ---"
python3.11 7_method_comparison.py curated
Rscript 8_method_comparison_table.R curated

echo "--- 9. Autism-only clustering sensitivity (Reviewer #4.6) ---"
python3.11 9_autism_only_clustering.py
Rscript 10_autism_only_table.R

echo "--- 11. Cluster scatter panels (IQ x SRS, 4 colourings) ---"
python3.11 11_cluster_scatter_panels.py

echo ""
echo "=== Done. Outputs in 1_clustering/outputs/ (cluster_autism/ for the autism-only sensitivity) ==="
echo "For the frozen/paper-reproduction pipeline: cd ../legacy/1_clustering && ./run_all.sh"
