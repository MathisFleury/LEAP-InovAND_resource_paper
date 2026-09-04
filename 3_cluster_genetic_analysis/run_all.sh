#!/usr/bin/env bash
# Cluster-level genetic analysis, gnomAD v4 + curated k-means clustering --
# the data used in the paper. Every other clustering-method variant
# (legacy/gmm/manual/spark_lof, hg19/hg38-v2) is superseded and lives in
# ../legacy/3_cluster_genetic_analysis/ with its own run_all.sh.
set -e
cd "$(dirname "$0")/scripts"

echo "=== 1. Carrier frequencies & odds ratios by clinical cluster (gnomAD v4) ==="
python3.11 1_carrier_freq_or_clusters.py

echo ""
echo "=== 2. Assemble gnomAD v4 interaction tables (curated clusters + v4 carriers + v4 PGS; PAN & EUR) ==="
python3.11 2_interaction_v4_inputs.py

echo ""
echo "=== 3. Common (PGS) x rare (carrier) interaction figures, by cluster/population (v4, PAN & EUR) ==="
ANCESTRY=PAN Rscript 3_interaction_common_rare.R
ANCESTRY=EUR Rscript 3_interaction_common_rare.R

echo ""
echo "=== 4. Cluster / LOEUF IQ x PGS-intelligence panel ==="
echo "    (needs ../2_genetic_analysis/scripts/2_iq_pgs_loeuf_figures.py run first)"
Rscript 4_plot_iq_pgs_cluster_loeuf.R

echo ""
echo "=== Done. Outputs in 3_cluster_genetic_analysis/outputs/ ==="
echo "For the legacy hg19/hg38-v2/non-kmeans pipeline: cd ../legacy/3_cluster_genetic_analysis && ./run_all.sh"
