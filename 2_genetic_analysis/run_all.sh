#!/bin/bash
# Run all genetic analysis scripts in order
set -e
cd "$(dirname "$0")/scripts"

echo "=== 1. Carrier annotation ==="
python3.11 01_carrier_annotation.py

echo "=== 1b. hg38 carrier annotation ==="
python3.11 01b_hg38_carrier_annotation.py

echo "=== 2. Carrier frequencies & odds ratios (population-level) ==="
python3.11 02_carrier_freq_or.py

echo "=== 3. Cluster-level analysis & PGS ==="
python3.11 03_carrier_freq_or_clusters.py

echo "=== 4. hg38 carrier frequencies & odds ratios ==="
python3.11 04_hg38_carrier_freq_or.py

echo "=== Done. Outputs in 2_genetic_analysis/outputs/ ==="
