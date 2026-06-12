#!/usr/bin/env bash
set -e
cd "$(dirname "$0")/scripts"

echo "=== 01. Psychomotor milestones by cluster ==="
python3.11 01_plot_psychomotor_milestones.py

echo ""
echo "=== 02. Verbal status by cluster ==="
python3.11 02_plot_verbal_status.py

echo ""
echo "=== Done. Outputs in 10_clinical_analysis/outputs/ ==="
