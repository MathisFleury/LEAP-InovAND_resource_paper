#!/usr/bin/env bash
# =============================================================================
# Run the anatomical MRI analysis pipeline -- QC+ComBat+regressed table, the
# data used in the paper (steps 1-21, see scripts/run_anatomical_analysis.py
# for exact ordering/dependencies). No frozen/deprecated-data reference
# anywhere in this tree; the legacy paper-reproduction pipeline (ComBat-only
# table) lives in ../legacy/4_anatomical_analysis/ with its own run_all.sh.
# =============================================================================
set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/scripts" && pwd)"

python3.11 "$SCRIPT_DIR/run_anatomical_analysis.py"

echo ""
echo "=== Done. Outputs in 4_anatomical_analysis/outputs/ ==="
echo "For the legacy paper-reproduction pipeline: cd ../legacy/4_anatomical_analysis && ./run_all.sh"
