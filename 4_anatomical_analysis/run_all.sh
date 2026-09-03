#!/usr/bin/env bash
# =============================================================================
# Run all anatomical analysis scripts.
# FROZEN / paper-reproduction pipeline (below) reproduces the original
# submission's figures from the ComBat-only table. CURATED (current, priority
# regime) reruns everything on the QC+ComBat+regressed table, plus the
# reviewer-driven analyses added since (Euler QC, clinical x MRI, Figure 7
# composite) -- see curated/scripts/run_anatomical_analysis_v2.py.
# =============================================================================
set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/scripts" && pwd)"
CURATED_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/curated/scripts" && pwd)"

echo "=== CURATED regime (current, priority) ==="
/usr/local/bin/python3.11 "$CURATED_DIR/run_anatomical_analysis_v2.py"

echo ""
echo "=== FROZEN / paper-reproduction regime (optional) ==="

echo "=== 1. Autism/TD anatomical MRI analysis ==="
/usr/local/bin/python3.11 "$SCRIPT_DIR/01_anatomical_mri_autism_td_analysis.py"

echo "=== 2. Structural MRI vs IQ/SRS ==="
/usr/local/bin/python3.11 "$SCRIPT_DIR/02_structural_mri_iq_srs_analysis.py"

echo "=== 3. Brain visualizations (ggseg, R) ==="
Rscript "$SCRIPT_DIR/03_plot_anatomical_mri_brain_visualizations.R"

echo "=== 3b. Combined single-page brain figure (4 features, FDR, ROI labels) ==="
Rscript "$SCRIPT_DIR/13_combined_brain_figure.R"

echo "=== 4. Subcortical yabplot ==="
/usr/local/bin/python3.11 "$SCRIPT_DIR/04_plot_subcortical_yabplot.py"

echo "=== 5. LOEUF vs MRI correlations (hg19) ==="
/usr/local/bin/python3.11 "$SCRIPT_DIR/05_loeuf_mri_correlations.py"

echo "=== 7. LOEUF vs MRI regression (hg19) ==="
/usr/local/bin/python3.11 "$SCRIPT_DIR/07_carrier_mri_analysis.py"

echo "=== 8. LOEUF vs MRI regression (hg38) ==="
/usr/local/bin/python3.11 "$SCRIPT_DIR/08_hg38_carrier_mri_analysis.py"

echo "=== 10. Standardized beta in significant regions ==="
/usr/local/bin/python3.11 "$SCRIPT_DIR/10_beta_coeff_sign_regions.py"

echo "=== 11. LOEUF vs MRI interaction — scatter, forest, IQ-SRS (hg19 & hg38) ==="
/usr/local/bin/python3.11 "$SCRIPT_DIR/11_loeuf_mri_interaction.py"

echo "=== 12. LOEUF vs MRI brain maps — ggseg, delloF + Miss, SynGO + CHROM (R) ==="
Rscript "$SCRIPT_DIR/12_loeuf_mri_brain_maps.R"

echo ""
echo "=== All anatomical analyses complete ==="
