#!/usr/bin/env bash
# =============================================================================
# Age x sex interaction analysis ("The Effects of Age and Sex" manuscript
# section) + age-binned sanity checks. This section had no entry point before
# -- run order below matches the internal importlib/Rscript chaining already
# in the scripts (01 is a hard dependency of 05/06/09/10 via dynamic import;
# 10 already shells out to 11, 12 does not shell out to 13).
# =============================================================================
set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/scripts" && pwd)"
cd "$SCRIPT_DIR"
PY=/usr/local/bin/python3.11
# Manuscript's primary rsfMRI regime -- matches 6_functional_analysis/concat's
# and 7_cluster_functional_analysis's run_all.sh. Scripts already default here;
# set explicitly so the data source is self-documenting, not implicit.
FMRI_CSV="$SCRIPT_DIR/../../6_functional_analysis/concat/preprocessing/outputs/nogsr_concat/df_conn_cohort_norm_nogsr_concat.csv"

echo "=== 01. Diagnosis x age / x sex interaction (anatomical + functional models) ==="
$PY 01_age_sex_interactions.py

echo ""
echo "=== 02. Functional NBS for the age/sex interaction ==="
AUTISM_TD_FMRI_CSV="$FMRI_CSV" $PY 02_func_nbs_interaction.py

echo ""
echo "=== 03. Functional interaction double network matrix ==="
Rscript 03_plot_func_interaction_double_matrix.R

echo ""
echo "=== 04. Anatomical interaction ggseg brain maps ==="
Rscript 04_plot_anat_interaction_ggseg.R

echo ""
echo "=== 08. Combined interaction brain figure (age + sex, one page) ==="
Rscript 08_combined_interaction_brain_figure.R

echo ""
echo "=== 05. Age trajectory, posterior cingulate thickness (autism vs NT) ==="
$PY 05_plot_frontal_thickness_age.py

echo ""
echo "=== 06. Age-composition sensitivity (5-22y range + matched sub-sample) ==="
$PY 06_age_range_and_matched_sensitivity.py

echo ""
echo "=== 07. Age by cluster (Kruskal-Wallis + pairwise Mann-Whitney) ==="
$PY 07_age_by_cluster.py

echo ""
echo "=== 09. Age-binned regional effect, headline interaction region ==="
$PY 09_age_bin_regional_effect.py

echo ""
echo "=== 10. Age-binned anatomical brain maps, all features (also runs 11's grid figure) ==="
$PY 10_age_bin_brain_maps.py

echo ""
echo "=== 12. Age-binned functional NBS ==="
AUTISM_TD_FMRI_CSV="$FMRI_CSV" $PY 12_age_bin_func_nbs.py

echo ""
echo "=== 13. Age-binned functional NBS grid figure ==="
Rscript 13_age_bin_func_nbs_grid.R

echo ""
echo "=== Done. Outputs in 11_age_sex_analysis/outputs/ ==="
