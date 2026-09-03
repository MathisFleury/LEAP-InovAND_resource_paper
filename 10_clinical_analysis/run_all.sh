#!/usr/bin/env bash
set -e
cd "$(dirname "$0")/scripts"

echo "=== 01. Psychomotor milestones by cluster ==="
python3.11 01_plot_psychomotor_milestones.py

echo ""
echo "=== 02. Verbal status by cluster ==="
python3.11 02_plot_verbal_status.py

echo ""
echo "=== 03. Clinical scores by cluster (panel C) ==="
python3.11 03_plot_clinical_scores.py

echo ""
echo "=== 04. IQ x SRS scatter (clusters + v4 genetics + adaptive) ==="
python3.11 04_plot_iq_srs_scatter.py

echo ""
echo "=== 04e. Cluster/LOEUF IQ x SRS panel, reframed (feeds 4_anatomical_analysis's Figure 7 composite) ==="
Rscript 04e_plot_iq_srs_cluster_reframe.R

echo ""
echo "=== 05. UpSet plot: data-modality availability by phenotype ==="
python3.11 05_upset_modality.py

echo ""
echo "=== 06. Cross-dataset IQ/SRS/Age by sex (LEAP-InovAND vs SPARK vs ABIDE) ==="
python3.11 06_datasets_iq_srs_age.py

echo ""
echo "=== 08. Cohort size, before (original submission) vs now (curated) ==="
python3.11 08_cohort_size_before_after.py

echo ""
echo "=== 09. Clinical scores, LEAP vs INOVAND, Autism vs NT (paper supp. figure update) ==="
python3.11 09_plot_clinical_scores_by_population.py

echo ""
echo "=== Done. Outputs in 10_clinical_analysis/outputs/ ==="
echo ""
echo "Optional (not run by default):"
echo "  07_legacy_vs_curated_diff.py           -- QA diagnostic, frozen vs curated clustering input"
echo "  10_presentation_iq_srs_build.R         -- oral-presentation slide build, not a paper figure"
echo "  11_presentation_iq_srs_hcndd.R         -- same, HCNDD genetics variant"
