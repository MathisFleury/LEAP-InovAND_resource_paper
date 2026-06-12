# LEAP-InovAND Resource

Analysis code for **LEAP-InovAND: a multiscale resource to explore genetics,
brain imaging and clinical data in autism**.

This repository is **code only** — no participant data or generated outputs
are versioned. Each section regenerates its own `outputs/` locally when
run against the appropriate data directory.

## Repository layout

```
1_clustering/                # Clinical clustering (IQ × SRS, k=3 K-means)
2_genetic_analysis/          # Rare-variant carrier annotation & ORs (population)
3_cluster_genetic_analysis/  # Rare-variant ORs and PGS by cluster
4_anatomical_analysis/       # Structural MRI: autism vs TD, IQ/SRS, LOEUF
4-1_anatomical_analysis/     # Anatomical analysis (demo-only QC variant)
5_cluster_anatomical_analysis/  # Structural MRI by cluster
6_functional_analysis/       # fMRI connectivity: autism vs TD (Schaefer + subcortical)
7_cluster_functional_analysis/  # fMRI connectivity by cluster
8_eeg_analysis/              # EEG alpha peak & power bands: autism vs TD
9_cluster_eeg_analysis/      # EEG by cluster
10_clinical_analysis/        # Psychomotor milestones & verbal status by cluster
_resources/                  # Custom subcortical atlas builder
```

Each section follows the same convention:

```
<section>/
├── run_all.sh        # entry point — runs every step in order
├── scripts/          # numbered R / Python scripts
└── outputs/          # generated figures/, tables/ (not tracked)
```

## Data

All scripts read participant-level data from a directory pointed to by
environment variables. Set these once in your shell:

```bash
export LEAP_INOVAND_DATA="/path/to/data/directory"
export LEAP_GENETIC_DATA="/path/to/genetic/files"     # optional, defaults to LEAP_INOVAND_DATA
export LEAP_GENELIST_FILE="/path/to/hgnc_genelist.txt"  # optional
```

Cluster-aware sections (3, 5, 7, 9, 10) read cluster labels written by
`1_clustering/scripts/04_run_clustering.py`. Run section 1 first.

## Dependencies

- **Python 3.11** (`/usr/local/bin/python3.11`):
  `pandas numpy scipy scikit-learn statsmodels matplotlib seaborn openpyxl`
- **R**:
  `readr FactoMineR factoextra NbClust mclust cluster fpc ggplot2 ggseg ggsegSchaefer corrplot dplyr gt`

---

## 1. Clustering

Clinical clustering based on IQ and SRS-2 dimensions.

- **PCA** on a broad set of clinical measures (IQ subscales, SRS-2,
  Vineland, RBS-R, SSP) to identify primary variance components.
- **Feature selection**: IQ and SRS maximise sample size (n=1,023 after
  cohort matching) and contribute most to PC1/PC2.
- **Cluster validation** (k=2–10): NbClust (26 indices), AIC/BIC/ICL,
  Silhouette, VRS, Davies-Bouldin, Pseudo-F, Gap statistic.
- **Final clustering**: k=3 on standardised IQ × SRS via **K-means**
  (`n_init=25`, `random_state=42`). The original paper used hierarchical
  Ward; the present resource adopts K-means as the primary method based
  on a reviewer-requested comparative analysis (see
  `06_method_comparison.py`) — k-means is more stable under bootstrap
  (per-subject modal probability 95 % vs 80 %; bootstrap ARI ⟨0.96⟩ vs
  ⟨0.72⟩) and produces better internal indices (silhouette 0.45 vs 0.41,
  Davies-Bouldin 0.83 vs 0.92, Calinski-Harabasz 1097 vs 929). Ward is
  retained as a sensitivity analysis.

```bash
cd 1_clustering/scripts
Rscript 01_pca_features.R                  # PCA on clinical features
Rscript 02_feature_selection_rationale.R   # sample sizes, PCA loadings
Rscript 03_cluster_validation.R            # NbClust, Gap, etc.
python3.11 04_run_clustering.py            # k=3, K-means + Ward sensitivity
python3.11 05_cluster_stability.py         # subsampling, noise, bootstrap
python3.11 06_method_comparison.py         # K-means vs Ward (reviewer-requested)
```

Curated variants (autism-only, by recruitment wave) live alongside as
`04_run_clustering_curated.py` plus `run_downstream_with_curated.py` /
`run_anatomical_with_other_methods.py` / `run_genetics_with_other_methods.py`.

## 2. Genetic analysis

Rare-variant carrier analysis across gene lists (HCNDD, SPARK/SFARI,
EAGLE, SynGO, ChromEpiTF) and ancestries.

- **Carrier annotation**: deletions (validated CNVs), LoF (LOFTEE HC),
  pathogenic missenses (AlphaMissense) on constrained genes.
- **Carrier frequencies & odds ratios**: Fisher exact vs TD, FDR-
  corrected, stratified by ancestry (EUR/non-EUR), cohort
  (InovAND/LEAP), and population grouping.
- **hg38 carrier annotation & ORs**: same pipeline rerun on hg38
  coordinates (`01b`, `04`).

```bash
cd 2_genetic_analysis && ./run_all.sh
```

See `2_genetic_analysis/README.md` for the full methodology and a
script-by-script output table.

## 3. Cluster-level genetic analysis

Carrier frequencies, odds ratios, and PGS by clinical cluster, across
multiple clustering choices (Ward / K-means / GMM / manual), NT pool
definitions (full TD vs C1-only), and genome builds (hg19 / hg38).
Includes external SPARK LoF replication.

```bash
cd 3_cluster_genetic_analysis && ./run_all.sh
```

## 4. Anatomical MRI

Structural MRI (FreeSurfer-derived) for autism vs TD, IQ/SRS correlates,
LOEUF carrier effects, and brain visualisations (ggseg + subcortical
yabplot).

```bash
cd 4_anatomical_analysis && ./run_all.sh
```

`4-1_anatomical_analysis/` is a parallel pipeline using the
`freesurfer_zscore_demo_only.tsv` QC variant (z-scored, demo-only). A
revision pipeline (`4_anatomical_analysis/revision/`) adds per-population
Euler-number checks and hg38 LOEUF maps.

## 5. Cluster anatomical MRI

Structural MRI contrasts by cluster (t-stats, Cohen's d) plus ggseg and
subcortical visualisations. Variants cover curated/curated-autism cohorts
across Ward/K-means/GMM clusterings.

```bash
cd 5_cluster_anatomical_analysis/scripts
python3.11 run_cluster_anatomical_analysis.py
# Optional: restrict NT pool to NT subjects with Cluster == 'C1'
python3.11 run_cluster_anatomical_analysis.py --nt-c1
```

## 6. Functional MRI

Functional connectivity, autism vs TD:

- `preprocessing/` — 7-step pipeline (load → QC → batch → connectivity →
  ComBat → regress → z-score → save).
- `scripts/` — connectivity analyses (full, 6-min sensitivity, covariate-
  adjusted), Schaefer atlas + subcortical visualisations, sensitivity
  comparison plots.

```bash
cd 6_functional_analysis/preprocessing
python3.11 run_pipeline.py
cd ../scripts
python3.11 run_functional_analysis.py
```

## 7. Cluster functional MRI

Cluster-wise connectivity (inputs, Schaefer plots, network matrices,
hyper-/hypo-connectivity counts, subcortical yabplot).

```bash
cd 7_cluster_functional_analysis/scripts
python3.11 run_cluster_functional_analysis.py
```

## 8. EEG

Resting-state EEG, autism vs TD: alpha peak frequency and multi-band
power.

```bash
cd 8_eeg_analysis/scripts
python3.11 run_eeg_analysis.py
```

## 9. Cluster EEG

Same alpha-peak / band-power analyses by cluster.

```bash
cd 9_cluster_eeg_analysis/scripts
python3.11 run_cluster_eeg_analysis.py
```

## 10. Clinical (developmental) analysis

Psychomotor milestone ages and verbal-status proportions by cluster.

```bash
cd 10_clinical_analysis && ./run_all.sh
```

## Citation

If you use this resource, please cite:

> LEAP-InovAND: A multiscale resource to explore genetics, brain imaging, and clinical data in autism.
> medRxiv preprint (2025). https://doi.org/10.1101/2025.11.24.25340858
