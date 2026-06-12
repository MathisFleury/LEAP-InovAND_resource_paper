# LEAP-InovAND Resource

This repository contains analysis code for our last paper: LEAP-InovAND a multiscale resource to explore genetics, brain imaging and clinical data in autism.

## Overview

Analysis is organized by paper section. Each analysis folder contains:
- `scripts/` — R and Python scripts
- `outputs/` — generated figures and tables
  - `figures/` — PDF figures for publication
  - `tables/` — CSV/TSV tables

## 1. Clustering Analysis

Clinical clustering based on IQ and SRS-2 dimensions, following the methodology described in the paper:

- **PCA** on a broad set of clinical measures (IQ subscales, SRS-2, Vineland, RBS-R, SSP) to identify primary variance components
- **Feature selection** rationale: IQ and SRS maximize sample size (n=1,023 after cohort matching) and contribute most to PC1/PC2
- **Cluster validation** (k=2–10): NbClust (26 indices), AIC/BIC/ICL, Silhouette, VRS, Davies-Bouldin, Pseudo-F, Gap statistic
- **Final clustering**: k=3 on standardised IQ × SRS via **K-means** (`n_init=25`, `random_state=42`). The original paper used hierarchical Ward; the present resource adopts K-means as the primary method based on a reviewer-requested comparative analysis (see `06_method_comparison.py`) — k-means is more stable under bootstrap (per-subject modal probability 95 % vs 80 %; bootstrap ARI ⟨0.96⟩ vs ⟨0.72⟩) and produces better internal indices (silhouette 0.45 vs 0.41, Davies-Bouldin 0.83 vs 0.92, Calinski-Harabasz 1097 vs 929). Ward is retained as a sensitivity analysis.

### Scripts (run in order)

```bash
cd 1_clustering/scripts

# 1. PCA on clinical features
Rscript 01_pca_features.R

# 2. Feature selection rationale (sample sizes, PCA loadings)
Rscript 02_feature_selection_rationale.R

# 3. Cluster number validation (NbClust, Gap, etc.)
Rscript 03_cluster_validation.R

# 4. Run final clustering (k=3, K-means primary + Ward sensitivity)
python 04_run_clustering.py

# 5. Cluster stability (subsampling, noise, bootstrap, sample-size)
python 05_cluster_stability.py

# 6. Method comparison (K-means vs Ward, reviewer-requested)
python 06_method_comparison.py
```

### Outputs

| Output | Description |
|-------|-------------|
| `figures/PCA_variance.pdf` | Scree plot |
| `figures/PCA_loadings.pdf` | Variable loadings (Supplementary Fig. 21 / S20) |
| `figures/PCA_contrib_PC1.pdf`, `PCA_contrib_PC2.pdf` | Variable contributions |
| `figures/sample_size_completeness.pdf` | Sample size by variable combination |
| `figures/gap_statistic.pdf` | Gap statistic for k selection |
| `figures/NbClust_consensus.pdf` | NbClust index consensus |
| `figures/silhouette_k3.pdf` | Silhouette plot for k=3 |
| `figures/clustering_scatter.pdf` | IQ × SRS scatter, K-means primary |
| `figures/clustering_scatter_ward.pdf` | IQ × SRS scatter, Ward sensitivity |
| `figures/cluster_stability.pdf` | 4-panel stability check (K-means) |
| `figures/method_comparison.pdf` | K-means vs Ward comparison |
| `tables/PCA_variable_importance.csv` | PCA loadings and contributions |
| `tables/sample_sizes_by_variable_combination.csv` | Supplementary Table 12 |
| `tables/cluster_validation_indices.csv` | Supplementary Table 13 |
| `tables/cluster_assignments.csv` | Individual cluster labels (K-means primary) |
| `tables/cluster_assignments_ward.csv` | Ward labels (sensitivity) |
| `tables/stability_subsampling.csv`, `stability_noise.csv`, `stability_per_individual.csv`, `stability_sample_size.csv` | Stability tables |
| `tables/method_comparison_*.csv` | Internal indices, concordance, bootstrap, Jaccard, modal-prob |

### Data

Scripts expect `individuals_metrics.tsv` in the data directory. Set the path via:

```bash
export LEAP_INOVAND_DATA="/path/to/data/directory"
```

Default path (relative to repository): `../imaging2genet/0_input/dataframes/`

### R Dependencies

```r
install.packages(c("readr", "FactoMineR", "factoextra", "NbClust", "mclust", "cluster", "fpc", "ggplot2", "corrplot", "dplyr", "gt"))
```

### Python Dependencies

```
pandas numpy matplotlib seaborn scikit-learn
```

## 2. Genetic Analysis

Rare variant carrier analysis across multiple gene lists (HCNDD, SPARK/SFARI, EAGLE, SynGO, ChromEpiTF) and clinical clusters:

- **Carrier annotation**: Deletions (validated CNVs), LoF (LOFTEE HC), and pathogenic missenses (AlphaMissense) on constrained genes
- **Carrier frequencies and odds ratios**: Fisher exact tests vs TD, FDR-corrected, stratified by ancestry (EUR/non-EUR), cohort (InovAND/LEAP), and population grouping
- **Cluster-level analysis**: Carrier frequencies and OR by clinical clusters (C1–C3, TD, IDD)
- **Polygenic scores**: PGS odds ratios (top 25% vs bottom 25%) by cluster for autism, ADHD, intelligence, and other traits

### Scripts (run in order)

```bash
cd 2_genetic_analysis/scripts

# 1. Carrier annotation (DEL, LoF, missense)
python 01_carrier_annotation.py

# 2. Population-level frequencies & odds ratios
python 02_carrier_freq_or.py

# 3. Cluster-level analysis & PGS
python 03_carrier_freq_or_clusters.py
```

### Data

Set data directory via environment variables:

```bash
export LEAP_INOVAND_DATA="/path/to/data/directory"
export LEAP_GENETIC_DATA="/path/to/genetic/files"    # optional, defaults to LEAP_INOVAND_DATA
export LEAP_GENELIST_FILE="/path/to/hgnc_genelist.txt"  # optional
```

### Python Dependencies

```
pandas numpy matplotlib seaborn scikit-learn scipy statsmodels openpyxl
```

## Citation

If you use this resource, please cite:

> LEAP-InovAND: A multiscale resource to explore genetics, brain imaging, and clinical data in autism.  
> medRxiv preprint (2025). https://doi.org/10.1101/2025.11.24.25340858
