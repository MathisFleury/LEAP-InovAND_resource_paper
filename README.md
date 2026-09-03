# LEAP-InovAND Resource

Analysis code for **LEAP-InovAND: a multiscale resource to explore genetics,
brain imaging and clinical data in autism**.

This repository is **code only** — no participant data or generated outputs
are versioned. Each section regenerates its own `outputs/` locally when
run against the appropriate data directory.

Two data regimes recur throughout: **curated** (current, priority regime —
what new work and most `run_all.sh` defaults target) and **frozen /
paper-reproduction** (the original submission's data snapshot, kept only to
reproduce the published figures, not the basis for new work). Where a
section keeps both, `run_all.sh` runs curated by default and the frozen
block separately, clearly labelled.

## Repository layout

```
1_clustering/                   # Clinical clustering (IQ × SRS, k=3 K-means)
2_genetic_analysis/             # Rare-variant carrier annotation & ORs (population)
3_cluster_genetic_analysis/     # Rare-variant ORs and PGS by cluster
4_anatomical_analysis/          # Structural MRI: autism vs NT, IQ/SRS, LOEUF
  └── curated/                  #   current pipeline: QC+ComBat+regressed table,
                                 #   Euler QC, clinical x MRI, Figure 7 composite
5_cluster_anatomical_analysis/  # Structural MRI by cluster
  └── curated/                  #   current pipeline (Figure 6a composite, gt tables)
6_functional_analysis/          # fMRI connectivity: autism vs NT (legacy XCP-D 0.8 pipeline)
6-1_functional_analysis_concat/ # fMRI sensitivity: concatenated INOVAND runs (>=6 min)
6-2_functional_analysis/        # fMRI connectivity, CURRENT pipeline (XCP-D 0.11, primary
                                 #   regime = nogsr_concat) -- Figures 3/6b, QC
7_cluster_functional_analysis/  # fMRI connectivity by cluster (incl. Figure 6b chain)
8_eeg_analysis/                 # EEG alpha peak & power bands: autism vs NT
9_cluster_eeg_analysis/         # EEG by cluster
10_clinical_analysis/           # Psychomotor milestones, IQ/SRS/cohort figures, verbal status
11_age_sex_analysis/            # Diagnosis x age / x sex interaction + age-binned checks
_resources/                     # Cross-section figure/table builders + atlas builder
```

Each numbered section follows the same convention:

```
<section>/
├── run_all.sh        # entry point — runs every step in order
├── scripts/          # numbered R / Python scripts
└── outputs/          # generated figures/, tables/ (not tracked)
```

Sections without a `run_all.sh` have a `run_*.py` orchestrator in `scripts/`
instead (noted per section below).

## Data

All scripts read participant-level data from a directory pointed to by
environment variables. Set these once in your shell:

```bash
export LEAP_INOVAND_DATA="/path/to/data/directory"
export LEAP_GENETIC_DATA="/path/to/genetic/files"     # optional, defaults to LEAP_INOVAND_DATA
export LEAP_GENELIST_FILE="/path/to/hgnc_genelist.txt"  # optional
```

Cluster-aware sections (3, 5, 7, 9, 10) read cluster labels written by
`1_clustering/scripts/04_run_clustering_curated.py` (curated) or
`04_run_clustering.py` (frozen). Run section 1 first.

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
- **Feature selection**: IQ and SRS maximise sample size and contribute
  most to PC1/PC2.
- **Cluster validation** (k=2–10): NbClust (26 indices), AIC/BIC/ICL,
  Silhouette, VRS, Davies-Bouldin, Pseudo-F, Gap statistic.
- **Final clustering**: k=3 on standardised IQ × SRS via **K-means**
  (`n_init=25`, `random_state=42`). The original paper used hierarchical
  Ward; the present resource adopts K-means as the primary method based
  on a reviewer-requested comparative analysis (see
  `06_method_comparison.py`) — k-means is more stable under bootstrap and
  produces better internal indices. Ward and GMM are retained as
  sensitivity analyses.
- **Autism-only clustering** (Reviewer #4.6) and a **reval**
  stability-based k-selection check are also included.

```bash
cd 1_clustering && ./run_all.sh
```

Runs the **curated** regime by default (clustering → PCA → stability →
method comparison → autism-only sensitivity → scatter panels), then a
separate, clearly labelled **frozen/paper-reproduction** block (original
01–06, 08–09, plus the reval check). Curated variants (`04_run_clustering_curated.py`,
`05_cluster_stability_curated.py`, …) can also be run individually; three
cross-section batch drivers (`run_anatomical_with_other_methods.py`,
`run_downstream_with_curated.py`, `run_genetics_with_other_methods.py`)
rerun the relevant downstream sections once per clustering method
(K-means / Ward / GMM).

## 2. Genetic analysis

Rare-variant carrier analysis across gene lists (HCNDD, SPARK/SFARI,
EAGLE, SynGO, ChromEpiTF) and ancestries.

- **Carrier annotation**: deletions (validated CNVs), LoF (LOFTEE HC),
  pathogenic missenses (AlphaMissense) on constrained genes — hg19/gnomAD
  v2 (`01`), hg38/gnomAD v2 (`01b`), and gnomAD v4 (`05`, current priority
  regime).
- **Carrier frequencies & odds ratios**: Fisher exact vs NT, FDR-
  corrected, stratified by ancestry (EUR/non-EUR), cohort
  (InovAND/LEAP), and population grouping.
- **IQ x PGS-intelligence x LOEUF** scatter (cluster-coloured).

```bash
cd 2_genetic_analysis && ./run_all.sh
```

See `2_genetic_analysis/README.md` for the full methodology and a
script-by-script output table.

## 3. Cluster-level genetic analysis

Carrier frequencies, odds ratios, and PGS by clinical cluster, across
multiple clustering choices (Ward / K-means / GMM / manual cluster-source
strategies), NT pool definitions (full NT vs C1-only), and genome builds
(hg19 / hg38 / gnomAD v4, the current priority regime — see `15`–`18`).
Includes external SPARK LoF replication.

```bash
cd 3_cluster_genetic_analysis && ./run_all.sh
```

## 4. Anatomical MRI

Structural MRI (FreeSurfer-derived) for autism vs NT, IQ/SRS correlates,
LOEUF carrier effects, and brain visualisations (ggseg + subcortical
yabplot).

```bash
cd 4_anatomical_analysis && ./run_all.sh
```

Runs the **curated** pipeline by default (`curated/scripts/run_anatomical_analysis_v2.py`
— rebuilt on the QC+ComBat+age/sex/eTIV-regressed FreeSurfer table, plus
Euler-number QC, LOEUF hg38 correlations/regression, clinical x MRI
correlations, per-site/age-imbalance robustness checks, and the Figure 7
composite), then the **frozen/paper-reproduction** pipeline (original
ComBat-only table, no QC/regression).

## 5. Cluster anatomical MRI

Structural MRI contrasts by cluster (t-stats, Cohen's d) plus ggseg and
subcortical visualisations, including the Figure 6a composite grid.

```bash
cd 5_cluster_anatomical_analysis/scripts
python3.11 run_cluster_anatomical_analysis.py       # frozen
# Optional: restrict NT pool to NT subjects with Cluster == 'C1'
python3.11 run_cluster_anatomical_analysis.py --nt-c1

cd ../curated/scripts
python3.11 run_cluster_anatomical_analysis_v2.py    # curated (current)
```

## 6. Functional MRI (legacy)

Functional connectivity, autism vs NT, on the original XCP-D 0.8
preprocessing. Superseded for reporting by `6-2_functional_analysis/`
below, but kept for the sensitivity/legacy-comparison checks it still
answers.

```bash
cd 6_functional_analysis/preprocessing
python3.11 run_pipeline.py
cd ../scripts
python3.11 run_functional_analysis.py
```

`run_functional_analysis.py` runs the primary chain (edge-level FDR,
network-block permutation, NBS, brain maps), then a sensitivity block
(6-min/mean-FD cutoff, covariate-adjusted, comparison figures), then a
legacy-vs-revised preprocessing QA comparison. See
`6_functional_analysis/METHODS_functional_6_vs_7.md` for how this
section's inference levels relate to section 7's per-cluster analysis.

## 6-1. Functional MRI, concatenated-run sensitivity

fMRI autism-vs-NT on runs concatenated per subject (≥6 min), for INOVAND
(+ LEAP). Requires the `/Volumes/Imaging5` mount.

```bash
cd 6-1_functional_analysis_concat && ./run_all.sh
```

## 6-2. Functional MRI (current pipeline)

Functional connectivity on XCP-D v0.11, self-contained (no dependency on
section 6). The manuscript's primary regime is `nogsr_concat`
(run-concatenated, no global-signal-regression, >6 min; n=393 autism / 327
NT) — this is what `11_age_sex_analysis` and the reviewer-response numbers
use.

```bash
cd 6-2_functional_analysis && ./run_all.sh
```

Runs preprocessing (build connectivity → ComBat → regress → z-score) on
`nogsr_concat`, the core edge/network-block/NBS analysis suite, the
Figure 6b/3 brain-map chain, and QC (motion confound, sequence/ComBat SVM,
FD x connectivity). `run_gsr_all_modes.sh` separately sweeps the other 8
sensitivity variants (gsr, mixed, regfirst, per-cohort, anatomical-QC-pass)
— run manually, not part of `run_all.sh` (it is a long sweep, not a single
pipeline run).

## 7. Cluster functional MRI

Cluster-wise connectivity (inputs, Schaefer plots, network matrices,
hyper-/hypo-connectivity counts, subcortical yabplot), including the
Figure 6b per-cluster NBS chain.

```bash
cd 7_cluster_functional_analysis/scripts
python3.11 run_cluster_functional_analysis.py
```

## 8. EEG

Resting-state EEG, autism vs NT: alpha peak frequency and multi-band
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

## 10. Clinical analysis

Psychomotor milestones and verbal status by cluster; cohort-composition
figures (UpSet plot, cross-dataset IQ/SRS/age); IQ x SRS x LOEUF scatter
(feeds section 4's Figure 7 composite); cohort-size before/after
accounting; clinical-score distributions by cohort and population.

```bash
cd 10_clinical_analysis && ./run_all.sh
```

## 11. Age x sex interaction analysis

Diagnosis x age and diagnosis x sex interaction models (anatomical +
functional), age-binned sanity checks, and age-composition sensitivity
analyses ("The Effects of Age and Sex" manuscript section). Depends on
`4_anatomical_analysis/curated/` directly (imports its config and script 01).

```bash
cd 11_age_sex_analysis && ./run_all.sh
```

## Citation

If you use this resource, please cite:

> LEAP-InovAND: A multiscale resource to explore genetics, brain imaging, and clinical data in autism.
> medRxiv preprint (2025). https://doi.org/10.1101/2025.11.24.25340858
