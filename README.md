# LEAP-InovAND Resource

Analysis code for **LEAP-InovAND: a multiscale resource to explore genetics,
brain imaging and clinical data in autism**.

This repository is **code only** — no participant data or generated outputs
are versioned. Each section regenerates its own `outputs/` locally when
run against the appropriate data directory.

Two data regimes recur throughout: **curated** (default) and **legacy**
(an earlier pipeline, kept alongside it). Where a section has both,
`run_all.sh` runs curated by default and the legacy block separately.

## Repository layout

```
1_clustering/                   # Clinical clustering (IQ × SRS, k=3 K-means)
2_genetic_analysis/             # Rare-variant carrier annotation & ORs (population)
3_cluster_genetic_analysis/     # Rare-variant ORs and PGS by cluster
4_anatomical_analysis/          # Structural MRI: autism vs NT, IQ/SRS, LOEUF
                                 #   (QC+ComBat+regressed table, Euler QC,
                                 #   clinical x MRI, Figure 7 composite)
5_cluster_anatomical_analysis/  # Structural MRI by cluster (Figure 6a composite, gt tables)
6_functional_analysis/          # fMRI connectivity: autism vs NT
                                 #   concat/     -- CURRENT pipeline (XCP-D 0.11, primary
                                 #                  regime = nogsr_concat) -- Figures 3/6b, QC
                                 #   non_concat/ -- legacy XCP-D 0.8 pipeline, secondary analysis
7_cluster_functional_analysis/  # fMRI connectivity by cluster (incl. Figure 6b chain)
8_eeg_analysis/                 # EEG alpha peak & power bands: autism vs NT
9_cluster_eeg_analysis/         # EEG by cluster
10_clinical_analysis/           # Psychomotor milestones, IQ/SRS/cohort figures, verbal status
11_age_sex_analysis/            # Diagnosis x age / x sex interaction + age-binned checks
legacy/                         # Superseded/frozen pipelines, moved out of their section
                                 #   (local only, not part of this public release)
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
`1_clustering/scripts/4_run_clustering.py` (curated) or the local, not
publicly released, frozen pipeline. Run section 1 first.

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
  on a comparative analysis (see `7_method_comparison.py`) — k-means is
  more stable under bootstrap and produces better internal indices. Ward
  and GMM are retained as sensitivity analyses.
- **Autism-only clustering** and a **reval** stability-based k-selection
  check are also included.

```bash
cd 1_clustering && ./run_all.sh
```

Runs only the **curated** regime (load cohort → PCA → clustering →
stability → method comparison → autism-only sensitivity → scatter panels)
— no frozen/deprecated-data reference anywhere in this tree. Scripts are
numbered in run order (`1`–`11`); PCA (`2`) runs before clustering (`4`)
since both read the cohort table `1_load_cohort.py` writes, so the
feature-justification step doesn't depend on the clustering step's output.
`7_method_comparison.py` still takes a `curated` arg since it's a
comparison and needs both regimes' data, but never touches the deprecated
path directly. Every legacy script lives in the local,
not publicly released, `legacy/1_clustering/` with its own, fully
self-contained `run_all.sh`. Three cross-section batch drivers
(`run_anatomical_with_other_methods.py`,
`run_downstream_with_curated.py`, `run_genetics_with_other_methods.py`)
rerun the relevant downstream sections once per clustering method
(K-means / Ward / GMM).

## 2. Genetic analysis

Rare-variant carrier analysis across gene lists (HCNDD, SPARK/SFARI,
EAGLE, SynGO, ChromEpiTF) and ancestries, on the **gnomAD v4** data used
in the paper. Population-level only — cluster-level analysis lives in
section 3 (matches the anatomical/functional section convention). The
earlier hg19/hg38-v2 pipeline is superseded and lives in
`legacy/2_genetic_analysis/` (local only, not part of this public release).

- **Carrier annotation, frequencies & odds ratios**: deletions (validated
  CNVs), LoF (LOFTEE HC), pathogenic missenses (AlphaMissense) on
  constrained genes; Fisher exact vs NT, FDR-corrected, stratified by
  ancestry (EUR/PAN), constraint definition, and variant combination.
- **IQ x PGS-intelligence x LOEUF** scatter (cluster-coloured).

```bash
cd 2_genetic_analysis && ./run_all.sh
```

See `2_genetic_analysis/README.md` for the full methodology and a
script-by-script output table.

## 3. Cluster-level genetic analysis

Carrier frequencies, odds ratios, and PGS by clinical cluster, on gnomAD v4
data stratified by the current curated k-means clustering. Every other
clustering choice explored previously (Ward/frozen k-means, GMM, manual
cluster-source strategies, NT pool definitions, hg19/hg38 genome builds)
plus external SPARK LoF replication is superseded and lives in
`legacy/3_cluster_genetic_analysis/` (local only, not part of this public
release).

```bash
cd 3_cluster_genetic_analysis && ./run_all.sh
```

## 4. Anatomical MRI

Structural MRI (FreeSurfer-derived) for autism vs NT, IQ/SRS correlates,
LOEUF carrier effects, and brain visualisations (ggseg + subcortical
yabplot), on the QC+ComBat+age/sex/eTIV-regressed FreeSurfer table used in
the paper. Also covers Euler-number QC, LOEUF hg38 correlations/regression,
clinical x MRI correlations, per-site/age-imbalance robustness checks, and
the Figure 7 composite (`scripts/run_anatomical_analysis.py`, 21 steps).
The earlier ComBat-only (no QC/regression) pipeline is superseded and lives
in `legacy/4_anatomical_analysis/` (local only, not part of this public
release).

```bash
cd 4_anatomical_analysis && ./run_all.sh
```

## 5. Cluster anatomical MRI

Structural MRI contrasts by cluster (t-stats, Cohen's d) plus ggseg and
subcortical visualisations, including the Figure 6a composite grid, on the
current curated k-means clustering. The earlier per-cluster pipeline
(frozen k-means/Ward) is superseded and lives in
`legacy/5_cluster_anatomical_analysis/` (local only, not part of this
public release).

```bash
cd 5_cluster_anatomical_analysis/scripts
python3.11 run_cluster_anatomical_analysis.py
# Optional: restrict NT pool to NT subjects with Cluster == 'C1'
python3.11 run_cluster_anatomical_analysis.py --nt-c1

# Legacy frozen pipeline (local only, not part of this public release)
cd ../../legacy/5_cluster_anatomical_analysis/scripts
python3.11 run_cluster_anatomical_analysis.py
```

## 6. Functional MRI

Functional connectivity, whole-group autism vs NT. Two self-contained,
parallel analyses under one folder, differing in preprocessing (XCP-D
version) and in whether runs are concatenated per subject before computing
connectivity — see `6_functional_analysis/README.md`.

**`6_functional_analysis/concat/`** (current) — XCP-D v0.11. The primary
variant is `nogsr_concat` (run-concatenated, no global-signal-regression,
>6 min; n=393 autism / 327 NT) — used directly by `11_age_sex_analysis`.

```bash
cd 6_functional_analysis/concat && ./run_all.sh
```

Runs preprocessing (build connectivity → ComBat → regress → z-score) on
`nogsr_concat`, the core edge/network-block/NBS analysis suite, the
Figure 6b/3 brain-map chain, and QC (motion confound, sequence/ComBat SVM,
FD x connectivity). `run_gsr_all_modes.sh` separately sweeps the other 8
sensitivity variants (gsr, mixed, regfirst, per-cohort, anatomical-QC-pass)
— run manually, not part of `run_all.sh` (it is a long sweep, not a single
pipeline run).

**`6_functional_analysis/non_concat/`** (secondary/historical analysis) —
the original XCP-D 0.8 preprocessing, single best run per subject. Kept for
the sensitivity/legacy-comparison checks it still answers.

```bash
cd 6_functional_analysis/non_concat/preprocessing
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

An earlier, now-superseded attempt at a concatenated-runs pipeline lives in
`legacy/6-1_functional_analysis_concat/` (local only, not part of this
public release).

## 7. Cluster functional MRI

Cluster-wise connectivity (inputs, Schaefer plots, network matrices,
hyper-/hypo-connectivity counts, subcortical yabplot), including the
Figure 6b per-cluster NBS chain. Reads the same primary-regime connectivity
table as section 6 (`6_functional_analysis/concat`'s `nogsr_concat`).

```bash
cd 7_cluster_functional_analysis && ./run_all.sh
```

See `7_cluster_functional_analysis/README.md`.

## 8. EEG

Resting-state EEG, autism vs NT: alpha peak frequency and multi-band
power. Alpha peak is rebuilt from curated raw features
(`8_eeg_analysis/preprocessing/`); power spectrum has no raw counterpart on
`$IMG5` and stays on the frozen sibling-repo table, with curated
demographics grafted on. See `8_eeg_analysis/README.md`.

```bash
cd 8_eeg_analysis/scripts
python3.11 run_eeg_analysis.py
```

## 9. Cluster EEG

Same alpha-peak / band-power analyses by cluster, sharing section 8's
preprocessing. See `9_cluster_eeg_analysis/README.md`.

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
analyses (corresponds to "The Effects of Age and Sex"). Depends on
`4_anatomical_analysis/` directly (imports its config and script 01).

```bash
cd 11_age_sex_analysis && ./run_all.sh
```

## Citation

If you use this resource, please cite:

> LEAP-InovAND: A multiscale resource to explore genetics, brain imaging, and clinical data in autism.
> medRxiv preprint (2025). https://doi.org/10.1101/2025.11.24.25340858
