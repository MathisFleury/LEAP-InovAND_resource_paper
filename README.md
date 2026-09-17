# LEAP-InovAND Resource


Analysis code for the paper **LEAP-InovAND: a multiscale resource to explore genetics,
brain imaging and clinical data in autism**
([doi.org/10.1101/2025.11.24.25340858](https://doi.org/10.1101/2025.11.24.25340858)).

[![Paper](https://img.shields.io/badge/Paper-Arxiv-red)](https://doi.org/10.1101/2025.11.24.25340858)

This repository is **code only** — no participant data or generated outputs
are versioned. Each section regenerates its own `outputs/` locally when
run against the appropriate data directory.

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
```

Each numbered section follows the same convention:

```
<section>/
├── run_all.sh        # entry point — runs every step in order
├── scripts/          # numbered R / Python scripts
└── outputs/          # generated figures/, tables/ (not tracked)
```

Sections without a `run_all.sh` have a `run_*.py` orchestrator in `scripts/`
instead.

## Data

All scripts read participant-level data from a directory pointed to by
environment variables. Set these once in your shell:

```bash
export LEAP_INOVAND_DATA="/path/to/data/directory"
export LEAP_GENETIC_DATA="/path/to/genetic/files"     # optional, defaults to LEAP_INOVAND_DATA
export LEAP_GENELIST_FILE="/path/to/hgnc_genelist.txt"  # optional
```

Cluster-aware sections (3, 5, 7, 9, 10) read cluster labels written by
`1_clustering/scripts/4_run_clustering.py`. Run section 1 first.

## Dependencies

- **Python 3.11** (`/usr/local/bin/python3.11`):
  `pandas numpy scipy scikit-learn statsmodels matplotlib seaborn openpyxl`
- **R**:
  `readr FactoMineR factoextra NbClust mclust cluster fpc ggplot2 ggseg ggsegSchaefer corrplot dplyr gt`

---

## Citation

If you use this resource, please cite:

> LEAP-InovAND: A multiscale resource to explore genetics, brain imaging, and clinical data in autism.
> medRxiv preprint (2025). https://doi.org/10.1101/2025.11.24.25340858
