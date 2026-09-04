# 2. Genetic Analysis

Rare variant carrier analysis (deletions, LoF, missenses) across gene lists and clinical clusters for the LEAP-InovAND Nature Neuroscience paper.

## Methodology Summary

1. **Carrier annotation**: Individuals are annotated as carriers of rare coding deletions (validated CNVs absent from gnomAD), LoF SNVs (LOFTEE HC, canonical), and pathogenic missenses (AlphaMissense) across gene lists: HCNDD dominant, HCNDD X-linked recessive, SPARK/SFARI1, EAGLE definitive/strong, SynGO, ChromEpiTF, and constrained genes
2. **LOEUF constraint scores**: Each carrier is assigned the lowest gnomAD LOEUF score among their affected genes
3. **Carrier frequencies**: Proportion of carriers per population group, with 95% CI
4. **Odds ratios**: Fisher exact test vs TD controls, with FDR correction (Benjamini-Hochberg)
5. **Stratifications**: Pan-ancestry / EUR / non-EUR, both cohorts / InovAND / LEAP, Population1 (Autism with/without IDD) / PopulationS1 (Autism grouped), DEL+LoF / DEL+LoF+Miss
6. **PGS**: Polygenic score odds ratios (top 25% vs bottom 25%) for autism, ADHD, intelligence, anxiety, MDD, PTSD, and other traits

This tree is **population-level only**. Cluster-level carrier freq/OR
analysis (C1, C2, C3 with TD and IDD) lives in `../3_cluster_genetic_analysis/`
instead, matching the `4_anatomical_analysis`/`5_cluster_anatomical_analysis`
and `6_functional_analysis`/`7_cluster_functional_analysis` convention.

This tree runs the **gnomAD v4** genetics data used in the paper. The
earlier hg19 and hg38/gnomAD-v2 pipeline is superseded and lives in
`../legacy/2_genetic_analysis/` (local only, not part of this public
release) with its own, fully self-contained `run_all.sh`.

## Scripts

| Script | Description |
|--------|-------------|
| `_config.py` | Shared configuration: paths, palettes, gene list labels |
| `_carrier_freq_or_shared.py` | Shared plotting/OR engine (`plot_carrier_frequencies`, `compute_and_plot_odds_ratios`). Dynamically imported by `1` below, by `3_cluster_genetic_analysis/1_carrier_freq_or_clusters.py`, *and* by the legacy hg19/hg38 scripts, so it stays in this tree rather than being duplicated. Its own standalone run (hg19 population-level frequencies) is invoked only from `legacy/2_genetic_analysis/run_all.sh`. |
| `1_carrier_freq_or.py` | Population-level carrier frequencies and odds ratios, gnomAD v4 (both constraint definitions, PAN + EUR) |
| `2_iq_pgs_loeuf_figures.py` | IQ x PGS-intelligence plot frame + IQ ~ -log10(LOEUF) beta coefficients (population-level). Writes `iq_pgs_plot_frame.csv`, read by the cluster-coloured panel below. |

The cluster-coloured version of the IQ x PGS-intelligence panel
(`4_plot_iq_pgs_cluster_loeuf.R`) lives in `../3_cluster_genetic_analysis/`
instead, matching the population-vs-cluster section split.

## Run

```bash
# From 2_genetic_analysis/
./run_all.sh

# Or individually from 2_genetic_analysis/scripts/
python3.11 1_carrier_freq_or.py
python3.11 2_iq_pgs_loeuf_figures.py

# Cluster-level genetic analysis (also runs the cluster-coloured IQ x PGS x LOEUF panel)
cd ../3_cluster_genetic_analysis && ./run_all.sh

# Legacy hg19/hg38-v2 pipeline (local only, not part of this public release)
cd ../legacy/2_genetic_analysis && ./run_all.sh
```

## Outputs

### Figures (`outputs/figures/`)

| Figure | Description |
|--------|-------------|
| `PAN_*_freq.pdf` / `PAN_*_or.pdf` | Population-level carrier freq/OR, PAN ancestry |
| `EUR_*_freq.pdf` / `EUR_*_or.pdf` | Population-level carrier freq/OR, EUR ancestry |
| `beta_coefficients_iq_loeuf.pdf` | IQ ~ -log10(LOEUF) beta coefficients by gene list |

### Tables (`outputs/tables/`)

| Table | Description |
|-------|-------------|
| `*_or.csv` | OR table per stratification/constraint definition |
| `iq_pgs_plot_frame.csv` | Plot-ready IQ x PGS-intelligence x LOEUF frame — read by `../3_cluster_genetic_analysis/scripts/4_plot_iq_pgs_cluster_loeuf.R` |

## Data

Scripts expect data files in the directory set by `LEAP_INOVAND_DATA` environment variable:

```bash
export LEAP_INOVAND_DATA="/path/to/data/directory"
```

The gnomAD v4 carrier matrix, LOEUF scores, and PGS files are read from
`$IMG5` (the mounted imaging/genetics volume) via `config.py` — see that
file for exact paths.

Additional paths:
- `LEAP_GENETIC_DATA` — alternative directory for genetic files (defaults to `LEAP_INOVAND_DATA`)
- `LEAP_GENELIST_FILE` — path to GenesTrek/HGNC gene list file

Cluster assignments are loaded from `1_clustering/outputs/tables/`.

## Python Dependencies

```
pandas numpy matplotlib seaborn scikit-learn scipy statsmodels openpyxl
```
