# 2. Genetic Analysis

Rare variant carrier analysis (deletions, LoF, missenses) across gene lists and clinical clusters for the LEAP-InovAND Nature Neuroscience paper.

## Methodology Summary

1. **Carrier annotation**: Individuals are annotated as carriers of rare coding deletions (validated CNVs absent from gnomAD), LoF SNVs (LOFTEE HC, canonical), and pathogenic missenses (AlphaMissense) across gene lists: HCNDD dominant, HCNDD X-linked recessive, SPARK/SFARI1, EAGLE definitive/strong, SynGO, ChromEpiTF, and constrained genes
2. **LOEUF constraint scores**: Each carrier is assigned the lowest gnomAD LOEUF score among their affected genes
3. **Carrier frequencies**: Proportion of carriers per population group, with 95% CI
4. **Odds ratios**: Fisher exact test vs TD controls, with FDR correction (Benjamini-Hochberg)
5. **Stratifications**: Pan-ancestry / EUR / non-EUR, both cohorts / InovAND / LEAP, Population1 (Autism with/without IDD) / PopulationS1 (Autism grouped), DEL+LoF / DEL+LoF+Miss
6. **Cluster analysis**: Same analyses stratified by clinical clusters (C1, C2, C3) with TD and IDD
7. **PGS**: Polygenic score odds ratios (top 25% vs bottom 25%) by cluster for autism, ADHD, intelligence, anxiety, MDD, PTSD, and other traits

## Scripts

| Script | Description |
|--------|-------------|
| `_config.py` | Shared configuration: paths, palettes, gene list labels |
| `01_carrier_annotation.py` | Annotate carriers (DEL, LoF, missense) across gene lists; assign LOEUF scores |
| `02_carrier_freq_or.py` | Population-level carrier frequencies and odds ratios (13 analysis blocks) |
| `03_carrier_freq_or_clusters.py` | Cluster-level carrier frequencies, odds ratios, and PGS analysis |

## Run

```bash
# From 2_genetic_analysis/
./run_all.sh

# Or individually from 2_genetic_analysis/scripts/
python 01_carrier_annotation.py
python 02_carrier_freq_or.py
python 03_carrier_freq_or_clusters.py
```

## Outputs

### Figures

| Figure | Description |
|--------|-------------|
| `PAN_constraint_freq_population1.pdf` | Carrier freq by population (Autism with/without IDD), constrained genes |
| `PAN_constraint_or_population1.pdf` | OR by population (Autism with/without IDD), constrained genes |
| `PAN_constraint_freq.pdf` | Carrier freq by population (Autism grouped), constrained genes |
| `PAN_constraint_or.pdf` | OR by population (Autism grouped), constrained genes |
| `all_ancestries_freq.pdf` | Carrier freq, all gene lists, pan-ancestry |
| `all_ancestries_or.pdf` | OR, all gene lists, pan-ancestry |
| `EUR_constraint_freq.pdf` | Carrier freq, constrained, EUR ancestry |
| `EUR_constraint_or.pdf` | OR, constrained, EUR ancestry |
| `EUR_freq.pdf` / `EUR_or.pdf` | Carrier freq/OR, all genes, EUR ancestry |
| `nonEUR_constraint_freq.pdf` / `nonEUR_constraint_or.pdf` | Non-EUR ancestry, constrained |
| `nonEUR_freq.pdf` / `nonEUR_or.pdf` | Non-EUR ancestry, all genes |
| `PAN_constraint_freq_inovand.pdf` / `PAN_constraint_or_inovand.pdf` | InovAND cohort, constrained |
| `PAN_freq_inovand.pdf` / `PAN_or_inovand.pdf` | InovAND cohort, all genes |
| `PAN_constraint_freq_leap.pdf` / `PAN_constraint_or_leap.pdf` | LEAP cohort, constrained |
| `PAN_freq_leap.pdf` / `PAN_or_leap.pdf` | LEAP cohort, all genes |
| `PAN_freq_dellofmiss_constraint.pdf` / `PAN_or_dellofmiss_constraint.pdf` | DEL+LoF+Miss, constrained |
| `PAN_freq_dellofmiss_allconstraint.pdf` / `PAN_or_dellofmiss_allconstraint.pdf` | DEL+LoF+Miss, all genes |
| `CLUSTER_TD_constraint_freq.pdf` | Carrier freq by cluster (C1-C3, TD, IDD) |
| `CLUSTER_TD_constraint_or.pdf` | OR by cluster |
| `CLUSTER_pgs_OR_top25_bottom25.pdf` | PGS OR forest plot by cluster |

### Tables

| Table | Description |
|-------|-------------|
| `carrier_annotations.tsv` | Full carrier annotation matrix (one row per individual) |
| `PAN_constraint_or_population1.csv` | OR table, pan-ancestry, constrained, Population1 |
| `PAN_constraint_or.csv` | OR table, pan-ancestry, constrained |
| `CLUSTER_TD_constraint_or.csv` | OR table, cluster-level |
| `pgs_odds_ratios_clusters.csv` | PGS OR by cluster |
| *(+ one CSV per OR figure)* | |

## Data

Scripts expect data files in the directory set by `LEAP_INOVAND_DATA` environment variable:

```bash
export LEAP_INOVAND_DATA="/path/to/data/directory"
```

Required files:
- `individuals_metrics.tsv` — clinical data
- `gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz` — gnomAD LOEUF scores
- `SV_DEL_annotation.tsv` — structural variant deletions
- `slivar_all_ensg_annotated.xlsx` — SNV/indels (LoF + missenses)
- `diag_listing.csv` — returnable causative variant flags
- `InovAND-LEAP_pgs_SBayesRC.tsv` — polygenic scores

Additional paths:
- `LEAP_GENETIC_DATA` — alternative directory for genetic files (defaults to `LEAP_INOVAND_DATA`)
- `LEAP_GENELIST_FILE` — path to GenesTrek/HGNC gene list file

Cluster assignments are loaded from `1_clustering/outputs/tables/`.

## Python Dependencies

```
pandas numpy matplotlib seaborn scikit-learn scipy statsmodels openpyxl
```
