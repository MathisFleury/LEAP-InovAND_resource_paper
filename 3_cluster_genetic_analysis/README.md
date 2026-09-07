# 3. Cluster-level genetic analysis

Rare-variant carrier frequencies, odds ratios, and PGS by clinical cluster,
on the **gnomAD v4** genetics data used in the paper, stratified by the
**current curated k-means clustering** (matches the
`4_anatomical_analysis`/`5_cluster_anatomical_analysis` and
`6_functional_analysis`/`7_cluster_functional_analysis` convention: this
section holds the cluster-level counterpart of `2_genetic_analysis`'s
population-level analysis).

Every other clustering-method variant explored previously
(legacy/frozen k-means-or-Ward, GMM, manual cutoffs, SPARK-LOF, and
hg19/hg38-v2 genetics) is superseded and lives in
`../legacy/3_cluster_genetic_analysis/` (local only, not part of this public
release) with its own, fully self-contained `run_all.sh`.

## Scripts

| Script | Description |
|--------|-------------|
| `_config.py` | Shared configuration: paths, palettes, gene list labels, output directories |
| `1_carrier_freq_or_clusters.py` | Cluster-level (C1/C2/C3 + NT + IDD) carrier freq/OR on the gnomAD v4 matrix. Moved here from `2_genetic_analysis`. Imports its plotting engine from `2_genetic_analysis/scripts/_carrier_freq_or_shared.py`. |
| `2_interaction_v4_inputs.py` | Assembles the gnomAD v4 input table (curated clusters + v4 carriers + v4 PGS) for the common(PGS) x rare(carrier) interaction figure. |
| `3_interaction_common_rare.R` | Common (PGS) x rare (carrier) interaction figures, by cluster/population (`ANCESTRY=PAN` or `EUR`). |
| `4_plot_iq_pgs_cluster_loeuf.R` | Cluster-coloured, LOEUF-sized IQ x PGS-intelligence panel. Moved here from `2_genetic_analysis` since it's cluster-coloured; reads the plot frame `2_genetic_analysis/scripts/2_iq_pgs_loeuf_figures.py` writes. |

## Run

```bash
cd 3_cluster_genetic_analysis && ./run_all.sh

# Legacy pipeline (local only, not part of this public release)
cd ../legacy/3_cluster_genetic_analysis && ./run_all.sh
```

## Outputs

Flat `outputs/figures/` and `outputs/tables/`, all filenames `PAN_`/`EUR_`-
prefixed for the ancestry-stratified analyses (matches
`1_carrier_freq_or_clusters.py`'s convention) — no `PAN/`/`EUR/`
subfolders.

## Legacy cluster-source glossary

The legacy tree repeats the same carrier-frequency/OR/PGS analysis across
several ways of defining "which cluster does each person belong to", from
when the project's clustering methodology was still being finalised (see
root `1_clustering/`). None of these names are used in the published
figures/tables — they're this repo's internal labels for where the cluster
column came from:

| Label | Cluster column | Source file |
|---|---|---|
| **Legacy** | `Cluster` | `1_clustering/outputs/tables/cluster_assignments.csv` (the original, frozen k-means/Ward run) |
| **GMM** | `Cluster_GMM` | `df_multi_dataset_with_clusters.csv` |
| **Manual** | `Cluster_Manual` | `df_multi_dataset_with_clusters.csv` |

Each of Legacy/GMM has an **NT-restricted** variant (`_nt_c1`) that pools
only NT subjects assigned `Cluster == 'C1'` as the reference group, instead
of all NT. Each of Legacy/GMM/Manual has an **hg38** counterpart re-running
the same analysis on GRCh38-coordinate carrier data.

**SPARK LoF** (`09`, `10`, `13`, `14`) is a separate axis: external SPARK-cohort
LoF-variant replication, cross-tabulated against GMM or Manual clusters (and,
in `13`/`14`, against LEAP-InovAND's own NT rather than SPARK's).

**gnomAD v4, all 5 strategies** (`15`, plus its curated-method wrapper
`run_curated_v4_cluster_genetic.py`): `15` reruns the paper-label cluster
strategies (legacy k-means, GMM, legacy_nt_c1, gmm_nt_c1, manual) on the v4
carrier matrix (both PLI definitions, PAN + EUR ancestry). Its k-means
branch is superseded by the current tree's `1_carrier_freq_or_clusters.py`,
whose own curated-method comparison (k-means/Ward/GMM) is handled by
`1_clustering/scripts/run_genetics_with_other_methods.py` (file-swap based)
rather than `run_curated_v4_cluster_genetic.py`'s env-var based approach.
