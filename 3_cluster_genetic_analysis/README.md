# 3. Cluster-level genetic analysis

Rare-variant carrier frequencies, odds ratios, and PGS by clinical cluster.
Run with `./run_all.sh` (18 scripts, all wired in, in order).

## Cluster-source glossary

The same carrier-frequency/OR/PGS analysis is repeated across four different
ways of defining "which cluster does each person belong to", because the
project's clustering methodology evolved during review (see root
`1_clustering/`). None of these names appear in the manuscript — they're
this repo's internal labels for where the cluster column came from:

| Label | Cluster column | Source file |
|---|---|---|
| **Legacy** | `Cluster` | `1_clustering/outputs/tables/cluster_assignments.csv` (the original, frozen k-means/Ward run) |
| **GMM** | `Cluster_GMM` | `df_multi_dataset_with_clusters.csv` |
| **Manual** | `Cluster_Manual` | `df_multi_dataset_with_clusters.csv` |
| **Curated** (15–18) | k-means / Ward / GMM, curated cohort | `1_clustering/outputs/curated/tables/*_curated*.csv` — the current priority regime |

Each of Legacy/GMM has an **NT-restricted** variant (`_nt_c1`) that pools
only NT subjects assigned `Cluster == 'C1'` as the reference group, instead
of all NT. Each of Legacy/GMM/Manual has an **hg38** counterpart re-running
the same analysis on GRCh38-coordinate carrier data.

**SPARK LoF** (`09`, `10`, `13`, `14`) is a separate axis: external SPARK-cohort
LoF-variant replication, cross-tabulated against GMM or Manual clusters (and,
in `13`/`14`, against LEAP-InovAND's own NT rather than SPARK's).

**gnomAD v4** (`15`–`18`) is the current priority regime: `15` reruns the
paper-label cluster strategies on the v4 carrier matrix (both PLI
definitions, PAN + EUR ancestry); `16` (`run_curated_v4_cluster_genetic.py`)
reruns it once per curated clustering method (k-means/Ward/GMM); `17`/`18`
assemble the common-variant (PGS) x rare-variant (carrier) interaction
figures by cluster/population.

## Run

```bash
cd 3_cluster_genetic_analysis && ./run_all.sh
```
