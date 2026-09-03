# 1. Clustering Analysis

Clinical clustering on IQ and SRS-2 dimensions for the LEAP-InovAND Nature Neuroscience paper.

## Methodology Summary

1. **Broad clinical measures**: IQ subscales (full-scale, verbal, non-verbal), SRS-2, Vineland (composite + domains), RBS-R, SSP
2. **PCA**: Exclude missing data and relatives; standardize features. First two components explain >75% variance; IQ and SRS contribute most.
3. **Feature selection**: IQ and SRS maximize sample size relative to the broader Vineland-inclusive set.
4. **Cluster validation** (k=2–10): NbClust (26 indices) plus AIC/BIC/ICL, Silhouette, VRS, Davies-Bouldin, Pseudo-F, Gap statistic — all favor k=3.
5. **Final clustering**: k=3 on standardised IQ × SRS via **K-means** (current priority regime — see root README for the K-means-vs-Ward rationale). Ward and GMM are retained as sensitivity analyses.

## Scripts

| Script | Description |
|--------|-------------|
| `01_pca_features.R` | PCA on all clinical measures; scree plot, loadings, variable contributions. **Curated-only** — reads the curated cluster table by default, no deprecated-data fallback in this file. For the frozen/paper-reproduction run, see `legacy/1_clustering/01_pca_features_frozen.R`. |
| `01b_sample_size_table.R` | Companion gt table for `01`'s sample-size completeness. `[curated]` arg. |
| `04_run_clustering.py` | K-means (+ Ward/GMM sensitivity) k=3 on IQ and SRS, frozen data. Its output (`cluster_assignments.csv`, `individuals_metrics_with_clusters.csv`) is still read by ~20 files across 8 other sections in their frozen/default mode, so — unlike `01`/`05` — this stays a genuinely separate file rather than folding into the curated version. |
| `04_run_clustering_curated.py` | Same clustering, curated per-cohort clinical TSVs (LEAP+INOVAND+INFOR) — a genuinely separate implementation (different multi-cohort loading logic), not a thin wrapper. Writes the curated cluster label tables most downstream sections read by default. |
| `05_cluster_stability.py` | Bootstrap/subsampling/noise stability checks, all 3 methods. **Curated-only**, same reasoning as `01`. Frozen run: `legacy/1_clustering/05_cluster_stability_frozen.py`. |
| `05b_stability_figure_merged.py` | Merges the 3 per-method stability figures into one (defaults to curated). |
| `06_method_comparison.py` / `06b_method_comparison_table.R` | K-means vs Ward vs GMM comparative metrics (internal indices, bootstrap ARI). `[curated]` arg — kept shared (not split) since it never touches the deprecated data path directly for either mode's core logic beyond the frozen sample it's given. |
| `07_autism_only_clustering.py` / `07b_autism_only_table.R` | Autism-only clustering sensitivity (Reviewer #4.6). |
| `10_cluster_scatter_panels.py` | Cluster scatter panels (IQ × SRS, 4 colourings); defaults to curated, `frozen` arg for the paper-reproduction table. |
| `run_anatomical_with_other_methods.py`, `run_downstream_with_curated.py`, `run_genetics_with_other_methods.py` | Cross-section batch drivers — rerun the relevant downstream sections once per clustering method (K-means/Ward/GMM) on curated data. |

`../legacy/1_clustering/` holds every frozen-only script that would otherwise
put a deprecated-data reference in this tree, or has zero downstream
dependents anywhere in the repo: `01_pca_features_frozen.R`,
`05_cluster_stability_frozen.py` (thin invokers for the shared engines
above), plus `02_feature_selection_rationale.R`, `03_cluster_validation.R`,
`08_cluster_heatmap.py`, `09_cluster_jointplot.py`, `12_reval_kselection.py`
(genuine leaves) — see that folder's own `run_all.sh`.

## Run

```bash
# From 1_clustering/ — curated (default) + frozen blocks (04, 06 only)
./run_all.sh

# From legacy/1_clustering/ — every frozen-only script (run after
# 1_clustering/run_all.sh's frozen block, which 08/09 read the output of)
cd ../legacy/1_clustering && ./run_all.sh
```

## Data

Curated: reads the curated per-cohort clinical TSVs (see root `CLAUDE.md`).
Frozen: requires `individuals_metrics.tsv` with columns `ID`,
`Relation_to_proposant`, `total_IQ`, `performance_IQ`, `SRS_tscore`, and
optionally `ssp_total`, `RBS-R_total`, `vabsdscoresc_dss`,
`vabsdscoresd_dss`, `vabsdscoress_dss`, `vabsabcabc_standard` for full PCA.

Set path: `export LEAP_INOVAND_DATA="/path/to/data/dir"`
