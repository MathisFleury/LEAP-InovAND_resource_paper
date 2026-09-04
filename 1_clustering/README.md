# 1. Clustering Analysis

Clinical clustering on IQ and SRS-2 dimensions for the LEAP-InovAND Nature Neuroscience paper.

## Methodology Summary

1. **Broad clinical measures**: IQ subscales (full-scale, verbal, non-verbal), SRS-2, Vineland (composite + domains), RBS-R, SSP
2. **PCA**: Exclude missing data and relatives; standardize features. First two components explain >75% variance; IQ and SRS contribute most.
3. **Feature selection**: IQ and SRS maximize sample size relative to the broader Vineland-inclusive set.
4. **Cluster validation** (k=2–10): NbClust (26 indices) plus AIC/BIC/ICL, Silhouette, VRS, Davies-Bouldin, Pseudo-F, Gap statistic — all favor k=3.
5. **Final clustering**: k=3 on standardised IQ × SRS via **K-means** (current priority regime — see root README for the K-means-vs-Ward rationale). Ward and GMM are retained as sensitivity analyses.

## Scripts

This tree is **curated-only** — no deprecated-data reference anywhere in it.
Every frozen/paper-reproduction script lives in `../legacy/1_clustering/`.

| Script | Description |
|--------|-------------|
| `04_run_clustering.py` | K-means (+ Ward/GMM sensitivity) k=3 on IQ and SRS, curated per-cohort clinical TSVs (LEAP+INOVAND+INFOR). Writes the curated cluster label tables most downstream sections read by default. Frozen counterpart: `legacy/1_clustering/04_run_clustering_frozen.py` — a genuinely separate implementation (different multi-cohort loading logic), not a thin wrapper, so it's a full separate file rather than an arg on this one. |
| `01_pca_features.R` | PCA on all clinical measures; scree plot, loadings, variable contributions. Reads the curated cluster table by default. Frozen run: `legacy/1_clustering/01_pca_features_frozen.R` (a thin invoker for this same engine). |
| `01b_sample_size_table.R` | Companion gt table for `01`'s sample-size completeness. `[curated]` arg. |
| `05_cluster_stability.py` | Bootstrap/subsampling/noise stability checks, all 3 methods. Frozen run: `legacy/1_clustering/05_cluster_stability_frozen.py` (thin invoker, same reasoning as `01`). |
| `05b_stability_figure_merged.py` | Merges the 3 per-method stability figures into one (defaults to curated). |
| `06_method_comparison.py` / `06b_method_comparison_table.R` | K-means vs Ward vs GMM comparative metrics (internal indices, bootstrap ARI). `[curated]` arg — kept shared (not split) since it's a comparison between both regimes' data, not a deprecated-path reference of its own; its frozen-mode invocation lives in `legacy/1_clustering/run_all.sh` since it needs frozen `04`'s output. |
| `07_autism_only_clustering.py` / `07b_autism_only_table.R` | Autism-only clustering sensitivity (Reviewer #4.6). |
| `10_cluster_scatter_panels.py` | Cluster scatter panels (IQ × SRS, 4 colourings); defaults to curated, `frozen` arg for the paper-reproduction table. |
| `run_anatomical_with_other_methods.py`, `run_downstream_with_curated.py`, `run_genetics_with_other_methods.py` | Cross-section batch drivers — rerun the relevant downstream sections once per clustering method (K-means/Ward/GMM) on curated data. |

`../legacy/1_clustering/` holds the entire frozen/paper-reproduction
pipeline, fully self-contained (own `run_all.sh`, starting with
`04_run_clustering_frozen.py`): the frozen `04`/`01`/`05` counterparts,
plus every genuinely frozen-only leaf script with no downstream dependents
anywhere in the repo (`02_feature_selection_rationale.R`,
`03_cluster_validation.R`, `08_cluster_heatmap.py`, `09_cluster_jointplot.py`,
`12_reval_kselection.py`).

## Run

```bash
# From 1_clustering/ — curated only (current, priority regime)
./run_all.sh

# From legacy/1_clustering/ — the entire frozen/paper-reproduction
# pipeline, fully self-contained (writes shared output to the real
# 1_clustering/outputs/, not legacy/outputs/)
cd ../legacy/1_clustering && ./run_all.sh
```

## Data

Curated: reads the curated per-cohort clinical TSVs (see root `CLAUDE.md`).
Frozen (`legacy/1_clustering/`): requires `individuals_metrics.tsv` with
columns `ID`, `Relation_to_proposant`, `total_IQ`, `performance_IQ`,
`SRS_tscore`, and optionally `ssp_total`, `RBS-R_total`,
`vabsdscoresc_dss`, `vabsdscoresd_dss`, `vabsdscoress_dss`,
`vabsabcabc_standard` for full PCA.

Set path: `export LEAP_INOVAND_DATA="/path/to/data/dir"`
