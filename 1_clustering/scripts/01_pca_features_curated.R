# PCA on clinical features, CURATED cohort (see 01_pca_features.R for the frozen/paper version).
# Runs 01_pca_features.R on individuals_metrics_with_clusters_curated.csv, writing to
# 1_clustering/outputs/curated/. The frozen-paper run (01_pca_features.R with no env) is unaffected.

input_csv <- file.path("..", "outputs", "curated", "tables",
                        "individuals_metrics_with_clusters_curated.csv")
if (!file.exists(input_csv)) {
  stop("Curated features not found: ", input_csv, "\nRun 04_run_clustering_curated.py first.")
}

Sys.setenv(
  PCA_INPUT_CSV = normalizePath(input_csv),
  PCA_OUTPUT_SUBDIR = "curated"
)
source("01_pca_features.R")
