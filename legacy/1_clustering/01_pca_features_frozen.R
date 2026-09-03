# =============================================================================
# 01 - PCA on clinical features, FROZEN / paper-reproduction
# =============================================================================
# Thin invoker for 1_clustering/scripts/01_pca_features.R's shared engine, in
# its original frozen mode: reads the DEPRECATED individuals_metrics.tsv (see
# root CLAUDE.md -- "retained only for reproducing the original paper").
# Curated is the default/only mode shown in 1_clustering/scripts/ itself; this
# wrapper is what keeps the deprecated-data reference out of that main tree.
#
# Run: cd legacy/1_clustering && Rscript 01_pca_features_frozen.R
# Output: 1_clustering/outputs/ (same place the frozen run always wrote to --
# self-contained under legacy/ isn't worth it here since nothing reads this
# script's own figures/tables, but 04_run_clustering.py's frozen output lives
# there and this keeps all frozen PCA/clustering artifacts together).
# =============================================================================
this_dir <- getwd()
engine_dir <- normalizePath(file.path(this_dir, "..", "..", "1_clustering", "scripts"))
frozen_csv <- normalizePath(
  file.path(this_dir, "..", "..", "..", "imaging2genet", "0_input", "dataframes",
            "individuals_metrics.tsv"),
  mustWork = FALSE
)

Sys.setenv(PCA_INPUT_CSV = frozen_csv, PCA_OUTPUT_SUBDIR = "")
setwd(engine_dir)
source("01_pca_features.R")
setwd(this_dir)
