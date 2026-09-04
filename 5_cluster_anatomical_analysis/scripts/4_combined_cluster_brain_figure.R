# =============================================================================
# Curated-pipeline wrapper — combined single-page cluster brain figure (clusters vs TD)
# =============================================================================
# Thin wrapper around the shared _combined_cluster_brain_figure_shared.R (also
# used directly, unparameterised, by the legacy frozen pipeline). Input/output
# dirs resolve to this tree automatically (the shared script derives them
# from this caller's location), so the per-cluster figures land in
# outputs/figures and are renamed per method by the curated runner.
# The curated pipeline reports Cohen's d, but both stats are emitted.
# =============================================================================

args        <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
script_dir  <- if (length(script_path) == 0) getwd() else dirname(normalizePath(script_path))
section_dir <- dirname(script_dir)                       # 5_cluster_anatomical_analysis/
shared      <- file.path(script_dir, "_combined_cluster_brain_figure_shared.R")
source(shared)
