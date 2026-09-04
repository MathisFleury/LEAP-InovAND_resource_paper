# =============================================================================
# Curated-pipeline wrapper — combined single-page cluster brain figure (clusters vs TD)
# =============================================================================
# Thin wrapper around ../../scripts/04_combined_cluster_brain_figure.R. Input /
# output dirs resolve to this curated/ tree automatically (the shared script
# derives them from this caller's location), so the per-cluster figures land in
# curated/outputs/figures and are renamed per method by the curated runner.
# The curated pipeline reports Cohen's d, but both stats are emitted.
# =============================================================================

args        <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
script_dir  <- if (length(script_path) == 0) getwd() else dirname(normalizePath(script_path))
section_dir <- dirname(script_dir)                       # .../curated
shared      <- normalizePath(file.path(section_dir, "..", "scripts",
                                       "04_combined_cluster_brain_figure.R"))
source(shared)
