# =============================================================================
# Curated-pipeline wrapper — combined single-page brain figure (Autism vs NT)
# =============================================================================
# Thin wrapper around the shared _combined_brain_figure_shared.R (also used
# directly, unparameterised, by the legacy frozen pipeline): sets the
# curated-pipeline parameters (Cohen's d fill, fixed +/-0.4 scale, Autism vs
# NT) then sources the shared builder. Input/output dirs resolve to this
# tree automatically (the shared script derives them from its caller's
# location, but we set them explicitly here anyway).
# =============================================================================

args        <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
script_dir  <- if (length(script_path) == 0) getwd() else dirname(normalizePath(script_path))
section_dir <- dirname(script_dir)                       # 4_anatomical_analysis/
shared      <- file.path(script_dir, "_combined_brain_figure_shared.R")

Sys.setenv(
  ANAT_INPUT_DIR  = file.path(section_dir, "outputs", "figures", "r_input_files"),
  ANAT_OUTPUT_DIR = file.path(section_dir, "outputs", "figures"),
  ANAT_FILL       = "cohensd",
  ANAT_SCALE_MAX  = "0.4",
  ANAT_COMPARISON = "Autism vs NT"
)
source(shared)
