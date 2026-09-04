# =============================================================================
# Cluster fMRI Connectivity — Schaefer Atlas Visualisation
# =============================================================================
#
# For each autism cluster (C1, C2, C3):
#   - Reads the R-ready input CSV produced by 01_generate_cluster_fmri_inputs.py
#   - Maps connectivity column names to Schaefer 7-network labels via the atlas TSV
#   - Sums t-statistics per network region
#   - Plots a Schaefer 7×100 brain map with a diverging RdBu palette
#   - Saves per-cluster PDF and a combined three-panel PDF
#
# Adapted from eeg_mri-pipeline — self-contained, no source() calls to that repo.
# =============================================================================

# ---- library loading order matters: plyr BEFORE dplyr ----------------------
library(plyr)
library(dplyr)
library(ggplot2)
library(colorspace)
library(ggseg)
library(ggsegSchaefer)
library(tidyverse)
library(conflicted)

conflicts_prefer(dplyr::rename)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::mutate)

# =============================================================================
# PATHS  (commandArgs-based, compatible with Rscript and interactive use)
# =============================================================================
args        <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
if (length(script_path) == 0) {
  script_dir <- getwd()
} else {
  script_dir <- dirname(normalizePath(script_path))
}
section_dir <- dirname(script_dir)
output_dir  <- file.path(section_dir, "outputs", "figures")
r_input_dir <- file.path(section_dir, "outputs", "tables", "r_input_files")
setwd(output_dir)

ATLAS_FILE <- "/Users/mfleury/POSTDOC/LIBRAIRY/eeg_mri-pipeline/ressources/atlases/SCHAEFER/atlas-4S156Parcels/atlas-4S156Parcels_dseg.tsv"


# =============================================================================
# PLOT FUNCTION
# =============================================================================

# Helper: map con_Source/Target region strings to atlas label_7network,
# counting both source and target appearances per region.
count_edges_per_region <- function(regions_vec, atlas_df) {
  if (length(regions_vec) == 0) return(data.frame(region = character(), n_edges = integer()))

  sources <- sub("^con_", "", sapply(strsplit(regions_vec, "/"), `[`, 1))
  targets <- sub(".*/", "", regions_vec)

  all_labels <- c(sources, targets)

  df_labels <- data.frame(label = all_labels, stringsAsFactors = FALSE) %>%
    dplyr::left_join(atlas_df[, c("label", "label_7network")], by = "label") %>%
    dplyr::filter(!is.na(label_7network))

  df_labels %>%
    dplyr::group_by(label_7network) %>%
    dplyr::summarise(n_edges = dplyr::n(), .groups = "drop") %>%
    dplyr::rename(region = label_7network)
}

# Helper: build one brain map for edge counts.
# color_low / color_high: sequential palette; 0 values are set to NA → grey.
make_edge_map <- function(region_data, title_str, color_low, color_high) {
  # Replace 0 with NA so those regions appear grey
  region_data$n_edges[region_data$n_edges == 0] <- NA

  region_data %>%
    ggseg(
      mapping = aes(fill = n_edges),
      atlas   = schaefer7_100,
      colour  = "white",
      size    = 0.3
    ) +
    scale_fill_gradientn(
      colours  = c(color_low, color_high),
      na.value = "grey80",
      name     = "No. of edges"
    ) +
    labs(title = title_str) +
    theme_void() +
    theme(
      plot.title   = element_text(hjust = 0.5, size = 11, face = "bold"),
      legend.title = element_text(size = 9),
      legend.text  = element_text(size = 8)
    )
}


plot_cluster_schaefer <- function(cluster_id, output_dir, r_input_dir, atlas_file) {

  cat(sprintf("\n--- Cluster %s ---\n", cluster_id))

  # ---- load R-ready input (has region, mean_value=t_stat, p_fdr) -----------
  csv_path <- file.path(r_input_dir, sprintf("cluster_%s_connectivity_for_r.csv", cluster_id))
  if (!file.exists(csv_path)) {
    cat(sprintf("Warning: input file not found, skipping %s: %s\n", cluster_id, csv_path))
    return(NULL)
  }
  df <- read.csv(csv_path)
  if (nrow(df) == 0) {
    cat(sprintf("Warning: empty input file for cluster %s\n", cluster_id))
    return(NULL)
  }
  cat(sprintf("Loaded %d connectivity features\n", nrow(df)))

  # ---- load atlas -----------------------------------------------------------
  if (!file.exists(atlas_file)) {
    cat(sprintf("Warning: atlas file not found: %s\n", atlas_file))
    return(NULL)
  }
  df_atlas <- read.table(atlas_file, sep = "\t", header = TRUE)
  df_atlas_schaefer <- df_atlas[df_atlas$atlas_name == "4S156", ]

  # ---- split: FDR-sig hyper / hypo, by Cohen's d direction (Reviewer 3) ----
  df_sig   <- df[!is.na(df$p_fdr) & df$p_fdr < 0.05, ]
  dir_col  <- if ("cohens_d" %in% names(df_sig)) "cohens_d" else "mean_value"
  df_hyper <- df_sig[df_sig[[dir_col]] > 0, ]
  df_hypo  <- df_sig[df_sig[[dir_col]] < 0, ]

  cat(sprintf("  FDR-sig: %d total, %d hyper, %d hypo\n",
              nrow(df_sig), nrow(df_hyper), nrow(df_hypo)))

  # ---- count edges per Schaefer region (source + target) -------------------
  hyper_counts <- count_edges_per_region(df_hyper$region, df_atlas_schaefer)
  hypo_counts  <- count_edges_per_region(df_hypo$region,  df_atlas_schaefer)

  # ---- build plots (NULL if no data) ---------------------------------------
  p_hyper <- NULL
  p_hypo  <- NULL

  if (nrow(hyper_counts) > 0) {
    p_hyper <- make_edge_map(
      hyper_counts,
      sprintf("Cluster %s — Hyperconnectivity vs NT (FDR < 0.05)", cluster_id),
      color_low  = "#fcae91",   # light red
      color_high = "#D62728"    # dark red
    )
    out_hyper <- file.path(output_dir,
                           sprintf("cluster_%s_connectivity_schaefer_hyper.pdf", cluster_id))
    ggsave(out_hyper, plot = p_hyper, width = 6, height = 4, device = "pdf")
    cat(sprintf("  Saved: %s\n", basename(out_hyper)))
  } else {
    cat(sprintf("  No hyper edges for %s — skipping hyper map\n", cluster_id))
  }

  if (nrow(hypo_counts) > 0) {
    p_hypo <- make_edge_map(
      hypo_counts,
      sprintf("Cluster %s — Hypoconnectivity vs NT (FDR < 0.05)", cluster_id),
      color_low  = "#9ecae1",   # light blue
      color_high = "#1F77B4"    # dark blue
    )
    out_hypo <- file.path(output_dir,
                          sprintf("cluster_%s_connectivity_schaefer_hypo.pdf", cluster_id))
    ggsave(out_hypo, plot = p_hypo, width = 6, height = 4, device = "pdf")
    cat(sprintf("  Saved: %s\n", basename(out_hypo)))
  } else {
    cat(sprintf("  No hypo edges for %s — skipping hypo map\n", cluster_id))
  }

  return(list(hyper = p_hyper, hypo = p_hypo))
}


# =============================================================================
# MAIN — loop over clusters
# =============================================================================

clusters <- c("C1", "C2", "C3")
cluster_results <- list()

for (cl in clusters) {
  res <- plot_cluster_schaefer(cl, output_dir, r_input_dir, ATLAS_FILE)
  if (!is.null(res)) {
    cluster_results[[cl]] <- res
  }
}


# =============================================================================
# COMBINED FIGURE — clusters (rows) × hyper/hypo (columns)
# =============================================================================

if (length(cluster_results) >= 1 && requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)

  # Collect hyper plots (left column) and hypo plots (right column)
  hyper_plots <- lapply(cluster_results, `[[`, "hyper")
  hypo_plots  <- lapply(cluster_results, `[[`, "hypo")

  # Replace NULLs with blank placeholder
  blank <- ggplot() + theme_void()
  hyper_plots <- lapply(hyper_plots, function(p) if (is.null(p)) blank else p)
  hypo_plots  <- lapply(hypo_plots,  function(p) if (is.null(p)) blank else p)

  # Interleave: hyper C1, hypo C1, hyper C2, hypo C2, ...
  all_plots <- unlist(
    mapply(list, hyper_plots, hypo_plots, SIMPLIFY = FALSE),
    recursive = FALSE
  )

  combined <- Reduce(`+`, all_plots) +
    plot_layout(ncol = 2, guides = "collect") +
    plot_annotation(
      title = "Autism Cluster Connectivity vs NT — Hyper (left) / Hypo (right)",
      theme = theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
    )

  out_combined <- file.path(output_dir, "clusters_combined_connectivity_schaefer.pdf")
  ggsave(out_combined, plot = combined,
         width = 12, height = 4 * length(cluster_results), device = "pdf")
  cat(sprintf("\nSaved combined figure: %s\n", basename(out_combined)))
}


cat("\n", paste(rep("=", 60), collapse = ""), "\n", sep = "")
cat("Cluster Schaefer atlas visualisation complete!\n")
cat(paste(rep("=", 60), collapse = ""), "\n", sep = "")
