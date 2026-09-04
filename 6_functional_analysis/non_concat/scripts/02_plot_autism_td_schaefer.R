# =============================================================================
# fMRI Autism vs TD Schaefer Atlas Visualization
# =============================================================================
#
# Adapted from eeg_mri-pipeline/analysis/figures_papers/fmri_autism_td_analysis/plot_autism_td_schaefer.R
# =============================================================================

library(ggseg)
library(plyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(colorspace)
library(tidyverse)
library(conflicted)
library(ggsegSchaefer)

conflicts_prefer(dplyr::rename)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::mutate)

# Paths
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
if (length(script_path) == 0) {
  script_dir <- getwd()
} else {
  script_dir <- dirname(normalizePath(script_path))
}
section_dir <- dirname(script_dir)
# Honor AUTISM_TD_OUTPUT_DIR so legacy/revised runs land in separate folders.
# Default: the refactored-pipeline output folder, mirroring scripts/01.
output_dir <- Sys.getenv(
  "AUTISM_TD_OUTPUT_DIR",
  unset = file.path(section_dir, "outputs", "figures", "revised")
)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
setwd(output_dir)

ATLAS_FILE <- "/Users/mfleury/POSTDOC/LIBRAIRY/eeg_mri-pipeline/ressources/atlases/SCHAEFER/atlas-4S156Parcels/atlas-4S156Parcels_dseg.tsv"

load_autism_td_connectivity_data <- function(connectivity_type) {
  if (connectivity_type == "all") {
    file_path <- "autism_vs_td_hyperconnectivity_full_data.csv"
  } else {
    file_path <- paste0("autism_vs_td_", connectivity_type, "connectivity_full_data.csv")
  }
  if (!file.exists(file_path)) {
    cat("Warning: File not found:", file_path, "\n")
    return(NULL)
  }
  data <- read.csv(file_path)
  if (nrow(data) == 0) {
    cat("Warning: No data found for connectivity type:", connectivity_type, "\n")
    return(NULL)
  }
  if (connectivity_type == "all") {
    hyper_file <- "autism_vs_td_hyperconnectivity_full_data.csv"
    hypo_file  <- "autism_vs_td_hypoconnectivity_full_data.csv"
    if (file.exists(hyper_file) && file.exists(hypo_file)) {
      data <- rbind(read.csv(hyper_file), read.csv(hypo_file))
    } else {
      cat("Warning: Could not find both hyper and hypo data files\n")
      return(NULL)
    }
  } else {
    if (connectivity_type == "hyper") data <- data[data$t_stat > 0, ]
    else if (connectivity_type == "hypo") data <- data[data$t_stat < 0, ]
  }
  if (!file.exists(ATLAS_FILE)) {
    cat("Warning: Atlas file not found:", ATLAS_FILE, "\n")
    return(NULL)
  }
  df_atlas <- read.table(ATLAS_FILE, sep = "\t", header = TRUE)
  df_atlas_schaefer <- df_atlas[df_atlas$atlas_name == "4S156", ]

  data <- data %>%
    dplyr::left_join(df_atlas_schaefer[c('label', 'label_7network')], by = c('source' = 'label')) %>%
    dplyr::filter(!is.na(label_7network)) %>%
    dplyr::select(-source) %>%
    dplyr::rename(region = label_7network) %>%
    dplyr::mutate(comparison = "Autism vs TD")

  region_data <- data %>%
    dplyr::group_by(region) %>%
    dplyr::summarise(mean_value = n(), n_connections = n(), .groups = 'drop') %>%
    dplyr::mutate(comparison = "Autism vs TD", connectivity_type = connectivity_type, cluster_label = "Autism vs TD")

  all_schaefer_regions <- df_atlas_schaefer$label_7network[!is.na(df_atlas_schaefer$label_7network)]
  missing_regions <- setdiff(all_schaefer_regions, region_data$region)
  if (length(missing_regions) > 0) {
    region_data <- rbind(region_data, data.frame(
      region = missing_regions, mean_value = 0, n_connections = 0,
      comparison = "Autism vs TD", connectivity_type = connectivity_type, cluster_label = "Autism vs TD"
    ))
  }
  cat("Loaded", nrow(data), "significant edges for Autism vs TD\n")
  return(region_data)
}

create_autism_td_plots <- function(connectivity_type = "all") {
  cat("Creating Autism vs TD Schaefer atlas plots for", connectivity_type, "connectivity...\n")
  data <- load_autism_td_connectivity_data(connectivity_type)
  output_file <- paste0("autism_vs_td_", connectivity_type, "_schaefer.pdf")
  if (is.null(data) || nrow(data) == 0) {
    cat("No data available for", connectivity_type, "connectivity — writing empty-atlas PDF for parity.\n")
    df_atlas <- read.table(ATLAS_FILE, sep = "\t", header = TRUE)
    df_atlas_schaefer <- df_atlas[df_atlas$atlas_name == "4S156", ]
    all_regs <- unique(df_atlas_schaefer$label_7network[!is.na(df_atlas_schaefer$label_7network)])
    empty <- data.frame(region = all_regs, mean_value = 0)
    empty <- empty[empty$region != "n/a", ]
    p <- empty %>%
      ggseg(mapping = aes(fill = mean_value), atlas = schaefer7_100, position = "stacked") +
      scale_fill_gradient(low = "#EEEEEE", high = "#EEEEEE", limits = c(0, 1)) +
      labs(title = paste0("Autism vs TD Connectivity (", connectivity_type, ") — no significant edges"),
           fill = "No. of edges") +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
    ggsave(output_file, plot = p, width = 10, height = 8, dpi = 300)
    cat("Saved:", output_file, "(empty)\n")
    return(empty)
  }
  data <- data[data$region != "n/a", ]
  if (connectivity_type == "hyper") {
    color_scale <- scale_fill_gradient(low = "#FFDDDD", high = "#B30000", limits = c(0, 10), oob = scales::squish)
  } else if (connectivity_type == "hypo") {
    color_scale <- scale_fill_gradient(low = "#DDEEFF", high = "#003366", limits = c(0, 30), oob = scales::squish)
  } else {
    color_scale <- scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0, rev = TRUE)
  }
  p <- data %>%
    ggseg(mapping = aes(fill = mean_value), atlas = schaefer7_100, position = "stacked") +
    color_scale +
    labs(title = paste("Autism vs TD Connectivity (", connectivity_type, ")"), fill = "No. of edges") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  output_file <- paste0("autism_vs_td_", connectivity_type, "_schaefer.pdf")
  ggsave(output_file, plot = p, width = 10, height = 8, dpi = 300)
  cat("Saved:", output_file, "\n")
  return(data)
}

# Run for all connectivity types
connectivity_types <- c("all", "hyper", "hypo")
all_autism_td_data <- list()
for (conn_type in connectivity_types) {
  cat("\n--- Processing", conn_type, "connectivity ---\n")
  autism_td_data <- create_autism_td_plots(conn_type)
  if (!is.null(autism_td_data)) {
    all_autism_td_data[[conn_type]] <- autism_td_data
  }
}

# Side-by-side combined plot
if ("hyper" %in% names(all_autism_td_data) && "hypo" %in% names(all_autism_td_data)) {
  library(patchwork)
  hyper_plot <- all_autism_td_data[["hyper"]] %>%
    ggseg(mapping = aes(fill = mean_value), atlas = schaefer7_100, position = "stacked") +
    scale_fill_gradient(low = "#FFDDDD", high = "#B30000", limits = c(0, 30), oob = scales::squish) +
    labs(fill = "No. of edges", title = "Autism vs TD - Hyperconnectivity") + theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  hypo_plot <- all_autism_td_data[["hypo"]] %>%
    ggseg(mapping = aes(fill = mean_value), atlas = schaefer7_100, position = "stacked") +
    scale_fill_gradient(low = "#DDEEFF", high = "#003366", limits = c(0, 10), oob = scales::squish) +
    labs(fill = "No. of edges", title = "Autism vs TD - Hypoconnectivity") + theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  combined_plot <- hyper_plot + hypo_plot +
    plot_layout(guides = "collect") +
    plot_annotation(title = "Autism vs TD - Connectivity Differences",
                    theme = theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")))
  ggsave("autism_vs_td_side_by_side_connectivity_schaefer.pdf", plot = combined_plot, width = 20, height = 8, dpi = 300)
  cat("Saved: autism_vs_td_side_by_side_connectivity_schaefer.pdf\n")
}

cat("\n", paste(rep("=", 60), collapse=""), "\n")
cat("Autism vs TD Schaefer atlas visualization complete!\n")
cat(paste(rep("=", 60), collapse=""), "\n")
