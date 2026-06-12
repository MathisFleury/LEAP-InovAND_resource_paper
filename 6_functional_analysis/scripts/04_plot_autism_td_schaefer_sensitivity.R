# =============================================================================
# Schaefer atlas visualisations for the Autism-vs-NT sensitivity analyses
# (6-min cutoff and covariate-adjusted).
#
# Mirrors 02_plot_autism_td_schaefer.R but reads the *_6min and *_covadj
# _full_data.csv files written under outputs/figures/sensitivity/.
# =============================================================================

suppressPackageStartupMessages({
  library(ggseg)
  library(plyr)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(colorspace)
  library(tidyverse)
  library(conflicted)
  library(ggsegSchaefer)
  library(patchwork)
})

conflicts_prefer(dplyr::rename)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::mutate)

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
script_dir <- if (length(script_path) == 0) getwd() else dirname(normalizePath(script_path))
section_dir <- dirname(script_dir)

INPUT_DIR  <- file.path(section_dir, "outputs", "figures", "sensitivity")
OUTPUT_DIR <- file.path(section_dir, "outputs", "figures", "sensitivity")
ATLAS_FILE <- "/Users/mfleury/POSTDOC/LIBRAIRY/eeg_mri-pipeline/ressources/atlases/SCHAEFER/atlas-4S156Parcels/atlas-4S156Parcels_dseg.tsv"

if (!file.exists(ATLAS_FILE)) {
  stop("Atlas file not found: ", ATLAS_FILE)
}
df_atlas <- read.table(ATLAS_FILE, sep = "\t", header = TRUE)
df_atlas_schaefer <- df_atlas[df_atlas$atlas_name == "4S156", ]
all_schaefer_regions <- df_atlas_schaefer$label_7network[!is.na(df_atlas_schaefer$label_7network)]

load_data <- function(suffix, connectivity_type) {
  fname <- sprintf("autism_vs_td_%sconnectivity_%s_full_data.csv", connectivity_type, suffix)
  fp <- file.path(INPUT_DIR, fname)
  if (!file.exists(fp)) {
    cat("Missing:", fname, "\n"); return(NULL)
  }
  data <- read.csv(fp)
  if (nrow(data) == 0) { cat("Empty:", fname, "\n"); return(NULL) }

  data <- data %>%
    dplyr::left_join(df_atlas_schaefer[c("label", "label_7network")], by = c("source" = "label")) %>%
    dplyr::filter(!is.na(label_7network)) %>%
    dplyr::rename(region = label_7network)

  region_data <- data %>%
    dplyr::group_by(region) %>%
    dplyr::summarise(mean_value = dplyr::n(), n_connections = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(connectivity_type = connectivity_type, sensitivity = suffix)

  missing_regions <- setdiff(all_schaefer_regions, region_data$region)
  if (length(missing_regions) > 0) {
    region_data <- rbind(region_data, data.frame(
      region = missing_regions, mean_value = 0, n_connections = 0,
      connectivity_type = connectivity_type, sensitivity = suffix
    ))
  }
  region_data <- region_data[region_data$region != "n/a", ]
  cat("Loaded", nrow(data), "edges for", suffix, connectivity_type, "\n")
  return(region_data)
}

plot_one <- function(region_data, title_str, color_high, lim_max) {
  ggseg(region_data, mapping = aes(fill = mean_value), atlas = schaefer7_100, position = "stacked") +
    scale_fill_gradient(low = "#EEEEEE", high = color_high, limits = c(0, lim_max), oob = scales::squish) +
    labs(title = title_str, fill = "No. of edges") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold"))
}

# Loop over sensitivity variants
variants <- c("6min", "covadj")
variant_labels <- c("6-min cutoff (+ FD < 1 mm)", "Covariate-adjusted (FD + min. quality)")
names(variant_labels) <- variants

for (variant in variants) {
  hyper_data <- load_data(variant, "hyper")
  hypo_data  <- load_data(variant, "hypo")

  if (!is.null(hyper_data)) {
    p_hyper <- plot_one(hyper_data,
                        paste("Autism vs NT (hyper) —", variant_labels[variant]),
                        "#B30000", max(10, max(hyper_data$mean_value, na.rm = TRUE)))
    out <- file.path(OUTPUT_DIR, paste0("autism_vs_td_hyper_schaefer_", variant, ".pdf"))
    ggsave(out, plot = p_hyper, width = 10, height = 8, dpi = 300)
    cat("Saved:", basename(out), "\n")
  }
  if (!is.null(hypo_data)) {
    p_hypo <- plot_one(hypo_data,
                       paste("Autism vs NT (hypo) —", variant_labels[variant]),
                       "#003366", max(20, max(hypo_data$mean_value, na.rm = TRUE)))
    out <- file.path(OUTPUT_DIR, paste0("autism_vs_td_hypo_schaefer_", variant, ".pdf"))
    ggsave(out, plot = p_hypo, width = 10, height = 8, dpi = 300)
    cat("Saved:", basename(out), "\n")
  }
  if (!is.null(hyper_data) && !is.null(hypo_data)) {
    p_combined <- p_hyper + p_hypo +
      plot_layout(guides = "collect") +
      plot_annotation(title = paste("Autism vs NT —", variant_labels[variant]),
                      theme = theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold")))
    out <- file.path(OUTPUT_DIR, paste0("autism_vs_td_side_by_side_schaefer_", variant, ".pdf"))
    ggsave(out, plot = p_combined, width = 20, height = 8, dpi = 300)
    cat("Saved:", basename(out), "\n")
  }
}

cat("\nSensitivity Schaefer visualisations complete.\n")
