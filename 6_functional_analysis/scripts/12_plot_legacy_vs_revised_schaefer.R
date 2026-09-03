# =============================================================================
# Side-by-side Schaefer atlas comparison: legacy vs revised preprocessing
# Reads per-edge `_full_data.csv` files from
#   outputs/figures/legacy/   (preprocessing_cohort.py output)
#   outputs/figures/revised/  (6_functional_analysis/preprocessing pipeline)
# and produces:
#   outputs/figures/legacy_vs_revised_schaefer_all.pdf
#   outputs/figures/legacy_vs_revised_schaefer_hyper.pdf
#   outputs/figures/legacy_vs_revised_schaefer_hypo.pdf
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
library(patchwork)

conflicts_prefer(dplyr::rename)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::mutate)
conflicts_prefer(dplyr::summarise)
conflicts_prefer(dplyr::summarize)

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
if (length(script_path) == 0) {
  script_dir <- getwd()
} else {
  script_dir <- dirname(normalizePath(script_path))
}
section_dir <- dirname(script_dir)
fig_root    <- file.path(section_dir, "outputs", "figures")

ATLAS_FILE <- file.path(
  "/Users/mfleury/POSTDOC/LIBRAIRY/eeg_mri-pipeline/ressources/atlases",
  "SCHAEFER/atlas-4S156Parcels/atlas-4S156Parcels_dseg.tsv"
)
df_atlas <- read.table(ATLAS_FILE, sep = "\t", header = TRUE)
df_atlas <- df_atlas[df_atlas$atlas_name == "4S156", ]
all_regions <- df_atlas$label_7network[!is.na(df_atlas$label_7network)]

load_region_counts <- function(variant_dir, conn_type) {
  if (conn_type == "all") {
    hyper_file <- file.path(variant_dir, "autism_vs_td_hyperconnectivity_full_data.csv")
    hypo_file  <- file.path(variant_dir, "autism_vs_td_hypoconnectivity_full_data.csv")
    if (!file.exists(hyper_file) || !file.exists(hypo_file)) return(NULL)
    data <- rbind(read.csv(hyper_file), read.csv(hypo_file))
  } else {
    fp <- file.path(
      variant_dir,
      paste0("autism_vs_td_", conn_type, "connectivity_full_data.csv")
    )
    if (!file.exists(fp)) return(NULL)
    data <- read.csv(fp)
  }
  if (nrow(data) == 0) return(NULL)

  data <- data %>%
    dplyr::left_join(
      df_atlas[c("label", "label_7network")],
      by = c("source" = "label")
    ) %>%
    dplyr::filter(!is.na(label_7network)) %>%
    dplyr::rename(region = label_7network)

  region_data <- data %>%
    dplyr::group_by(region) %>%
    dplyr::summarise(mean_value = dplyr::n(), .groups = "drop")

  missing <- setdiff(all_regions, region_data$region)
  if (length(missing) > 0) {
    region_data <- rbind(
      region_data,
      data.frame(region = missing, mean_value = 0)
    )
  }
  region_data <- region_data[region_data$region != "n/a", ]
  return(region_data)
}

panel_palette <- function(conn_type, ceiling_value) {
  if (conn_type == "hyper") {
    scale_fill_gradient(
      low = "#FFDDDD", high = "#B30000",
      limits = c(0, ceiling_value), oob = scales::squish
    )
  } else if (conn_type == "hypo") {
    scale_fill_gradient(
      low = "#DDEEFF", high = "#003366",
      limits = c(0, ceiling_value), oob = scales::squish
    )
  } else {
    scale_fill_continuous_divergingx(
      palette = "RdBu", mid = 0, rev = TRUE,
      limits = c(0, ceiling_value), oob = scales::squish
    )
  }
}

make_panel <- function(region_data, conn_type, ceiling_value, title) {
  if (is.null(region_data) || all(region_data$mean_value == 0)) {
    return(
      ggplot() +
        annotate("text", x = 0, y = 0,
                 label = paste0(title, "\n(no significant edges)"),
                 size = 5, lineheight = 1.2) +
        theme_void()
    )
  }
  region_data %>%
    ggseg(mapping = aes(fill = mean_value),
          atlas = schaefer7_100, position = "stacked") +
    panel_palette(conn_type, ceiling_value) +
    labs(title = title, fill = "No. of edges") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"))
}

build_comparison <- function(conn_type, out_path) {
  cat(sprintf("\n--- %s ---\n", conn_type))
  legacy_dir  <- file.path(fig_root, "legacy")
  revised_dir <- file.path(fig_root, "revised")

  legacy  <- load_region_counts(legacy_dir,  conn_type)
  revised <- load_region_counts(revised_dir, conn_type)

  ceiling_value <- max(
    if (!is.null(legacy))  max(legacy$mean_value)  else 0,
    if (!is.null(revised)) max(revised$mean_value) else 0,
    1
  )

  panel_legacy  <- make_panel(legacy,  conn_type, ceiling_value,
                              paste0("Legacy — ", conn_type))
  panel_revised <- make_panel(revised, conn_type, ceiling_value,
                              paste0("Revised — ", conn_type))

  combined <- panel_legacy + panel_revised +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = paste0("Autism vs NT (Schaefer 7-network) — ", conn_type),
      subtitle = "Left: legacy pipeline   |   Right: revised pipeline",
      theme = theme(
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10)
      )
    )

  ggsave(out_path, plot = combined, width = 20, height = 8, dpi = 300)
  cat(sprintf("  wrote  %s\n", basename(out_path)))
}

for (conn_type in c("all", "hyper", "hypo")) {
  build_comparison(
    conn_type,
    file.path(fig_root, paste0("legacy_vs_revised_schaefer_", conn_type, ".pdf"))
  )
}

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("Legacy vs revised Schaefer comparison complete\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
