# =============================================================================
# Cluster MRI Brain Visualizations (REVISION — Reviewer 3)
# =============================================================================
#
# v2 of ../../scripts/02_plot_cluster_mri_brain.R, re-rendered on the
# QC+ComBat+age/sex/eTIV-regressed FreeSurfer dataset produced by
# 01_generate_cluster_mri_inputs_v2.py (see that script for the new dataset
# path and LEAP_W1 / INOVAND_T1 wave selection).
#
# Fill is Cohen's d (Reviewer 3); FDR-significance shown by black outlines.
# Reads input from: ../outputs/figures/r_input_files/
# Writes output to: ../outputs/figures/
# =============================================================================

library(plyr)      # must be loaded BEFORE dplyr
library(dplyr)
library(ggplot2)
library(colorspace)
library(ggseg)
library(scales)    # for scales::squish (clip out-of-range values to scale ends)

# =============================================================================
# PATH SETUP
# =============================================================================
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
if (length(script_path) == 0) {
  script_dir <- getwd()
} else {
  script_dir <- dirname(normalizePath(script_path))
}
section_dir <- dirname(script_dir)
input_dir   <- file.path(section_dir, "outputs", "figures", "r_input_files")
output_dir  <- file.path(section_dir, "outputs", "figures")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

#' Convert region labels that use lh_/rh_ prefix to _left/_right suffix form
#' (ggseg dk atlas expects lowercase labels like "bankssts_left").
transform_region_column <- function(df) {
  df <- df %>%
    plyr::mutate(label = case_when(
      endsWith(label, "_left")  ~ paste0("lh_", sub("_left$",  "", label)),
      endsWith(label, "_right") ~ paste0("rh_", sub("_right$", "", label)),
      TRUE ~ label
    ))
  return(df)
}

#' Rename aseg subcortical labels from lowercase ggseg defaults to FreeSurfer
#' naming convention expected by the data files.
rename_subcortical_labels <- function(df) {
  df <- df %>%
    dplyr::mutate(label = case_when(
      label == "left-thalamus"                 ~ "Left-Thalamus-Proper",
      label == "right-thalamus"                ~ "Right-Thalamus-Proper",
      label == "left-caudate"                  ~ "Left-Caudate",
      label == "right-caudate"                 ~ "Right-Caudate",
      label == "left-putamen"                  ~ "Left-Putamen",
      label == "right-putamen"                 ~ "Right-Putamen",
      label == "left-pallidum"                 ~ "Left-Pallidum",
      label == "right-pallidum"                ~ "Right-Pallidum",
      label == "left-hippocampus"              ~ "Left-Hippocampus",
      label == "right-hippocampus"             ~ "Right-Hippocampus",
      label == "left-amygdala"                 ~ "Left-Amygdala",
      label == "right-amygdala"                ~ "Right-Amygdala",
      label == "left-accumbens"                ~ "Left-Accumbens",
      label == "right-accumbens"               ~ "Right-Accumbens",
      label == "left-ventraldc"                ~ "Left-VentralDC",
      label == "right-ventraldc"               ~ "Right-VentralDC",
      label == "left-choroid-plexus"           ~ "Left-Choroid-Plexus",
      label == "right-choroid-plexus"          ~ "Right-Choroid-Plexus",
      label == "cc_posterior"                  ~ "CC_Posterior",
      label == "cc_mid_posterior"              ~ "CC_Mid_Posterior",
      label == "cc_central"                    ~ "CC_Central",
      label == "cc_mid_anterior"               ~ "CC_Mid_Anterior",
      label == "cc_anterior"                   ~ "CC_Anterior",
      label == "right-cerebellum-white-matter" ~ "Right-Cerebellum-White-Matter",
      label == "left-cerebellum-white-matter"  ~ "Left-Cerebellum-White-Matter",
      label == "right-cerebellum-cortex"       ~ "Right-Cerebellum-Cortex",
      label == "left-cerebellum-cortex"        ~ "Left-Cerebellum-Cortex",
      label == "3rd-ventricle"                 ~ "x3rd-ventricle",
      label == "4th-ventricle"                 ~ "x4th-ventricle",
      label == "brain-stem"                    ~ "Brain-Stem",
      label == "left-lateral-ventricle"        ~ "Left-Lateral-Ventricle",
      label == "right-lateral-ventricle"       ~ "Right-Lateral-Ventricle",
      TRUE ~ label
    ))
  return(df)
}

#' Load a cluster-vs-TD CSV and return a data frame, or NULL if missing.
#'
#' The Python generator (01_generate_cluster_mri_inputs.py) writes
#'   label, t_stat, p_val, cohens_d, cohens_d_lo, cohens_d_hi, n_a, n_b, p_fdr
#' but for backward compatibility we also accept the legacy 5-column form.
load_cluster_data <- function(cluster, metric, atlas_type) {
  fname <- paste0("t_stat_cluster_", cluster, "_", metric, "_", atlas_type,
                  "_mri_cluster_vs_td.csv")
  fpath <- file.path(input_dir, fname)
  if (!file.exists(fpath)) {
    warning(paste("File not found:", fpath))
    return(NULL)
  }
  data <- read.csv(fpath, stringsAsFactors = FALSE)
  required <- c("label", "t_stat", "p_val", "cohens_d", "p_fdr")
  if (!all(required %in% names(data))) {
    data <- read.csv(fpath, header = FALSE,
                     col.names = required,
                     stringsAsFactors = FALSE)
  }
  # Optional CI columns — populate with NA if absent
  for (cc in c("cohens_d_lo", "cohens_d_hi")) {
    if (!cc %in% names(data)) data[[cc]] <- NA_real_
  }
  data <- data %>% dplyr::mutate(fdr_corrected = p_fdr, significant = p_fdr < 0.05)
  return(data)
}

# =============================================================================
# PLOTTING FUNCTIONS
# =============================================================================

#' Create and save a pair of brain plots (plain + FDR-outlined) for one
#' cluster/metric combination using the dk cortical atlas.
#'
#' Per Reviewer 3 (revision): the fill is now **Cohen's d** (not the t-stat),
#' since differential cluster sample sizes inflate statistical power
#' differentially. FDR significance is shown by black region outlines.
plot_cortical_cluster <- function(data, cluster, metric, atlas_type, output_dir,
                                  d_max_override = NULL) {
  data <- transform_region_column(data)
  data <- data %>% plyr::mutate(label = tolower(label))
  d_max <- if (!is.null(d_max_override)) d_max_override else
            max(abs(range(data$cohens_d, na.rm = TRUE)))
  n_cluster <- if ("n_a" %in% names(data)) data$n_a[1] else NA
  n_label <- if (!is.na(n_cluster)) sprintf(" (n = %s)", n_cluster) else ""

  # Plain Cohen's d plot
  p1 <- data %>%
    ggseg(mapping = aes(fill = cohens_d), atlas = dk, position = "stacked") +
    scale_fill_continuous_divergingx(
      palette = 'RdBu', mid = 0, rev = TRUE,
      limits = c(-d_max, d_max),
      oob = scales::squish,    # clip |d| > d_max to scale endpoints
      na.value = "grey80"
    ) +
    labs(title = paste0("Cluster ", cluster, n_label, " - ", metric, " vs NT"),
         fill = "Cohen's d") +
    theme_minimal()

  fname_base <- paste0("cluster_", cluster, "_", metric, "_", atlas_type, "_vs_td")
  ggsave(file.path(output_dir, paste0(fname_base, ".pdf")), plot = p1,
         width = 10, height = 10)

  # FDR-outlined plot
  sig_labels <- data %>% dplyr::filter(fdr_corrected < 0.05) %>% dplyr::pull(label)
  data <- data %>%
    plyr::mutate(
      outline_color = ifelse(label %in% sig_labels, "black", "white"),
      outline_size  = ifelse(label %in% sig_labels, 1.2,    0.2)
    )

  p2 <- data %>%
    ggseg(mapping = aes(fill = cohens_d, col = outline_color, size = outline_size),
          atlas = "dk", position = "stacked") +
    scale_fill_continuous_divergingx(
      palette = 'RdBu', mid = 0, rev = TRUE,
      limits = c(-d_max, d_max),
      oob = scales::squish,    # clip |d| > d_max to scale endpoints
      na.value = "grey80"
    ) +
    scale_colour_identity() +
    scale_size_identity() +
    labs(title = paste0("Cluster ", cluster, n_label, " - ", metric, " vs NT (FDR p<0.05 outlined)"),
         fill = "Cohen's d") +
    theme_minimal()

  ggsave(file.path(output_dir, paste0(fname_base, "_fdr.pdf")), plot = p2,
         width = 10, height = 10)

  cat(paste0("  Saved: ", fname_base, ".pdf  +  ", fname_base, "_fdr.pdf\n"))
  invisible(list(plain = p1, fdr = p2))
}

#' Create and save a pair of brain plots for subcortical volume using aseg atlas.
#' Fill is **Cohen's d** (Reviewer 3 revision).
plot_subcortical_cluster <- function(data, cluster, metric, atlas_type, output_dir,
                                     d_max_override = NULL) {
  data <- data %>% dplyr::mutate(label = tolower(label))
  data <- rename_subcortical_labels(data)
  d_max <- if (!is.null(d_max_override)) d_max_override else
            max(abs(range(data$cohens_d, na.rm = TRUE)))
  n_cluster <- if ("n_a" %in% names(data)) data$n_a[1] else NA
  n_label <- if (!is.na(n_cluster)) sprintf(" (n = %s)", n_cluster) else ""

  p1 <- data %>%
    ggseg(mapping = aes(fill = cohens_d), atlas = aseg) +
    scale_fill_continuous_divergingx(
      palette = 'RdBu', mid = 0, rev = TRUE,
      limits = c(-d_max, d_max),
      oob = scales::squish,    # clip |d| > d_max to scale endpoints
      na.value = "grey80"
    ) +
    labs(title = paste0("Cluster ", cluster, n_label, " - ", metric, " vs NT"),
         fill = "Cohen's d") +
    theme_minimal()

  fname_base <- paste0("cluster_", cluster, "_", metric, "_", atlas_type, "_vs_td")
  ggsave(file.path(output_dir, paste0(fname_base, ".pdf")), plot = p1,
         width = 10, height = 10)

  sig_labels <- data %>% dplyr::filter(fdr_corrected < 0.05) %>% dplyr::pull(label)
  data <- data %>%
    dplyr::mutate(
      outline_color = ifelse(label %in% sig_labels, "black", "white"),
      outline_size  = ifelse(label %in% sig_labels, 1.7,    0.007)
    )

  p2 <- data %>%
    ggseg(mapping = aes(fill = cohens_d, colour = outline_color, size = outline_size),
          atlas = aseg) +
    scale_fill_continuous_divergingx(
      palette = 'RdBu', mid = 0, rev = TRUE,
      limits = c(-d_max, d_max),
      oob = scales::squish,    # clip |d| > d_max to scale endpoints
      na.value = "grey80"
    ) +
    scale_colour_identity() +
    scale_size_identity() +
    labs(title = paste0("Cluster ", cluster, n_label, " - ", metric, " vs NT (FDR p<0.05 outlined)"),
         fill = "Cohen's d") +
    theme_minimal()

  ggsave(file.path(output_dir, paste0(fname_base, "_fdr.pdf")), plot = p2,
         width = 10, height = 10)

  cat(paste0("  Saved: ", fname_base, ".pdf  +  ", fname_base, "_fdr.pdf\n"))
  invisible(list(plain = p1, fdr = p2))
}

# =============================================================================
# COMBINED 3-PANEL FIGURES
# =============================================================================

#' Attempt to combine three per-cluster plots into a single side-by-side PDF.
#' Uses patchwork if available, then cowplot, then falls back to individual files.
save_combined_figure <- function(plots_list, metric, atlas_type, output_dir,
                                 suffix = "") {
  if (length(plots_list) < 3) {
    cat("  Not enough plots for combined figure — skipping\n")
    return(invisible(NULL))
  }
  p1 <- plots_list[["C1"]]
  p2 <- plots_list[["C2"]]
  p3 <- plots_list[["C3"]]

  fname <- file.path(output_dir,
                     paste0("combined_clusters_", metric, "_", atlas_type,
                            "_vs_td", suffix, ".pdf"))

  if (requireNamespace("patchwork", quietly = TRUE)) {
    combined <- p1 + p2 + p3 + patchwork::plot_layout(ncol = 3)
    ggplot2::ggsave(fname, plot = combined, width = 30, height = 10)
    cat(paste0("  Saved combined (patchwork): ", basename(fname), "\n"))
  } else if (requireNamespace("cowplot", quietly = TRUE)) {
    combined <- cowplot::plot_grid(p1, p2, p3, ncol = 3,
                                   labels = c("C1", "C2", "C3"))
    ggplot2::ggsave(fname, plot = combined, width = 30, height = 10)
    cat(paste0("  Saved combined (cowplot): ", basename(fname), "\n"))
  } else {
    cat("  patchwork and cowplot not available — combined figure skipped\n")
    cat("  Individual per-cluster PDFs already saved.\n")
  }
}

# =============================================================================
# MAIN LOOP
# =============================================================================

clusters <- c("C1", "C2", "C3")

metrics_config <- list(
  list(metric = "thickness", atlas = "dk",   is_cortical = TRUE),
  list(metric = "area",      atlas = "dk",   is_cortical = TRUE),
  list(metric = "grayvol",   atlas = "dk",   is_cortical = TRUE),
  list(metric = "volume",    atlas = "aseg", is_cortical = FALSE)
)

for (cfg in metrics_config) {
  metric      <- cfg$metric
  atlas_type  <- cfg$atlas
  is_cortical <- cfg$is_cortical

  cat(paste0("\n============================================================\n"))
  cat(paste0("Processing: ", metric, " (", atlas_type, ")\n"))
  cat(paste0("============================================================\n"))

  plain_plots <- list()
  fdr_plots   <- list()

  # Pre-load all three clusters so we can share a single colorbar scale
  # (= global max |d| across the three).  Avoids out-of-range clipping while
  # keeping the per-metric cross-cluster comparison valid.
  per_cluster_data <- list()
  for (cluster in clusters) {
    d <- load_cluster_data(cluster, metric, atlas_type)
    if (!is.null(d)) per_cluster_data[[cluster]] <- d
  }
  d_max_metric <- if (length(per_cluster_data) > 0) {
    max(sapply(per_cluster_data, function(d)
      max(abs(d$cohens_d), na.rm = TRUE)), na.rm = TRUE)
  } else NULL
  if (!is.null(d_max_metric)) {
    d_max_metric <- max(d_max_metric, 0.1)
    cat(sprintf("  Shared d_max for %s: %.3f\n", metric, d_max_metric))
  }

  for (cluster in clusters) {
    cat(paste0("\n  Cluster ", cluster, ":\n"))
    data <- per_cluster_data[[cluster]]
    if (is.null(data)) {
      cat(paste0("    No data for cluster ", cluster, " — skipping\n"))
      next
    }

    if (is_cortical) {
      pair <- plot_cortical_cluster(data, cluster, metric, atlas_type, output_dir,
                                    d_max_override = d_max_metric)
    } else {
      pair <- plot_subcortical_cluster(data, cluster, metric, atlas_type, output_dir,
                                       d_max_override = d_max_metric)
    }

    plain_plots[[cluster]] <- pair$plain
    fdr_plots[[cluster]]   <- pair$fdr
  }

  # Combined 3-panel figures
  cat("\n  Creating combined 3-panel figures...\n")
  save_combined_figure(plain_plots, metric, atlas_type, output_dir, suffix = "")
  save_combined_figure(fdr_plots,   metric, atlas_type, output_dir, suffix = "_fdr")
}

# =============================================================================
# COHEN'S d FOREST PLOT (Reviewer 3): per-cluster FDR-significant ROIs
# =============================================================================
#
# For each metric × atlas, gather every FDR-significant ROI from the three
# clusters and render a single forest plot showing Cohen's d (with 95 % CI
# where available). Provides explicit effect-size comparison across clusters
# despite their different sample sizes.

forest_cluster_d <- function(metric, atlas_type, output_dir) {
  clusters <- c("C1", "C2", "C3")
  all_rows <- list()
  for (cluster in clusters) {
    data <- load_cluster_data(cluster, metric, atlas_type)
    if (is.null(data)) next
    sig <- dplyr::filter(data, p_fdr < 0.05)
    if (nrow(sig) == 0) next
    sig$cluster <- cluster
    all_rows[[cluster]] <- sig[, c("cluster", "label", "cohens_d",
                                   "cohens_d_lo", "cohens_d_hi",
                                   "p_fdr", "n_a", "n_b"), drop = FALSE]
  }
  if (length(all_rows) == 0) {
    cat(sprintf("  No FDR-significant ROIs for %s — skipping forest\n", metric))
    return(invisible(NULL))
  }
  forest_df <- dplyr::bind_rows(all_rows)
  # Order ROIs by absolute d (max across clusters)
  ord <- forest_df %>% dplyr::group_by(label) %>%
    dplyr::summarise(abs_d = max(abs(cohens_d), na.rm = TRUE)) %>%
    dplyr::arrange(dplyr::desc(abs_d)) %>% dplyr::pull(label)
  forest_df$label <- factor(forest_df$label, levels = rev(ord))

  cluster_palette <- c("C1" = "#324095", "C2" = "#5CAEE1", "C3" = "#C1C2BC")

  p <- ggplot(forest_df, aes(x = cohens_d, y = label, colour = cluster)) +
    geom_vline(xintercept = 0, colour = "black", linewidth = 0.3) +
    geom_errorbarh(aes(xmin = cohens_d_lo, xmax = cohens_d_hi),
                   height = 0, linewidth = 0.7,
                   position = position_dodge(width = 0.6)) +
    geom_point(size = 2.4, position = position_dodge(width = 0.6)) +
    scale_colour_manual(values = cluster_palette) +
    labs(title = sprintf("Cluster effect sizes (FDR<0.05) — %s [%s]",
                         metric, atlas_type),
         x = "Cohen's d (95 % CI), cluster vs NT",
         y = NULL, colour = NULL) +
    theme_minimal() +
    theme(legend.position = "top",
          plot.title = element_text(size = 12, face = "bold"),
          axis.text.y = element_text(size = 7))

  fname <- file.path(output_dir,
                     sprintf("cluster_%s_%s_cohens_d_forest.pdf", metric, atlas_type))
  ggsave(fname, plot = p, width = 7,
         height = max(3, 0.18 * length(unique(forest_df$label)) + 1.8))
  cat(sprintf("  Saved: %s\n", basename(fname)))

  # Also write the table for the response letter.
  csv_out <- file.path(output_dir,
                       sprintf("cluster_%s_%s_cohens_d_forest.csv", metric, atlas_type))
  write.csv(forest_df, csv_out, row.names = FALSE)
  cat(sprintf("  Saved: %s\n", basename(csv_out)))
}

cat("\n============================================================\n")
cat("Cohen's d forest plots (Reviewer 3)\n")
cat("============================================================\n")
for (cfg in metrics_config) {
  forest_cluster_d(cfg$metric, cfg$atlas, output_dir)
}

cat("\n============================================================\n")
cat("Cluster brain visualizations complete!\n")
cat(paste0("Output directory: ", output_dir, "\n"))
cat("============================================================\n")
