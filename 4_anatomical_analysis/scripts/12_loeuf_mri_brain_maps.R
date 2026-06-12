#!/usr/bin/env Rscript
# =============================================================================
# 12 - LOEUF vs MRI Brain Maps (ggseg) — delloF + Miss, SynGO + CHROM
# =============================================================================
# Reads genetics_mri_correlation_statistics.csv produced by script 05 (hg19)
# and loeuf_mri_interaction_statistics_hg38.csv produced by script 11 (hg38).
# Creates ggseg brain maps (Desikan-Killiany atlas) for:
#   variant types: dellof_constrained, miss_constrained, dellof_all, miss_all
#   gene lists:    SynGO, CHROM (excluding combined syngo_chromepitf)
#   metrics:       thickness, area
#
# Significant regions (FDR-corrected p < 0.05) highlighted with thick black
# border.  FDR is applied across all CHROM + SynGO regions combined within
# each variant-type x metric combination.
#
# Outputs:
#   figures_genetics/{variant_type}_{metric}_brain_map.pdf      (hg19)
#   figures_genetics_hg38/{variant_type}_{metric}_brain_map.pdf (hg38)
# =============================================================================

library(ggseg)
library(plyr)
library(dplyr)
library(ggplot2)
library(colorspace)
library(tidyverse)
library(gridExtra)
library(grid)

args       <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("--file=", args, value = TRUE)
if (length(script_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("--file=", "", script_arg)))
} else {
  script_dir <- getwd()
}

output_dir_hg19 <- file.path(script_dir, "..", "outputs", "figures_genetics")
output_dir_hg38 <- file.path(script_dir, "..", "outputs", "figures_genetics_hg38")
tables_hg19     <- file.path(script_dir, "..", "outputs", "tables_genetics")
tables_hg38     <- file.path(script_dir, "..", "outputs", "tables_genetics_hg38")
dir.create(output_dir_hg19, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir_hg38, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Helper: transform MRI column names to ggseg label format
# ---------------------------------------------------------------------------
prepare_ggseg_data <- function(feature_data, metric) {
  # Build label: "left" -> "lh_regionname", "right" -> "rh_regionname"
  ggseg_data <- feature_data %>%
    mutate(
      label = ifelse(hemisphere == "left",
                     paste0("lh_", gsub(" ", "", region)),
                     paste0("rh_", gsub(" ", "", region)))
    ) %>%
    select(label, value = correlation, p_unc = p_value)

  return(ggseg_data)
}

# ---------------------------------------------------------------------------
# Helper: short display label from genetic_feature column name
# ---------------------------------------------------------------------------
gene_list_name <- function(feat) {
  if (grepl("chromepitf", feat, ignore.case = TRUE)) "CHROM" else "SYNGO"
}

# ---------------------------------------------------------------------------
# Create and save brain maps for one variant type x metric combination
# ---------------------------------------------------------------------------
create_brain_plot <- function(stats_df, genetic_features, metric, variant_type,
                              out_dir) {
  cat(sprintf("  Creating %s | %s brain map...\n", variant_type, metric))

  variant_data <- stats_df %>%
    dplyr::filter(mri_type == metric, genetic_feature %in% genetic_features)

  if (nrow(variant_data) == 0) {
    cat(sprintf("    No data for %s %s — skipping\n", variant_type, metric))
    return(invisible(NULL))
  }

  # Build per-gene-list ggseg data frames
  ggseg_list <- list()
  for (feat in genetic_features) {
    feat_data <- variant_data %>% dplyr::filter(genetic_feature == !!feat)
    if (nrow(feat_data) == 0) next

    gl_name   <- gene_list_name(feat)
    gd        <- prepare_ggseg_data(feat_data, metric)
    gd$gene_list         <- gl_name
    gd$carrier_count     <- feat_data$sample_size[1]
    gd$gene_list_label   <- paste0(gl_name, "\n(n=", feat_data$sample_size[1], ")")

    ggseg_list[[gl_name]] <- gd
  }

  if (length(ggseg_list) == 0) return(invisible(NULL))

  combined <- do.call(rbind, ggseg_list)
  combined$gene_list       <- factor(combined$gene_list, levels = c("CHROM", "SYNGO"))
  combined$gene_list_label <- factor(combined$gene_list_label)

  # FDR correction across all regions and gene lists for this subset
  valid_idx <- !is.na(combined$p_unc)
  if (sum(valid_idx) > 0) {
    combined$p_fdr <- NA_real_
    combined$p_fdr[valid_idx] <- p.adjust(combined$p_unc[valid_idx], method = "fdr")
  } else {
    combined$p_fdr <- combined$p_unc
  }

  # Build plot data: keep grouped by gene_list_label so ggseg processes each
  # group's atlas join separately — prevents a phantom NA facet for unmatched
  # atlas regions (known ggseg + facet_wrap interaction)
  plot_data <- combined %>%
    dplyr::filter(!is.na(value), !is.na(gene_list_label)) %>%
    group_by(gene_list_label) %>%
    mutate(
      border_color = ifelse(is.na(p_fdr) | p_fdr >= 0.05, "transparent", "black"),
      border_size  = ifelse(is.na(p_fdr) | p_fdr >= 0.05, 0.5, 1.5),
      significance = ifelse(is.na(p_fdr) | p_fdr >= 0.05, "Non-significant", "Significant")
    ) %>%
    select(label, value, border_color, border_size, gene_list_label)

  p <- plot_data %>%
    ggseg(mapping = aes(fill = value, colour = border_color, size = border_size),
          atlas = dk, position = "stacked") +
    scale_fill_continuous_divergingx(palette = "RdBu", mid = 0, rev = TRUE,
                                     name = "Pearson Correlation") +
    scale_colour_identity() +
    scale_size_identity() +
    facet_wrap(~gene_list_label, ncol = 2) +
    theme_minimal() +
    theme(
      plot.title       = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle    = element_text(size = 14, hjust = 0.5),
      legend.title     = element_text(size = 12, face = "bold"),
      legend.text      = element_text(size = 10),
      legend.position  = "bottom",
      panel.grid       = element_blank(),
      axis.text        = element_blank(),
      axis.ticks       = element_blank(),
      strip.text       = element_text(size = 11, face = "bold"),
      strip.background = element_rect(fill = "lightgray", color = "black")
    ) +
    labs(title    = sprintf("-log(LOEUF) %s vs %s Pearson Correlations",
                            gsub("_", " ", toupper(variant_type)),
                            gsub("_", " ", toupper(metric))),
         subtitle = "Pearson correlations for Chrom and Syngo gene lists (significant regions highlighted with thick black border, FDR-corrected p < 0.05)")

  fname <- file.path(out_dir, sprintf("%s_%s_brain_map.pdf", variant_type, metric))
  ggsave(fname, plot = p, width = 12, height = 8, dpi = 300)
  cat(sprintf("    Saved: %s\n", basename(fname)))

  return(invisible(p))
}

# ---------------------------------------------------------------------------
# Main analysis loop
# ---------------------------------------------------------------------------
run_analysis <- function(stats_file, output_dir, variant_types, label = "hg19") {
  if (!file.exists(stats_file)) {
    cat(sprintf("  Skipping %s: file not found (%s)\n", label, stats_file))
    return(invisible(NULL))
  }

  cat(sprintf("\nLoading %s correlation statistics...\n", label))
  stats_df <- read.csv(stats_file)
  cat(sprintf("  Loaded %d rows\n", nrow(stats_df)))

  # Filter: LOEUF only, CHROM and SynGO only, exclude combined list
  loeuf_df <- stats_df %>%
    dplyr::filter(genetic_type == "LOEUF") %>%
    dplyr::filter(grepl("chromepitf|syngo", genetic_feature, ignore.case = TRUE)) %>%
    dplyr::filter(!grepl("syngo_chromepitf", genetic_feature, ignore.case = TRUE))

  cat(sprintf("  %d rows after filtering (CHROM + SynGO only)\n", nrow(loeuf_df)))

  for (vtype in variant_types) {
    cat(sprintf("\n=== %s | %s ===\n", label, vtype))

    # Identify genetic features belonging to this variant type
    if (vtype == "dellof_constrained") {
      features <- unique(loeuf_df$genetic_feature[
        grepl("^dellof_.*constraint.*_score$", loeuf_df$genetic_feature)])
    } else if (vtype == "miss_constrained") {
      features <- unique(loeuf_df$genetic_feature[
        grepl("^miss_.*constraint.*_score$", loeuf_df$genetic_feature)])
    } else if (vtype == "dellof_all") {
      features <- unique(loeuf_df$genetic_feature[
        grepl("^dellof_.*_score$", loeuf_df$genetic_feature) &
        !grepl("constraint", loeuf_df$genetic_feature)])
    } else if (vtype == "miss_all") {
      features <- unique(loeuf_df$genetic_feature[
        grepl("^miss_.*_score$", loeuf_df$genetic_feature) &
        !grepl("constraint", loeuf_df$genetic_feature)])
    } else {
      # hg38: match by prefix
      features <- unique(loeuf_df$genetic_feature[
        grepl(vtype, loeuf_df$genetic_feature, fixed = FALSE)])
    }

    if (length(features) == 0) {
      cat(sprintf("    No features found for %s — skipping\n", vtype))
      next
    }
    cat(sprintf("    Features: %s\n", paste(features, collapse = ", ")))

    for (metric in c("thickness", "area")) {
      create_brain_plot(loeuf_df, features, metric, vtype, output_dir)
    }
  }
}

# ---------------------------------------------------------------------------
# hg19: reads from script 05's genetics_mri_correlation_statistics.csv
# ---------------------------------------------------------------------------
hg19_stats  <- file.path(tables_hg19, "genetics_mri_correlation_statistics.csv")
hg19_vtypes <- c("dellof_constrained", "miss_constrained", "dellof_all", "miss_all")
run_analysis(hg19_stats, output_dir_hg19, hg19_vtypes, label = "hg19")

# ---------------------------------------------------------------------------
# hg38: reads from script 11's loeuf_mri_interaction_statistics_hg38.csv
# (hg38 has lof_constrained + del_constrained only, no miss_)
# ---------------------------------------------------------------------------
hg38_stats  <- file.path(tables_hg38, "loeuf_mri_interaction_statistics_hg38.csv")

if (file.exists(hg38_stats)) {
  cat(sprintf("\nLoading hg38 statistics...\n"))
  hg38_df <- read.csv(hg38_stats)
  cat(sprintf("  Loaded %d rows\n", nrow(hg38_df)))

  # hg38 stats have gene_list_label instead of variant_type; add genetic_type
  if (!"genetic_type" %in% colnames(hg38_df)) {
    hg38_df$genetic_type <- "LOEUF"
  }

  # Filter for CHROM and SynGO only
  hg38_filtered <- hg38_df %>%
    dplyr::filter(grepl("chromepitf|syngo", genetic_feature, ignore.case = TRUE)) %>%
    dplyr::filter(!grepl("syngo_chromepitf", genetic_feature, ignore.case = TRUE))

  # dellof combined (lof + del, most constrained gene) — constrained
  dellof_c_features <- unique(hg38_filtered$genetic_feature[
    grepl("^dellof_.*constraint.*_score.*_hg38$", hg38_filtered$genetic_feature)])
  # dellof combined — all carriers
  dellof_a_features <- unique(hg38_filtered$genetic_feature[
    grepl("^dellof_.*_score.*_hg38$", hg38_filtered$genetic_feature) &
    !grepl("constraint", hg38_filtered$genetic_feature)])

  hg38_subsets <- list(
    dellof_constrained = dellof_c_features,
    dellof_all         = dellof_a_features
  )
  hg38_subsets <- Filter(function(x) length(x) > 0, hg38_subsets)

  for (vtype in names(hg38_subsets)) {
    features <- hg38_subsets[[vtype]]
    if (length(features) == 0) next
    cat(sprintf("\n=== hg38 | %s ===\n", vtype))
    cat(sprintf("    Features: %s\n", paste(features, collapse = ", ")))
    for (metric in c("thickness", "area")) {
      create_brain_plot(hg38_filtered, features, metric, vtype, output_dir_hg38)
    }
  }
} else {
  cat(sprintf("\n  Skipping hg38: %s not found (run script 11 first)\n", hg38_stats))
}

cat("\n=== Brain maps complete ===\n")
cat(sprintf("  hg19: %s\n", output_dir_hg19))
cat(sprintf("  hg38: %s\n", output_dir_hg38))
