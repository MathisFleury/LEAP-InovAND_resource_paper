#!/usr/bin/env Rscript
# =============================================================================
# LOEUF (hg38) vs MRI Brain Maps (ggseg) — REVISION
# =============================================================================
# v2 of ../../4-1_anatomical_analysis/scripts/12_loeuf_mri_brain_maps.R, hg38
# branch only, consuming the stats CSV produced by
# 05_loeuf_mri_correlations_hg38_v2.py on the QC+ComBat+regression MRI table.
#
# Input  : ../outputs/tables_genetics_hg38/loeuf_mri_interaction_statistics_hg38.csv
# Output : ../outputs/figures_genetics_hg38/dellof_{all,constrained}_{thickness,area}_brain_map.pdf
# =============================================================================

library(ggseg)
library(plyr)
library(dplyr)
library(ggplot2)
library(colorspace)
library(tidyverse)

args       <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("--file=", args, value = TRUE)
if (length(script_arg) > 0) {
  script_dir <- dirname(normalizePath(sub("--file=", "", script_arg)))
} else {
  script_dir <- getwd()
}

section_dir <- dirname(script_dir)
tables_dir  <- file.path(section_dir, "outputs", "tables_genetics_hg38")
output_dir  <- file.path(section_dir, "outputs", "figures_genetics_hg38")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

stats_file <- file.path(tables_dir, "loeuf_mri_interaction_statistics_hg38.csv")

# ---------------------------------------------------------------------------
prepare_ggseg_data <- function(feature_data) {
  feature_data %>%
    mutate(label = ifelse(hemisphere == "left",
                          paste0("lh_", gsub(" ", "", region)),
                          paste0("rh_", gsub(" ", "", region)))) %>%
    select(label, value = correlation, p_unc = p_value)
}

gene_list_name <- function(feat) {
  if (grepl("chromepitf", feat, ignore.case = TRUE)) "CHROM" else "SYNGO"
}

# Map raw subcortical column names to ggseg `aseg` labels. The canonical
# mapping lives in 02_plot_anatomical_mri_brain_v2.R::rename_subcortical_labels;
# we keep an identical copy here so this script is self-contained but stays in
# lock-step with the rest of the revision pipeline. Update both if either
# changes.
rename_subcortical_labels <- function(df) {
  df %>% dplyr::mutate(label = dplyr::case_when(
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
    label == "cc_posterior"                  ~ "CC_Posterior",
    label == "cc_mid_posterior"              ~ "CC_Mid_Posterior",
    label == "cc_central"                    ~ "CC_Central",
    label == "cc_mid_anterior"               ~ "CC_Mid_Anterior",
    label == "cc_anterior"                   ~ "CC_Anterior",
    label == "left-cerebellum-cortex"        ~ "Left-Cerebellum-Cortex",
    label == "right-cerebellum-cortex"       ~ "Right-Cerebellum-Cortex",
    label == "left-cerebellum-white-matter"  ~ "Left-Cerebellum-White-Matter",
    label == "right-cerebellum-white-matter" ~ "Right-Cerebellum-White-Matter",
    label == "brain-stem"                    ~ "Brain-Stem",
    label == "left-lateral-ventricle"        ~ "Left-Lateral-Ventricle",
    label == "right-lateral-ventricle"       ~ "Right-Lateral-Ventricle",
    label == "3rd-ventricle"                 ~ "x3rd-ventricle",
    label == "4th-ventricle"                 ~ "x4th-ventricle",
    TRUE ~ label))
}

prepare_aseg_data <- function(feature_data) {
  feature_data %>%
    dplyr::mutate(label = tolower(mri_col)) %>%
    rename_subcortical_labels() %>%
    dplyr::select(label, value = correlation, p_unc = p_value)
}

create_brain_plot <- function(stats_df, genetic_features, metric, variant_type,
                              out_dir, atlas_obj = dk, prep_fn = prepare_ggseg_data) {
  cat(sprintf("  Creating %s | %s brain map...\n", variant_type, metric))

  variant_data <- stats_df %>%
    dplyr::filter(mri_type == metric, genetic_feature %in% genetic_features)
  if (nrow(variant_data) == 0) {
    cat(sprintf("    No data for %s %s — skipping\n", variant_type, metric))
    return(invisible(NULL))
  }

  ggseg_list <- list()
  for (feat in genetic_features) {
    feat_data <- variant_data %>% dplyr::filter(genetic_feature == !!feat)
    if (nrow(feat_data) == 0) next
    gl_name <- gene_list_name(feat)
    gd <- prep_fn(feat_data)
    gd$gene_list       <- gl_name
    gd$carrier_count   <- feat_data$sample_size[1]
    gd$gene_list_label <- paste0(gl_name, "\n(n=", feat_data$sample_size[1], ")")
    ggseg_list[[gl_name]] <- gd
  }
  if (length(ggseg_list) == 0) return(invisible(NULL))

  combined <- do.call(rbind, ggseg_list)
  combined$gene_list       <- factor(combined$gene_list, levels = c("CHROM", "SYNGO"))
  combined$gene_list_label <- factor(combined$gene_list_label)

  # FDR per panel (one ggseg facet = one gene list × this metric).
  # p.adjust passes NA through and excludes them from the denominator.
  combined <- combined %>%
    dplyr::group_by(gene_list_label) %>%
    dplyr::mutate(p_fdr = p.adjust(p_unc, method = "fdr")) %>%
    dplyr::ungroup()

  # Outline threshold: per-panel FDR q < 0.05.
  plot_data <- combined %>%
    dplyr::filter(!is.na(value), !is.na(gene_list_label)) %>%
    group_by(gene_list_label) %>%
    mutate(border_color = ifelse(is.na(p_fdr) | p_fdr >= 0.05, "transparent", "black"),
           border_size  = ifelse(is.na(p_fdr) | p_fdr >= 0.05, 0.5, 1.5),
           significance = ifelse(is.na(p_fdr) | p_fdr >= 0.05, "Non-significant", "Significant")) %>%
    select(label, value, border_color, border_size, gene_list_label)

  p <- plot_data %>%
    ggseg(mapping = aes(fill = value, colour = border_color, size = border_size),
          atlas = atlas_obj, position = "stacked") +
    scale_fill_continuous_divergingx(palette = "RdBu", mid = 0, rev = TRUE,
                                     name = "Pearson Correlation") +
    scale_colour_identity() +
    scale_size_identity() +
    facet_wrap(~gene_list_label, ncol = 2) +
    theme_minimal() +
    theme(plot.title       = element_text(size = 16, face = "bold", hjust = 0.5),
          plot.subtitle    = element_text(size = 13, hjust = 0.5),
          legend.title     = element_text(size = 12, face = "bold"),
          legend.text      = element_text(size = 10),
          legend.position  = "bottom",
          panel.grid       = element_blank(),
          axis.text        = element_blank(),
          axis.ticks       = element_blank(),
          strip.text       = element_text(size = 11, face = "bold"),
          strip.background = element_rect(fill = "lightgray", color = "black")) +
    labs(title    = sprintf("-log(LOEUF) %s vs %s Pearson Correlations (revision)",
                            gsub("_", " ", toupper(variant_type)),
                            gsub("_", " ", toupper(metric))),
         subtitle = "QC+ComBat+regression MRI · FDR per panel (q<0.05 outlined in black)")

  fname <- file.path(out_dir, sprintf("%s_%s_brain_map.pdf", variant_type, metric))
  ggsave(fname, plot = p, width = 12, height = 8, dpi = 300)
  cat(sprintf("    Saved: %s\n", basename(fname)))
  invisible(p)
}

# ---------------------------------------------------------------------------
if (!file.exists(stats_file)) {
  stop(sprintf("Stats file not found: %s\nRun 05_loeuf_mri_correlations_hg38_v2.py first.",
               stats_file))
}

cat(sprintf("\nLoading hg38 statistics: %s\n", stats_file))
hg38_df <- read.csv(stats_file)
cat(sprintf("  Loaded %d rows\n", nrow(hg38_df)))
if (!"genetic_type" %in% colnames(hg38_df)) hg38_df$genetic_type <- "LOEUF"

hg38_filtered <- hg38_df %>%
  dplyr::filter(grepl("chromepitf|syngo", genetic_feature, ignore.case = TRUE)) %>%
  dplyr::filter(!grepl("syngo_chromepitf", genetic_feature, ignore.case = TRUE))

# All-carriers only — the constrained subset is dropped upstream (in the
# Python stats step) to avoid inflating FDR when both subsets share a pool.
dellof_a_features <- unique(hg38_filtered$genetic_feature[
  grepl("^dellof_.*_score.*_hg38$", hg38_filtered$genetic_feature) &
  !grepl("constraint", hg38_filtered$genetic_feature)])

hg38_subsets <- list(dellof_all = dellof_a_features)
hg38_subsets <- Filter(function(x) length(x) > 0, hg38_subsets)

for (vtype in names(hg38_subsets)) {
  features <- hg38_subsets[[vtype]]
  cat(sprintf("\n=== %s ===\n", vtype))
  cat(sprintf("    Features: %s\n", paste(features, collapse = ", ")))
  for (metric in c("thickness", "area")) {
    create_brain_plot(hg38_filtered, features, metric, vtype, output_dir,
                      atlas_obj = dk, prep_fn = prepare_ggseg_data)
  }
  create_brain_plot(hg38_filtered, features, "subcortical", vtype, output_dir,
                    atlas_obj = aseg, prep_fn = prepare_aseg_data)
}

cat(sprintf("\nBrain maps written to: %s\n", output_dir))
