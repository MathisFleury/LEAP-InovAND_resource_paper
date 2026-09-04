# =============================================================================
# Anatomical MRI Brain Visualizations — Autism vs NT (REVISION, Reviewer 3)
# =============================================================================
#
# v2 of ../../scripts/03_plot_anatomical_mri_brain_visualizations.R, using
# the QC+ComBat+regression FreeSurfer dataset (LEAP_W1 / INOVAND_T1 priority)
# produced by 1_anatomical_mri_autism_nt_v2.py.
#
# Fill: Cohen's d (same convention as the cluster v2 plots).
# Scale: fixed [-0.4, 0.4] across thickness / area / volume for direct
#        comparability with the per-cluster maps.
# FDR-significant regions outlined in black.
# Also emits a Cohen's d forest plot (95 % CI) of FDR-significant ROIs.
#
# Input : ../outputs/figures/r_input_files/t_stat_anat_*_mri_autism_vs_control.csv
# Output: ../outputs/figures/
# =============================================================================

library(plyr)
library(dplyr)
library(ggplot2)
library(colorspace)
library(ggseg)
library(scales)    # provides scales::squish for out-of-bounds clipping

# Colorbar limits are computed per-metric from the loaded data
# (symmetric, = max |d|).  No clipping, no NA / grey regions.

# --- Paths (script-relative) -------------------------------------------------
args        <- commandArgs(trailingOnly = FALSE)
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

brain_map_theme <- theme_minimal() +
  theme(plot.title       = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle    = element_text(size = 13, hjust = 0.5),
        legend.title     = element_text(size = 12, face = "bold"),
        legend.text      = element_text(size = 10),
        legend.position  = "bottom",
        panel.grid       = element_blank(),
        axis.text        = element_blank(),
        axis.ticks       = element_blank(),
        strip.text       = element_text(size = 11, face = "bold"),
        strip.background = element_rect(fill = "lightgray", color = "black"))


transform_region_column <- function(df) {
  # Accept either "<region>_left"/"<region>_right" or "lh_<region>"/"rh_<region>"
  df %>% plyr::mutate(label = case_when(
    endsWith(label, "_left")  ~ paste0("lh_", sub("_left$",  "", label)),
    endsWith(label, "_right") ~ paste0("rh_", sub("_right$", "", label)),
    TRUE ~ label
  ))
}


rename_subcortical_labels <- function(df) {
  df %>% dplyr::mutate(label = case_when(
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
    TRUE ~ label
  ))
}


load_t_stat_data <- function(metric, atlas_type) {
  fname <- paste0("t_stat_anat_", atlas_type, "_", metric,
                  "_mri_autism_vs_control.csv")
  fpath <- file.path(input_dir, fname)
  if (!file.exists(fpath)) {
    warning(paste("File not found:", fpath))
    return(NULL)
  }
  data <- read.csv(fpath, stringsAsFactors = FALSE)
  required <- c("label", "t_stat", "p_val", "cohens_d", "p_fdr")
  if (!all(required %in% names(data))) {
    data <- read.csv(fpath, header = FALSE,
                     col.names = required, stringsAsFactors = FALSE)
  }
  for (cc in c("cohens_d_lo", "cohens_d_hi", "n_a", "n_b", "p_bonf")) {
    if (!cc %in% names(data)) data[[cc]] <- NA
  }
  data %>% dplyr::mutate(
    fdr_corrected  = p_fdr,
    bonf_corrected = p_bonf,
    significant    = p_fdr < 0.05
  )
}


# =============================================================================
# CORTICAL  (thickness + area)
# =============================================================================
#' Symmetric ggseg brain map for any continuous fill column ("cohens_d" or
#' "t_stat").  Emits three variants per fill column:
#'   - plain (no outlines)
#'   - FDR-outlined        (suffix "_fdr")
#'   - Bonferroni-outlined (suffix "_bonf")
.plot_cortical_one <- function(data, metric, fill_col, prefix, label,
                               sig_fdr_lbl, sig_bonf_lbl) {
  vmax <- max(abs(range(data[[fill_col]], na.rm = TRUE)), 0.1, na.rm = TRUE)
  cat(sprintf("  %s_max for %s (dk): %.3f\n", fill_col, metric, vmax))
  sub <- sprintf("%s (Autism - NT); scale = max = %.2f", label, vmax)

  # --- plain map ---
  p1 <- data %>%
    ggseg(mapping = aes(fill = .data[[fill_col]]), atlas = dk,
          position = "stacked") +
    scale_fill_continuous_divergingx(palette = "RdBu", mid = 0, rev = TRUE,
                                     limits = c(-vmax, vmax),
                                     oob = scales::squish,
                                     na.value = "grey80",
                                     name = label) +
    brain_map_theme +
    labs(title = paste0("Autism vs NT - ", metric, " (Desikan)"),
         subtitle = sub)
  out1 <- file.path(output_dir,
                    paste0(prefix, "_", metric, "_dk.pdf"))
  ggsave(out1, plot = p1, width = 10, height = 10)
  cat(paste0("  Saved: ", basename(out1), "\n"))

  # --- outline variants ---
  render_outlined <- function(sig_lbl, suffix, badge) {
    if (length(sig_lbl) == 0) {
      cat(paste0("  [skip] no ", badge, "-significant ROI for ",
                 metric, " (", fill_col, ")\n"))
      return(invisible(NULL))
    }
    d2 <- data %>% plyr::mutate(
      outline_color = ifelse(label %in% sig_lbl, "black", "white"),
      outline_size  = ifelse(label %in% sig_lbl, 1.2, 0.2))
    p <- d2 %>%
      ggseg(mapping = aes(fill = .data[[fill_col]],
                          col = outline_color,
                          size = outline_size),
            atlas = "dk", position = "stacked") +
      scale_fill_continuous_divergingx(palette = "RdBu", mid = 0, rev = TRUE,
                                       limits = c(-vmax, vmax),
                                       oob = scales::squish,
                                       na.value = "grey80",
                                       name = label) +
      scale_colour_identity() + scale_size_identity() +
      brain_map_theme +
      labs(title = paste0("Autism vs NT - ", metric, " (Desikan)"),
           subtitle = paste0(sub, " (", badge, "<0.05 outlined)"))
    out <- file.path(output_dir,
                     paste0(prefix, "_", metric, "_dk_", suffix, ".pdf"))
    ggsave(out, plot = p, width = 10, height = 10)
    cat(paste0("  Saved: ", basename(out), "\n"))
  }
  render_outlined(sig_fdr_lbl,  "fdr",  "FDR")
  render_outlined(sig_bonf_lbl, "bonf", "Bonferroni")
}


plot_cortical <- function(metric) {
  data <- load_t_stat_data(metric, "dk")
  if (is.null(data)) return(NULL)
  data <- transform_region_column(data) %>% plyr::mutate(label = tolower(label))

  sig_fdr  <- data %>% dplyr::filter(fdr_corrected  < 0.05) %>%
    dplyr::pull(label)
  sig_bonf <- data %>% dplyr::filter(bonf_corrected < 0.05) %>%
    dplyr::pull(label)
  cat(sprintf("  %s: FDR<0.05 = %d ROIs · Bonferroni<0.05 = %d ROIs\n",
              metric, length(sig_fdr), length(sig_bonf)))

  # --- Cohen's d (primary) -------------------------------------------------
  .plot_cortical_one(data, metric,
                     fill_col = "cohens_d",
                     prefix   = "autism_vs_nt",
                     label    = "Cohen's d",
                     sig_fdr_lbl  = sig_fdr,
                     sig_bonf_lbl = sig_bonf)
  # --- t-statistic (companion) --------------------------------------------
  .plot_cortical_one(data, metric,
                     fill_col = "t_stat",
                     prefix   = "autism_vs_nt_tstat",
                     label    = "t-statistic",
                     sig_fdr_lbl  = sig_fdr,
                     sig_bonf_lbl = sig_bonf)
  invisible(data)
}


# =============================================================================
# SUBCORTICAL
# =============================================================================
.plot_subcortical_one <- function(data, fill_col, prefix, label,
                                  sig_fdr_lbl, sig_bonf_lbl) {
  vmax <- max(abs(range(data[[fill_col]], na.rm = TRUE)), 0.1, na.rm = TRUE)
  cat(sprintf("  %s_max for volume (aseg): %.3f\n", fill_col, vmax))
  sub <- sprintf("%s (Autism - NT); scale = max = %.2f", label, vmax)

  p1 <- data %>%
    ggseg(mapping = aes(fill = .data[[fill_col]]), atlas = aseg) +
    scale_fill_continuous_divergingx(palette = "RdBu", mid = 0, rev = TRUE,
                                     limits = c(-vmax, vmax),
                                     oob = scales::squish,
                                     na.value = "grey80",
                                     name = label) +
    brain_map_theme +
    labs(title = "Autism vs NT - Subcortical Volume (aseg)",
         subtitle = sub)
  out1 <- file.path(output_dir, paste0(prefix, "_volume_aseg.pdf"))
  ggsave(out1, plot = p1, width = 10, height = 10)
  cat(paste0("  Saved: ", basename(out1), "\n"))

  render_outlined <- function(sig_lbl, suffix, badge) {
    if (length(sig_lbl) == 0) {
      cat(paste0("  [skip] no ", badge, "-significant subcortical ROI\n"))
      return(invisible(NULL))
    }
    d2 <- data %>% dplyr::mutate(
      outline_color = ifelse(label %in% sig_lbl, "black", "white"),
      outline_size  = ifelse(label %in% sig_lbl, 1.7, 0.007))
    p <- d2 %>%
      ggseg(mapping = aes(fill = .data[[fill_col]],
                          colour = outline_color,
                          size = outline_size),
            atlas = aseg) +
      scale_fill_continuous_divergingx(palette = "RdBu", mid = 0, rev = TRUE,
                                       limits = c(-vmax, vmax),
                                       oob = scales::squish,
                                       na.value = "grey80",
                                       name = label) +
      scale_colour_identity() + scale_size_identity() +
      brain_map_theme +
      labs(title = "Autism vs NT - Subcortical Volume (aseg)",
           subtitle = paste0(sub, " (", badge, "<0.05 outlined)"))
    out <- file.path(output_dir,
                     paste0(prefix, "_volume_aseg_", suffix, ".pdf"))
    ggsave(out, plot = p, width = 10, height = 10)
    cat(paste0("  Saved: ", basename(out), "\n"))
  }
  render_outlined(sig_fdr_lbl,  "fdr",  "FDR")
  render_outlined(sig_bonf_lbl, "bonf", "Bonferroni")
}


plot_subcortical <- function() {
  data <- load_t_stat_data("volume", "aseg")
  if (is.null(data)) return(NULL)
  data <- data %>% dplyr::mutate(label = tolower(label))
  data <- rename_subcortical_labels(data)
  sig_fdr  <- data %>% dplyr::filter(fdr_corrected  < 0.05) %>%
    dplyr::pull(label)
  sig_bonf <- data %>% dplyr::filter(bonf_corrected < 0.05) %>%
    dplyr::pull(label)
  cat(sprintf("  volume: FDR<0.05 = %d ROIs · Bonferroni<0.05 = %d ROIs\n",
              length(sig_fdr), length(sig_bonf)))
  .plot_subcortical_one(data, "cohens_d", "autism_vs_nt",
                        "Cohen's d", sig_fdr, sig_bonf)
  .plot_subcortical_one(data, "t_stat", "autism_vs_nt_tstat",
                        "t-statistic", sig_fdr, sig_bonf)
  invisible(data)
}


# =============================================================================
# FOREST  (Cohen's d + 95 % CI of FDR-significant ROIs across metrics)
# =============================================================================
build_forest_table <- function() {
  rows <- list()
  for (cfg in list(
        list(metric = "thickness", atlas = "dk"),
        list(metric = "area",      atlas = "dk"),
        list(metric = "grayvol",   atlas = "dk"),
        list(metric = "volume",    atlas = "aseg"))) {
    data <- load_t_stat_data(cfg$metric, cfg$atlas)
    if (is.null(data)) next
    sig <- dplyr::filter(data, p_fdr < 0.05)
    if (nrow(sig) == 0) next
    sig$metric <- cfg$metric
    cols_keep <- c("metric", "label", "cohens_d", "cohens_d_lo",
                   "cohens_d_hi", "p_fdr", "n_a", "n_b")
    rows[[cfg$metric]] <- sig[, cols_keep, drop = FALSE]
  }
  if (length(rows) == 0) return(NULL)
  dplyr::bind_rows(rows)
}


plot_forest <- function() {
  df <- build_forest_table()
  if (is.null(df) || nrow(df) == 0) {
    cat("  No FDR-significant ROI — forest skipped\n")
    return(invisible(NULL))
  }
  # Make y-axis labels unique across metrics (some region names — e.g.
  # Brain-Stem or lateral ventricle — can appear in more than one input).
  df$y_label <- paste0(df$label, "  [", df$metric, "]")
  # Order by |d|, descending
  df <- df %>% dplyr::arrange(dplyr::desc(abs(cohens_d)))
  df$y_label <- factor(df$y_label, levels = rev(df$y_label))

  palette <- c("thickness" = "#324095", "area" = "#5CAEE1", "volume" = "#C1382F")

  p <- ggplot(df, aes(x = cohens_d, y = y_label, colour = metric)) +
    geom_vline(xintercept = 0, colour = "black", linewidth = 0.3) +
    geom_errorbarh(aes(xmin = cohens_d_lo, xmax = cohens_d_hi),
                   height = 0, linewidth = 0.6) +
    geom_point(size = 2.2) +
    scale_colour_manual(values = palette) +
    labs(title = "Autism vs NT — FDR-significant Cohen's d (95 % CI)",
         x = "Cohen's d (Autism − NT)",
         y = NULL, colour = NULL) +
    theme_minimal() +
    theme(legend.position = "top",
          plot.title = element_text(size = 12, face = "bold"),
          axis.text.y = element_text(size = 7))

  fname <- file.path(output_dir, "autism_vs_nt_cohens_d_forest.pdf")
  ggsave(fname, plot = p, width = 7,
         height = max(3, 0.18 * nrow(df) + 1.8))
  csv_out <- file.path(output_dir, "autism_vs_nt_cohens_d_forest.csv")
  write.csv(df, csv_out, row.names = FALSE)
  cat(paste0("  Saved: ", basename(fname), "  +  ", basename(csv_out), "\n"))
}


# =============================================================================
# MAIN
# =============================================================================
cat("============================================================\n")
cat("Autism vs NT brain visualisations (REVISION)\n")
cat("============================================================\n")
plot_cortical("thickness")
plot_cortical("area")
plot_cortical("grayvol")    # per-Desikan cortical gray-matter volume
plot_subcortical()
cat("\nForest plot of FDR-significant Cohen's d:\n")
plot_forest()
cat("\nDone.\n")
