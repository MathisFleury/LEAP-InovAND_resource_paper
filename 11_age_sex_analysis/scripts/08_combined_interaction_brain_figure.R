#!/usr/bin/env Rscript
# =============================================================================
# Combined single-page brain figure (ggseg) — diagnosis x age / x sex
# interaction, full sample.
# =============================================================================
# One column, four stacked ggseg maps per modulator (age, sex):
#   Cortical Thickness | Surface Area | Cortical Volume | Subcortical Volumes
# Fill = interaction t-statistic (RdBu diverging); FDR<0.05 ROIs outlined in
# black. Adapted from ../../4_anatomical_analysis/scripts/13_combined_brain_figure.R
# for the interaction CSVs written by 01_age_sex_interactions.py
#   outputs/figures/r_input_files/t_stat_anat_<atlas>_<metric>_interaction_<model>.csv
#   (headerless: label, t_stat, p_val, cohens_d, p_fdr)
# Output: outputs/figures/combined_interaction_brain_maps_<model>.pdf
# =============================================================================
suppressMessages({
  library(plyr); library(dplyr); library(ggplot2)
  library(colorspace); library(ggseg); library(ggrepel)
  library(patchwork)
})

args        <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
script_dir  <- if (length(script_path) == 0) getwd() else dirname(normalizePath(script_path))
section_dir <- dirname(script_dir)

input_dir  <- file.path(section_dir, "outputs", "figures", "r_input_files")
output_dir <- file.path(section_dir, "outputs", "figures")

n_suffix <- ""
cp <- file.path(input_dir, "group_counts.csv")
if (file.exists(cp)) {
  gc <- read.csv(cp)
  n_suffix <- sprintf("  (Autism n = %s, NT n = %s)", gc$n_autism[1], gc$n_nt[1])
}

PANELS <- list(
  list(title = "Cortical Thickness",  metric = "thickness", atlas = "dk"),
  list(title = "Surface Area",        metric = "area",      atlas = "dk"),
  list(title = "Cortical Volume",     metric = "volume",    atlas = "dk"),
  list(title = "Subcortical Volumes", metric = "volume",    atlas = "aseg")
)

transform_region_column <- function(df) df %>% plyr::mutate(label = case_when(
  endsWith(label, "_left")  ~ paste0("lh_", sub("_left$",  "", label)),
  endsWith(label, "_right") ~ paste0("rh_", sub("_right$", "", label)),
  TRUE ~ label))

aseg_label <- function(df) df %>% dplyr::mutate(label = dplyr::recode(label,
  "Left-Thalamus" = "Left-Thalamus-Proper", "Right-Thalamus" = "Right-Thalamus-Proper",
  "3rd-Ventricle" = "x3rd-ventricle", "4th-Ventricle" = "x4th-ventricle"))

load_metric <- function(model, metric, atlas_type) {
  fp <- file.path(input_dir, sprintf("t_stat_anat_%s_%s_interaction_%s.csv", atlas_type, metric, model))
  if (!file.exists(fp)) { warning(paste("missing", fp)); return(NULL) }
  d <- read.csv(fp, header = FALSE,
                col.names = c("label", "t_stat", "p_val", "cohens_d", "p_fdr"),
                stringsAsFactors = FALSE)
  if (atlas_type == "dk") transform_region_column(d) else aseg_label(d)
}

panel_theme <- theme_void(base_family = "Helvetica") +
  theme(plot.title    = element_text(size = 23, face = "plain", hjust = 0.5,
                                     margin = margin(b = 2, t = 6)),
        legend.position = "right",
        legend.key.height = unit(20, "pt"), legend.key.width = unit(25, "pt"),
        legend.title  = element_text(size = 11), legend.text = element_text(size = 9))

make_panel <- function(spec, model) {
  d <- load_metric(model, spec$metric, spec$atlas)
  if (is.null(d)) return(NULL)
  atlas_obj <- if (spec$atlas == "dk") dk else aseg
  atlas_labels <- unique(as.data.frame(atlas_obj)$label)
  dropped <- setdiff(d$label, atlas_labels)
  if (length(dropped)) cat("  no ggseg region for:", paste(dropped, collapse = ", "), "\n")
  d <- d %>% dplyr::filter(label %in% atlas_labels)
  lim <- max(abs(d$t_stat), na.rm = TRUE)
  d <- d %>% plyr::mutate(
    outline_color = ifelse(!is.na(p_fdr) & p_fdr < 0.05, "black", "white"),
    outline_size  = ifelse(!is.na(p_fdr) & p_fdr < 0.05, 0.8, 0.15))
  p <- ggseg(d, mapping = aes(fill = t_stat, colour = outline_color, size = outline_size),
             atlas = atlas_obj, position = "stacked")
  p + scale_fill_continuous_divergingx(palette = "RdBu", mid = 0, rev = TRUE,
                                       limits = c(-lim, lim), name = "interaction t") +
    scale_colour_identity() + scale_size_identity() + panel_theme + labs(title = spec$title)
}

for (model in c("age", "sex")) {
  panels <- Filter(Negate(is.null), lapply(PANELS, make_panel, model = model))
  if (length(panels) == 0) { warning("No panels for ", model); next }
  combined <- wrap_plots(panels, ncol = 1) +
    plot_annotation(
      title = sprintf("Diagnosis x %s interaction — t-statistic", model),
      subtitle = paste0("FDR p < 0.05 outlined", n_suffix),
      theme = theme(plot.title = element_text(size = 24, face = "plain", hjust = 0.5, family = "Helvetica"),
                    plot.subtitle = element_text(size = 12, hjust = 0.5, family = "Helvetica")))
  out_pdf <- file.path(output_dir, sprintf("combined_interaction_brain_maps_%s.pdf", model))
  ggsave(out_pdf, plot = combined, units = "cm",
         width = 18, height = 18 * (5.2 * length(panels) + 1) / 8, limitsize = FALSE)
  cat("Saved:", out_pdf, "\n")
}
