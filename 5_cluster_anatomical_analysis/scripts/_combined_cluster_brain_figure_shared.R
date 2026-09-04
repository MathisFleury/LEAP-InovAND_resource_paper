# =============================================================================
# Combined single-page brain figure (ggseg) — clusters vs TD
# =============================================================================
# Cluster analogue of 4_anatomical_analysis/scripts/13_combined_brain_figure.R.
# For each cluster (C1, C2, C3) builds ONE single-column figure stacking the
# available feature maps:
#   Cortical Thickness | Surface Area | Cortical Volume | Subcortical Volumes
# Each panel: divergingx RdBu fill, FDR-significant ROIs outlined in black, a
# per-feature title. Two figures per cluster: t-statistic and Cohen's d.
#
# Input files (one per cluster x metric x atlas):
#   t_stat_cluster_<C>_<metric>_<atlas>_mri_cluster_vs_td.csv
# Cortical volume (grayvol_dk) is only present in the curated pipeline; if
# absent the Cortical Volume panel is simply skipped.
#
# Env overrides (default to this section's outputs):
#   ANAT_INPUT_DIR, ANAT_OUTPUT_DIR, ANAT_SCALE_MAX (fixed |Cohen's d| limit)
# =============================================================================

suppressMessages({
  library(plyr); library(dplyr); library(ggplot2)
  library(colorspace); library(ggseg); library(patchwork); library(scales)
})

# --- Resolve paths / params --------------------------------------------------
args        <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
script_dir  <- if (length(script_path) == 0) getwd() else dirname(normalizePath(script_path))
section_dir <- dirname(script_dir)

input_dir  <- Sys.getenv("ANAT_INPUT_DIR",  file.path(section_dir, "outputs", "tables", "r_input_files"))
output_dir <- Sys.getenv("ANAT_OUTPUT_DIR", file.path(section_dir, "outputs", "figures"))
scale_max  <- suppressWarnings(as.numeric(Sys.getenv("ANAT_SCALE_MAX", "")))
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

CLUSTERS <- c("C1", "C2", "C3")
# (display title, metric token, atlas type). Cortical volume = grayvol (dk).
PANELS <- list(
  list(title = "Cortical Thickness",  metric = "thickness", atlas = "dk"),
  list(title = "Surface Area",        metric = "area",      atlas = "dk"),
  list(title = "Cortical Volume",     metric = "grayvol",   atlas = "dk"),
  list(title = "Subcortical Volumes", metric = "volume",    atlas = "aseg")
)

# --- Helpers (shared with 13_combined_brain_figure.R) ------------------------
transform_region_column <- function(df) {
  df %>% plyr::mutate(label = case_when(
    endsWith(label, "_left")  ~ paste0("lh_", sub("_left$",  "", label)),
    endsWith(label, "_right") ~ paste0("rh_", sub("_right$", "", label)),
    TRUE ~ label))
}
rename_subcortical_labels <- function(df) {
  df %>% dplyr::mutate(label = case_when(
    label == "left-thalamus" ~ "Left-Thalamus-Proper",   label == "right-thalamus" ~ "Right-Thalamus-Proper",
    label == "left-caudate" ~ "Left-Caudate",            label == "right-caudate" ~ "Right-Caudate",
    label == "left-putamen" ~ "Left-Putamen",            label == "right-putamen" ~ "Right-Putamen",
    label == "left-pallidum" ~ "Left-Pallidum",          label == "right-pallidum" ~ "Right-Pallidum",
    label == "left-hippocampus" ~ "Left-Hippocampus",    label == "right-hippocampus" ~ "Right-Hippocampus",
    label == "left-amygdala" ~ "Left-Amygdala",          label == "right-amygdala" ~ "Right-Amygdala",
    label == "left-accumbens" ~ "Left-Accumbens",        label == "right-accumbens" ~ "Right-Accumbens",
    label == "left-ventraldc" ~ "Left-VentralDC",        label == "right-ventraldc" ~ "Right-VentralDC",
    label == "left-choroid-plexus" ~ "Left-Choroid-Plexus",
    label == "right-choroid-plexus" ~ "Right-Choroid-Plexus",
    label == "cc_posterior" ~ "CC_Posterior",            label == "cc_mid_posterior" ~ "CC_Mid_Posterior",
    label == "cc_central" ~ "CC_Central",                label == "cc_mid_anterior" ~ "CC_Mid_Anterior",
    label == "cc_anterior" ~ "CC_Anterior",
    label == "right-cerebellum-white-matter" ~ "Right-Cerebellum-White-Matter",
    label == "left-cerebellum-white-matter" ~ "Left-Cerebellum-White-Matter",
    label == "right-cerebellum-cortex" ~ "Right-Cerebellum-Cortex",
    label == "left-cerebellum-cortex" ~ "Left-Cerebellum-Cortex",
    label == "3rd-ventricle" ~ "x3rd-ventricle",         label == "4th-ventricle" ~ "x4th-ventricle",
    label == "brain-stem" ~ "Brain-Stem",
    label == "left-lateral-ventricle" ~ "Left-Lateral-Ventricle",
    label == "right-lateral-ventricle" ~ "Right-Lateral-Ventricle",
    TRUE ~ label))
}
load_metric <- function(cluster, metric, atlas_type) {
  fp <- file.path(input_dir, paste0("t_stat_cluster_", cluster, "_", metric, "_",
                                    atlas_type, "_mri_cluster_vs_td.csv"))
  if (!file.exists(fp)) return(NULL)
  has_header <- grepl("label", readLines(fp, n = 1), fixed = TRUE)
  d <- if (has_header) read.csv(fp, header = TRUE, stringsAsFactors = FALSE)
       else read.csv(fp, header = FALSE,
                     col.names = c("label", "t_stat", "p_val", "cohens_d", "p_fdr"),
                     stringsAsFactors = FALSE)
  if (atlas_type == "dk") { d <- transform_region_column(d); d$label <- tolower(d$label) }
  else                    { d$label <- tolower(d$label); d <- rename_subcortical_labels(d) }
  d$significant <- !is.na(d$p_fdr) & d$p_fdr < 0.05
  d
}

panel_theme <- theme_void(base_family = "Helvetica") +
  theme(plot.title    = element_text(size = 23, face = "plain", hjust = 0.5,
                                     margin = margin(b = 2, t = 6)),
        legend.position = "right",
        legend.key.height = unit(15, "pt"), legend.key.width = unit(34, "pt"),
        legend.title  = element_text(size = 11), legend.text = element_text(size = 9))

make_panel <- function(spec, cluster, fill_kind) {
  d <- load_metric(cluster, spec$metric, spec$atlas)
  if (is.null(d)) return(NULL)
  fill_col  <- if (fill_kind == "cohensd") "cohens_d" else "t_stat"
  fill_name <- if (fill_kind == "cohensd") "Cohen's d" else "t-statistic"
  d$fillval <- d[[fill_col]]
  lim <- if (fill_kind == "cohensd" && is.finite(scale_max)) scale_max
         else max(abs(d$fillval), na.rm = TRUE)
  d <- d %>% plyr::mutate(outline_color = ifelse(significant, "black", "grey70"),
                          outline_size  = ifelse(significant, 1.1, 0.15))
  if (spec$atlas == "dk") {
    p <- ggseg(d, mapping = aes(fill = fillval, col = outline_color, size = outline_size),
               atlas = "dk", position = "stacked")
  } else {
    p <- ggseg(d, mapping = aes(fill = fillval, colour = outline_color, size = outline_size),
               atlas = aseg)
  }
  p +
    scale_fill_continuous_divergingx(palette = "RdBu", mid = 0, rev = TRUE,
                                     limits = c(-lim, lim), oob = scales::squish,
                                     name = fill_name) +
    scale_colour_identity() + scale_size_identity() +
    panel_theme + labs(title = spec$title)
}

for (cl in CLUSTERS) {
  for (fk in c("tstat", "cohensd")) {
    panels <- Filter(Negate(is.null), lapply(PANELS, make_panel, cluster = cl, fill_kind = fk))
    if (length(panels) == 0) next
    stat_lbl <- if (fk == "cohensd") "Cohen's d" else "t-statistic"
    combined <- wrap_plots(panels, ncol = 1) +
      plot_annotation(title    = paste0(cl, " vs TD — ", stat_lbl),
                      subtitle = "FDR p < 0.05 outlined",
                      theme = theme(
                        plot.title    = element_text(size = 24, face = "plain", hjust = 0.5, family = "Helvetica"),
                        plot.subtitle = element_text(size = 12, hjust = 0.5, family = "Helvetica")))
    out_pdf <- file.path(output_dir, paste0("combined_cluster_brain_maps_", cl, "_fdr_", fk, ".pdf"))
    ggsave(out_pdf, plot = combined, width = 8, height = 5.2 * length(panels) + 1, limitsize = FALSE)
    cat("Saved:", out_pdf, "\n")
  }
}
