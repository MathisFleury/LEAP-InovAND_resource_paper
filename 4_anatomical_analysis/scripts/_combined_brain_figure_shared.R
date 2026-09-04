# =============================================================================
# Combined single-page brain figure (ggseg) — Autism vs control
# =============================================================================
# One column, four stacked FDR-corrected ggseg maps:
#   Cortical Thickness | Surface Area | Cortical Volume | Subcortical Volumes
# Each panel: divergingx RdBu fill, FDR-significant ROIs outlined in black and
# labelled with a leader line (ggrepel). A per-feature title sits above each map.
#
# Parameterised via environment variables so the same script serves the main
# (t-statistic) and the /curated (Cohen's d) pipelines:
#   ANAT_INPUT_DIR   input r_input_files dir   (default: main outputs)
#   ANAT_OUTPUT_DIR  output dir                (default: main outputs)
#   ANAT_FILL        "tstat" | "cohensd"       (default: tstat)
#   ANAT_SCALE_MAX   fixed symmetric limit, e.g. 0.4 (default: per-panel max|x|)
#   ANAT_COMPARISON  label for the page header (default: "Autism vs TD")
# The cortical-volume metric file ("volume" or "grayvol") is auto-detected.
# =============================================================================

suppressMessages({
  library(plyr); library(dplyr); library(ggplot2)
  library(colorspace); library(ggseg); library(ggrepel)
  library(patchwork); library(scales)
})

# --- Resolve paths / params --------------------------------------------------
args        <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
script_dir  <- if (length(script_path) == 0) getwd() else dirname(normalizePath(script_path))
section_dir <- dirname(script_dir)

input_dir  <- Sys.getenv("ANAT_INPUT_DIR",  file.path(section_dir, "outputs", "tables", "r_input_files"))
output_dir <- Sys.getenv("ANAT_OUTPUT_DIR", file.path(section_dir, "outputs", "figures"))
# ANAT_SCALE_MAX (if set) gives a fixed symmetric limit for the Cohen's d figure
# only; t-stats always use a per-panel data-driven limit.
scale_max  <- suppressWarnings(as.numeric(Sys.getenv("ANAT_SCALE_MAX", "")))
comparison <- Sys.getenv("ANAT_COMPARISON", "Autism vs TD")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Sample sizes for the page header. The curated-pipeline CSVs carry n_a (autism / group A)
# and n_b (NT) columns; the main headerless pipeline does not, so n_suffix stays "".
n_suffix <- ""
for (mm in list(c("thickness", "dk"), c("area", "dk"), c("volume", "aseg"))) {
  fp <- file.path(input_dir, paste0("t_stat_anat_", mm[2], "_", mm[1], "_mri_autism_vs_control.csv"))
  if (file.exists(fp) && grepl("label", readLines(fp, n = 1), fixed = TRUE)) {
    d0 <- read.csv(fp, header = TRUE, stringsAsFactors = FALSE)
    if (all(c("n_a", "n_b") %in% names(d0))) {
      n_suffix <- sprintf("  (Autism n = %s, NT n = %s)", d0$n_a[1], d0$n_b[1]); break
    }
  }
}

# Auto-detect the cortical-volume file metric token (main: volume, curated: grayvol)
cortvol_metric <- NA
for (m in c("volume", "grayvol")) {
  if (file.exists(file.path(input_dir, paste0("t_stat_anat_dk_", m, "_mri_autism_vs_control.csv")))) {
    cortvol_metric <- m; break
  }
}

# (display title, file metric token, atlas type)
PANELS <- list(
  list(title = "Cortical Thickness",  metric = "thickness",      atlas = "dk"),
  list(title = "Surface Area",        metric = "area",           atlas = "dk"),
  list(title = "Cortical Volume",     metric = cortvol_metric,   atlas = "dk"),
  list(title = "Subcortical Volumes", metric = "volume",         atlas = "aseg")
)

# --- Helpers -----------------------------------------------------------------
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
load_metric <- function(metric, atlas_type) {
  if (is.na(metric)) return(NULL)
  fp <- file.path(input_dir, paste0("t_stat_anat_", atlas_type, "_", metric, "_mri_autism_vs_control.csv"))
  if (!file.exists(fp)) { warning(paste("missing", fp)); return(NULL) }
  # main pipeline: headerless 5-col; curated: header + extra CI/n columns.
  has_header <- grepl("label", readLines(fp, n = 1), fixed = TRUE)
  if (has_header) {
    d <- read.csv(fp, header = TRUE, stringsAsFactors = FALSE)
  } else {
    d <- read.csv(fp, header = FALSE,
                  col.names = c("label", "t_stat", "p_val", "cohens_d", "p_fdr"),
                  stringsAsFactors = FALSE)
  }
  if (atlas_type == "dk") {
    d <- transform_region_column(d); d$label <- tolower(d$label)
  } else {
    d$label <- tolower(d$label); d <- rename_subcortical_labels(d)
  }
  # Ensure a Bonferroni column. Revision CSVs already carry p_bonf (per metric);
  # the main headerless form only has p_val, so derive it there.
  if (!"p_bonf" %in% names(d) && "p_val" %in% names(d)) {
    d$p_bonf <- p.adjust(d$p_val, method = "bonferroni")
  }
  d
}

panel_theme <- theme_void(base_family = "Helvetica") +
  theme(plot.title    = element_text(size = 23, face = "plain", hjust = 0.5,
                                     margin = margin(b = 2, t = 6)),
        legend.position = "right",
        legend.key.height = unit(20, "pt"), legend.key.width = unit(25, "pt"),
        legend.title  = element_text(size = 11), legend.text = element_text(size = 9))

make_panel <- function(spec, fill_kind, sig_col) {
  d <- load_metric(spec$metric, spec$atlas)
  if (is.null(d)) return(NULL)
  fill_col  <- if (fill_kind == "cohensd") "cohens_d" else "t_stat"
  fill_name <- if (fill_kind == "cohensd") "Cohen's d" else "t-statistic"
  d$fillval <- d[[fill_col]]
  d$significant <- sig_col %in% names(d) & !is.na(d[[sig_col]]) & d[[sig_col]] < 0.05
  # fixed limit only for Cohen's d (when ANAT_SCALE_MAX set); else per-panel max
  lim <- if (fill_kind == "cohensd" && is.finite(scale_max)) scale_max
         else max(abs(d$fillval), na.rm = TRUE)

  d <- d %>% plyr::mutate(
    outline_color = ifelse(significant, "black", "white"),
    outline_size  = ifelse(significant, 0.5, 0.15))

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

# Emit one figure per (multiple-comparison method × fill statistic):
#   FDR + Bonferroni  ×  t-statistic + Cohen's d.
SIG_METHODS <- list(
  list(key = "fdr",        col = "p_fdr",  label = "FDR"),
  list(key = "bonferroni", col = "p_bonf", label = "Bonferroni")
)
for (sm in SIG_METHODS) {
  for (fk in c("tstat", "cohensd")) {
    panels <- Filter(Negate(is.null),
                     lapply(PANELS, make_panel, fill_kind = fk, sig_col = sm$col))
    if (length(panels) == 0) { warning("No panels for ", sm$key, " ", fk); next }
    stat_lbl <- if (fk == "cohensd") "Cohen's d" else "t-statistic"
    combined <- wrap_plots(panels, ncol = 1) +
      plot_annotation(title    = paste0(comparison, " — ", stat_lbl),
                      subtitle = paste0(sm$label, " p < 0.05 outlined", n_suffix),
                      theme = theme(
                        plot.title    = element_text(size = 24, face = "plain", hjust = 0.5, family = "Helvetica"),
                        plot.subtitle = element_text(size = 12, hjust = 0.5, family = "Helvetica")))
    out_pdf <- file.path(output_dir, paste0("combined_brain_maps_", sm$key, "_", fk, ".pdf"))
    # Width capped at 18 cm; height keeps the original 8:(5.2n+1) aspect ratio.
    ggsave(out_pdf, plot = combined, units = "cm",
           width = 18, height = 18 * (5.2 * length(panels) + 1) / 8,
           limitsize = FALSE)
    cat("Saved:", out_pdf, "\n")
  }
}
