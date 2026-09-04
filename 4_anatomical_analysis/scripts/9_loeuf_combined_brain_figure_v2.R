# =============================================================================
# Combined single-page LOEUF brain figure (ggseg) — per pathway
# =============================================================================
# Template adapted from 07/13_combined_brain_figure.R (Autism vs NT) for the
# LOEUF regression + permutation outputs of
# 08_loeuf_mri_regression_permutation_v2.py.
#
# One figure per gene set (Protein coding / SynGO / ChromEpiTF); each is a
# single column of four stacked ggseg maps:
#   Cortical Thickness | Surface Area | Cortical Volume | Subcortical Volumes
# Fill = OLS beta (metric ~ -log10 LOEUF); ROIs with permutation p < 0.05 are
# outlined in black. Colour limit is a per-figure symmetric max|beta|.
#
# Input  : ../outputs/tables_genetics_hg38_regperm/loeuf_reg_<pathway>_<feat>.csv
# Output : ../outputs/figures_genetics_hg38_regperm/combined_loeuf_<pathway>_brain.pdf
# =============================================================================

suppressMessages({
  library(plyr); library(dplyr); library(ggplot2)
  library(colorspace); library(ggseg); library(patchwork); library(scales)
})

# --- Resolve paths -----------------------------------------------------------
args        <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
script_dir  <- if (length(script_path) == 0) getwd() else dirname(normalizePath(script_path))
section_dir <- dirname(script_dir)

# Output base dir (default "outputs"; override with LOEUF_OUT_DIR if needed).
out_base   <- Sys.getenv("LOEUF_OUT_DIR", "outputs")
input_dir  <- file.path(section_dir, out_base, "tables_genetics_hg38_regperm")
output_dir <- file.path(section_dir, out_base, "figures_genetics_hg38_regperm")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

PATHWAYS <- list(
  list(short = "proteincoding", title = "Protein coding"),
  list(short = "syngo",         title = "SynGO"),
  list(short = "chromepitf",    title = "ChromEpiTF")
)
# (display title, file feature token, atlas)
PANELS <- list(
  list(title = "Cortical Thickness",  feat = "thickness",  atlas = "dk"),
  list(title = "Surface Area",        feat = "area",       atlas = "dk"),
  list(title = "Cortical Volume",     feat = "grayvol",    atlas = "dk"),
  list(title = "Subcortical Volumes", feat = "subcortical", atlas = "aseg")
)

# --- Helpers (same conventions as 13_/06_) -----------------------------------
rename_subcortical_labels <- function(df) {
  df %>% dplyr::mutate(label = dplyr::case_when(
    label == "left-thalamus" ~ "Left-Thalamus-Proper",   label == "right-thalamus" ~ "Right-Thalamus-Proper",
    label == "left-caudate" ~ "Left-Caudate",            label == "right-caudate" ~ "Right-Caudate",
    label == "left-putamen" ~ "Left-Putamen",            label == "right-putamen" ~ "Right-Putamen",
    label == "left-pallidum" ~ "Left-Pallidum",          label == "right-pallidum" ~ "Right-Pallidum",
    label == "left-hippocampus" ~ "Left-Hippocampus",    label == "right-hippocampus" ~ "Right-Hippocampus",
    label == "left-amygdala" ~ "Left-Amygdala",          label == "right-amygdala" ~ "Right-Amygdala",
    label == "left-accumbens" ~ "Left-Accumbens",        label == "right-accumbens" ~ "Right-Accumbens",
    label == "left-ventraldc" ~ "Left-VentralDC",        label == "right-ventraldc" ~ "Right-VentralDC",
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

load_panel <- function(short, feat, atlas) {
  fp <- file.path(input_dir, sprintf("loeuf_reg_%s_%s.csv", short, feat))
  if (!file.exists(fp)) { warning(paste("missing", fp)); return(NULL) }
  d <- read.csv(fp, stringsAsFactors = FALSE)
  if (atlas == "dk") {
    d$label <- tolower(ifelse(d$hemisphere == "left",
                              paste0("lh_", gsub(" ", "", d$region)),
                              paste0("rh_", gsub(" ", "", d$region))))
  } else {
    d$label <- tolower(d$mri_col)
    d <- rename_subcortical_labels(d)
  }
  d$significant <- !is.na(d$p_permutation) & d$p_permutation < 0.05
  d
}

panel_theme <- theme_void(base_family = "Helvetica") +
  theme(plot.title    = element_text(size = 20, face = "plain", hjust = 0.5,
                                     margin = margin(b = 2, t = 6)),
        legend.position = "right",
        legend.key.height = unit(20, "pt"), legend.key.width = unit(22, "pt"),
        legend.title  = element_text(size = 11), legend.text = element_text(size = 9))

make_panel <- function(d, title, atlas, lim) {
  d <- d %>% plyr::mutate(
    outline_color = ifelse(significant, "black", "grey70"),
    outline_size  = ifelse(significant, 1.1, 0.15))
  # Keep ONLY the join key + aesthetics. ggseg joins on every shared column
  # name, so leftover `region`/`hemisphere`/`mri_col` cols clash with the atlas
  # (harmless for dk where region names align, but blanks the aseg map).
  d <- dplyr::select(d, label, beta, outline_color, outline_size)
  if (atlas == "dk") {
    p <- ggseg(d, mapping = aes(fill = beta, col = outline_color, size = outline_size),
               atlas = "dk", position = "stacked")
  } else {
    p <- ggseg(d, mapping = aes(fill = beta, colour = outline_color, size = outline_size),
               atlas = aseg)
  }
  p +
    scale_fill_continuous_divergingx(palette = "RdBu", mid = 0, rev = TRUE,
                                     limits = c(-lim, lim), oob = scales::squish,
                                     na.value = "grey90", name = "OLS beta") +
    scale_colour_identity() + scale_size_identity() +
    panel_theme + labs(title = title)
}

# --- Shared colour scale across ALL pathways/features ------------------------
# Pin to a number for a fixed scale, or NA for one global data-driven limit
# (so every figure/panel is directly comparable).
SCALE_MAX <- 1.0

cache <- lapply(PATHWAYS, function(pw)
  lapply(PANELS, function(s) load_panel(pw$short, s$feat, s$atlas)))
names(cache) <- vapply(PATHWAYS, function(pw) pw$short, character(1))

global_lim <- if (!is.na(SCALE_MAX)) SCALE_MAX else {
  vals <- unlist(lapply(cache, function(pd) lapply(pd, function(d)
    if (is.null(d)) NA_real_ else max(abs(d$beta), na.rm = TRUE))))
  max(0.05, max(vals, na.rm = TRUE))
}
cat(sprintf("Shared colour limit: +/-%.3f\n", global_lim))

# --- Build one combined figure per pathway (shared scale) --------------------
for (pw in PATHWAYS) {
  panel_data <- cache[[pw$short]]
  keep <- !vapply(panel_data, is.null, logical(1))
  if (!any(keep)) { cat(sprintf("No data for %s — skipping\n", pw$short)); next }

  n_carriers <- panel_data[keep][[1]]$n_total_carriers[1]

  panels <- Map(function(d, s) if (is.null(d)) NULL else make_panel(d, s$title, s$atlas, global_lim),
                panel_data, PANELS)
  panels <- Filter(Negate(is.null), panels)

  combined <- wrap_plots(panels, ncol = 1) +
    plot_annotation(
      title    = sprintf("-log10(LOEUF) x MRI  —  %s  (n = %s carriers)", pw$title, n_carriers),
      subtitle = "OLS beta fill · permutation p < 0.05 outlined in black",
      theme = theme(plot.title    = element_text(size = 22, hjust = 0.5, family = "Helvetica"),
                    plot.subtitle = element_text(size = 12, hjust = 0.5, family = "Helvetica")))

  out_pdf <- file.path(output_dir, sprintf("combined_loeuf_%s_brain.pdf", pw$short))
  # Width capped at 18 cm; height keeps the original 8:(5.2n+1) aspect ratio.
  ggsave(out_pdf, plot = combined, units = "cm",
         width = 18, height = 18 * (5.2 * length(panels) + 1) / 8,
         limitsize = FALSE)
  cat(sprintf("Saved: %s  (lim=%.3f)\n", basename(out_pdf), global_lim))
}

cat(sprintf("\nCombined LOEUF figures written to: %s\n", output_dir))
