#!/usr/bin/env Rscript
# =============================================================================
# Clinical x MRI brain maps (ggseg) — autism (revised anat, curated clinical)
# =============================================================================
# Region-wise Pearson r between MRI measures and four clinical features, drawn
# as ggseg brain maps. 4 rows (modalities: CT | SA | cortical volume |
# subcortical volume) x 4 columns (SRS-2, RBS-R, FSIQ, VABS-II). Shared
# diverging r scale; BH-FDR q<0.05 outlined black.
#
# Consumes ../outputs/tables_clinical_mri/clinmri_<clin>_<feature>.csv written by
# 12_clinical_mri_correlations.py. Adapted from 6_loeuf_mri_brain_maps_hg38.R.
# Output: ../outputs/figures_clinical_mri/clinical_mri_brain_maps.pdf
# =============================================================================
suppressMessages({
  library(dplyr); library(ggplot2); library(colorspace)
  library(ggseg); library(patchwork); library(scales)
})

args <- commandArgs(trailingOnly = FALSE)
sp <- sub("--file=", "", args[grep("--file=", args)])
script_dir <- if (length(sp) == 0) getwd() else dirname(normalizePath(sp))
section_dir <- dirname(script_dir)
out_base <- Sys.getenv("LOEUF_OUT_DIR", "outputs")
input_dir <- file.path(section_dir, out_base, "tables_clinical_mri")
output_dir <- file.path(section_dir, out_base, "figures_clinical_mri")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

CLIN <- list(list(tok = "srs", title = "SRS-2 t-score"), list(tok = "rbsr", title = "RBS-R"),
             list(tok = "fsiq", title = "FSIQ"),          list(tok = "vabs", title = "VABS-II"))
MODS <- list(list(title = "CT",                 feat = "thickness",   atlas = "dk"),
             list(title = "SA",                 feat = "area",        atlas = "dk"),
             list(title = "Cortical volume",    feat = "grayvol",     atlas = "dk"),
             list(title = "Subcortical volume", feat = "subcortical", atlas = "aseg"))
SCALE_MAX <- 0.22   # fixed shared colour limit (+/-)

# aseg labels are CamelCase and mri_col already matches, except: Thalamus needs
# the "-Proper" suffix and the ventricles are prefixed "x" in ggseg.
aseg_label <- function(mri_col) {
  dplyr::recode(mri_col,
    "Left-Thalamus" = "Left-Thalamus-Proper", "Right-Thalamus" = "Right-Thalamus-Proper",
    "3rd-Ventricle" = "x3rd-ventricle", "4th-Ventricle" = "x4th-ventricle",
    .default = mri_col)
}

load_cell <- function(tok, feat, atlas) {
  fp <- file.path(input_dir, sprintf("clinmri_%s_%s.csv", tok, feat))
  if (!file.exists(fp)) return(NULL)
  d <- read.csv(fp, stringsAsFactors = FALSE)
  if (atlas == "dk")
    d$label <- tolower(ifelse(d$hemisphere == "left",
                              paste0("lh_", gsub(" ", "", d$region)),
                              paste0("rh_", gsub(" ", "", d$region))))
  else
    d$label <- aseg_label(d$mri_col)
  d$significant <- !is.na(d$p_fdr) & d$p_fdr < 0.05
  d
}

make_panel <- function(d, atlas, lim, title = NULL) {
  d <- d %>% mutate(outline_color = ifelse(significant, "black", "grey70"),
                    outline_size  = ifelse(significant, 1.0, 0.15)) %>%
    dplyr::select(label, correlation, outline_color, outline_size)
  p <- if (atlas == "dk")
    ggseg(d, mapping = aes(fill = correlation, col = outline_color, size = outline_size),
          atlas = "dk", position = "stacked")
  else
    ggseg(d, mapping = aes(fill = correlation, colour = outline_color, size = outline_size),
          atlas = aseg)
  p + scale_fill_continuous_divergingx(palette = "RdBu", mid = 0, rev = TRUE,
        limits = c(-lim, lim), oob = scales::squish, na.value = "grey90", name = "Pearson r",
        guide = guide_colourbar(direction = "horizontal", title.position = "left",
                                title.vjust = 0.9, barwidth = unit(9, "cm"),
                                barheight = unit(0.45, "cm"))) +
    scale_colour_identity() + scale_size_identity() +
    theme_void(base_family = "Helvetica") +
    theme(plot.title = element_text(size = 15, hjust = 0.5)) +
    labs(title = title)
}

cells <- list()
for (m in MODS) for (cl in CLIN)
  cells[[paste(m$feat, cl$tok, sep = "|")]] <- load_cell(cl$tok, m$feat, m$atlas)

lim <- if (!is.na(SCALE_MAX)) SCALE_MAX else
  max(0.15, max(vapply(cells, function(d) if (is.null(d)) NA_real_ else
    max(abs(d$correlation), na.rm = TRUE), numeric(1)), na.rm = TRUE))
cat(sprintf("Shared colour limit: +/-%.3f\n", lim))

row_label <- function(txt) ggplot() + annotate("text", 1, 1, label = txt, angle = 90,
  size = 5.5, fontface = "bold") + theme_void()
rows <- lapply(MODS, function(m) {
  maps <- lapply(CLIN, function(cl) {
    d <- cells[[paste(m$feat, cl$tok, sep = "|")]]
    ttl <- if (identical(m$title, MODS[[1]]$title)) cl$title else NULL
    if (is.null(d)) ggplot() + theme_void() else make_panel(d, m$atlas, lim, ttl)
  })
  wrap_plots(c(list(row_label(m$title)), maps), nrow = 1, widths = c(0.14, 1, 1, 1, 1))
})
combined <- wrap_plots(rows, ncol = 1) +
  plot_layout(guides = "collect") +          # ONE shared colourbar for the whole figure
  plot_annotation(
    title = "Cortical & subcortical MRI vs clinical features - autism",
    subtitle = "Region-wise Pearson r (fill); BH-FDR q<0.05 outlined black. CT: cortical thickness; SA: surface area.",
    theme = theme(plot.title = element_text(size = 20, hjust = 0.5, family = "Helvetica"),
                  plot.subtitle = element_text(size = 12, hjust = 0.5, family = "Helvetica")))
combined <- combined & theme(legend.position = "bottom")   # at the bottom

out <- file.path(output_dir, "clinical_mri_brain_maps.pdf")
ggsave(out, combined, width = 20, height = 15, limitsize = FALSE)
cat("Saved:", out, "\n")
