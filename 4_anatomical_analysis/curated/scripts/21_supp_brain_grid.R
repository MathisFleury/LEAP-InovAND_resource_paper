#!/usr/bin/env Rscript
# =============================================================================
# Supplementary brain grid: 4 features (rows) x 3 gene lists (cols).
# Columns: ChromEpiTF | SynGO | Protein coding.
# Rows   : Cortical Thickness | Surface Area | Cortical Volume | Subcortical.
# Feature names appear ONCE as left-side row labels; gene-list names ONCE as
# column headers; ONE shared colourbar at the bottom. Fill = OLS beta
# (metric ~ -log10 LOEUF); perm p<0.05 outlined in black.
# Output: outputs/figures_genetics_hg38_regperm/composite/supp_brain_all_features.pdf
# =============================================================================
suppressMessages({
  library(plyr); library(dplyr); library(ggplot2); library(grid)
  library(colorspace); library(ggseg); library(patchwork); library(scales)
})

args        <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
script_dir  <- if (length(script_path) == 0) getwd() else dirname(normalizePath(script_path))
section_dir <- dirname(script_dir)
out_base    <- Sys.getenv("LOEUF_OUT_DIR", "outputs")
input_dir   <- file.path(section_dir, out_base, "tables_genetics_hg38_regperm")
output_dir  <- file.path(section_dir, out_base, "figures_genetics_hg38_regperm", "composite")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ponytail: n hardcoded (stable regperm rosters); reads from table if changed.
PATHWAYS <- list(list(short = "chromepitf", title = "ChromEpiTF (n = 246)"),
                 list(short = "syngo",      title = "SynGO (n = 293)"),
                 list(short = "proteincoding", title = "Protein coding (n = 661)"))
PANELS <- list(list(row = "Cortical\nThickness", feat = "thickness",  atlas = "dk"),
               list(row = "Surface\nArea",        feat = "area",       atlas = "dk"),
               list(row = "Cortical\nVolume",     feat = "grayvol",    atlas = "dk"),
               list(row = "Subcortical\nVolumes", feat = "subcortical", atlas = "aseg"))
SCALE_MAX <- 1.0

rename_subcortical_labels <- function(df) df %>% dplyr::mutate(label = dplyr::case_when(
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

load_panel <- function(short, feat, atlas) {
  d <- read.csv(file.path(input_dir, sprintf("loeuf_reg_%s_%s.csv", short, feat)),
                stringsAsFactors = FALSE)
  if (atlas == "dk") {
    d$label <- tolower(ifelse(d$hemisphere == "left",
                              paste0("lh_", gsub(" ", "", d$region)),
                              paste0("rh_", gsub(" ", "", d$region))))
  } else {
    d$label <- tolower(d$mri_col); d <- rename_subcortical_labels(d)
  }
  d$significant <- !is.na(d$p_permutation) & d$p_permutation < 0.05
  d %>% plyr::mutate(outline_color = ifelse(significant, "black", "grey70"),
                     outline_size  = ifelse(significant, 1.1, 0.15)) %>%
    dplyr::select(label, beta, outline_color, outline_size)
}

make_map <- function(short, feat, atlas, title = NULL) {
  d <- load_panel(short, feat, atlas)
  p <- if (atlas == "dk")
    ggseg(d, mapping = aes(fill = beta, col = outline_color, size = outline_size),
          atlas = "dk", position = "stacked")
  else
    ggseg(d, mapping = aes(fill = beta, colour = outline_color, size = outline_size), atlas = aseg)
  p + scale_fill_continuous_divergingx(palette = "RdBu", mid = 0, rev = TRUE,
                                       limits = c(-SCALE_MAX, SCALE_MAX), oob = scales::squish,
                                       na.value = "grey90",
                                       name = "OLS beta  (metric ~ -log10 LOEUF)") +
    scale_colour_identity() + scale_size_identity() +
    theme_void(base_family = "Helvetica") +
    theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5, margin = margin(b = 2)),
          legend.title = element_text(size = 11), legend.text = element_text(size = 9),
          legend.key.height = unit(8, "pt"), legend.key.width = unit(70, "pt")) +
    labs(title = title)
}

featlab <- function(txt) wrap_elements(grid::textGrob(
  txt, rot = 90, gp = grid::gpar(fontsize = 14, fontface = "bold", fontfamily = "Helvetica")))

rows <- lapply(seq_along(PANELS), function(i) {
  s <- PANELS[[i]]
  maps <- lapply(PATHWAYS, function(pw)
    make_map(pw$short, s$feat, s$atlas, title = if (i == 1) pw$title else NULL))
  wrap_plots(c(list(featlab(s$row)), maps), nrow = 1, widths = c(0.10, 1, 1, 1))
})

combined <- wrap_plots(rows, ncol = 1) +
  plot_layout(guides = "collect") &
  guides(fill = guide_colourbar(direction = "horizontal", title.position = "top",
                                title.hjust = 0.5)) &
  theme(legend.position = "bottom")

out <- file.path(output_dir, "supp_brain_all_features.pdf")
ggsave(out, combined, units = "cm", width = 26, height = 31, limitsize = FALSE)
cat("Saved:", out, "\n")
