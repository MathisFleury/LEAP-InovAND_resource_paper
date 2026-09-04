#!/usr/bin/env Rscript
# =============================================================================
# Age-binned Autism-vs-NT brain maps, ALL metrics in ONE figure: rows =
# metric (Cortical Thickness / Surface Area / Cortical Volume / Subcortical
# Volume), columns = age bin. Fill = Cohen's d (group effect, cohort[+sex]-
# adjusted), FDR<0.05 outlined; one shared colour scale + legend for the
# whole grid (pattern matches 4_anatomical_analysis's
# 21_supp_brain_grid.R). Bin/N labels appear once, on the top row.
#
# Reads the per-bin ggseg-input CSVs written by 10_age_bin_brain_maps.py:
#   outputs/tables/r_input_bin_<bin>/t_stat_anat_<atlas>_<metric>_interaction_groupdiff.csv
#   outputs/tables/r_input_bin_<bin>/group_counts.csv
# Output: outputs/figures/composite/age_bin_groupdiff_all_features.pdf
# =============================================================================
suppressPackageStartupMessages({
  library(plyr); library(dplyr); library(ggplot2); library(grid)
  library(colorspace); library(ggseg); library(patchwork)
})

args <- commandArgs(trailingOnly = FALSE)
sp <- sub("--file=", "", args[grep("--file=", args)])
script_dir <- if (length(sp) == 0) getwd() else dirname(normalizePath(sp))
fig_dir <- file.path(dirname(script_dir), "outputs", "figures")
tables_dir <- file.path(dirname(script_dir), "outputs", "tables")
output_dir <- file.path(fig_dir, "composite")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# order = chronological, left-to-right columns
# <6 kept separate from 6-9 (2026-08-31): merging created a within-bin age
# gap between groups (NT<6 is rare, n=6) that produced a false-positive
# signal -- see 10_age_bin_brain_maps.py. <6's low N shows in its own column
# header; read it before trusting that column.
BINS <- c("lt6" = "<6", "6_9" = "6-9", "10_13" = "10-13", "14_17" = "14-17",
          "18_24" = "18-24", "25plus" = "25+")
ROWS <- list(list(label = "Cortical\nThickness", atlas = "dk",   metric = "thickness"),
             list(label = "Surface\nArea",        atlas = "dk",   metric = "area"),
             list(label = "Cortical\nVolume",     atlas = "dk",   metric = "volume"),
             list(label = "Subcortical\nVolume",  atlas = "aseg", metric = "volume"))

dk_label <- function(df) df %>% plyr::mutate(label = tolower(case_when(
  endsWith(label, "_left")  ~ paste0("lh_", sub("_left$", "", label)),
  endsWith(label, "_right") ~ paste0("rh_", sub("_right$", "", label)),
  TRUE ~ label)))
aseg_label <- function(df) df %>% dplyr::mutate(label = dplyr::recode(label,
  "Left-Thalamus" = "Left-Thalamus-Proper", "Right-Thalamus" = "Right-Thalamus-Proper",
  "3rd-Ventricle" = "x3rd-ventricle", "4th-Ventricle" = "x4th-ventricle"))

load_bin <- function(bin_key, atlas_type, metric) {
  bindir <- file.path(tables_dir, sprintf("r_input_bin_%s", bin_key))
  fp <- file.path(bindir, sprintf("t_stat_anat_%s_%s_interaction_groupdiff.csv", atlas_type, metric))
  if (!file.exists(fp)) { warning(paste("missing", fp)); return(NULL) }
  d <- read.csv(fp, header = FALSE, col.names = c("label", "t_stat", "p_val", "cohens_d", "p_fdr"),
                stringsAsFactors = FALSE)
  d <- if (atlas_type == "dk") dk_label(d) else aseg_label(d)
  gc <- read.csv(file.path(bindir, "group_counts.csv"))
  d$col_title <- sprintf("%s\n(A=%d, NT=%d)", BINS[[bin_key]], gc$n_autism[1], gc$n_nt[1])
  d
}

make_panel <- function(d, atlas_obj, dmax, title = NULL) {
  atlas_labels <- unique(as.data.frame(atlas_obj)$label)
  dropped <- setdiff(d$label, atlas_labels)
  if (length(dropped)) cat("  no ggseg region for:", paste(dropped, collapse = ", "), "\n")
  d <- d %>% dplyr::filter(label %in% atlas_labels)
  sig <- d %>% dplyr::filter(p_fdr < 0.05)
  d <- d %>% plyr::mutate(
    outline_color = ifelse(label %in% sig$label, "black", "white"),
    outline_size  = ifelse(label %in% sig$label, 1.0, 0.15))
  ggseg(d, mapping = aes(fill = cohens_d, col = outline_color, size = outline_size),
        atlas = atlas_obj, position = "stacked") +
    scale_fill_continuous_divergingx(palette = "RdBu", mid = 0, rev = TRUE,
                                     limits = c(-dmax, dmax), name = "Cohen's d (Autism-NT)") +
    scale_colour_identity() + scale_size_identity() +
    theme_void(base_family = "Helvetica") +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
    labs(title = title)
}

featlab <- function(txt) wrap_elements(grid::textGrob(
  txt, rot = 90, gp = grid::gpar(fontsize = 12, fontface = "bold", fontfamily = "Helvetica")))

# load once: needed both for the shared colour scale and for building the rows
all_d <- lapply(ROWS, function(r) {
  ds <- lapply(names(BINS), load_bin, atlas_type = r$atlas, metric = r$metric)
  ds[!sapply(ds, is.null)]
})
dmax <- max(abs(unlist(lapply(all_d, function(ds)
  unlist(lapply(ds, function(d) range(d$cohens_d, na.rm = TRUE)))))))

rows <- lapply(seq_along(ROWS), function(i) {
  r <- ROWS[[i]]; ds <- all_d[[i]]
  if (length(ds) == 0) { cat("  skip", r$label, "(no input)\n"); return(NULL) }
  atlas_obj <- if (r$atlas == "dk") dk else aseg
  panels <- lapply(ds, function(d)
    make_panel(d, atlas_obj, dmax, title = if (i == 1) d$col_title[1] else NULL))
  wrap_plots(c(list(featlab(r$label)), panels), nrow = 1, widths = c(0.12, rep(1, length(panels))))
})
rows <- rows[!sapply(rows, is.null)]

combined <- wrap_plots(rows, ncol = 1) +
  plot_layout(guides = "collect") &
  guides(fill = guide_colourbar(direction = "horizontal", title.position = "top",
                                title.hjust = 0.5, barwidth = unit(140, "pt"))) &
  theme(legend.position = "bottom")
combined <- combined + plot_annotation(
  title = "Autism vs NT across age bins (Cohen's d, FDR<0.05 outlined)",
  theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)))

out <- file.path(output_dir, "age_bin_groupdiff_all_features.pdf")
ggsave(out, combined, units = "cm", width = 34, height = 26, limitsize = FALSE)
cat("Saved:", basename(out), "\n")
