#!/usr/bin/env Rscript
# =============================================================================
# ggseg brain maps of the diagnosis x age / x sex interaction t-statistics.
# Adapts 4_anatomical_analysis/scripts/03_plot_anatomical_mri_brain_visualizations.R
# Reads the per-metric interaction CSVs written by 01_age_sex_interactions.py
#   outputs/figures/r_input_files/t_stat_anat_<atlas>_<metric>_interaction_<model>.csv
# One ggseg map per (model, metric); FDR<0.05 regions outlined (none expected).
# =============================================================================
suppressPackageStartupMessages({
  library(plyr); library(dplyr); library(ggplot2); library(colorspace); library(ggseg)
})

brain_map_theme <- theme_minimal() +
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 13, hjust = 0.5),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10), legend.position = "bottom",
        panel.grid = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(),
        strip.text = element_text(size = 11, face = "bold"),
        strip.background = element_rect(fill = "lightgray", color = "black"))

args <- commandArgs(trailingOnly = FALSE)
sp <- sub("--file=", "", args[grep("--file=", args)])
script_dir <- if (length(sp) == 0) getwd() else dirname(normalizePath(sp))
# R_INPUT_DIR / OUT_TAG let a caller point at a subset's inputs and tag outputs
# (e.g. the age-5-22 / age-matched sensitivity subsets).
input_dir <- Sys.getenv("R_INPUT_DIR",
  file.path(dirname(script_dir), "outputs", "figures", "r_input_files"))
output_dir <- file.path(dirname(script_dir), "outputs", "figures")
out_tag <- Sys.getenv("OUT_TAG", "")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

n_suffix <- ""
cp <- file.path(input_dir, "group_counts.csv")
if (file.exists(cp)) { gc <- read.csv(cp); n_suffix <- sprintf("  |  Autism n=%s, NT n=%s", gc$n_autism[1], gc$n_nt[1]) }

dk_label <- function(df) df %>% plyr::mutate(label = tolower(case_when(
  endsWith(label, "_left")  ~ paste0("lh_", sub("_left$", "", label)),
  endsWith(label, "_right") ~ paste0("rh_", sub("_right$", "", label)),
  TRUE ~ label)))

# input aseg labels are already FreeSurfer-cased; only 4 differ from ggseg's aseg
aseg_label <- function(df) df %>% dplyr::mutate(label = dplyr::recode(label,
  "Left-Thalamus" = "Left-Thalamus-Proper", "Right-Thalamus" = "Right-Thalamus-Proper",
  "3rd-Ventricle" = "x3rd-ventricle", "4th-Ventricle" = "x4th-ventricle"))

load_one <- function(model, metric, atlas_type) {
  fp <- file.path(input_dir, sprintf("t_stat_anat_%s_%s_interaction_%s.csv", atlas_type, metric, model))
  if (!file.exists(fp)) { warning(paste("missing", fp)); return(NULL) }
  d <- read.csv(fp, header = FALSE, col.names = c("label", "t_stat", "p_val", "cohens_d", "p_fdr"),
                stringsAsFactors = FALSE)
  d$significant <- d$p_fdr < 0.05
  d
}

plot_map <- function(model, metric, atlas_type, atlas_obj, pretty) {
  d <- load_one(model, metric, atlas_type); if (is.null(d)) return(invisible())
  d <- if (atlas_type == "dk") dk_label(d) else aseg_label(d)
  # ggseg's aseg geometry has only the right cerebellum; drop labels with no
  # plottable region so the map renders without "not merged" gaps.
  atlas_labels <- unique(as.data.frame(atlas_obj)$label)
  dropped <- setdiff(d$label, atlas_labels)
  if (length(dropped)) cat("  no ggseg region for:", paste(dropped, collapse = ", "), "\n")
  d <- d %>% dplyr::filter(label %in% atlas_labels)
  tmax <- max(abs(range(d$t_stat, na.rm = TRUE)))
  sig <- d %>% dplyr::filter(p_fdr < 0.05)
  d <- d %>% plyr::mutate(
    outline_color = ifelse(label %in% sig$label, "black", "white"),
    outline_size  = ifelse(label %in% sig$label, 1.2, 0.2))
  # "groupdiff" (e.g. per-age-bin Autism-vs-NT maps) reads differently from an
  # interaction model; everything else keeps the original interaction wording.
  is_groupdiff <- model == "groupdiff"
  stat_name <- if (is_groupdiff) "group (Autism-NT) t" else "interaction t"
  title_txt <- if (is_groupdiff) sprintf("%s: Autism vs NT", pretty)
               else sprintf("%s: diagnosis x %s interaction", pretty, model)
  subtitle_txt <- paste0(if (is_groupdiff) "group t-statistic" else "interaction t-statistic",
                         " (FDR<0.05 outlined)", n_suffix)
  p <- d %>% ggseg(mapping = aes(fill = t_stat, col = outline_color, size = outline_size),
                   atlas = atlas_obj, position = "stacked") +
    scale_fill_continuous_divergingx(palette = "RdBu", mid = 0, rev = TRUE,
                                     limits = c(-tmax, tmax), name = stat_name) +
    scale_colour_identity() + scale_size_identity() + brain_map_theme +
    labs(title = title_txt, subtitle = subtitle_txt)
  out <- file.path(output_dir, sprintf("anat_interaction_%s_%s_%s%s.pdf", atlas_type, metric, model, out_tag))
  ggsave(out, p, width = 10, height = 10)
  write.csv(d, sub("\\.pdf$", ".csv", out), row.names = FALSE)
  cat("Saved:", basename(out), "\n")
}

specs <- list(c("thickness", "dk"), c("area", "dk"), c("volume", "dk"), c("volume", "aseg"))
pretty <- c(thickness = "Cortical Thickness", area = "Surface Area", volume = "Volume")
# Detect which model(s) this input_dir actually has files for (age/sex
# interaction, or groupdiff) instead of hardcoding -- avoids noisy missing-file
# warnings when a caller (e.g. per-age-bin maps) only ever writes one model.
avail <- list.files(input_dir, pattern = "^t_stat_anat_.*_interaction_.+\\.csv$")
models <- unique(sub(".*_interaction_(.+)\\.csv$", "\\1", avail))
if (length(models) == 0) models <- c("age", "sex")
for (model in models) for (s in specs) {
  metric <- s[1]; atlas_type <- s[2]
  atlas_obj <- if (atlas_type == "dk") dk else aseg
  lab <- if (atlas_type == "aseg") "Subcortical Volume" else pretty[[metric]]
  plot_map(model, metric, atlas_type, atlas_obj, lab)
}
cat("Done.\n")
