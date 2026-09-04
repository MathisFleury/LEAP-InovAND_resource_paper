# =============================================================================
# Figure 6a — anatomical results grid (clusters vs NT), ggseg t-stat maps
# =============================================================================
# Reproduces panel (a) of the Nature Neuroscience figure_6 layout:
#   3 modality rows  x  3 cluster columns
#     rows : Cortical Thickness | Surface Area | Subcortical Volumes
#     cols : Cluster C1 | Cluster C2 | Cluster C3
# Fill = Cohen's d (cluster autism vs pooled NT), diverging RdBu centred at 0,
# FDR q<0.05 ROIs outlined black. One shared "Cohen's d" colourbar per row (its
# limit = max|d| across the three clusters for that modality).
#
# Input : <CLUSTER_OUT_DIR>/tables/r_input_files/t_stat_cluster_<C>_<metric>_<atlas>_mri_cluster_vs_td.csv
#         (CLUSTER_OUT_DIR default: "outputs" -- the staging dir the curated
#         wrapper stages into before renaming to outputs_curated_<method>/;
#         must match the other steps' default so this step's output survives
#         that rename instead of being written straight into a dir the
#         wrapper later deletes)
# Output: <CLUSTER_OUT_DIR>/figures/figure6a_anatomical_clusters.{pdf,png}
#
# Env overrides: ANAT_INPUT_DIR, ANAT_OUTPUT_DIR (take precedence over
# CLUSTER_OUT_DIR when set, e.g. for a standalone run against an
# already-renamed outputs/).
# =============================================================================
suppressMessages({
  library(plyr); library(dplyr); library(ggplot2)
  library(colorspace); library(ggseg); library(patchwork); library(scales)
})

args        <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
script_dir  <- if (length(script_path) == 0) getwd() else dirname(normalizePath(script_path))
section_dir <- dirname(script_dir)                                  # 5_cluster_anatomical_analysis/

out_base           <- file.path(section_dir, Sys.getenv("CLUSTER_OUT_DIR", "outputs"))
default_fig_dir    <- file.path(out_base, "figures")
default_tables_dir <- file.path(out_base, "tables")
input_dir  <- Sys.getenv("ANAT_INPUT_DIR",  file.path(default_tables_dir, "r_input_files"))
output_dir <- Sys.getenv("ANAT_OUTPUT_DIR", default_fig_dir)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

CLUSTERS <- c("C1", "C2", "C3")
ROWS <- list(
  list(title = "Cortical Thickness",  metric = "thickness", atlas = "dk"),
  list(title = "Surface Area",        metric = "area",      atlas = "dk"),
  list(title = "Subcortical Volumes", metric = "volume",    atlas = "aseg")
)

# --- label mapping (shared with 04_combined_cluster_brain_figure.R) ----------
transform_region_column <- function(df) df %>% plyr::mutate(label = case_when(
  endsWith(label, "_left")  ~ paste0("lh_", sub("_left$",  "", label)),
  endsWith(label, "_right") ~ paste0("rh_", sub("_right$", "", label)),
  TRUE ~ label))
rename_subcortical_labels <- function(df) df %>% dplyr::mutate(label = case_when(
  label == "left-thalamus" ~ "Left-Thalamus-Proper",   label == "right-thalamus" ~ "Right-Thalamus-Proper",
  label == "left-caudate" ~ "Left-Caudate",            label == "right-caudate" ~ "Right-Caudate",
  label == "left-putamen" ~ "Left-Putamen",            label == "right-putamen" ~ "Right-Putamen",
  label == "left-pallidum" ~ "Left-Pallidum",          label == "right-pallidum" ~ "Right-Pallidum",
  label == "left-hippocampus" ~ "Left-Hippocampus",    label == "right-hippocampus" ~ "Right-Hippocampus",
  label == "left-amygdala" ~ "Left-Amygdala",          label == "right-amygdala" ~ "Right-Amygdala",
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

load_metric <- function(cluster, metric, atlas) {
  fp <- file.path(input_dir, paste0("t_stat_cluster_", cluster, "_", metric, "_",
                                    atlas, "_mri_cluster_vs_td.csv"))
  if (!file.exists(fp)) return(NULL)
  d <- read.csv(fp, stringsAsFactors = FALSE)
  if (atlas == "dk") { d <- transform_region_column(d); d$label <- tolower(d$label) }
  else               { d$label <- tolower(d$label); d <- rename_subcortical_labels(d) }
  d$significant <- !is.na(d$p_fdr) & d$p_fdr < 0.05
  d
}

# --- panel + row-label builders ----------------------------------------------
panel_theme <- theme_void(base_family = "Helvetica") +
  theme(plot.title = element_text(size = 13, hjust = 0.5, margin = margin(b = 1, t = 4)),
        legend.key.height = unit(34, "pt"), legend.key.width = unit(10, "pt"),
        legend.title = element_text(size = 10), legend.text = element_text(size = 8))

make_panel <- function(spec, cluster, lim) {
  d <- load_metric(cluster, spec$metric, spec$atlas)
  if (is.null(d)) return(patchwork::plot_spacer())
  d <- d %>% plyr::mutate(outline_color = ifelse(significant, "black", "grey75"),
                          outline_size  = ifelse(significant, 1.1, 0.12))
  p <- if (spec$atlas == "dk")
    ggseg(d, mapping = aes(fill = cohens_d, col = outline_color, size = outline_size),
          atlas = "dk", position = "stacked")
  else
    ggseg(d, mapping = aes(fill = cohens_d, colour = outline_color, size = outline_size),
          atlas = aseg)
  p + scale_fill_continuous_divergingx(palette = "RdBu", mid = 0, rev = TRUE,
        limits = c(-lim, lim), oob = scales::squish, name = "Cohen's d") +
    scale_colour_identity() + scale_size_identity() +
    panel_theme + labs(title = paste("Cluster", cluster))
}

# Horizontal modality header, left-aligned, sitting above each row of brains.
row_header <- function(txt) ggplot() +
  annotate("text", 0, 0.5, label = txt, hjust = 0, vjust = 0.5, size = 5.4,
           fontface = "bold", family = "Helvetica") +
  xlim(0, 1) + ylim(0, 1) + theme_void()

# --- assemble ----------------------------------------------------------------
# ONE shared Cohen's d scale across every panel (all modalities, all clusters).
GLIM <- max(unlist(lapply(ROWS, function(spec)
  lapply(CLUSTERS, function(cl) {
    d <- load_metric(cl, spec$metric, spec$atlas)
    if (is.null(d)) NA_real_ else max(abs(d$cohens_d), na.rm = TRUE)
  }))), na.rm = TRUE)
GLIM <- ceiling(GLIM * 10) / 10                     # round up to nearest 0.1
cat(sprintf("shared Cohen's d limit (all panels) = +/-%.1f\n", GLIM))

# Identical scale (GLIM) everywhere; one colourbar per row; horizontal header
# above each modality's row of brains.
blocks <- lapply(ROWS, function(spec) {
  panels <- lapply(CLUSTERS, function(cl) make_panel(spec, cl, GLIM))
  prow <- wrap_plots(panels, nrow = 1) +
    plot_layout(guides = "collect") & theme(legend.position = "right")
  wrap_plots(list(row_header(spec$title), prow), ncol = 1, heights = c(0.08, 1))
})

combined <- wrap_plots(blocks, ncol = 1) +
  plot_annotation(
    title = "Anatomical alterations per cluster (vs neurotypicals)",
    subtitle = "Region-wise Cohen's d (fill); FDR q<0.05 outlined black. Full-sample-fit, curated k-means.",
    theme = theme(plot.title    = element_text(size = 17, hjust = 0.5, family = "Helvetica"),
                  plot.subtitle = element_text(size = 10, hjust = 0.5, family = "Helvetica")))

out_pdf <- file.path(output_dir, "figure6a_anatomical_clusters.pdf")
ggsave(out_pdf, combined, width = 13, height = 11, limitsize = FALSE)
ggsave(sub("\\.pdf$", ".png", out_pdf), combined, width = 13, height = 11, dpi = 200, limitsize = FALSE)
cat("Saved:", out_pdf, "\n")
