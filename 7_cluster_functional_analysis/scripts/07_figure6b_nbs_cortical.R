# =============================================================================
# Figure 6b — functional NBS results (clusters vs NT), full composite
# =============================================================================
# Panel (b) analogue of figure_6, one page. Per direction block
# (Hypoconnectivity blue | Hyperconnectivity red) x 3 cluster columns (C1/C2/C3):
#   - cortical row : Schaefer-7 surface, fill = node-wise NBS-component edge count
#   - subcortical  : yabplot inset (08_*.png) embedded beneath each cluster
# Panel titles carry the component's NBS FWER p. One shared "No. of edges"
# colourbar per block (cortical); the subcortical insets share their own 1..max.
#
# Input : outputs/tables/nbs/cluster_<C>_nbs_edges.csv      (06_cluster_nbs.py)
#         outputs/tables/nbs/nbs_summary.csv                (FWER per cluster/dir)
#         outputs/figures/nbs/cluster_<C>_nbs_subcortical_<dir>.png  (08_*.py)
# Output: outputs/figures/nbs/figure6b_nbs.{pdf,png}
# =============================================================================
suppressMessages({
  library(plyr); library(dplyr); library(ggplot2)
  library(ggseg); library(ggsegSchaefer); library(patchwork); library(tidyr)
  library(magick); library(grid)
})

args        <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
script_dir  <- if (length(script_path) == 0) getwd() else dirname(normalizePath(script_path))
section_dir <- dirname(script_dir)
nbs_sub     <- Sys.getenv("NBS_SUBDIR", "nbs")     # per-variant output folder
tables_dir  <- file.path(section_dir, "outputs", "tables", nbs_sub)
nbs_dir     <- file.path(section_dir, "outputs", "figures", nbs_sub)
dir.create(nbs_dir, showWarnings = FALSE, recursive = TRUE)
ATLAS_FILE  <- "/Users/mfleury/POSTDOC/LIBRAIRY/eeg_mri-pipeline/ressources/atlases/SCHAEFER/atlas-4S156Parcels/atlas-4S156Parcels_dseg.tsv"

atlas <- read.table(ATLAS_FILE, sep = "\t", header = TRUE)
lab2net <- setNames(atlas$label_7network, atlas$label)   # node -> 7Networks_* parcel

# thresholding description for the caption (from the edge-mode marker in summary)
summ0 <- read.csv(file.path(tables_dir, "nbs_summary.csv"), stringsAsFactors = FALSE)
edge_mode <- if ("mode" %in% names(summ0)) summ0$mode[1] else "uncorrected"
thr_p     <- if ("thresh_p" %in% names(summ0)) summ0$thresh_p[1] else 0.05
CAPTION <- switch(edge_mode,
  nbs = sprintf("largest NBS component per direction (edge threshold p<%.2f)", thr_p),
  fdr = sprintf("edges at q<%.2f (BH-FDR)", thr_p),
  sprintf("edges at p<%.2f uncorrected", thr_p))

CLUSTERS <- c("C1", "C2", "C3")
DIRS <- list(list(key = "hypo",  title = "Hypoconnectivity",  low = "#deebf7", high = "#08519c"),
             list(key = "hyper", title = "Hyperconnectivity", low = "#fee0d2", high = "#a50f15"))

# node-wise edge count for one cluster x direction, mapped to Schaefer parcels
count_cortical <- function(cluster, direction) {
  fp <- file.path(tables_dir, sprintf("cluster_%s_nbs_edges.csv", cluster))
  if (!file.exists(fp)) return(NULL)
  d <- read.csv(fp, stringsAsFactors = FALSE)
  d <- d[d$direction == direction, , drop = FALSE]
  if (nrow(d) == 0) return(data.frame(region = character(), n_edges = integer()))
  ends <- sub("^con_", "", d$region)
  src <- sub("/.*$", "", ends); tgt <- sub("^.*/", "", ends)
  nodes <- c(src, tgt)
  parcels <- lab2net[nodes]                       # cortical Schaefer parcels
  # subcortical nodes carry label_7network == "n/a" (literal) — drop them; they
  # are shown by the yabplot subcortical insets, not on the cortical surface.
  parcels <- parcels[!is.na(parcels) & parcels != "n/a" & parcels != ""]
  if (length(parcels) == 0) return(data.frame(region = character(), n_edges = integer()))
  as.data.frame(table(region = parcels), stringsAsFactors = FALSE) |>
    dplyr::rename(n_edges = Freq)
}

# global max edge-count (shared scale across everything)
GMAX <- 0
grid <- list()
for (dl in DIRS) for (cl in CLUSTERS) {
  cc <- count_cortical(cl, dl$key)
  grid[[paste(dl$key, cl)]] <- cc
  if (!is.null(cc) && nrow(cc)) GMAX <- max(GMAX, cc$n_edges)
}
cat(sprintf("shared No. of edges limit = %d\n", GMAX))

make_panel <- function(cc, cluster, dl) {
  if (is.null(cc) || nrow(cc) == 0)
    cc <- data.frame(region = "7Networks_LH_Vis_1", n_edges = NA_real_)
  ttl <- paste("Cluster", cluster)
  p <- ggseg(cc, mapping = aes(fill = n_edges), atlas = schaefer7_100,
             position = "stacked", colour = "grey85", size = 0.1)
  p + scale_fill_gradient(low = dl$low, high = dl$high, na.value = "grey92",
                          limits = c(1, GMAX), name = "No. of edges") +
    theme_void(base_family = "Helvetica") +
    theme(plot.title = element_text(size = 12, hjust = 0.5, margin = margin(b = 1, t = 4)),
          legend.key.height = unit(30, "pt"), legend.key.width = unit(9, "pt"),
          legend.title = element_text(size = 9), legend.text = element_text(size = 8)) +
    labs(title = ttl)
}

# subcortical yabplot inset (PNG) -> trimmed rasterGrob wrapped as a plot
sub_inset <- function(cluster, direction) {
  fp <- file.path(nbs_dir, sprintf("cluster_%s_nbs_subcortical_%s.png", cluster, direction))
  if (!file.exists(fp)) return(patchwork::plot_spacer())
  img <- magick::image_trim(magick::image_read(fp))
  g <- grid::rasterGrob(as.raster(img), interpolate = TRUE)
  wrap_elements(full = g)
}

row_header <- function(txt) ggplot() +
  annotate("text", 0, 0.5, label = txt, hjust = 0, vjust = 0.5, size = 5.4,
           fontface = "bold", family = "Helvetica") +
  xlim(0, 1) + ylim(0, 1) + theme_void()

blocks <- lapply(DIRS, function(dl) {
  cort <- lapply(CLUSTERS, function(cl) make_panel(grid[[paste(dl$key, cl)]], cl, dl))
  crow <- wrap_plots(cort, nrow = 1) +
    plot_layout(guides = "collect") & theme(legend.position = "right")
  srow <- wrap_plots(lapply(CLUSTERS, sub_inset, direction = dl$key), nrow = 1)
  # header, cortical row (tall), subcortical inset row (short)
  wrap_plots(list(row_header(dl$title), crow, srow),
             ncol = 1, heights = c(0.07, 1, 0.42))
})

combined <- wrap_plots(blocks, ncol = 1) +
  plot_annotation(
    title = "Functional connectivity per cluster (vs neurotypicals)",
    subtitle = paste0("Node-wise edge count — ", CAPTION,
                      " (cortical surface + subcortical insets). ",
                      "Curated connectivity, k-means clusters."),
    theme = theme(plot.title    = element_text(size = 17, hjust = 0.5, family = "Helvetica"),
                  plot.subtitle = element_text(size = 10, hjust = 0.5, family = "Helvetica")))

out <- file.path(nbs_dir, "figure6b_nbs.pdf")
ggsave(out, combined, width = 13, height = 13, limitsize = FALSE)
ggsave(sub("\\.pdf$", ".png", out), combined, width = 13, height = 13, dpi = 200, limitsize = FALSE)
cat("Saved:", out, "\n")
