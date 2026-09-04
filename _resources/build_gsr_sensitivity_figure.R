#!/usr/bin/env Rscript
# =============================================================================
# GSR sensitivity figure -- functional Autism-vs-NT connectivity (R1.4 / R4.5a).
#
# DATA SOURCE: the manuscript's exact whole-group result (p.12: "393 autistic
# and 327 NT... hypoconnected component (122 edges, P_FWER=0.023)...
# hyperconnected component (139 edges, P_FWER=0.014)") is reproduced EXACTLY
# by 6_functional_analysis/concat/outputs/tables/nbs_wholegroup_nogsr_concat/
# nbs_summary.csv (na=393, nb=327, hyper=139 p=0.0138, hypo=122 p=0.0228) --
# NOT 6_functional_analysis/concat/outputs/nogsr_concat/ (a separate, lenient-QC,
# stale rebuild). nbs_wholegroup_gsr_concat/ is the matching GSR-toggle
# sensitivity run (same family, run-concatenated held fixed). Built from each
# dir's nbs_edges.csv (per-edge region/direction/t_stat list) -- no pre-built
# matrix exists there -- via the same atlas network mapping and count/
# proportion logic as 06_plot_double_network_matrix.R.
#
# 2 rows (hyper, hypo) x 2 conditions (Primary no-GSR, Sensitivity GSR), each
# condition a touching count|proportion pair (upper-left triangle = count,
# lower-right = proportion, matching Fig. 3d/f's convention and network
# order: rows run opposite to columns, Visual->Amyg.&Hippoc. top-to-bottom).
# Colour is always the full red/blue scale -- significance is stated in the
# subtitle (n edges, p_FWER), not colour-coded.
#
# Reads: 6_functional_analysis/concat/outputs/tables/nbs_wholegroup_{nogsr_concat,gsr_concat}/
#        nbs_edges.csv    (region [con_<src>/<tgt>], direction, t_stat)
#        nbs_summary.csv  (direction, comp_size_edges, p_fwer, na, nb)
# Output: _resources/Figure_GSR_Sensitivity.pdf
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(reshape2); library(ggplot2)
  library(grid); library(patchwork)
})

args <- commandArgs(trailingOnly = FALSE)
sp <- sub("--file=", "", args[grep("--file=", args)])
script_dir <- if (length(sp) == 0) getwd() else dirname(normalizePath(sp))
FUNC <- file.path(dirname(script_dir), "6_functional_analysis", "concat", "outputs", "tables")
ATLAS_FILE <- "/Users/mfleury/POSTDOC/LIBRAIRY/eeg_mri-pipeline/ressources/atlases/SCHAEFER/atlas-4S156Parcels/atlas-4S156Parcels_dseg.tsv"
OUT <- file.path(script_dir, "Figure_GSR_Sensitivity.pdf")

NET_ORDER <- c("Amyg. & Hippoc.", "Striatum", "Cerebellum", "Thalamus",
               "Default", "Cont", "Limbic", "SalVentAttn", "DorsAttn", "SomMot", "Vis")
READABLE <- c("Amyg. & Hippoc." = "Amyg.&Hip.", "Striatum" = "Striat.",
              "Cerebellum" = "Cereb.", "Thalamus" = "Thal.",
              "Default" = "DMN", "Cont" = "FPN", "Limbic" = "Limbic",
              "SalVentAttn" = "VAN", "DorsAttn" = "DAN", "SomMot" = "SomMot", "Vis" = "Vis")

atlas <- read_tsv(ATLAS_FILE, show_col_types = FALSE)
net_of <- setNames(atlas$network_label, atlas$label)
nets_all <- intersect(NET_ORDER, unique(atlas$network_label))
counts_per_net <- atlas %>% count(network_label, name = "n") %>% tibble::deframe()
total_mat <- outer(nets_all, nets_all, Vectorize(function(a, b) {
  na_ <- counts_per_net[[a]]; nb_ <- counts_per_net[[b]]
  if (a == b) na_ * (na_ - 1) / 2 else na_ * nb_
}))
dimnames(total_mat) <- list(nets_all, nets_all)

DIRECTIONS <- list(list(label = "Hypo",  direction = "hypo",  hue = "#003366"),
                    list(label = "Hyper", direction = "hyper", hue = "#B30000"))
CONDITIONS <- list(list(mode = "nogsr_concat", header = "PRIMARY (no GSR)"),
                    list(mode = "gsr_concat",  header = "SENSITIVITY (GSR)"))

featlab <- function(txt) wrap_elements(grid::textGrob(
  txt, rot = 90, gp = grid::gpar(fontsize = 9, fontface = "bold", fontfamily = "Helvetica")))

load_summary <- function(mode, direction) {
  fp <- file.path(FUNC, sprintf("nbs_wholegroup_%s", mode), "nbs_summary.csv")
  r <- read_csv(fp, show_col_types = FALSE)
  row <- r[r$direction == direction, ]
  list(n_edges = row$comp_size_edges[1], p_fwer = row$p_fwer[1])
}

# builds BOTH the count matrix and the proportion matrix for one (mode,
# direction) from the raw edge list (no pre-built matrix exists for this
# verified data source, unlike the per-condition analysis-suite outputs).
build_matrices <- function(mode, direction) {
  fp <- file.path(FUNC, sprintf("nbs_wholegroup_%s", mode), "nbs_edges.csv")
  edges <- read_csv(fp, show_col_types = FALSE) %>% filter(direction == !!direction)
  pairs <- strsplit(sub("^con_", "", edges$region), "/", fixed = TRUE)
  src <- vapply(pairs, `[`, character(1), 1); tgt <- vapply(pairs, `[`, character(1), 2)
  d <- data.frame(src = net_of[src], tgt = net_of[tgt]) %>%
    filter(src %in% nets_all, tgt %in% nets_all)
  cnt <- matrix(0, length(nets_all), length(nets_all), dimnames = list(nets_all, nets_all))
  for (k in seq_len(nrow(d))) {
    i <- d$src[k]; j <- d$tgt[k]
    cnt[i, j] <- cnt[i, j] + 1
    if (i != j) cnt[j, i] <- cnt[j, i] + 1
  }
  prp <- cnt / total_mat; prp[is.nan(prp)] <- 0
  list(count = cnt, proportion = prp)
}

sig_lab <- function(p) if (p < 0.05) sprintf("p_FWER=%.3f", p) else sprintf("n.s. (p=%.2f)", p)

# is_pair_start: TRUE for the count (left) panel of a Primary/Sensitivity
# pair, FALSE for its proportion (right) partner -- controls which side gets
# the near-zero margin (the two touch) vs. the small gap that separates one
# condition's pair from the next.
make_tri_panel <- function(mats, metric, hue, title_txt, subtitle_txt) {
  m <- mats[[metric]]
  nets <- intersect(NET_ORDER, rownames(m))
  m <- m[nets, nets]
  idx <- setNames(seq_along(nets), nets)
  labs <- READABLE[nets]; N <- length(nets)
  is_count <- metric == "count"

  long <- reshape2::melt(m, varnames = c("row", "col"), value.name = "value")
  long$xi <- idx[as.character(long$col)]
  long$yi <- idx[as.character(long$row)]
  # row order runs OPPOSITE to column order (Visual at the top row,
  # Amyg.&Hippoc. at the bottom; columns run Amyg.&Hippoc.->Visual
  # left-to-right), matching the reference figure. Both triangles INCLUDE the
  # diagonal (self-loop network pairs) -- confirmed against the reference:
  # its proportion panel shows a Visual-Visual cell, not just off-diagonal
  # pairs. count = upper-left+diag (yi>=xi); proportion = lower-right+diag
  # (yi<=xi) -- each cell just carries a different statistic in each panel.
  long <- if (is_count) long[long$yi >= long$xi, ] else long[long$yi <= long$xi, ]
  fmt <- if (is_count) "%.0f" else "%.2f"
  # near-zero margin on the touching side (within a pair); a small but
  # nonzero margin on the outer side (separates one condition's pair from
  # the next, or from the page edge).
  panel_margin <- if (is_count) margin(2, 0, 2, 3) else margin(2, 3, 2, 0)

  ggplot(long, aes(xi, yi, fill = value)) +
    geom_tile(color = "grey70", linewidth = 0.2) +
    scale_fill_gradient(low = "white", high = hue,
                        name = if (is_count) "No. of edges" else "Normalised\nproportion",
                        labels = function(v) sprintf(fmt, v),
                        guide = guide_colorbar(barwidth = unit(42, "pt"), barheight = unit(5, "pt"),
                                              title.position = "top", title.hjust = 0.5)) +
    scale_x_continuous(breaks = seq_len(N), labels = labs,
                       position = if (is_count) "top" else "bottom",
                       limits = c(0.5, N + 0.5), expand = c(0, 0)) +
    scale_y_continuous(breaks = seq_len(N), labels = labs,
                       position = if (is_count) "left" else "right",
                       limits = c(0.5, N + 0.5), expand = c(0, 0)) +
    coord_cartesian(clip = "off") + theme_minimal(base_family = "Helvetica") +
    theme(panel.grid = element_blank(), axis.title = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = if (is_count) 0 else 1,
                                     vjust = 0.5, size = 6),
          axis.text.y = element_text(size = 6),
          plot.title = element_text(size = 8.5, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 7, hjust = 0.5),
          plot.margin = panel_margin,
          legend.position = "bottom", legend.key.height = unit(5, "pt"),
          legend.key.width = unit(11, "pt"),
          legend.title = element_text(size = 6.5), legend.text = element_text(size = 5.5)) +
    labs(title = title_txt, subtitle = subtitle_txt)
}

rows <- lapply(DIRECTIONS, function(dir) {
  panels <- list()
  for (cond in CONDITIONS) {
    summ <- load_summary(cond$mode, dir$direction)
    mats <- build_matrices(cond$mode, dir$direction)
    panels[[length(panels) + 1]] <- make_tri_panel(
      mats, "count", dir$hue, title_txt = cond$header,
      subtitle_txt = sprintf("%d edges, %s", summ$n_edges, sig_lab(summ$p_fwer)))
    panels[[length(panels) + 1]] <- make_tri_panel(
      mats, "proportion", dir$hue, title_txt = NULL, subtitle_txt = NULL)
  }
  wrap_plots(c(list(featlab(dir$label)), panels), nrow = 1,
             widths = c(0.10, 1, 1, 1, 1))
})

combined <- wrap_plots(rows, ncol = 1) +
  plot_annotation(
    title = "GSR sensitivity: Autism vs NT functional connectivity",
    subtitle = "Both conditions run-concatenated (>6 min); GSR is the only difference between Primary and Sensitivity.",
    theme = theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(size = 8.5, hjust = 0.5)))

ggsave(OUT, combined, width = 18, height = 20, units = "cm", limitsize = FALSE)
cat("Saved:", OUT, "\n")
