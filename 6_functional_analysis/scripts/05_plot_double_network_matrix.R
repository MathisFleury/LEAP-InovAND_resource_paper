#!/usr/bin/env Rscript
# =============================================================================
# Combined double network matrix (Autism vs NT), hyper & hypo.
# =============================================================================
# Single figure per connectivity type:
#   upper-left triangle (+ diagonal) = raw edge COUNT     (left colorbar)
#   lower-right triangle             = normalised PROPORTION (right colorbar)
#       proportion = count / total possible network-pair connections
# Reproduces eeg_mri-pipeline .../autism_vs_td_<type>_double_matrix_final.pdf
# but renders the combined figure directly in R (the original emitted the two
# triangles as separate PDFs assembled by hand).
#
# Improvements over the original:
#   - one combined figure via ggnewscale (two fill scales / colorbars)
#   - square tiles (coord_fixed), Helvetica, axis labels on all four sides
#   - parameterised input/output via AUTISM_TD_OUTPUT_DIR (defaults to revised/)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(reshape2)
  library(ggplot2); library(ggnewscale); library(cowplot)
})

ATLAS_FILE <- "/Users/mfleury/POSTDOC/LIBRAIRY/eeg_mri-pipeline/ressources/atlases/SCHAEFER/atlas-4S156Parcels/atlas-4S156Parcels_dseg.tsv"

args <- commandArgs(trailingOnly = FALSE)
sp <- sub("--file=", "", args[grep("--file=", args)])
script_dir <- if (length(sp) == 0) getwd() else dirname(normalizePath(sp))
section_dir <- dirname(script_dir)
IO_DIR <- Sys.getenv("AUTISM_TD_OUTPUT_DIR",
                     file.path(section_dir, "outputs", "figures", "revised"))

NET_ORDER <- c("Amyg. & Hippoc.", "Striatum", "Cerebellum", "Thalamus",
               "Default", "Cont", "Limbic", "SalVentAttn", "DorsAttn",
               "SomMot", "Vis")
READABLE <- c("Amyg. & Hippoc." = "Amyg. & Hippoc.", "Striatum" = "Striatum",
              "Cerebellum" = "Cerebellum", "Thalamus" = "Thalamus",
              "Default" = "Default Mode", "Cont" = "Frontoparietal",
              "Limbic" = "Limbic", "SalVentAttn" = "Ventral Attention",
              "DorsAttn" = "Dorsal Attention", "SomMot" = "Somatomotor", "Vis" = "Visual")

atlas <- read_tsv(ATLAS_FILE, show_col_types = FALSE)

# total possible connections per network pair (denominator for proportion)
counts_per_net <- atlas %>% count(network_label, name = "n") %>% tibble::deframe()
nets <- intersect(NET_ORDER, names(counts_per_net))
total_mat <- outer(nets, nets, Vectorize(function(a, b) {
  na <- counts_per_net[[a]]; nb <- counts_per_net[[b]]
  if (a == b) na * (na - 1) / 2 else na * nb
}))
dimnames(total_mat) <- list(nets, nets)

build_count_matrix <- function(df) {
  d <- df %>%
    left_join(atlas %>% select(label, network_label), by = c("source" = "label")) %>%
    rename(src = network_label) %>%
    left_join(atlas %>% select(label, network_label), by = c("target" = "label")) %>%
    rename(tgt = network_label) %>%
    filter(src %in% nets, tgt %in% nets)
  m <- matrix(0, length(nets), length(nets), dimnames = list(nets, nets))
  for (k in seq_len(nrow(d))) {
    i <- d$src[k]; j <- d$tgt[k]
    m[i, j] <- m[i, j] + 1
    if (i != j) m[j, i] <- m[j, i] + 1
  }
  m
}

plot_double <- function(conn_type) {
  fp <- file.path(IO_DIR, sprintf("autism_vs_td_%sconnectivity_full_data.csv", conn_type))
  if (!file.exists(fp)) { cat("Missing:", fp, "\n"); return(invisible()) }
  df <- read_csv(fp, show_col_types = FALSE)
  if (nrow(df) == 0) { cat("Empty:", basename(fp), "\n"); return(invisible()) }

  cnt <- build_count_matrix(df)
  prop <- cnt / total_mat; prop[is.nan(prop)] <- 0

  idx <- setNames(seq_along(nets), nets)          # Amyg=1 ... Visual=N
  long <- reshape2::melt(cnt, varnames = c("row", "col"), value.name = "count")
  long$prop <- reshape2::melt(prop)$value
  long$xi <- idx[as.character(long$col)]          # x: col, Amyg left
  long$yi <- idx[as.character(long$row)]          # y: row, Amyg bottom
  hue <- if (conn_type == "hyper") "#B30000" else "#003366"
  # Separate panels, so BOTH triangles keep the diagonal (within-network) cells.
  cdat <- long[long$yi >= long$xi, ]   # count triangle + diagonal
  pdat <- long[long$yi <= long$xi, ]   # proportion triangle + diagonal

  labs <- READABLE[nets]
  base_family <- "Helvetica"
  N <- length(nets)

  # One self-contained TRIANGLE panel (only its half is drawn → triangular shape).
  # x_pos / y_pos place the network-name axes on the requested sides.
  tri_panel <- function(dat, fillvar, legname, x_pos = "bottom", y_pos = "left") {
    ggplot(dat, aes(xi, yi, fill = .data[[fillvar]])) +
      geom_tile(color = "grey55", linewidth = 0.4) +
      scale_fill_gradient(low = "white", high = hue, name = legname,
                          guide = guide_colorbar(barheight = grid::unit(3, "cm"))) +
      scale_x_continuous(breaks = seq_len(N), labels = labs, position = x_pos,
                         limits = c(0.5, N + 0.5), expand = c(0, 0)) +
      scale_y_continuous(breaks = seq_len(N), labels = labs, position = y_pos,
                         limits = c(0.5, N + 0.5), expand = c(0, 0)) +
      coord_fixed() +
      theme_minimal(base_family = base_family) +
      theme(
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = if (x_pos == "top") 0 else 1, size = 9),
        axis.text.y = element_text(size = 9),
        legend.title = element_text(size = 9), legend.text = element_text(size = 8)
      )
  }

  # Left panel: network names on TOP; right panel: network names on the RIGHT.
  p_count <- tri_panel(cdat, "count", "No. of edges", x_pos = "top",    y_pos = "left")
  p_prop  <- tri_panel(pdat, "prop", "Normalised\nproportion", x_pos = "bottom", y_pos = "right")
  combined <- plot_grid(p_count, p_prop, nrow = 1, align = "h")

  out <- file.path(IO_DIR, sprintf("autism_vs_td_%s_double_matrix_final.pdf", conn_type))
  ggsave(out, combined, width = 15, height = 7)  # base pdf device (Helvetica core font)
  write.csv(cnt, file.path(IO_DIR, sprintf("autism_vs_td_%s_raw_count_matrix.csv", conn_type)))
  write.csv(prop, file.path(IO_DIR, sprintf("autism_vs_td_%s_proportion_matrix.csv", conn_type)))
  cat("Saved:", basename(out), sprintf("(%d edges)\n", nrow(df)))
}

cat("Double network matrix — input/output dir:", IO_DIR, "\n")
for (t in c("hyper", "hypo")) plot_double(t)
