#!/usr/bin/env Rscript
# =============================================================================
# Double network matrix of the NBS interaction component (functional).
# Adapts 6_functional_analysis/scripts/05_plot_double_network_matrix.R to the
# diagnosis x age / x sex interaction components written by 02_func_nbs_interaction.py.
#   upper-left triangle (+diag) = raw edge COUNT      (left colorbar)
#   lower-right triangle        = normalised PROPORTION (right colorbar)
# One figure per model x sign (default: age/pos = the Autism>NT contrast).
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(reshape2)
  library(ggplot2); library(ggnewscale); library(cowplot)
})

ATLAS_FILE <- "/Users/mfleury/POSTDOC/LIBRAIRY/eeg_mri-pipeline/ressources/atlases/SCHAEFER/atlas-4S156Parcels/atlas-4S156Parcels_dseg.tsv"

args <- commandArgs(trailingOnly = FALSE)
sp <- sub("--file=", "", args[grep("--file=", args)])
script_dir <- if (length(sp) == 0) getwd() else dirname(normalizePath(sp))
IO_DIR <- file.path(dirname(script_dir), "outputs", "nbs")

NET_ORDER <- c("Amyg. & Hippoc.", "Striatum", "Cerebellum", "Thalamus",
               "Default", "Cont", "Limbic", "SalVentAttn", "DorsAttn", "SomMot", "Vis")
READABLE <- c("Amyg. & Hippoc." = "Amyg. & Hippoc.", "Striatum" = "Striatum",
              "Cerebellum" = "Cerebellum", "Thalamus" = "Thalamus",
              "Default" = "Default Mode", "Cont" = "Frontoparietal",
              "Limbic" = "Limbic", "SalVentAttn" = "Ventral Attention",
              "DorsAttn" = "Dorsal Attention", "SomMot" = "Somatomotor", "Vis" = "Visual")
SIGN_LAB <- c(pos = "Autism > NT", neg = "Autism < NT")

atlas <- read_tsv(ATLAS_FILE, show_col_types = FALSE)
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

plot_double <- function(model, sign) {
  fp <- file.path(IO_DIR, sprintf("interaction_edges_%s_%s.csv", model, sign))
  if (!file.exists(fp)) { cat("Missing:", fp, "\n"); return(invisible()) }
  df <- read_csv(fp, show_col_types = FALSE)
  if (nrow(df) == 0) { cat("Empty (no component):", basename(fp), "\n"); return(invisible()) }

  cnt <- build_count_matrix(df)
  prop <- cnt / total_mat; prop[is.nan(prop)] <- 0
  idx <- setNames(seq_along(nets), nets)
  long <- reshape2::melt(cnt, varnames = c("row", "col"), value.name = "count")
  long$prop <- reshape2::melt(prop)$value
  long$xi <- idx[as.character(long$col)]
  long$yi <- idx[as.character(long$row)]
  hue <- if (sign == "pos") "#B30000" else "#003366"
  cdat <- long[long$yi >= long$xi, ]
  pdat <- long[long$yi <= long$xi, ]
  labs <- READABLE[nets]; N <- length(nets)

  tri_panel <- function(dat, fillvar, legname, x_pos, y_pos) {
    ggplot(dat, aes(xi, yi, fill = .data[[fillvar]])) +
      geom_tile(color = "grey55", linewidth = 0.4) +
      scale_fill_gradient(low = "white", high = hue, name = legname,
                          guide = guide_colorbar(barheight = grid::unit(3, "cm"))) +
      scale_x_continuous(breaks = seq_len(N), labels = labs, position = x_pos,
                         limits = c(0.5, N + 0.5), expand = c(0, 0)) +
      scale_y_continuous(breaks = seq_len(N), labels = labs, position = y_pos,
                         limits = c(0.5, N + 0.5), expand = c(0, 0)) +
      coord_fixed() + theme_minimal(base_family = "Helvetica") +
      theme(panel.grid = element_blank(), axis.title = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = if (x_pos == "top") 0 else 1, size = 9),
            axis.text.y = element_text(size = 9),
            legend.title = element_text(size = 9), legend.text = element_text(size = 8))
  }

  p_count <- tri_panel(cdat, "count", "No. of edges", "top", "left")
  p_prop  <- tri_panel(pdat, "prop", "Normalised\nproportion", "bottom", "right")
  title <- ggdraw() + draw_label(
    sprintf("Diagnosis x %s interaction — %s component (%d edges)",
            model, SIGN_LAB[[sign]], nrow(df)), fontface = "bold", size = 13)
  combined <- plot_grid(title, plot_grid(p_count, p_prop, nrow = 1, align = "h"),
                        ncol = 1, rel_heights = c(0.08, 1))

  out <- file.path(IO_DIR, sprintf("interaction_double_matrix_%s_%s.pdf", model, sign))
  ggsave(out, combined, width = 15, height = 7.5)
  write.csv(cnt, file.path(IO_DIR, sprintf("interaction_count_matrix_%s_%s.csv", model, sign)))
  cat("Saved:", basename(out), sprintf("(%d edges)\n", nrow(df)))
}

cat("Interaction double matrix — dir:", IO_DIR, "\n")
for (m in c("age", "sex")) for (s in c("pos", "neg")) plot_double(m, s)
