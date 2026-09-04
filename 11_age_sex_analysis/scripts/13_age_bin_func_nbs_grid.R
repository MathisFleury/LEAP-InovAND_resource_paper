#!/usr/bin/env Rscript
# =============================================================================
# Age-binned functional NBS results, ONE figure: rows = direction (Autism>NT
# hyper / Autism<NT hypo), columns = age bin. Both directions are always
# shown -- picking only the "winning" sign per bin would hide, e.g., 10-13's
# hypo component or 14-17's hyper component even when non-significant, which
# matters for judging the overall developmental pattern. Each panel is that
# bin+sign's top (largest) component's network x network edge-count matrix;
# shown grey + "n.s." when it doesn't survive p_FWER < 0.05. Adapts
# 03_plot_func_interaction_double_matrix.R's network-matrix builder to
# 12_age_bin_func_nbs.py's per-bin outputs (scanner[+motion][+sex]-adjusted).
#
# Reads: outputs/tables/nbs/age_bin_nbs_edges_<bin>_<sign>.csv
#        outputs/tables/nbs/age_bin_nbs_<bin>_<sign>.csv          (for p_FWER)
# Output: outputs/figures/nbs/composite/age_bin_nbs_grid.pdf
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(reshape2); library(ggplot2)
  library(grid); library(patchwork)
})

ATLAS_FILE <- "/Users/mfleury/POSTDOC/LIBRAIRY/eeg_mri-pipeline/ressources/atlases/SCHAEFER/atlas-4S156Parcels/atlas-4S156Parcels_dseg.tsv"

args <- commandArgs(trailingOnly = FALSE)
sp <- sub("--file=", "", args[grep("--file=", args)])
script_dir <- if (length(sp) == 0) getwd() else dirname(normalizePath(sp))
IO_DIR <- file.path(dirname(script_dir), "outputs", "tables", "nbs")
OUT_DIR <- file.path(dirname(script_dir), "outputs", "figures", "nbs", "composite")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

NET_ORDER <- c("Amyg. & Hippoc.", "Striatum", "Cerebellum", "Thalamus",
               "Default", "Cont", "Limbic", "SalVentAttn", "DorsAttn", "SomMot", "Vis")
READABLE <- c("Amyg. & Hippoc." = "Amyg.&Hip.", "Striatum" = "Striat.",
              "Cerebellum" = "Cereb.", "Thalamus" = "Thal.",
              "Default" = "DMN", "Cont" = "FPN", "Limbic" = "Limbic",
              "SalVentAttn" = "VAN", "DorsAttn" = "DAN", "SomMot" = "SomMot", "Vis" = "Vis")

# order = chronological, left-to-right columns
# <6 kept separate from 6-9 (2026-08-31, matching the anatomical bins): <6 has
# zero NT subjects with usable resting-state connectivity, so it's untestable,
# not just sparse -- its column header shows "(A=n, NT=0)" via group_counts.
BINS <- c("lt6" = "<6", "6_9" = "6-9", "10_13" = "10-13", "14_17" = "14-17",
          "18_24" = "18-24", "25plus" = "25+")
SIGNS <- list(pos = list(label = "Autism > NT (hyper)", hue = "#B30000"),
              neg = list(label = "Autism < NT (hypo)",  hue = "#003366"))

GROUP_COUNTS <- read_csv(file.path(IO_DIR, "age_bin_nbs_group_counts.csv"), show_col_types = FALSE)
bin_header <- function(bin_key) {
  lbl <- BINS[[bin_key]]
  gc <- GROUP_COUNTS[GROUP_COUNTS$bin == lbl, ]
  if (nrow(gc) == 0) return(lbl)
  sprintf("%s\n(A=%d, NT=%d)", lbl, gc$n_autism[1], gc$n_nt[1])
}

atlas <- read_tsv(ATLAS_FILE, show_col_types = FALSE)
nets <- intersect(NET_ORDER, unique(atlas$network_label))

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

featlab <- function(txt) wrap_elements(grid::textGrob(
  txt, rot = 90, gp = grid::gpar(fontsize = 10, fontface = "bold", fontfamily = "Helvetica")))

load_top <- function(bin_key, sgn) {
  fp <- file.path(IO_DIR, sprintf("age_bin_nbs_%s_%s.csv", bin_key, sgn))
  if (!file.exists(fp)) return(NULL)
  r <- read_csv(fp, show_col_types = FALSE)
  if (nrow(r) == 0) return(NULL)
  r[1, ]
}

make_panel <- function(bin_key, bin_label, sgn, show_title) {
  top <- load_top(bin_key, sgn)
  title_txt <- if (show_title) bin_label else NULL
  if (is.null(top)) {
    return(ggplot() + theme_void() +
             labs(title = title_txt, subtitle = "no component") +
             theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
                   plot.subtitle = element_text(size = 7.5, hjust = 0.5)))
  }
  ep <- file.path(IO_DIR, sprintf("age_bin_nbs_edges_%s_%s.csv", bin_key, sgn))
  edges <- read_csv(ep, show_col_types = FALSE)
  cnt <- build_count_matrix(edges)
  sig <- top$p_fwer < 0.05
  hue <- if (sig) SIGNS[[sgn]]$hue else "grey40"
  idx <- setNames(seq_along(nets), nets)
  labs <- READABLE[nets]; N <- length(nets)
  long <- reshape2::melt(cnt, varnames = c("row", "col"), value.name = "count")
  long$xi <- idx[as.character(long$col)]
  long$yi <- N + 1 - idx[as.character(long$row)]  # flip directly -- avoid trans="reverse"
  sig_lab <- if (sig) sprintf("p=%.3f", top$p_fwer) else "n.s."

  ggplot(long, aes(xi, yi, fill = count)) +
    geom_tile(color = "grey70", linewidth = 0.3) +
    scale_fill_gradient(low = "white", high = hue, name = "edges") +
    scale_x_continuous(breaks = seq_len(N), labels = labs, position = "top",
                       limits = c(0.5, N + 0.5), expand = c(0, 0)) +
    scale_y_continuous(breaks = seq_len(N), labels = rev(labs),
                       limits = c(0.5, N + 0.5), expand = c(0, 0)) +
    coord_fixed(clip = "off") + theme_minimal(base_family = "Helvetica") +
    theme(panel.grid = element_blank(), axis.title = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 6.5),
          axis.text.y = element_text(size = 6.5),
          plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 7.5, hjust = 0.5),
          plot.margin = margin(4, 8, 4, 8),
          legend.position = "bottom", legend.key.height = unit(5, "pt"),
          legend.title = element_text(size = 7), legend.text = element_text(size = 6)) +
    labs(title = title_txt, subtitle = sprintf("%d edges, %s", top$n_edges, sig_lab))
}

rows <- lapply(names(SIGNS), function(sgn) {
  panels <- lapply(seq_along(BINS), function(i)
    make_panel(names(BINS)[i], bin_header(names(BINS)[i]), sgn, show_title = (sgn == names(SIGNS)[1])))
  wrap_plots(c(list(featlab(SIGNS[[sgn]]$label)), panels), nrow = 1,
             widths = c(0.12, rep(1, length(panels))))
})

combined <- wrap_plots(rows, ncol = 1) +
  plot_annotation(
    title = "Age-binned functional connectivity: Autism vs NT",
    subtitle = "Each panel = top component for that bin x direction; grey = not FWER-significant",
    theme = theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(size = 10, hjust = 0.5)))

out <- file.path(OUT_DIR, "age_bin_nbs_grid.pdf")
ggsave(out, combined, width = 40, height = 16, units = "cm", limitsize = FALSE)
cat("Saved:", out, "\n")
