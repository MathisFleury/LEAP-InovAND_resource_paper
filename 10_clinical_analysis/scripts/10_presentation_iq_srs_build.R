#!/usr/bin/env Rscript
# =============================================================================
# 08 - IQ x SRS cluster figure, built step-by-step for an oral presentation
# =============================================================================
# Progressive reveal of the cluster/LOEUF IQ x SRS panel (same style as 04d).
#
# IMPORTANT: to keep the PLOTTING PANEL IN THE EXACT SAME POSITION across every
# slide, all slides carry the FULL slide-6 legends (Cluster + LOEUF + Carrier)
# — so nothing shifts/shrinks between frames. Objects are revealed by toggling
# alpha (invisible = alpha 0) rather than by adding/removing legend-bearing
# layers. Layers/scales are identical on every slide; only visibility changes.
#
#   1  axes + full (empty) legends
#   2  + clinical thresholds (IQ 70/130, SRS-2 60/75)
#   3  + NT (grey) vs Cluster C1
#   4  + Cluster C2
#   5  + Cluster C3
#   6  + carrier rings (SynGO/ChromEpiTF) + gene names
#
# Reads outputs/tables/iq_srs_plot_frame.csv (built by 04_plot_iq_srs_scatter.py).
# Writes outputs/figures/presentation/slide_1..6.pdf + presentation_iq_srs.pdf.
# =============================================================================
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2)
  library(ggrepel); library(ggnewscale)
})

here <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
if (length(here) == 0) here <- "."
FIG_DIR <- normalizePath(file.path(here, "..", "outputs", "figures", "presentation"), mustWork = FALSE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
FRAME <- normalizePath(file.path(here, "..", "outputs", "tables", "iq_srs_plot_frame.csv"))

W_CM <- 22; H_CM <- 15
AXIS_TXT <- 12; AXIS_TITLE <- 13; LEG_TXT <- 10; LEG_TITLE <- 11; GENE_TXT <- 7; TITLE_TXT <- 15
CIRCLE_MULT <- 0.30; DASH_LW <- 0.5

CLUS_LEVELS <- c("NT", "C1", "C2", "C3")
PAL_CLUSTER <- c(NT="#C1C2BC", C1="#7A8B47", C2="#ff9fa0", C3="#e7ba52")
LOEUF_SIZE  <- c(`1`=60, `2`=120, `3`=220, `4`=450, `5`=850)
LOEUF_LAB   <- c(`5`="LOEUF < 0.10", `4`="0.10 <= LOEUF < 0.20", `3`="0.20 <= LOEUF < 0.36",
                 `2`="0.36 <= LOEUF < 0.45", `1`="LOEUF >= 0.45")
LOEUF_ANNOTATE_BELOW <- 0.36
CARRIERS <- list(
  list(col = "gs_syngo", label = "SynGO",              colour = "#4a90e2"),
  list(col = "gs_chrom", label = "ChromEpiTF",         colour = "#8e44ad"),
  list(col = "gs_both",  label = "SynGO + ChromEpiTF", colour = "#008000"))
RING_COLS <- setNames(sapply(CARRIERS, `[[`, "colour"), sapply(CARRIERS, `[[`, "label"))

area_to_mm <- function(area) sqrt(area) * CIRCLE_MULT
size_vals <- area_to_mm(LOEUF_SIZE); names(size_vals) <- names(LOEUF_SIZE)

df <- read_csv(FRAME, show_col_types = FALSE) %>% mutate(ID = as.character(ID))
df <- df[df$Cluster %in% CLUS_LEVELS, ]
df$band <- factor(df$LOEUF, levels = names(LOEUF_SIZE))
df$hue  <- factor(df$Cluster, levels = CLUS_LEVELS)
known <- df[!df$was_nan_LOEUF, ]; miss <- df[df$was_nan_LOEUF, ]
car <- do.call(rbind, lapply(CARRIERS, function(s) {
  x <- df[df[[s$col]] == 1 & !df$was_nan_LOEUF, ]; if (nrow(x)) transform(x, ring = s$label) else NULL }))
car$ring <- factor(car$ring, levels = names(RING_COLS))
lab_all <- df[Reduce(`|`, lapply(CARRIERS, function(s) df[[s$col]] == 1)) &
              !is.na(df$gene) & df$LOEUF_raw < LOEUF_ANNOTATE_BELOW, ]

XLIM <- range(df$SRS_tscore, na.rm = TRUE) + c(-3, 3)
YLIM <- range(df$IQ, na.rm = TRUE) + c(-5, 5)

# reveal spec per slide: which clusters are visible, whether rings/labels/lines show
UNIFORM_PT <- 2.0  # fixed marker size on slides 3-5 (before the LOEUF reveal)

make_slide <- function(title, clusters = character(0),
                       rings = FALSE, labels = FALSE, thresholds = FALSE, size_reveal = FALSE) {
  a_known <- ifelse(known$Cluster %in% clusters, 0.70, 0)
  a_miss  <- ifelse(miss$Cluster  %in% clusters, 0.70, 0)
  a_ring  <- if (rings) 0.9 else 0

  # invisible layer that ALWAYS renders the LOEUF size legend, so the panel
  # stays put even before the circle-size reveal on slide 6.
  g <- ggplot() +
    geom_point(data = known, aes(SRS_tscore, IQ, size = band), colour = "grey40", alpha = 0) +
    scale_size_manual(values = size_vals, labels = LOEUF_LAB, name = "LOEUF",
                      guide = guide_legend(override.aes = list(colour = "grey40", alpha = 1)))
  # visible cluster points: uniform size until slide 6, then LOEUF-sized
  if (size_reveal) {
    # LOEUF reveal: known = sized circles, missing-LOEUF = small diamonds
    g <- g +
      geom_point(data = known, aes(SRS_tscore, IQ, colour = hue, size = band, alpha = a_known)) +
      geom_point(data = miss, aes(SRS_tscore, IQ, colour = hue, alpha = a_miss), shape = 18, size = 1.6)
  } else {
    # before the reveal: ALL individuals are identical uniform circles
    a_all <- ifelse(df$Cluster %in% clusters, 0.70, 0)
    g <- g + geom_point(data = df, aes(SRS_tscore, IQ, colour = hue, alpha = a_all), size = UNIFORM_PT)
  }
  g <- g +
    scale_colour_manual(values = PAL_CLUSTER, name = "Cluster", drop = FALSE,
                        guide = guide_legend(override.aes = list(alpha = 1, size = 2.5))) +
    scale_alpha_identity(guide = "none") +
    new_scale_colour() +
    # carrier rings — always present (keeps the Carrier legend), revealed via alpha
    geom_point(data = car, aes(SRS_tscore, IQ, size = band, colour = ring),
               shape = 21, fill = NA, stroke = 0.6, alpha = a_ring, show.legend = TRUE) +
    scale_colour_manual(values = RING_COLS, name = "Carrier",
                        guide = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    coord_cartesian(xlim = XLIM, ylim = YLIM) +
    labs(x = "SRS-2 total T-score  (higher = more social-communication difficulties)",
         y = "Measured IQ  (higher = better cognition)", title = title) +
    scale_x_continuous(breaks = seq(40, 90, 10)) +
    scale_y_continuous(breaks = seq(40, 160, 20)) +
    theme_classic(base_size = AXIS_TITLE) +
    theme(plot.title = element_text(size = TITLE_TXT, face = "bold"),
          legend.text = element_text(size = LEG_TXT), legend.title = element_text(size = LEG_TITLE),
          legend.box = "vertical", legend.key.size = unit(0.9, "lines"),
          axis.text = element_text(size = AXIS_TXT, colour = "black"),
          axis.title = element_text(size = AXIS_TITLE, colour = "black"))

  if (thresholds)
    g <- g +
      geom_vline(xintercept = c(60, 75), linetype = "dashed", linewidth = DASH_LW) +
      geom_hline(yintercept = c(70, 130), linetype = "dashed", linewidth = DASH_LW) +
      annotate("text", x = c(60, 75), y = YLIM[2], label = c("SRS 60", "SRS 75"),
               vjust = 1.2, size = GENE_TXT / .pt, colour = "grey30") +
      annotate("text", x = XLIM[1], y = c(70, 130), label = c("IQ 70 (ID)", "IQ 130"),
               hjust = 0, vjust = -0.4, size = GENE_TXT / .pt, colour = "grey30")
  if (labels)
    g <- g + geom_text_repel(data = lab_all, aes(SRS_tscore, IQ, label = gene),
                             fontface = "italic", size = GENE_TXT / .pt, colour = "black",
                             segment.colour = "grey40", segment.size = 0.2, min.segment.length = 0,
                             box.padding = 0.3, max.overlaps = Inf, force = 3, max.iter = 20000, seed = 42)
  g
}

slides <- list(
  make_slide("1  The clinical space"),
  make_slide("2  Clinical thresholds", thresholds = TRUE),
  make_slide("3  Neurotypical (NT) vs Cluster C1", c("NT", "C1"), thresholds = TRUE),
  make_slide("4  + Cluster C2", c("NT", "C1", "C2"), thresholds = TRUE),
  make_slide("5  + Cluster C3", c("NT", "C1", "C2", "C3"), thresholds = TRUE),
  make_slide("6  Genetic constraint: LOEUF circle size",
             c("NT", "C1", "C2", "C3"), thresholds = TRUE, size_reveal = TRUE),
  make_slide("7  + carrier rings + gene names",
             c("NT", "C1", "C2", "C3"), rings = TRUE, labels = TRUE, thresholds = TRUE, size_reveal = TRUE)
)

for (i in seq_along(slides)) {
  ggsave(file.path(FIG_DIR, sprintf("slide_%d.pdf", i)), slides[[i]],
         width = W_CM, height = H_CM, units = "cm")
  cat("  Saved slide", i, "\n")
}
pdf(file.path(FIG_DIR, "presentation_iq_srs.pdf"), width = W_CM / 2.54, height = H_CM / 2.54)
for (g in slides) print(g)
invisible(dev.off())
cat("Saved: presentation_iq_srs.pdf (", length(slides), "pages ) ->", FIG_DIR, "\n")
