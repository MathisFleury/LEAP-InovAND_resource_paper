#!/usr/bin/env Rscript
# =============================================================================
# 08b - IQ x SRS presentation build — HCNDD genetics variant
# =============================================================================
# Same step-by-step, fixed-panel build as 10_presentation_iq_srs_build.R, but
# the genetic layer (slide 6) highlights HIGH-CONFIDENCE NDD (HCNDD, constrained
# dominant + rec-X-linked) carriers instead of SynGO / ChromEpiTF.
#
# All slides carry the full legend (Cluster + LOEUF + Carrier) so the plotting
# panel is in the IDENTICAL position on every frame; objects reveal via alpha.
#
# Reads outputs/tables/iq_srs_plot_frame.csv (must include the HCNDD best-gene
# columns — rerun 04_plot_iq_srs_scatter.py first).
# Writes outputs/figures/presentation/hcndd_slide_1..6.pdf + presentation_iq_srs_hcndd.pdf.
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
# single carrier set: HCNDD (constrained dominant + rec-X-linked)
RING_COLS <- c("HCNDD" = "#1b7837")
CARRIER_COL <- "dellof_hcndddomv6_contraint_carrier"
HCNDD_GENE  <- "dellof_hcndddomv6_constraint_genes_best_gene"
HCNDD_SCORE <- "dellof_hcndddomv6_constraint_genes_best_score"

area_to_mm <- function(area) sqrt(area) * CIRCLE_MULT
size_vals <- area_to_mm(LOEUF_SIZE); names(size_vals) <- names(LOEUF_SIZE)

df <- read_csv(FRAME, show_col_types = FALSE) %>% mutate(ID = as.character(ID))
stopifnot(all(c(CARRIER_COL, HCNDD_GENE, HCNDD_SCORE) %in% names(df)))
df <- df[df$Cluster %in% CLUS_LEVELS, ]
df$band <- factor(df$LOEUF, levels = names(LOEUF_SIZE))
df$hue  <- factor(df$Cluster, levels = CLUS_LEVELS)
known <- df[!df$was_nan_LOEUF, ]; miss <- df[df$was_nan_LOEUF, ]
car <- df[df[[CARRIER_COL]] == 1 & !df$was_nan_LOEUF, ]
if (nrow(car)) car$ring <- factor("HCNDD", levels = names(RING_COLS))
lab_all <- df[df[[CARRIER_COL]] == 1 & !is.na(df[[HCNDD_GENE]]) & df[[HCNDD_SCORE]] < LOEUF_ANNOTATE_BELOW, ]

XLIM <- range(df$SRS_tscore, na.rm = TRUE) + c(-3, 3)
YLIM <- range(df$IQ, na.rm = TRUE) + c(-5, 5)

UNIFORM_PT <- 2.0  # fixed marker size on slides 3-5 (before the LOEUF reveal)

make_slide <- function(title, clusters = character(0),
                       rings = FALSE, labels = FALSE, thresholds = FALSE, size_reveal = FALSE) {
  a_known <- ifelse(known$Cluster %in% clusters, 0.70, 0)
  a_miss  <- ifelse(miss$Cluster  %in% clusters, 0.70, 0)
  a_ring  <- if (rings) 0.9 else 0

  # invisible layer that ALWAYS renders the LOEUF size legend (panel stays fixed)
  g <- ggplot() +
    geom_point(data = known, aes(SRS_tscore, IQ, size = band), colour = "grey40", alpha = 0) +
    scale_size_manual(values = size_vals, labels = LOEUF_LAB, name = "LOEUF",
                      guide = guide_legend(override.aes = list(colour = "grey40", alpha = 1)))
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
    g <- g + geom_text_repel(data = lab_all, aes(SRS_tscore, IQ, label = .data[[HCNDD_GENE]]),
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
  make_slide("7  + HCNDD carrier rings + gene names",
             c("NT", "C1", "C2", "C3"), rings = TRUE, labels = TRUE, thresholds = TRUE, size_reveal = TRUE)
)

for (i in seq_along(slides)) {
  ggsave(file.path(FIG_DIR, sprintf("hcndd_slide_%d.pdf", i)), slides[[i]],
         width = W_CM, height = H_CM, units = "cm")
  cat("  Saved hcndd slide", i, "\n")
}
pdf(file.path(FIG_DIR, "presentation_iq_srs_hcndd.pdf"), width = W_CM / 2.54, height = H_CM / 2.54)
for (g in slides) print(g)
invisible(dev.off())
cat("Saved: presentation_iq_srs_hcndd.pdf (", length(slides), "pages ) ->", FIG_DIR, "\n")
