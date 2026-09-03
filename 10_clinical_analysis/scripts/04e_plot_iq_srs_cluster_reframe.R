#!/usr/bin/env Rscript
# =============================================================================
# 04e - Cluster / LOEUF IQ×SRS panel (reframed), authored at 16 cm print width
# =============================================================================
# Companion to 04c (gene-set panel @ 18 cm): the cluster-coloured, LOEUF-sized
# IQ×SRS scatter (SynGO / ChromEpiTF / both rings + gene names), authored at
# EXACTLY 16 cm wide. Reads the frame dumped by 04_plot_iq_srs_scatter.py, so it
# reflects the current gnomAD v4 genetics.
# Output: figures/figure_IQ_SRS_cluster_loeuf_reframe.{pdf,svg}
#
# NOTE: the original generator for this filename was not in the repo; this
# reconstructs it from the 04c template at 16 cm. Adjust content if the prior
# version differed.
# =============================================================================
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2)
  library(ggrepel); library(ggnewscale)
})

here <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
if (length(here) == 0) here <- "."
FIG_DIR <- normalizePath(file.path(here, "..", "outputs", "figures"), mustWork = FALSE)
FRAME   <- normalizePath(file.path(here, "..", "outputs", "tables", "iq_srs_plot_frame.csv"))

# ---- print-size knobs (16 cm wide) -----------------------------------------
W_CM <- 16; H_CM <- 11.5
AXIS_TXT <- 7; AXIS_TITLE <- 8; LEG_TXT <- 6; LEG_TITLE <- 6.5; GENE_TXT <- 5
CIRCLE_MULT <- 0.21; DIAMOND_SZ <- 1.3
DASH_LW <- 0.35; AXIS_LW <- 0.4; SEG_LW <- 0.2

PAL_CLUSTER <- c(NT="#C1C2BC", C1="#7A8B47", C2="#ff9fa0", C3="#e7ba52", IDD="#D8A4CB")
LOEUF_SIZE <- c(`1`=60, `2`=120, `3`=220, `4`=450, `5`=850)
LOEUF_LAB  <- c(`5`="LOEUF < 0.10", `4`="0.10 <= LOEUF < 0.20", `3`="0.20 <= LOEUF < 0.36",
                `2`="0.36 <= LOEUF < 0.45", `1`="LOEUF >= 0.45")
LOEUF_ANNOTATE_BELOW <- 0.36
CARRIERS <- list(
  list(col = "gs_syngo", label = "SynGO",              colour = "#4a90e2"),
  list(col = "gs_chrom", label = "ChromEpiTF",         colour = "#8e44ad"),
  list(col = "gs_both",  label = "SynGO + ChromEpiTF", colour = "#008000"))

area_to_mm <- function(area) sqrt(area) * CIRCLE_MULT
df <- read_csv(FRAME, show_col_types = FALSE) %>% mutate(ID = as.character(ID))

df$band <- factor(df$LOEUF, levels = names(LOEUF_SIZE))
size_vals <- area_to_mm(LOEUF_SIZE); names(size_vals) <- names(LOEUF_SIZE)
df$hue <- factor(df$Cluster, levels = names(PAL_CLUSTER))
known <- df[!df$was_nan_LOEUF, ]; miss <- df[df$was_nan_LOEUF, ]

ring_cols <- setNames(sapply(CARRIERS, `[[`, "colour"), sapply(CARRIERS, `[[`, "label"))
car <- do.call(rbind, lapply(CARRIERS, function(s) {
  sub <- df[df[[s$col]] == 1 & !df$was_nan_LOEUF, ]
  if (nrow(sub)) transform(sub, ring = s$label) else NULL
}))
if (!is.null(car)) car$ring <- factor(car$ring, levels = names(ring_cols))
cc <- Reduce(`|`, lapply(CARRIERS, function(s) df[[s$col]] == 1))
lab <- df[cc & !is.na(df$gene) & df$LOEUF_raw < LOEUF_ANNOTATE_BELOW, ]

g <- ggplot() +
  geom_point(data = known, aes(SRS_tscore, IQ, colour = hue, size = band), alpha = 0.7) +
  geom_point(data = miss, aes(SRS_tscore, IQ, colour = hue), shape = 18, size = DIAMOND_SZ, alpha = 0.7) +
  scale_colour_manual(values = PAL_CLUSTER, name = "Cluster", drop = TRUE) +
  scale_size_manual(values = size_vals, labels = LOEUF_LAB, name = "LOEUF",
                    guide = guide_legend(override.aes = list(colour = "grey40"))) +
  new_scale_colour() +
  geom_point(data = car, aes(SRS_tscore, IQ, size = band, colour = ring),
             shape = 21, fill = NA, stroke = 0.35) +
  scale_colour_manual(values = ring_cols, name = "Carrier",
                      guide = guide_legend(override.aes = list(size = 2.5))) +
  geom_vline(xintercept = c(60, 75), linetype = "dashed", linewidth = DASH_LW) +
  geom_hline(yintercept = c(70, 130), linetype = "dashed", linewidth = DASH_LW) +
  geom_text_repel(data = lab, aes(SRS_tscore, IQ, label = gene), fontface = "italic",
                  size = GENE_TXT / .pt, colour = "black", segment.colour = "grey40",
                  segment.size = SEG_LW, min.segment.length = 0, box.padding = 0.3,
                  point.padding = 0.2, max.overlaps = Inf, force = 3, max.iter = 20000, seed = 42) +
  labs(x = "SRS-2 total T-score", y = "Measured IQ") +
  scale_x_continuous(breaks = seq(40, 90, 10)) +
  scale_y_continuous(breaks = seq(40, 160, 20)) +
  theme_classic(base_size = AXIS_TITLE) +
  theme(aspect.ratio = 1 / 1.35,   # panel x-axis 1.35x wider than y-axis
        legend.key.size = unit(0.5, "lines"),
        legend.text = element_text(size = LEG_TXT),
        legend.title = element_text(size = LEG_TITLE),
        plot.margin = margin(3, 3, 3, 3),
        axis.line = element_line(linewidth = AXIS_LW),
        axis.ticks = element_line(linewidth = AXIS_LW),
        axis.text = element_text(size = AXIS_TXT, colour = "black"),
        axis.title = element_text(size = AXIS_TITLE, colour = "black"))

for (ext in c("pdf", "svg")) {
  out <- file.path(FIG_DIR, paste0("figure_IQ_SRS_cluster_loeuf_reframe.", ext))
  ggsave(out, g, width = W_CM, height = H_CM, units = "cm",
         device = if (ext == "svg") svglite::svglite else pdf)
  cat("  Saved:", out, "\n")
}
