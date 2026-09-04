#!/usr/bin/env Rscript
# =============================================================================
# 4 - Cluster / LOEUF IQ x PGS-intelligence panel (16 cm)
# =============================================================================
# Cluster-level analogue of 10_clinical_analysis/04d_plot_iq_srs_cluster_16cm.R:
# the cluster-coloured, LOEUF-sized scatter, but x = intelligence PGS instead of
# SRS-2, plus an OLS regression line with R^2 / p annotated on the figure.
# Reads the frame written by 2_genetic_analysis/scripts/2_iq_pgs_loeuf_figures.py
# (moved here from 2_genetic_analysis since this figure is cluster-coloured,
# matching the population-vs-cluster section split).
# Output: figures/figure_IQ_PGSint_cluster_loeuf_16cm.{pdf,svg}
# =============================================================================
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(ggplot2)
  library(ggrepel); library(ggnewscale)
})

here <- dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))
if (length(here) == 0) here <- "."
FIG_DIR <- normalizePath(file.path(here, "..", "outputs", "figures"), mustWork = FALSE)
FRAME   <- normalizePath(file.path(here, "..", "..", "2_genetic_analysis", "outputs", "tables", "iq_pgs_plot_frame.csv"))
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

W_CM <- 20; H_CM <- 11.5
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

# --- OLS IQ ~ PGS-intelligence (regression line + R^2 / p annotation) --------
fit <- lm(IQ ~ PGS_int, data = df)
sm  <- summary(fit)
r   <- as.numeric(sign(coef(fit)["PGS_int"]) * sqrt(sm$r.squared))
pv  <- sm$coefficients["PGS_int", "Pr(>|t|)"]
p_txt <- formatC(pv, format = "g", digits = 3)  # real p-value, not "< 0.001"
ttl <- bquote("PGS intelligence" ~ italic("vs.") ~
              .(sprintf("IQ scores (r = %.3f, p = %s)", r, p_txt)))

g <- ggplot() +
  geom_point(data = known, aes(PGS_int, IQ, colour = hue, size = band), alpha = 0.7) +
  geom_point(data = miss, aes(PGS_int, IQ, colour = hue), shape = 18, size = DIAMOND_SZ, alpha = 0.7) +
  scale_colour_manual(values = PAL_CLUSTER, name = "Cluster", drop = TRUE) +
  scale_size_manual(values = size_vals, labels = LOEUF_LAB, name = "LOEUF",
                    guide = guide_legend(override.aes = list(colour = "grey40"))) +
  new_scale_colour() +
  geom_point(data = car, aes(PGS_int, IQ, size = band, colour = ring),
             shape = 21, fill = NA, stroke = 0.35) +
  scale_colour_manual(values = ring_cols, name = "Carrier",
                      guide = guide_legend(override.aes = list(size = 2.5))) +
  geom_smooth(data = df, aes(PGS_int, IQ), method = "lm", formula = y ~ x,
              se = TRUE, colour = "black", linewidth = 0.5, fill = "grey70", alpha = 0.3) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = DASH_LW) +
  geom_hline(yintercept = c(70, 130), linetype = "dashed", linewidth = DASH_LW) +
  geom_text_repel(data = lab, aes(PGS_int, IQ, label = gene), fontface = "italic",
                  size = GENE_TXT / .pt, colour = "black", segment.colour = "grey40",
                  segment.size = SEG_LW, min.segment.length = 0, box.padding = 0.3,
                  point.padding = 0.2, max.overlaps = Inf, force = 3, max.iter = 20000, seed = 42) +
  labs(x = "Intelligence PGS (z-score)", y = "Measured IQ", title = ttl) +
  scale_x_continuous(limits = c(-3, 3), breaks = seq(-3, 3, 1)) +
  scale_y_continuous(breaks = seq(40, 160, 20)) +
  theme_classic(base_size = AXIS_TITLE) +
  theme(plot.title = element_text(hjust = 0.5, face = "plain", size = AXIS_TITLE),
        legend.key.size = unit(0.5, "lines"),
        legend.text = element_text(size = LEG_TXT),
        legend.title = element_text(size = LEG_TITLE),
        plot.margin = margin(3, 3, 3, 3),
        axis.line = element_line(linewidth = AXIS_LW),
        axis.ticks = element_line(linewidth = AXIS_LW),
        axis.text = element_text(size = AXIS_TXT, colour = "black"),
        axis.title = element_text(size = AXIS_TITLE, colour = "black"))

for (ext in c("pdf", "svg")) {
  out <- file.path(FIG_DIR, paste0("figure_IQ_PGSint_cluster_loeuf_16cm.", ext))
  ggsave(out, g, width = W_CM, height = H_CM, units = "cm",
         device = if (ext == "svg") svglite::svglite else pdf)
  cat("  Saved:", out, "\n")
}
