#!/usr/bin/env Rscript
# =============================================================================
# Panel a for the composite: ChromEpiTF | SynGO thickness ggseg maps, side by
# side with ONE shared colourbar (patchwork guides="collect"). Matches the
# published Figure 7 layout. Fill = OLS beta (metric ~ -log10 LOEUF);
# perm p<0.05 outlined in black.
# Output: outputs/figures_genetics_hg38_regperm/composite/panel_a_thickness.pdf
# =============================================================================
suppressMessages({
  library(plyr); library(dplyr); library(ggplot2)
  library(colorspace); library(ggseg); library(patchwork); library(scales)
})

args        <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
script_dir  <- if (length(script_path) == 0) getwd() else dirname(normalizePath(script_path))
section_dir <- dirname(script_dir)
out_base    <- Sys.getenv("LOEUF_OUT_DIR", "outputs")
input_dir   <- file.path(section_dir, out_base, "tables_genetics_hg38_regperm")
output_dir  <- file.path(section_dir, out_base, "figures_genetics_hg38_regperm", "composite")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# left -> right, matching the reference figure
PATHWAYS <- list(list(short = "chromepitf", title = "ChromEpiTF"),
                 list(short = "syngo",      title = "SynGO"))
SCALE_MAX <- 1.0

load_thickness <- function(short) {
  d <- read.csv(file.path(input_dir, sprintf("loeuf_reg_%s_thickness.csv", short)),
                stringsAsFactors = FALSE)
  d$label <- tolower(ifelse(d$hemisphere == "left",
                            paste0("lh_", gsub(" ", "", d$region)),
                            paste0("rh_", gsub(" ", "", d$region))))
  d$significant <- !is.na(d$p_permutation) & d$p_permutation < 0.05
  d
}

panel_theme <- theme_void(base_family = "Helvetica") +
  theme(plot.title = element_text(size = 15, hjust = 0.5, margin = margin(b = 1, t = 3)),
        legend.title = element_text(size = 10), legend.text = element_text(size = 8),
        legend.key.height = unit(16, "pt"), legend.key.width = unit(16, "pt"))

make_map <- function(short, title) {
  d <- load_thickness(short)
  d <- d %>% plyr::mutate(outline_color = ifelse(significant, "black", "grey70"),
                          outline_size  = ifelse(significant, 1.1, 0.15))
  d <- dplyr::select(d, label, beta, outline_color, outline_size)
  ggseg(d, mapping = aes(fill = beta, col = outline_color, size = outline_size),
        atlas = "dk", position = "stacked") +
    scale_fill_continuous_divergingx(palette = "RdBu", mid = 0, rev = TRUE,
                                     limits = c(-SCALE_MAX, SCALE_MAX), oob = scales::squish,
                                     na.value = "grey90",
                                     name = "OLS beta\nCortical Thickness\nvs -log(LOEUF)") +
    scale_colour_identity() + scale_size_identity() +
    panel_theme + labs(title = title)
}

maps <- lapply(PATHWAYS, function(pw) make_map(pw$short, pw$title))
combined <- wrap_plots(maps, nrow = 1) +
  plot_layout(guides = "collect") +
  plot_annotation(caption = "perm p < 0.05 outlined in black",
                  theme = theme(plot.caption = element_text(size = 8, hjust = 0.5,
                                                            family = "Helvetica"))) &
  # horizontal colourbar BELOW the maps -> frees horizontal space (ref layout)
  guides(fill = guide_colourbar(direction = "horizontal", title.position = "left",
                                title.vjust = 0.9, barwidth = unit(70, "pt"),
                                barheight = unit(7, "pt"))) &
  theme(legend.position = "bottom", legend.box = "horizontal")

out <- file.path(output_dir, "panel_a_thickness.pdf")
ggsave(out, combined, units = "cm", width = 15, height = 9, limitsize = FALSE)
cat("Saved:", out, "\n")
