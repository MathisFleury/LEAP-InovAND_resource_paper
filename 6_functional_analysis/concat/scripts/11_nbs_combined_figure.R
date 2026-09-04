#!/usr/bin/env Rscript
# =============================================================================
# Combined NBS figure (Autism vs NT), per variant — section 6-2
# =============================================================================
# Two stacked panels:
#   TOP    "Hypoconnected network in autism"  — left: double network matrix (hypo)
#                                                right: cortical brain map (top)
#                                                       + subcortical inset (below)
#   BOTTOM "Hyperconnected network in autism"  — same, hyper.
#
# All three views (matrix / cortical / subcortical) derive from the SAME NBS
# largest-component edges (04's full_data source,target), so they are consistent.
#
# Inputs per variant (outputs/<v>/nbs_double/):
#   autism_vs_td_{hypo,hyper}connectivity_full_data.csv      (04_nbs_component_edges.py)
#   autism_vs_nt_nbs_subcortical_{hypo,hyper}.png            (10_nbs_brainmap_subcortical.py)
# Output: outputs/<v>/nbs_double/figure_nbs_combined_<v>.{pdf,png}
#
# Variants via env VARIANTS (comma-sep) or auto-detected outputs/*/nbs_double.
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(reshape2)
  library(ggplot2); library(ggseg); library(ggsegSchaefer)
  library(patchwork); library(magick); library(grid)
})

ATLAS_FILE <- "/Users/mfleury/POSTDOC/LIBRAIRY/eeg_mri-pipeline/ressources/atlases/SCHAEFER/atlas-4S156Parcels/atlas-4S156Parcels_dseg.tsv"
args <- commandArgs(trailingOnly = FALSE)
sp <- sub("--file=", "", args[grep("--file=", args)])
script_dir <- if (length(sp) == 0) getwd() else dirname(normalizePath(sp))
section_dir <- dirname(script_dir)
out_root <- file.path(section_dir, "outputs")

atlas   <- read_tsv(ATLAS_FILE, show_col_types = FALSE)
lab2net <- setNames(atlas$label_7network, atlas$label)   # node -> 7Networks_* parcel

NET_ORDER <- c("Amyg. & Hippoc.", "Striatum", "Cerebellum", "Thalamus",
               "Default", "Cont", "Limbic", "SalVentAttn", "DorsAttn", "SomMot", "Vis")
READABLE <- c("Amyg. & Hippoc." = "Amyg. & Hippoc.", "Striatum" = "Striatum",
              "Cerebellum" = "Cerebellum", "Thalamus" = "Thalamus",
              "Default" = "Default Mode", "Cont" = "Frontoparietal", "Limbic" = "Limbic",
              "SalVentAttn" = "Ventral Attention", "DorsAttn" = "Dorsal Attention",
              "SomMot" = "Somatomotor", "Vis" = "Visual")
counts_per_net <- atlas %>% count(network_label, name = "n") %>% tibble::deframe()
nets <- intersect(NET_ORDER, names(counts_per_net))
total_mat <- outer(nets, nets, Vectorize(function(a, b) {
  na <- counts_per_net[[a]]; nb <- counts_per_net[[b]]
  if (a == b) na * (na - 1) / 2 else na * nb
})); dimnames(total_mat) <- list(nets, nets)

DIRS <- list(hypo  = list(title = "Hypoconnected network in autism",  hue = "#003366",
                          low = "#deebf7", high = "#08519c"),
             hyper = list(title = "Hyperconnected network in autism", hue = "#B30000",
                          low = "#fee0d2", high = "#a50f15"))

# ---- double network matrix (two triangles), from 06's logic ------------------
build_count_matrix <- function(df) {
  d <- df %>%
    left_join(atlas %>% select(label, network_label), by = c("source" = "label")) %>% rename(src = network_label) %>%
    left_join(atlas %>% select(label, network_label), by = c("target" = "label")) %>% rename(tgt = network_label) %>%
    filter(src %in% nets, tgt %in% nets)
  m <- matrix(0, length(nets), length(nets), dimnames = list(nets, nets))
  for (k in seq_len(nrow(d))) { i <- d$src[k]; j <- d$tgt[k]
    m[i, j] <- m[i, j] + 1; if (i != j) m[j, i] <- m[j, i] + 1 }
  m
}

matrix_panel <- function(df, dl) {
  cnt <- build_count_matrix(df); prop <- cnt / total_mat; prop[is.nan(prop)] <- 0
  idx <- setNames(seq_along(nets), nets)
  long <- reshape2::melt(cnt, varnames = c("row", "col"), value.name = "count")
  long$prop <- reshape2::melt(prop)$value
  long$xi <- idx[as.character(long$col)]; long$yi <- idx[as.character(long$row)]
  labs <- READABLE[nets]; N <- length(nets)
  tri <- function(dat, fillvar, legname, x_pos, y_pos) {
    ggplot(dat, aes(xi, yi, fill = .data[[fillvar]])) +
      geom_tile(color = "grey55", linewidth = 0.4) +
      scale_fill_gradient(low = "white", high = dl$hue, name = legname,
                          guide = guide_colorbar(barheight = grid::unit(2.2, "cm"))) +
      scale_x_continuous(breaks = seq_len(N), labels = labs, position = x_pos,
                         limits = c(0.5, N + 0.5), expand = c(0, 0)) +
      scale_y_continuous(breaks = seq_len(N), labels = labs, position = y_pos,
                         limits = c(0.5, N + 0.5), expand = c(0, 0)) +
      coord_fixed() + theme_minimal(base_family = "Helvetica") +
      theme(panel.grid = element_blank(), axis.title = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = if (x_pos == "top") 0 else 1, size = 8),
            axis.text.y = element_text(size = 8),
            legend.title = element_text(size = 8), legend.text = element_text(size = 7))
  }
  p_count <- tri(long[long$yi >= long$xi, ], "count", "No. of edges", "top", "left")
  p_prop  <- tri(long[long$yi <= long$xi, ], "prop", "Normalised\nproportion", "bottom", "right")
  p_count + p_prop
}

# ---- cortical Schaefer map (node-wise NBS-edge count) ------------------------
cortical_counts <- function(df) {
  nodes <- c(df$source, df$target)
  parcels <- lab2net[nodes]
  parcels <- parcels[!is.na(parcels) & parcels != "n/a" & parcels != ""]
  if (length(parcels) == 0) return(data.frame(region = character(), n_edges = integer()))
  as.data.frame(table(region = parcels), stringsAsFactors = FALSE) |> dplyr::rename(n_edges = Freq)
}
cortical_panel <- function(cc, dl, lim) {
  if (nrow(cc) == 0) cc <- data.frame(region = "7Networks_LH_Vis_1", n_edges = NA_real_)
  ggseg(cc, mapping = aes(fill = n_edges), atlas = schaefer7_100,
        position = "stacked", colour = "grey85", size = 0.1) +
    scale_fill_gradient(low = dl$low, high = dl$high, na.value = "grey92",
                        limits = c(1, lim), name = "No. of edges") +
    theme_void(base_family = "Helvetica") +
    theme(legend.key.height = unit(26, "pt"), legend.key.width = unit(8, "pt"),
          legend.title = element_text(size = 8), legend.text = element_text(size = 7))
}
sub_inset <- function(png_fp) {
  if (!file.exists(png_fp)) return(patchwork::plot_spacer())
  img <- magick::image_trim(magick::image_read(png_fp))
  wrap_elements(full = grid::rasterGrob(as.raster(img), interpolate = TRUE))
}

# ---- assemble one variant ----------------------------------------------------
build_variant <- function(v) {
  io <- file.path(out_root, v, "nbs_double")
  fp <- function(d) file.path(io, sprintf("autism_vs_td_%sconnectivity_full_data.csv", d))
  if (!all(file.exists(fp("hypo")), file.exists(fp("hyper")))) {
    cat("skip", v, "- missing full_data\n"); return(invisible())
  }
  dfs <- lapply(c(hypo = "hypo", hyper = "hyper"), function(d) read_csv(fp(d), show_col_types = FALSE))
  ccs <- lapply(dfs, cortical_counts)
  lim <- max(1, unlist(lapply(ccs, function(c) if (nrow(c)) max(c$n_edges) else 0)))

  header <- function(txt) ggplot() +
    annotate("text", 0, 0.5, label = txt, hjust = 0, size = 5.6, fontface = "bold",
             family = "Helvetica") + xlim(0, 1) + ylim(0, 1) + theme_void()
  blocks <- lapply(names(DIRS), function(d) {
    dl <- DIRS[[d]]
    mat <- matrix_panel(dfs[[d]], dl)
    cort <- cortical_panel(ccs[[d]], dl, lim)
    sub  <- sub_inset(file.path(io, sprintf("autism_vs_nt_nbs_subcortical_%s.png", d)))
    right <- wrap_plots(cort, sub, ncol = 1, heights = c(1, 0.75))
    body  <- wrap_plots(mat, right, nrow = 1, widths = c(1.7, 1))
    wrap_plots(header(dl$title), body, ncol = 1, heights = c(0.06, 1))
  })
  fig <- wrap_plots(blocks, ncol = 1) +
    plot_annotation(
      title = sprintf("Functional connectivity — Autism vs NT (NBS, %s)", v),
      subtitle = sprintf("Largest NBS component, edge threshold p<%s. XCP-D v0.11, merged thalamus.",
                         Sys.getenv("NBS_P", "0.05")),
      theme = theme(plot.title = element_text(size = 16, hjust = 0.5, family = "Helvetica"),
                    plot.subtitle = element_text(size = 9, hjust = 0.5, family = "Helvetica")))
  fig_dir <- file.path(out_root, "figures")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  out <- file.path(fig_dir, sprintf("figure_nbs_combined_%s.pdf", v))
  ggsave(out, fig, width = 17, height = 15, limitsize = FALSE)
  ggsave(sub("\\.pdf$", ".png", out), fig, width = 17, height = 15, dpi = 170, limitsize = FALSE)
  cat("Saved:", out, "\n")
}

env_v <- Sys.getenv("VARIANTS", "")
VARIANTS <- if (nzchar(env_v)) {
  trimws(strsplit(env_v, ",")[[1]])
} else {
  basename(dirname(Sys.glob(file.path(out_root, "*", "nbs_double"))))
}
VARIANTS <- unique(VARIANTS)
cat("variants:", paste(VARIANTS, collapse = ", "), "\n")
for (v in VARIANTS) build_variant(v)
