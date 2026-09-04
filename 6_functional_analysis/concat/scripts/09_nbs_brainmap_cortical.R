# =============================================================================
# Whole-group Autism-vs-NT NBS brain map (figure_6b style) — section 6-2
# =============================================================================
# Analogue of 7_cluster_functional_analysis/scripts/07_figure6b_nbs_cortical.R,
# but the contrast is whole-group Autism-vs-NT (not per cluster) and the COLUMNS
# are the XCP-D v0.11 GSR VARIANTS (nogsr | nogsr_concat | gsr | gsr_concat).
#
# Per direction block (Hypoconnectivity blue | Hyperconnectivity red):
#   - cortical row : Schaefer-7 surface, fill = node-wise NBS-component edge count
#   - subcortical  : yabplot inset (10_*.png) embedded beneath each variant
#
# Produces:
#   - combined:   figures/figure_nbs_autism_vs_nt_all_variants.{pdf,png}
#   - per-variant: figures/nbs_wholegroup_<v>/figure_nbs_autism_vs_nt_<v>.{pdf,png}
#
# Inputs: tables/nbs_wholegroup_<v>/{nbs_edges.csv,nbs_summary.csv} (08's output),
#         <v>/nbs_double/autism_vs_nt_nbs_subcortical_<dir>.png (10's output,
#         same copy 11_nbs_combined_figure.R uses -- single source of truth,
#         no separate figures/nbs_wholegroup_<v>/ copy to go stale).
# =============================================================================
suppressMessages({
  library(plyr); library(dplyr); library(ggplot2)
  library(ggseg); library(ggsegSchaefer); library(patchwork)
  library(magick); library(grid)
})

args        <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
script_dir  <- if (length(script_path) == 0) getwd() else dirname(normalizePath(script_path))
section_dir <- dirname(script_dir)
out_root    <- file.path(section_dir, "outputs")
fig_dir     <- file.path(out_root, "figures")
tables_dir  <- file.path(out_root, "tables")
ATLAS_FILE  <- "/Users/mfleury/POSTDOC/LIBRAIRY/eeg_mri-pipeline/ressources/atlases/SCHAEFER/atlas-4S156Parcels/atlas-4S156Parcels_dseg.tsv"

atlas   <- read.table(ATLAS_FILE, sep = "\t", header = TRUE)
lab2net <- setNames(atlas$label_7network, atlas$label)

# variants = env VARIANTS (comma-sep) or auto-detected nbs_wholegroup_* table dirs.
env_v <- Sys.getenv("VARIANTS", "")
if (nzchar(env_v)) {
  VARIANTS <- trimws(strsplit(env_v, ",")[[1]])
} else {
  VARIANTS <- sub("^nbs_wholegroup_", "",
                  basename(Sys.glob(file.path(tables_dir, "nbs_wholegroup_*"))))
}
VARIANTS <- VARIANTS[vapply(VARIANTS, function(v)
  file.exists(file.path(tables_dir, paste0("nbs_wholegroup_", v), "nbs_edges.csv")), logical(1))]
cat("variants:", paste(VARIANTS, collapse = ", "), "\n")

DIRS <- list(list(key = "hypo",  title = "Hypoconnectivity",  low = "#deebf7", high = "#08519c"),
             list(key = "hyper", title = "Hyperconnectivity", low = "#fee0d2", high = "#a50f15"))

# figures: per-variant figure output. tables: nbs_edges.csv input.
# nbs_double: subcortical inset input -- same dir 11_nbs_combined_figure.R reads.
# The primary regime (nogsr_concat) lives at outputs/<v>/; the other sweep
# variants (gsr, gsr_concat, nogsr) live under outputs/sensitivity/<v>/.
vdir <- function(v) file.path(fig_dir, paste0("nbs_wholegroup_", v))
tdir <- function(v) file.path(tables_dir, paste0("nbs_wholegroup_", v))
ndir <- function(v) {
  direct <- file.path(out_root, v, "nbs_double")
  if (dir.exists(direct)) return(direct)
  file.path(out_root, "sensitivity", v, "nbs_double")
}

# node-wise edge count for one variant x direction, mapped to Schaefer parcels
count_cortical <- function(variant, direction) {
  fp <- file.path(tdir(variant), "nbs_edges.csv")
  if (!file.exists(fp)) return(data.frame(region = character(), n_edges = integer()))
  d <- read.csv(fp, stringsAsFactors = FALSE)
  d <- d[d$direction == direction, , drop = FALSE]
  if (nrow(d) == 0) return(data.frame(region = character(), n_edges = integer()))
  ends <- sub("^con_", "", d$region)
  nodes <- c(sub("/.*$", "", ends), sub("^.*/", "", ends))
  parcels <- lab2net[nodes]
  parcels <- parcels[!is.na(parcels) & parcels != "n/a" & parcels != ""]
  if (length(parcels) == 0) return(data.frame(region = character(), n_edges = integer()))
  as.data.frame(table(region = parcels), stringsAsFactors = FALSE) |>
    dplyr::rename(n_edges = Freq)
}

# shared scale across everything (fair variant comparison)
GMAX <- 0
grid <- list()
for (dl in DIRS) for (v in VARIANTS) {
  cc <- count_cortical(v, dl$key); grid[[paste(dl$key, v)]] <- cc
  if (nrow(cc)) GMAX <- max(GMAX, cc$n_edges)
}
GMAX <- max(GMAX, 1)
cat(sprintf("shared No. of edges limit = %d\n", GMAX))

make_panel <- function(cc, title_str, dl) {
  if (is.null(cc) || nrow(cc) == 0)
    cc <- data.frame(region = "7Networks_LH_Vis_1", n_edges = NA_real_)
  p <- ggseg(cc, mapping = aes(fill = n_edges), atlas = schaefer7_100,
             position = "stacked", colour = "grey85", size = 0.1)
  p + scale_fill_gradient(low = dl$low, high = dl$high, na.value = "grey92",
                          limits = c(1, GMAX), name = "No. of edges") +
    theme_void(base_family = "Helvetica") +
    theme(plot.title = element_text(size = 12, hjust = 0.5, margin = margin(b = 1, t = 4)),
          legend.key.height = unit(30, "pt"), legend.key.width = unit(9, "pt"),
          legend.title = element_text(size = 9), legend.text = element_text(size = 8)) +
    labs(title = title_str)
}

sub_inset <- function(variant, direction) {
  fp <- file.path(ndir(variant), sprintf("autism_vs_nt_nbs_subcortical_%s.png", direction))
  if (!file.exists(fp)) return(patchwork::plot_spacer())
  img <- magick::image_trim(magick::image_read(fp))
  wrap_elements(full = grid::rasterGrob(as.raster(img), interpolate = TRUE))
}

row_header <- function(txt) ggplot() +
  annotate("text", 0, 0.5, label = txt, hjust = 0, vjust = 0.5, size = 5.4,
           fontface = "bold", family = "Helvetica") +
  xlim(0, 1) + ylim(0, 1) + theme_void()

# --- build one composite over an arbitrary set of column variants -------------
build_fig <- function(cols, title, subtitle) {
  blocks <- lapply(DIRS, function(dl) {
    cort <- lapply(cols, function(v) make_panel(grid[[paste(dl$key, v)]], v, dl))
    crow <- wrap_plots(cort, nrow = 1) +
      plot_layout(guides = "collect") & theme(legend.position = "right")
    srow <- wrap_plots(lapply(cols, sub_inset, direction = dl$key), nrow = 1)
    wrap_plots(list(row_header(dl$title), crow, srow), ncol = 1, heights = c(0.07, 1, 0.42))
  })
  wrap_plots(blocks, ncol = 1) +
    plot_annotation(title = title, subtitle = subtitle,
      theme = theme(plot.title    = element_text(size = 17, hjust = 0.5, family = "Helvetica"),
                    plot.subtitle = element_text(size = 10, hjust = 0.5, family = "Helvetica")))
}

SUB <- paste0("Node-wise edge count — largest NBS component per direction ",
              "(edge threshold p<0.01). XCP-D v0.11, merged thalamus.")

# 1) combined figure: variants as columns
if (length(VARIANTS) >= 1) {
  g <- build_fig(VARIANTS, "Functional connectivity — Autism vs NT (NBS, GSR variants)", SUB)
  w <- max(6, 3.2 * length(VARIANTS))
  out <- file.path(fig_dir, "figure_nbs_autism_vs_nt_all_variants.pdf")
  ggsave(out, g, width = w, height = 13, limitsize = FALSE)
  ggsave(sub("\\.pdf$", ".png", out), g, width = w, height = 13, dpi = 200, limitsize = FALSE)
  cat("Saved:", out, "\n")
}

# 2) per-variant figures (single Autism-vs-NT column)
for (v in VARIANTS) {
  g <- build_fig(v, sprintf("Functional connectivity — Autism vs NT (NBS, %s)", v), SUB)
  out <- file.path(vdir(v), sprintf("figure_nbs_autism_vs_nt_%s.pdf", v))
  ggsave(out, g, width = 5.5, height = 13, limitsize = FALSE)
  ggsave(sub("\\.pdf$", ".png", out), g, width = 5.5, height = 13, dpi = 200, limitsize = FALSE)
  cat("Saved:", out, "\n")
}
