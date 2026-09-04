#!/usr/bin/env Rscript
#
# Cluster vs NT — Double Network Matrix Visualization
# ====================================================
#
# Creates per-cluster (C1/C2/C3) network connectivity matrices where:
#   - Upper triangle: Proportion of significant connections (normalised by total possible)
#   - Lower triangle: Raw count of significant edges
#
# Reads: cluster_{C}_vs_NT_significant_results.csv (feature, t_stat)
# Splits by t_stat sign → Hyper (>0) and Hypo (<0) matrices.
#
# Outputs (per cluster × connectivity type):
#   cluster_{C}_vs_nt_{type}_raw_count_lower.pdf
#   cluster_{C}_vs_nt_{type}_proportion_upper.pdf
#   cluster_{C}_vs_nt_{type}_raw_count_matrix.csv
#   cluster_{C}_vs_nt_{type}_normalized_proportion_matrix.csv

suppressPackageStartupMessages({
  library(plyr)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(readr)
  library(reshape2)
})

# =============================================================================
# PATHS (commandArgs-based detection for Rscript compatibility)
# =============================================================================
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
if (length(script_path) == 0) {
  script_dir <- getwd()
} else {
  script_dir <- dirname(normalizePath(script_path))
}
section_dir <- dirname(script_dir)

INPUT_DIR  <- file.path(section_dir, "outputs", "tables")    # reads 01's per-cluster CSVs
OUTPUT_DIR <- file.path(section_dir, "outputs", "figures")    # the 2 matrix plots
TABLES_DIR <- file.path(section_dir, "outputs", "tables")     # the 2 matrix CSVs
ATLAS_FILE <- "/Users/mfleury/POSTDOC/LIBRAIRY/eeg_mri-pipeline/ressources/atlases/SCHAEFER/atlas-4S156Parcels/atlas-4S156Parcels_dseg.tsv"

# =============================================================================
# NETWORK ORDER & LABELS
# =============================================================================
NETWORK_ORDER <- c(
  "Amyg. & Hippoc.", "Striatum", "Subcortical", "Cerebellum", "Thalamus",
  "Default", "Cont", "Limbic", "SalVentAttn", "DorsAttn", "SomMot", "Vis"
)

READABLE_LABELS <- c(
  "Amyg. & Hippoc."  = "Amyg. & Hippoc.",
  "Striatum"         = "Striatum",
  "Subcortical"      = "Subcortical",
  "Cerebellum"       = "Cerebellum",
  "Thalamus"         = "Thalamus",
  "Default"          = "Default Mode",
  "Cont"             = "Frontoparietal",
  "Limbic"           = "Limbic",
  "SalVentAttn"      = "Ventral Attention",
  "DorsAttn"         = "Dorsal Attention",
  "SomMot"           = "Somatomotor",
  "Vis"              = "Visual"
)

# =============================================================================
# ATLAS LOADING
# =============================================================================
load_atlas_info <- function() {
  if (!file.exists(ATLAS_FILE)) {
    stop("Atlas file not found: ", ATLAS_FILE)
  }
  read_tsv(ATLAS_FILE, col_names = TRUE, show_col_types = FALSE)
}

calculate_total_possible_connections <- function(atlas_data) {
  network_counts <- atlas_data %>%
    filter(network_label %in% NETWORK_ORDER) %>%
    dplyr::count(network_label, name = "region_count") %>%
    deframe()

  available <- intersect(NETWORK_ORDER, names(network_counts))

  mat <- matrix(0,
    nrow = length(available), ncol = length(available),
    dimnames = list(available, available)
  )

  for (i in seq_along(available)) {
    for (j in seq_along(available)) {
      n1 <- network_counts[[available[i]]]
      n2 <- network_counts[[available[j]]]
      if (i == j) {
        mat[i, j] <- n1 * (n1 - 1) / 2
      } else {
        mat[i, j] <- n1 * n2
      }
    }
  }
  mat
}

# =============================================================================
# DATA LOADING
# =============================================================================
load_cluster_data <- function(cluster, connectivity_type) {
  file_path <- file.path(INPUT_DIR, paste0("cluster_", cluster, "_vs_NT_significant_results.csv"))
  if (!file.exists(file_path)) {
    cat("  WARNING: File not found:", file_path, "\n")
    return(NULL)
  }
  data <- read_csv(file_path, show_col_types = FALSE)
  if (!all(c("t_stat", "feature") %in% names(data))) {
    cat("  WARNING: Unexpected columns in:", basename(file_path), "\n")
    return(NULL)
  }
  # Use Cohen's d for direction when available (Reviewer 3); fall back to
  # t_stat for backward compatibility with legacy CSVs.
  dir_col <- if ("cohens_d" %in% names(data)) "cohens_d" else "t_stat"
  if (connectivity_type == "hyper") {
    data <- data %>% filter(.data[[dir_col]] > 0)
  } else if (connectivity_type == "hypo") {
    data <- data %>% filter(.data[[dir_col]] < 0)
  }
  if (nrow(data) == 0) return(NULL)

  data %>% mutate(
    source = gsub("con_", "", gsub("/.*", "", feature)),
    target = gsub(".*/", "", feature)
  )
}

# =============================================================================
# MATRIX CREATION
# =============================================================================
create_network_matrix <- function(connectivity_data, atlas_data) {
  if (is.null(connectivity_data) || nrow(connectivity_data) == 0) return(NULL)

  atlas_lookup <- atlas_data %>%
    filter(network_label %in% NETWORK_ORDER) %>%
    select(label, network_label)

  data_with_networks <- connectivity_data %>%
    left_join(atlas_lookup, by = c("source" = "label")) %>%
    dplyr::rename(net_source = network_label) %>%
    left_join(atlas_lookup, by = c("target" = "label")) %>%
    dplyr::rename(net_target = network_label) %>%
    filter(!is.na(net_source), !is.na(net_target))

  if (nrow(data_with_networks) == 0) return(NULL)

  available <- intersect(NETWORK_ORDER, unique(c(data_with_networks$net_source,
                                                  data_with_networks$net_target)))
  # Use full NETWORK_ORDER intersected with what's in atlas
  all_nets <- intersect(NETWORK_ORDER,
                        atlas_data %>% filter(network_label %in% NETWORK_ORDER) %>%
                          pull(network_label) %>% unique())

  mat <- matrix(0,
    nrow = length(all_nets), ncol = length(all_nets),
    dimnames = list(all_nets, all_nets)
  )

  for (i in seq_len(nrow(data_with_networks))) {
    s <- data_with_networks$net_source[i]
    t <- data_with_networks$net_target[i]
    if (s %in% rownames(mat) && t %in% colnames(mat)) {
      mat[s, t] <- mat[s, t] + 1
      if (s != t) mat[t, s] <- mat[t, s] + 1
    }
  }
  mat
}

# =============================================================================
# PLOTTING
# =============================================================================
create_matrix_plots <- function(count_matrix, proportion_matrix, title_suffix, output_prefix) {
  if (is.null(count_matrix) || is.null(proportion_matrix)) return(NULL)

  # Apply readable labels
  rl <- READABLE_LABELS[rownames(count_matrix)]
  rownames(count_matrix)       <- rl
  colnames(count_matrix)       <- rl
  rownames(proportion_matrix)  <- rl
  colnames(proportion_matrix)  <- rl

  lower_mask       <- lower.tri(count_matrix, diag = TRUE)
  upper_diag_mask  <- upper.tri(proportion_matrix, diag = TRUE)

  count_disp <- matrix(NA, nrow = nrow(count_matrix), ncol = ncol(count_matrix),
                       dimnames = list(rownames(count_matrix), colnames(count_matrix)))
  count_disp[lower_mask] <- count_matrix[lower_mask]

  prop_disp <- matrix(NA, nrow = nrow(proportion_matrix), ncol = ncol(proportion_matrix),
                      dimnames = list(rownames(proportion_matrix), colnames(proportion_matrix)))
  prop_disp[upper_diag_mask] <- proportion_matrix[upper_diag_mask]

  colors <- c("white", "#666666")

  count_df <- reshape2::melt(count_disp)
  prop_df  <- reshape2::melt(prop_disp)
  count_df <- count_df[!is.na(count_df$value), ]
  prop_df  <- prop_df[!is.na(prop_df$value), ]

  count_plot <- ggplot(count_df, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile(color = "grey50", linewidth = 0.5) +
    scale_fill_gradientn(colours = colors, name = "No. of edges",
                         guide = guide_colorbar(title.position = "right",
                                                title.vjust = 1, direction = "vertical")) +
    theme_minimal() +
    theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y  = element_text(size = 10),
          axis.title   = element_blank(),
          panel.grid   = element_blank(),
          plot.title   = element_text(size = 13, hjust = 0.5),
          legend.position = "right") +
    labs(title = paste0("Raw Count — ", title_suffix)) +
    scale_x_discrete(position = "bottom") +
    scale_y_discrete(position = "left")

  prop_plot <- ggplot(prop_df, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile(color = "grey50", linewidth = 0.5) +
    scale_fill_gradientn(colours = colors, name = "Normalised\nProportion",
                         guide = guide_colorbar(title.position = "right",
                                                title.vjust = 1, direction = "vertical")) +
    theme_minimal() +
    theme(axis.text.x  = element_text(angle = 45, hjust = 0, size = 10),
          axis.text.y  = element_text(size = 10),
          axis.title   = element_blank(),
          panel.grid   = element_blank(),
          plot.title   = element_text(size = 13, hjust = 0.5),
          legend.position = "right") +
    labs(title = paste0("Normalised Proportion — ", title_suffix)) +
    scale_x_discrete(position = "top") +
    scale_y_discrete(position = "right")

  count_file <- file.path(OUTPUT_DIR, paste0(output_prefix, "_raw_count_lower.pdf"))
  prop_file  <- file.path(OUTPUT_DIR, paste0(output_prefix, "_proportion_upper.pdf"))
  count_csv  <- file.path(TABLES_DIR, paste0(output_prefix, "_raw_count_matrix.csv"))
  prop_csv   <- file.path(TABLES_DIR, paste0(output_prefix, "_normalized_proportion_matrix.csv"))

  ggsave(count_file, count_plot, width = 10, height = 8, dpi = 300)
  ggsave(prop_file,  prop_plot,  width = 10, height = 8, dpi = 300)
  write.csv(count_matrix,      count_csv, row.names = TRUE)
  write.csv(proportion_matrix, prop_csv,  row.names = TRUE)

  cat("  Saved:", basename(count_file), "\n")
  cat("  Saved:", basename(prop_file), "\n")
  cat("  Saved:", basename(count_csv), "\n")
  cat("  Saved:", basename(prop_csv), "\n")

  list(count = count_plot, proportion = prop_plot)
}

# =============================================================================
# MAIN
# =============================================================================
main <- function() {
  cat("============================================================\n")
  cat("Cluster vs NT — Network Connectivity Matrix\n")
  cat("============================================================\n")
  cat("Output dir:", OUTPUT_DIR, "\n\n")

  if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
  if (!dir.exists(TABLES_DIR)) dir.create(TABLES_DIR, recursive = TRUE)

  atlas_data <- load_atlas_info()
  total_possible <- calculate_total_possible_connections(atlas_data)

  clusters <- c("C1", "C2", "C3")
  types    <- c("hyper", "hypo")

  for (cl in clusters) {
    for (tp in types) {
      cat("\nProcessing", cl, "vs NT —", tp, "...\n")
      conn_data <- load_cluster_data(cl, tp)
      if (is.null(conn_data)) {
        cat("  No data — skipping\n")
        next
      }
      count_mat <- create_network_matrix(conn_data, atlas_data)
      if (is.null(count_mat)) {
        cat("  No valid network pairs — skipping\n")
        next
      }

      # Align proportion matrix dimensions with count matrix
      shared_nets <- intersect(rownames(count_mat), rownames(total_possible))
      count_sub   <- count_mat[shared_nets, shared_nets, drop = FALSE]
      total_sub   <- total_possible[shared_nets, shared_nets, drop = FALSE]
      prop_mat    <- count_sub / total_sub
      prop_mat[is.nan(prop_mat) | is.infinite(prop_mat)] <- 0

      title_suffix  <- paste0("Cluster ", cl, " vs NT — ",
                               ifelse(tp == "hyper", "Hyperconnectivity", "Hypoconnectivity"))
      output_prefix <- paste0("cluster_", tolower(cl), "_vs_nt_", tp)

      create_matrix_plots(count_sub, prop_mat, title_suffix, output_prefix)
    }
  }

  cat("\n============================================================\n")
  cat("Network Matrix Analysis Complete\n")
  cat("============================================================\n")
  invisible(TRUE)
}

if (!interactive()) {
  main()
}
