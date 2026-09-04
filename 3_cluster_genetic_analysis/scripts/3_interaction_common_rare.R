#!/usr/bin/env Rscript
# =============================================================================
# Common (PGS) x rare (carrier) interaction figures, by cluster / population.
# =============================================================================
# Port of eeg_mri-pipeline/analysis/figures_papers/figure_interaction_commun_rare.R
# into this repo, with these fixes:
#   - input = gnomad v4 interaction table (current priority regime),
#     assembled by 2_interaction_v4_inputs.py, not the Ward frozen df
#   - plyr::mutate / ggpubr::mutate -> dplyr::mutate (those verbs don't exist)
#   - IDD/NT population spelling + v4 HCNDD "Dom rec-X" carrier column
#   - input/output parameterised via DF_CLUSTERS / CLUSTER_GENETIC_FIG_DIR env
#
# To reproduce the frozen v2 paper figure instead, set:
#   DF_CLUSTERS=.../df_clusters_complete_kmeans.csv HCNDD_CARRIER_COL=dellof_hcndddomv6_contraint_carrier
# (the superset palette + NT/TD handling below cover both regimes).
#
# Ancestry version chosen by env ANCESTRY (PAN default | EUR).
# Emits, under CLUSTER_GENETIC_FIG_DIR (default outputs/figures/), flat and
# <ANCESTRY>_-prefixed (matches 1_carrier_freq_or_clusters.py's convention):
#   <ANCESTRY>_carrier_stacked_barplots.pdf          % carriers per cluster x phenotype
#   <ANCESTRY>_pgs_by_carrier_type.pdf               PGS x cluster, carrier vs non, per gene set
#   <ANCESTRY>_pgs_by_carrier_type_Population1.pdf   PGS x population, carrier vs non
#   <ANCESTRY>_pgs_by_cluster_<gene set>.pdf         PGS x (NT,C1,C2,C3), carrier vs non
# read.csv (not readr) is deliberate: it make.names-mangles "a/b" PGS columns to
# "a.b", which the column names below rely on.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(forcats); library(tidyr); library(patchwork)
})
if (!requireNamespace("Hmisc", quietly = TRUE))
  stop("Hmisc is required for stat_summary(mean_cl_normal). Install it first.")

args <- commandArgs(trailingOnly = FALSE)
sp <- sub("--file=", "", args[grep("--file=", args)])
script_dir <- if (length(sp) == 0) getwd() else dirname(normalizePath(sp))
section_dir <- dirname(script_dir)

# Ancestry version: PAN (pan-ancestry PGS, whole cohort) or EUR (European-only PGS,
# EUR_ancestry==True). The assembler (2_interaction_v4_inputs.py) already applies
# the cohort filter, so this script does NOT re-filter on ancestry.
ANCESTRY <- Sys.getenv("ANCESTRY", "PAN")
DF_CLUSTERS <- Sys.getenv(
  "DF_CLUSTERS",
  file.path(section_dir, "outputs", "tables", sprintf("df_v4_interaction_inputs_%s.csv", ANCESTRY)))
output_dir <- Sys.getenv("CLUSTER_GENETIC_FIG_DIR",
                         file.path(section_dir, "outputs", "figures"))
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
# Flat, ANCESTRY-prefixed filenames (matches 1_carrier_freq_or_clusters.py's
# PAN_/EUR_ convention) instead of a PAN/EUR subfolder.
fig_path <- function(name) file.path(output_dir, paste0(ANCESTRY, "_", name))

cat("Ancestry:", ANCESTRY, "\nInput :", DF_CLUSTERS, "\nOutput:", output_dir, "\n")
df_genetic <- read.csv(DF_CLUSTERS)

# v4 HCNDD constrained "Dom rec-X" set; override for the v2 frozen figure.
hcndd_col <- Sys.getenv("HCNDD_CARRIER_COL", "dellof_hcndddom_xlinked_boyz_contraint_carrier")
carrier_types <- setNames(
  list("HCNDD Dom rec-X carriers", "SPARKS-FARI1 carriers", " EAGLE carriers",
       "SYNGO carriers", "CHROMATIN/EPI/TF carriers"),
  c(hcndd_col, "dellof_sparksfari1_contraint_carrier", "dellof_eagle_contraint_carrier",
    "dellof_syngo_contraint_carrier", "dellof_chromepitf_contraint_carrier"))

# Superset palette: v4 (IDD/NT) + v2 (ID/TD) keys, so either regime renders.
population_palette <- c(
  "NT" = "#C1C2BC", "TD" = "#C1C2BC",
  "Autism with IDD" = "#324095", "Autism with ID" = "#8991FA",
  "IDD" = "#D8A4CB", "ID" = "#D8A4CB",
  "Autism without IDD" = "#5CAEE1", "Autism without ID" = "#5CAEE1",
  "Relatives" = "#9AD5D3", "Autism to exclude" = "#D9D9D9"
)

###### % CARRIERS PER CLUSTER x PHENOTYPE ######
plots <- list()
for (carrier_col in names(carrier_types)) {

  df_carriers <- df_genetic %>%
    dplyr::filter(!is.na(Cluster), !!sym(carrier_col) == "True")

  total_counts <- df_genetic %>%
    dplyr::filter(!is.na(Cluster)) %>%
    dplyr::count(Cluster, name = "total_n")

  stacked_data <- df_carriers %>%
    dplyr::count(Cluster, Population1, name = "n_carriers") %>%
    left_join(total_counts, by = "Cluster") %>%
    dplyr::mutate(pct = (n_carriers / total_n) * 100) %>%
    dplyr::filter(Cluster != "")

  p <- ggplot(stacked_data, aes(x = Cluster, y = pct, fill = Population1)) +
    geom_bar(stat = "identity", color = "black", width = 0.7) +
    scale_fill_manual(values = population_palette) +
    labs(title = carrier_types[[carrier_col]], y = "% of Carriers", x = NULL, fill = "Phenotype") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold", size = 13),
          legend.position = "right") +
    ylim(0, 60)

  plots[[carrier_col]] <- p
}
combined_plot <- wrap_plots(plots, nrow = 1)
ggsave(fig_path("carrier_stacked_barplots.pdf"), plot = combined_plot, width = 20, height = 8)

###### PGS x CLUSTER, CARRIER vs NON ######
pgs_columns <- c("autism_grove2019", "adhd_demontis2023", "int_savage2018",
                 "anxiety_purves2020", "ptsd_nievergelt2024", "ukbiobank_watanabe2019.felt_loved")
pdf(fig_path("pgs_by_carrier_type.pdf"), width = 10, height = 8)
for (carrier_col in names(carrier_types)) {

  df_carrier_plot <- df_genetic %>%
    dplyr::filter(!is.na(Cluster), !is.na(!!sym(carrier_col))) %>%
    dplyr::mutate(Carrier = as.character(!!sym(carrier_col))) %>%
    dplyr::filter(Carrier %in% c("True", "False"))

  plots <- list()
  for (pgs in pgs_columns) {
    p <- ggplot(df_carrier_plot, aes(x = Cluster, y = .data[[pgs]], color = Carrier)) +
      stat_summary(fun.data = mean_cl_normal, geom = "pointrange",
                   position = position_dodge(width = 0.6), linewidth = 0.6) +
      stat_summary(fun = mean, geom = "point", shape = 21, fill = "white", size = 2.5,
                   position = position_dodge(width = 0.6)) +
      coord_cartesian(ylim = c(-1.2, 1.2)) +
      labs(title = paste(pgs), y = "PGS Score", x = "Cluster", color = "Carrier") +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(face = "bold", size = 14))
    plots[[pgs]] <- p + plot_annotation(title = carrier_types[[carrier_col]])
  }
  print(plots[[pgs_columns[1]]] + plots[[pgs_columns[2]]] + plots[[pgs_columns[3]]] +
        plots[[pgs_columns[4]]] + plots[[pgs_columns[5]]] + plots[[pgs_columns[6]]])
}
dev.off()

###### PGS x POPULATION, CARRIER vs NON ######
pgs_columns <- c("autism_grove2019", "int_savage2018")
pgs_labels <- c(autism_grove2019 = "PGS of Autism", int_savage2018 = "PGS of Intel.")
pdf_path <- fig_path("pgs_by_carrier_type_Population1.pdf")

pdf(pdf_path, width = 4, height = 8)
for (carrier_col in names(carrier_types)) {

  df_carrier_plot <- df_genetic %>%
    dplyr::filter(!is.na(PopulationS1), !is.na(!!sym(carrier_col))) %>%
    dplyr::mutate(Carrier = as.character(!!sym(carrier_col))) %>%
    dplyr::filter(Carrier %in% c("True", "False")) %>%
    dplyr::filter(PopulationS1 != "Autism to exclude") %>%
    dplyr::filter(Relation_to_proposant == "participant")

  plots <- list()
  for (pgs in pgs_columns) {
    p <- ggplot(df_carrier_plot, aes(x = PopulationS1, y = .data[[pgs]], color = Carrier)) +
      stat_summary(fun.data = mean_cl_normal, geom = "pointrange",
                   position = position_dodge(width = 0.6), linewidth = 0.6) +
      stat_summary(fun = mean, geom = "point", shape = 21, fill = "white", size = 2.5,
                   position = position_dodge(width = 0.6)) +
      coord_cartesian(ylim = c(-1.2, 1.2)) +
      labs(title = paste(pgs), y = "PGS Score", x = "Population1", color = "Carrier") +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(face = "bold", size = 14))
    plots[[pgs]] <- p + plot_annotation(title = carrier_types[[carrier_col]])
  }
  print(plots[[pgs_columns[1]]] / plots[[pgs_columns[2]]])
}
dev.off()

# same, titled per-gene-set annotation, y-axis = readable PGS label
pdf(pdf_path, width = 4, height = 8)
for (carrier_col in names(carrier_types)) {

  df_carrier_plot <- df_genetic %>%
    dplyr::filter(!is.na(PopulationS1), !is.na(.data[[carrier_col]]),
                  PopulationS1 != "Autism to exclude",
                  Relation_to_proposant == "participant") %>%
    dplyr::mutate(Carrier = as.character(.data[[carrier_col]])) %>%
    dplyr::filter(Carrier %in% c("True", "False"))

  plots <- list()
  for (pgs in pgs_columns) {
    p <- ggplot(df_carrier_plot, aes(x = PopulationS1, y = .data[[pgs]], color = Carrier)) +
      stat_summary(fun.data = mean_cl_normal, geom = "pointrange",
                   position = position_dodge(width = 0.6), linewidth = 0.6) +
      stat_summary(fun = mean, geom = "point", shape = 21, fill = "white", size = 2.5,
                   position = position_dodge(width = 0.6)) +
      coord_cartesian(ylim = c(-1.2, 1.2)) +
      labs(title = NULL, y = pgs_labels[[pgs]], x = "Population", color = "Carrier") +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    plots[[pgs]] <- p
  }
  final_plot <- (plots[[pgs_columns[1]]] / plots[[pgs_columns[2]]]) +
    plot_annotation(title = carrier_types[[carrier_col]]) &
    theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
  print(final_plot)
}
dev.off()

###### PGS x CLUSTER WITH NT ######
pgs_columns <- c("autism_grove2019", "adhd_demontis2023", "int_savage2018",
                 "anxiety_purves2020", "ptsd_nievergelt2024", "ukbiobank_watanabe2019.felt_loved")
pgs_labels <- list(
  autism_grove2019 = "Autism \n (Grove et al. 2019)",
  adhd_demontis2023 = "ADHD \n (Demontis et al. 2023)",
  int_savage2018 = "Intelligence \n (Savage et al. 2018)",
  anxiety_purves2020 = "Anxiety \n (Purves et al. 2020)",
  ptsd_nievergelt2024 = "PTSD \n (Nievergelt et al. 2024)",
  "ukbiobank_watanabe2019.felt_loved" = "Felt Loved as a child \n (Watanabe et al. 2019)"
)
for (carrier_col in names(carrier_types)) {

  df_carrier_plot <- df_genetic %>%
    dplyr::filter(!is.na(!!sym(carrier_col))) %>%
    dplyr::mutate(
      Carrier = as.character(!!sym(carrier_col)),
      Cluster_modified = dplyr::case_when(
        PopulationS1 == "Autism" & !is.na(Cluster) & Cluster %in% c("C1", "C2", "C3") ~ Cluster,
        PopulationS1 %in% c("NT", "TD") ~ "NT",
        TRUE ~ NA_character_)) %>%
    dplyr::filter(!is.na(Cluster_modified), Carrier %in% c("True", "False")) %>%
    dplyr::mutate(Cluster_modified = forcats::fct_relevel(Cluster_modified, "NT", "C1", "C2", "C3"))

  plots <- list()
  for (pgs in pgs_columns) {
    p <- ggplot(df_carrier_plot, aes(x = Cluster_modified, y = .data[[pgs]], color = Carrier)) +
      stat_summary(fun.data = mean_cl_normal, geom = "pointrange",
                   position = position_dodge(width = 0.6), linewidth = 0.6) +
      stat_summary(fun = mean, geom = "point", shape = 21, fill = "white", size = 2.5,
                   position = position_dodge(width = 0.6)) +
      coord_cartesian(ylim = c(-1.2, 1.2)) +
      labs(title = pgs_labels[[pgs]], y = "z-score", x = NULL, color = "Carrier") +
      theme_minimal(base_size = 12) +
      theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14),
            axis.title.y = element_text(size = 14),
            plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
    plots[[pgs]] <- p
  }
  combined_plot <- wrap_plots(plots, ncol = 3, nrow = 2) +
    plot_annotation(title = carrier_types[[carrier_col]],
                    theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5)))

  carrier_name_clean <- gsub("[^A-Za-z0-9]", "_", carrier_types[[carrier_col]])
  carrier_name_clean <- gsub("_+", "_", carrier_name_clean)
  carrier_name_clean <- gsub("^_|_$", "", carrier_name_clean)
  pdf_path <- fig_path(paste0("pgs_by_cluster_", carrier_name_clean, ".pdf"))
  pdf(pdf_path, width = 15, height = 10)
  print(combined_plot)
  dev.off()
  cat("Saved:", pdf_path, "\n")
}
