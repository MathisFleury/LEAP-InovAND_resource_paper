# =============================================================================
# Anatomical MRI Autism vs TD Brain Visualizations using ggseg
# =============================================================================
#
# Adapted from eeg_mri-pipeline/analysis/figures_papers/mri_autism_td_analysis/plot_anatomical_mri_brain_visualizations.R
# =============================================================================

library(plyr)
library(dplyr)
library(ggplot2)
library(colorspace)
library(ggseg)

# Shared theme for all ggseg brain maps
brain_map_theme <- theme_minimal() +
  theme(
    plot.title       = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle    = element_text(size = 14, hjust = 0.5),
    legend.title     = element_text(size = 12, face = "bold"),
    legend.text      = element_text(size = 10),
    legend.position  = "bottom",
    panel.grid       = element_blank(),
    axis.text        = element_blank(),
    axis.ticks       = element_blank(),
    strip.text       = element_text(size = 11, face = "bold"),
    strip.background = element_rect(fill = "lightgray", color = "black")
  )

# Paths
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
if (length(script_path) == 0) {
  script_dir <- getwd()
} else {
  script_dir <- dirname(normalizePath(script_path))
}
section_dir <- dirname(script_dir)
input_dir  <- file.path(section_dir, "outputs", "figures", "r_input_files")
output_dir <- file.path(section_dir, "outputs", "figures")

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

transform_region_column <- function(df) {
  df <- df %>%
    plyr::mutate(label = case_when(
      endsWith(label, "_left")  ~ paste0("lh_", sub("_left$", "", label)),
      endsWith(label, "_right") ~ paste0("rh_", sub("_right$", "", label)),
      TRUE ~ label
    ))
  return(df)
}

rename_subcortical_labels <- function(df) {
  df <- df %>%
    dplyr::mutate(label = case_when(
      label == "left-thalamus"                    ~ "Left-Thalamus-Proper",
      label == "right-thalamus"                   ~ "Right-Thalamus-Proper",
      label == "left-caudate"                     ~ "Left-Caudate",
      label == "right-caudate"                    ~ "Right-Caudate",
      label == "left-putamen"                     ~ "Left-Putamen",
      label == "right-putamen"                    ~ "Right-Putamen",
      label == "left-pallidum"                    ~ "Left-Pallidum",
      label == "right-pallidum"                   ~ "Right-Pallidum",
      label == "left-hippocampus"                 ~ "Left-Hippocampus",
      label == "right-hippocampus"                ~ "Right-Hippocampus",
      label == "left-amygdala"                    ~ "Left-Amygdala",
      label == "right-amygdala"                   ~ "Right-Amygdala",
      label == "left-accumbens"                   ~ "Left-Accumbens",
      label == "right-accumbens"                  ~ "Right-Accumbens",
      label == "left-ventraldc"                   ~ "Left-VentralDC",
      label == "right-ventraldc"                  ~ "Right-VentralDC",
      label == "left-choroid-plexus"              ~ "Left-Choroid-Plexus",
      label == "right-choroid-plexus"             ~ "Right-Choroid-Plexus",
      label == "cc_posterior"                     ~ "CC_Posterior",
      label == "cc_mid_posterior"                 ~ "CC_Mid_Posterior",
      label == "cc_central"                       ~ "CC_Central",
      label == "cc_mid_anterior"                  ~ "CC_Mid_Anterior",
      label == "cc_anterior"                      ~ "CC_Anterior",
      label == "right-cerebellum-white-matter"    ~ "Right-Cerebellum-White-Matter",
      label == "left-cerebellum-white-matter"     ~ "Left-Cerebellum-White-Matter",
      label == "right-cerebellum-cortex"          ~ "Right-Cerebellum-Cortex",
      label == "left-cerebellum-cortex"           ~ "Left-Cerebellum-Cortex",
      label == "3rd-ventricle"                    ~ "x3rd-ventricle",
      label == "4th-ventricle"                    ~ "x4th-ventricle",
      label == "brain-stem"                       ~ "Brain-Stem",
      label == "left-lateral-ventricle"           ~ "Left-Lateral-Ventricle",
      label == "right-lateral-ventricle"          ~ "Right-Lateral-Ventricle",
      TRUE ~ label
    ))
  return(df)
}

load_t_stat_data <- function(metric, atlas_type) {
  file_path <- file.path(input_dir, paste0("t_stat_anat_", atlas_type, "_", metric, "_mri_autism_vs_control.csv"))
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(NULL)
  }
  data <- read.csv(file_path, header = FALSE, col.names = c("label", "t_stat", "p_val", "cohens_d", "p_fdr"),
                   stringsAsFactors = FALSE)
  data <- data %>%
    dplyr::mutate(fdr_corrected = p_fdr, significant = fdr_corrected < 0.05)
  return(data)
}

# =============================================================================
# CORTICAL THICKNESS
# =============================================================================
print("Processing cortical thickness...")
cortical_thickness_data <- load_t_stat_data("thickness", "dk")
if (!is.null(cortical_thickness_data)) {
  cortical_thickness_data <- transform_region_column(cortical_thickness_data)
  cortical_thickness_data <- cortical_thickness_data %>% plyr::mutate(label = tolower(label))
  t_stat_max <- max(abs(range(cortical_thickness_data$t_stat, na.rm = TRUE)))

  p1 <- cortical_thickness_data %>%
    ggseg(mapping = aes(fill = t_stat), atlas = dk) +
    scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0, rev = TRUE,
                                     limits = c(-t_stat_max, t_stat_max), name = "t-statistic") +
    brain_map_theme +
    labs(title = "Cortical Thickness: Autism vs TD", subtitle = "t-statistics (FDR-corrected p < 0.05 highlighted)")
  ggsave(file.path(output_dir, "cortical_thickness_autism_control.pdf"), plot = p1, width = 10, height = 10)
  write.csv(cortical_thickness_data, file.path(output_dir, "cortical_thickness_autism_control_all.csv"), row.names = FALSE)

  cortical_significant <- cortical_thickness_data %>% dplyr::filter(fdr_corrected < 0.05)
  cortical_thickness_data <- cortical_thickness_data %>%
    plyr::mutate(outline_color = ifelse(label %in% cortical_significant$label, "black", "white"),
                 outline_size  = ifelse(label %in% cortical_significant$label, 1.2, 0.2))

  p2 <- cortical_thickness_data %>%
    ggseg(mapping = aes(fill = t_stat, col = outline_color, size = outline_size), atlas = "dk", position = "stacked") +
    scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0, rev = TRUE,
                                     limits = c(-t_stat_max, t_stat_max), name = "t-statistic") +
    scale_colour_identity() + scale_size_identity() +
    brain_map_theme +
    labs(title = "Cortical Thickness: Autism vs TD", subtitle = "t-statistics (FDR-corrected p < 0.05 highlighted with thick black border)")
  ggsave(file.path(output_dir, "cortical_thickness_autism_control_fdr.pdf"), plot = p2, width = 10, height = 10)
  write.csv(cortical_thickness_data, file.path(output_dir, "cortical_thickness_autism_control_fdr.csv"), row.names = FALSE)
  print("Saved: cortical_thickness_autism_control.pdf + cortical_thickness_autism_control_fdr.pdf")
}

# =============================================================================
# SURFACE AREA
# =============================================================================
print("Processing surface area...")
cortical_area_data <- load_t_stat_data("area", "dk")
if (!is.null(cortical_area_data)) {
  cortical_area_data <- transform_region_column(cortical_area_data)
  cortical_area_data <- cortical_area_data %>% plyr::mutate(label = tolower(label))
  t_stat_max <- max(abs(range(cortical_area_data$t_stat, na.rm = TRUE)))

  p3 <- cortical_area_data %>%
    ggseg(mapping = aes(fill = t_stat), atlas = dk) +
    scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0, rev = TRUE,
                                     limits = c(-t_stat_max, t_stat_max), name = "t-statistic") +
    brain_map_theme +
    labs(title = "Surface Area: Autism vs TD", subtitle = "t-statistics (FDR-corrected p < 0.05 highlighted)")
  ggsave(file.path(output_dir, "surface_area_autism_control.pdf"), plot = p3, width = 10, height = 10)
  write.csv(cortical_area_data, file.path(output_dir, "surface_area_autism_control_all.csv"), row.names = FALSE)

  cortical_area_significant <- cortical_area_data %>% dplyr::filter(fdr_corrected < 0.05)
  cortical_area_data <- cortical_area_data %>%
    plyr::mutate(outline_color = ifelse(label %in% cortical_area_significant$label, "black", "white"),
                 outline_size  = ifelse(label %in% cortical_area_significant$label, 2, 0.2))

  p4 <- cortical_area_data %>%
    ggseg(mapping = aes(fill = t_stat, col = outline_color, size = outline_size), atlas = "dk", position = "stacked") +
    scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0, rev = TRUE,
                                     limits = c(-t_stat_max, t_stat_max), name = "t-statistic") +
    scale_colour_identity() + scale_size_identity() +
    brain_map_theme +
    labs(title = "Surface Area: Autism vs TD", subtitle = "t-statistics (FDR-corrected p < 0.05 highlighted with thick black border)")
  ggsave(file.path(output_dir, "surface_area_autism_control_fdr.pdf"), plot = p4, width = 10, height = 10)
  write.csv(cortical_area_data, file.path(output_dir, "surface_area_autism_control_fdr.csv"), row.names = FALSE)
  print("Saved: surface_area_autism_control.pdf + surface_area_autism_control_fdr.pdf")
}

# =============================================================================
# SUBCORTICAL VOLUMES
# =============================================================================
print("Processing subcortical volumes...")
subcortical_data <- load_t_stat_data("volume", "aseg")
if (!is.null(subcortical_data)) {
  subcortical_data <- subcortical_data %>% dplyr::mutate(label = tolower(label))
  subcortical_data <- rename_subcortical_labels(subcortical_data)
  t_stat_max <- max(abs(range(subcortical_data$t_stat, na.rm = TRUE)))

  p5 <- subcortical_data %>%
    ggseg(mapping = aes(fill = t_stat), atlas = aseg) +
    scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0, rev = TRUE,
                                     limits = c(-t_stat_max, t_stat_max), name = "t-statistic") +
    brain_map_theme +
    labs(title = "Subcortical Volumes: Autism vs TD", subtitle = "t-statistics (FDR-corrected p < 0.05 highlighted)")
  ggsave(file.path(output_dir, "subcortical_volume_autism_control.pdf"), plot = p5, width = 10, height = 10)
  write.csv(subcortical_data, file.path(output_dir, "subcortical_volume_autism_control_all.csv"), row.names = FALSE)

  sub_cortical_significant <- subcortical_data %>% dplyr::filter(fdr_corrected < 0.05)
  subcortical_data <- subcortical_data %>%
    dplyr::mutate(outline_color = ifelse(label %in% sub_cortical_significant$label, "black", "white"),
                  outline_size  = ifelse(label %in% sub_cortical_significant$label, 1.7, 0.007))

  p6 <- subcortical_data %>%
    ggseg(mapping = aes(fill = t_stat, colour = outline_color, size = outline_size), atlas = aseg) +
    scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0, rev = TRUE,
                                     limits = c(-t_stat_max, t_stat_max), name = "t-statistic") +
    scale_colour_identity() + scale_size_identity() +
    brain_map_theme +
    labs(title = "Subcortical Volumes: Autism vs TD", subtitle = "t-statistics (FDR-corrected p < 0.05 highlighted with thick black border)")
  ggsave(file.path(output_dir, "subcortical_volume_autism_control_fdr.pdf"), plot = p6, width = 10, height = 10)
  write.csv(subcortical_data, file.path(output_dir, "subcortical_volume_autism_control_fdr.csv"), row.names = FALSE)
  print("Saved: subcortical_volume_autism_control.pdf + subcortical_volume_autism_control_fdr.pdf")
}

cat("============================================================\n")
cat("Brain visualizations complete!\n")
cat(paste("Output directory:", output_dir, "\n"))
cat("============================================================\n")
