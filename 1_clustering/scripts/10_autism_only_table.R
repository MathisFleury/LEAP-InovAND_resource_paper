# =============================================================================
# Supplementary table — autism-only clustering profile & concordance (R4.6)
# =============================================================================
# Companion to 07_autism_only_clustering.py: assembles the autism-only cluster
# profile (IQ/SRS-2), its crosstab against the full-cohort assignment, and the
# concordance/silhouette numbers into one publication-ready gt table (CSV + PDF).
# Same rendering convention as 06b_method_comparison_table.R.
#
# Run:  Rscript 07b_autism_only_table.R
# =============================================================================

library(readr)
library(dplyr)
library(gt)

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
SCRIPT_DIR <- if (length(script_path)) dirname(normalizePath(script_path)) else getwd()
TABLES_DIR <- normalizePath(file.path(SCRIPT_DIR, "..", "outputs", "curated", "tables"))

profile     <- read_csv(file.path(TABLES_DIR, "autism_only_cluster_profile.csv"), show_col_types = FALSE)
crosstab    <- read_csv(file.path(TABLES_DIR, "autism_only_vs_fullcohort_crosstab.csv"), show_col_types = FALSE)
concordance <- read_csv(file.path(TABLES_DIR, "autism_only_concordance.csv"), show_col_types = FALSE)
silhouette  <- read_csv(file.path(TABLES_DIR, "autism_only_k_selection_silhouette.csv"), show_col_types = FALSE)

display <- profile %>%
  transmute(
    Cluster = ClusterAutOnly,
    n = n,
    `Full-scale IQ` = sprintf("%.1f (%.1f)", IQ_mean, IQ_sd),
    `SRS-2 t-score` = sprintf("%.1f (%.1f)", SRS_mean, SRS_sd)
  ) %>%
  left_join(crosstab, by = c("Cluster" = "ClusterAutOnly"))

write.csv(display, file.path(TABLES_DIR, "autism_only_summary_table.csv"), row.names = FALSE)

# modal (largest) full-cohort overlap per row, for bolding
modal_col <- apply(display[, c("C1", "C2", "C3")], 1, function(r) names(r)[which.max(r)])
sil_txt <- paste(sprintf("k=%d: %.3f", silhouette$k, silhouette$silhouette), collapse = ", ")

gt_table <- display %>%
  gt() %>%
  tab_header(
    title = md("**Autism-only clustering profile and concordance with the full-cohort partition**"),
    subtitle = md(sprintf(
      "Curated cohort, autistic participants only (n = %d), k = 3 (Reviewer #4.6). C1-C3 columns show how many members of each autism-only cluster were assigned to each full-cohort cluster; modal (largest) overlap per row in **bold**. Concordance with the full-cohort assignment: ARI = %.2f, NMI = %.2f. Silhouette width by k: %s.",
      concordance$n_autistic, concordance$ari, concordance$nmi, sil_txt))
  ) %>%
  tab_spanner(label = "Full-cohort assignment", columns = c(C1, C2, C3)) %>%
  tab_options(table.font.size = px(11), column_labels.font.weight = "bold")

for (i in seq_len(nrow(display))) {
  gt_table <- gt_table %>%
    tab_style(style = cell_text(weight = "bold"),
              locations = cells_body(columns = all_of(modal_col[i]), rows = i))
}

gtsave(gt_table, file.path(TABLES_DIR, "autism_only_summary_table.pdf"))

cat("Saved:\n  ", file.path(TABLES_DIR, "autism_only_summary_table.csv"), "\n  ",
    file.path(TABLES_DIR, "autism_only_summary_table.pdf"), "\n")
print(display)
