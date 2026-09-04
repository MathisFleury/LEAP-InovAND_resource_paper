# =============================================================================
# Supplementary table — sample-size completeness across variable combinations
# =============================================================================
# Companion to 2_pca_features.R: renders the PCA feature-selection sample-size
# table it already computes (sample_sizes_by_variable_combination.csv) as a
# publication-ready gt table (CSV + PDF). Table structure/style mirrors
# eeg_mri-pipeline/analysis/figures_papers/feature_selection_rationale.R's
# "Sample Size Completeness" table (cols_label + fmt_number + bold-max row).
#
# Run:  Rscript 3_sample_size_table.R [curated]
#   curated  → write *_curated.{csv,pdf} (reads the flat outputs/tables/,
#              matching 2_pca_features.R's default output location)
# =============================================================================

library(readr)
library(dplyr)
library(gt)

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
SCRIPT_DIR <- if (length(script_path)) dirname(normalizePath(script_path)) else getwd()

user_args <- commandArgs(trailingOnly = TRUE)
CURATED <- any(user_args %in% c("curated", "--curated"))

TABLES_DIR <- normalizePath(file.path(SCRIPT_DIR, "..", "outputs", "tables"))
OUT_SUFFIX <- if (CURATED) "_curated" else ""

sizes <- read_csv(file.path(TABLES_DIR, "sample_sizes_by_variable_combination.csv"),
                   show_col_types = FALSE)
pool_n <- sizes$Pool_N[1]

display <- sizes %>%
  transmute(
    Variable_Combination = Variable_Combination,
    N_Variables = N_Variables,
    N_Participants = Sample_Size,
    Percentage_Complete = Percent_of_Pool
  ) %>%
  arrange(desc(N_Participants))

write.csv(display, file.path(TABLES_DIR, paste0("sample_size_table", OUT_SUFFIX, ".csv")),
          row.names = FALSE)

sample_size_table <- display %>%
  gt() %>%
  tab_header(
    title = md("**Sample size completeness across variable combinations**"),
    subtitle = md(sprintf(
      "%sNumber of participants with complete data for each candidate clinical-variable combination, out of the full participant pool (N = %d). Higher sample sizes provide greater statistical power for cluster validation.",
      if (CURATED) "Curated cohort. " else "", pool_n))
  ) %>%
  cols_label(
    Variable_Combination = "Variable Combination",
    N_Variables = "N Variables",
    N_Participants = "N Participants",
    Percentage_Complete = "% of Pool"
  ) %>%
  fmt_number(columns = c(N_Variables, N_Participants), decimals = 0) %>%
  fmt_number(columns = Percentage_Complete, decimals = 1) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = N_Participants,
      rows = N_Participants == max(N_Participants)
    )
  ) %>%
  tab_options(table.width = pct(100), table.font.size = px(10))

gtsave(sample_size_table, file.path(TABLES_DIR, paste0("sample_size_table", OUT_SUFFIX, ".pdf")))

cat("Saved:\n  ", file.path(TABLES_DIR, paste0("sample_size_table", OUT_SUFFIX, ".csv")), "\n  ",
    file.path(TABLES_DIR, paste0("sample_size_table", OUT_SUFFIX, ".pdf")), "\n")
print(display)
