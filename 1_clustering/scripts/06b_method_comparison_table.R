# =============================================================================
# Supplementary table — clustering method comparison (k-means vs GMM vs Ward)
# =============================================================================
# Companion to 06_method_comparison.py / 05_cluster_stability.py: assembles the
# per-method internal-validity and resampling-stability metrics they wrote into
# outputs/tables/ and renders a publication-ready gt table (CSV + PDF), with the
# best value in each column highlighted in bold.
#
#   Internal validity : silhouette, Davies-Bouldin, Calinski-Harabasz
#   Stability         : bootstrap ARI (80% subsample), % subjects with modal
#                       assignment >= 0.90, mean per-cluster Jaccard
#
# Run:  Rscript 06b_method_comparison_table.R [curated]
#   curated  → read outputs/curated/tables/ and write *_curated.{csv,pdf}
# =============================================================================

library(readr)
library(dplyr)
library(gt)

# --- args / paths ------------------------------------------------------------
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
SCRIPT_DIR <- if (length(script_path)) dirname(normalizePath(script_path)) else getwd()

user_args <- commandArgs(trailingOnly = TRUE)
CURATED <- any(user_args %in% c("curated", "--curated"))

# curated is the priority regime; its metrics live under outputs/curated/.
TABLES_DIR <- normalizePath(file.path(
  SCRIPT_DIR, "..", "outputs", if (CURATED) "curated" else ".", "tables"))
OUT_SUFFIX <- if (CURATED) "_curated" else ""

BOOT_FRACTION <- 0.8   # subsample fraction reported for bootstrap ARI

# --- read the metric tables 06/05 already produced ---------------------------
internal <- read_csv(file.path(TABLES_DIR, "method_comparison_internal_indices.csv"),
                     show_col_types = FALSE)
stab     <- read_csv(file.path(TABLES_DIR, "stability_summary_all_methods.csv"),
                     show_col_types = FALSE)
jacc     <- read_csv(file.path(TABLES_DIR, "method_comparison_jaccard.csv"),
                     show_col_types = FALSE)

# bootstrap ARI at the chosen subsample fraction, per method
boot <- lapply(c("kmeans", "ward", "gmm"), function(m) {
  f <- file.path(TABLES_DIR, sprintf("stability_subsampling_%s.csv", m))
  d <- read_csv(f, show_col_types = FALSE)
  r <- d[which.min(abs(d$fraction - BOOT_FRACTION)), ]
  data.frame(method = m, ARI_mean = r$ARI_mean, ARI_sd = r$ARI_sd)
}) %>% bind_rows()

jacc_mean <- jacc %>% group_by(method) %>%
  summarise(jaccard = mean(jaccard), .groups = "drop")

# --- assemble one row per method ---------------------------------------------
labels <- c(kmeans = "K-means", gmm = "GMM", ward = "Ward")
order  <- c("kmeans", "gmm", "ward")   # primary method first

tab <- data.frame(method = order) %>%
  left_join(internal, by = "method") %>%
  left_join(boot,     by = "method") %>%
  left_join(stab %>% select(method, modal_prob_ge_0.90), by = "method") %>%
  left_join(jacc_mean, by = "method") %>%
  transmute(
    Method            = labels[method],
    Silhouette        = round(silhouette, 3),
    `Davies-Bouldin`  = round(davies_bouldin, 3),
    `Calinski-Harabasz` = round(calinski_harabasz, 0),
    Bootstrap_ARI_mean = ARI_mean,                      # numeric, for bolding
    Bootstrap_ARI_sd   = ARI_sd,
    `Subjects stable (%)` = round(100 * modal_prob_ge_0.90, 1),
    `Mean Jaccard`     = round(jaccard, 3)
  )

# display string "mean (sd)" for the bootstrap ARI column
tab$`Bootstrap ARI` <- sprintf("%.3f (%.3f)", tab$Bootstrap_ARI_mean, tab$Bootstrap_ARI_sd)

display <- tab %>%
  select(Method, Silhouette, `Davies-Bouldin`, `Calinski-Harabasz`,
         `Bootstrap ARI`, `Subjects stable (%)`, `Mean Jaccard`)

# --- best value per column (bold) --------------------------------------------
# higher is better for all except Davies-Bouldin
best <- list(
  Silhouette            = which.max(tab$Silhouette),
  `Davies-Bouldin`      = which.min(tab$`Davies-Bouldin`),
  `Calinski-Harabasz`   = which.max(tab$`Calinski-Harabasz`),
  `Bootstrap ARI`       = which.max(tab$Bootstrap_ARI_mean),
  `Subjects stable (%)` = which.max(tab$`Subjects stable (%)`),
  `Mean Jaccard`        = which.max(tab$`Mean Jaccard`)
)

# --- CSV ---------------------------------------------------------------------
write.csv(display, file.path(TABLES_DIR, paste0("method_comparison_table", OUT_SUFFIX, ".csv")),
          row.names = FALSE)

# --- gt table ----------------------------------------------------------------
gt_table <- display %>%
  gt() %>%
  tab_header(
    title = md("**Comparison of clustering algorithms: internal validity and stability**"),
    subtitle = md(sprintf(
      "%sK-means (primary) vs. GMM and Ward (sensitivity), on standardised IQ x SRS (k = 3). Higher is better for all indices except Davies-Bouldin (lower is better). Bootstrap ARI and per-cluster Jaccard are means over %d%%-subsample resamples; \"subjects stable\" = %% with modal cluster assignment >= 0.90. Best value per column in **bold**.",
      if (CURATED) "Curated cohort. " else "", round(BOOT_FRACTION * 100)))
  ) %>%
  tab_options(table.font.size = px(11), column_labels.font.weight = "bold")

for (col in names(best)) {
  gt_table <- gt_table %>%
    tab_style(style = cell_text(weight = "bold"),
              locations = cells_body(columns = all_of(col), rows = best[[col]]))
}

gtsave(gt_table, file.path(TABLES_DIR, paste0("method_comparison_table", OUT_SUFFIX, ".pdf")))

cat("Saved:\n  ", file.path(TABLES_DIR, paste0("method_comparison_table", OUT_SUFFIX, ".csv")), "\n  ",
    file.path(TABLES_DIR, paste0("method_comparison_table", OUT_SUFFIX, ".pdf")), "\n")
print(display)
