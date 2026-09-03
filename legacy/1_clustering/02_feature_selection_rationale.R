# =============================================================================
# 02 - Feature Selection Rationale
# =============================================================================
# Generates figures and tables supporting the selection of IQ and SRS for
# clustering: (1) sample size completeness across variable combinations,
# (2) PCA loadings showing IQ and SRS as primary contributors.
#
# Paper: IQ and SRS maximize sample size (n=1,026 vs ~500 with VABS) and
# contribute most to PC1/PC2; VABS contributes to PC1 but minimally to PC2.
# =============================================================================

# --- Configuration ---
DATA_PATH <- Sys.getenv("LEAP_INOVAND_DATA", unset = NA)
if (is.na(DATA_PATH) || DATA_PATH == "") {
  DATA_PATH <- file.path(getwd(), "..", "..", "..", "imaging2genet", "0_input", "dataframes")
}
INDIVIDUALS_METRICS <- file.path(DATA_PATH, "individuals_metrics.tsv")
OUTPUT_BASE <- normalizePath(file.path(getwd(), "outputs"), mustWork = FALSE)
FIGURES_DIR <- file.path(OUTPUT_BASE, "figures")
TABLES_DIR <- file.path(OUTPUT_BASE, "tables")
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)

# --- Libraries ---
library(readr)
library(dplyr)
library(ggplot2)
library(FactoMineR)
library(factoextra)
library(gt)

# --- Data preparation ---
cat("Loading data from:", INDIVIDUALS_METRICS, "\n")
df <- readr::read_tsv(INDIVIDUALS_METRICS)
df <- df[!duplicated(df$ID), ]
df[df == 999] <- NA
df <- df[df$Relation_to_proposant == "participant", ]
df$total_IQ <- ifelse(is.na(df$total_IQ), df$performance_IQ, df$total_IQ)

clinical_columns <- c(
  "performance_IQ", "SRS_tscore", "ssp_total", "RBS-R_total", "verbal_IQ",
  "vabsdscoresc_dss", "vabsdscoresd_dss", "vabsdscoress_dss", "vabsabcabc_standard", "total_IQ"
)
clinical_columns <- intersect(clinical_columns, colnames(df))

# Transform (as in original pipeline)
df$ssp_total <- log1p(190 + 38 - df$ssp_total)
df$"RBS-R_total" <- log1p(df$"RBS-R_total")

# --- Table 1: Sample size completeness ---
variable_combinations <- list(
  "IQ + SRS" = c("total_IQ", "SRS_tscore"),
  "IQ + VABS Composite" = c("total_IQ", "vabsabcabc_standard"),
  "IQ + SRS + VABS Composite" = c("total_IQ", "SRS_tscore", "vabsabcabc_standard"),
  "IQ + SRS + SSP" = c("total_IQ", "SRS_tscore", "ssp_total"),
  "IQ + SRS + RBS-R" = c("total_IQ", "SRS_tscore", "RBS-R_total"),
  "IQ + SRS + VABS + SSP" = c("total_IQ", "SRS_tscore", "vabsabcabc_standard", "ssp_total"),
  "IQ + SRS + VABS + RBS-R" = c("total_IQ", "SRS_tscore", "vabsabcabc_standard", "RBS-R_total"),
  "IQ + SRS + VABS + SSP + RBS-R" = c("total_IQ", "SRS_tscore", "vabsabcabc_standard", "ssp_total", "RBS-R_total")
)

total_participants <- nrow(df)
completeness_results <- data.frame(
  Combination = character(),
  N_Variables = integer(),
  N_Participants = integer(),
  Percentage_Complete = numeric(),
  stringsAsFactors = FALSE
)
for (i in seq_along(variable_combinations)) {
  vars <- variable_combinations[[i]][variable_combinations[[i]] %in% colnames(df)]
  if (length(vars) > 0) {
    n_complete <- sum(complete.cases(df[, vars, drop = FALSE]))
    completeness_results <- rbind(completeness_results, data.frame(
      Combination = names(variable_combinations)[i],
      N_Variables = length(vars),
      N_Participants = n_complete,
      Percentage_Complete = round(100 * n_complete / total_participants, 1),
      stringsAsFactors = FALSE
    ))
  }
}
completeness_results <- completeness_results[order(-completeness_results$N_Participants), ]

write.csv(completeness_results, file.path(TABLES_DIR, "sample_size_completeness.csv"), row.names = FALSE)
cat("Sample size completeness saved.\n")

# Bar plot
p1 <- ggplot(completeness_results, aes(x = reorder(Combination, N_Participants), y = N_Participants)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
  geom_text(aes(label = paste0(N_Participants, "\n(", Percentage_Complete, "%)")),
    hjust = 0.5, vjust = -0.3, size = 3) +
  coord_flip() +
  labs(
    title = "Sample Size Completeness Across Variable Combinations",
    x = "Variable Combination",
    y = "Number of Participants with Complete Data"
  ) +
  theme_minimal()
ggsave(file.path(FIGURES_DIR, "sample_size_completeness.pdf"), plot = p1, width = 10, height = 8)

# --- Table 2: PCA results ---
df_pca <- df[complete.cases(df[, clinical_columns]), ]
cat("PCA sample size:", nrow(df_pca), "\n")

res.pca <- PCA(df_pca[, clinical_columns], scale.unit = TRUE, graph = FALSE)
eig.val <- get_eigenvalue(res.pca)

pca_variance_table <- data.frame(
  Component = paste0("PC", seq_len(nrow(eig.val))),
  Eigenvalue = round(eig.val[, "eigenvalue"], 2),
  Variance_Explained_Percent = round(eig.val[, "variance.percent"], 1),
  Cumulative_Variance_Percent = round(eig.val[, "cumulative.variance.percent"], 1),
  stringsAsFactors = FALSE
)
write.csv(pca_variance_table, file.path(TABLES_DIR, "pca_variance_explained.csv"), row.names = FALSE)

var <- get_pca_var(res.pca)
pca_loadings_table <- data.frame(
  Variable = rownames(var$contrib),
  Contribution_PC1_Percent = round(var$contrib[, 1], 2),
  Correlation_PC1 = round(var$cor[, 1], 3),
  Contribution_PC2_Percent = round(var$contrib[, 2], 2),
  Correlation_PC2 = round(var$cor[, 2], 3),
  stringsAsFactors = FALSE
)
pca_loadings_table <- pca_loadings_table[order(-pca_loadings_table$Contribution_PC1_Percent), ]
write.csv(pca_loadings_table, file.path(TABLES_DIR, "pca_variable_loadings.csv"), row.names = FALSE)

# --- Figures: PCA ---
p2 <- fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 60),
  title = "PCA Scree Plot - Variance Explained by Each Component") + theme_minimal()
ggsave(file.path(FIGURES_DIR, "pca_scree_plot.pdf"), plot = p2, width = 8, height = 6)

p3 <- fviz_contrib(res.pca, choice = "var", axes = 1, top = 10,
  title = paste0("Contribution of Variables to PC1 (", round(eig.val[1, 2], 1), "% variance)")) + theme_minimal()
ggsave(file.path(FIGURES_DIR, "pca_contrib_pc1.pdf"), plot = p3, width = 8, height = 6)

p4 <- fviz_contrib(res.pca, choice = "var", axes = 2, top = 10,
  title = paste0("Contribution of Variables to PC2 (", round(eig.val[2, 2], 1), "% variance)")) + theme_minimal()
ggsave(file.path(FIGURES_DIR, "pca_contrib_pc2.pdf"), plot = p4, width = 8, height = 6)

p5 <- fviz_pca_var(res.pca, col.var = "contrib",
  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
  title = "PCA Variable Plot - Contributions to First Two Components",
  repel = TRUE) + theme_minimal()
ggsave(file.path(FIGURES_DIR, "pca_variables_plot.pdf"), plot = p5, width = 10, height = 8)

cat("Done. Outputs in:", OUTPUT_BASE, "\n")
