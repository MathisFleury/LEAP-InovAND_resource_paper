# =============================================================================
# 2 - PCA on Clinical Features
# =============================================================================
# Principal component analysis to identify features explaining most variance.
# Used to justify focusing on IQ and SRS for clustering (see paper Methods).
#
# Outputs: PCA scree plot, variable loadings plot, contribution plots, tables
# =============================================================================

# --- Configuration ---
# Curated (current, priority regime) by default -- reads the cohort table
# written by 1_load_cohort.py. No frozen/deprecated-data fallback here: the
# paper-reproduction run lives in the local (not publicly released) frozen
# pipeline, which points PCA_INPUT_CSV at its own reference table and
# invokes this same engine.
DEFAULT_INPUT_CSV <- file.path(getwd(), "..", "outputs", "tables", "cohort_curated.csv")
INDIVIDUALS_METRICS <- Sys.getenv("PCA_INPUT_CSV", unset = DEFAULT_INPUT_CSV)
if (!file.exists(INDIVIDUALS_METRICS)) {
  stop("Input not found: ", INDIVIDUALS_METRICS, "\nRun 1_load_cohort.py first",
       " (or set PCA_INPUT_CSV to point elsewhere).")
}
# FIGURES_DIR / TABLES_DIR default to the flat, unqualified outputs/ -- curated
# is the current pipeline, so that's simply its home. The frozen wrapper
# (legacy/1_clustering/1_pca_features_frozen.R) overrides both to keep its
# own figures/tables self-contained instead of colliding with these.
FIGURES_DIR <- Sys.getenv("PCA_FIGURES_DIR",
                          unset = normalizePath(file.path(getwd(), "..", "outputs", "figures"), mustWork = FALSE))
TABLES_DIR <- Sys.getenv("PCA_TABLES_DIR",
                         unset = normalizePath(file.path(getwd(), "..", "outputs", "tables"), mustWork = FALSE))
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLES_DIR, recursive = TRUE, showWarnings = FALSE)

# --- Libraries ---
library(readr)
library(FactoMineR)
library(factoextra)
library(corrplot)

# --- Data preparation ---
cat("Loading data from:", INDIVIDUALS_METRICS, "\n")
df <- if (grepl("\\.csv$", INDIVIDUALS_METRICS)) {
  readr::read_csv(INDIVIDUALS_METRICS, show_col_types = FALSE)
} else {
  readr::read_tsv(INDIVIDUALS_METRICS)
}
df <- df[!duplicated(df$ID), ]
df[df == 999] <- NA
relation_col <- if ("relation_to_proposant" %in% colnames(df)) "relation_to_proposant" else "Relation_to_proposant"
df <- df[df[[relation_col]] == "participant", ]
df$total_IQ <- ifelse(is.na(df$total_IQ), df$performance_IQ, df$total_IQ)

# Broad set of clinical measures (paper: IQ subscales, SRS-2, Vineland, RBS-R, SSP)
clinical_columns <- c(
  "performance_IQ", "verbal_IQ", "total_IQ",
  "SRS_tscore",
  "ssp_total", "RBS-R_total",
  "vabsdscoresc_dss", "vabsdscoresd_dss", "vabsdscoress_dss", "vabsabcabc_standard"
)
# Use LEAP columns if INOVAND columns missing
col_map <- c(
  "LEAP_t1_ssp_total" = "ssp_total",
  "LEAP_t1_vabsdscoresc_dss" = "vabsdscoresc_dss",
  "LEAP_t1_vabsdscoresd_dss" = "vabsdscoresd_dss",
  "LEAP_t1_vabsdscoress_dss" = "vabsdscoress_dss",
  "LEAP_t1_vabsabcabc_standard" = "vabsabcabc_standard"
)
for (lep in names(col_map)) {
  if (lep %in% colnames(df) && !col_map[lep] %in% colnames(df))
    df[[col_map[lep]]] <- df[[lep]]
}
clinical_columns <- intersect(clinical_columns, colnames(df))

# Exclude missing data and relatives (paper: n = 380 or 374 excluded)
df_pca <- df[complete.cases(df[, clinical_columns]), ]
cat("Sample size for PCA (complete cases, participants only):", nrow(df_pca), "\n")

# Standardize and PCA
res.pca <- PCA(df_pca[, clinical_columns], scale.unit = TRUE, graph = FALSE)
eig.val <- get_eigenvalue(res.pca)

# --- Variance explained ---
cat("\n=== PCA VARIANCE EXPLAINED ===\n")
cat("First two components explain", round(sum(eig.val[1:2, 2]), 1), "% of variance\n")
cat("Dim1:", round(eig.val[1, 2], 1), "%, Dim2:", round(eig.val[2, 2], 1), "%\n\n")

# --- Figures ---
# Scree plot
p_scree <- fviz_eig(res.pca, addlabels = TRUE, ylim = c(0, 60))
ggsave(file.path(FIGURES_DIR, "PCA_variance.pdf"), plot = p_scree, width = 10, height = 10)

# Variable loadings (contributions)
p_loadings <- fviz_pca_var(res.pca, col.var = "black")
ggsave(file.path(FIGURES_DIR, "PCA_loadings.pdf"), plot = p_loadings, width = 10, height = 10)

# Combined scree + biplot, panel-labelled (a/b) — for Supplementary Fig. 25
library(patchwork)
combined <- (p_scree | p_loadings) + plot_annotation(tag_levels = "a")
ggsave(file.path(FIGURES_DIR, "PCA_scree_and_biplot.pdf"), plot = combined, width = 14, height = 7)

# Cos2 correlation plot
var <- get_pca_var(res.pca)
pdf(file.path(FIGURES_DIR, "PCA_cos2.pdf"), width = 8, height = 6)
corrplot(var$cos2, is.corr = FALSE,
  col = colorRampPalette(c("white", "red"))(200),
  tl.col = "black", tl.srt = 45)
dev.off()

# Contributions to PC1, PC2, PC3
fviz_contrib(res.pca, choice = "var", axes = 1, top = 10)
ggsave(file.path(FIGURES_DIR, "PCA_contrib_PC1.pdf"), width = 10, height = 8)
fviz_contrib(res.pca, choice = "var", axes = 2, top = 10)
ggsave(file.path(FIGURES_DIR, "PCA_contrib_PC2.pdf"), width = 10, height = 8)
fviz_contrib(res.pca, choice = "var", axes = 3, top = 10)
ggsave(file.path(FIGURES_DIR, "PCA_contrib_PC3.pdf"), width = 10, height = 8)

# --- Summary table ---
loadings_pc1 <- res.pca$var$coord[, 1]
loadings_pc2 <- res.pca$var$coord[, 2]
contrib_pc1 <- res.pca$var$contrib[, 1]
contrib_pc2 <- res.pca$var$contrib[, 2]
cos2_pc1 <- res.pca$var$cos2[, 1]
cos2_pc2 <- res.pca$var$cos2[, 2]

pca_summary <- data.frame(
  Variable = names(loadings_pc1),
  Loading_PC1 = round(loadings_pc1, 3),
  Loading_PC2 = round(loadings_pc2, 3),
  Contribution_PC1 = round(contrib_pc1, 2),
  Contribution_PC2 = round(contrib_pc2, 2),
  Cos2_PC1 = round(cos2_pc1, 3),
  Cos2_PC2 = round(cos2_pc2, 3),
  Cos2_Total = round(cos2_pc1 + cos2_pc2, 3)
)
pca_summary <- pca_summary[order(abs(pca_summary$Loading_PC1), decreasing = TRUE), ]

write.csv(pca_summary, file.path(TABLES_DIR, "PCA_variable_importance.csv"), row.names = FALSE)
cat("\nTable saved: PCA_variable_importance.csv\n")

# --- Sample sizes by variable combination (Supplementary Table 9) ---
df_orig <- df
combinations <- list(
  "All clinical variables" = clinical_columns,
  "IQ and SRS only" = c("total_IQ", "SRS_tscore"),
  "IQ subscales and SRS" = c("total_IQ", "verbal_IQ", "performance_IQ", "SRS_tscore"),
  "IQ, SRS, and VABS composite" = c("total_IQ", "SRS_tscore", "vabsabcabc_standard"),
  "IQ, SRS, and all VABS" = c("total_IQ", "SRS_tscore", "vabsdscoresc_dss", "vabsdscoresd_dss", "vabsdscoress_dss", "vabsabcabc_standard"),
  "IQ, SRS, SSP, RBS-R" = c("total_IQ", "SRS_tscore", "ssp_total", "RBS-R_total")
)
sample_sizes <- data.frame(
  Variable_Combination = character(),
  Sample_Size = integer(),
  stringsAsFactors = FALSE
)
POOL_N <- nrow(df_orig)  # denominator for Percent_of_Pool: full participant
                          # pool (deduped, participants only), before the
                          # per-combination completeness filter below.
for (i in seq_along(combinations)) {
  vars <- combinations[[i]][combinations[[i]] %in% colnames(df_orig)]
  if (length(vars) > 0) {
    n_complete <- sum(complete.cases(df_orig[, vars, drop = FALSE]))
    sample_sizes <- rbind(sample_sizes, data.frame(
      Variable_Combination = names(combinations)[i],
      N_Variables = length(vars),
      Sample_Size = n_complete,
      Percent_of_Pool = round(100 * n_complete / POOL_N, 1),
      Pool_N = POOL_N,
      stringsAsFactors = FALSE
    ))
    cat(sprintf("%s: n = %d\n", names(combinations)[i], n_complete))
  }
}
write.csv(sample_sizes, file.path(TABLES_DIR, "sample_sizes_by_variable_combination.csv"), row.names = FALSE)
cat("\nDone. Tables in:", TABLES_DIR, "\n      Figures in:", FIGURES_DIR, "\n")
