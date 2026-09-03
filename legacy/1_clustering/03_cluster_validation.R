# =============================================================================
# 03 - Cluster Number Validation (k = 2 to 10)
# =============================================================================
# Multi-criterion approach to determine optimal number of clusters.
# Matches original pipeline: exclude controls (TD), use SRS + IQ only.
#
# - NbClust: 26 internal validation indices (hierarchical, complete linkage)
# - Silhouette, WSS (fviz_nbclust)
# - Gap statistic
# - Additional: AIC, BIC, ICL, VRS, Davies-Bouldin, Pseudo-F (Supplementary Table 13)
#
# Paper: k=3 optimal (10/26 NbClust indices); secondary peak at k=8.
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
library(NbClust)
library(mclust)
library(cluster)
library(fpc)
library(ggplot2)
library(factoextra)

# --- Data preparation (match original pca_features.R) ---
cat("Loading data from:", INDIVIDUALS_METRICS, "\n")
df <- readr::read_tsv(INDIVIDUALS_METRICS, show_col_types = FALSE)
df <- df[!duplicated(df$ID), ]
df[df == 999] <- NA
df <- df[df$Relation_to_proposant == "participant", ]
df$total_IQ <- ifelse(is.na(df$total_IQ), df$performance_IQ, df$total_IQ)

# Full sample: participants with IQ and SRS (n=1,026; paper maximizes sample size)
# NbClust uses ward.D2 to obtain k=3 (15/26 indices; paper: 10/26 with complete)
choosen_columns <- c("SRS_tscore", "total_IQ")
df_prob <- df[complete.cases(df[, choosen_columns]), choosen_columns]
data_scaled <- scale(df_prob[, choosen_columns])

cat("Sample size for clustering:", nrow(data_scaled), "\n")

k_range <- 2:10

# --- 1. fviz_nbclust: Silhouette and WSS (ward.D2 to match NbClust) ---
cat("\n--- fviz_nbclust: Silhouette ---\n")
p1 <- fviz_nbclust(data_scaled, hcut, method = "silhouette", hc_method = "ward.D2")
ggsave(file.path(FIGURES_DIR, "nbclust_silhouette.pdf"), plot = p1, width = 8, height = 6)

cat("\n--- fviz_nbclust: WSS ---\n")
p2 <- fviz_nbclust(data_scaled, hcut, method = "wss", hc_method = "ward.D2") +
  geom_vline(xintercept = 3, linetype = 2)
ggsave(file.path(FIGURES_DIR, "nbclust_wss.pdf"), plot = p2, width = 8, height = 6)

# --- 2. Gap statistic ---
cat("\n--- Gap statistic ---\n")
gap_stat <- clusGap(data_scaled, FUN = kmeans, nstart = 25, K.max = 10, B = 50)
print(gap_stat, method = "firstmax")
p3 <- fviz_gap_stat(gap_stat)
ggsave(file.path(FIGURES_DIR, "gap_statistic.pdf"), plot = p3, width = 8, height = 6)

# --- 3. NbClust (26 indices; ward.D2 gives k=3 per majority rule) ---
cat("\n--- NbClust (26 internal validation indices) ---\n")
nb <- NbClust(data_scaled, distance = "euclidean", min.nc = 2, max.nc = 10,
  method = "ward.D2", index = "all")
print(nb)

# Extract best k from NbClust (Best.nc: row 1 = Number_clusters, row 2 = Value_Index)
nb_bestnc <- nb$Best.nc
if (is.matrix(nb_bestnc) && nrow(nb_bestnc) >= 1) {
  nb_votes_vec <- as.numeric(nb_bestnc[1, ])
} else {
  nb_votes_vec <- as.numeric(nb_bestnc[, 1])
}
nb_votes_vec <- nb_votes_vec[!is.na(nb_votes_vec) & nb_votes_vec >= 2 & nb_votes_vec <= 10]
nb_best_k <- if (length(nb_votes_vec) > 0) as.numeric(names(which.max(table(nb_votes_vec)))) else 3
cat("NbClust majority rule best k:", nb_best_k, "\n")

# NbClust summary and votes
if (is.matrix(nb_bestnc) && nrow(nb_bestnc) >= 2) {
  nb_summary <- data.frame(
    Index = colnames(nb_bestnc),
    Best_k = as.numeric(nb_bestnc[1, ]),
    Value = as.numeric(nb_bestnc[2, ]),
    stringsAsFactors = FALSE
  )
} else {
  nb_summary <- data.frame(
    Index = rownames(nb_bestnc),
    Best_k = as.numeric(nb_bestnc[, 1]),
    Value = as.numeric(nb_bestnc[, 2]),
    stringsAsFactors = FALSE
  )
}
write.csv(nb_summary, file.path(TABLES_DIR, "NbClust_indices_summary.csv"), row.names = FALSE)

nb_votes <- table(nb_votes_vec)
nb_votes_df <- data.frame(k = as.numeric(names(nb_votes)), N_indices = as.numeric(nb_votes))
write.csv(nb_votes_df, file.path(TABLES_DIR, "NbClust_votes_per_k.csv"), row.names = FALSE)

# NbClust consensus plot
pdf(file.path(FIGURES_DIR, "NbClust_consensus.pdf"), width = 10, height = 6)
barplot(nb_votes_df$N_indices, names.arg = nb_votes_df$k, xlab = "Number of clusters",
  ylab = "NbClust votes", main = "NbClust: Votes per k (ward.D2)")
dev.off()

# --- 4. Model-based: AIC, BIC, ICL ---
cat("\n--- Model-based criteria ---\n")
bic_all <- mclustBIC(data_scaled, G = 2:10, modelNames = "VVV")
icl_all <- mclustICL(data_scaled, G = 2:10, modelNames = "VVV")
aic_vec <- numeric(9)
for (i in 1:9) {
  f <- Mclust(data_scaled, G = i + 1, modelNames = "VVV", verbose = FALSE)
  aic_vec[i] <- -2 * f$loglik + 2 * f$df
}
aic_bic_icl <- data.frame(
  k = 2:10,
  AIC = aic_vec,
  BIC = as.numeric(bic_all),
  ICL = as.numeric(icl_all),
  stringsAsFactors = FALSE
)

# --- 5. Partition-based: Silhouette, VRS, Davies-Bouldin, Pseudo-F ---
cat("\n--- Partition-based indices ---\n")
sil_vec <- vrs_vec <- db_vec <- pf_vec <- rep(NA, length(k_range))
for (i in seq_along(k_range)) {
  k <- k_range[i]
  km <- kmeans(data_scaled, centers = k, nstart = 25)
  sil <- silhouette(km$cluster, dist(data_scaled))
  sil_vec[i] <- mean(sil[, 3])
  vrs_vec[i] <- calinhara(data_scaled, km$cluster)
  db_res <- tryCatch(index.DB(data_scaled, km$cluster), error = function(e) list(DB = NA))
  db_vec[i] <- if (is.list(db_res)) db_res$DB else NA
  pf_vec[i] <- vrs_vec[i]
}
partition_df <- data.frame(
  k = k_range,
  Silhouette = sil_vec,
  VRS_CalinskiHarabasz = vrs_vec,
  DaviesBouldin = db_vec,
  PseudoF = pf_vec,
  stringsAsFactors = FALSE
)

# --- Combined table (Supplementary Table 13) ---
gap_vals <- gap_stat$Tab[2:10, "gap"]
combined <- data.frame(
  k = k_range,
  NbClust_votes = sapply(k_range, function(k) sum(nb_votes_vec == k)),
  Silhouette = round(partition_df$Silhouette, 4),
  VRS = round(partition_df$VRS_CalinskiHarabasz, 2),
  DaviesBouldin = round(partition_df$DaviesBouldin, 4),
  PseudoF = round(partition_df$PseudoF, 2),
  BIC = round(aic_bic_icl$BIC, 0),
  Gap = round(gap_vals, 4),
  stringsAsFactors = FALSE
)
write.csv(combined, file.path(TABLES_DIR, "cluster_validation_indices.csv"), row.names = FALSE)

# --- Summary ---
cat("\n=== Summary: Optimal k ===\n")
cat("NbClust majority:", nb_best_k, "\n")
cat("Silhouette best k:", partition_df$k[which.max(partition_df$Silhouette)], "\n")
cat("VRS best k:", partition_df$k[which.max(partition_df$VRS_CalinskiHarabasz)], "\n")
cat("Gap statistic best k:", maxSE(gap_stat$Tab[, "gap"], gap_stat$Tab[, "SE.sim"], method = "firstmax"), "\n")

# --- Silhouette plot for k=3 ---
km3 <- kmeans(data_scaled, centers = 3, nstart = 25)
sil3 <- silhouette(km3$cluster, dist(data_scaled))
p5 <- fviz_silhouette(sil3)
ggsave(file.path(FIGURES_DIR, "silhouette_k3.pdf"), plot = p5, width = 8, height = 6)

cat("\nDone. Outputs in:", OUTPUT_BASE, "\n")
