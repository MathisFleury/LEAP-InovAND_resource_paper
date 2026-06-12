#!/usr/bin/env python3
# =============================================================================
# 04 - Run Clustering (k=3)
# =============================================================================
# Primary clustering method: KMeans on standardised IQ × SRS (k=3, n_init=25).
# Sensitivity method: hierarchical Ward as used in the original paper
# (op.heatmap convention — see _ward_assignments() below).
#
# The k-means method was adopted as primary after a reviewer-requested
# comparative analysis (see 06_method_comparison.py): k-means is more
# stable under bootstrap (per-subject modal probability 95 % vs 80 %;
# bootstrap ARI ⟨0.96⟩ vs ⟨0.72⟩), produces better internal indices
# (silhouette 0.45 vs 0.41, Davies-Bouldin 0.83 vs 0.92, Calinski-
# Harabasz 1097 vs 929), and is robust to small input perturbations.
# Ward is retained as a sensitivity analysis to enable direct
# reproduction of the original paper's figures.
#
# Outputs (1_clustering/outputs/):
#   tables/cluster_assignments.csv                 — k-means (primary)
#   tables/individuals_metrics_with_clusters.csv   — full table + k-means
#   tables/cluster_assignments_ward.csv            — Ward (sensitivity)
#   tables/individuals_metrics_with_clusters_ward.csv
#   figures/clustering_scatter.pdf                 — k-means scatter
#   figures/clustering_scatter_ward.pdf            — Ward scatter
#
# Additionally writes the k-means labels into the sibling repo as
#   eeg_mri-pipeline/results/dataset_paper/dataframes/df_clusters_complete_kmeans.csv
# so that downstream cluster-* analyses across the project pick up the
# new primary method via a one-line filename swap.
# =============================================================================

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist, squareform
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler


# -----------------------------------------------------------------------------
# Clustering helpers
# -----------------------------------------------------------------------------

def _label_by_size_desc(raw_labels, k):
    """Relabel cluster IDs by descending cluster size: largest → 0 (= C1).
    Stable label naming across runs and methods."""
    sizes = pd.Series(raw_labels).value_counts()
    order = list(sizes.index[:k])
    remap = {old: new for new, old in enumerate(order)}
    return np.array([remap[l] for l in raw_labels])


def _kmeans_assignments(X_scaled, k, random_state=42):
    km = KMeans(n_clusters=k, n_init=25, random_state=random_state).fit(X_scaled)
    return _label_by_size_desc(km.labels_, k)


def _ward_assignments(X_scaled, k):
    """Hierarchical Ward (op.heatmap convention) — retained as sensitivity.
    omniplot.plot.heatmap runs linkage on squareform(pdist(X)), i.e. treats
    the n×n distance matrix as observations. This is NOT plain Ward on the
    features; it reproduces the original paper's df_clusters_complete.csv
    exactly when combined with the EXCLUDED_IDS filter below."""
    D = squareform(pdist(X_scaled, metric="euclidean"))
    Z = linkage(D, method="ward")
    raw = fcluster(Z, t=k, criterion="maxclust") - 1
    return _label_by_size_desc(raw, k)


def _gmm_assignments(X_scaled, k, random_state=42):
    """Gaussian Mixture (full covariance, matches mclust 'VVV' in
    03_cluster_validation.R) — retained as sensitivity."""
    gmm = GaussianMixture(
        n_components=k, covariance_type="full", n_init=10, random_state=random_state,
    ).fit(X_scaled)
    return _label_by_size_desc(gmm.predict(X_scaled), k)


# IDs present in our input but not labelled in the paper's
# df_clusters_complete.csv (excluded at QC time when the paper's dataframe
# was frozen). Removing them keeps the cohort identical to the paper.
EXCLUDED_IDS = {
    "429385763020",
    "C0733-011-137-001",
    "C0733-011-155-001",
}


# --- Configuration ---
_script_dir = os.path.dirname(os.path.abspath(__file__))
DATA_PATH = os.environ.get(
    "LEAP_INOVAND_DATA",
    os.path.join(_script_dir, "..", "..", "..", "imaging2genet", "0_input", "dataframes")
)
INDIVIDUALS_METRICS = os.path.join(DATA_PATH, "individuals_metrics.tsv")
OUTPUT_BASE = os.path.normpath(os.path.join(_script_dir, "..", "outputs"))
FIGURES_DIR = os.path.join(OUTPUT_BASE, "figures")
TABLES_DIR = os.path.join(OUTPUT_BASE, "tables")
os.makedirs(FIGURES_DIR, exist_ok=True)
os.makedirs(TABLES_DIR, exist_ok=True)

# Cluster-complete dataframe propagated to the sibling repo so Group A
# downstream consumers can switch with a single filename change.
SIBLING_CLUSTERS_DIR = os.path.normpath(os.path.join(
    _script_dir, "..", "..", "..", "eeg_mri-pipeline",
    "results", "dataset_paper", "dataframes",
))
SIBLING_REF_FILE     = os.path.join(SIBLING_CLUSTERS_DIR, "df_clusters_complete.csv")
SIBLING_KMEANS_FILE  = os.path.join(SIBLING_CLUSTERS_DIR, "df_clusters_complete_kmeans.csv")

N_CLUSTERS   = 3
RANDOM_STATE = 42

# --- Data preparation ---
print(f"Loading data from: {INDIVIDUALS_METRICS}")
df = pd.read_csv(INDIVIDUALS_METRICS, sep="\t", low_memory=False)
df = df.drop_duplicates(subset=["ID"])
df = df.replace(999, np.nan)
df = df.replace(998, np.nan)
df = df[df["Relation_to_proposant"] == "participant"]

df["IQ"] = df["total_IQ"].fillna(df["performance_IQ"])

clinical_features = ["IQ", "SRS_tscore"]
df_clust = df[clinical_features + ["ID", "Population1", "PopulationS1"]].dropna()
df_clust = df_clust[~df_clust["ID"].astype(str).isin(EXCLUDED_IDS)]
X = df_clust[clinical_features].values
X_scaled = StandardScaler().fit_transform(X)

print(f"Sample size for clustering: {len(df_clust)}")

# --- Clustering (primary = k-means; sensitivity = Ward + GMM) ---
labels_kmeans = _kmeans_assignments(X_scaled, N_CLUSTERS, random_state=RANDOM_STATE)
labels_ward   = _ward_assignments(X_scaled, N_CLUSTERS)
labels_gmm    = _gmm_assignments(X_scaled, N_CLUSTERS, random_state=RANDOM_STATE)

df_clust = df_clust.copy()
df_clust["Cluster"]      = [f"C{c + 1}" for c in labels_kmeans]
df_clust["Cluster_ward"] = [f"C{c + 1}" for c in labels_ward]
df_clust["Cluster_gmm"]  = [f"C{c + 1}" for c in labels_gmm]


def _emit(method_name, label_col, scatter_pdf, assign_csv, full_csv):
    """Write scatter + per-method tables for one labelling method."""
    # Same palette as 2_genetic_analysis/scripts/_config.py::PALETTE_CLUSTERS
    palette = {"C1": "#7A8B47", "C2": "#ff9fa0", "C3": "#e7ba52"}
    fig, ax = plt.subplots(figsize=(10, 8))
    sns.scatterplot(
        data=df_clust, x="SRS_tscore", y="IQ",
        hue=label_col, palette=palette, alpha=0.7, s=50, ax=ax,
    )
    ax.set_xlabel("SRS-2 t-score")
    ax.set_ylabel("Full-scale IQ")
    ax.set_title(f"Clinical Clustering: IQ and SRS (k=3, {method_name})")
    ax.legend(title="Cluster")
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, scatter_pdf), dpi=300, bbox_inches="tight")
    plt.close()

    out_cols = ["ID", "IQ", "SRS_tscore", label_col, "Population1", "PopulationS1"]
    df_clust[out_cols].rename(columns={label_col: "Cluster"}).to_csv(
        os.path.join(TABLES_DIR, assign_csv), index=False,
    )
    full = df.merge(df_clust[["ID", label_col]].rename(columns={label_col: "Cluster"}),
                    on="ID", how="left")
    full.to_csv(os.path.join(TABLES_DIR, full_csv), index=False)

    print(f"  [{method_name}] sizes:", df_clust[label_col].value_counts().sort_index().to_dict())
    print(f"  [{method_name}] saved: {scatter_pdf}, {assign_csv}, {full_csv}")


print("\n--- Writing per-method outputs ---")
_emit("k-means", "Cluster",       "clustering_scatter.pdf",
      "cluster_assignments.csv", "individuals_metrics_with_clusters.csv")
_emit("Ward (sensitivity)", "Cluster_ward", "clustering_scatter_ward.pdf",
      "cluster_assignments_ward.csv", "individuals_metrics_with_clusters_ward.csv")
_emit("GMM (sensitivity)", "Cluster_gmm", "clustering_scatter_gmm.pdf",
      "cluster_assignments_gmm.csv", "individuals_metrics_with_clusters_gmm.csv")


# --- Propagate k-means labels to the sibling repo ---
# All Group A scripts read df_clusters_complete.csv from the sibling repo.
# We write a parallel df_clusters_complete_kmeans.csv with the same row
# structure (same IDs, same metadata columns) but with the Cluster column
# replaced by k-means assignments. Downstream scripts that want k-means
# just swap the filename suffix.
print("\n--- Writing sibling-repo df_clusters_complete_kmeans.csv ---")
if os.path.exists(SIBLING_REF_FILE):
    ref_df = pd.read_csv(SIBLING_REF_FILE, low_memory=False)
    ref_df["ID"] = ref_df["ID"].astype(str)
    ref_df = ref_df.drop(columns=["Cluster"], errors="ignore")
    new_labels = df_clust[["ID", "Cluster"]].copy()
    new_labels["ID"] = new_labels["ID"].astype(str)
    out = ref_df.merge(new_labels, on="ID", how="left")
    out.to_csv(SIBLING_KMEANS_FILE, index=False)
    print(f"  wrote: {SIBLING_KMEANS_FILE}")
    print(f"  rows: {len(out):,}  labelled: {out['Cluster'].notna().sum():,}")
else:
    print(f"  WARNING: {SIBLING_REF_FILE} not found — skipping sibling-repo propagation.")

print("\nDone. Outputs in:", OUTPUT_BASE)
