#!/usr/bin/env python3
# =============================================================================
# 06 - Method Comparison: hierarchical Ward (paper / op.heatmap) vs K-means
#                        vs Gaussian Mixture Model (GMM)
# =============================================================================
# Reports objective comparative metrics between three candidate clustering
# methods on the same IQ × SRS feature set (N=1023, cohort matched to the
# paper's df_clusters_complete.csv).
#
# Reference partitions (all k=3, labels sorted by descending cluster size):
#   - "ward"    : omniplot.plot.heatmap convention
#                 Z = linkage(squareform(pdist(X, "euclidean")), method="ward")
#                 (i.e. Ward on the distance matrix treated as observations
#                  — matches the paper).
#   - "kmeans"  : sklearn KMeans(n_clusters=3, n_init=25, random_state=42).
#   - "gmm"     : sklearn GaussianMixture(n_components=3,
#                 covariance_type="full", n_init=10, random_state=42).
#                 Matches the model-based block in 03_cluster_validation.R
#                 (mclust VVV, full covariance per component).
#
# Outputs (under ../outputs/):
#   tables/method_comparison_internal_indices.csv
#       Per-method silhouette, Davies-Bouldin, Calinski-Harabasz on the
#       reference partition.
#   tables/method_comparison_concordance.csv
#       Pairwise ARI / NMI / exact-after-alignment between every pair of
#       reference partitions.
#   tables/method_comparison_bootstrap.csv
#       Per-method bootstrap distribution of ARI vs reference (1000 draws).
#   tables/method_comparison_jaccard.csv
#       Per-method, per-reference-cluster Jaccard stability across 1000
#       bootstraps (Hungarian-aligned).
#   tables/method_comparison_modal_prob.csv
#       Per-subject modal-cluster probability under bootstrap, for each
#       method.
#   figures/method_comparison.pdf
#       3-panel figure: (a) bootstrap-ARI violins for the three methods,
#                       (b) per-cluster Jaccard bars,
#                       (c) per-subject modal-probability histograms.
#
# Why this exists: a reviewer asked for objective, quantitative comparison
# beyond visual inspection (ARI between bootstrap solutions, silhouette,
# Davies-Bouldin, Jaccard stability).
# =============================================================================

import os
import sys
import warnings
from itertools import combinations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist, pdist, squareform
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture
from sklearn.metrics import (
    adjusted_rand_score,
    calinski_harabasz_score,
    davies_bouldin_score,
    normalized_mutual_info_score,
    silhouette_score,
)
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings(
    "ignore",
    message=".*looks suspiciously like an uncondensed distance matrix.*",
)

# --- Configuration ---
# curated is the priority regime; pass "curated" to compare methods on the
# curated cluster inputs (features read from the curated assignment table,
# outputs written flat under outputs/, matching 1_run_clustering.py).
CURATED = "curated" in sys.argv[1:] or "--curated" in sys.argv[1:]

_script_dir = os.path.dirname(os.path.abspath(__file__))
DATA_PATH = os.environ.get(
    "LEAP_INOVAND_DATA",
    os.path.join(_script_dir, "..", "..", "..", "imaging2genet", "0_input", "dataframes"),
)
INDIVIDUALS_METRICS = os.path.join(DATA_PATH, "individuals_metrics.tsv")
_SECTION_OUTPUTS = os.path.normpath(os.path.join(_script_dir, "..", "outputs"))
# Curated is flat in the main tree (current pipeline); frozen writes to its
# own self-contained legacy/ location so the unsuffixed filenames below
# (method_comparison.pdf, method_comparison_*.csv) never collide.
_LEGACY_OUTPUTS = os.path.normpath(os.path.join(_script_dir, "..", "..", "legacy", "1_clustering", "outputs"))
OUTPUT_BASE = _SECTION_OUTPUTS if CURATED else _LEGACY_OUTPUTS
FIGURES_DIR = os.path.join(OUTPUT_BASE, "figures")
TABLES_DIR = os.path.join(OUTPUT_BASE, "tables")
# curated features come from the assignment table 1_run_clustering.py wrote.
CURATED_FEATURES = os.path.join(_SECTION_OUTPUTS, "tables", "cluster_assignments_curated.csv")
os.makedirs(FIGURES_DIR, exist_ok=True)
os.makedirs(TABLES_DIR, exist_ok=True)

N_CLUSTERS = 3
RANDOM_STATE = 42
N_ITERS = 1000

EXCLUDED_IDS = {
    "429385763020",
    "C0733-011-137-001",
    "C0733-011-155-001",
}

METHODS = ["ward", "kmeans", "gmm"]
METHOD_LABELS = {"ward": "Ward (paper)", "kmeans": "K-means", "gmm": "GMM"}
METHOD_COLORS = {"ward": "#324095", "kmeans": "#C1C2BC", "gmm": "#5CAEE1"}

rng = np.random.default_rng(RANDOM_STATE)


# --- Helpers ---
def _label_by_size_desc(raw_labels, k):
    sizes = pd.Series(raw_labels).value_counts()
    order = list(sizes.index[:k])
    remap = {old: new for new, old in enumerate(order)}
    return np.array([remap[l] for l in raw_labels])


def fit_ward(X, k=N_CLUSTERS):
    """op.heatmap-equivalent Ward (see 1_run_clustering.py)."""
    D = squareform(pdist(X, metric="euclidean"))
    Z = linkage(D, method="ward")
    raw = fcluster(Z, t=k, criterion="maxclust") - 1
    return _label_by_size_desc(raw, k)


def fit_kmeans(X, k=N_CLUSTERS, seed=RANDOM_STATE):
    km = KMeans(n_clusters=k, n_init=25, random_state=seed).fit(X)
    return _label_by_size_desc(km.labels_, k), km.cluster_centers_


def fit_gmm(X, k=N_CLUSTERS, seed=RANDOM_STATE):
    """Gaussian Mixture (full covariance, matches mclust 'VVV')."""
    gmm = GaussianMixture(
        n_components=k, covariance_type="full", n_init=10, random_state=seed,
    ).fit(X)
    return _label_by_size_desc(gmm.predict(X), k), gmm


def centroids_from_labels(X, labels, k=N_CLUSTERS):
    return np.vstack([X[labels == c].mean(axis=0) for c in range(k)])


def assign_nearest(X, centroids):
    return cdist(X, centroids).argmin(axis=1)


def align_labels(ref, pred, k=N_CLUSTERS):
    cm = np.zeros((k, k), dtype=int)
    for r, p in zip(ref, pred):
        cm[r, p] += 1
    row, col = linear_sum_assignment(-cm)
    mapping = {int(p): int(r) for r, p in zip(row, col)}
    return np.array([mapping[int(p)] for p in pred])


def jaccard_per_cluster(ref, pred_aligned, k=N_CLUSTERS):
    out = np.zeros(k)
    for c in range(k):
        ref_set = ref == c
        pred_set = pred_aligned == c
        inter = np.logical_and(ref_set, pred_set).sum()
        union = np.logical_or(ref_set, pred_set).sum()
        out[c] = inter / union if union > 0 else np.nan
    return out


# --- Data preparation ---
clinical_features = ["IQ", "SRS_tscore"]
if CURATED:
    # Features from the curated assignment table (1_run_clustering.py).
    df_clust = pd.read_csv(CURATED_FEATURES, low_memory=False)
    df_clust = df_clust[clinical_features + ["ID"]].dropna().reset_index(drop=True)
    print(f"[curated] features from {os.path.basename(CURATED_FEATURES)}")
else:
    # Frozen / paper regime (identical to 04 / 05).
    df = pd.read_csv(INDIVIDUALS_METRICS, sep="\t", low_memory=False)
    df = df.drop_duplicates(subset=["ID"]).replace({999: np.nan, 998: np.nan})
    df = df[df["Relation_to_proposant"] == "participant"]
    df["IQ"] = df["total_IQ"].fillna(df["performance_IQ"])
    df_clust = df[clinical_features + ["ID"]].dropna().reset_index(drop=True)
    df_clust = df_clust[~df_clust["ID"].astype(str).isin(EXCLUDED_IDS)].reset_index(drop=True)
X = df_clust[clinical_features].values
X_scaled = StandardScaler().fit_transform(X)
N = len(df_clust)
print(f"N for method comparison{' (curated)' if CURATED else ''}: {N}")

# --- Reference partitions ---
ref_ward = fit_ward(X_scaled, k=N_CLUSTERS)
ref_km, _ = fit_kmeans(X_scaled, k=N_CLUSTERS, seed=RANDOM_STATE)
ref_gmm, _ = fit_gmm(X_scaled, k=N_CLUSTERS, seed=RANDOM_STATE)
REFERENCES = {"ward": ref_ward, "kmeans": ref_km, "gmm": ref_gmm}

print("\nReference cluster sizes:")
for m, lab in REFERENCES.items():
    print(f"  {METHOD_LABELS[m]:14s}:", dict(zip(*np.unique(lab, return_counts=True))))

# --- 1. Internal validation indices on reference partitions ---
internal = []
for m, lab in REFERENCES.items():
    internal.append({
        "method": m,
        "silhouette":         silhouette_score(X_scaled, lab, metric="euclidean"),
        "davies_bouldin":     davies_bouldin_score(X_scaled, lab),
        "calinski_harabasz":  calinski_harabasz_score(X_scaled, lab),
    })
internal_df = pd.DataFrame(internal)
internal_df.to_csv(os.path.join(TABLES_DIR, "method_comparison_internal_indices.csv"), index=False)
print("\n[1/4] Internal validation indices:")
print(internal_df.to_string(index=False))

# --- 2. Pairwise concordance between reference partitions ---
concord = []
for a, b in combinations(METHODS, 2):
    la, lb = REFERENCES[a], REFERENCES[b]
    lb_aligned = align_labels(la, lb, k=N_CLUSTERS)
    concord.append({
        "method_A": a, "method_B": b,
        "ARI":   adjusted_rand_score(la, lb),
        "NMI":   normalized_mutual_info_score(la, lb),
        "exact_match_after_alignment": float((la == lb_aligned).mean()),
    })
concord_df = pd.DataFrame(concord)
concord_df.to_csv(os.path.join(TABLES_DIR, "method_comparison_concordance.csv"), index=False)
print("\n[2/4] Pairwise concordance between reference partitions:")
print(concord_df.round(4).to_string(index=False))

# --- 3 & 4. Bootstrap stability for each method ---
boot_records = []
jacc_records = []
modal = {m: np.zeros((N, N_CLUSTERS), dtype=float) for m in METHODS}

print(f"\n[3/4] Bootstrap stability ({N_ITERS} iterations × {len(METHODS)} methods)...")
for it in range(N_ITERS):
    boot_idx = rng.choice(N, size=N, replace=True)
    Xb = X_scaled[boot_idx]

    # Ward (centroid-based out-of-sample assignment)
    lab_b_w = fit_ward(Xb, k=N_CLUSTERS)
    ward_c_b = centroids_from_labels(Xb, lab_b_w, k=N_CLUSTERS)
    pred_w_all = align_labels(ref_ward, assign_nearest(X_scaled, ward_c_b), k=N_CLUSTERS)

    # K-means (uses fitted centroids)
    _, km_c_b = fit_kmeans(Xb, k=N_CLUSTERS, seed=it)
    pred_k_all = align_labels(ref_km, assign_nearest(X_scaled, km_c_b), k=N_CLUSTERS)

    # GMM (uses .predict on the fitted mixture)
    _, gmm_b = fit_gmm(Xb, k=N_CLUSTERS, seed=it)
    pred_g_all = align_labels(ref_gmm, gmm_b.predict(X_scaled), k=N_CLUSTERS)

    for m, ref, pred in [("ward", ref_ward, pred_w_all),
                         ("kmeans", ref_km, pred_k_all),
                         ("gmm", ref_gmm, pred_g_all)]:
        ari = adjusted_rand_score(ref, pred)
        jac = jaccard_per_cluster(ref, pred, k=N_CLUSTERS)
        modal[m][np.arange(N), pred] += 1
        boot_records.append({"iter": it, "method": m, "ARI": ari})
        for c in range(N_CLUSTERS):
            jacc_records.append({"iter": it, "method": m, "cluster": f"C{c+1}", "jaccard": jac[c]})

boot_df = pd.DataFrame(boot_records)
jacc_df = pd.DataFrame(jacc_records)
boot_df.to_csv(os.path.join(TABLES_DIR, "method_comparison_bootstrap.csv"), index=False)
jacc_df.to_csv(os.path.join(TABLES_DIR, "method_comparison_jaccard.csv"), index=False)

print("\n[4/4] Bootstrap summary (ARI vs reference):")
summary_ari = boot_df.groupby("method")["ARI"].agg(
    mean="mean", sd="std",
    p2_5=lambda s: np.percentile(s, 2.5),
    p97_5=lambda s: np.percentile(s, 97.5),
    share_ge_075=lambda s: float((s >= 0.75).mean()),
).round(4)
print(summary_ari)

print("\nBootstrap summary (Jaccard per cluster):")
summary_j = (jacc_df.groupby(["method", "cluster"])["jaccard"]
             .agg(mean="mean", sd="std",
                  p2_5=lambda s: np.percentile(s, 2.5),
                  p97_5=lambda s: np.percentile(s, 97.5))
             .round(4))
print(summary_j)

# Per-subject modal probability
prob = {m: modal[m] / N_ITERS for m in METHODS}
modal_df = pd.DataFrame({"ID": df_clust["ID"]})
for m in METHODS:
    modal_df[f"ref_{m}"] = [f"C{c+1}" for c in REFERENCES[m]]
    modal_df[f"modal_prob_{m}"] = prob[m].max(axis=1)
modal_df.to_csv(os.path.join(TABLES_DIR, "method_comparison_modal_prob.csv"), index=False)

print("\nPer-subject modal probability (share ≥ 0.90):")
for m in METHODS:
    print(f"  {METHOD_LABELS[m]:14s}: {(modal_df[f'modal_prob_{m}'] >= 0.90).mean():.3f}")

# --- Figure ---
fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))

# (a) Bootstrap ARI distributions
ax = axes[0]
positions = np.arange(len(METHODS))
data = [boot_df[boot_df["method"] == m]["ARI"].values for m in METHODS]
parts = ax.violinplot(data, positions=positions, widths=0.7, showmedians=True, showextrema=False)
for body, m in zip(parts["bodies"], METHODS):
    body.set_facecolor(METHOD_COLORS[m]); body.set_alpha(0.7); body.set_edgecolor("black")
ax.axhline(0.75, ls="--", color="grey", lw=0.8)
ax.set_xticks(positions); ax.set_xticklabels([METHOD_LABELS[m] for m in METHODS])
ax.set_ylabel("ARI vs reference")
ax.set_title("(a) Bootstrap ARI distribution")
ax.set_ylim(0, 1.02)

# (b) Per-cluster Jaccard
ax = axes[1]
clusters = ["C1", "C2", "C3"]
xs = np.arange(len(clusters))
width = 0.27
for i, m in enumerate(METHODS):
    vals = [jacc_df[(jacc_df["method"] == m) & (jacc_df["cluster"] == c)]["jaccard"].values for c in clusters]
    means = np.array([v.mean() for v in vals])
    los   = np.array([np.percentile(v, 2.5)  for v in vals])
    his   = np.array([np.percentile(v, 97.5) for v in vals])
    yerr = np.vstack([means - los, his - means])
    offset = (i - (len(METHODS) - 1) / 2) * width
    ax.bar(xs + offset, means, width=width, yerr=yerr, capsize=3,
           color=METHOD_COLORS[m], edgecolor="black", alpha=0.85,
           label=METHOD_LABELS[m])
ax.set_xticks(xs); ax.set_xticklabels(clusters)
ax.set_ylabel("Jaccard vs reference cluster")
ax.set_title("(b) Per-cluster Jaccard stability")
ax.set_ylim(0, 1.02)
ax.legend(frameon=False)

# (c) Modal-probability histograms
ax = axes[2]
bins = np.linspace(0, 1, 31)
for m in METHODS:
    ax.hist(modal_df[f"modal_prob_{m}"], bins=bins,
            color=METHOD_COLORS[m], alpha=0.45,
            label=METHOD_LABELS[m], edgecolor="black")
ax.axvline(0.90, ls="--", color="grey", lw=0.8)
ax.set_xlabel("Modal cluster assignment probability")
ax.set_ylabel("Individuals")
ax.set_title("(c) Per-individual reassignment stability")
ax.legend(frameon=False)

plt.tight_layout()
out_pdf = os.path.join(FIGURES_DIR, "method_comparison.pdf")
plt.savefig(out_pdf, dpi=300, bbox_inches="tight")
plt.close()
print(f"\nFigure saved: {out_pdf}")
print("Done.")
