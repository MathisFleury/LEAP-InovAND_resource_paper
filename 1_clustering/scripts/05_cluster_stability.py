#!/usr/bin/env python3
# =============================================================================
# 05 - Cluster Stability Analysis
# =============================================================================
# Robustness checks for the k=3 IQ × SRS partition, run for all three
# candidate methods (K-means = primary, Ward = sensitivity, GMM = sensitivity).
# See 04_run_clustering.py for the methods; see 06_method_comparison.py for
# the head-to-head metrics that motivated K-means as primary.
#
# For each method, the four standard stability analyses are run:
#  (1) Subsampling stability — ARI vs reference at several resample fractions.
#  (2) Noise perturbation    — ARI under Gaussian jitter on standardised features.
#  (3) Per-individual reassignment — bootstrap co-assignment matrix.
#  (4) Sample-size curve     — N required for stable cluster structure.
#
# Outputs (1_clustering/outputs/):
#   tables/stability_<analysis>_<method>.csv         (4 × 3 = 12 tables)
#   tables/stability_per_individual_<method>.csv     (3 tables)
#   figures/cluster_stability_<method>.pdf           (3 figures)
# where <method> ∈ {kmeans, ward, gmm} and <analysis> ∈ {subsampling, noise,
# sample_size}.
# =============================================================================

import os
import sys
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist, pdist, squareform
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings(
    "ignore",
    message=".*looks suspiciously like an uncondensed distance matrix.*",
)

# --- Configuration (mirrors 04_run_clustering.py) ---
# Curated (current, priority regime) by default -- reads the ready feature CSV
# (ID, IQ, SRS_tscore) written by 04_run_clustering_curated.py. No frozen/
# deprecated-data fallback here: for the paper-reproduction run, see
# legacy/1_clustering/05_cluster_stability_frozen.py, which derives the same
# ready-CSV shape from the deprecated individuals_metrics.tsv (see root
# CLAUDE.md) and invokes this same engine -- kept out of this file so that
# reference isn't visible in the main tree.
_script_dir = os.path.dirname(os.path.abspath(__file__))
_DEFAULT_FEATURES_CSV = os.path.join(_script_dir, "..", "outputs", "curated", "tables",
                                     "cluster_assignments_curated.csv")
FEATURES_CSV = os.environ.get("STABILITY_FEATURES_CSV", _DEFAULT_FEATURES_CSV)
if not os.path.exists(FEATURES_CSV):
    sys.exit(f"Input not found: {FEATURES_CSV}\nRun 04_run_clustering_curated.py first"
              " (or set STABILITY_FEATURES_CSV to point elsewhere).")
_OUT_SUB = os.environ.get("STABILITY_OUTPUT_SUBDIR", "curated")
OUTPUT_BASE = os.path.normpath(os.path.join(_script_dir, "..", "outputs", _OUT_SUB))
FIGURES_DIR = os.path.join(OUTPUT_BASE, "figures")
TABLES_DIR = os.path.join(OUTPUT_BASE, "tables")
os.makedirs(FIGURES_DIR, exist_ok=True)
os.makedirs(TABLES_DIR, exist_ok=True)

N_CLUSTERS = 3
RANDOM_STATE = 42
N_ITERS = 1000
SUBSAMPLE_FRACS = [0.50, 0.60, 0.70, 0.80, 0.90]
NOISE_SDS = [0.05, 0.10, 0.20, 0.30, 0.50]
SAMPLE_SIZE_GRID = [50, 100, 150, 200, 300, 400, 500, 600, 800, 1000]
STABLE_ARI = 0.75

METHODS = ["kmeans", "ward", "gmm"]
METHOD_LABELS = {"kmeans": "K-means (primary)", "ward": "Ward (sensitivity)", "gmm": "GMM (sensitivity)"}
METHOD_COLORS = {"kmeans": "#324095", "ward": "#5CAEE1", "gmm": "#7a7a7a"}


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
def _label_by_size_desc(raw_labels, k):
    sizes = pd.Series(raw_labels).value_counts()
    order = list(sizes.index[:k])
    remap = {old: new for new, old in enumerate(order)}
    return np.array([remap[l] for l in raw_labels])


def fit_kmeans(X, k=N_CLUSTERS, seed=RANDOM_STATE):
    km = KMeans(n_clusters=k, n_init=25, random_state=seed).fit(X)
    return _label_by_size_desc(km.labels_, k)


def fit_kmeans_with_predict_full(X_fit, X_full, k=N_CLUSTERS, seed=RANDOM_STATE):
    """Fit on X_fit, return labels for X_full (via nearest centroid)."""
    km = KMeans(n_clusters=k, n_init=25, random_state=seed).fit(X_fit)
    return km.predict(X_full)


def fit_ward(X, k=N_CLUSTERS, seed=RANDOM_STATE):
    """op.heatmap convention: linkage on the squareform distance matrix."""
    D = squareform(pdist(X, metric="euclidean"))
    Z = linkage(D, method="ward")
    raw = fcluster(Z, t=k, criterion="maxclust") - 1
    return _label_by_size_desc(raw, k)


def fit_ward_with_predict_full(X_fit, X_full, k=N_CLUSTERS, seed=RANDOM_STATE):
    """Ward has no native predictor — fit, compute per-cluster centroids on
    the fitting sample, then assign full-sample points to nearest centroid."""
    labels = fit_ward(X_fit, k=k)
    centroids = np.vstack([X_fit[labels == c].mean(axis=0) for c in range(k)])
    return cdist(X_full, centroids).argmin(axis=1)


def fit_gmm(X, k=N_CLUSTERS, seed=RANDOM_STATE):
    gmm = GaussianMixture(
        n_components=k, covariance_type="full", n_init=10, random_state=seed,
    ).fit(X)
    return _label_by_size_desc(gmm.predict(X), k)


def fit_gmm_with_predict_full(X_fit, X_full, k=N_CLUSTERS, seed=RANDOM_STATE):
    gmm = GaussianMixture(
        n_components=k, covariance_type="full", n_init=10, random_state=seed,
    ).fit(X_fit)
    return gmm.predict(X_full)


def align_labels(ref, pred, k=N_CLUSTERS):
    cm = np.zeros((k, k), dtype=int)
    for r, p in zip(ref, pred):
        cm[r, p] += 1
    row, col = linear_sum_assignment(-cm)
    mapping = {int(p): int(r) for r, p in zip(row, col)}
    return np.array([mapping[int(p)] for p in pred])


FIT_FNS = {"kmeans": fit_kmeans, "ward": fit_ward, "gmm": fit_gmm}
FIT_FULL_FNS = {
    "kmeans": fit_kmeans_with_predict_full,
    "ward":   fit_ward_with_predict_full,
    "gmm":    fit_gmm_with_predict_full,
}


# -----------------------------------------------------------------------------
# Data preparation
# -----------------------------------------------------------------------------
clinical_features = ["IQ", "SRS_tscore"]
print(f"Loading features from ready CSV: {FEATURES_CSV}")
df = pd.read_csv(FEATURES_CSV, low_memory=False)
for c in clinical_features:
    df[c] = pd.to_numeric(df.get(c), errors="coerce")
df_clust = df[clinical_features + ["ID"]].dropna().reset_index(drop=True)
X = df_clust[clinical_features].values
X_scaled = StandardScaler().fit_transform(X)
N = len(df_clust)
print(f"N for stability analysis: {N}")


# -----------------------------------------------------------------------------
# Per-method runner
# -----------------------------------------------------------------------------
def run_stability_for_method(method):
    """Run all four stability analyses for one method, write tables + figure."""
    rng = np.random.default_rng(RANDOM_STATE)
    fit_fn = FIT_FNS[method]
    fit_full_fn = FIT_FULL_FNS[method]
    suffix = f"_{method}"
    label = METHOD_LABELS[method]
    color = METHOD_COLORS[method]

    print(f"\n=== {label} ===")
    ref_labels = fit_fn(X_scaled, k=N_CLUSTERS, seed=RANDOM_STATE)

    # (1) Subsampling
    print("  [1/4] Subsampling stability...")
    sub_records = []
    for frac in SUBSAMPLE_FRACS:
        aris = np.empty(N_ITERS)
        for it in range(N_ITERS):
            idx = rng.choice(N, size=int(frac * N), replace=False)
            sub_labels = fit_fn(X_scaled[idx], k=N_CLUSTERS, seed=it)
            aris[it] = adjusted_rand_score(ref_labels[idx], sub_labels)
        sub_records.append({
            "fraction": frac, "n_subset": int(frac * N),
            "ARI_mean": aris.mean(), "ARI_sd": aris.std(),
            "ARI_p2.5": np.percentile(aris, 2.5),
            "ARI_p97.5": np.percentile(aris, 97.5),
            "share_ARI_ge_0.75": float((aris >= STABLE_ARI).mean()),
        })
    sub_df = pd.DataFrame(sub_records)
    sub_df.to_csv(os.path.join(TABLES_DIR, f"stability_subsampling{suffix}.csv"), index=False)

    # (2) Noise
    print("  [2/4] Noise perturbation...")
    noise_records = []
    for sd in NOISE_SDS:
        aris = np.empty(N_ITERS)
        for it in range(N_ITERS):
            Xn = X_scaled + rng.normal(0, sd, size=X_scaled.shape)
            noise_labels = fit_fn(Xn, k=N_CLUSTERS, seed=it)
            aris[it] = adjusted_rand_score(ref_labels, noise_labels)
        noise_records.append({
            "noise_sd": sd,
            "ARI_mean": aris.mean(), "ARI_sd": aris.std(),
            "ARI_p2.5": np.percentile(aris, 2.5),
            "ARI_p97.5": np.percentile(aris, 97.5),
        })
    noise_df = pd.DataFrame(noise_records)
    noise_df.to_csv(os.path.join(TABLES_DIR, f"stability_noise{suffix}.csv"), index=False)

    # (3) Per-individual reassignment
    print("  [3/4] Per-individual reassignment...")
    co_assign = np.zeros((N, N_CLUSTERS), dtype=float)
    for it in range(N_ITERS):
        boot_idx = rng.choice(N, size=N, replace=True)
        pred_raw = fit_full_fn(X_scaled[boot_idx], X_scaled, k=N_CLUSTERS, seed=it)
        pred_aligned = align_labels(ref_labels, pred_raw)
        co_assign[np.arange(N), pred_aligned] += 1
    prob = co_assign / N_ITERS
    modal_prob = prob.max(axis=1)
    df_individual = pd.DataFrame({
        "ID": df_clust["ID"],
        "ref_cluster": [f"C{c + 1}" for c in ref_labels],
        "modal_assignment_prob": modal_prob,
        "P_C1": prob[:, 0], "P_C2": prob[:, 1], "P_C3": prob[:, 2],
    })
    df_individual.to_csv(os.path.join(TABLES_DIR, f"stability_per_individual{suffix}.csv"), index=False)
    print(f"    share with modal prob >= 0.90: {(modal_prob >= 0.90).mean():.3f}")
    print(f"    share with modal prob >= 0.80: {(modal_prob >= 0.80).mean():.3f}")
    print(f"    median modal prob: {np.median(modal_prob):.3f}")

    # (4) Sample-size curve
    print("  [4/4] Sample-size curve...")
    size_records = []
    for n in SAMPLE_SIZE_GRID:
        if n > N:
            continue
        aris = np.empty(N_ITERS)
        for it in range(N_ITERS):
            idx = rng.choice(N, size=n, replace=False)
            size_labels = fit_fn(X_scaled[idx], k=N_CLUSTERS, seed=it)
            aris[it] = adjusted_rand_score(ref_labels[idx], size_labels)
        size_records.append({
            "n": n,
            "ARI_mean": aris.mean(), "ARI_sd": aris.std(),
            "ARI_p2.5": np.percentile(aris, 2.5),
            "ARI_p97.5": np.percentile(aris, 97.5),
            "share_ARI_ge_0.75": float((aris >= STABLE_ARI).mean()),
        })
    size_df = pd.DataFrame(size_records)
    size_df.to_csv(os.path.join(TABLES_DIR, f"stability_sample_size{suffix}.csv"), index=False)

    n_stable = size_df.loc[size_df["ARI_mean"] >= STABLE_ARI, "n"]
    n_star = int(n_stable.min()) if len(n_stable) else None
    print(f"    smallest N with mean ARI >= {STABLE_ARI}: {n_star}")

    # --- Figure ---
    fig, axes = plt.subplots(1, 4, figsize=(20, 4.5))
    axes[0].hist(modal_prob, bins=30, color=color, alpha=0.85)
    axes[0].axvline(0.90, ls="--", color="grey", lw=0.8)
    axes[0].set_xlabel("Modal cluster assignment probability")
    axes[0].set_ylabel("Individuals")
    axes[0].set_title(f"(a) Per-individual stability — {label}")

    for ax, d, x, xlabel, title in [
        (axes[1], sub_df,   "n_subset", "Subsample size",            "(b) Subsampling"),
        (axes[2], noise_df, "noise_sd", "Noise SD (std-space units)", "(c) Noise"),
        (axes[3], size_df,  "n",        "N (training subset)",        "(d) Sample-size curve"),
    ]:
        ax.errorbar(
            d[x], d["ARI_mean"],
            yerr=[d["ARI_mean"] - d["ARI_p2.5"], d["ARI_p97.5"] - d["ARI_mean"]],
            marker="o", capsize=3, color=color,
        )
        ax.axhline(STABLE_ARI, ls="--", color="grey", lw=0.8)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("ARI vs reference")
        ax.set_title(title)
        ax.set_ylim(0, 1.02)

    plt.tight_layout()
    plt.savefig(os.path.join(FIGURES_DIR, f"cluster_stability{suffix}.pdf"),
                dpi=300, bbox_inches="tight")
    plt.close()

    return {
        "method": method,
        "modal_prob_ge_0.90": float((modal_prob >= 0.90).mean()),
        "modal_prob_ge_0.80": float((modal_prob >= 0.80).mean()),
        "modal_prob_median":  float(np.median(modal_prob)),
        "smallest_n_stable":  n_star,
    }


# -----------------------------------------------------------------------------
# Run all three methods
# -----------------------------------------------------------------------------
summary_rows = [run_stability_for_method(m) for m in METHODS]
summary = pd.DataFrame(summary_rows)
summary.to_csv(os.path.join(TABLES_DIR, "stability_summary_all_methods.csv"), index=False)

print("\n========================================")
print("Summary across methods")
print("========================================")
print(summary.to_string(index=False))
print("\nDone. Tables in:", TABLES_DIR)
print("Figures in:", FIGURES_DIR)
