#!/usr/bin/env python3
# =============================================================================
# 07 - reval k-selection (stability-based RELATIVE validation)
# =============================================================================
# Complements 05_cluster_stability.py. Where 05 asks "given k=3, how robust is
# the partition?", this asks "which k generalises out-of-sample?" using the
# stability-based relative validation of Lange et al. (2004), as implemented in
# the reval package (IIT-LAND/reval_clustering).
#
# Method: clustering is recast as supervised classification. Data is split
# train/test; for each candidate k the training set is clustered, a classifier
# (KNN) learns those labels, predicts test labels, and the misclassification vs
# an independent test-set clustering gives a normalised stability. Repeated
# 10x2 CV; the k minimising validation stability is "best".
#
# Outputs (1_clustering/outputs/):
#   tables/reval_kselection.csv          — per-k train/val stability
#   figures/reval_kselection.pdf         — stability vs k curve
#
# Data prep mirrors 04_run_clustering.py exactly (same cohort, same features,
# same EXCLUDED_IDS) so the k it endorses is comparable to the k=3 used there.
# =============================================================================

import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.model_selection import train_test_split
from sklearn.neighbors import KNeighborsClassifier
from sklearn.preprocessing import StandardScaler

from reval.best_nclust_cv import FindBestClustCV

# --- Configuration (mirrors 04_run_clustering.py) ---
_script_dir = os.path.dirname(os.path.abspath(__file__))
DATA_PATH = os.environ.get(
    "LEAP_INOVAND_DATA",
    os.path.join(_script_dir, "..", "..", "..", "imaging2genet", "0_input", "dataframes"),
)
INDIVIDUALS_METRICS = os.path.join(DATA_PATH, "individuals_metrics.tsv")
# Self-contained under legacy/1_clustering/ (moved here; script sits directly
# in this dir now, not in a nested scripts/ subfolder, so no ".." needed).
OUTPUT_BASE = os.path.normpath(os.path.join(_script_dir, "outputs"))
FIGURES_DIR = os.path.join(OUTPUT_BASE, "figures")
TABLES_DIR = os.path.join(OUTPUT_BASE, "tables")
os.makedirs(FIGURES_DIR, exist_ok=True)
os.makedirs(TABLES_DIR, exist_ok=True)

RANDOM_STATE = 42
NCLUST_RANGE = list(range(2, 7))   # reval treats this as an explicit list: {2,3,4,5,6}
N_RAND = 10             # random-labelling baseline iterations
N_FOLD = 2              # 2-fold inner CV
ITER_CV = 10            # repeated CV (=> 10x2)

EXCLUDED_IDS = {
    "429385763020",
    "C0733-011-137-001",
    "C0733-011-155-001",
}

# --- Data preparation (identical to 04_run_clustering.py) ---
df = pd.read_csv(INDIVIDUALS_METRICS, sep="\t", low_memory=False)
df = df.drop_duplicates(subset=["ID"]).replace({999: np.nan, 998: np.nan})
df = df[df["Relation_to_proposant"] == "participant"]
df["IQ"] = df["total_IQ"].fillna(df["performance_IQ"])
clinical_features = ["IQ", "SRS_tscore"]
df_clust = df[clinical_features + ["ID"]].dropna()
df_clust = df_clust[~df_clust["ID"].astype(str).isin(EXCLUDED_IDS)].reset_index(drop=True)
X_scaled = StandardScaler().fit_transform(df_clust[clinical_features].values)
print(f"N for reval k-selection: {len(X_scaled)}")

# --- reval: train/test split, then stability-based relative validation ---
X_tr, X_ts = train_test_split(
    X_scaled, test_size=0.40, random_state=RANDOM_STATE,
)

classifier = KNeighborsClassifier(n_neighbors=5)
clustering = KMeans(n_init=25, random_state=RANDOM_STATE)  # reval sets n_clusters per k

findbest = FindBestClustCV(
    s=classifier, c=clustering,
    nrand=N_RAND, nfold=N_FOLD, nclust_range=NCLUST_RANGE,
)
metrics, nbest = findbest.best_nclust(X_tr, iter_cv=ITER_CV)
out = findbest.evaluate(X_tr, X_ts, nclust=nbest)

# metrics['train'|'val'][k] = [mean_normalised_stability, (sem, error)]
ks = sorted(metrics["val"].keys())
rows = []
for k in ks:
    val = metrics["val"][k]
    tr = metrics["train"][k]
    rows.append({
        "k": k,
        "val_stability": float(val[0]),
        # reval stores (mean_stability, (mean_stability, CI_half_width)); take the CI.
        "val_stability_sem": float(val[1][1]),
        "train_stability": float(tr[0]),
        "is_best": k == nbest,
    })
res = pd.DataFrame(rows)
res.to_csv(os.path.join(TABLES_DIR, "reval_kselection.csv"), index=False)

print("\nNormalised validation stability (lower = more stable / generalisable):")
print(res.to_string(index=False))
print(f"\nreval best k = {nbest}")
print(f"held-out test accuracy at k={nbest}: {out.test_acc:.3f}")

# --- Figure: stability vs k ---
fig, ax = plt.subplots(figsize=(7, 5))
ax.errorbar(res["k"], res["val_stability"], yerr=res["val_stability_sem"],
            marker="o", capsize=3, color="#324095", label="Validation")
ax.plot(res["k"], res["train_stability"], marker="s", ls="--",
        color="#5CAEE1", label="Training")
ax.axvline(nbest, ls=":", color="grey", lw=1)
ax.set_xlabel("Number of clusters (k)")
ax.set_ylabel("Normalised stability (lower = better)")
ax.set_title(f"reval relative validation — best k = {nbest}")
ax.set_xticks(ks)
ax.legend()
plt.tight_layout()
plt.savefig(os.path.join(FIGURES_DIR, "reval_kselection.pdf"), dpi=300, bbox_inches="tight")
plt.close()

print("\nDone. Table: reval_kselection.csv | Figure: reval_kselection.pdf")
