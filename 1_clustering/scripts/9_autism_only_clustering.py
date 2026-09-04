"""Autism-only clustering sensitivity analysis (Reviewer #4.6).

Question: if we cluster ONLY autistic participants on standardised IQ x SRS-2
(removing NT and IDD-only), does a similar three-cluster within-spectrum
structure emerge, and how concordant is it with the full-cohort assignment?

Method mirrors 04_run_clustering.py: K-means, k=3, n_init=25,
random_state=42, features standardised within the autistic subset, labels
relabelled by descending size (A1/A2/A3). Reads the same feature values the
full clustering used (the curated label file) so the two partitions are
directly comparable on the shared autistic individuals.

# ponytail: reuses the label file's IQ/SRS rather than re-loading raw TSVs,
# so this can't drift from the full-cohort run it is compared against.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (adjusted_rand_score, normalized_mutual_info_score,
                             silhouette_score)

ROOT = Path(__file__).resolve().parents[2]
LABELS = ROOT / "1_clustering/outputs/curated/tables/individuals_metrics_with_clusters_curated.csv"
OUT = ROOT / "1_clustering/outputs/curated/tables"
FIG_OUT = ROOT / "1_clustering/outputs/curated_autism/figures"
AUT = {"Autism without IDD", "Autism with IDD", "Autism to exclude"}
RS = 42

# Same palette + reference lines + axis limits as the main clustering figure
# (10_cluster_scatter_panels.py::draw_panel / config.PALETTE_CLUSTERS), so the
# autism-only figure reads as the same kind of plot, not a different style.
PALETTE = {"A1": "#7A8B47", "A2": "#ff9fa0", "A3": "#e7ba52"}
SRS_LINES = [60, 75]
IQ_LINES = [70, 130]
XLIM = (32, 94)   # matches the full-cohort curated SRS_tscore range (36-4, 90+4)
YLIM = (34, 161)  # matches the full-cohort curated IQ range (40-6, 155+6)


def plot_autism_only_scatter(aut, out_path):
    fig, ax = plt.subplots(figsize=(3.0, 3.0))
    for g in ["A1", "A2", "A3"]:
        sub = aut[aut["ClusterAutOnly"] == g]
        ax.scatter(sub["SRS_feat"], sub["IQ_feat"], s=9, color=PALETTE[g],
                   alpha=0.65, edgecolor="none", zorder=2, label=f"{g} (n={len(sub)})")
        ax.scatter(sub["SRS_feat"].mean(), sub["IQ_feat"].mean(), marker="o", s=90,
                   color=PALETTE[g], edgecolor="black", linewidth=1.5, zorder=5)
    for xv in SRS_LINES:
        ax.axvline(xv, color="#555555", ls=(0, (4, 3)), lw=0.9, zorder=4)
    for yv in IQ_LINES:
        ax.axhline(yv, color="#555555", ls=(0, (4, 3)), lw=0.9, zorder=4)
    ax.set_xlim(*XLIM); ax.set_ylim(*YLIM)
    ax.set_xticks(SRS_LINES); ax.set_yticks(IQ_LINES)
    ax.tick_params(labelsize=7)
    ax.set_aspect(abs((XLIM[1] - XLIM[0]) / (YLIM[1] - YLIM[0])))
    for s in ax.spines.values():
        s.set_edgecolor("#999999")
    ax.set_title("Autism-only k-means clusters", fontsize=11, fontweight="bold")
    ax.set_xlabel("SRS-2 t-score", fontsize=8)
    ax.set_ylabel("Measured IQ", fontsize=8)
    ax.legend(loc="upper right", fontsize=6.5, frameon=False,
              handletextpad=0.2, borderaxespad=0.2, labelspacing=0.25)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def relabel_by_size(lbl, k):
    order = list(pd.Series(lbl).value_counts().index[:k])
    remap = {old: new for new, old in enumerate(order)}
    return np.array([remap[x] for x in lbl])


def kmeans(X, k):
    km = KMeans(n_clusters=k, n_init=25, random_state=RS).fit(X)
    return relabel_by_size(km.labels_, k)


def main():
    df = pd.read_csv(LABELS)
    df = df[df["Cluster"].notna()].copy()
    # same feature definition as the full clustering: IQ x SRS-2 t-score
    df["IQ_feat"] = pd.to_numeric(df["total_IQ"], errors="coerce").fillna(
        pd.to_numeric(df["performance_IQ"], errors="coerce"))
    df["SRS_feat"] = pd.to_numeric(df["SRS_tscore"], errors="coerce")

    aut = df[df["population"].isin(AUT)].dropna(subset=["IQ_feat", "SRS_feat"]).copy()
    print(f"Autistic participants with IQ+SRS: {len(aut)}")

    X = StandardScaler().fit_transform(aut[["IQ_feat", "SRS_feat"]].values)

    # k selection among autistic-only (silhouette)
    print("\nk selection (autism-only, silhouette):")
    k_sil = []
    for k in range(2, 7):
        lab = kmeans(X, k)
        sil = silhouette_score(X, lab)
        print(f"  k={k}: silhouette={sil:.3f}")
        k_sil.append({"k": k, "silhouette": sil})

    # primary: k=3
    aut["ClusterAutOnly"] = [f"A{c+1}" for c in kmeans(X, 3)]

    print("\nAutism-only clusters (k=3) — IQ/SRS profile:")
    prof = aut.groupby("ClusterAutOnly").agg(
        n=("IQ_feat", "size"), IQ_mean=("IQ_feat", "mean"),
        IQ_sd=("IQ_feat", "std"), SRS_mean=("SRS_feat", "mean"),
        SRS_sd=("SRS_feat", "std")).round(1)
    print(prof.to_string())

    # concordance vs full-cohort assignment on the SAME autistic individuals
    ari = adjusted_rand_score(aut["Cluster"], aut["ClusterAutOnly"])
    nmi = normalized_mutual_info_score(aut["Cluster"], aut["ClusterAutOnly"])
    print(f"\nConcordance with full-cohort labels (autistic individuals only):")
    print(f"  ARI = {ari:.3f}   NMI = {nmi:.3f}")

    ct = pd.crosstab(aut["ClusterAutOnly"], aut["Cluster"])
    print("\nCross-tab (rows=autism-only A1-3, cols=full-cohort C1-3):")
    print(ct.to_string())

    # save
    pd.DataFrame(k_sil).to_csv(OUT / "autism_only_k_selection_silhouette.csv", index=False)
    prof.to_csv(OUT / "autism_only_cluster_profile.csv")
    ct.to_csv(OUT / "autism_only_vs_fullcohort_crosstab.csv")
    aut[["Cluster", "ClusterAutOnly", "IQ_feat", "SRS_feat", "population"]].to_csv(
        OUT / "autism_only_cluster_labels.csv", index=False)
    pd.DataFrame([{"ari": ari, "nmi": nmi, "n_autistic": len(aut)}]).to_csv(
        OUT / "autism_only_concordance.csv", index=False)
    print(f"\nwrote -> {OUT}/autism_only_*.csv")

    fig_path = FIG_OUT / "clustering_scatter_curated_autism.pdf"
    plot_autism_only_scatter(aut, fig_path)
    print(f"wrote -> {fig_path}")


if __name__ == "__main__":
    main()
