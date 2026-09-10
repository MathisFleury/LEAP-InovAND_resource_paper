#!/usr/bin/env python3
# =============================================================================
# 12 - Project SPARK (Litman et al. GFMM groups) into our IQ x SRS clusters
# =============================================================================
# Fits our primary clustering (k-means, k=3, IQ x SRS-2, n_init=25, seed 42) on
# the curated LEAP/INOVAND cohort, then assigns SPARK individuals to the nearest
# of OUR centroids using OUR scaler -- i.e. a pure projection, no refitting.
# The crosstab GFMM group x our cluster answers "where does each Litman group
# land in our space?".
#
# Coverage caveat: SPARK measures IQ on a small subset, so of the 5,392
# GFMM-labelled subjects only ~235 have both IQ and SRS-2 and can be projected.
# The projected subset is IQ-enriched and NOT representative of SPARK.
# =============================================================================

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import chi2_contingency
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

_HERE = os.path.dirname(os.path.abspath(__file__))
_OUT = os.path.normpath(os.path.join(_HERE, "..", "outputs"))
FIGURES_DIR = os.path.join(_OUT, "figures")
TABLES_DIR = os.path.join(_OUT, "tables")

SPARK_BASE = "/Volumes/Imaging5/EEG_MRI-MF/SPARK/_clinical_data"
SPARK_GFMM = f"{SPARK_BASE}/raw/cluster_asignment/GFMMgroupassignments_5392cohort.tsv"
SPARK_CLIN = f"{SPARK_BASE}/curated/SPARK_clinical_curated.tsv"

PALETTE = {"C1": "#7A8B47", "C2": "#ff9fa0", "C3": "#e7ba52"}
# GFMM group order = Litman et al.'s ordering (least -> most affected)
GFMM_ORDER = ["Social/Behavioral", "Moderate Challenges",
              "Mixed ASD with DD", "Broadly Impacted"]
GFMM_COLORS = dict(zip(GFMM_ORDER, ["#4c72b0", "#55a868", "#c44e52", "#8172b3"]))


# -----------------------------------------------------------------------------
# Our cohort -- refit the primary clustering to recover scaler + centroids
# -----------------------------------------------------------------------------
ours = pd.read_csv(os.path.join(TABLES_DIR, "cohort_curated.csv"), low_memory=False)
ours = ours.dropna(subset=["IQ", "SRS_tscore"])
Xo = ours[["SRS_tscore", "IQ"]].values  # plot order: x = SRS, y = IQ

scaler = StandardScaler().fit(Xo)
km = KMeans(n_clusters=3, n_init=25, random_state=42).fit(scaler.transform(Xo))
# Relabel by descending cluster size (largest -> C1), same as 4_run_clustering.py
order = list(pd.Series(km.labels_).value_counts().index)
remap = {old: f"C{new + 1}" for new, old in enumerate(order)}
ours["Cluster"] = [remap[l] for l in km.labels_]
print(f"Our cohort: n={len(ours)}, sizes={ours.Cluster.value_counts().sort_index().to_dict()}")


# -----------------------------------------------------------------------------
# SPARK -- merge GFMM labels with curated clinical, project onto our centroids
# -----------------------------------------------------------------------------
gfmm = pd.read_csv(SPARK_GFMM, sep="\t")
clin = pd.read_csv(SPARK_CLIN, sep="\t", low_memory=False)
sp = gfmm.merge(clin, left_on="subject_sp_id", right_on="ID", how="left")
sp["IQ"] = sp["total_IQ"].fillna(sp["performance_IQ"])  # same IQ rule as ours
sp = sp.dropna(subset=["IQ", "SRS_tscore"])
sp["Cluster"] = [remap[l] for l in km.predict(scaler.transform(sp[["SRS_tscore", "IQ"]].values))]
print(f"SPARK: {len(gfmm)} GFMM-labelled -> {len(sp)} projectable (IQ x SRS complete)")

ct = pd.crosstab(sp["mixed_pred"], sp["Cluster"]).reindex(GFMM_ORDER).fillna(0).astype(int)
pct_row = ct.div(ct.sum(axis=1), axis=0) * 100   # where does each GFMM group go?
pct_col = ct.div(ct.sum(axis=0), axis=1) * 100   # what is inside each of our clusters?
chi2, p, dof, _ = chi2_contingency(ct)
print(f"\n{ct}\n\nrow % (GFMM -> ours):\n{pct_row.round(1)}"
      f"\n\ncol % (inside our clusters):\n{pct_col.round(1)}"
      f"\n\nchi2={chi2:.1f}, dof={dof}, p={p:.2e}")

ct.assign(**{f"pct_of_gfmm_group_{c}": pct_row[c].round(1) for c in ct.columns},
          **{f"pct_of_cluster_{c}": pct_col[c].round(1) for c in ct.columns}).to_csv(
    os.path.join(TABLES_DIR, "spark_projection_crosstab.csv"))
sp[["subject_sp_id", "mixed_pred", "IQ", "SRS_tscore", "Cluster", "age_yrs", "Sex"]].to_csv(
    os.path.join(TABLES_DIR, "spark_projection_assignments.csv"), index=False)


# -----------------------------------------------------------------------------
# Figure: (A) our clusters + boundaries, (B) SPARK projected, (C) row-% bars
# -----------------------------------------------------------------------------
# Two rows (A|B on top, C full-width below) so the figure stays legible when
# scaled to text width in the response letter. ponytail: hand-tuned layout.
plt.rcParams.update({"font.size": 13, "axes.titlesize": 14, "axes.labelsize": 13,
                     "xtick.labelsize": 11, "ytick.labelsize": 11})
fig = plt.figure(figsize=(12.5, 11))
gs = fig.add_gridspec(2, 2, height_ratios=[1.35, 1], hspace=0.28, wspace=0.22)
axes = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, :])]

# Decision regions of our k-means, shared by panels A and B
xx, yy = np.meshgrid(np.linspace(25, 100, 300), np.linspace(15, 150, 300))
grid = np.c_[xx.ravel(), yy.ravel()]
zz = np.array([int(remap[l][1]) for l in km.predict(scaler.transform(grid))]).reshape(xx.shape)
region_cmap = [PALETTE["C1"], PALETTE["C2"], PALETTE["C3"]]
cent = scaler.inverse_transform(km.cluster_centers_)

for ax in axes[:2]:
    ax.contourf(xx, yy, zz, levels=[0.5, 1.5, 2.5, 3.5], colors=region_cmap, alpha=0.10)
    ax.contour(xx, yy, zz, levels=[1.5, 2.5], colors="0.5", linewidths=0.8)
    ax.set_xlabel("SRS-2 t-score")
    ax.set_ylabel("Full-scale IQ")
    ax.set_xlim(25, 100)
    ax.set_ylim(15, 150)

for c, g in ours.groupby("Cluster"):
    axes[0].scatter(g["SRS_tscore"], g["IQ"], c=PALETTE[c], s=18, alpha=0.65,
                    lw=0, label=f"{c} (n={len(g)})")
axes[0].scatter(cent[:, 0], cent[:, 1], marker="X", s=200, c="k", zorder=5)
axes[0].legend(title="Our cluster", frameon=False, loc="lower left", fontsize=11)
axes[0].set_title(f"A. LEAP/INOVAND clusters (k-means, n={len(ours)})")

for g_name in GFMM_ORDER:
    g = sp[sp["mixed_pred"] == g_name]
    axes[1].scatter(g["SRS_tscore"], g["IQ"], c=GFMM_COLORS[g_name], s=28,
                    alpha=0.8, lw=0.3, edgecolor="w", label=f"{g_name} (n={len(g)})")
axes[1].scatter(cent[:, 0], cent[:, 1], marker="X", s=200, c="k", zorder=5)
axes[1].legend(title="SPARK GFMM group", frameon=False, loc="lower left", fontsize=10)
axes[1].set_title(f"B. SPARK projected into our space (n={len(sp)})")

CLUSTERS = ["C1", "C2", "C3"]
bottom = np.zeros(len(CLUSTERS))
ypos = np.arange(len(CLUSTERS))
for g_name in GFMM_ORDER:
    vals = np.array([pct_col.loc[g_name, c] if c in pct_col else 0.0 for c in CLUSTERS])
    axes[2].barh(ypos, vals, left=bottom, color=GFMM_COLORS[g_name],
                 label=g_name, height=0.55)
    for y, (v, b) in enumerate(zip(vals, bottom)):
        if v >= 7:
            axes[2].text(b + v / 2, y, f"{v:.0f}%", ha="center", va="center",
                         fontsize=12, color="w")
    bottom += vals
axes[2].set_yticks(ypos)
axes[2].set_yticklabels([f"{c}\n(n={int(ct[c].sum()) if c in ct else 0})" for c in CLUSTERS])
axes[2].invert_yaxis()
axes[2].set_xlim(0, 100)
axes[2].set_xlabel("% of projected SPARK subjects in our cluster")
axes[2].legend(title="SPARK GFMM group", frameon=False, fontsize=11,
               loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=4)
axes[2].set_title(f"C. GFMM composition of our clusters (chi2={chi2:.1f}, p={p:.1e})")


out_pdf = os.path.join(FIGURES_DIR, "spark_projection.pdf")
plt.savefig(out_pdf, dpi=300, bbox_inches="tight")
plt.savefig(out_pdf.replace(".pdf", ".png"), dpi=200, bbox_inches="tight")
plt.close()
print(f"\nSaved: {out_pdf}")
