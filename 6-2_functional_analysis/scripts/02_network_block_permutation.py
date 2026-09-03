#!/usr/bin/env python3
"""
Network-block permutation test for Autism vs NT connectivity.

Edge-level FDR is underpowered for distributed effects and ignores network
dependence. This adds a block-level inference that respects it:

  1. average each subject's edge connectivity within every network-pair block
     (11 networks -> 66 blocks, incl. within-network diagonal),
  2. observed Welch t (Autism vs NT) and mean difference per block,
  3. permutation null by shuffling the group labels (default 5000 perms),
  4. two-sided permutation p per block, BH-FDR across the 66 blocks.

So each network-block cell carries a real p-value instead of a raw count of
edge-level survivors.

Adds (never overwrites the edge-level double_matrix outputs):
    network_block_permutation_results.csv
    network_block_matrix_permutation.pdf
into AUTISM_TD_OUTPUT_DIR (default outputs/figures/revised).
Input connectivity via AUTISM_TD_FMRI_CSV (same env as 01); clusters/atlas fixed.
"""
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _sensitivity_utils import filter_connectivity_cols, ATLAS_FILE, DF_CLUSTERS_FILE  # noqa: E402

_SECTION = Path(__file__).resolve().parent.parent
FMRI_FILE = Path(os.environ.get(
    "AUTISM_TD_FMRI_CSV", _SECTION / "preprocessing" / "outputs" / "df_conn_cohort_norm.csv"))
OUT_DIR = Path(os.environ.get("AUTISM_TD_OUTPUT_DIR", _SECTION / "outputs" / "figures" / "revised"))
OUT_DIR.mkdir(parents=True, exist_ok=True)
N_PERM = int(os.environ.get("N_PERM", "5000"))
RNG = np.random.default_rng(42)

NET_ORDER = ["Amyg. & Hippoc.", "Striatum", "Cerebellum", "Thalamus", "Default",
             "Cont", "Limbic", "SalVentAttn", "DorsAttn", "SomMot", "Vis"]
READABLE = {"Amyg. & Hippoc.": "Amyg. & Hippoc.", "Striatum": "Striatum",
            "Cerebellum": "Cerebellum", "Thalamus": "Thalamus", "Default": "Default Mode",
            "Cont": "Frontoparietal", "Limbic": "Limbic", "SalVentAttn": "Ventral Attention",
            "DorsAttn": "Dorsal Attention", "SomMot": "Somatomotor", "Vis": "Visual"}


def welch_t(block_means, mask_a):
    """Per-block Welch t (group A vs B) across columns. block_means: (n_subj, n_block)."""
    a, b = block_means[mask_a], block_means[~mask_a]
    na, nb = a.shape[0], b.shape[0]
    va, vb = a.var(0, ddof=1), b.var(0, ddof=1)
    denom = np.sqrt(va / na + vb / nb)
    denom[denom == 0] = np.nan
    return (a.mean(0) - b.mean(0)) / denom


def main() -> int:
    print(f"Connectivity: {FMRI_FILE}")
    df = pd.read_csv(FMRI_FILE, low_memory=False)
    df["ID"] = df["ID"].astype(str)
    # Diagnosis from the curated clinical merge (population_group) carried in the
    # preprocessing output, not the paper PopulationS1 from the cluster file.
    df["PopulationS1"] = df["population_group"].replace("TD", "NT")
    df = df[df["PopulationS1"].isin(["Autism", "NT"])].copy()
    print(f"Autism/NT N: {df['PopulationS1'].value_counts().to_dict()}")

    conn_cols = filter_connectivity_cols(df.columns.tolist())
    atlas = pd.read_csv(ATLAS_FILE, sep="\t")
    net_of = dict(zip(atlas["label"], atlas["network_label"]))

    # map each edge column -> network-pair block (unordered), keep blocks in NET_ORDER
    order_idx = {n: i for i, n in enumerate(NET_ORDER)}
    block_cols: dict[tuple, list] = {}
    for c in conn_cols:
        s, t = c[len("con_"):].split("/", 1)
        ns, nt = net_of.get(s), net_of.get(t)
        if ns not in order_idx or nt not in order_idx:
            continue
        key = tuple(sorted((ns, nt), key=lambda n: order_idx[n]))
        block_cols.setdefault(key, []).append(c)

    blocks = sorted(block_cols, key=lambda k: (order_idx[k[0]], order_idx[k[1]]))
    # subject x block mean connectivity
    B = np.column_stack([df[block_cols[b]].mean(axis=1).values for b in blocks])
    mask_a = (df["PopulationS1"] == "Autism").values

    t_obs = welch_t(B, mask_a)
    mean_a = B[mask_a].mean(0)
    mean_b = B[~mask_a].mean(0)
    mean_diff = mean_a - mean_b

    # permutation null
    null_ge = np.zeros(len(blocks))
    for _ in range(N_PERM):
        t_perm = welch_t(B, RNG.permutation(mask_a))
        null_ge += (np.abs(t_perm) >= np.abs(t_obs)).astype(float)
    p_perm = (1.0 + null_ge) / (N_PERM + 1.0)
    p_fdr = multipletests(p_perm, method="fdr_bh")[1]

    res = pd.DataFrame({
        "network_1": [READABLE[b[0]] for b in blocks],
        "network_2": [READABLE[b[1]] for b in blocks],
        "n_edges": [len(block_cols[b]) for b in blocks],
        "mean_autism": mean_a, "mean_nt": mean_b, "mean_diff": mean_diff,
        "t_obs": t_obs, "p_perm": p_perm, "p_fdr": p_fdr,
        "n_autism": int(mask_a.sum()), "n_nt": int((~mask_a).sum()),
    }).sort_values("p_fdr")
    csv = OUT_DIR / "network_block_permutation_results.csv"
    res.to_csv(csv, index=False)
    print(f"Wrote {csv}  ({(p_fdr < 0.05).sum()} blocks FDR<0.05, {N_PERM} perms)")

    # ---- matrix figure: signed effect, FDR-significant blocks asterisked ----
    N = len(NET_ORDER)
    M = np.full((N, N), np.nan)
    sig = np.zeros((N, N), bool)
    for k, b in enumerate(blocks):
        i, j = order_idx[b[0]], order_idx[b[1]]
        M[i, j] = M[j, i] = mean_diff[k]
        if p_fdr[k] < 0.05:
            sig[i, j] = sig[j, i] = True

    plt.rcParams["font.family"] = "Helvetica"
    vmax = np.nanmax(np.abs(M))
    fig, ax = plt.subplots(figsize=(8.5, 7.5))
    im = ax.imshow(M, cmap="RdBu_r", vmin=-vmax, vmax=vmax, origin="lower")
    labs = [READABLE[n] for n in NET_ORDER]
    ax.set_xticks(range(N)); ax.set_xticklabels(labs, rotation=45, ha="right", fontsize=9)
    ax.set_yticks(range(N)); ax.set_yticklabels(labs, fontsize=9)
    for i in range(N):
        for j in range(N):
            if sig[i, j]:
                ax.text(j, i, "*", ha="center", va="center", fontsize=14, fontweight="bold")
    ax.set_title("Autism vs NT — network-block connectivity\n"
                 f"(mean difference; * permutation FDR<0.05, {N_PERM} perms)", fontsize=12)
    cb = fig.colorbar(im, ax=ax, shrink=0.7)
    cb.set_label("Mean connectivity difference (Autism − NT)", fontsize=9)
    fig.tight_layout()
    out_pdf = OUT_DIR / "network_block_matrix_permutation.pdf"
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {out_pdf}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
