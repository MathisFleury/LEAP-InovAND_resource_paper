#!/usr/bin/env python3
"""
Network-Based Statistic (NBS; Zalesky et al. 2010) for Autism vs NT connectivity.

Complements the edge-level FDR and the network-block permutation test with a
component-level, FWER-controlled inference that exploits the network topology:

  1. Welch t per edge (Autism vs NT); fix a primary edge threshold t* from p<NBS_P.
  2. Among supra-threshold edges (separately for Autism>NT and Autism<NT), find
     connected components in the 156-node graph; component size = #edges (extent).
  3. Permute group labels (default 5000); each permutation keeps the largest
     component size per direction -> null max distribution.
  4. FWER p per observed component = P(null max >= observed size).

Adds (never overwrites): per AUTISM_TD_OUTPUT_DIR (default outputs/figures/revised)
    nbs_components_<dir>.csv          significant components + edges/networks
    nbs_block_counts_<dir>.pdf        network-block edge counts of sig components
Input connectivity via AUTISM_TD_FMRI_CSV (same env as 01).
"""
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import t as tdist
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _sensitivity_utils import filter_connectivity_cols, ATLAS_FILE, DF_CLUSTERS_FILE  # noqa: E402

_SECTION = Path(__file__).resolve().parent.parent
FMRI_FILE = Path(os.environ.get(
    "AUTISM_TD_FMRI_CSV", _SECTION / "preprocessing" / "outputs" / "df_conn_cohort_norm.csv"))
_BASE_DIR = Path(os.environ.get("AUTISM_TD_OUTPUT_DIR", _SECTION / "outputs" / "figures" / "revised"))
TABLES_DIR = _BASE_DIR / "tables"
FIGURES_DIR = _BASE_DIR / "figures"
TABLES_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)
N_PERM = int(os.environ.get("N_PERM", "5000"))
NBS_P = float(os.environ.get("NBS_P", "0.01"))   # primary edge-forming threshold
RNG = np.random.default_rng(42)

NET_ORDER = ["Amyg. & Hippoc.", "Striatum", "Cerebellum", "Thalamus", "Default",
             "Cont", "Limbic", "SalVentAttn", "DorsAttn", "SomMot", "Vis"]
READABLE = {"Amyg. & Hippoc.": "Amyg. & Hippoc.", "Striatum": "Striatum",
            "Cerebellum": "Cerebellum", "Thalamus": "Thalamus", "Default": "Default Mode",
            "Cont": "Frontoparietal", "Limbic": "Limbic", "SalVentAttn": "Ventral Attention",
            "DorsAttn": "Dorsal Attention", "SomMot": "Somatomotor", "Vis": "Visual"}


def welch_t(B, mask_a):
    a, b = B[mask_a], B[~mask_a]
    na, nb = a.shape[0], b.shape[0]
    va, vb = a.var(0, ddof=1), b.var(0, ddof=1)
    denom = np.sqrt(va / na + vb / nb)
    denom[denom == 0] = np.nan
    return (a.mean(0) - b.mean(0)) / denom


def largest_component_size(supra, ei, ej, n_nodes):
    """Max #edges in any connected component among supra-threshold edges."""
    if not supra.any():
        return 0
    si, sj = ei[supra], ej[supra]
    adj = coo_matrix((np.ones(si.size), (si, sj)), shape=(n_nodes, n_nodes))
    adj = adj + adj.T
    _, lab = connected_components(adj, directed=False)
    # component of each supra edge = lab[endpoint]; count edges per component
    return int(np.bincount(lab[si]).max())


def components(supra, ei, ej, n_nodes):
    si, sj = ei[supra], ej[supra]
    adj = coo_matrix((np.ones(si.size), (si, sj)), shape=(n_nodes, n_nodes))
    adj = adj + adj.T
    ncomp, lab = connected_components(adj, directed=False)
    out = []
    for c in np.unique(lab[si]):
        em = lab[si] == c
        out.append({"nodes": np.unique(np.r_[si[em], sj[em]]),
                    "edge_idx": np.where(supra)[0][em], "size": int(em.sum())})
    return sorted(out, key=lambda d: -d["size"])


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

    # node indexing + edge endpoint arrays
    pairs = [c[len("con_"):].split("/", 1) for c in conn_cols]
    nodes = sorted({n for p in pairs for n in p})
    nidx = {n: i for i, n in enumerate(nodes)}
    ei = np.array([nidx[s] for s, _ in pairs])
    ej = np.array([nidx[t] for _, t in pairs])
    n_nodes = len(nodes)

    B = df[conn_cols].values.astype(float)
    mask_a = (df["PopulationS1"] == "Autism").values
    na, nb = int(mask_a.sum()), int((~mask_a).sum())
    t_thresh = tdist.ppf(1 - NBS_P / 2, na + nb - 2)
    print(f"{len(conn_cols)} edges, {n_nodes} nodes, primary t*={t_thresh:.3f} (p<{NBS_P})")

    t_obs = welch_t(B, mask_a)
    directions = {"hypo": t_obs < 0, "hyper": t_obs > 0}   # Autism<NT / Autism>NT

    # observed components per direction
    obs = {d: components((np.abs(t_obs) > t_thresh) & sign_mask, ei, ej, n_nodes)
           for d, sign_mask in directions.items()}

    # permutation null: max component size per direction
    null_max = {"hypo": np.zeros(N_PERM), "hyper": np.zeros(N_PERM)}
    for k in range(N_PERM):
        tp = welch_t(B, RNG.permutation(mask_a))
        supra = np.abs(tp) > t_thresh
        null_max["hypo"][k] = largest_component_size(supra & (tp < 0), ei, ej, n_nodes)
        null_max["hyper"][k] = largest_component_size(supra & (tp > 0), ei, ej, n_nodes)
    np.savez(TABLES_DIR / "nbs_null_max.npz", **null_max)

    READABLE_DIR = {"hypo": "Autism < NT (hypo)", "hyper": "Autism > NT (hyper)"}

    nodes_arr = np.array(nodes)
    for d in ("hypo", "hyper"):
        rows = []
        for ci, comp in enumerate(obs[d], 1):
            p_fwer = (1.0 + (null_max[d] >= comp["size"]).sum()) / (N_PERM + 1.0)
            comp_nodes = nodes_arr[comp["nodes"]]
            nets = pd.Series([net_of.get(x) for x in comp_nodes]).value_counts().to_dict()
            rows.append({"direction": d, "component": ci, "n_edges": comp["size"],
                         "n_nodes": len(comp["nodes"]), "p_fwer": p_fwer,
                         "networks": ";".join(f"{k}:{v}" for k, v in nets.items())})
        res = pd.DataFrame(rows)
        res.to_csv(TABLES_DIR / f"nbs_components_{d}.csv", index=False)
        sig = res[res["p_fwer"] < 0.05] if len(res) else res
        print(f"[{d}] {len(res)} components; largest={res['n_edges'].max() if len(res) else 0} edges; "
              f"FWER<0.05: {len(sig)}")

        # --- NBS diagnostic plots (always, regardless of significance) ---
        if not len(res):
            print(f"    [{d}] no supra-threshold components — nothing to plot")
            continue
        comp = obs[d][0]
        p_top = float(res.iloc[0]["p_fwer"])
        sig_txt = "significant" if p_top < 0.05 else "n.s."
        col = "firebrick" if d == "hyper" else "navy"
        plt.rcParams["font.family"] = "Helvetica"

        # (1) permutation null of the largest component size vs observed
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.hist(null_max[d], bins=40, color="0.75", edgecolor="white")
        ax.axvline(comp["size"], color=col, lw=2,
                   label=f"observed largest = {comp['size']} edges")
        ax.set_xlabel("Largest component size under H₀ (no. of edges)")
        ax.set_ylabel(f"Permutations (N={N_PERM})")
        ax.set_title(f"NBS null distribution — {READABLE_DIR[d]}\n"
                     f"p_FWER = {p_top:.4f} ({sig_txt})", fontsize=11)
        ax.legend(fontsize=8, frameon=False)
        fig.tight_layout()
        fig.savefig(FIGURES_DIR / f"nbs_null_distribution_{d}.pdf", dpi=300, bbox_inches="tight")
        plt.close(fig)

        # (2) network-block edge-count matrix of the largest observed component
        M = np.zeros((len(NET_ORDER), len(NET_ORDER)))
        oi = {n: i for i, n in enumerate(NET_ORDER)}
        for e in comp["edge_idx"]:
            ns, nt = net_of.get(nodes[ei[e]]), net_of.get(nodes[ej[e]])
            if ns in oi and nt in oi:
                M[oi[ns], oi[nt]] += 1
                M[oi[nt], oi[ns]] += 1
        fig, ax = plt.subplots(figsize=(8, 7))
        im = ax.imshow(np.where(M > 0, M, np.nan), cmap="Reds" if d == "hyper" else "Blues",
                       origin="lower")
        labs = [READABLE[n] for n in NET_ORDER]
        ax.set_xticks(range(len(NET_ORDER))); ax.set_xticklabels(labs, rotation=45, ha="right", fontsize=9)
        ax.set_yticks(range(len(NET_ORDER))); ax.set_yticklabels(labs, fontsize=9)
        ax.set_title(f"NBS {d}-connectivity largest component (Autism vs NT)\n"
                     f"{comp['size']} edges, p_FWER={p_top:.4f} ({sig_txt})", fontsize=12)
        cb = fig.colorbar(im, ax=ax, shrink=0.7); cb.set_label("No. of component edges", fontsize=9)
        fig.tight_layout()
        fig.savefig(FIGURES_DIR / f"nbs_block_counts_{d}.pdf", dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"    saved nbs_null_distribution_{d}.pdf + nbs_block_counts_{d}.pdf")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
