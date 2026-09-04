#!/usr/bin/env python3
"""
Per-cluster Network-Based Statistic (NBS; Zalesky 2010) — cluster autism vs NT.

Cluster analogue of 6_functional_analysis/non_concat/scripts/07_nbs_test.py. For each
autism cluster (C1/C2/C3 vs the pooled NT group) and each direction
(hypo = cluster<NT, hyper = cluster>NT):
  1. Welch t per edge; primary edge threshold t* from p<NBS_P.
  2. Connected components among supra-threshold edges (156-node graph).
  3. Permute cluster/NT labels (N_PERM) -> null max component size -> FWER p.
  4. Keep the LARGEST observed component (the NBS subnetwork) and write its edges.

Outputs (outputs/tables/nbs/):
  nbs_summary.csv                          cluster x direction: size, n_nodes, FWER
  cluster_<C>_nbs_edges.csv                largest-component edges (region, direction, t)

These edge lists feed the brain-map plotters (cortical Schaefer + subcortical
yabplot). Uses the CURATED connectivity (df_conn_cohort_norm) and the k-means
cluster labels, diagnosis from population_group (curated).

Env: N_PERM (default 5000), NBS_P (default 0.01).
"""
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t as tdist
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from statsmodels.stats.multitest import multipletests

# reuse the section-6 loaders (connectivity + cluster merge + col filter)
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent
                       / "6_functional_analysis" / "non_concat" / "scripts"))
from _sensitivity_utils import load_connectivity_and_clusters, filter_connectivity_cols  # noqa: E402

_SECTION = Path(__file__).resolve().parent.parent
# NBS_SUBDIR lets each connectivity variant (e.g. the 6-2 v0.11 GSR modes) write
# to its own folder, e.g. NBS_SUBDIR=nbs_nogsr_concat. Connectivity source is
# selected upstream via FMRI_CONN_FILE (honoured by _sensitivity_utils).
OUT_DIR = _SECTION / "outputs" / "tables" / os.environ.get("NBS_SUBDIR", "nbs")
OUT_DIR.mkdir(parents=True, exist_ok=True)

N_PERM = int(os.environ.get("N_PERM", "5000"))
NBS_P  = float(os.environ.get("NBS_P", "0.05"))   # edge-forming / thresholding p
CLUSTERS = ["C1", "C2", "C3"]
RNG = np.random.default_rng(42)

# EDGE_MODE controls what gets written to the edge lists that feed the maps:
#   "uncorrected" (default) — ALL edges with uncorrected p<NBS_P, per direction.
#       No component restriction, no permutation. The figure caption states
#       "p<0.05 uncorrected". Fast.
#   "fdr" — edges with BH-FDR q<NBS_P (across all edges), split by direction.
#   "nbs" — Zalesky NBS: keep only the LARGEST connected component per direction
#       and estimate its FWER p by label permutation (N_PERM).
EDGE_MODE = os.environ.get("EDGE_MODE", "uncorrected")


def welch_t(B, mask_a):
    a, b = B[mask_a], B[~mask_a]
    va, vb = a.var(0, ddof=1), b.var(0, ddof=1)
    d = np.sqrt(va / a.shape[0] + vb / b.shape[0])
    d[d == 0] = np.nan
    return (a.mean(0) - b.mean(0)) / d


def _adj_labels(supra, ei, ej, n_nodes):
    si, sj = ei[supra], ej[supra]
    adj = coo_matrix((np.ones(si.size), (si, sj)), shape=(n_nodes, n_nodes))
    adj = adj + adj.T
    _, lab = connected_components(adj, directed=False)
    return lab, si, sj


def largest_component_size(supra, ei, ej, n_nodes):
    if not supra.any():
        return 0
    lab, si, _ = _adj_labels(supra, ei, ej, n_nodes)
    return int(np.bincount(lab[si]).max())


def largest_component_edges(supra, ei, ej, n_nodes):
    """Return boolean mask over edges belonging to the largest component."""
    if not supra.any():
        return np.zeros(supra.size, dtype=bool)
    lab, si, _ = _adj_labels(supra, ei, ej, n_nodes)
    biggest = np.bincount(lab[si]).argmax()
    edge_idx = np.where(supra)[0]
    keep = edge_idx[lab[si] == biggest]
    out = np.zeros(supra.size, dtype=bool)
    out[keep] = True
    return out


def main() -> int:
    df, conn_cols = load_connectivity_and_clusters()
    df = df.copy()
    df["dx"] = df["population_group"].replace("TD", "NT")

    pairs = [c[len("con_"):].split("/", 1) for c in conn_cols]
    nodes = sorted({n for p in pairs for n in p})
    nidx = {n: i for i, n in enumerate(nodes)}
    ei = np.array([nidx[s] for s, _ in pairs])
    ej = np.array([nidx[t] for _, t in pairs])
    n_nodes = len(nodes)

    nt = df["dx"] == "NT"
    print(f"EDGE_MODE={EDGE_MODE}, NBS_P={NBS_P}"
          + (f", N_PERM={N_PERM}" if EDGE_MODE == "nbs" else "")
          + f", NT n={int(nt.sum())}")
    summary = []

    for C in CLUSTERS:
        sub = df[((df["dx"] == "Autism") & (df["Cluster"] == C)) | nt]
        B = sub[conn_cols].values.astype(float)
        mask_a = (sub["dx"] == "Autism").values
        na, nb = int(mask_a.sum()), int((~mask_a).sum())
        t_thresh = tdist.ppf(1 - NBS_P / 2, na + nb - 2)
        t_obs = welch_t(B, mask_a)
        dirs = {"hypo": t_obs < 0, "hyper": t_obs > 0}

        # FDR mode: two-sided p per edge, BH across all edges, keep q<NBS_P
        fdr_pass = None
        if EDGE_MODE == "fdr":
            p_edge = 2.0 * tdist.sf(np.abs(t_obs), na + nb - 2)
            ok = np.isfinite(p_edge)
            fdr_pass = np.zeros(p_edge.size, dtype=bool)
            if ok.any():
                fdr_pass[ok] = multipletests(p_edge[ok], method="fdr_bh")[1] < NBS_P

        # NBS mode only: permutation null of the largest component size
        null_max = None
        if EDGE_MODE == "nbs":
            null_max = {"hypo": np.zeros(N_PERM), "hyper": np.zeros(N_PERM)}
            for k in range(N_PERM):
                tp = welch_t(B, RNG.permutation(mask_a))
                s = np.abs(tp) > t_thresh
                null_max["hypo"][k] = largest_component_size(s & (tp < 0), ei, ej, n_nodes)
                null_max["hyper"][k] = largest_component_size(s & (tp > 0), ei, ej, n_nodes)

        rows = []
        for d, sign in dirs.items():
            supra = (np.abs(t_obs) > t_thresh) & sign
            if EDGE_MODE == "nbs":
                keep = largest_component_edges(supra, ei, ej, n_nodes)
                size = int(keep.sum())
                fwer = ((1.0 + (null_max[d] >= size).sum()) / (N_PERM + 1.0)
                        if size else np.nan)
            elif EDGE_MODE == "fdr":  # BH q<NBS_P, split by direction
                keep = fdr_pass & sign
                size = int(keep.sum())
                fwer = np.nan
            else:  # uncorrected: every edge with p<NBS_P
                keep = supra
                size = int(keep.sum())
                fwer = np.nan
            comp_nodes = sorted(set(np.r_[ei[keep], ej[keep]].tolist()))
            for e in np.where(keep)[0]:
                rows.append({"region": conn_cols[e], "direction": d,
                             "t_stat": float(t_obs[e])})
            summary.append({"cluster": C, "direction": d, "comp_size_edges": size,
                            "n_nodes": len(comp_nodes), "p_fwer": fwer,
                            "mode": EDGE_MODE, "thresh_p": NBS_P, "na": na, "nb": nb})
            tag = {"nbs": f"FWER={fwer:.3f}", "fdr": "q<%.2f FDR" % NBS_P}.get(
                EDGE_MODE, "p<%.2f uncorr" % NBS_P)
            print(f"  {C} {d:5s}: {size} edges, {len(comp_nodes)} nodes, {tag}")
        pd.DataFrame(rows).to_csv(OUT_DIR / f"cluster_{C}_nbs_edges.csv", index=False)

    pd.DataFrame(summary).to_csv(OUT_DIR / "nbs_summary.csv", index=False)
    print(f"\nSaved edge lists + nbs_summary.csv to {OUT_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
