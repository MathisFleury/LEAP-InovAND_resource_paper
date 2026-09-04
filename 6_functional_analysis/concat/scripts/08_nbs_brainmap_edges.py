#!/usr/bin/env python3
"""
Whole-group Autism-vs-NT NBS edges — for the figure_6b-style brain map.

Section-6-2 (XCP-D v0.11, merged thalamus → 2 ROIs) analogue of
7_cluster_functional_analysis/scripts/06_cluster_nbs.py, but the case side is
ALL autistic participants pooled (not split by cluster).

For each direction (hypo = Autism<NT, hyper = Autism>NT):
  1. Welch t per edge (Autism vs NT), primary edge threshold t* from p<NBS_P.
  2. Largest connected component among supra-threshold edges (Zalesky NBS).
  3. FWER p = P(null max >= observed size) over N_PERM label permutations.
  4. Write that component's edges.

Connectivity source selected upstream via FMRI_CONN_FILE (a 6-2 GSR variant);
output folder via NBS_SUBDIR (e.g. nbs_wholegroup_nogsr_concat).

Outputs (outputs/tables/<NBS_SUBDIR>/):
  nbs_edges.csv    region,direction,t_stat   (largest-component edges)
  nbs_summary.csv  direction,comp_size_edges,n_nodes,p_fwer,mode,thresh_p,na,nb

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

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _sensitivity_utils import filter_connectivity_cols  # noqa: E402

# Whole-group sample: read the connectivity file DIRECTLY (all Autism+NT), the
# same way 04_nbs_component_edges.py does — NOT the cluster/QC-merged loader,
# which drops subjects meant only for the cluster analysis. This keeps the FWER
# (here) on the same sample as the figures (04 -> 06/05/10/11).
_FMRI_FILE = os.environ.get(
    "FMRI_CONN_FILE",
    str(Path(__file__).resolve().parent.parent / "preprocessing" / "outputs"
        / "df_conn_cohort_norm.csv"))


def load_connectivity_and_clusters():
    df = pd.read_csv(_FMRI_FILE, low_memory=False)
    df["ID"] = df["ID"].astype(str)
    conn_cols = filter_connectivity_cols(df.columns.tolist())
    return df, conn_cols

_SECTION = Path(__file__).resolve().parent.parent
# Pure-table output (nbs_edges.csv, nbs_summary.csv) -- lives under outputs/tables/,
# a sibling of outputs/figures/ where 09/11 write the actual brain-map figures.
OUT_DIR = _SECTION / "outputs" / "tables" / os.environ.get("NBS_SUBDIR", "nbs_wholegroup")
OUT_DIR.mkdir(parents=True, exist_ok=True)

N_PERM = int(os.environ.get("N_PERM", "5000"))
NBS_P  = float(os.environ.get("NBS_P", "0.01"))
RNG = np.random.default_rng(42)


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
    df = df[df["dx"].isin(["Autism", "NT"])]

    pairs = [c[len("con_"):].split("/", 1) for c in conn_cols]
    nodes = sorted({n for p in pairs for n in p})
    nidx = {n: i for i, n in enumerate(nodes)}
    ei = np.array([nidx[s] for s, _ in pairs])
    ej = np.array([nidx[t] for _, t in pairs])
    n_nodes = len(nodes)

    B = df[conn_cols].values.astype(float)
    mask_a = (df["dx"] == "Autism").values
    na, nb = int(mask_a.sum()), int((~mask_a).sum())
    t_thresh = tdist.ppf(1 - NBS_P / 2, na + nb - 2)
    t_obs = welch_t(B, mask_a)
    dirs = {"hypo": t_obs < 0, "hyper": t_obs > 0}
    print(f"NBS whole-group (Autism {na} vs NT {nb}), {n_nodes} nodes, "
          f"{len(conn_cols)} edges, N_PERM={N_PERM}, primary p<{NBS_P} (t*={t_thresh:.3f})")

    null_max = {"hypo": np.zeros(N_PERM), "hyper": np.zeros(N_PERM)}
    for k in range(N_PERM):
        tp = welch_t(B, RNG.permutation(mask_a))
        s = np.abs(tp) > t_thresh
        null_max["hypo"][k] = largest_component_size(s & (tp < 0), ei, ej, n_nodes)
        null_max["hyper"][k] = largest_component_size(s & (tp > 0), ei, ej, n_nodes)

    rows, summary = [], []
    for d, sign in dirs.items():
        supra = (np.abs(t_obs) > t_thresh) & sign
        keep = largest_component_edges(supra, ei, ej, n_nodes)
        size = int(keep.sum())
        fwer = (1.0 + (null_max[d] >= size).sum()) / (N_PERM + 1.0) if size else np.nan
        comp_nodes = sorted(set(np.r_[ei[keep], ej[keep]].tolist()))
        for e in np.where(keep)[0]:
            rows.append({"region": conn_cols[e], "direction": d, "t_stat": float(t_obs[e])})
        summary.append({"direction": d, "comp_size_edges": size, "n_nodes": len(comp_nodes),
                        "p_fwer": fwer, "mode": "nbs", "thresh_p": NBS_P, "na": na, "nb": nb})
        print(f"  {d:5s}: largest component = {size} edges, {len(comp_nodes)} nodes, FWER={fwer:.4f}")

    pd.DataFrame(rows).to_csv(OUT_DIR / "nbs_edges.csv", index=False)
    pd.DataFrame(summary).to_csv(OUT_DIR / "nbs_summary.csv", index=False)
    print(f"Saved nbs_edges.csv + nbs_summary.csv -> {OUT_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
