#!/usr/bin/env python3
"""
Export the largest NBS component's edges (per direction) as a source,target CSV
in the exact format 06_plot_double_network_matrix.R consumes, so the NBS result
can be rendered in the same double-triangle (count + proportion) style as
autism_vs_td_<dir>_double_matrix_final.pdf.

Recomputes the *observed* largest component only (Welch t + edge threshold +
connected components) — instant, no permutations. FWER significance comes from
the already-computed nbs_components_<dir>.csv, so 07 need not be rerun.

Writes, under NBS_DOUBLE_DIR (default: <AUTISM_TD_OUTPUT_DIR>/nbs_double):
    autism_vs_td_<dir>connectivity_full_data.csv   (source,target — for 05.R)
Then run 05 with AUTISM_TD_OUTPUT_DIR=<NBS_DOUBLE_DIR> to emit the PDF.
"""
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t as tdist
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _sensitivity_utils import filter_connectivity_cols, ATLAS_FILE  # noqa: E402

_SECTION = Path(__file__).resolve().parent.parent
FMRI_FILE = Path(os.environ.get(
    "AUTISM_TD_FMRI_CSV", _SECTION / "preprocessing" / "outputs" / "df_conn_cohort_norm.csv"))
OUT_DIR = Path(os.environ.get("AUTISM_TD_OUTPUT_DIR", _SECTION / "outputs" / "figures" / "revised"))
NBS_DOUBLE_DIR = Path(os.environ.get("NBS_DOUBLE_DIR", OUT_DIR / "nbs_double"))
NBS_DOUBLE_DIR.mkdir(parents=True, exist_ok=True)
NBS_P = float(os.environ.get("NBS_P", "0.01"))  # must match 07's edge-forming threshold


def welch_t(B, mask_a):
    a, b = B[mask_a], B[~mask_a]
    va, vb = a.var(0, ddof=1), b.var(0, ddof=1)
    denom = np.sqrt(va / a.shape[0] + vb / b.shape[0])
    denom[denom == 0] = np.nan
    return (a.mean(0) - b.mean(0)) / denom


def largest_component_edges(supra, ei, ej, n_nodes):
    """edge_idx of the largest connected component among supra-threshold edges."""
    if not supra.any():
        return np.array([], dtype=int)
    si, sj = ei[supra], ej[supra]
    adj = coo_matrix((np.ones(si.size), (si, sj)), shape=(n_nodes, n_nodes))
    adj = adj + adj.T
    _, lab = connected_components(adj, directed=False)
    top = np.bincount(lab[si]).argmax()
    return np.where(supra)[0][lab[si] == top]


def main() -> int:
    df = pd.read_csv(FMRI_FILE, low_memory=False)
    df["ID"] = df["ID"].astype(str)
    df["PopulationS1"] = df["population_group"].replace("TD", "NT")
    df = df[df["PopulationS1"].isin(["Autism", "NT"])].copy()

    conn_cols = filter_connectivity_cols(df.columns.tolist())
    pairs = [c[len("con_"):].split("/", 1) for c in conn_cols]
    nodes = sorted({n for p in pairs for n in p})
    nidx = {n: i for i, n in enumerate(nodes)}
    ei = np.array([nidx[s] for s, _ in pairs])
    ej = np.array([nidx[t] for _, t in pairs])

    B = df[conn_cols].values.astype(float)
    mask_a = (df["PopulationS1"] == "Autism").values
    na, nb = int(mask_a.sum()), int((~mask_a).sum())
    t_thresh = tdist.ppf(1 - NBS_P / 2, na + nb - 2)
    t_obs = welch_t(B, mask_a)
    directions = {"hypo": t_obs < 0, "hyper": t_obs > 0}  # Autism<NT / Autism>NT

    for d, sign_mask in directions.items():
        supra = (np.abs(t_obs) > t_thresh) & sign_mask
        idx = largest_component_edges(supra, ei, ej, len(nodes))
        edges = pd.DataFrame({"source": [nodes[ei[e]] for e in idx],
                              "target": [nodes[ej[e]] for e in idx]})
        out = NBS_DOUBLE_DIR / f"autism_vs_td_{d}connectivity_full_data.csv"
        edges.to_csv(out, index=False)
        print(f"[{d}] largest component: {len(edges)} edges -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
