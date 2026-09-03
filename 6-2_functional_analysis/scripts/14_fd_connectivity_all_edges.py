#!/usr/bin/env python3
"""
Mean-FD vs. connectivity, ALL edges (not just the 2 significant NBS
subnetworks) -- functional counterpart of 4_anatomical_analysis/curated/
scripts/18_euler_phenotype_all_features.py, for a combined sMRI+fMRI
image-quality effect-size figure (per user request, 2026-09-01).

For every edge (9,043), correlates mean_fd against connectivity, pooled
(whole Autism+NT sample, nogsr_concat, n=720) -- reported as an EQUIVALENT
COHEN'S D (d = 2r/sqrt(1-r^2)), same convention as the structural check, so
both modalities' confound checks are expressed in the same units. Aggregated
to a network x network mean-d matrix (same NET_ORDER/atlas as
03_plot_func_interaction_double_matrix.R) for a network-level summary; also
flags which network pairs contain >=1 edge from the significant whole-cohort
NBS component (hyper or hypo), mirroring the structural figure's Bonferroni-
significant-region outline.

Output: outputs/qc_motion/fd_conn_all_edges_network_matrix.csv
        outputs/qc_motion/fd_conn_all_edges_sig_pairs.csv
"""
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _sensitivity_utils import filter_connectivity_cols, ATLAS_FILE  # noqa: E402

_SECTION = Path(__file__).resolve().parent.parent
VARIANT = os.environ.get("XCPD_VARIANT_TAG", "nogsr_concat")
CONN = _SECTION / "preprocessing" / "outputs" / VARIANT / f"df_conn_cohort_norm_{VARIANT}.csv"
NBS_DIR = _SECTION / "outputs" / VARIANT / "nbs_double"
OUT = _SECTION / "outputs" / "qc_motion"
OUT.mkdir(parents=True, exist_ok=True)

NET_ORDER = ["Amyg. & Hippoc.", "Striatum", "Cerebellum", "Thalamus",
             "Default", "Cont", "Limbic", "SalVentAttn", "DorsAttn", "SomMot", "Vis"]


def main():
    df = pd.read_csv(CONN, low_memory=False)
    df["dx"] = df["population_group"].replace("TD", "NT")
    df = df[df["dx"].isin(["Autism", "NT"])].copy()
    conn_cols = filter_connectivity_cols(df.columns.tolist())
    fd = df["mean_fd"].values

    atlas = pd.read_csv(ATLAS_FILE, sep="\t")
    net_of = dict(zip(atlas["label"], atlas["network_label"]))

    # significant whole-cohort NBS edges (both directions)
    sig_edges = set()
    for direction in ("hyper", "hypo"):
        fp = NBS_DIR / f"autism_vs_td_{direction}connectivity_full_data.csv"
        if fp.exists():
            e = pd.read_csv(fp)
            for s, t in zip(e["source"], e["target"]):
                sig_edges.add(frozenset((s, t)))
    print(f"n={len(df)}  edges={len(conn_cols)}  significant edges (hyper+hypo): {len(sig_edges)}")

    rows = []
    nets = [n for n in NET_ORDER if n in set(net_of.values())]
    net_idx = {n: i for i, n in enumerate(nets)}
    sum_d = np.zeros((len(nets), len(nets)))
    n_edges_pp = np.zeros((len(nets), len(nets)))
    has_sig = np.zeros((len(nets), len(nets)), dtype=bool)

    for col in conn_cols:
        src, tgt = col[len("con_"):].split("/", 1)
        ns, nt = net_of.get(src), net_of.get(tgt)
        if ns not in net_idx or nt not in net_idx:
            continue
        y = df[col].values
        ok = np.isfinite(fd) & np.isfinite(y)
        if ok.sum() < 10:
            continue
        r = np.corrcoef(fd[ok], y[ok])[0, 1]
        d_equiv = 2 * r / np.sqrt(1 - r ** 2) if abs(r) < 1 else np.nan
        i, j = net_idx[ns], net_idx[nt]
        sum_d[i, j] += d_equiv; n_edges_pp[i, j] += 1
        if i != j:
            sum_d[j, i] += d_equiv; n_edges_pp[j, i] += 1
        if frozenset((src, tgt)) in sig_edges:
            has_sig[i, j] = has_sig[j, i] = True
        rows.append({"source": src, "target": tgt, "net_source": ns, "net_target": nt,
                     "r": round(r, 4), "cohens_d": round(d_equiv, 4)})

    edge_df = pd.DataFrame(rows)
    edge_df.to_csv(OUT / "fd_conn_all_edges.csv", index=False)

    mean_d = np.divide(sum_d, n_edges_pp, out=np.full_like(sum_d, np.nan), where=n_edges_pp > 0)
    mat = pd.DataFrame(mean_d, index=nets, columns=nets)
    mat.to_csv(OUT / "fd_conn_all_edges_network_matrix.csv")
    sig_mat = pd.DataFrame(has_sig, index=nets, columns=nets)
    sig_mat.to_csv(OUT / "fd_conn_all_edges_sig_pairs.csv")

    print(f"Overall: mean|d_equiv|={edge_df['cohens_d'].abs().mean():.3f}  "
          f"max|d_equiv|={edge_df['cohens_d'].abs().max():.3f}")
    print(f"Saved: {OUT}/fd_conn_all_edges_network_matrix.csv, fd_conn_all_edges_sig_pairs.csv")


if __name__ == "__main__":
    main()
