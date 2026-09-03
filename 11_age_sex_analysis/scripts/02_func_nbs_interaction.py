#!/usr/bin/env python3
"""
NBS (Network-Based Statistic) for diagnosis x age and diagnosis x sex
interactions in functional connectivity.

Edge-level FDR over 10,731 edges is underpowered for an interaction, so this
runs a component-level, FWER-controlled test instead — the same NBS machinery
as 6_functional_analysis/scripts/07_nbs_test.py, but the per-edge statistic is
the *interaction* t from the model

    y ~ diag + age_c + Sex + cohort + diag:<modulator>        (modulator = age_c or Sex)

computed via Frisch-Waugh-Lovell: the full-model interaction t equals the t of
Y_perp on z_perp, where both the edges (Y) and the interaction regressor (z) are
residualised on the reduced design (everything except the interaction term).
Inference is by permuting the orthogonalised regressor z_perp (Smith/randomise
scheme) — exchangeable under H0 of no partial association and fully vectorised
across all edges, so N_PERM permutations cost one matrix-vector product each.

  1. t_e per edge (interaction); primary edge threshold t* from p<NBS_P.
  2. supra-threshold edges -> connected components in the 156-node graph;
     component size = #edges, separately for positive / negative interaction t.
  3. permute z_perp (N_PERM); keep largest component size per sign -> null max.
  4. FWER p per observed component = P(null max >= observed size).

Outputs (11_age_sex_analysis/outputs/nbs/):
  nbs_interaction_<model>_<sign>.csv     components + p_fwer + networks
  nbs_interaction_summary.csv            largest component + p_fwer per model/sign
"""
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t as tdist
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

_FUNC_SCRIPTS = '/Users/mfleury/POSTDOC/LIBRAIRY/LEAP-InovAND_resource/6_functional_analysis/scripts'
sys.path.insert(0, _FUNC_SCRIPTS)
from _sensitivity_utils import filter_connectivity_cols, ATLAS_FILE  # noqa: E402

_SECTION = Path(__file__).resolve().parent.parent
# Must match the manuscript's primary rsfMRI pipeline (verified against
# 6-2_functional_analysis/scripts/03_nbs_test.py and
# 7_cluster_functional_analysis/scripts/06_cluster_nbs.py, both run against
# this exact file for the reported whole-cohort/cluster NBS results: n=393
# autism/327 NT, matching the manuscript). Previously pointed at
# 6_functional_analysis's df_conn_cohort_norm.csv -- an older, non-run-
# concatenated table with a different sample (355/295) -- which was wrong.
FMRI_FILE = Path(os.environ.get(
    'AUTISM_TD_FMRI_CSV',
    '/Users/mfleury/POSTDOC/LIBRAIRY/LEAP-InovAND_resource/6-2_functional_analysis'
    '/preprocessing/outputs/nogsr_concat/df_conn_cohort_norm_nogsr_concat.csv'))
OUT_DIR = _SECTION / 'outputs' / 'nbs'
OUT_DIR.mkdir(parents=True, exist_ok=True)
N_PERM = int(os.environ.get('N_PERM', '5000'))
NBS_P = float(os.environ.get('NBS_P', '0.01'))   # primary edge-forming threshold
RNG = np.random.default_rng(42)


def residualise(M, X):
    """Residuals of columns of M (n x k) on design X (n x p). Least squares."""
    beta, _, _, _ = np.linalg.lstsq(X, M, rcond=None)
    return M - X @ beta


def largest_component_size(supra, ei, ej, n_nodes):
    if not supra.any():
        return 0
    si, sj = ei[supra], ej[supra]
    adj = coo_matrix((np.ones(si.size), (si, sj)), shape=(n_nodes, n_nodes))
    adj = adj + adj.T
    _, lab = connected_components(adj, directed=False)
    return int(np.bincount(lab[si]).max())


def components(supra, ei, ej, n_nodes):
    si, sj = ei[supra], ej[supra]
    adj = coo_matrix((np.ones(si.size), (si, sj)), shape=(n_nodes, n_nodes))
    adj = adj + adj.T
    _, lab = connected_components(adj, directed=False)
    out = []
    for c in np.unique(lab[si]):
        em = lab[si] == c
        out.append({'nodes': np.unique(np.r_[si[em], sj[em]]),
                    'edge_idx': np.where(supra)[0][em], 'size': int(em.sum())})
    return sorted(out, key=lambda d: -d['size'])


def partial_t(zp, Yp, df):
    """FWL interaction t per edge from orthogonalised z (n,) and Y (n x E)."""
    num = zp @ Yp                              # (E,)
    zz = zp @ zp
    beta = num / zz                            # (E,)
    resid = Yp - np.outer(zp, beta)            # (n x E)
    rss = (resid ** 2).sum(0)                  # (E,)
    se = np.sqrt(rss / df / zz)
    return beta / se


def run_model(df, conn_cols, ei, ej, n_nodes, net_of, nodes, model_name, modulator):
    n = len(df)
    diag = (df['diag'].values == 'Autism').astype(float)
    age_c = df['age_c'].values
    sex = df['Sex'].values
    # reduced design: intercept + diag + age_c + sex + cohort dummies (no interaction)
    cohort_d = pd.get_dummies(df['cohort'], drop_first=True).values.astype(float)
    Xred = np.column_stack([np.ones(n), diag, age_c, sex, cohort_d])
    zint = diag * (age_c if modulator == 'age' else sex)  # interaction regressor
    # full design adds zint -> df_resid = n - (Xred cols + 1)
    df_resid = n - (Xred.shape[1] + 1)

    Y = df[conn_cols].values.astype(float)
    zp = residualise(zint.reshape(-1, 1), Xred).ravel()
    Yp = residualise(Y, Xred)

    t_obs = partial_t(zp, Yp, df_resid)
    t_thresh = tdist.ppf(1 - NBS_P / 2, df_resid)
    print(f"  [{model_name}] {len(conn_cols)} edges, n={n}, df={df_resid}, "
          f"t*={t_thresh:.3f} (p<{NBS_P})")

    directions = {'pos': t_obs > 0, 'neg': t_obs < 0}
    obs = {d: components((np.abs(t_obs) > t_thresh) & m, ei, ej, n_nodes)
           for d, m in directions.items()}

    # export edges of the largest component per sign -> double-matrix plotter
    # (source,target node labels; p two-sided from the edge interaction t).
    pairs = [c[len('con_'):].split('/', 1) for c in conn_cols]
    for d in ('pos', 'neg'):
        comp = obs[d][0] if obs[d] else None
        rows = []
        if comp is not None:
            for e in comp['edge_idx']:
                s, t = pairs[e]
                pe = 2 * tdist.sf(abs(t_obs[e]), df_resid)
                rows.append({'source': s, 'target': t, 't_stat': t_obs[e],
                             'p_value': pe, 'p_value_fdr': np.nan})
        pd.DataFrame(rows, columns=['source', 'target', 't_stat', 'p_value', 'p_value_fdr']) \
          .to_csv(OUT_DIR / f'interaction_edges_{model_name}_{d}.csv', index=False)

    # permutation null: permute orthogonalised regressor -> recompute t per edge
    null_max = {'pos': np.zeros(N_PERM), 'neg': np.zeros(N_PERM)}
    for k in range(N_PERM):
        tp = partial_t(zp[RNG.permutation(n)], Yp, df_resid)
        supra = np.abs(tp) > t_thresh
        null_max['pos'][k] = largest_component_size(supra & (tp > 0), ei, ej, n_nodes)
        null_max['neg'][k] = largest_component_size(supra & (tp < 0), ei, ej, n_nodes)

    nodes_arr = np.array(nodes)
    summary = []
    for d in ('pos', 'neg'):
        rows = []
        for ci, comp in enumerate(obs[d], 1):
            p_fwer = (1.0 + (null_max[d] >= comp['size']).sum()) / (N_PERM + 1.0)
            comp_nodes = nodes_arr[comp['nodes']]
            nets = pd.Series([net_of.get(x) for x in comp_nodes]).value_counts().to_dict()
            rows.append({'model': f'diagnosis x {model_name}', 'sign': d,
                         'component': ci, 'n_edges': comp['size'],
                         'n_nodes': len(comp['nodes']), 'p_fwer': p_fwer,
                         'networks': ';'.join(f'{k}:{v}' for k, v in nets.items())})
        res = pd.DataFrame(rows)
        res.to_csv(OUT_DIR / f'nbs_interaction_{model_name}_{d}.csv', index=False)
        largest = int(res['n_edges'].max()) if len(res) else 0
        top_p = float(res.iloc[0]['p_fwer']) if len(res) else float('nan')
        n_sig = int((res['p_fwer'] < 0.05).sum()) if len(res) else 0
        print(f"    [{model_name}/{d}] {len(res)} components; largest={largest} edges; "
              f"top p_FWER={top_p:.4f}; FWER<0.05: {n_sig}")
        summary.append({'model': f'diagnosis x {model_name}', 'sign': d,
                        'n_components': len(res), 'largest_n_edges': largest,
                        'top_p_fwer': top_p, 'n_sig_fwer_05': n_sig})
    return summary


def main():
    print('=' * 60)
    print(f'NBS interaction analysis (functional)  N_PERM={N_PERM}  p*<{NBS_P}')
    print('=' * 60)
    df = pd.read_csv(FMRI_FILE, low_memory=False)
    df['diag'] = df['population_group'].replace('TD', 'NT')
    df = df[df['diag'].isin(['Autism', 'NT'])].copy()
    df = df[df['age_yrs'].notna() & df['Sex'].notna() & df['cohort'].notna()].copy()
    df['age_c'] = df['age_yrs'].astype(float) - df['age_yrs'].astype(float).mean()
    df['Sex'] = df['Sex'].astype(float)
    df['cohort'] = df['cohort'].map(str)
    print(f"n={len(df)}  {df['diag'].value_counts().to_dict()}")

    conn_cols = filter_connectivity_cols(df.columns.tolist())
    atlas = pd.read_csv(ATLAS_FILE, sep='\t')
    net_of = dict(zip(atlas['label'], atlas['network_label']))

    pairs = [c[len('con_'):].split('/', 1) for c in conn_cols]
    nodes = sorted({n for p in pairs for n in p})
    nidx = {n: i for i, n in enumerate(nodes)}
    ei = np.array([nidx[s] for s, _ in pairs])
    ej = np.array([nidx[t] for _, t in pairs])
    n_nodes = len(nodes)

    summary = []
    summary += run_model(df, conn_cols, ei, ej, n_nodes, net_of, nodes, 'age', 'age')
    summary += run_model(df, conn_cols, ei, ej, n_nodes, net_of, nodes, 'sex', 'sex')

    df_sum = pd.DataFrame(summary)
    df_sum.to_csv(OUT_DIR / 'nbs_interaction_summary.csv', index=False)
    print('\n' + '=' * 60)
    print('SUMMARY')
    print('=' * 60)
    print(df_sum.to_string(index=False))
    print(f"\nSaved: {OUT_DIR / 'nbs_interaction_summary.csv'}")


if __name__ == '__main__':
    raise SystemExit(main())
