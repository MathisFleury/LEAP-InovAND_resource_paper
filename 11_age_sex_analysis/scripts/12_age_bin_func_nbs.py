#!/usr/bin/env python3
"""
Age-binned NBS for functional connectivity Autism-vs-NT differences --
functional analogue of 10_age_bin_brain_maps.py, and a direct answer to R1's
actual question (Uddin et al. 2013, "developmental perspective") for the
modality that paper is about. Replaces the linear diagnosis x age NBS
interaction test (02_func_nbs_interaction.py) with a per-age-bin main-effect
NBS: for each bin, is there an Autism-vs-NT connectivity difference, adjusted
for scanner [+ motion] [+ sex]? A stratified view catches non-monotonic/
stage-specific patterns a single linear interaction term would miss.

Same bins as the structural analysis (<6, 6-9, 10-13, 14-17, 18-24, 25+;
2026-08-31, for consistency across modalities). <6 has zero NT subjects with
usable resting-state connectivity (fMRI is essentially unobtainable below
age 6 in this cohort) -- untestable, not just sparse -- so it is skipped and
reported as "no NT data" rather than run through NBS.
Adjusts for `machine_batch` (the individual-scanner ComBat batch variable,
strictly nested within cohort -- e.g. in the 18-24 bin Autism is dominated by
LEAP batches 7/13 while NT skews toward batch 2/3/8; 25+ similarly skews by
scanner) rather than the coarser 2-level cohort dummy, plus `mean_fd` (autism
has higher motion than NT specifically in the 10-13/14-17/18-24 bins, echoing
the whole-cohort meanFD confound already flagged for reviewer #4). Per-edge
model is

    edge ~ diag + machine_batch + mean_fd [+ sex]

fit via Frisch-Waugh-Lovell, reusing 02_func_nbs_interaction.py's
residualise/partial_t/components/largest_component_size machinery with
z = diag instead of a diag x modulator interaction regressor -- NOT the
plain unadjusted Welch-t of 6_functional_analysis/non_concat/scripts/07_nbs_test.py.

Outputs (11_age_sex_analysis/outputs/tables/nbs/):
  age_bin_nbs_<bin>_<sign>.csv         components + p_fwer + networks, per bin
  age_bin_nbs_edges_<bin>_<sign>.csv   edges of the top component, per bin
  age_bin_nbs_summary.csv              largest component + p_fwer per bin/sign
"""
import os
import sys
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t as tdist

_SCRIPT_DIR = Path(__file__).resolve().parent
_SECTION = _SCRIPT_DIR.parent
OUT_DIR = _SECTION / 'outputs' / 'tables' / 'nbs'  # pure-table step: no plots
OUT_DIR.mkdir(parents=True, exist_ok=True)

_spec = importlib.util.spec_from_file_location('nbs02', _SCRIPT_DIR / '02_func_nbs_interaction.py')
nbs = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(nbs)

sys.path.insert(0, str(_SECTION.parent / '6_functional_analysis' / 'non_concat' / 'scripts'))
from _sensitivity_utils import filter_connectivity_cols, ATLAS_FILE  # noqa: E402

N_PERM = int(os.environ.get('N_PERM', '5000'))
NBS_P = float(os.environ.get('NBS_P', '0.01'))   # primary edge-forming threshold
RNG = np.random.default_rng(42)

BIN_EDGES = [0, 6, 10, 14, 18, 25, 200]
BIN_LABELS = ['<6', '6-9', '10-13', '14-17', '18-24', '25+']
SAFE = {'<6': 'lt6', '6-9': '6_9', '10-13': '10_13', '14-17': '14_17',
        '18-24': '18_24', '25+': '25plus'}
MIN_N = 2  # below this in either group, the bin is untestable (not just sparse)


def run_bin(df, conn_cols, pairs, ei, ej, n_nodes, net_of, nodes, bin_label):
    n = len(df)
    diag = (df['diag'].values == 'Autism').astype(float)
    cols = [np.ones(n)]
    if df['machine_batch'].nunique() > 1:
        cols.append(pd.get_dummies(df['machine_batch'], drop_first=True).values.astype(float))
    cols.append(df['mean_fd'].values.reshape(-1, 1).astype(float))
    if df['Sex'].nunique() > 1:
        cols.append(df['Sex'].values.reshape(-1, 1))
    Xred = np.column_stack(cols)
    df_resid = n - Xred.shape[1] - 1

    Y = df[conn_cols].values.astype(float)
    zp = nbs.residualise(diag.reshape(-1, 1), Xred).ravel()
    Yp = nbs.residualise(Y, Xred)

    t_obs = nbs.partial_t(zp, Yp, df_resid)
    t_thresh = tdist.ppf(1 - NBS_P / 2, df_resid)
    n_a, n_n = int(diag.sum()), int((1 - diag).sum())
    print(f"  [{bin_label}] n={n} (A={n_a}, NT={n_n}), df_resid={df_resid}, t*={t_thresh:.3f}")

    directions = {'pos': t_obs > 0, 'neg': t_obs < 0}
    obs = {d: nbs.components((np.abs(t_obs) > t_thresh) & m, ei, ej, n_nodes)
           for d, m in directions.items()}

    for d in ('pos', 'neg'):
        comp = obs[d][0] if obs[d] else None
        rows = []
        if comp is not None:
            for e in comp['edge_idx']:
                s, t = pairs[e]
                pe = 2 * tdist.sf(abs(t_obs[e]), df_resid)
                rows.append({'source': s, 'target': t, 't_stat': t_obs[e], 'p_value': pe})
        pd.DataFrame(rows, columns=['source', 'target', 't_stat', 'p_value']).to_csv(
            OUT_DIR / f"age_bin_nbs_edges_{SAFE[bin_label]}_{d}.csv", index=False)

    null_max = {'pos': np.zeros(N_PERM), 'neg': np.zeros(N_PERM)}
    for k in range(N_PERM):
        tp = nbs.partial_t(zp[RNG.permutation(n)], Yp, df_resid)
        supra = np.abs(tp) > t_thresh
        null_max['pos'][k] = nbs.largest_component_size(supra & (tp > 0), ei, ej, n_nodes)
        null_max['neg'][k] = nbs.largest_component_size(supra & (tp < 0), ei, ej, n_nodes)

    nodes_arr = np.array(nodes)
    summary = []
    for d in ('pos', 'neg'):
        rows = []
        for ci, comp in enumerate(obs[d], 1):
            p_fwer = (1.0 + (null_max[d] >= comp['size']).sum()) / (N_PERM + 1.0)
            comp_nodes = nodes_arr[comp['nodes']]
            nets = pd.Series([net_of.get(x) for x in comp_nodes]).value_counts().to_dict()
            rows.append({'bin': bin_label, 'sign': d, 'component': ci, 'n_edges': comp['size'],
                         'n_nodes': len(comp['nodes']), 'p_fwer': p_fwer,
                         'networks': ';'.join(f'{k}:{v}' for k, v in nets.items())})
        res = pd.DataFrame(rows)
        res.to_csv(OUT_DIR / f"age_bin_nbs_{SAFE[bin_label]}_{d}.csv", index=False)
        largest = int(res['n_edges'].max()) if len(res) else 0
        top_p = float(res.iloc[0]['p_fwer']) if len(res) else float('nan')
        n_sig = int((res['p_fwer'] < 0.05).sum()) if len(res) else 0
        print(f"    [{bin_label}/{d}] {len(res)} components; largest={largest} edges; "
              f"top p_FWER={top_p:.4f}; FWER<0.05: {n_sig}")
        summary.append({'bin': bin_label, 'sign': d, 'n_autism': n_a, 'n_nt': n_n,
                         'n_components': len(res), 'largest_n_edges': largest,
                         'top_p_fwer': top_p, 'n_sig_fwer_05': n_sig})
    return summary


def main():
    df = pd.read_csv(nbs.FMRI_FILE, low_memory=False)
    df['diag'] = df['population_group'].replace('TD', 'NT')
    df = df[df['diag'].isin(['Autism', 'NT'])].copy()
    df = df[df['age_yrs'].notna() & df['Sex'].notna() & df['machine_batch'].notna()
            & df['mean_fd'].notna()].copy()
    df['age_yrs'] = df['age_yrs'].astype(float)
    df['Sex'] = df['Sex'].astype(float)
    df['machine_batch'] = df['machine_batch'].map(str)
    df['bin'] = pd.cut(df['age_yrs'], bins=BIN_EDGES, labels=BIN_LABELS, right=False)
    print(f"n={len(df)}  {df['diag'].value_counts().to_dict()}  N_PERM={N_PERM}")

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
    counts = []
    for b in BIN_LABELS:
        g = df[df['bin'] == b]
        n_a = int((g['diag'] == 'Autism').sum())
        n_n = int((g['diag'] == 'NT').sum())
        counts.append({'bin': b, 'n_autism': n_a, 'n_nt': n_n})
        if n_a < MIN_N or n_n < MIN_N:
            print(f"  [{b}] SKIPPED -- untestable (A={n_a}, NT={n_n}, need >={MIN_N} each)")
            for d in ('pos', 'neg'):
                pd.DataFrame(columns=['bin', 'sign', 'component', 'n_edges', 'n_nodes',
                                       'p_fwer', 'networks']).to_csv(
                    OUT_DIR / f"age_bin_nbs_{SAFE[b]}_{d}.csv", index=False)
                pd.DataFrame(columns=['source', 'target', 't_stat', 'p_value']).to_csv(
                    OUT_DIR / f"age_bin_nbs_edges_{SAFE[b]}_{d}.csv", index=False)
                summary.append({'bin': b, 'sign': d, 'n_autism': n_a, 'n_nt': n_n,
                                 'n_components': 0, 'largest_n_edges': np.nan,
                                 'top_p_fwer': np.nan, 'n_sig_fwer_05': 0})
            continue
        summary += run_bin(g, conn_cols, pairs, ei, ej, n_nodes, net_of, nodes, b)

    pd.DataFrame(counts).to_csv(OUT_DIR / 'age_bin_nbs_group_counts.csv', index=False)
    out = pd.DataFrame(summary)
    out.to_csv(OUT_DIR / 'age_bin_nbs_summary.csv', index=False)
    print('\n' + '=' * 70)
    print('AGE-BINNED FUNCTIONAL NBS (Autism vs NT, cohort[+sex]-adjusted)')
    print('=' * 70)
    print(out.to_string(index=False))
    print(f"\nSaved: {OUT_DIR / 'age_bin_nbs_summary.csv'}")


if __name__ == '__main__':
    raise SystemExit(main())
