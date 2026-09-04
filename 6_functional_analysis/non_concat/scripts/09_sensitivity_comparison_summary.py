#!/usr/bin/env python3
"""
Compare primary vs sensitivity analyses (6-min cutoff & covariate-adjusted)
for the Autism-vs-NT and the cluster-stratified fMRI connectivity contrasts.

For each contrast we report:
  - n_subjects per group used by each model
  - n significant connections at FDR < 0.05
  - hyper / hypo split
  - Jaccard overlap of the significant feature set against the primary set

Writes:
  6_functional_analysis/non_concat/outputs/tables/sensitivity/primary_vs_sensitivity_summary.csv
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

REPO_ROOT = Path('/Users/mfleury/POSTDOC/LIBRAIRY/LEAP-InovAND_resource')
SEC6 = REPO_ROOT / '6_functional_analysis' / 'non_concat'
SEC7 = REPO_ROOT / '7_cluster_functional_analysis'
OUT_DIR = SEC6 / 'outputs' / 'tables' / 'sensitivity'
OUT_DIR.mkdir(parents=True, exist_ok=True)


def _read_csv_if_exists(p: Path) -> pd.DataFrame | None:
    return pd.read_csv(p) if p.is_file() else None


def _sig_set(df: pd.DataFrame, alpha: float = 0.05,
             feature_col: str = 'feature', p_col: str = 'p_value_fdr') -> set[str]:
    """Return the set of feature names with FDR < alpha."""
    if df is None or len(df) == 0:
        return set()
    if p_col not in df.columns:
        for alt in ('p_fdr', 'p_value_fdr'):
            if alt in df.columns:
                p_col = alt
                break
    if feature_col not in df.columns:
        for alt in ('feature', 'region', 'source_target'):
            if alt in df.columns:
                feature_col = alt
                break
    if p_col not in df.columns or feature_col not in df.columns:
        return set()
    return set(df.loc[df[p_col] < alpha, feature_col].astype(str).tolist())


def _jaccard(a: set, b: set) -> float:
    if not a and not b:
        return float('nan')
    return len(a & b) / len(a | b)


# -----------------------------------------------------------------------------
# Section 6 — Autism vs NT
# -----------------------------------------------------------------------------

def autism_vs_nt_rows() -> list[dict]:
    """Build rows comparing primary, 6-min, and covariate-adjusted Autism-vs-NT."""
    rows = []

    # ---- Primary ----
    # The primary script only saves split full-data CSVs (no combined all_results).
    p_hyper = _read_csv_if_exists(SEC6 / 'outputs' / 'figures' / 'autism_vs_td_hyperconnectivity_full_data.csv')
    p_hypo = _read_csv_if_exists(SEC6 / 'outputs' / 'figures' / 'autism_vs_td_hypoconnectivity_full_data.csv')

    def _feat_set_split(df: pd.DataFrame | None) -> set[str]:
        if df is None or len(df) == 0:
            return set()
        # reconstruct 'feature' as con_<source>/<target>
        return set('con_' + df['source'].astype(str) + '/' + df['target'].astype(str))

    primary_hyper = _feat_set_split(p_hyper)
    primary_hypo = _feat_set_split(p_hypo)
    primary_sig = primary_hyper | primary_hypo

    rows.append({
        'analysis': 'autism_vs_NT',
        'variant': 'primary',
        'model': 'Welch_t_test',
        'min_minutes': None,
        'max_mean_fd': None,
        'covariates': '',
        'n_significant_fdr05': len(primary_sig),
        'n_hyperconnectivity': len(primary_hyper),
        'n_hypoconnectivity': len(primary_hypo),
        'jaccard_vs_primary': 1.0,
    })

    # ---- 6-min sensitivity ----
    s6 = _read_csv_if_exists(SEC6 / 'outputs' / 'tables' / 'sensitivity' / 'autism_vs_td_6min_all_results.csv')
    sum6 = _read_csv_if_exists(SEC6 / 'outputs' / 'tables' / 'sensitivity' / 'autism_vs_td_6min_summary.csv')
    sig6 = _sig_set(s6)
    if s6 is not None:
        hyper6 = set(s6.loc[(s6['p_value_fdr'] < 0.05) & (s6['t_stat'] > 0), 'feature'].astype(str))
        hypo6 = set(s6.loc[(s6['p_value_fdr'] < 0.05) & (s6['t_stat'] < 0), 'feature'].astype(str))
    else:
        hyper6, hypo6 = set(), set()
    rows.append({
        'analysis': 'autism_vs_NT',
        'variant': '6min_cutoff',
        'model': 'Welch_t_test',
        'min_minutes': float(sum6['min_minutes_quality_data'].iloc[0]) if sum6 is not None else None,
        'max_mean_fd': float(sum6['max_mean_fd'].iloc[0]) if sum6 is not None else None,
        'covariates': '',
        'n_significant_fdr05': len(sig6),
        'n_hyperconnectivity': len(hyper6),
        'n_hypoconnectivity': len(hypo6),
        'jaccard_vs_primary': _jaccard(primary_sig, sig6),
    })

    # ---- Covariate-adjusted ----
    sc = _read_csv_if_exists(SEC6 / 'outputs' / 'tables' / 'sensitivity' / 'autism_vs_td_covadj_all_results.csv')
    sumc = _read_csv_if_exists(SEC6 / 'outputs' / 'tables' / 'sensitivity' / 'autism_vs_td_covadj_summary.csv')
    sigc = _sig_set(sc)
    if sc is not None:
        hyperc = set(sc.loc[(sc['p_value_fdr'] < 0.05) & (sc['t_stat'] > 0), 'feature'].astype(str))
        hypoc = set(sc.loc[(sc['p_value_fdr'] < 0.05) & (sc['t_stat'] < 0), 'feature'].astype(str))
    else:
        hyperc, hypoc = set(), set()
    rows.append({
        'analysis': 'autism_vs_NT',
        'variant': 'covariate_adjusted',
        'model': 'OLS_group+mean_fd+minutes_quality',
        'min_minutes': None,
        'max_mean_fd': float(sumc['max_mean_fd'].iloc[0]) if sumc is not None else None,
        'covariates': str(sumc['covariates'].iloc[0]) if sumc is not None else '',
        'n_significant_fdr05': len(sigc),
        'n_hyperconnectivity': len(hyperc),
        'n_hypoconnectivity': len(hypoc),
        'jaccard_vs_primary': _jaccard(primary_sig, sigc),
    })

    return rows


# -----------------------------------------------------------------------------
# Section 7 — cluster comparisons
# -----------------------------------------------------------------------------

def cluster_rows() -> list[dict]:
    rows = []

    # ---- Primary cluster outputs ----
    sec7_fig = SEC7 / 'outputs' / 'figures'
    sec7_sens = SEC7 / 'outputs' / 'figures' / 'sensitivity'

    # Cluster vs NT
    for cluster in ['C1', 'C2', 'C3']:
        primary_path = sec7_fig / f'cluster_{cluster}_connectivity_autism_vs_td.csv'
        s6_path = sec7_sens / f'cluster_{cluster}_connectivity_autism_vs_td_6min.csv'
        sc_path = sec7_sens / f'cluster_{cluster}_connectivity_autism_vs_td_covadj.csv'

        primary = _read_csv_if_exists(primary_path)
        s6 = _read_csv_if_exists(s6_path)
        sc = _read_csv_if_exists(sc_path)

        primary_sig = _sig_set(primary, feature_col='region', p_col='p_fdr')
        s6_sig = _sig_set(s6, feature_col='region', p_col='p_fdr')
        sc_sig = _sig_set(sc, feature_col='feature', p_col='p_fdr')

        if primary is not None:
            n_aut = int(primary['autism_n'].iloc[0]) if 'autism_n' in primary.columns and len(primary) else None
            n_nt = int(primary['td_n'].iloc[0]) if 'td_n' in primary.columns and len(primary) else None
        else:
            n_aut = n_nt = None
        rows.append({
            'analysis': f'cluster_{cluster}_vs_NT',
            'variant': 'primary',
            'model': 'Welch_t_test',
            'n_group_a': n_aut, 'n_group_b': n_nt,
            'n_significant_fdr05': len(primary_sig),
            'jaccard_vs_primary': 1.0,
        })
        if s6 is not None:
            n_aut6 = int(s6['autism_n'].iloc[0]) if 'autism_n' in s6.columns and len(s6) else None
            n_nt6 = int(s6['td_n'].iloc[0]) if 'td_n' in s6.columns and len(s6) else None
            rows.append({
                'analysis': f'cluster_{cluster}_vs_NT',
                'variant': '6min_cutoff',
                'model': 'Welch_t_test',
                'n_group_a': n_aut6, 'n_group_b': n_nt6,
                'n_significant_fdr05': len(s6_sig),
                'jaccard_vs_primary': _jaccard(primary_sig, s6_sig),
            })
        if sc is not None:
            n_autc = int(sc['autism_n'].iloc[0]) if 'autism_n' in sc.columns and len(sc) else None
            n_ntc = int(sc['td_n'].iloc[0]) if 'td_n' in sc.columns and len(sc) else None
            rows.append({
                'analysis': f'cluster_{cluster}_vs_NT',
                'variant': 'covariate_adjusted',
                'model': 'OLS_group+mean_fd+minutes_quality',
                'n_group_a': n_autc, 'n_group_b': n_ntc,
                'n_significant_fdr05': len(sc_sig),
                'jaccard_vs_primary': _jaccard(primary_sig, sc_sig),
            })

    # Pairwise
    for c1, c2 in [('C1', 'C2'), ('C1', 'C3'), ('C2', 'C3')]:
        primary_path = sec7_fig / f'pairwise_{c1}_vs_{c2}_connectivity.csv'
        s6_path = sec7_sens / f'pairwise_{c1}_vs_{c2}_connectivity_6min.csv'
        sc_path = sec7_sens / f'pairwise_{c1}_vs_{c2}_connectivity_covadj.csv'

        primary = _read_csv_if_exists(primary_path)
        s6 = _read_csv_if_exists(s6_path)
        sc = _read_csv_if_exists(sc_path)

        primary_sig = _sig_set(primary, feature_col='region', p_col='p_fdr')
        s6_sig = _sig_set(s6, feature_col='region', p_col='p_fdr')
        sc_sig = _sig_set(sc, feature_col='feature', p_col='p_fdr')

        rows.append({
            'analysis': f'pairwise_{c1}_vs_{c2}',
            'variant': 'primary',
            'model': 'Welch_t_test',
            'n_group_a': int(primary[f'{c1}_n'].iloc[0]) if primary is not None and f'{c1}_n' in primary.columns else None,
            'n_group_b': int(primary[f'{c2}_n'].iloc[0]) if primary is not None and f'{c2}_n' in primary.columns else None,
            'n_significant_fdr05': len(primary_sig),
            'jaccard_vs_primary': 1.0,
        })
        if s6 is not None:
            rows.append({
                'analysis': f'pairwise_{c1}_vs_{c2}',
                'variant': '6min_cutoff',
                'model': 'Welch_t_test',
                'n_group_a': int(s6[f'{c1}_n'].iloc[0]) if f'{c1}_n' in s6.columns else None,
                'n_group_b': int(s6[f'{c2}_n'].iloc[0]) if f'{c2}_n' in s6.columns else None,
                'n_significant_fdr05': len(s6_sig),
                'jaccard_vs_primary': _jaccard(primary_sig, s6_sig),
            })
        if sc is not None:
            rows.append({
                'analysis': f'pairwise_{c1}_vs_{c2}',
                'variant': 'covariate_adjusted',
                'model': 'OLS_group+mean_fd+minutes_quality',
                'n_group_a': int(sc[f'n_{c1}'].iloc[0]) if f'n_{c1}' in sc.columns else None,
                'n_group_b': int(sc[f'n_{c2}'].iloc[0]) if f'n_{c2}' in sc.columns else None,
                'n_significant_fdr05': len(sc_sig),
                'jaccard_vs_primary': _jaccard(primary_sig, sc_sig),
            })

    return rows


def main() -> int:
    sec6_rows = autism_vs_nt_rows()
    sec7_rows = cluster_rows()

    df = pd.DataFrame(sec6_rows + sec7_rows)
    out_path = OUT_DIR / 'primary_vs_sensitivity_summary.csv'
    df.to_csv(out_path, index=False)
    print(f'Saved comparison summary: {out_path}\n')
    print(df.to_string(index=False))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
