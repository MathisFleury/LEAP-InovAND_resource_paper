#!/usr/bin/env python3
"""
Visual summary of the primary vs sensitivity contrasts.

Produces (in 6_functional_analysis/outputs/figures/sensitivity/):
  - sensitivity_n_significant_bar.pdf
      Bar chart of n significant connections (FDR < 0.05) per contrast,
      for primary / 6-min cutoff / covariate-adjusted.
  - sensitivity_sample_sizes_bar.pdf
      Group sample sizes per contrast for the three variants.
  - sensitivity_jaccard_heatmap.pdf
      Jaccard overlap of the significant-feature set, sensitivity vs primary.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path('/Users/mfleury/POSTDOC/LIBRAIRY/LEAP-InovAND_resource')
SUMMARY_CSV = REPO / '6_functional_analysis' / 'outputs' / 'tables' / 'sensitivity' / 'primary_vs_sensitivity_summary.csv'
OUT_DIR = REPO / '6_functional_analysis' / 'outputs' / 'figures' / 'sensitivity'
OUT_DIR.mkdir(parents=True, exist_ok=True)

VARIANT_ORDER = ['primary', '6min_cutoff', 'covariate_adjusted']
VARIANT_LABELS = {
    'primary': 'Primary\n(Welch, all subj.)',
    '6min_cutoff': '6-min cutoff\n+ FD < 1 mm',
    'covariate_adjusted': 'OLS adj.\nFD + min. quality',
}
VARIANT_COLORS = {'primary': '#4C72B0', '6min_cutoff': '#DD8452', 'covariate_adjusted': '#55A868'}


def plot_n_significant(df: pd.DataFrame, out_path: Path) -> None:
    analyses = df['analysis'].drop_duplicates().tolist()
    n = len(analyses)
    x = np.arange(n)
    width = 0.25

    fig, ax = plt.subplots(figsize=(max(8, 1.4 * n + 2), 5))
    for i, variant in enumerate(VARIANT_ORDER):
        vals = []
        for a in analyses:
            row = df[(df['analysis'] == a) & (df['variant'] == variant)]
            vals.append(int(row['n_significant_fdr05'].iloc[0]) if len(row) else 0)
        offset = (i - 1) * width
        bars = ax.bar(x + offset, vals, width, label=VARIANT_LABELS[variant],
                      color=VARIANT_COLORS[variant], edgecolor='black', linewidth=0.4)
        for rect, v in zip(bars, vals):
            ax.text(rect.get_x() + rect.get_width() / 2, rect.get_height() + max(vals) * 0.01,
                    str(v), ha='center', va='bottom', fontsize=8)

    ax.set_xticks(x)
    ax.set_xticklabels(analyses, rotation=30, ha='right')
    ax.set_ylabel('Significant connections (FDR < 0.05)')
    ax.set_title('Primary vs sensitivity: number of significant connections')
    ax.legend(loc='upper right', frameon=False, fontsize=9)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f'Saved: {out_path.name}')


def plot_sample_sizes(df: pd.DataFrame, out_path: Path) -> None:
    # For autism_vs_NT, n_group_a / n_group_b are stored separately as n_autism/n_nt in
    # the summary csv via the column-naming choices in 09_sensitivity_comparison_summary.py.
    # For cluster_*/pairwise we already have n_group_a, n_group_b.
    df = df.copy()
    # Harmonise sample-size columns: section-6 summary writes nothing for n_group_*
    # but section-7 does. For section-6 we read the standalone summary CSVs.
    def _autism_nt_sizes(analysis: str, variant: str) -> tuple[float, float]:
        if analysis != 'autism_vs_NT':
            return (np.nan, np.nan)
        if variant == 'primary':
            # The primary script writes 933-N counts inside autism_td_connectivity_analysis_summary.csv
            # but not group sizes. Pull from the merged sample (937 → 728 Autism/NT).
            return (427.0, 301.0)  # from script run earlier
        if variant == '6min_cutoff':
            p = REPO / '6_functional_analysis' / 'outputs' / 'tables' / 'sensitivity' / 'autism_vs_td_6min_summary.csv'
            if p.exists():
                s = pd.read_csv(p)
                return (float(s['n_autism'].iloc[0]), float(s['n_nt'].iloc[0]))
        if variant == 'covariate_adjusted':
            p = REPO / '6_functional_analysis' / 'outputs' / 'tables' / 'sensitivity' / 'autism_vs_td_covadj_summary.csv'
            if p.exists():
                s = pd.read_csv(p)
                return (float(s['n_autism'].iloc[0]), float(s['n_nt'].iloc[0]))
        return (np.nan, np.nan)

    for i, row in df.iterrows():
        a, v = row['analysis'], row['variant']
        if pd.isna(row.get('n_group_a')) or pd.isna(row.get('n_group_b')):
            n_a, n_b = _autism_nt_sizes(a, v)
            df.at[i, 'n_group_a'] = n_a
            df.at[i, 'n_group_b'] = n_b

    analyses = df['analysis'].drop_duplicates().tolist()
    n = len(analyses)
    x = np.arange(n)
    width = 0.25

    fig, ax = plt.subplots(figsize=(max(8, 1.4 * n + 2), 5))
    for i, variant in enumerate(VARIANT_ORDER):
        vals = []
        for a in analyses:
            row = df[(df['analysis'] == a) & (df['variant'] == variant)]
            if len(row) and not pd.isna(row['n_group_a'].iloc[0]) and not pd.isna(row['n_group_b'].iloc[0]):
                vals.append(int(row['n_group_a'].iloc[0]) + int(row['n_group_b'].iloc[0]))
            else:
                vals.append(0)
        offset = (i - 1) * width
        bars = ax.bar(x + offset, vals, width, label=VARIANT_LABELS[variant],
                      color=VARIANT_COLORS[variant], edgecolor='black', linewidth=0.4)
        for rect, v_ in zip(bars, vals):
            ax.text(rect.get_x() + rect.get_width() / 2, rect.get_height() + max(vals) * 0.01,
                    str(v_), ha='center', va='bottom', fontsize=8)

    ax.set_xticks(x)
    ax.set_xticklabels(analyses, rotation=30, ha='right')
    ax.set_ylabel('Total participants in contrast')
    ax.set_title('Primary vs sensitivity: sample sizes per contrast')
    ax.legend(loc='upper right', frameon=False, fontsize=9)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f'Saved: {out_path.name}')


def plot_jaccard(df: pd.DataFrame, out_path: Path) -> None:
    analyses = df['analysis'].drop_duplicates().tolist()
    variants_to_show = ['6min_cutoff', 'covariate_adjusted']

    M = np.full((len(analyses), len(variants_to_show)), np.nan)
    for i, a in enumerate(analyses):
        for j, v in enumerate(variants_to_show):
            row = df[(df['analysis'] == a) & (df['variant'] == v)]
            if len(row):
                val = row['jaccard_vs_primary'].iloc[0]
                if not pd.isna(val):
                    M[i, j] = float(val)

    fig, ax = plt.subplots(figsize=(6, max(4, 0.45 * len(analyses) + 2)))
    im = ax.imshow(M, cmap='viridis', vmin=0, vmax=1, aspect='auto')
    ax.set_xticks(np.arange(len(variants_to_show)))
    ax.set_xticklabels([VARIANT_LABELS[v] for v in variants_to_show])
    ax.set_yticks(np.arange(len(analyses)))
    ax.set_yticklabels(analyses)
    ax.set_title('Jaccard overlap of significant-feature set\n(sensitivity vs primary)')
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            txt = '—' if np.isnan(M[i, j]) else f'{M[i, j]:.2f}'
            ax.text(j, i, txt, ha='center', va='center',
                    color='white' if not np.isnan(M[i, j]) and M[i, j] < 0.5 else 'black',
                    fontsize=9)
    fig.colorbar(im, ax=ax, label='Jaccard')
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    print(f'Saved: {out_path.name}')


def main() -> int:
    if not SUMMARY_CSV.exists():
        raise FileNotFoundError(
            f'Comparison summary not found: {SUMMARY_CSV}. '
            'Run 09_sensitivity_comparison_summary.py first.'
        )
    df = pd.read_csv(SUMMARY_CSV)

    plot_n_significant(df, OUT_DIR / 'sensitivity_n_significant_bar.pdf')
    plot_sample_sizes(df, OUT_DIR / 'sensitivity_sample_sizes_bar.pdf')
    plot_jaccard(df, OUT_DIR / 'sensitivity_jaccard_heatmap.pdf')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
