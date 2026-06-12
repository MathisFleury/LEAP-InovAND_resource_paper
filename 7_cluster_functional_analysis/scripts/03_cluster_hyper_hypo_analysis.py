#!/usr/bin/env python3
"""
Cluster Functional Connectivity — Hyper/Hypo Edge Count + Effect Size

For each autism cluster vs NT, counts the number of:
- Hyperconnected edges (Cluster > NT, FDR-corrected p < 0.05, cohens_d > 0)
- Hypoconnected edges (Cluster < NT, FDR-corrected p < 0.05, cohens_d < 0)

Per Reviewer 3 (revision): also reports the distribution of Cohen's d
within each hyper / hypo bucket, so the reader can distinguish between
"many small effects" and "few large effects" — which differential cluster
sample sizes can obscure if only significance is reported.
"""

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# =============================================================================
# PATHS
# =============================================================================
_SCRIPT_DIR  = Path(__file__).parent
_SECTION_DIR = _SCRIPT_DIR.parent
OUTPUT_DIR   = _SECTION_DIR / 'outputs' / 'figures'

CLUSTERS = ['C1', 'C2', 'C3']

CLUSTER_COLORS = {
    'C1': '#393b79',
    'C2': '#e7ba52',
    'C3': '#ff9fa0',
}

HYPER_COLOR = '#D62728'   # red
HYPO_COLOR  = '#1F77B4'   # blue


def load_cluster_results():
    """Load cluster vs NT connectivity CSVs and classify hyper/hypo by Cohen's d."""
    results = {}
    for cluster in CLUSTERS:
        path = OUTPUT_DIR / f'cluster_{cluster}_connectivity_autism_vs_td.csv'
        if not path.exists():
            print(f"  WARNING: {path.name} not found, skipping {cluster}.")
            continue
        df = pd.read_csv(path)
        # Hyper / hypo classified by Cohen's d (same sign as t_stat) — keeps
        # the direction unambiguous and surfaces the effect size used downstream.
        direction = df['cohens_d'] if 'cohens_d' in df.columns else df['t_stat']
        sig = df[df['p_fdr'] < 0.05]
        sig_dir = direction[df['p_fdr'] < 0.05]
        hyper = sig[sig_dir > 0]
        hypo  = sig[sig_dir < 0]
        results[cluster] = {
            'total_sig':  len(sig),
            'n_hyper':    len(hyper),
            'n_hypo':     len(hypo),
            'n_total':    len(df),
            'mean_abs_d_hyper':   float(hyper['cohens_d'].abs().mean())  if len(hyper) and 'cohens_d' in hyper else float('nan'),
            'mean_abs_d_hypo':    float(hypo['cohens_d'].abs().mean())   if len(hypo)  and 'cohens_d' in hypo  else float('nan'),
            'median_abs_d_hyper': float(hyper['cohens_d'].abs().median()) if len(hyper) and 'cohens_d' in hyper else float('nan'),
            'median_abs_d_hypo':  float(hypo['cohens_d'].abs().median())  if len(hypo)  and 'cohens_d' in hypo  else float('nan'),
            'df':         df,
            'df_hyper':   hyper,
            'df_hypo':    hypo,
        }
        print(f"  {cluster}: {len(hyper)} hyper, {len(hypo)} hypo  "
              f"(total sig: {len(sig)} / {len(df)}) · "
              f"mean |d| hyper={results[cluster]['mean_abs_d_hyper']:.2f} · "
              f"hypo={results[cluster]['mean_abs_d_hypo']:.2f}")
    return results


def plot_edge_counts(results):
    """Grouped bar chart: hyper and hypo edge counts per cluster."""
    clusters = list(results.keys())
    n_hyper  = [results[c]['n_hyper'] for c in clusters]
    n_hypo   = [results[c]['n_hypo']  for c in clusters]

    x = np.arange(len(clusters))
    width = 0.35

    fig, ax = plt.subplots(figsize=(7, 5))

    bars_hyper = ax.bar(x - width/2, n_hyper, width,
                        color=HYPER_COLOR, label='Hyperconnectivity (Cluster > NT)',
                        alpha=0.85)
    bars_hypo  = ax.bar(x + width/2, n_hypo,  width,
                        color=HYPO_COLOR,  label='Hypoconnectivity (Cluster < NT)',
                        alpha=0.85)

    # Value labels on bars
    for bar in bars_hyper:
        h = bar.get_height()
        if h > 0:
            ax.text(bar.get_x() + bar.get_width()/2, h + 0.5, str(int(h)),
                    ha='center', va='bottom', fontsize=11, fontweight='bold')
    for bar in bars_hypo:
        h = bar.get_height()
        if h > 0:
            ax.text(bar.get_x() + bar.get_width()/2, h + 0.5, str(int(h)),
                    ha='center', va='bottom', fontsize=11, fontweight='bold')

    ax.set_xticks(x)
    ax.set_xticklabels([f'Cluster {c}' for c in clusters], fontsize=13)
    ax.set_ylabel('Number of significant edges (FDR < 0.05)', fontsize=12)
    ax.set_title('Functional Connectivity — Clusters vs NT\nHyper / Hypoconnectivity', fontsize=13)
    ax.legend(fontsize=10)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    out = OUTPUT_DIR / 'cluster_hyper_hypo_edge_counts.pdf'
    fig.savefig(out, bbox_inches='tight', dpi=300)
    plt.close(fig)
    print(f"Saved: {out.name}")


def plot_hypo_hyper_comparison(results):
    """Stacked horizontal bar chart: percentage hyper vs hypo per cluster."""
    clusters = [c for c in CLUSTERS if c in results]
    totals   = [results[c]['n_hyper'] + results[c]['n_hypo'] for c in clusters]

    fig, ax = plt.subplots(figsize=(6, 4))

    y = np.arange(len(clusters))
    for i, cluster in enumerate(clusters):
        total = totals[i]
        if total == 0:
            continue
        n_h = results[cluster]['n_hyper']
        n_y = results[cluster]['n_hypo']
        ax.barh(i, n_h, color=HYPER_COLOR, alpha=0.85)
        ax.barh(i, n_y, left=n_h, color=HYPO_COLOR, alpha=0.85)
        ax.text(n_h / 2, i, str(n_h), ha='center', va='center',
                fontsize=10, color='white', fontweight='bold')
        ax.text(n_h + n_y / 2, i, str(n_y), ha='center', va='center',
                fontsize=10, color='white', fontweight='bold')

    ax.set_yticks(y)
    ax.set_yticklabels([f'Cluster {c}' for c in clusters], fontsize=12)
    ax.set_xlabel('Number of significant edges (FDR < 0.05)', fontsize=11)
    ax.set_title('Hyper vs Hypo Connectivity per Cluster vs NT', fontsize=12)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    patches = [mpatches.Patch(color=HYPER_COLOR, label='Hyper (Cluster > NT)'),
               mpatches.Patch(color=HYPO_COLOR,  label='Hypo (Cluster < NT)')]
    ax.legend(handles=patches, fontsize=9)

    out = OUTPUT_DIR / 'cluster_hyper_hypo_stacked.pdf'
    fig.savefig(out, bbox_inches='tight', dpi=300)
    plt.close(fig)
    print(f"Saved: {out.name}")


def save_summary_csv(results):
    """Save summary CSV with edge counts and Cohen's d summaries per cluster."""
    rows = []
    for cluster, r in results.items():
        rows.append({
            'cluster':             cluster,
            'n_total_edges':       r['n_total'],
            'n_significant':       r['total_sig'],
            'n_hyper':             r['n_hyper'],
            'n_hypo':              r['n_hypo'],
            'pct_hyper':           round(100 * r['n_hyper'] / r['total_sig'], 1) if r['total_sig'] > 0 else 0,
            'pct_hypo':            round(100 * r['n_hypo']  / r['total_sig'], 1) if r['total_sig'] > 0 else 0,
            'mean_abs_d_hyper':    r.get('mean_abs_d_hyper'),
            'mean_abs_d_hypo':     r.get('mean_abs_d_hypo'),
            'median_abs_d_hyper':  r.get('median_abs_d_hyper'),
            'median_abs_d_hypo':   r.get('median_abs_d_hypo'),
        })
    df = pd.DataFrame(rows)
    out = OUTPUT_DIR / 'cluster_hyper_hypo_summary.csv'
    df.to_csv(out, index=False)
    print(f"Saved: {out.name}")
    print(df.to_string(index=False))


def plot_cohens_d_distribution(results):
    """Violin / strip distribution of Cohen's d per cluster × direction.

    Reviewer 3: shows the magnitude of effect sizes alongside counts so that
    'many small effects' and 'few large effects' can be distinguished —
    important when cluster sample sizes differ.
    """
    rows = []
    for cluster, r in results.items():
        if 'cohens_d' in r['df_hyper'].columns:
            for v in r['df_hyper']['cohens_d'].dropna():
                rows.append({'cluster': cluster, 'direction': 'Hyper', 'd': float(v)})
        if 'cohens_d' in r['df_hypo'].columns:
            for v in r['df_hypo']['cohens_d'].dropna():
                rows.append({'cluster': cluster, 'direction': 'Hypo', 'd': float(v)})
    if not rows:
        print("  No Cohen's d to plot.")
        return
    import seaborn as sns
    df = pd.DataFrame(rows)
    palette = {'Hyper': HYPER_COLOR, 'Hypo': HYPO_COLOR}

    fig, ax = plt.subplots(figsize=(7, 4.5))
    sns.violinplot(data=df, x='cluster', y='d', hue='direction',
                   order=[c for c in CLUSTERS if c in results],
                   hue_order=['Hyper', 'Hypo'],
                   palette=palette, inner=None, linewidth=0,
                   split=True, cut=0, ax=ax)
    sns.stripplot(data=df, x='cluster', y='d', hue='direction',
                  order=[c for c in CLUSTERS if c in results],
                  hue_order=['Hyper', 'Hypo'],
                  palette={'Hyper': 'black', 'Hypo': 'black'},
                  size=1.6, alpha=0.45, dodge=True, jitter=0.25, ax=ax)
    ax.axhline(0, color='black', linewidth=0.4)
    ax.set_xlabel('')
    ax.set_ylabel("Cohen's d (cluster vs NT)", fontsize=11)
    ax.set_title('Effect-size distribution of FDR-significant edges per cluster',
                 fontsize=12)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    # Trim duplicate seaborn legend (one from violin + one from strip)
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend([by_label.get('Hyper'), by_label.get('Hypo')],
              ['Hyper (Cluster > NT)', 'Hypo (Cluster < NT)'],
              frameon=False, loc='upper right', fontsize=9)
    plt.tight_layout()
    out = OUTPUT_DIR / 'cluster_hyper_hypo_cohens_d_distribution.pdf'
    fig.savefig(out, bbox_inches='tight', dpi=300)
    plt.close(fig)
    print(f"Saved: {out.name}")


def main():
    print("=" * 60)
    print("CLUSTER HYPER/HYPO CONNECTIVITY ANALYSIS")
    print("=" * 60)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    results = load_cluster_results()
    if not results:
        print("No data found. Run 01_generate_cluster_fmri_inputs.py first.")
        return 1

    print()
    save_summary_csv(results)
    plot_edge_counts(results)
    plot_hypo_hyper_comparison(results)
    plot_cohens_d_distribution(results)

    print(f"\nAll outputs saved to: {OUTPUT_DIR}")
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())
