#!/usr/bin/env python3
"""
Subcortical Volume Brain Visualization — Clusters vs NT (REVISION, R3)

v2 of ../../scripts/03_plot_subcortical_cluster_yabplot.py, re-rendered on
the QC+ComBat+regression FreeSurfer dataset (LEAP_W1 / INOVAND_T1 priority).
Reads ../outputs/tables/r_input_files/ produced by
1_generate_cluster_mri_inputs.py.
"""

from pathlib import Path
import pandas as pd
import pyvista as pv

pv.global_theme.notebook = True

from yabplot import plot_subcortical, get_atlas_regions

# =============================================================================
# PATHS
# =============================================================================
_SCRIPT_DIR  = Path(__file__).parent
_SECTION_DIR = _SCRIPT_DIR.parent
INPUT_DIR    = _SECTION_DIR / 'outputs' / 'tables' / 'r_input_files'
OUTPUT_DIR   = _SECTION_DIR / 'outputs' / 'figures'

CLUSTERS = ['C1', 'C2', 'C3']

LABEL_MAP = {
    'Left-Thalamus':    'Left_Thalamus',
    'Right-Thalamus':   'Right_Thalamus',
    'Left-Caudate':     'Left_Caudate',
    'Right-Caudate':    'Right_Caudate',
    'Left-Putamen':     'Left_Putamen',
    'Right-Putamen':    'Right_Putamen',
    'Left-Pallidum':    'Left_Pallidum',
    'Right-Pallidum':   'Right_Pallidum',
    'Left-Hippocampus': 'Left_Hippocampus',
    'Right-Hippocampus':'Right_Hippocampus',
    'Left-Amygdala':    'Left_Amygdala',
    'Right-Amygdala':   'Right_Amygdala',
}

ALL_REGIONS = get_atlas_regions('aseg', 'subcortical')


def load_data(cluster):
    path = INPUT_DIR / f't_stat_cluster_{cluster}_volume_aseg_mri_cluster_vs_td.csv'
    if not path.exists():
        print(f"  WARNING: {path.name} not found, skipping.")
        return None
    df = pd.read_csv(path)
    return df


def make_yabplot_data(df, fdr_only=False, value_col='cohens_d'):
    data = {}
    for _, row in df.iterrows():
        yab_label = LABEL_MAP.get(row['label'])
        if yab_label:
            if fdr_only and row.get('p_fdr', 1.0) >= 0.05:
                data[yab_label] = float('nan')
            else:
                data[yab_label] = float(row[value_col])
    # Fill missing regions with NaN
    for r in ALL_REGIONS:
        if r not in data:
            data[r] = float('nan')
    return data


def plot_and_save(data, vmax, out_path, title='', cbar_label="Cohen's d"):
    plotter = plot_subcortical(
        data=data,
        atlas='aseg',
        cmap='RdBu_r',
        vminmax=[-vmax, vmax],
        display_type='object',
        figsize=(1200, 500),
    )
    if title:
        plotter.add_title(title, font_size=12)
    plotter.screenshot(str(out_path))
    plotter.close()
    print(f"  Saved: {out_path.name}")


def main():
    print("=" * 60)
    print("SUBCORTICAL VOLUME YABPLOT — Clusters vs NT  (Cohen's d, R3 revision)")
    print("=" * 60)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Determine global Cohen's d range across all clusters for consistent scale
    all_d = []
    dfs = {}
    for cluster in CLUSTERS:
        df = load_data(cluster)
        if df is not None:
            dfs[cluster] = df
            for _, row in df.iterrows():
                if LABEL_MAP.get(row['label']) and 'cohens_d' in row:
                    val = row['cohens_d']
                    if pd.notna(val):
                        all_d.append(abs(float(val)))

    if not all_d:
        print("ERROR: No data found for any cluster.")
        return 1

    d_max = max(all_d)
    print(f"Global |Cohen's d| max: {d_max:.3f}")

    for cluster, df in dfs.items():
        print(f"\nCluster {cluster}:")
        data_all = make_yabplot_data(df, fdr_only=False, value_col='cohens_d')
        plot_and_save(
            data_all, d_max,
            OUTPUT_DIR / f'cluster_{cluster}_volume_aseg_vs_nt_yabplot.png',
            title=f'Cluster {cluster} vs NT — Subcortical Volume (Cohen\'s d)',
        )
        data_fdr = make_yabplot_data(df, fdr_only=True, value_col='cohens_d')
        plot_and_save(
            data_fdr, d_max,
            OUTPUT_DIR / f'cluster_{cluster}_volume_aseg_vs_nt_yabplot_fdr.png',
            title=f'Cluster {cluster} vs NT — Subcortical Volume (Cohen\'s d, FDR<0.05)',
        )

    print(f"\nAll outputs saved to: {OUTPUT_DIR}")
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())
