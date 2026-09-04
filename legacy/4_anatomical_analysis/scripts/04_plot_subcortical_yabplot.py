#!/usr/bin/env python3
"""
Subcortical Volume Brain Visualization — Autism vs NT
Uses yabplot for 3D subcortical structure plots.
"""

from pathlib import Path
import pandas as pd
import pyvista as pv

# Offscreen rendering required outside Jupyter
pv.global_theme.notebook = True

from yabplot import plot_subcortical, get_atlas_regions

# =============================================================================
# PATHS
# =============================================================================
_SCRIPT_DIR  = Path(__file__).parent
_SECTION_DIR = _SCRIPT_DIR.parent
INPUT_DIR    = _SECTION_DIR / 'outputs' / 'figures' / 'r_input_files'
OUTPUT_DIR   = _SECTION_DIR / 'outputs' / 'figures'

INPUT_FILE = INPUT_DIR / 't_stat_anat_aseg_volume_mri_autism_vs_control.csv'

# Mapping: our CSV label  →  yabplot aseg region name
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


def load_data(filepath):
    df = pd.read_csv(filepath, header=None,
                     names=['label', 't_stat', 'p_val', 'cohens_d', 'p_fdr'])
    return df


def make_yabplot_data(df):
    """Map CSV labels to yabplot region names; return dict {region: t_stat}."""
    data = {}
    for _, row in df.iterrows():
        yab_label = LABEL_MAP.get(row['label'])
        if yab_label:
            data[yab_label] = float(row['t_stat'])
    return data


def plot_and_save(data, t_stat_max, out_path, title=''):
    plotter = plot_subcortical(
        data=data,
        atlas='aseg',
        cmap='RdBu_r',
        vminmax=[-t_stat_max, t_stat_max],
        display_type='object',
        figsize=(1200, 500),
    )
    if title:
        plotter.add_title(title, font_size=12)
    plotter.screenshot(str(out_path))
    plotter.close()
    print(f"Saved: {out_path.name}")


def main():
    print("=" * 60)
    print("SUBCORTICAL VOLUME YABPLOT — Autism vs NT")
    print("=" * 60)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    if not INPUT_FILE.exists():
        print(f"ERROR: Input file not found: {INPUT_FILE}")
        return 1

    df = load_data(INPUT_FILE)
    data = make_yabplot_data(df)

    if not data:
        print("ERROR: No matching regions found.")
        return 1

    t_stat_max = max(abs(v) for v in data.values())
    print(f"Regions mapped: {len(data)}, t-stat range: ±{t_stat_max:.2f}")

    # All regions (no FDR filter)
    plot_and_save(data, t_stat_max,
                  OUTPUT_DIR / 'subcortical_volume_autism_vs_nt_yabplot.png',
                  title='Autism vs NT — Subcortical Volume (t-statistic)')

    # FDR-significant regions only (others set to NaN → shown in grey)
    sig_labels = set(df[df['p_fdr'] < 0.05]['label'])
    data_fdr = {k: v for k, v in data.items()
                if any(LABEL_MAP.get(lbl) == k for lbl in sig_labels)}
    # Pad with NaN for non-significant regions
    all_regions = get_atlas_regions('aseg', 'subcortical')
    data_fdr_full = {r: data_fdr.get(r, float('nan')) for r in all_regions}
    plot_and_save(data_fdr_full, t_stat_max,
                  OUTPUT_DIR / 'subcortical_volume_autism_vs_nt_yabplot_fdr.png',
                  title='Autism vs NT — Subcortical Volume (FDR p<0.05)')

    print(f"\nAll outputs saved to: {OUTPUT_DIR}")
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(main())
