#!/usr/bin/env python3
"""
Age trajectory of a significant-interaction thickness ROI (posterior cingulate),
autism vs NT.

Uses the SAME data as the interaction models (01_age_sex_interactions.py):
the regressed z-score anat file + curated-concat clinical, via that script's
load_anat(). NOTE: age (+ age^2) is regressed out of this file, so the
thickness values are age-detrended and the diagnosis x age interaction is ~null
here — the trajectories are expected to be near-flat. This is the honest
picture on the analysis substrate (contrast the ComBat file, where a fragile
effect appears that does not survive age>=10 restriction).

Per-group OLS fit + 95% CI on thickness adjusted for sex + cohort; age<10 shaded
(NT sparse) with a per-group age rug; interaction stat annotated per panel.
"""
import os
import sys
import importlib.util
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.formula.api as smf
import warnings
warnings.filterwarnings('ignore')

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_SECTION_DIR = os.path.dirname(_SCRIPT_DIR)
OUT = os.path.join(_SECTION_DIR, 'outputs', 'figures',
                   'posteriorcingulate_thickness_age_autism_vs_nt.pdf')

# reuse the analysis script's loaders (identical anat z-score + curated sample)
_spec = importlib.util.spec_from_file_location(
    'agesex', os.path.join(_SCRIPT_DIR, '01_age_sex_interactions.py'))
m = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m)

# lh posterior cingulate is the one thickness ROI with a significant
# diagnosis x age interaction on this substrate (p_fdr=0.025); rh is a
# non-significant contralateral contrast.
REGIONS = [('lh_posteriorcingulate_thickness', 'Left posterior cingulate (sig.)'),
           ('rh_posteriorcingulate_thickness', 'Right posterior cingulate (n.s.)')]
COLORS = {'Autism': '#d1495b', 'NT': '#2e6f95'}


def adjust(d, col):
    """Residualise on sex + cohort (add grand mean) so shown slopes match a
    sex+cohort-adjusted model."""
    sub = d[[col, 'Sex', 'cohort']].dropna()
    X = np.column_stack([np.ones(len(sub)), sub['Sex'].astype(float).values,
                         pd.get_dummies(sub['cohort'], drop_first=True).values.astype(float)])
    y = sub[col].values
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return pd.Series(y - X @ beta + y.mean(), index=sub.index)


def interaction_stat(d, col):
    s = d[[col, 'group', 'age_c', 'Sex', 'cohort']].dropna().rename(columns={col: 'y'})
    mod = smf.ols("y ~ C(group, Treatment('NT'))*age_c + C(Sex) + C(cohort)", s).fit()
    name = [p for p in mod.params.index if ':' in p][0]
    return mod.tvalues[name], mod.pvalues[name], 2 * mod.tvalues[name] / np.sqrt(mod.df_resid)


df = m.prep_common(m.load_anat(), 'population_group')
df['group'] = pd.Categorical(df['diag'].map(str))
df['age_yrs'] = df['age_yrs'].astype(float)
print(f"n={len(df)}  " + str(df['group'].value_counts().to_dict()))

fig, axes = plt.subplots(1, 2, figsize=(13, 5.2))
for ax, (col, title) in zip(axes, REGIONS):
    df['_adj'] = adjust(df, col)
    d = df[['age_yrs', 'group', '_adj']].dropna().rename(columns={'_adj': col})
    t, p, dval = interaction_stat(df, col)
    ax.axvspan(d['age_yrs'].min(), 10, color='0.5', alpha=0.10, lw=0)
    for grp in ('NT', 'Autism'):
        g = d[d['group'] == grp]
        sns.regplot(x='age_yrs', y=col, data=g, ax=ax, color=COLORS[grp], ci=95,
                    scatter_kws=dict(s=10, alpha=0.30, edgecolor='none'),
                    line_kws=dict(lw=2.2), label=f'{grp} (n={len(g)})')
        sns.rugplot(x=g['age_yrs'], ax=ax, color=COLORS[grp], height=0.04, alpha=0.5, lw=0.6)
    ax.set_title(title, fontsize=12, fontweight='bold')
    ax.set_xlabel('Age (years)   — grey: age<10 (NT sparse)')
    ax.set_ylabel('Thickness (regressed z-score, adj. sex+cohort)')
    ax.annotate(f"diagnosis × age:\n t = {t:.2f},  p = {p:.3f},  d = {dval:.2f}",
                xy=(0.97, 0.97), xycoords='axes fraction', ha='right', va='top',
                fontsize=9, bbox=dict(boxstyle='round', fc='white', ec='0.7', alpha=0.85))
    ax.legend(frameon=False, fontsize=9, loc='lower left')
    sns.despine(ax=ax)

fig.suptitle('Posterior cingulate thickness × age (autism vs NT) — significant '
             'diagnosis×age ROI, regressed z-score substrate',
             fontsize=12, fontweight='bold')
fig.tight_layout(rect=(0, 0, 1, 0.96))
fig.savefig(OUT, dpi=300, bbox_inches='tight')
print(f'Saved: {OUT}')
