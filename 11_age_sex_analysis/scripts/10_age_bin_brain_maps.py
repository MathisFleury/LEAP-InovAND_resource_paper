#!/usr/bin/env python3
"""
Autism-vs-NT anatomical brain maps computed separately within each age bin --
extends 09_age_bin_regional_effect.py's single-region check to all DK cortical
(thickness/area/volume) + ASEG subcortical regions, rendered with ggseg.

Same bins as 09 (<6, 6-9, 10-13, 14-17, 18-24, 25+). <6 and 6-9 are kept
SEPARATE, not merged: merging them made "<10" compare mostly-young autism
cases (mean age 5.8y, spread 2-10y) against mostly-older NT cases (mean 8.0y,
concentrated 6-10y, since NT<6 is rare, n=6) -- a ~2.2-year within-bin age
gap between groups that produced a broad, thickness-heavy false-positive
signal (40/232 features) driven by normal developmental cortical thinning,
not diagnosis (2026-08-31 correction). No bin is suppressed: a bin's
group_counts.csv carries its N, so the R figure's column header shows it
automatically (e.g. "NT n=6" for <6) -- read it before trusting a sparse-bin
column; <6 (NT=6) is expected to be uninterpretable by inspection, not
hidden.

Per bin: Welch's two-sample t-test (autism vs NT) per feature, identical to
the primary whole-cohort test in
4_anatomical_analysis/scripts/1_anatomical_mri_autism_nt.py
(ttest_ind(a, b, equal_var=False); pooled-SD Cohen's d) -- per user decision,
2026-08-30: no covariates in the test itself, for consistency with the
primary analysis's test type. FDR-BH within each metric family (thickness_dk
/ area_dk / volume_dk / volume_aseg), matching the primary analysis's
correction convention too.
Uses load_anat_noage() (not load_anat()): the standard z-score table has age
regressed out globally, which would make an age-bin check close to
meaningless; this input is ComBat + regressed for sex/eTIV/euler only, same
as the primary group analysis, with age preserved as the stratifying
variable -- so sex/eTIV/euler are already adjusted for upstream, same as the
primary test's z-score table.
Known, accepted trade-off (per user decision): a bin's cohort composition can
be imbalanced between groups -- e.g. the 18-24 bin is 95% LEAP for Autism vs
65% LEAP for NT -- so a raw t-test there can pick up cohort along with
diagnosis. An earlier version of this script added +cohort [+sex] as
covariates in an OLS model to guard against this; that was dropped in favour
of matching the primary analysis's test exactly. Writes ggseg-input CSVs,
then runs 11_age_bin_brain_grid.R, which composes every (metric x bin) panel
into ONE figure -- rows = metric, columns = age bin, one shared Cohen's d
colour scale.

Output: outputs/tables/r_input_bin_<bin>/*.csv
        outputs/figures/composite/age_bin_groupdiff_all_features.pdf
"""
import os
import subprocess
import importlib.util
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings('ignore')

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_SECTION_DIR = os.path.dirname(_SCRIPT_DIR)
TABLES_DIR = os.path.join(_SECTION_DIR, 'outputs', 'tables')  # pure-table step: r_input_bin_* only
R_GRID_SCRIPT = os.path.join(_SCRIPT_DIR, '11_age_bin_brain_grid.R')

_spec = importlib.util.spec_from_file_location(
    'agesex', os.path.join(_SCRIPT_DIR, '01_age_sex_interactions.py'))
m = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m)

BIN_EDGES = [0, 6, 10, 14, 18, 25, 200]
BIN_LABELS = ['<6', '6-9', '10-13', '14-17', '18-24', '25+']
SAFE = {'<6': 'lt6', '6-9': '6_9', '10-13': '10_13', '14-17': '14_17',
        '18-24': '18_24', '25+': '25plus'}


def _cohens_d(a, b):
    """Pooled-SD Cohen's d -- identical formula to
    1_anatomical_mri_autism_nt.py's _cohens_d."""
    n_a, n_b = len(a), len(b)
    sd = np.sqrt(((n_a - 1) * np.var(a, ddof=1) + (n_b - 1) * np.var(b, ddof=1))
                 / (n_a + n_b - 2))
    if sd == 0:
        return np.nan
    return (np.mean(a) - np.mean(b)) / sd


def _fit_diag(sub):
    """Welch's two-sample t-test (autism vs NT), no covariates -- identical
    test to the primary whole-cohort analysis. Returns (t, p, cohens_d)."""
    a = sub.loc[sub['diag'] == 'Autism', 'y'].values
    b = sub.loc[sub['diag'] == 'NT', 'y'].values
    if len(a) < 2 or len(b) < 2:
        return None
    t, p = ttest_ind(a, b, equal_var=False)
    return t, p, _cohens_d(a, b)


def _selftest():
    rng = np.random.default_rng(0)
    n = 100
    sub = pd.DataFrame({
        'y': np.concatenate([rng.normal(1, 1, n), rng.normal(0, 1, n)]),
        'diag': pd.Categorical(['Autism'] * n + ['NT'] * n),
    })
    t, p, d = _fit_diag(sub)
    assert 0.6 < d < 1.4


def run_bin(g, feats):
    """Welch's t-test per feature, FDR-BH within (atlas, metric)."""
    rows = []
    for f in feats:
        sub = g[[f, 'diag']].dropna().rename(columns={f: 'y'})
        if sub['diag'].nunique() < 2 or len(sub) < 10:
            continue
        res_f = _fit_diag(sub)
        if res_f is None:
            continue
        t, p, d = res_f
        rows.append({'feature': f, 't': t, 'p_value': p, 'cohens_d': d})
    res = pd.DataFrame(rows)
    if res.empty:
        return res
    lab = res['feature'].map(m._ggseg_label)
    res['label'] = [x[0] for x in lab]
    res['atlas'] = [x[1] for x in lab]
    res['metric'] = [x[2] for x in lab]
    res['p_fdr'] = np.nan
    for _, idx in res.groupby(['atlas', 'metric']).groups.items():
        res.loc[idx, 'p_fdr'] = multipletests(res.loc[idx, 'p_value'], method='fdr_bh')[1]
    return res


def export(res, n_aut, n_nt, outdir):
    os.makedirs(outdir, exist_ok=True)
    pd.DataFrame([{'n_autism': n_aut, 'n_nt': n_nt}]).to_csv(
        os.path.join(outdir, 'group_counts.csv'), index=False)
    for (atlas, metric), grp in res.groupby(['atlas', 'metric']):
        grp[['label', 't', 'p_value', 'cohens_d', 'p_fdr']].to_csv(
            os.path.join(outdir, f't_stat_anat_{atlas}_{metric}_interaction_groupdiff.csv'),
            index=False, header=False)


def main():
    _selftest()
    df = m.prep_common(m.load_anat_noage(), 'population_group')
    df['age_yrs'] = df['age_yrs'].astype(float)
    feats = m.anat_feature_cols(df)
    df['bin'] = pd.cut(df['age_yrs'], bins=BIN_EDGES, labels=BIN_LABELS, right=False)

    for b in BIN_LABELS:
        g = df[df['bin'] == b]
        n_aut, n_nt = int((g['diag'] == 'Autism').sum()), int((g['diag'] == 'NT').sum())
        print(f"\n[{b}] Autism={n_aut}  NT={n_nt}")
        res = run_bin(g, feats)
        outdir = os.path.join(TABLES_DIR, f'r_input_bin_{SAFE[b]}')
        export(res, n_aut, n_nt, outdir)

    print("\nComposing per-metric grids (age bins as columns)...")
    r = subprocess.run(['Rscript', R_GRID_SCRIPT], capture_output=True, text=True)
    print(r.stdout.strip())
    if r.returncode != 0:
        print(f"  R FAILED:\n{r.stderr[-2000:]}")
    print("\nDone.")


if __name__ == '__main__':
    main()
