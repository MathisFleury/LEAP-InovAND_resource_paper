#!/usr/bin/env python3
"""
Two age-composition sensitivity analyses for the age/sex interaction models,
reusing 01_age_sex_interactions.py's loaders + models unchanged:

  A) AGE 5-22   : keep only subjects aged 5-22 years.
  B) AGE-MATCHED: frequency-match NT and Autism within 2-year age bins (1:1 per
                  bin), so the two groups' age distributions overlap and the
                  NT:Autism ratio is constant across age.

For each subsample we re-run the anat (per-metric FDR + global-size adjustment
for area/volume) and func interaction models and compare to the full sample.

Output: ../outputs/tables/age_range_matched_summary.csv
"""
import os
import importlib.util
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_FIG = os.path.join(os.path.dirname(_SCRIPT_DIR), 'outputs', 'figures')
TABLES = os.path.join(os.path.dirname(_SCRIPT_DIR), 'outputs', 'tables')
os.makedirs(TABLES, exist_ok=True)

_spec = importlib.util.spec_from_file_location(
    'm', os.path.join(_SCRIPT_DIR, '01_age_sex_interactions.py'))
m = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(m)

SEED = 42
BIN = 2  # years


def recenter(df):
    df = df.copy()
    df['age_c'] = df['age_yrs'].astype(float) - df['age_yrs'].astype(float).mean()
    return df


def age_range(df, lo, hi):
    return recenter(df[(df['age_yrs'] >= lo) & (df['age_yrs'] <= hi)])


def age_matched(df):
    """1:1 NT/Autism within 2-year bins -> overlapping age distributions."""
    df = df.copy()
    df['_bin'] = (df['age_yrs'].astype(float) // BIN).astype(int)
    keep = []
    for _, g in df.groupby('_bin'):
        a = g[g['diag'] == 'Autism']; n = g[g['diag'] == 'NT']
        k = min(len(a), len(n))
        if k == 0:
            continue
        keep.append(a.sample(k, random_state=SEED))
        keep.append(n.sample(k, random_state=SEED))
    return recenter(pd.concat(keep).drop(columns='_bin'))


def anat_setup(df):
    """Attach the global-size covariates + return (feats, family_of, cov_of)."""
    feats = m.anat_feature_cols(df)
    area = [c for c in feats if c.endswith('_area')]
    vol = [c for c in feats if c.endswith('_volume') and c[:3] in ('lh_', 'rh_')]
    df['glob_area'] = df[area].mean(axis=1)
    df['glob_volume'] = df[vol].mean(axis=1)

    def cov(f):
        atlas, metric = m._ggseg_label(f)[1:]
        if metric == 'area':
            return 'glob_area'
        if metric == 'volume' and atlas == 'dk':
            return 'glob_volume'
        return None
    fam = lambda f: '_'.join(m._ggseg_label(f)[1:])
    return feats, fam, cov


def export_ggseg(results, outdir, n_aut, n_nt):
    """Write per-metric ggseg-input CSVs (as 04_plot_*.R expects) for a subset."""
    os.makedirs(outdir, exist_ok=True)
    pd.DataFrame([{'n_autism': n_aut, 'n_nt': n_nt}]).to_csv(
        os.path.join(outdir, 'group_counts.csv'), index=False)
    for model_name, res in results.items():
        lab = res['feature'].map(m._ggseg_label)
        d = res.assign(label=[x[0] for x in lab], atlas=[x[1] for x in lab],
                       metric=[x[2] for x in lab])
        for (atlas, metric), grp in d.groupby(['atlas', 'metric']):
            grp = grp.copy()
            grp['p_fdr'] = multipletests(grp['p_value'], method='fdr_bh')[1]
            grp[['label', 't', 'p_value', 'cohens_d', 'p_fdr']].to_csv(
                os.path.join(outdir, f't_stat_anat_{atlas}_{metric}_interaction_{model_name}.csv'),
                index=False, header=False)


def report_sig(results, tag):
    for model_name, res in results.items():
        sig = res[res['p_fdr'] < 0.05].sort_values('p_fdr')
        print(f"  {tag}/{model_name}: {len(sig)} region(s) FDR<0.05")
        if len(sig):
            print(sig[['feature', 't', 'cohens_d', 'p_value', 'p_fdr']]
                  .round(4).to_string(index=False))


def describe(df, tag):
    a = df[df['diag'] == 'Autism']['age_yrs'].astype(float)
    n = df[df['diag'] == 'NT']['age_yrs'].astype(float)
    print(f"  [{tag}] Autism n={len(a)} (med {a.median():.1f})  "
          f"NT n={len(n)} (med {n.median():.1f})  ratio A:N={len(a)/max(len(n),1):.2f}")


def main():
    anat = m.prep_common(m.load_anat(), 'population_group')
    func = m.prep_common(pd.read_csv(m.FUNC_FILE), 'population_group')

    subsets = {
        'full': lambda d: recenter(d),
        'age_5_22': lambda d: age_range(d, 5, 22),
        'age_matched': age_matched,
    }

    rows = []
    for name, fn in subsets.items():
        print(f"\n=== {name} ===")
        a_sub = fn(anat); f_sub = fn(func)
        describe(a_sub, f'anat/{name}'); describe(f_sub, f'func/{name}')

        feats, fam, cov = anat_setup(a_sub)
        a_summ, a_res = m.run_models(a_sub, feats, 'anat', family_of=fam, cov_of=cov)
        f_summ, _ = m.run_models(f_sub, m.func_feature_cols(f_sub), 'func')
        for s in a_summ + f_summ:
            s = dict(s); s['subset'] = name
            rows.append(s)

        # precision: per-region tables + surviving regions + ggseg inputs
        report_sig(a_res, f'anat/{name}')
        na = int((a_sub['diag'] == 'Autism').sum()); nn = int((a_sub['diag'] == 'NT').sum())
        for model_name, res in a_res.items():
            res.to_csv(os.path.join(TABLES, f'interaction_anat_{model_name}_{name}.csv'), index=False)
        export_ggseg(a_res, os.path.join(TABLES, f'r_input_{name}'), na, nn)

    df = pd.DataFrame(rows)[
        ['subset', 'modality', 'interaction', 'n_features',
         'n_sig_uncorrected_p05', 'n_sig_fdr_05', 'min_p_fdr', 'max_abs_cohens_d']]
    out = os.path.join(TABLES, 'age_range_matched_summary.csv')
    df.to_csv(out, index=False)
    print("\n" + "=" * 70 + "\nAGE-RANGE / AGE-MATCHED SENSITIVITY\n" + "=" * 70)
    print(df.to_string(index=False))
    print(f"\nSaved: {out}")


if __name__ == '__main__':
    main()
