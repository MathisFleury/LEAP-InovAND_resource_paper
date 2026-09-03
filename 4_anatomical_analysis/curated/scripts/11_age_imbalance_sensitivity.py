#!/usr/bin/env python3
"""
Age-imbalance sensitivity for the Autism-vs-NT anatomical analysis.

Run on the ComBat file (retains age + un-regressed region values) — NOT the
regressed z-score file, whose pooled age+age^2 regression would erase the very
age signal this check needs.

NT is sparse at young ages while autism dominates there. For every region we
compare three per-region Autism-vs-NT estimates:
  1. naive   — Welch t-test (no age adjustment)                 [what a raw t-test gives]
  2. adjusted— OLS  y ~ diagnosis + age + age^2 + sex + cohort  [age-adjusted group effect]
  3. age>=10 — Welch t-test restricted to the common-support age range
FDR-BH within each metric. A finding robust to the imbalance should keep its
effect size / significance across (1)-(3).

Outputs: ../outputs/tables/age_imbalance_sensitivity_{metric}.csv
         ../outputs/tables/age_imbalance_summary.csv
"""
import os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
import warnings
warnings.filterwarnings('ignore')

LIB = '/Users/mfleury/POSTDOC/LIBRAIRY'
# Canonical anat input (see memory: anat-func-canonical-inputs).
ANAT_FILE = os.path.join(LIB, 'imaging2genet', '0_input', 'dataframes',
                         'MRI_ANAT_INOVAND_LEAP_COMBAT.tsv')
TABLES_DIR = Path(__file__).parent.parent / 'outputs' / 'tables'
TABLES_DIR.mkdir(parents=True, exist_ok=True)

ASEG_REGIONS = [
    "Left-Lateral-Ventricle", "Left-Thalamus", "Left-Caudate", "Left-Putamen",
    "Left-Pallidum", "3rd-Ventricle", "4th-Ventricle", "Brain-Stem",
    "Left-Hippocampus", "Left-Amygdala", "Left-VentralDC",
    "Right-Lateral-Ventricle", "Right-Thalamus", "Right-Caudate", "Right-Putamen",
    "Right-Pallidum", "Right-Hippocampus", "Right-Amygdala", "Right-VentralDC",
    "CC_Posterior", "CC_Mid_Posterior", "CC_Central", "CC_Mid_Anterior",
    "CC_Anterior", "Right-Cerebellum-Cortex", "Right-Cerebellum-White-Matter",
    "Left-Cerebellum-Cortex", "Left-Cerebellum-White-Matter",
]


def cohens_d(a, b):
    na, nb = len(a), len(b)
    if na < 2 or nb < 2:
        return np.nan
    sd = np.sqrt(((na - 1) * np.var(a, ddof=1) + (nb - 1) * np.var(b, ddof=1)) / (na + nb - 2))
    return (np.mean(a) - np.mean(b)) / sd if sd else np.nan


def feature_cols(df, metric, atlas):
    if atlas == 'dk':
        return [c for c in df.columns
                if c.startswith(('lh_', 'rh_')) and c.endswith('_' + metric)
                and 'WhiteSurfArea' not in c and 'MeanThickness' not in c]
    return [c for c in ASEG_REGIONS if c in df.columns]


def run_metric(df, metric, atlas):
    rows = []
    for c in feature_cols(df, metric, atlas):
        sub = df[[c, 'diag', 'age_c', 'age_c2', 'Sex', 'cohort']].dropna()
        a = sub[sub.diag == 'Autism'][c].values
        b = sub[sub.diag == 'NT'][c].values
        if len(a) < 3 or len(b) < 3:
            continue
        # 1. naive
        t_n, p_n = ttest_ind(a, b, equal_var=False)
        d_n = cohens_d(a, b)
        # 2. age + age^2 + sex + cohort adjusted
        s = sub.rename(columns={c: 'y'})
        m = smf.ols("y ~ C(diag, Treatment('NT')) + age_c + age_c2 + C(Sex) + C(cohort)", s).fit()
        dx = [p for p in m.params.index if 'diag' in p][0]
        t_a, p_a = m.tvalues[dx], m.pvalues[dx]
        d_a = 2 * t_a / np.sqrt(m.df_resid)
        # 3. age >= 10 naive
        s10 = sub[sub.age_c + AGE_MEAN >= 10]
        a10 = s10[s10.diag == 'Autism'][c].values
        b10 = s10[s10.diag == 'NT'][c].values
        if len(a10) >= 3 and len(b10) >= 3:
            t_10, p_10 = ttest_ind(a10, b10, equal_var=False)
            d_10 = cohens_d(a10, b10)
        else:
            t_10 = p_10 = d_10 = np.nan
        rows.append({'label': c, 'd_naive': d_n, 'p_naive': p_n,
                     'd_adj': d_a, 'p_adj': p_a, 'd_age10': d_10, 'p_age10': p_10,
                     'n_a': len(a), 'n_b': len(b)})
    res = pd.DataFrame(rows)
    for col in ('naive', 'adj', 'age10'):
        pcol = f'p_{col}'
        mask = res[pcol].notna()
        res[f'pfdr_{col}'] = np.nan
        if mask.any():
            res.loc[mask, f'pfdr_{col}'] = multipletests(res.loc[mask, pcol], method='fdr_bh')[1]
    return res


def main():
    print(f"Loading {ANAT_FILE}")
    df = pd.read_csv(ANAT_FILE, sep='\t', low_memory=False)
    df = df[df['PopulationS1'].isin(['Autism', 'TD'])].copy()
    df['diag'] = df['PopulationS1'].replace('TD', 'NT')
    df['age_yrs'] = df['age_yrs'].astype(float)
    df = df[df['age_yrs'].notna() & df['Sex'].notna() & df['cohort'].notna()].copy()
    global AGE_MEAN
    AGE_MEAN = df['age_yrs'].mean()
    df['age_c'] = df['age_yrs'] - AGE_MEAN
    df['age_c2'] = df['age_c'] ** 2
    df['Sex'] = df['Sex'].astype(float)
    # patsy can't read pandas StringDtype; Categorical is safe
    df['diag'] = pd.Categorical(df['diag'].map(str))
    df['cohort'] = pd.Categorical(df['cohort'].map(str))

    print(f"\nAutism={int((df.diag=='Autism').sum())}  NT={int((df.diag=='NT').sum())}")
    print("=== age distribution by group ===")
    for grp in ('NT', 'Autism'):
        g = df[df.diag == grp]['age_yrs']
        print(f"  {grp:6s} n={len(g):4d} median={g.median():.1f}  "
              f"<6:{int((g<6).sum())}  6-10:{int(((g>=6)&(g<10)).sum())}  >=10:{int((g>=10).sum())}")

    summary = []
    for metric, atlas in [('thickness', 'dk'), ('area', 'dk'), ('volume', 'dk'), ('volume', 'aseg')]:
        res = run_metric(df, metric, atlas)
        if res.empty:
            continue
        tag = f'{metric}_{atlas}'
        res.to_csv(TABLES_DIR / f'age_imbalance_sensitivity_{tag}.csv', index=False)
        headline = res.iloc[res['d_naive'].abs().argmax()]
        hl = res[res.label == headline['label']].iloc[0]
        x = res[['d_naive', 'd_adj', 'd_age10']].dropna()
        summary.append({
            'metric': tag, 'headline': headline['label'],
            'n_fdr_naive': int((res.pfdr_naive < .05).sum()),
            'n_fdr_adj': int((res.pfdr_adj < .05).sum()),
            'n_fdr_age10': int((res.pfdr_age10 < .05).sum()),
            'hl_d_naive': round(hl.d_naive, 3), 'hl_d_adj': round(hl.d_adj, 3),
            'hl_d_age10': round(hl.d_age10, 3),
            'hl_pfdr_naive': f'{hl.pfdr_naive:.2g}', 'hl_pfdr_adj': f'{hl.pfdr_adj:.2g}',
            'hl_pfdr_age10': f'{hl.pfdr_age10:.2g}',
            'r_naive_adj': round(np.corrcoef(x.d_naive, x.d_adj)[0, 1], 2),
            'r_naive_age10': round(np.corrcoef(x.d_naive, x.d_age10)[0, 1], 2),
        })
    sm = pd.DataFrame(summary)
    sm.to_csv(TABLES_DIR / 'age_imbalance_summary.csv', index=False)
    print("\n" + "=" * 60 + "\nAGE-IMBALANCE SENSITIVITY (ComBat input)\n" + "=" * 60)
    print(sm.to_string(index=False))
    print(f"\nWrote {TABLES_DIR / 'age_imbalance_summary.csv'}")


if __name__ == '__main__':
    main()
