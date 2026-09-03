#!/usr/bin/env python3
"""
Age-binned sanity check for the diagnosis x age interaction: N per bin (by
group) alongside the region's Autism-vs-NT effect per bin, for the headline
region from 01_age_sex_interactions.py (lowest p_fdr in
outputs/tables/interaction_anat_age_full.csv).

Same idea as 06_age_range_and_matched_sensitivity.py (age >= 10 / age-matched
restriction) but shown bin-by-bin: if the full-sample effect is being pulled
by a handful of subjects in a sparse bin, the per-bin Cohen's d should be
unstable/large exactly where N is small.

Per bin: raw Cohen's d (pooled SD) between Autism and NT, no covariates in
the comparison itself -- per user decision, 2026-08-30, to match the primary
whole-cohort test's test type exactly (Welch t-test / pooled-SD d on the
z-score table, no per-test covariates, since sex/eTIV/euler are already
regressed out upstream in the NOAGE table -- see load_anat_noage()). An
earlier version residualised on sex+cohort per bin before computing d; that
extra adjustment was dropped for consistency with the primary analysis.

Output: ../outputs/figures/age_bin_regional_effect_<region>.pdf (+ .png, .csv)
"""
import os
import importlib.util
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_SECTION_DIR = os.path.dirname(_SCRIPT_DIR)
TABLES = os.path.join(_SECTION_DIR, 'outputs', 'tables')
FIG_DIR = os.path.join(_SECTION_DIR, 'outputs', 'figures')
os.makedirs(FIG_DIR, exist_ok=True)

_spec = importlib.util.spec_from_file_location(
    'agesex', os.path.join(_SCRIPT_DIR, '01_age_sex_interactions.py'))
m = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m)

COLORS = {'Autism': '#d1495b', 'NT': '#2e6f95'}
# <6 and 6-9 kept SEPARATE (not merged into one <10 bin, 2026-08-31): merging
# made "<10" compare mostly-young autism cases against mostly-older-within-
# bracket NT cases (NT<6 is rare, n=6), a within-bin age gap that produced a
# false-positive signal in 10_age_bin_brain_maps.py. <6 is expected to be
# flagged unstable below (NT=6 < MIN_N) rather than shown as a point.
BIN_EDGES = [0, 6, 10, 14, 18, 25, 200]
BIN_LABELS = ['<6', '6-9', '10-13', '14-17', '18-24', '25+']
MIN_N = 10  # below this per-group N in a bin, d is flagged unstable
            # (matches the repo's age>=10 "adequately represented" cutoff, e.g.
            # 06_age_range_and_matched_sensitivity.py / age_imbalance_NOTE.md)


def headline_region():
    """Region with the lowest p_fdr in the full-sample age-interaction table;
    falls back to the last-known headline if the table is missing."""
    path = os.path.join(TABLES, 'interaction_anat_age_full.csv')
    if not os.path.exists(path):
        return 'lh_posteriorcingulate_thickness'
    t = pd.read_csv(path)
    return t.sort_values('p_fdr').iloc[0]['feature']


def cohens_d_ci(a, b, alpha=0.05):
    """Cohen's d (pooled SD) + Hedges-Olkin normal-approx 95% CI."""
    na, nb = len(a), len(b)
    if na < 2 or nb < 2:
        return np.nan, np.nan, np.nan
    sd = np.sqrt(((na - 1) * np.var(a, ddof=1) + (nb - 1) * np.var(b, ddof=1)) / (na + nb - 2))
    d = (np.mean(a) - np.mean(b)) / sd if sd else np.nan
    se = np.sqrt((na + nb) / (na * nb) + d ** 2 / (2 * (na + nb)))
    z = 1.959963984540054
    return d, d - z * se, d + z * se


def _selftest():
    rng = np.random.default_rng(0)
    a = rng.normal(1, 1, 200)
    b = rng.normal(0, 1, 200)
    d, lo, hi = cohens_d_ci(a, b)
    assert lo < d < hi
    assert 0.6 < d < 1.4  # true d=1, generous tolerance
    bins = pd.cut([1, 5, 9, 30], bins=BIN_EDGES, labels=BIN_LABELS, right=False)
    assert list(bins) == ['<6', '<6', '6-9', '25+']


def main():
    _selftest()
    region = headline_region()
    print(f"Headline region: {region}")

    df = m.prep_common(m.load_anat_noage(), 'population_group')
    df['age_yrs'] = df['age_yrs'].astype(float)
    df['bin'] = pd.cut(df['age_yrs'], bins=BIN_EDGES, labels=BIN_LABELS, right=False)
    df = df.dropna(subset=['bin', region])

    rows = []
    for b in BIN_LABELS:
        g = df[df['bin'] == b]
        a = g.loc[g['diag'] == 'Autism', region].values
        n = g.loc[g['diag'] == 'NT', region].values
        d, lo, hi = cohens_d_ci(a, n)
        rows.append({'bin': b, 'n_autism': len(a), 'n_nt': len(n),
                     'cohens_d': d, 'ci_lo': lo, 'ci_hi': hi,
                     'unstable': len(a) < MIN_N or len(n) < MIN_N})
    res = pd.DataFrame(rows)
    csv_path = os.path.join(FIG_DIR, f'age_bin_regional_effect_{region}.csv')
    res.to_csv(csv_path, index=False)
    print(res.round(3).to_string(index=False))

    fig, (ax_n, ax_d) = plt.subplots(2, 1, figsize=(9, 7), sharex=True,
                                      gridspec_kw={'height_ratios': [1, 1.3]})
    x = np.arange(len(BIN_LABELS))
    w = 0.35
    ax_n.bar(x - w / 2, res['n_autism'], width=w, color=COLORS['Autism'], label='Autism')
    ax_n.bar(x + w / 2, res['n_nt'], width=w, color=COLORS['NT'], label='NT')
    ax_n.set_ylabel('N')
    ax_n.set_title('Subjects per age bin', fontsize=11)
    ax_n.legend(frameon=False, fontsize=9)
    for s in ('top', 'right'):
        ax_n.spines[s].set_visible(False)

    # Bins with n < MIN_N in either group are omitted from the plot rather than
    # shown as a point: the Hedges-Olkin CI is a large-sample approximation and
    # the "pooled" SD is dominated by the other group's variance at this n, so a
    # displayed d/CI would look more precise than it is (see review discussion).
    stable = ~res['unstable']
    ax_d.errorbar(x[stable], res.loc[stable, 'cohens_d'],
                  yerr=[res.loc[stable, 'cohens_d'] - res.loc[stable, 'ci_lo'],
                        res.loc[stable, 'ci_hi'] - res.loc[stable, 'cohens_d']],
                  fmt='o', color='0.15', capsize=3, ms=7, zorder=3)
    ax_d.axhline(0, color='0.6', lw=1, ls='--')
    for i in x[~stable]:
        r = res.iloc[i]
        ax_d.text(i, 0, f"insufficient N\n(n_A={r.n_autism}, n_NT={r.n_nt})",
                  ha='center', va='center', fontsize=8, color='0.5', style='italic')
    ax_d.set_xticks(x)
    ax_d.set_xticklabels(BIN_LABELS)
    ax_d.set_xlabel('Age bin (years)')
    ax_d.set_ylabel("Cohen's d (Autism − NT)")
    ax_d.set_title(f'{region} — per-bin group effect (Welch t / pooled-SD d, no covariates)', fontsize=11)
    for s in ('top', 'right'):
        ax_d.spines[s].set_visible(False)

    fig.suptitle('Age-binned sanity check: N imbalance vs. regional effect', fontweight='bold')
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out_base = os.path.join(FIG_DIR, f'age_bin_regional_effect_{region}')
    fig.savefig(out_base + '.pdf', bbox_inches='tight')
    fig.savefig(out_base + '.png', dpi=300, bbox_inches='tight')
    print(f"Saved: {out_base}.pdf/.png, {csv_path}")


if __name__ == '__main__':
    main()
