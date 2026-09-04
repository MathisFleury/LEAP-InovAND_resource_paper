#!/usr/bin/env python3
"""
LOEUF (v4) × MRI — OLS regression + permutation test (REVISION, paper method)

Reproduces the first-iteration paper statistic (distinct from
5_loeuf_mri_correlations_hg38_v2.py, which uses Pearson r + parametric p):

  For each gene set and each MRI feature/ROI:
    - OLS regression of the MRI metric on -log10(LOEUF) best-score:
        metric ~ const + (-log10 LOEUF)   ->  beta, SE, parametric p
    - Permutation test (1000 iters, two-tailed, seed 43) with a
      GENE-LIST-SPECIFIC null:
        * Protein coding : shuffle X (random variant assignment).
        * SynGO / CHROM  : draw X from the protein-coding background
                           distribution (resample with replacement).
      p_perm = (#|perm beta| as extreme as observed + 1) / (n_perm + 1) * 2.

Gene sets : Protein coding (also the background), SynGO, ChromEpiTF.
Carriers  : autistic only when AUTISM_ONLY (reuses 05's loaders / roster).

Outputs (NEW dirs — nothing in 05's tables_genetics_hg38 is touched):
  ../outputs/tables_genetics_hg38_regperm/loeuf_reg_<pathway>_<feature>.csv
"""

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats as sps

# Reuse the loaders / roster / carrier logic from script 05.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
_m05 = import_module("5_loeuf_mri_correlations_hg38_v2")  # noqa: E402
import _config  # noqa: E402  shared paths + hyperparameters (see _config.py)

neg_log10        = _m05.neg_log10
SUBCORTICAL_COLS = list(_m05.SUBCORTICAL_COLS)

# Centralised in _config.py (shared with 05 so the two LOEUF analyses stay in
# sync). LOEUF_OUT_DIR overrides the output base dir if a copy is wanted.
AUTISM_ONLY = _config.AUTISM_ONLY
N_PERM      = _config.N_PERM
SEED        = _config.SEED
TABLES_DIR  = Path(_config.OUTPUT_BASE) / "tables_genetics_hg38_regperm"

# (output feature name, lof best-score col, del best-score col, label, short)
GENE_SETS = [
    ("dellof_proteincoding_best_score", "lof_proteincoding_genes_best_score",
     "del_Gene_name_best_score",        "Protein coding", "proteincoding"),
    ("dellof_syngo_best_score",         "lof_syngo_genes_best_score",
     "del_syngo_genes_best_score",      "SynGO",          "syngo"),
    ("dellof_chromepitf_best_score",    "lof_chromepitf_genes_best_score",
     "del_chromepitf_genes_best_score", "CHROM",          "chromepitf"),
]
PROTEIN_FEATURE = "dellof_proteincoding_best_score"  # background for pathway nulls


def build_dataset():
    mri  = _m05.load_curated_mri_with_canonical_id()
    gens = _m05.load_carriers_hg38()
    mri = mri.drop(columns=[c for c in ("cohort", "PopulationS1") if c in mri.columns])
    df = gens.merge(mri, on="ID", how="inner")
    print(f"  merged: {len(df)} individuals")

    if AUTISM_ONLY:
        before = len(df)
        df = df[df["PopulationS1"] == "Autism"].copy()
        print(f"  AUTISM-ONLY filter: {before} -> {len(df)} autistic individuals")

    for new_col, lof_col, del_col, _, _ in GENE_SETS:
        lof_s = df[lof_col] if lof_col in df.columns else pd.Series(np.nan, index=df.index)
        del_s = df[del_col] if del_col in df.columns else pd.Series(np.nan, index=df.index)
        df[new_col] = np.fmin(lof_s.values, del_s.values)

    ct_cols = [c for c in df.columns
               if (c.endswith("_thickness") or c.endswith("_area")
                   or c.endswith("_grayvol"))
               and (c.startswith("lh_") or c.startswith("rh_"))
               and "MeanThickness" not in c and "WhiteSurfArea" not in c]
    sc_cols = [c for c in SUBCORTICAL_COLS if c in df.columns]
    print(f"  cortical columns: {len(ct_cols)} · subcortical columns: {len(sc_cols)}")
    return df, ct_cols, sc_cols


def _ols_beta_se_p(x, y):
    """OLS y ~ const + x. Returns (beta, se, p)."""
    model = sm.OLS(y, sm.add_constant(x)).fit()
    return float(model.params[1]), float(model.bse[1]), float(model.pvalues[1])


def _residuals_not_normal(resid):
    """Four complementary normality tests (paper method); flag if any p<0.05."""
    flags = []
    if len(resid) < 5000:
        flags.append(sps.shapiro(resid).pvalue < 0.05)
    flags.append(sps.normaltest(resid).pvalue < 0.05)          # D'Agostino
    flags.append(sps.jarque_bera(resid).pvalue < 0.05)
    ad = sps.anderson(resid, dist="norm")                      # 5% crit = index 2
    flags.append(ad.statistic > ad.critical_values[2])
    return any(flags)


def _bootstrap_median_beta(X, y, n=1000):
    """Robust beta: median of 1000 resample-with-replacement OLS slopes
    (seed 42, paper method / plot_final_figure.py)."""
    np.random.seed(42)
    betas = []
    for _ in range(n):
        idx = np.random.choice(len(X), size=len(X), replace=True)
        try:
            betas.append(sm.OLS(y[idx], sm.add_constant(X[idx])).fit().params[1])
        except Exception:
            continue
    return float(np.median(betas)) if betas else None


def _perm_slopes(Xp, y):
    """Vectorised simple-regression slopes for a (n_perm, n) matrix Xp vs y."""
    xc = Xp - Xp.mean(axis=1, keepdims=True)
    yc = y - y.mean()
    denom = (xc ** 2).sum(axis=1)
    denom[denom == 0] = np.nan
    return (xc * yc).sum(axis=1) / denom


def run(df, ct_cols, sc_cols):
    rng = np.random.default_rng(SEED)
    background = neg_log10(df[PROTEIN_FEATURE]).dropna().values
    rows = []
    for feat, _lof, _del, gl_label, short in GENE_SETS:
        if feat not in df.columns:
            print(f"  skip {feat}: missing"); continue
        x_full = neg_log10(df[feat])
        n_carriers = len(x_full)
        is_protein = (feat == PROTEIN_FEATURE)
        for mri_col in list(ct_cols) + list(sc_cols):
            mri = df[mri_col].dropna()
            common = x_full.index.intersection(mri.index)
            if len(common) < 10:
                continue
            X = x_full.loc[common].values.astype(float)
            y = mri.loc[common].values.astype(float)
            model = sm.OLS(y, sm.add_constant(X)).fit()
            beta, se, p_param = float(model.params[1]), float(model.bse[1]), float(model.pvalues[1])

            # Robust beta: bootstrap median when OLS residuals fail normality
            # (paper method). `beta` stays OLS so the brain maps in 09 keep their
            # "OLS beta fill"; 10 plots `beta_robust`.
            not_normal = _residuals_not_normal(model.resid)
            boot = _bootstrap_median_beta(X, y) if not_normal else None
            beta_robust = boot if boot is not None else beta

            # Gene-list-specific null
            if is_protein:
                Xp = np.array([rng.permutation(X) for _ in range(N_PERM)])
            else:
                Xp = rng.choice(background, size=(N_PERM, len(X)), replace=True)
            perm_betas = _perm_slopes(Xp, y)
            perm_betas = perm_betas[np.isfinite(perm_betas)]

            if perm_betas.size:
                if beta >= 0:
                    n_ext = int(np.sum(perm_betas >= beta))
                else:
                    n_ext = int(np.sum(perm_betas <= beta))
                p_perm = min((n_ext + 1) / (N_PERM + 1) * 2, 1.0)
                z = ((beta - perm_betas.mean()) / perm_betas.std()
                     if perm_betas.std() > 0 else 0.0)
            else:
                p_perm, z = np.nan, np.nan

            m = re.match(r"^(lh|rh)_(.+)_(thickness|area|grayvol)$", mri_col)
            if m:
                hemi = "left" if m.group(1) == "lh" else "right"
                region, mtype = m.group(2), m.group(3)
            else:
                hemi = ("left" if mri_col.startswith("Left-")
                        else "right" if mri_col.startswith("Right-") else "")
                region, mtype = mri_col, "subcortical"
            rows.append({
                "genetic_type": "LOEUF", "genetic_feature": feat,
                "gene_list_label": gl_label, "pathway_short": short,
                "mri_col": mri_col, "mri_type": mtype,
                "hemisphere": hemi, "region": region,
                "sample_size": int(len(common)), "n_total_carriers": int(n_carriers),
                "beta": beta, "beta_robust": beta_robust, "se": se,
                "residuals_normal": (not not_normal),
                "p_value": p_param, "p_permutation": float(p_perm),
                "perm_z": float(z),
            })
        print(f"  {gl_label}: {n_carriers} carriers"
              f"{' (background+self via label-shuffle)' if is_protein else ''}")
    return pd.DataFrame(rows)


def main():
    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    print("=" * 60)
    print("LOEUF x MRI — OLS + permutation (paper method)")
    print(f"Output: {TABLES_DIR}")
    print(f"AUTISM_ONLY={AUTISM_ONLY} · N_PERM={N_PERM} · seed={SEED}")
    print("=" * 60)

    df, ct_cols, sc_cols = build_dataset()
    stats = run(df, ct_cols, sc_cols)
    if stats.empty:
        print("No stats produced."); return

    print("\nWriting per-pathway × per-feature files:")
    for (short, mtype), sub in stats.groupby(["pathway_short", "mri_type"]):
        out = TABLES_DIR / f"loeuf_reg_{short}_{mtype}.csv"
        sub.to_csv(out, index=False)
        n_sig = int((sub["p_permutation"] < 0.05).sum())
        print(f"  {out.name}: {len(sub)} ROIs · perm p<0.05 {n_sig}")


if __name__ == "__main__":
    main()
