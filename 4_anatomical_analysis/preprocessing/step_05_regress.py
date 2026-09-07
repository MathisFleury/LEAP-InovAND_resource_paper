"""
Step 5 — Per-feature linear regression of nuisance covariates.

For each MRI feature, fit:
    y = b0 + b1*age + b2*age^2 + b3*sex + b4*eTIV + b5*mean_euler + e
on the fitting subset (a reference subgroup if FIT_REFERENCE_COL is set, else
the whole sample),
then REPLACE the feature column with its residuals.

When the feature being residualised IS eTIV, the design matrix drops
eTIV (to avoid trivial self-regression).

The loop is per-feature (rather than a single lstsq across all features)
so we can skip features containing any NaN and use the right design per
feature.
"""
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression

from config import FIT_REFERENCE_COL, FIT_REFERENCE_VALUE, REGRESS_COVARS


def _fitting_subset(df: pd.DataFrame) -> pd.DataFrame:
    if FIT_REFERENCE_COL is None or FIT_REFERENCE_VALUE is None:
        print(f"[05] fitting on FULL sample (FIT_REFERENCE_COL=None)  N={len(df):,}")
        return df
    sub = df[df[FIT_REFERENCE_COL] == FIT_REFERENCE_VALUE]
    print(f"[05] fitting on reference subgroup "
          f"({FIT_REFERENCE_COL}={FIT_REFERENCE_VALUE})  N={len(sub):,}")
    return sub


def fit_and_residualise(
    df: pd.DataFrame, features: list[str]
) -> tuple[pd.DataFrame, list[str], np.ndarray]:
    """
    Returns:
        df          — input df, with feature columns replaced by residuals
        kept        — list of feature names that were residualised
                      (skipped features are not in this list)
        betas_mat   — (n_covars+1, n_kept_feats) matrix for reporting/persistence
    """
    fit_sub = _fitting_subset(df)

    kept: list[str] = []
    skipped_nan = []
    betas_list = []

    for feature in features:
        if df[feature].isna().any():
            skipped_nan.append(feature)
            continue

        # When the feature IS eTIV itself, drop eTIV from the design matrix to
        # avoid trivial self-regression (residuals would otherwise be ~0).
        covars  = [c for c in REGRESS_COVARS if c != feature]
        X_fit   = fit_sub[covars].to_numpy(dtype=float)
        X_full  = df[covars].to_numpy(dtype=float)

        y_fit = fit_sub[feature].to_numpy(dtype=float)
        reg   = LinearRegression().fit(X_fit, y_fit)
        residuals = df[feature].to_numpy(dtype=float) - reg.predict(X_full)

        df[feature] = residuals            # replace in place — no rename
        kept.append(feature)
        betas_list.append((feature, reg.intercept_, reg.coef_, covars))

    if skipped_nan:
        print(f"[05] skipped {len(skipped_nan)} feature(s) containing NaN")
    print(f"[05] residualised {len(kept):,} features in place")

    # Betas matrix: pad each column to width = 1 + |REGRESS_COVARS|. The eTIV
    # row's b_eTIV cell is NaN since eTIV was not in its own design matrix.
    cov_order = list(REGRESS_COVARS)
    betas = np.full((1 + len(cov_order), len(betas_list)), np.nan)
    for j, (_, intercept, coef, used_covars) in enumerate(betas_list):
        betas[0, j] = intercept
        for c, b in zip(used_covars, coef):
            betas[1 + cov_order.index(c), j] = b
    return df, kept, betas


if __name__ == "__main__":
    print("Run via run_pipeline.py — depends on ComBat output.")
