"""
Step 5 — Per-feature linear regression of nuisance covariates.

For every connectivity feature, fit on the full sample:
    y = β0 + β1·age + β2·age² + β3·sex + β4·cohort + ε
then replace the feature column with the residuals. The structural
pipeline_ndd uses the same approach (with eTIV / mean_euler added).
"""
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression

from config import REGRESS_COVARS, AGE_COL, AGE2_COL, COHORT_COL


def _prepare_design(df_meta: pd.DataFrame) -> pd.DataFrame:
    df = df_meta.copy()
    df[AGE2_COL] = pd.to_numeric(df[AGE_COL], errors="coerce") ** 2
    cohort_codes = {"LEAP": 1, "INOVAND": 2}
    df[COHORT_COL] = df[COHORT_COL].map(cohort_codes).astype(float)
    return df


def fit_and_residualise(
    df_conn: pd.DataFrame, df_meta: pd.DataFrame
) -> tuple[pd.DataFrame, list[str], np.ndarray]:
    df_meta = _prepare_design(df_meta)
    assert (df_meta["ID"].astype(str).values == df_conn.index.values).all()

    X = df_meta[REGRESS_COVARS].to_numpy(dtype=float)
    if np.isnan(X).any():
        n_bad = int(np.isnan(X).any(axis=1).sum())
        raise ValueError(f"design matrix has NaNs ({n_bad} rows) — check covariates")

    features = list(df_conn.columns)
    Y = df_conn.to_numpy(dtype=float)

    reg = LinearRegression().fit(X, Y)
    residuals = Y - reg.predict(X)

    df_out = pd.DataFrame(residuals, index=df_conn.index, columns=features)

    # Betas matrix (n_covars+1, n_features) for persistence / reporting.
    betas = np.vstack([reg.intercept_.reshape(1, -1), reg.coef_.T])
    print(f"[05] residualised {len(features):,} features against {REGRESS_COVARS}")
    return df_out, features, betas


if __name__ == "__main__":
    print("Run via run_pipeline.py — depends on ComBat output.")
