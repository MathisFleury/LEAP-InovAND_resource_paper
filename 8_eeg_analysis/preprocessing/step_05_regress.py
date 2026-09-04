"""
Step 5 — Regress nuisance covariates out of each EEG feature.

Per feature, fit on the FULL sample (DX_COL=None):
    y = β0 + β1·age + β2·age² + β3·sex + ε
then replace the feature column with the residuals. Mirrors the anat pipeline
step_05_regress (which additionally carries eTIV / euler for structural MRI).

Fit is per-feature so that features with differing missing-data patterns each
use their own valid rows.
"""
import numpy as np
from sklearn.linear_model import LinearRegression

from config import REGRESS_COVARS


def regress(df, features, covars=REGRESS_COVARS, min_n=10):
    df = df.copy()
    if "age_yrs_sq" in covars and "age_yrs_sq" not in df.columns:
        df["age_yrs_sq"] = df["age_yrs"] ** 2

    X = df[covars].copy()
    if "Sex" in X.columns:
        X["Sex"] = X["Sex"].fillna(0)

    for feat in features:
        mask = df[feat].notna() & X.notna().all(axis=1)
        if int(mask.sum()) < min_n:
            continue
        model = LinearRegression().fit(X[mask], df.loc[mask, feat])
        df.loc[mask, feat] = df.loc[mask, feat] - model.predict(X[mask])

    print(f"[05] residualised {len(features)} feature(s) vs {covars} (full sample)")
    return df


if __name__ == "__main__":
    # ponytail self-check: a covariate-driven signal is removed, mean → ~0
    import pandas as pd
    rng_age = np.linspace(5, 40, 200)
    demo = pd.DataFrame({
        "age_yrs": rng_age,
        "Sex": [0, 1] * 100,
        "feat": 2.0 * rng_age + 3.0,  # pure age effect
    })
    out = regress(demo, ["feat"])
    assert abs(out["feat"].mean()) < 1e-6, out["feat"].mean()
    assert out["feat"].std() < 1e-6, "linear age effect should be fully removed"
    print("ok")
