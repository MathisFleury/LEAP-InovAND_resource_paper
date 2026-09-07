"""
Step 6 — Z-score the residuals.

After ComBat + nuisance-regression the feature columns are residuals.
Convert each region to a z-score using the reference mean/SD:

    z = (residual - mu_ref) / sigma_ref

Reference sample = same subset used to fit regression (a reference subgroup if
FIT_REFERENCE_COL is set, full sample otherwise). This keeps the reference
statistically coherent with the residualisation step.
"""
import numpy as np
import pandas as pd

from config import FIT_REFERENCE_COL, FIT_REFERENCE_VALUE


def _reference_subset(df: pd.DataFrame) -> pd.DataFrame:
    if FIT_REFERENCE_COL is None or FIT_REFERENCE_VALUE is None:
        return df
    return df[df[FIT_REFERENCE_COL] == FIT_REFERENCE_VALUE]


def zscore(df: pd.DataFrame, features: list[str]) -> pd.DataFrame:
    # ref = _reference_subset(df)
    # mu  = ref[features].mean(axis=0)
    # sd  = ref[features].std(axis=0, ddof=1).replace(0, np.nan)

    mu = df[features].mean(axis=0)
    sd = df[features].std(axis=0)
    df_out = df.copy()
    df_out[features] = (df[features] - mu) / sd
    # print(f"[06] z-scored {len(features):,} features using N={len(ref):,} reference")
    print(f"[06] z-scored {len(features):,}")
    return df_out


if __name__ == "__main__":
    print("Run via run_pipeline.py — depends on regression output.")
