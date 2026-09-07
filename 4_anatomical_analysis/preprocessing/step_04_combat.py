"""
Step 4 — ComBat harmonization across cohorts.

We use neuroHarmonize (preferred over neuroCombat: handles NaNs better,
modern API). Biological/nuisance covariates are PRESERVED so the batch
model does not absorb age/sex/euler signal.

    pip install neuroHarmonize
"""
import numpy as np
import pandas as pd
from neuroHarmonize import harmonizationLearn

from config import BATCH_COL, COMBAT_COVARS


def run_combat(
    df: pd.DataFrame,
    features: list[str],
) -> tuple[pd.DataFrame, dict]:
    """
    Returns the input df with feature columns replaced by harmonized values,
    plus the fitted ComBat model dictionary (for reuse / reporting).
    """
    covars = df[[BATCH_COL] + COMBAT_COVARS].rename(columns={BATCH_COL: "SITE"})
    data   = df[features].to_numpy(dtype=float)

    print(f"[04] ComBat  N={len(df):,}  features={len(features):,}  "
          f"batches={covars['SITE'].nunique()}  covars={COMBAT_COVARS}")

    model, harmonized = harmonizationLearn(data, covars)

    df_out = df.copy()
    df_out[features] = harmonized
    return df_out, model


if __name__ == "__main__":
    from step_01_load import load_dataframe
    from step_02_filter_qc import filter_qc, add_euler, drop_missing_covariates
    from step_03_select_features import select_features, restrict_to_numeric

    df = drop_missing_covariates(add_euler(filter_qc(load_dataframe())))
    feats = restrict_to_numeric(df, select_features(df))
    df_h, model = run_combat(df, feats)
    print(f"[04] done  ->  harmonized matrix shape {df_h[feats].shape}")
