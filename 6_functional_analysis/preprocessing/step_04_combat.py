"""
Step 4 — ComBat harmonization across scanners (batch = `machine_batch`).

Uses neuroHarmonize (modern API, handles NaNs, preserves biological
covariates) to match the structural pipeline_ndd. Biological/nuisance
covariates are passed in so the batch model does not absorb their signal.
COMBAT_COVARS includes age, age², sex, cohort.
"""
import numpy as np
import pandas as pd
from neuroHarmonize import harmonizationLearn

from config import BATCH_COL, COMBAT_COVARS, AGE2_COL, AGE_COL, COHORT_COL


def _prepare_covars(df_meta: pd.DataFrame) -> pd.DataFrame:
    df = df_meta.copy()
    df[AGE2_COL] = pd.to_numeric(df[AGE_COL], errors="coerce") ** 2

    # Cohort → integer code so neuroHarmonize can use it as a covariate.
    cohort_codes = {"LEAP": 1, "INOVAND": 2}
    df[COHORT_COL] = df[COHORT_COL].map(cohort_codes).astype(float)
    return df


def run_combat(
    df_conn: pd.DataFrame, df_meta: pd.DataFrame
) -> tuple[pd.DataFrame, dict]:
    df_meta = _prepare_covars(df_meta)
    assert (df_meta["ID"].astype(str).values == df_conn.index.values).all(), \
        "row order mismatch between metadata and connectivity matrix"

    covars = df_meta[[BATCH_COL] + COMBAT_COVARS].rename(columns={BATCH_COL: "SITE"})
    data = df_conn.to_numpy(dtype=float)

    print(f"[04] ComBat  N={len(df_conn):,}  features={df_conn.shape[1]:,}  "
          f"batches={covars['SITE'].nunique()}  covars={COMBAT_COVARS}")

    model, harmonized = harmonizationLearn(data, covars)
    df_out = pd.DataFrame(harmonized, index=df_conn.index, columns=df_conn.columns)
    return df_out, model


if __name__ == "__main__":
    print("Run via run_pipeline.py — depends on connectivity loader output.")
