"""
Step 2b — Build the ComBat batch column (scanner ID).

INOVAND rows already carry a `machine` integer (1, 2, 3 → Intera, Ingenia,
Ingenia 3T). LEAP rows have `LEAP_t1_site` instead, and we offset those
codes by LEAP_SITE_OFFSET so LEAP and INOVAND batches do not collide.
This mirrors the fallback used in the original preprocessing_cohort.py:

    df['machine'] = df['machine'].fillna(df['LEAP_t1_site'] + LEAP_SITE_OFFSET)
"""
import numpy as np
import pandas as pd

from config import BATCH_COL, LEAP_SITE_OFFSET


def assign_batch(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()

    machine = pd.to_numeric(df.get("machine"), errors="coerce")
    leap_site = pd.to_numeric(df.get("LEAP_t1_site"), errors="coerce")
    leap_site_offset = leap_site + LEAP_SITE_OFFSET

    batch = machine.fillna(leap_site_offset)
    n_missing = int(batch.isna().sum())
    if n_missing:
        per_cohort = df.loc[batch.isna(), "cohort"].value_counts().to_dict()
        print(f"[02b] dropping {n_missing} subjects without scanner info  ({per_cohort})")
        df = df.loc[batch.notna()].copy()
        batch = batch.loc[batch.notna()]

    df[BATCH_COL] = batch.astype(int).astype(str).values
    n_batch = df[BATCH_COL].nunique()
    print(f"[02b] {BATCH_COL} assigned  →  {n_batch} batches across {len(df):,} subjects")
    print(f"[02b] batch counts: {df[BATCH_COL].value_counts().to_dict()}")
    return df


if __name__ == "__main__":
    from step_01_load import load_dataframe
    from step_02_filter_qc import filter_qc, restrict_to_groups, keep_best_session
    df = keep_best_session(restrict_to_groups(filter_qc(load_dataframe())))
    df = assign_batch(df)
    print(df.groupby([BATCH_COL, "cohort"]).size())
