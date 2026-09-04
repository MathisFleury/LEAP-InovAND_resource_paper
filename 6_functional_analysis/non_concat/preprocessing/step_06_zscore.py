"""
Step 6 — Z-score each connectivity feature across subjects.

Reference sample = full sample (DX-stratified z-scoring is not needed for
the autism-vs-TD t-test that consumes this CSV; downstream scripts only
require harmonized, comparable units across subjects).
"""
import numpy as np
import pandas as pd


def zscore(df_conn: pd.DataFrame) -> pd.DataFrame:
    mu = df_conn.mean(axis=0)
    sd = df_conn.std(axis=0, ddof=1).replace(0, np.nan)
    df_out = (df_conn - mu) / sd
    print(f"[06] z-scored {df_conn.shape[1]:,} features across N={len(df_conn):,}")
    return df_out


if __name__ == "__main__":
    print("Run via run_pipeline.py — depends on regression output.")
