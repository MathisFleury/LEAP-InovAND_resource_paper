"""
Step 2 — QC filter, cohort restriction, and deduplication.

  • filter_qc: drop scans with mean_fd above FD_THRESHOLD_MM (mm).
  • restrict_to_groups: keep only ASD-bucket + control-bucket subjects
    (matches the original preprocessing_cohort.py grouping).
  • keep_best_session: collapse to one row per ID, preferring the run with
    the lowest mean_fd (deterministic; the original used drop_duplicates).
"""
import numpy as np
import pandas as pd

from config import (
    FD_THRESHOLD_MM, DEDUP_STRATEGY,
    ASD_STATUSES, CONTROL_STATUSES, KEEP_NAN_STATUS_AS_ASD,
)


def filter_qc(df: pd.DataFrame) -> pd.DataFrame:
    n0 = len(df)
    df = df[df["mean_fd"] < FD_THRESHOLD_MM].copy()
    print(f"[02] FD < {FD_THRESHOLD_MM} mm  →  kept {len(df):,} / {n0:,}")
    return df


def _group_label(status: object) -> str:
    if status in ASD_STATUSES:
        return "ASD"
    if status in CONTROL_STATUSES:
        return "Control"
    if pd.isna(status) and KEEP_NAN_STATUS_AS_ASD:
        return "ASD"
    return "Other"


def restrict_to_groups(df: pd.DataFrame) -> pd.DataFrame:
    n0 = len(df)
    df = df.copy()
    df["group"] = df["control_status"].map(_group_label)
    df = df[df["group"].isin(["ASD", "Control"])].copy()
    counts = df["group"].value_counts().to_dict()
    print(f"[02] ASD/Control restriction  →  kept {len(df):,} / {n0:,}  ({counts})")
    return df


def keep_best_session(df: pd.DataFrame) -> pd.DataFrame:
    """Collapse multi-session/run subjects to a single row.

    DEDUP_STRATEGY = 'lowest_fd' (default) keeps the run with the smallest
    mean_fd per ID — deterministic and clinically motivated. 'first' falls
    back on lexical (session, run) ordering.
    """
    n0 = len(df)
    if DEDUP_STRATEGY == "lowest_fd":
        df = (df.sort_values("mean_fd", kind="stable")
                .drop_duplicates(subset="ID", keep="first"))
    elif DEDUP_STRATEGY == "first":
        df = (df.sort_values(["session", "run"], kind="stable")
                .drop_duplicates(subset="ID", keep="first"))
    else:
        raise ValueError(f"Unknown DEDUP_STRATEGY: {DEDUP_STRATEGY}")
    print(f"[02] dedup ({DEDUP_STRATEGY}) per ID  →  kept {len(df):,} / {n0:,}  "
          f"(dropped {n0 - len(df):,})")
    return df.reset_index(drop=True)


if __name__ == "__main__":
    from step_01_load import load_dataframe
    df = keep_best_session(restrict_to_groups(filter_qc(load_dataframe())))
    print(df["cohort"].value_counts())
    print(df["group"].value_counts())
