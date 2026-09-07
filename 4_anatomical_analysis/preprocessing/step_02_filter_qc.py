"""
Step 2 — Cohort restriction + QC filter + Euler covariate (NDD pipeline).
"""
import re
import numpy as np
import pandas as pd
from config import (
    QC_COL, QC_KEEP,
    AGE_COL, SEX_COL, ETIV_COL, BATCH_COL, MIN_BATCH_N,
    COHORT_KEEP, KEEP_BEST_RUN, USE_EULER_COVARIATE,
    EULER_TSV, EULER_COL,
)


def restrict_to_ndd(df: pd.DataFrame) -> pd.DataFrame:
    n0 = len(df)
    df = df[df["cohort"].isin(COHORT_KEEP)].copy()
    print(f"[02] NDD-cohort restriction  ->  kept {len(df):,} / {n0:,}  "
          f"({sorted(df['cohort'].unique())})")
    return df


def filter_qc(df: pd.DataFrame) -> pd.DataFrame:
    n0 = len(df)
    df = df[df[QC_COL].isin(QC_KEEP)].copy()
    print(f"[02] QC drop  grade 3/4/NaN  ->  kept {len(df):,} / {n0:,}")
    return df


def _norm_id(s: object) -> str:
    s = str(s)
    if s.endswith(".0"):
        s = s[:-2]
    m = re.match(r"^sub-?0*(\d+)$", s)
    if m:
        return m.group(1)
    return s.lstrip("0") or s


def add_euler(df: pd.DataFrame) -> pd.DataFrame:
    """Merge per-subject mean Euler number from recon-all logs.

    Source TSV is produced by build_euler_table.py and lives at EULER_TSV.
    Match keys: (MRI_ID, session_id) <-> (MRI_ID_n, ses_id) — both normalized.
    """
    if not USE_EULER_COVARIATE:
        print(f"[02] {EULER_COL} covariate disabled (USE_EULER_COVARIATE=False)")
        return df
    if not EULER_TSV.exists():
        raise FileNotFoundError(
            f"Euler TSV not found: {EULER_TSV}\n"
            f"Run build_euler_table.py first."
        )
    e = pd.read_csv(EULER_TSV, sep="\t", dtype={"MRI_ID_n": str, "ses_id": str})
    # We want one row per (MRI_ID_n, ses_id). The recon-all walk already
    # gives that, but defensive dedup just in case.
    e = e.drop_duplicates(["MRI_ID_n", "ses_id"], keep="first")
    e_idx = e.set_index(["MRI_ID_n", "ses_id"])["mean_euler_nofix"]

    df = df.copy()
    keys = list(zip(df["MRI_ID"].map(_norm_id), df["session_id"].map(_norm_id)))
    df[EULER_COL] = pd.Series(keys, index=df.index).map(e_idx)
    matched = df[EULER_COL].notna().sum()
    print(f"[02] {EULER_COL} from recon-all  ->  matched {matched:,}/{len(df):,}  "
          f"(median={df[EULER_COL].median():.1f}, "
          f"range=[{df[EULER_COL].min():.1f}, {df[EULER_COL].max():.1f}])")
    return df


def drop_missing_covariates(df: pd.DataFrame) -> pd.DataFrame:
    need = [AGE_COL, SEX_COL, ETIV_COL]
    if USE_EULER_COVARIATE:
        need.append(EULER_COL)
    n0 = len(df)
    df = df.dropna(subset=need).copy()
    print(f"[02] dropped missing covariates  ->  kept {len(df):,} / {n0:,}")
    return df


def keep_best_run(df: pd.DataFrame) -> pd.DataFrame:
    """When a subject (MRI_ID) has multiple runs/waves, keep the single BEST run.

    Ranking, best first:
      1. complete covariates (age / Sex / eTIV [/ mean_euler]) — so we never
         drop a subject by picking a run that lacks a covariate another run has
      2. lowest QC_seg grade (1 best ... 4 worst; NaN treated as worst)
      3. highest mean_euler (fewer surface defects = cleaner recon)
      4. earliest wave — stable tie-break on the cohort tag (T1<T2, W1<W2<W3)
    """
    if not KEEP_BEST_RUN:
        return df
    n0 = len(df)
    need = [AGE_COL, SEX_COL, ETIV_COL] + ([EULER_COL] if USE_EULER_COVARIATE else [])
    d = df.copy()
    d["_complete"] = d[need].notna().all(axis=1)
    d["_qc"] = pd.to_numeric(d[QC_COL], errors="coerce").fillna(99)
    d["_euler"] = (pd.to_numeric(d[EULER_COL], errors="coerce").fillna(-np.inf)
                   if USE_EULER_COVARIATE and EULER_COL in d.columns else 0.0)
    d = (d.sort_values(["_complete", "_qc", "_euler", "cohort"],
                       ascending=[False, True, False, True], kind="stable")
           .drop_duplicates(subset="MRI_ID", keep="first")
           .drop(columns=["_complete", "_qc", "_euler"]))
    print(f"[02] keep best-QC run per MRI_ID  ->  kept {len(d):,} / {n0:,}  "
          f"(dropped {n0 - len(d):,} duplicate runs)")
    return d


def drop_tiny_cohorts(df: pd.DataFrame) -> pd.DataFrame:
    """Drop batches below ComBat-minimum size (MIN_BATCH_N=2 by default -> keep all)."""
    n0 = len(df)
    counts = df[BATCH_COL].value_counts()
    keep   = counts[counts >= MIN_BATCH_N].index
    drop   = counts[counts < MIN_BATCH_N]
    if len(drop):
        print(f"[02] dropping batches with N < {MIN_BATCH_N}: "
              f"{dict(drop)}  ->  -{int(drop.sum())} subjects")
    df = df[df[BATCH_COL].isin(keep)].copy()
    print(f"[02] batches kept: {len(keep)}  ->  N={len(df):,} / {n0:,}")
    return df


if __name__ == "__main__":
    from step_01_load import load_dataframe
    df = restrict_to_ndd(load_dataframe())
    df = drop_missing_covariates(add_euler(filter_qc(df)))
    print(df["cohort"].value_counts())
