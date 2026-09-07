"""
Step 3 — Select FreeSurfer regional measures to harmonize.

FreeSurfer summary columns have NO prefix in the input TSV. A column is
treated as a FS measure iff it ends with one of FS_SUFFIXES (e.g.
_Volume_mm3, _ThickAvg, ...) or matches a known global name. Columns used
as covariates (eTIV, surface holes) are excluded here.
"""
import pandas as pd
from config import FEATURE_EXCLUDE, FS_SUFFIXES, FS_GLOBALS, BATCH_COL


def _is_fs_measure(colname: str) -> bool:
    if colname in FS_GLOBALS:
        return True
    return any(colname.endswith(suf) for suf in FS_SUFFIXES)


def select_features(df: pd.DataFrame) -> list[str]:
    feats = [
        c for c in df.columns
        if c not in FEATURE_EXCLUDE
        and not c.startswith("corr-")        # ML-corrected variant — not used
        and _is_fs_measure(c)
    ]
    print(f"[03] selected  {len(feats):,} FreeSurfer features (raw, non corr-*)")
    return feats


def restrict_to_numeric(df: pd.DataFrame, feats: list[str]) -> list[str]:
    numeric = df[feats].apply(pd.to_numeric, errors="coerce")
    keep = [c for c in feats if numeric[c].notna().any()]
    dropped = set(feats) - set(keep)
    if dropped:
        print(f"[03] dropped {len(dropped)} non-numeric feature(s)")
    return keep


def drop_subjects_missing_features(
    df: pd.DataFrame, feats: list[str]
) -> pd.DataFrame:
    """ComBat requires complete data — drop subjects with any FS NaN."""
    n0 = len(df)
    df = df.dropna(subset=feats).copy()
    print(f"[03] dropped subjects missing any FS feature  ->  kept {len(df):,} / {n0:,}")
    return df


def drop_zero_variance_per_batch(
    df: pd.DataFrame, feats: list[str]
) -> list[str]:
    """Drop features that are constant within ANY batch — ComBat divides by
    per-batch SD, so these blow up the whole harmonization (NaN propagation).
    """
    bad = set()
    for _, sub in df.groupby(BATCH_COL):
        v = sub[feats].var(numeric_only=True)
        bad |= set(v[v == 0].index)
    feats_kept = [f for f in feats if f not in bad]
    if bad:
        print(f"[03] dropped {len(bad):,} features constant in some batch  "
              f"->  kept {len(feats_kept):,}")
    return feats_kept


if __name__ == "__main__":
    from step_01_load     import load_dataframe
    from step_02_filter_qc import filter_qc, add_euler, drop_missing_covariates
    df = drop_missing_covariates(add_euler(filter_qc(load_dataframe())))
    feats = restrict_to_numeric(df, select_features(df))
    df = drop_subjects_missing_features(df, feats)
    print(f"[03] final: {len(df):,} subjects x {len(feats):,} features")
