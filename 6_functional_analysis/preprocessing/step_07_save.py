"""
Step 7 — Save outputs.

Writes a drop-in replacement for the legacy
`results/dataframes/fmri/df_conn_cohort_norm.csv` so the downstream
01_*.py scripts can be repointed to OUT_DIR without further changes:

    df_conn_cohort_norm.csv     wide table — `con_*` features + metadata
    n_at_each_stage_<tag>.tsv   audit trail (N kept at each stage)
    features_used_<tag>.txt     list of connectivity feature columns
    combat_model_<tag>.pkl      fitted ComBat model (neuroHarmonize dict)
    regression_betas_<tag>.npy  per-feature regression coefficients
"""
from __future__ import annotations

import pickle
import numpy as np
import pandas as pd

from config import (
    OUT_NORM_CSV, OUT_AUDIT, OUT_FEATURES,
    OUT_COMBAT_MODEL, OUT_REGRESS_BETAS,
    BATCH_COL, COHORT_COL, AGE_COL, AGE2_COL, SEX_COL,
)


# Columns from the metadata frame that we keep alongside `con_*` features.
_META_KEEP = (
    "ID", "cohort", "session", "run", "path",
    "MRI_ID", "Sex", "age_yrs", "age_yrs_sq",
    "control_status", "group", "mean_fd",
    "PopulationS1", "Population1",
    "machine_batch",
    "LEAP_t1_site",
    "SRS_tscore", "total_IQ", "verbal_IQ", "performance_IQ",
)


def _assemble_output(
    df_conn: pd.DataFrame, df_meta: pd.DataFrame
) -> pd.DataFrame:
    df_conn = df_conn.copy()
    df_conn.insert(0, "ID", df_conn.index.astype(str))
    df_conn = df_conn.reset_index(drop=True)

    meta_cols = [c for c in _META_KEEP if c in df_meta.columns]
    meta = df_meta[meta_cols].copy()
    meta["ID"] = meta["ID"].astype(str)

    out = meta.merge(df_conn, on="ID", how="inner")
    return out


def save_zscored_csv(df_conn: pd.DataFrame, df_meta: pd.DataFrame) -> None:
    out = _assemble_output(df_conn, df_meta)
    out.to_csv(OUT_NORM_CSV, index=False)
    n_feat = sum(c.startswith("con_") for c in out.columns)
    print(f"[07] wrote  {OUT_NORM_CSV.name}  ({len(out):,} rows × {n_feat:,} features)")


def save_audit(n_log: dict[str, int]) -> None:
    pd.DataFrame(list(n_log.items()), columns=["stage", "N"]).to_csv(
        OUT_AUDIT, sep="\t", index=False
    )
    print(f"[07] wrote  {OUT_AUDIT.name}")


def save_features(features: list[str]) -> None:
    OUT_FEATURES.write_text("\n".join(features))
    print(f"[07] wrote  {OUT_FEATURES.name}  ({len(features):,} features)")


def save_models(combat_model: dict, betas: np.ndarray) -> None:
    with open(OUT_COMBAT_MODEL, "wb") as f:
        pickle.dump(combat_model, f)
    np.save(OUT_REGRESS_BETAS, betas)
    print(f"[07] wrote  {OUT_COMBAT_MODEL.name}  +  {OUT_REGRESS_BETAS.name}")


if __name__ == "__main__":
    print("Run via run_pipeline.py — depends on full pipeline output.")
