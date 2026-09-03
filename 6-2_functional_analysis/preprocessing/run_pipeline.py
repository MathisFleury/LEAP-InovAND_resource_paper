#!/usr/bin/env python3
"""
XCP-D v0.11 fMRI post-processing pipeline (6-2_functional_analysis).

Self-contained: runs entirely on 6-2's own modules — no dependency on
6_functional_analysis. The cohort is built from CURATED clinical data.

XCP-D 0.11 no longer writes connectivity matrices (only parcel mean timeseries),
so connectivity is computed here from the timeseries, then harmonised.

Pipeline, step by step
----------------------
  Step 1  step_01_build_connectivity.py     XCP-D timeseries -> raw connectivity CSV
          (or step_01_build_connectivity_    [run SEPARATELY, before this script]
           concat / _mixed.py for those modes)
  Step 2  step_02_build_cohort.build_cohort  curated clinical (diagnosis, scanner,
                                             age, Sex) + FD filter + lowest-FD dedup
                                             + ComBat batch assignment
  Step 3  step_03_combat.run_combat          ComBat scanner harmonisation
  Step 4  step_04_regress.fit_and_residualise  regress out age/age²/Sex/cohort/FD
  Step 5  step_05_zscore.zscore              per-feature z-score across subjects
  Step 6  _save                              df_conn_cohort_norm_<variant><suffix>.csv

Run modes (env vars) each add a suffix so outputs never overwrite each other —
see the block below. Example:

    XCPD_VARIANT=gsr /usr/local/bin/python3.11 run_pipeline.py
    XCPD_VARIANT=gsr CONCAT=1 COHORT_FILTER=LEAP /usr/local/bin/python3.11 run_pipeline.py
"""
from __future__ import annotations

import os
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))          # 6-2's own config + step modules

import config as cfg                    # noqa: E402
import step_02_build_cohort             # noqa: E402  (Step 2)
import step_03_combat                   # noqa: E402  (module ref for cohort-filter tweak)
import step_04_regress                  # noqa: E402
from step_03_combat import run_combat            # noqa: E402  (Step 3)
from step_04_regress import fit_and_residualise  # noqa: E402  (Step 4)
from step_05_zscore import zscore                # noqa: E402  (Step 5)

# --------------------------------------------------------------- run modes
# Each env var maps to a piece of the output-name suffix, so every variant/mode
# combination writes to its own file and directory.
#   CONCAT=1         within-session concatenated runs (approach B)
#   MIXED=1          INOVAND concat + LEAP single best run (stitched)
#   CONN_METHOD=lw   Ledoit-Wolf connectivity raw (else empirical correlation)
#   REGRESS_FIRST=1  regress -> ComBat (else ComBat -> regress)
#   COHORT_FILTER    'LEAP' | 'INOVAND'  -> single-cohort run
#   MIN_MINUTES      retained-time floor for concat/mixed (default 6)
CONCAT        = os.environ.get("CONCAT", "0") == "1"
MIXED         = os.environ.get("MIXED", "0") == "1"
REGRESS_FIRST = os.environ.get("REGRESS_FIRST", "0") == "1"
COHORT_FILTER = os.environ.get("COHORT_FILTER") or None
MIN_MINUTES   = float(os.environ.get("MIN_MINUTES", "6"))
_EST   = "_lw" if os.environ.get("CONN_METHOD", "emp").lower() == "lw" else ""
_MODE  = ("_mixed" if MIXED else "_concat" if CONCAT else "") + _EST
_ORDER = "_regfirst" if REGRESS_FIRST else ""
SUFFIX = _MODE + _ORDER + (f"_{COHORT_FILTER}" if COHORT_FILTER else "")

# Input raw (per variant + combination) and output dir (per variant + full suffix).
_RAW_DIR = cfg.OUT_BASE / f"{cfg.VARIANT}{_MODE}" if _MODE else cfg.OUT_BASE / cfg.VARIANT
cfg.CONN_RAW_CSV = _RAW_DIR / f"df_conn_raw_{cfg.VARIANT}{_MODE}.csv"
cfg.OUT_DIR = cfg.OUT_BASE / f"{cfg.VARIANT}{SUFFIX}"
cfg.OUT_DIR.mkdir(parents=True, exist_ok=True)

# Single-cohort run: `cohort` is constant, so drop it from ComBat & regression
# (harmonise across that cohort's scanner batches only, avoid a collinear covar).
if COHORT_FILTER:
    step_03_combat.COMBAT_COVARS = [c for c in step_03_combat.COMBAT_COVARS if c != "cohort"]
    step_04_regress.REGRESS_COVARS = [c for c in step_04_regress.REGRESS_COVARS if c != "cohort"]


def _save(df_conn, df_meta, n_log, combat_model, betas, features):
    """Step 6 — write the harmonised connectivity table + provenance sidecars."""
    out_csv = cfg.OUT_DIR / f"df_conn_cohort_norm_{cfg.VARIANT}{SUFFIX}.csv"
    conn = df_conn.copy()
    conn.insert(0, "ID", conn.index.astype(str))
    conn = conn.reset_index(drop=True)
    meta = df_meta.copy()
    meta["ID"] = meta["ID"].astype(str)
    out = meta.merge(conn, on="ID", how="inner")          # curated metadata + edges
    out.to_csv(out_csv, index=False)
    n_feat = sum(c.startswith("con_") for c in out.columns)
    print(f"[step 6] wrote {out_csv.name}  ({len(out):,} rows × {n_feat:,} features)")

    pd.DataFrame(list(n_log.items()), columns=["stage", "N"]).to_csv(
        cfg.OUT_DIR / f"n_at_each_stage_{cfg.VARIANT}{SUFFIX}.tsv", sep="\t", index=False)
    (cfg.OUT_DIR / f"features_used_{cfg.VARIANT}{SUFFIX}.txt").write_text("\n".join(features))
    with open(cfg.OUT_DIR / f"combat_model_{cfg.VARIANT}{SUFFIX}.pkl", "wb") as f:
        pickle.dump(combat_model, f)
    np.save(cfg.OUT_DIR / f"regression_betas_{cfg.VARIANT}{SUFFIX}.npy", betas)


def main() -> None:
    n_log = {}

    # Step 1 — load the raw connectivity built from XCP-D timeseries (built earlier).
    raw = pd.read_csv(cfg.CONN_RAW_CSV, low_memory=False)
    n_log["v11_raw_rows"] = len(raw)
    print(f"[step 1] raw connectivity: {len(raw)} run-rows  ({cfg.CONN_RAW_CSV.name})")
    if COHORT_FILTER:
        raw = raw[raw["cohort"] == COHORT_FILTER].reset_index(drop=True)
        n_log[f"cohort={COHORT_FILTER}"] = len(raw)
    if (CONCAT or MIXED) and "minutes" in raw.columns:
        n0 = len(raw)
        raw = raw[~(raw["minutes"] < MIN_MINUTES)].copy()
        print(f"[step 1] ≥{MIN_MINUTES:g} min retained: {len(raw)} / {n0} rows (NaN kept)")

    # Step 2 — build the analysis cohort from CURATED clinical data.
    # step_02_build_cohort merges the curated per-cohort clinical TSVs (diagnosis/
    # scanner/age/Sex), applies FD filter + lowest-FD dedup, and assigns the ComBat
    # scanner batch. This REPLACES the old frozen df_xcp_subjects load / QC / batch.
    min_batch = int(os.environ.get("MIN_BATCH", "10"))
    df_conn, df_meta = step_02_build_cohort.build_cohort(raw, min_batch=min_batch)
    n_log["final"] = len(df_conn)
    print(f"[step 2] curated cohort: {len(df_conn)} subjects × {df_conn.shape[1]:,} edges")

    # Steps 3–5 — harmonise, residualise, z-score (order toggled by REGRESS_FIRST).
    if REGRESS_FIRST:
        df_conn, features, betas = fit_and_residualise(df_conn, df_meta)   # Step 4 first
        df_conn, combat_model = run_combat(df_conn, df_meta)               # then Step 3
        print("[order] regress -> ComBat -> z-score")
    else:
        df_conn, combat_model = run_combat(df_conn, df_meta)               # Step 3
        df_conn, features, betas = fit_and_residualise(df_conn, df_meta)   # Step 4
        print("[order] ComBat -> regress -> z-score")
    df_conn = zscore(df_conn)                                              # Step 5

    # Step 6 — save.
    _save(df_conn, df_meta, n_log, combat_model, betas, features)
    print(f"\n✓ XCP-D v0.11 [{cfg.VARIANT}{SUFFIX}] preprocessing complete")


if __name__ == "__main__":
    main()
