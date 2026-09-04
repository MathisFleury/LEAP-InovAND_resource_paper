#!/usr/bin/env python3
"""
Mixed variant: INOVAND = concat (precision-weighted Fisher-z avg of within-session
runs, approach B), LEAP = single best run. This is the original 6-1 intent
("concat INOVAND + the one run from LEAP"), now on v0.11.

Just stitches the two existing raw tables (no recompute):
  LEAP rows    <- per-run raw      outputs/<variant>/df_conn_raw_<variant>.csv
  INOVAND rows <- concat raw       outputs/<variant>_concat/df_conn_raw_<variant>_concat.csv

LEAP rows keep one row per run (the pipeline's dedup picks the best run);
INOVAND rows are one row per (subject, session) with `minutes` for the ≥6-min
filter. Columns are unioned (LEAP gets NaN minutes/n_runs/n_retained), which the
pipeline's minutes filter tolerates (NaN-minutes rows are kept).

Output: outputs/<variant>_mixed/df_conn_raw_<variant>_mixed.csv
Run   : XCPD_VARIANT=nogsr python3.11 step_01_build_connectivity_mixed.py
"""
import os
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import OUT_BASE, VARIANT  # noqa: E402

PERRUN = OUT_BASE / VARIANT / f"df_conn_raw_{VARIANT}.csv"
CONCAT = OUT_BASE / f"{VARIANT}_concat" / f"df_conn_raw_{VARIANT}_concat.csv"
OUT_DIR = OUT_BASE / f"{VARIANT}_mixed"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_CSV = OUT_DIR / f"df_conn_raw_{VARIANT}_mixed.csv"


def main() -> int:
    perrun = pd.read_csv(PERRUN, low_memory=False)
    concat = pd.read_csv(CONCAT, low_memory=False)

    leap = perrun[perrun["cohort"] == "LEAP"].copy()                 # per-run (best picked downstream)
    ino = concat[concat["cohort"].isin(["INOVAND", "INFOR"])].copy()  # session-concat (INFOR folds in)

    out = pd.concat([leap, ino], ignore_index=True, sort=False)
    out.to_csv(OUT_CSV, index=False)

    print(f"Wrote {OUT_CSV}")
    print(f"  LEAP per-run rows: {len(leap)}  (subjects {leap['ID'].nunique()})")
    print(f"  INOVAND concat rows: {len(ino)}  (subjects {ino['ID'].nunique()})")
    print(f"  total rows={len(out)}  edges={sum(c.startswith('con_') for c in out.columns)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
