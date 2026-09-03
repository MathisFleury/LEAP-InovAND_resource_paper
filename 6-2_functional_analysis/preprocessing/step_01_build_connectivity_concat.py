#!/usr/bin/env python3
"""
Concatenation variant (approach B): precision-weighted Fisher-z averaging of
the PER-RUN v0.11 connectivity, to combine INOVAND's short runs into a single
>6-min FC estimate (LEAP is single-run -> unchanged).

Why B (not XCP-D --combine-runs timeseries concat): runs differ in length and
censoring, so we treat each run's Fisher-z FC as an independent estimate and
combine by inverse-variance weighting. For Fisher-z, Var ≈ 1/(n_retained-3),
so the precision weight is w_r = n_retained_r - 3:

    z_bar(edge) = Σ_r w_r·z_r / Σ_r w_r        (over runs with a finite z_r)

This stays in Fisher-z space (the per-run raw already is), so it feeds the same
ComBat -> regress -> z-score -> analysis machinery unchanged.

Input : outputs/<variant>/df_conn_raw_<variant>.csv      (per-run, from step_01_build_connectivity)
Output: outputs/<variant>_concat/df_conn_raw_<variant>_concat.csv
        one row per (subject, session), run='concat', + minutes / n_runs / n_retained.
Run   : XCPD_VARIANT=nogsr python3.11 step_01_build_connectivity_concat.py
"""
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import CONN_RAW_CSV, OUT_BASE, VARIANT  # noqa: E402

OUT_DIR = OUT_BASE / f"{VARIANT}_concat"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_CSV = OUT_DIR / f"df_conn_raw_{VARIANT}_concat.csv"


def _linc_qc_path(ts_path: str) -> Path:
    return Path(str(ts_path).replace(
        "_seg-4S156Parcels_stat-mean_timeseries", "_desc-linc_qc"))


def _retained(ts_path: str) -> float:
    try:
        return float(pd.read_csv(_linc_qc_path(ts_path), sep="\t")["num_retained_volumes"].iloc[0])
    except Exception:
        return np.nan


def _tr(func_dir: Path) -> float:
    for j in sorted(func_dir.glob("*_desc-denoised_bold.json")):
        try:
            v = json.loads(j.read_text()).get("RepetitionTime")
            if v:
                return float(v)
        except Exception:
            continue
    return np.nan


def main() -> int:
    raw = pd.read_csv(CONN_RAW_CSV, low_memory=False)
    con = [c for c in raw.columns if c.startswith("con_")]
    raw["n_retained"] = raw["path"].map(_retained)
    print(f"[concat] {len(raw)} per-run rows; retained read for "
          f"{int(raw['n_retained'].notna().sum())}")

    out_rows, out_feats = [], []
    for (sub, cohort, ses), g in raw.groupby(["ID", "cohort", "session"], sort=False):
        Z = g[con].to_numpy(float)                      # (n_runs, n_edges), may hold NaN
        ret = g["n_retained"].to_numpy(float)
        w = np.clip(ret - 3.0, 1.0, None)               # Fisher-z precision weight
        w = np.where(np.isfinite(w), w, 1.0)
        finite = np.isfinite(Z)
        W = finite * w[:, None]
        den = W.sum(axis=0)
        num = np.nansum(np.where(finite, Z, 0.0) * W, axis=0)
        zbar = np.where(den > 0, num / den, np.nan)

        tr = _tr(Path(str(g["path"].iloc[0])).parent)
        ret_sum = float(np.nansum(ret))
        out_feats.append(zbar)
        out_rows.append({
            "ID": sub, "cohort": cohort, "session": ses, "run": "concat",
            "n_runs": int(len(g)), "n_retained": ret_sum,
            "minutes": (ret_sum * tr / 60.0) if np.isfinite(tr) else np.nan,
            "mean_fd_v11": float(np.average(g["mean_fd_v11"],
                                            weights=np.where(ret > 0, ret, 1.0))),
            "n_nan_edges": int(np.isnan(zbar).sum()),
            "path": str(g["path"].iloc[0]),
        })

    meta = pd.DataFrame(out_rows)
    feats = pd.DataFrame(np.vstack(out_feats), columns=con)
    out = pd.concat([meta.reset_index(drop=True), feats], axis=1)
    out.to_csv(OUT_CSV, index=False)

    mins = meta["minutes"].dropna()
    print(f"Wrote {OUT_CSV}")
    print(f"  subject-sessions={len(out)}  edges={len(con)}  "
          f"cohorts={meta['cohort'].value_counts().to_dict()}")
    print(f"  n_runs per session: {meta['n_runs'].value_counts().sort_index().to_dict()}")
    print(f"  minutes: median={mins.median():.2f} min={mins.min():.2f} "
          f"max={mins.max():.2f}  pct>=6={100*(mins>=6).mean():.0f}%")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
