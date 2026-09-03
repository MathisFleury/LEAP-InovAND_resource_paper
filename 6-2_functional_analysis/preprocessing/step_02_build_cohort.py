#!/usr/bin/env python3
"""
Curated-clinical cohort builder for the v0.11 fMRI pipeline.

Replaces the frozen df_xcp_subjects_* / qc_xcp_d_* files entirely. Demographics,
diagnosis (population_group) and SCANNER (mri_machine) come from the curated
clinical TSVs; connectivity + motion (mean_fd_v11) come from the v0.11 raw built
by step_01_build_connectivity.py. This recovers every v0.11 subject (e.g. INOVAND
573, incl. scanners/sites missing from the old 420-subject list) and lets ComBat
batch on the true scanner.

Keys
  INOVAND : v0.11 sub-id == curated MRI_ID ('sub-XXXX'); analysis ID = curated `ID`
  LEAP    : v0.11 sub-id 'sub-<n>' -> ID <n>; wave by session (ses-0k -> t{k}, fallback t1)

ComBat batch = mri_machine (string, unique across cohorts -> factorised). Scanners
with < MIN_BATCH subjects in the final sample are merged into 'OTHER_<cohort>'.
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd

_CUR = Path(os.environ.get("IMG5_ROOT", "/Volumes/Imaging5/EEG_MRI-MF"))
INOVAND_MERGED = _CUR / "INOVAND/_clinical_data/curated/INOVAND_clinical_curated_merged.tsv"
INFOR_CURATED = _CUR / "INFOR/_clinical_data/curated/INFOR_clinical_curated.tsv"
LEAP_WAVES = {f"ses-0{k}": _CUR / f"LEAP/_clinical_data/curated/LEAP_clinical_curated_t{k}.tsv"
              for k in (1, 2, 3)}
_COLS = ["age_yrs", "Sex", "population_group", "mri_machine", "site_id"]
VALID_DX = {"Autism", "NT", "Relatives", "IDD", "Potential autism"}


def _sex(s):
    sl = str(s).strip().lower()
    if sl in ("m", "male", "1", "1.0"):
        return 1.0
    if sl in ("f", "female", "0", "0.0"):
        return 0.0
    try:
        v = float(s)
        return v if v in (0.0, 1.0) else np.nan
    except Exception:
        return np.nan


def _mri_keyed_lookup(path: Path) -> pd.DataFrame:
    """Curated TSV keyed by BIDS MRI_ID ('sub-XXXX'); analysis id = its `ID`.
    Used for INOVAND and INFOR (both one row per scan)."""
    d = pd.read_csv(path, sep="\t", low_memory=False)
    d = d[d["MRI_ID"].notna()].copy()
    d["sub_id"] = d["MRI_ID"].astype(str).str.strip()
    d["_h"] = d["mri_machine"].notna().astype(int)
    d = d.sort_values("_h", ascending=False).drop_duplicates("sub_id", keep="first")
    return d[["sub_id", "ID", *_COLS]].rename(columns={"ID": "ana_id"})


def _leap_lookup() -> pd.DataFrame:
    frames = []
    for ses, path in LEAP_WAVES.items():
        if not path.exists():
            continue
        d = pd.read_csv(path, sep="\t", low_memory=False)
        d["ID"] = d["ID"].astype(str).str.strip()
        d = d.drop_duplicates("ID", keep="first")
        d["session"] = ses
        frames.append(d[["ID", "session", *_COLS]])
    leap = pd.concat(frames, ignore_index=True)
    leap = leap.rename(columns={"ID": "ana_id"})
    leap["cohort"] = "LEAP"
    return leap


def build_cohort(raw: pd.DataFrame, fd_thresh: float = 0.5, min_batch: int = 10
                 ) -> tuple[pd.DataFrame, pd.DataFrame]:
    con = [c for c in raw.columns if c.startswith("con_")]
    raw = raw.copy()
    raw["sub_id"] = raw["ID"].astype(str)

    # INOVAND + INFOR: keyed by BIDS MRI_ID. INFOR is a separate acquisition
    # SOURCE (distinct scanners) but folded into the INOVAND biological cohort.
    ino_raw = raw[raw["cohort"] == "INOVAND"].merge(
        _mri_keyed_lookup(INOVAND_MERGED), on="sub_id", how="left")
    mri_frames = [ino_raw]
    if (raw["cohort"] == "INFOR").any() and INFOR_CURATED.exists():
        infor_raw = raw[raw["cohort"] == "INFOR"].merge(
            _mri_keyed_lookup(INFOR_CURATED), on="sub_id", how="left")
        mri_frames.append(infor_raw)
        print(f"[curated] INFOR rows bridged: {len(infor_raw)}")

    leap = raw[raw["cohort"] == "LEAP"].copy()
    leap["ana_id"] = leap["sub_id"].str.replace("sub-", "", regex=False)
    L = _leap_lookup()
    # session-matched wave, then fallback to any wave (earliest) for that ID
    leap_m = leap.merge(L.drop(columns="cohort"), on=["ana_id", "session"], how="left")
    miss = leap_m["age_yrs"].isna()
    if miss.any():
        fb = L.sort_values("session").drop_duplicates("ana_id", keep="first")
        fb = fb.drop(columns=["cohort", "session"])
        leap_m = leap_m.merge(fb, on="ana_id", how="left", suffixes=("", "_fb"))
        for c in _COLS:
            leap_m[c] = leap_m[c].fillna(leap_m[f"{c}_fb"])
        leap_m = leap_m.drop(columns=[f"{c}_fb" for c in _COLS])

    df = pd.concat([*mri_frames, leap_m], ignore_index=True, sort=False)

    # scanner SOURCE (LEAP/INOVAND/INFOR) drives batching; biological COHORT
    # (LEAP/INOVAND) is the ComBat/regression covariate — INFOR folds into INOVAND.
    df["src"] = df["cohort"]
    df["cohort"] = df["cohort"].replace({"INFOR": "INOVAND"})

    # --- clean covariates / restrict ---
    df["Sex"] = df["Sex"].map(_sex)
    df["age_yrs"] = pd.to_numeric(df["age_yrs"], errors="coerce")
    df["mean_fd"] = pd.to_numeric(df["mean_fd_v11"], errors="coerce")
    n0 = len(df)
    df = df[df["age_yrs"].notna() & df["Sex"].notna() & df["mri_machine"].notna()
            & df["population_group"].isin(VALID_DX)].copy()
    print(f"[curated] valid age/Sex/machine/dx: {len(df)} / {n0} rows")
    df = df[df["mean_fd"] < fd_thresh].copy()
    print(f"[curated] FD<{fd_thresh}: {len(df)} rows  ({df['cohort'].value_counts().to_dict()})")

    # --- one row per subject: prefer clean run, then lowest FD ---
    df["is_clean"] = (df["n_nan_edges"] == 0).astype(int)
    sort_cols, asc = ["is_clean", "mean_fd"], [False, True]
    if "minutes" in df.columns:           # concat/mixed: also prefer more retained time
        sort_cols, asc = ["is_clean", "minutes", "mean_fd"], [False, False, True]
    df = (df.sort_values(sort_cols, ascending=asc)
            .drop_duplicates("ana_id", keep="first").reset_index(drop=True))
    print(f"[curated] one row per subject: {len(df)}  ({df['cohort'].value_counts().to_dict()})")

    # drop residual NaN-edge subjects -> keep full edge grid
    nbad = int((df["n_nan_edges"] > 0).sum())
    if nbad:
        df = df[df["n_nan_edges"] == 0].reset_index(drop=True)
        print(f"[curated] dropped {nbad} subjects with NaN edges")

    # --- batch = source × mri_machine (INFOR scanners distinct from INOVAND's
    #     even if same model); merge tiny scanners into per-source OTHER ---
    df["mri_machine"] = df["mri_machine"].astype(str).str.strip()
    n_by = df.groupby(["src", "mri_machine"])["ana_id"].transform("size")
    df["batch_label"] = np.where(n_by < min_batch,
                                 "OTHER_" + df["src"], df["mri_machine"])
    df["machine_batch"] = pd.factorize(df["src"] + "|" + df["batch_label"])[0].astype(str)
    merged_small = df.loc[n_by < min_batch].groupby(["src", "mri_machine"]).size().to_dict()
    if merged_small:
        print(f"[curated] scanners <{min_batch} merged into OTHER_<src>: {merged_small}")
    print(f"[curated] ComBat batches: {df.groupby('src')['machine_batch'].nunique().to_dict()}  "
          f"({df['machine_batch'].nunique()} total)")
    print(f"[curated] scanner counts:\n{df.groupby(['src','mri_machine']).size().to_string()}")

    df_conn = df[con].copy()
    df_conn.index = df["ana_id"].astype(str).values
    meta_cols = ["ana_id", "sub_id", "cohort", "src", "session", "run", "mean_fd", "age_yrs",
                 "Sex", "population_group", "mri_machine", "site_id", "machine_batch",
                 "n_nan_edges", "path_v11" if "path_v11" in df.columns else "path"]
    meta_cols += [c for c in ("minutes", "n_runs", "n_retained") if c in df.columns]
    df_meta = df[[c for c in meta_cols if c in df.columns]].copy()
    df_meta = df_meta.rename(columns={"ana_id": "ID"})
    return df_conn, df_meta


if __name__ == "__main__":
    import sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from config import CONN_RAW_CSV
    raw = pd.read_csv(CONN_RAW_CSV, low_memory=False)
    dc, dm = build_cohort(raw)
    print("\nfinal:", dc.shape, "| Autism/NT:",
          dm["population_group"].value_counts().to_dict())
