#!/usr/bin/env python3
"""
Clinical x MRI region-wise Pearson correlations (REVISION anat + CURATED clinical)

For the autism cohort, correlate every regional MRI measure (cortical thickness /
surface area / cortical volume / subcortical volume) with four clinical features
— SRS-2 t-score, RBS-R total, full-scale IQ, VABS-II ABC — using Pearson r.

MRI      : revised FreeSurfer z-score+ComBat+regression table, loaded (with the
           tidy region columns + canonical ID) via 05's loader.
Clinical : the CURATED cluster table (curated SRS/RBS-R/IQ/VABS + PopulationS1),
           bridged to the MRI by the same canonical join key.

Output (feeds 13_clinical_mri_brain_maps_v2.R):
  ../outputs/tables_clinical_mri/clinmri_<clin>_<feature>.csv
  columns: mri_col, mri_type, hemisphere, region, correlation, p_value, p_fdr, sample_size
"""

import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, zscore

sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
_m05 = import_module("5_loeuf_mri_correlations_hg38_v2")
_mod01 = import_module("1_anatomical_mri_autism_nt_v2")
import _config  # noqa: E402
sys.path.insert(0, os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "..")))
import config as root_config  # noqa: E402  (for the curated cluster table)

SUBCORTICAL_COLS = list(_m05.SUBCORTICAL_COLS)
TABLES_DIR = Path(_config.OUTPUT_BASE) / "tables_clinical_mri"

# clinical: output token -> (column in curated table, display label)
CLIN = [("srs",  "SRS_tscore",          "SRS-2 t-score"),
        ("rbsr", "RBS-R_total",         "RBS-R"),
        ("fsiq", "total_IQ",            "FSIQ"),
        ("vabs", "vabsabcabc_standard", "VABS-II")]


def build_dataset():
    mri = _m05.load_curated_mri_with_canonical_id()          # region cols + _join_key
    mri = mri.drop(columns=[c for c in ("PopulationS1", "cohort") if c in mri.columns])

    cur = pd.read_csv(root_config.CLUSTERS_CURATED, low_memory=False)
    cur["_join_key"] = cur.apply(_mod01._canonical_join_key_clusters, axis=1)
    keep = ["_join_key", "PopulationS1"] + [c[1] for c in CLIN]
    cur = (cur[[c for c in keep if c in cur.columns]]
           .dropna(subset=["_join_key"]).drop_duplicates("_join_key"))

    df = cur.merge(mri, on="_join_key", how="inner")
    print(f"  MRI x curated-clinical merged: {len(df)}")
    df = df[df["PopulationS1"] == "Autism"].copy()
    print(f"  autism only: {len(df)}")

    ct_cols = [c for c in df.columns
               if (c.endswith("_thickness") or c.endswith("_area") or c.endswith("_grayvol"))
               and (c.startswith("lh_") or c.startswith("rh_"))
               and "MeanThickness" not in c and "WhiteSurfArea" not in c]
    sc_cols = [c for c in SUBCORTICAL_COLS if c in df.columns]
    print(f"  cortical columns: {len(ct_cols)} · subcortical columns: {len(sc_cols)}")
    return df, ct_cols, sc_cols


def correlate(df, clin_col, mri_cols):
    rows = []
    y_full = pd.to_numeric(df[clin_col], errors="coerce")
    for mri_col in mri_cols:
        x = pd.to_numeric(df[mri_col], errors="coerce")
        ok = x.notna() & y_full.notna()
        if ok.sum() < 10:
            continue
        r, p = pearsonr(zscore(x[ok].values), zscore(y_full[ok].values))
        m = re.match(r"^(lh|rh)_(.+)_(thickness|area|grayvol)$", mri_col)
        if m:
            hemi = "left" if m.group(1) == "lh" else "right"
            region, mtype = m.group(2), m.group(3)
        else:
            hemi = ("left" if mri_col.startswith("Left-") else
                    "right" if mri_col.startswith("Right-") else "")
            region, mtype = mri_col, "subcortical"
        rows.append({"mri_col": mri_col, "mri_type": mtype, "hemisphere": hemi,
                     "region": region, "correlation": float(r), "p_value": float(p),
                     "sample_size": int(ok.sum())})
    return pd.DataFrame(rows)


def main():
    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    df, ct_cols, sc_cols = build_dataset()
    all_cols = list(ct_cols) + list(sc_cols)
    for tok, col, lab in CLIN:
        if col not in df.columns:
            print(f"  skip {lab}: {col} missing"); continue
        s = correlate(df, col, all_cols)
        if s.empty:
            print(f"  {lab}: no correlations"); continue
        # BH-FDR within each mri_type (one ggseg facet), matching 06's pooling
        s["p_fdr"] = np.nan
        for mt, sub in s.groupby("mri_type"):
            from statsmodels.stats.multitest import multipletests
            s.loc[sub.index, "p_fdr"] = multipletests(sub["p_value"], method="fdr_bh")[1]
        for mt, sub in s.groupby("mri_type"):
            out = TABLES_DIR / f"clinmri_{tok}_{mt}.csv"
            sub.to_csv(out, index=False)
        n_sig = int((s["p_fdr"] < 0.05).sum())
        print(f"  {lab}: {len(s)} regions · FDR q<0.05 {n_sig} · n≈{int(s['sample_size'].median())}")
    print(f"\nWritten to: {TABLES_DIR}")


if __name__ == "__main__":
    main()
