#!/usr/bin/env python3
"""
Cross-reference the fMRI subjects (df_conn_cohort_norm_<variant>.csv) against
the anatomical FreeSurfer QC rating (QC_seg, 1=best … 4=fail) and flag those
that did NOT pass anat QC (QC_seg > 1).

Anat QC source: the qc1234 z-score table (all QC levels retained), which still
carries the per-scan QC_seg column. A subject is taken at its WORST (max)
QC_seg across anat sessions (conservative for exclusion).

ID bridge (cohort-specific):
  LEAP     : fMRI ID (numeric)        ↔ anat MRI_ID (numeric)
  INOVAND  : fMRI MRI_ID ('sub-0290') ↔ anat MRI_ID ('290')   [strip sub-, leading zeros]

Output: qc/anat_qc_flags_<variant>.csv   (one row per fMRI subject + QC_seg, pass flag)
Run:    XCPD_VARIANT=nogsr python3.11 anat_qc_crossref.py
"""
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

VARIANT = os.environ.get("XCPD_VARIANT", "nogsr")
HERE = Path(__file__).resolve().parent
CONN = HERE.parent / "preprocessing" / "outputs" / VARIANT / f"df_conn_cohort_norm_{VARIANT}.csv"

_PROJECT_DIR = HERE.parent.parent.parent  # = LEAP-InovAND_resource/
if str(_PROJECT_DIR) not in sys.path:
    sys.path.insert(0, str(_PROJECT_DIR))
from config import MRI_CURATED  # noqa: E402 -- single source of truth for the anat dir

# qc1234 (all QC levels retained, unlike the qc1-only MRI_CURATED default) --
# still carries the per-scan QC_seg column needed here. Same dir as MRI_CURATED
# (built in-repo by 4_anatomical_analysis/preprocessing/run_pipeline.py).
ANAT_QC = Path(os.environ.get(
    "ANAT_QC_FILE",
    os.path.join(os.path.dirname(MRI_CURATED), "freesurfer_zscore_qc1234_combat_regress.tsv")))


def _to_int(x):
    try:
        return str(int(float(re.sub(r"^sub-", "", str(x)))))
    except Exception:
        return None


def main() -> int:
    anat = pd.read_csv(ANAT_QC, sep="\t")
    anat["base"] = anat["cohort"].astype(str).str.split("_").str[0]
    anat["key"] = anat["base"] + "_" + anat["MRI_ID"].map(_to_int)
    qc = (anat.groupby("key")
              .agg(QC_seg=("QC_seg", "max"),
                   n_anat_sessions=("QC_seg", "size"),
                   mean_euler=("mean_euler", "min"))
              .reset_index())

    fmri = pd.read_csv(CONN, usecols=["ID", "MRI_ID", "cohort", "population_group"],
                       low_memory=False)
    fmri["key"] = np.where(fmri["cohort"].eq("INOVAND"),
                           "INOVAND_" + fmri["MRI_ID"].map(_to_int),
                           "LEAP_" + fmri["ID"].map(_to_int))

    m = fmri.merge(qc, on="key", how="left")
    m["anat_qc_matched"] = m["QC_seg"].notna()
    m["passed_anat_qc"] = m["QC_seg"].le(1)        # QC_seg==1 passes; >1 fails; NaN -> False

    out = HERE / f"anat_qc_flags_{VARIANT}.csv"
    m.to_csv(out, index=False)

    fail = m[m["QC_seg"] > 1]
    print(f"fMRI subjects: {len(m)}  |  matched to anat QC: {int(m['anat_qc_matched'].sum())}  "
          f"|  unmatched: {int((~m['anat_qc_matched']).sum())}")
    print(f"did NOT pass anat QC (QC_seg>1): {len(fail)}")
    print(f"  QC_seg: {fail['QC_seg'].value_counts().sort_index().to_dict()}")
    print(f"  cohort: {fail['cohort'].value_counts().to_dict()}")
    print(f"  group : {fail['population_group'].value_counts().to_dict()}")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
