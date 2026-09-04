#!/usr/bin/env python3
"""
Rebuild the combined corrected alpha peak from RAW $IMG5 features (retake of
correct_alpha_peak_LEAP_INOVAND.py), curated + full-sample.

Flow:
  1. read raw LEAP + INOVAND Alpha_peak.csv (different separators / ID schemes)
  2. harmonise IDs to the curated clinical ID
       LEAP     : Subject_ID is already the clinical ID
       INOVAND  : Subject_ID -> sub-XXXX -> clinical ID via INOVAND_BIDS_MAP
  3. merge curated demographics + k-means Cluster (curated file only)
  4. regress age/age²/sex then z-score, both on the FULL sample
Writes config.ALPHA_PEAK_CORRECTED, consumed by 8_ and 9_ analysis.
The final z-score lives in column `alpha_peak_corrected` (downstream does NOT
re-correct); `alpha_peak` keeps the raw peak.
"""
import sys
from pathlib import Path

import pandas as pd

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))
from config import (LEAP_FEATURES, INOVAND_FEATURES,  # noqa: E402
                    INOVAND_BIDS_MAP, OUT_DIR, ALPHA_PEAK_CORRECTED)
from curated_map import load_curated_demographics, load_curated_clusters  # noqa: E402
from step_05_regress import regress  # noqa: E402
from step_06_zscore import zscore    # noqa: E402


def _read_alpha(path: Path, sep: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep=sep)
    df = df[[c for c in df.columns if not str(c).lower().startswith("unnamed")]]
    return df[["Subject_ID", "alpha_peak"]].copy()


def load_raw_alpha() -> pd.DataFrame:
    leap = _read_alpha(LEAP_FEATURES / "Alpha_peak.csv", ",")
    leap["ID"] = leap["Subject_ID"].astype(str).str.replace("sub-", "", regex=False)
    leap["cohort"] = "LEAP"

    ino = _read_alpha(INOVAND_FEATURES / "Alpha_peak.csv", ";")
    ino["bids"] = ino["Subject_ID"].apply(lambda x: f"sub-{int(x):04d}")
    xwalk = (pd.read_csv(INOVAND_BIDS_MAP, low_memory=False)[["ID", "INOVAND_EEG_bids"]]
             .dropna().rename(columns={"INOVAND_EEG_bids": "bids"}))
    xwalk["ID"] = xwalk["ID"].astype(str)
    ino = ino.merge(xwalk, on="bids", how="inner")
    ino["cohort"] = "INOVAND"

    cols = ["ID", "cohort", "alpha_peak"]
    combined = pd.concat([leap[cols], ino[cols]], ignore_index=True)
    print(f"[load] raw alpha: LEAP={len(leap)}, INOVAND mapped={len(ino)}, total={len(combined)}")
    return combined


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    alpha = load_raw_alpha()
    alpha["ID"] = alpha["ID"].astype(str)
    alpha = alpha[alpha["alpha_peak"].notna() & (alpha["alpha_peak"] != 0)]

    # demographics from the BROAD curated clinical TSVs; Cluster left-joined
    alpha = alpha.merge(load_curated_demographics(), on="ID", how="inner")
    alpha = alpha.merge(load_curated_clusters(), on="ID", how="left")
    alpha = alpha.dropna(subset=["age_yrs"])
    print(f"[merge] after curated demographics merge: n={len(alpha)}")

    # full-sample regression + z-score (see config); result in alpha_peak_corrected
    alpha["alpha_peak_corrected"] = alpha["alpha_peak"]
    alpha = regress(alpha, ["alpha_peak_corrected"])
    alpha = zscore(alpha, ["alpha_peak_corrected"])

    alpha.to_csv(ALPHA_PEAK_CORRECTED, index=False)
    print(f"Wrote {ALPHA_PEAK_CORRECTED}  n={len(alpha)}  cohort={alpha['cohort'].value_counts().to_dict()}")
    print("PopulationS1:", alpha["PopulationS1"].value_counts(dropna=False).to_dict())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
