"""
Curated clinical merge for the fMRI preprocessing pipeline.

Replaces the old `control_status`-based grouping with the up-to-date curated
diagnosis (`population_group`), joining each cohort on the key that actually
matches the fMRI identifiers:

  LEAP    -> int(ID)                    (numeric subject id)
  INOVAND -> MRI_ID string ("sub-XXXX") via the *merged* curated TSV

Why the merged INOVAND file: the fMRI INOVAND MRI_ID is the BIDS `sub-XXXX`
id, which the merged curated TSV uses verbatim (420/420 match). The non-merged
TSV (used by 4_anatomical for the FreeSurfer `894`-scheme) does NOT match the
BIDS ids — that earlier gave only ~180/420 coincidental overlaps.

LEAP is taken from curated_clinical.load_all_cohorts (T1->T2->T3 backfill,
participant filter, population_group passthrough); only the INOVAND join key /
source file differs here.
"""
import sys
from pathlib import Path
import numpy as np
import pandas as pd

_CUR = Path("/Volumes/Imaging5/EEG_MRI-MF")
LEAP_WAVES = [_CUR / "LEAP/_clinical_data/curated" / f"LEAP_clinical_curated_{w}.tsv"
              for w in ("t1", "t2", "t3")]
INOVAND_MERGED = _CUR / "INOVAND/_clinical_data/curated/INOVAND_clinical_curated_merged.tsv"
_CARRY = ("population_group", "site_id")


def _leap_key(v):
    try:
        return "L" + str(int(float(v)))
    except (TypeError, ValueError):
        return None


def _fmri_key(row):
    cohort = str(row.get("cohort", ""))
    if cohort.startswith("LEAP"):
        return _leap_key(row.get("ID"))
    if cohort.startswith("INOVAND"):
        mri = str(row.get("MRI_ID", "")).strip()
        return "I" + mri if mri and mri.lower() != "nan" else None
    return None


def _curated_keyed() -> pd.DataFrame:
    """One row per canonical key carrying population_group (+ site_id).

    No participant filter — every subject with a scan gets a population_group
    (relatives tagged 'Relatives'); the Autism/NT restriction downstream excludes
    them. LEAP keyed by int(ID) over T1->T2->T3 (earliest wins); INOVAND keyed by
    the BIDS MRI_ID string via the *merged* curated TSV (matches fMRI sub-XXXX).
    """
    # LEAP: T1 -> T2 -> T3, earliest non-duplicate wins.
    leap_frames = []
    for order, path in enumerate(LEAP_WAVES):
        if path.exists():
            d = pd.read_csv(path, sep="\t", low_memory=False)
            d["_ord"] = order
            leap_frames.append(d)
    leap = pd.concat(leap_frames, ignore_index=True)
    leap["_ckey"] = leap["ID"].map(_leap_key)
    leap = (leap.dropna(subset=["_ckey"]).sort_values("_ord")
                .drop_duplicates("_ckey", keep="first"))

    # INOVAND: merged TSV keyed by BIDS MRI_ID (string).
    ino = pd.read_csv(INOVAND_MERGED, sep="\t", low_memory=False)
    ino = ino[ino["MRI_ID"].notna()].copy()
    ino["_ckey"] = "I" + ino["MRI_ID"].astype(str).str.strip()

    keep = ["_ckey", *_CARRY]
    frames = []
    for d in (leap, ino):
        for c in _CARRY:
            if c not in d.columns:
                d[c] = np.nan
        frames.append(d[keep])
    clin = pd.concat(frames, ignore_index=True).dropna(subset=["_ckey"])
    # Prefer a row with a known population_group when a key repeats.
    clin = (clin.sort_values("population_group", na_position="last")
                .drop_duplicates("_ckey", keep="first"))
    return clin


def merge_curated_clinical(df: pd.DataFrame, **_) -> pd.DataFrame:
    df = df.copy()
    df["_ckey"] = df.apply(_fmri_key, axis=1)
    clin = _curated_keyed()
    merged = df.merge(clin, on="_ckey", how="left", suffixes=("", "_curated"))
    merged["clinical_matched"] = merged["population_group"].notna()

    print("[clinical] curated-clinical merge coverage:")
    for coh, sub in merged.groupby("cohort"):
        print(f"    {coh:8s}: matched {int(sub['clinical_matched'].sum()):4d} / {len(sub):4d}")
    vc = merged.loc[merged.clinical_matched, "population_group"].value_counts().to_dict()
    print(f"    curated population_group: {vc}")
    return merged.drop(columns=["_ckey"])


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from step_01_load import load_dataframe
    m = load_dataframe()  # already runs merge_curated_clinical internally
    print(m[["ID", "cohort", "control_status", "population_group",
             "clinical_matched"]].head())
