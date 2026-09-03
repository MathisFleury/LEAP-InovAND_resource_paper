#!/usr/bin/env python3
# =============================================================================
# 08 - Cohort size, before (original submission) vs now (current curated)
# =============================================================================
# Built for the reviewer-response cover letter's "Addition of new
# participants" paragraph (2026-09-01/02). Answers two questions per group
# (Autism / NT / IDD / Relatives):
#   1. How many individuals were there at the time of the original
#      submission, vs now?
#   2. Of the ORIGINAL individuals, how many can still be found by ID in the
#      CURRENT curated data at all (regardless of group), and how many
#      cannot be found anywhere -- i.e. real attrition, not just a group
#      relabelling.
#
# "Before" = the frozen, paper-aligned snapshot used for the original
# submission (single file, one row per person, cohort tag not always
# present -- see the NO_COHORT_TAG check below).
# "After"  = the current curated regime: one curated TSV per cohort
# (LEAP has 3 waves; INFOR is fully contained within INOVAND's IDs -- see
# the INFOR_SUBSET_OF_INOVAND check below, so it is NOT added a second time).
#
# All paths below are absolute and spelled out on purpose (per project
# convention, this script is meant to be read top-to-bottom by a human, not
# just run) -- nothing is hidden behind a shared config import.
#
# Outputs (10_clinical_analysis/outputs/tables/):
#   cohort_size_before_after.csv     one row per group: before / after / delta
#   cohort_size_missing_ids.csv      one row per original individual not
#                                     findable (by ID) anywhere in current data
#   cohort_size_new_ids.csv          the opposite: one row per current
#                                     individual not findable (by ID) in the
#                                     original data -- i.e. genuinely new
#                                     participants added since the original
#                                     submission
# =============================================================================
import os

import pandas as pd

# --- "Before": original-submission frozen snapshot ---------------------------
# One person per row. Group column is `PopulationS1`. Cohort tag (`cohort`)
# is present for LEAP/INOVAND but is MISSING for a large fraction of rows
# (see NO_COHORT_TAG_N below) -- those rows still count in the "before" total,
# they just cannot be cross-checked against a specific current-cohort file.
FROZEN_INDIVIDUALS_METRICS = (
    "/Users/mfleury/POSTDOC/LIBRAIRY/LEAP-InovAND_resource/script_zakaria/individuals_metrics_pull_9juillet2025.tsv"
)

# --- "After": current curated per-cohort snapshots ---------------------------
# LEAP has three waves (T1/T2/T3); a person can appear in more than one, so
# we deduplicate by ID across all three before counting.
CURRENT_LEAP_WAVES = [
    "/Volumes/Imaging5/EEG_MRI-MF/LEAP/_clinical_data/curated/LEAP_clinical_curated_t1.tsv",
    "/Volumes/Imaging5/EEG_MRI-MF/LEAP/_clinical_data/curated/LEAP_clinical_curated_t2.tsv",
    "/Volumes/Imaging5/EEG_MRI-MF/LEAP/_clinical_data/curated/LEAP_clinical_curated_t3.tsv",
]
CURRENT_INOVAND = "/Volumes/Imaging5/EEG_MRI-MF/INOVAND/_clinical_data/curated/INOVAND_clinical_curated.tsv"
# INFOR's IDs are a COMPLETE subset of INOVAND's (verified below at runtime) --
# every INFOR individual already has a row in the INOVAND file, so INFOR is
# read only for that verification, never added into the "after" total.
# Naively summing all three files (LEAP + INOVAND + INFOR) double-counts the
# INFOR/INOVAND overlap -- this was an actual bug in an earlier hand-computed
# version of this comparison; this script exists so that mistake can't recur.
CURRENT_INFOR = "/Volumes/Imaging5/EEG_MRI-MF/INFOR/_clinical_data/curated/INFOR_clinical_curated.tsv"

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_SECTION_DIR = os.path.dirname(_SCRIPT_DIR)
TABLES_DIR = os.path.join(_SECTION_DIR, "outputs", "tables")

GROUPS = ["Autism", "NT", "IDD", "Relatives"]
# The frozen file spells NT "TD" and IDD "ID"; normalise so both regimes use
# the same labels as the current curated `population_group` column.
FROZEN_GROUP_MAP = {"TD": "NT", "ID": "IDD"}


def _normalize_id(raw):
    """"sub-0001" -> "1"; "894.0" -> "894"; barcode-style IDs (contain "-"
    beyond a leading "sub-", e.g. "C0733-011-190-001") pass through
    unchanged. Same convention as curated_clinical.py's `_norm_mri` --
    without this, a person recorded as an MRI-style "sub-XXXX" ID in one
    snapshot and a plain numeric ID in the other looks like two different
    people and shows up as a false "missing" hit."""
    s = str(raw).strip()
    if s.lower() in ("", "nan"):
        return s
    if s.startswith("sub-"):
        s = s[4:]
    try:
        return str(int(float(s)))
    except (TypeError, ValueError):
        return s


def load_before():
    print(f"[before] reading {FROZEN_INDIVIDUALS_METRICS}")
    df = pd.read_csv(FROZEN_INDIVIDUALS_METRICS, sep="\t", low_memory=False)
    df["ID_raw"] = df["ID"].astype(str)
    df["ID"] = df["ID_raw"].map(_normalize_id)
    df["group"] = df["PopulationS1"].replace(FROZEN_GROUP_MAP)
    n_no_cohort = df["cohort"].isna().sum()
    n_renamed = (df["ID"] != df["ID_raw"]).sum()
    print(f"[before] {len(df)} rows total; {n_no_cohort} have no cohort tag "
          f"(cannot be cross-checked against a specific current-cohort file below, "
          f"but still count in the totals); {n_renamed} IDs were normalised "
          f"(e.g. 'sub-0001' -> '1')")
    return df[["ID", "ID_raw", "group", "cohort"]]


def load_after():
    leap_frames = []
    for wave_path in CURRENT_LEAP_WAVES:
        print(f"[after]  reading {wave_path}")
        d = pd.read_csv(wave_path, sep="\t", low_memory=False, usecols=["ID", "population_group"])
        d["ID"] = d["ID"].astype(str).map(_normalize_id)
        leap_frames.append(d)
    leap = pd.concat(leap_frames, ignore_index=True).drop_duplicates(subset=["ID"])
    leap["cohort"] = "LEAP"
    print(f"[after]  LEAP, deduplicated across 3 waves: {len(leap)} unique IDs")

    print(f"[after]  reading {CURRENT_INOVAND}")
    inovand = pd.read_csv(CURRENT_INOVAND, sep="\t", low_memory=False, usecols=["ID", "population_group"])
    inovand["ID"] = inovand["ID"].astype(str).map(_normalize_id)
    inovand = inovand.drop_duplicates(subset=["ID"])
    inovand["cohort"] = "INOVAND"
    print(f"[after]  INOVAND: {len(inovand)} unique IDs")

    print(f"[after]  reading {CURRENT_INFOR}  (verification only, see header note)")
    infor = pd.read_csv(CURRENT_INFOR, sep="\t", low_memory=False, usecols=["ID"])
    infor["ID"] = infor["ID"].astype(str).map(_normalize_id)
    n_infor_only = (~infor["ID"].isin(inovand["ID"])).sum()
    if n_infor_only:
        print(f"[after]  WARNING: {n_infor_only} INFOR IDs are NOT in INOVAND -- "
              f"the 'INFOR is a subset of INOVAND' assumption no longer holds, "
              f"they are being added separately")
        infor_extra = infor.loc[~infor["ID"].isin(inovand["ID"])].copy()
        infor_extra["population_group"] = pd.NA
        infor_extra["cohort"] = "INFOR"
    else:
        print(f"[after]  confirmed: all {len(infor)} INFOR IDs already present in INOVAND -- "
              f"not added again")
        infor_extra = infor.iloc[0:0]

    df = pd.concat([leap, inovand, infor_extra], ignore_index=True)
    df = df.rename(columns={"population_group": "group"})
    return df[["ID", "group", "cohort"]]


def main():
    os.makedirs(TABLES_DIR, exist_ok=True)
    before = load_before()
    after = load_after()

    before_ids = set(before["ID"])
    after_ids = set(after["ID"])

    print("\n" + "=" * 70)
    print("COHORT SIZE: before (original submission) vs after (current curated)")
    print("=" * 70)
    rows = []
    for g in GROUPS:
        n_before = int((before["group"] == g).sum())
        n_after = int((after["group"] == g).sum())
        rows.append({"group": g, "before": n_before, "after": n_after, "delta": n_after - n_before})
    rows.append({"group": "TOTAL", "before": len(before), "after": len(after),
                 "delta": len(after) - len(before)})
    summary = pd.DataFrame(rows)
    print(summary.to_string(index=False))

    summary_path = os.path.join(TABLES_DIR, "cohort_size_before_after.csv")
    summary.to_csv(summary_path, index=False)

    # --- ID-level attrition: original individuals not findable at all now ---
    print("\n" + "=" * 70)
    print("ATTRITION: original individuals not found by ID anywhere in current data")
    print("=" * 70)
    missing = before[~before["ID"].isin(after_ids)].copy()
    for g in GROUPS:
        n_g = int((before["group"] == g).sum())
        n_missing_g = int((missing["group"] == g).sum())
        tag_counts = missing.loc[missing["group"] == g, "cohort"].value_counts(dropna=False).to_dict()
        print(f"  {g:10s}: {n_missing_g:4d} / {n_g:4d} missing  (by original cohort tag: {tag_counts})")

    missing_path = os.path.join(TABLES_DIR, "cohort_size_missing_ids.csv")
    missing.to_csv(missing_path, index=False)

    # --- ID-level growth: current individuals not present originally --------
    # The opposite of the attrition check above: who is genuinely new.
    print("\n" + "=" * 70)
    print("NEW PARTICIPANTS: current individuals not found by ID in the original data")
    print("=" * 70)
    new = after[~after["ID"].isin(before_ids)].copy()
    for g in GROUPS:
        n_g = int((after["group"] == g).sum())
        n_new_g = int((new["group"] == g).sum())
        tag_counts = new.loc[new["group"] == g, "cohort"].value_counts(dropna=False).to_dict()
        print(f"  {g:10s}: {n_new_g:4d} / {n_g:4d} new  (by current cohort tag: {tag_counts})")

    new_path = os.path.join(TABLES_DIR, "cohort_size_new_ids.csv")
    new.to_csv(new_path, index=False)

    print(f"\nSaved: {summary_path}\n       {missing_path}\n       {new_path}")
    print("\nNOTE: 'missing' means not found by ID in the current curated LEAP/INOVAND/INFOR")
    print("files; 'new' is the reverse, not found by ID in the original frozen snapshot.")
    print("Neither establishes WHY on its own -- consent withdrawal, a later QC exclusion,")
    print("an ID renumbering between pulls, a genuine data-migration gap, and real new")
    print("recruitment are all possible and were not distinguished here. See")
    print("cohort_size_missing_ids.csv / cohort_size_new_ids.csv for the specific IDs to")
    print("trace against the source database.")


def _selftest():
    """Sanity check the group-merge logic on a tiny synthetic example."""
    d = pd.DataFrame({"ID": ["1", "2", "3"], "group": ["Autism", "TD", "Relatives"]})
    d["group"] = d["group"].replace(FROZEN_GROUP_MAP)
    assert list(d["group"]) == ["Autism", "NT", "Relatives"]


if __name__ == "__main__":
    _selftest()
    main()
