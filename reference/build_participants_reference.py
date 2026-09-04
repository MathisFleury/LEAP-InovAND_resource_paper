#!/usr/bin/env python3
"""
Build a participant-level REFERENCE table for external reuse (curated base).

One row per individual PER TIMEPOINT (participants AND relatives). LEAP T1 and
LEAP T2 are kept as SEPARATE rows because their clinical values and their
MRI/EEG data differ by visit; WGS is per-person (identical across a subject's
timepoints). INOVAND is a single timepoint here.

Availability flags (MRI is BEFORE QC = acquired/processed, not QC-passed):
    has_MRI        any MRI (anatomical or functional) at that timepoint
    has_MRI_ANAT   anatomical MRI  (FreeSurfer processed) at that timepoint
    has_MRI_FUNC   functional MRI  (xcp_d fMRI)          at that timepoint
    has_EEG        EEG data available (subject-level; the EEG table has no
                   timepoint, so both LEAP waves get the same value)
    has_WGS        whole-genome sequencing available (per-person)

Timepoint matching:
  - anat roster distinguishes visits by its `cohort` (LEAP_T1 / LEAP_T2 / ...)
  - func roster distinguishes visits by `session` (ses-01 -> T1, ses-02 -> T2)
Ids differ by cohort (LEAP barcodes vs INOVAND sub-ids), so a subject matches a
roster if ANY of its identifiers (ID or MRI_ID, normalised) is in that roster.

Run:  python3.11 build_participants_reference.py
Output: reference/participants_reference_curated.tsv
"""

import pandas as pd

# ---------------------------------------------------------------------------
# PATHS (hard-coded on purpose — self-contained)
# ---------------------------------------------------------------------------
IMG5 = "/Volumes/Imaging5/EEG_MRI-MF"
LIB = "/Users/mfleury/POSTDOC/LIBRAIRY"

# curated clinical base, per cohort/timepoint (participants + relatives)
LEAP_T1_CLINICAL = f"{IMG5}/LEAP/_clinical_data/curated/LEAP_clinical_curated_t1.tsv"
LEAP_T2_CLINICAL = f"{IMG5}/LEAP/_clinical_data/curated/LEAP_clinical_curated_t2.tsv"
INOVAND_CLINICAL = f"{IMG5}/INOVAND/_clinical_data/curated/INOVAND_clinical_curated.tsv"

# modality rosters (MRI = before QC)
ANAT_ROSTER = f"{IMG5}/ALL/results/tabular/anat/freesurfer/8.1.0/output/curated_freesurfer_summary_all_cohorts.tsv"
# func roster = the CURRENT XCP-D v0.11 per-run connectivity (6-2). One combined
# file for both cohorts, with `cohort` and `session` (ses-01->T1, ses-02->T2).
# (Replaces the stale v0.x df_xcp_subjects_*.csv rosters.)
FUNC_ROSTER_V11 = f"{LIB}/LEAP-InovAND_resource/6_functional_analysis/concat/preprocessing/outputs/nogsr/df_conn_raw_nogsr.csv"
EEG_ROSTER = f"{LIB}/eeg_mri-pipeline/results/dataset_paper/dataframes/Alpha_peak_combined_corrected.csv"
WGS_ROSTER = f"{IMG5}/ALL/results/tabular/genetics/df_carrier_genelist_DEL_LOF_MISS_withalphamissense_LOEUF_gnomadv4.tsv"

OUTPUT = f"{LIB}/LEAP-InovAND_resource/reference/participants_reference_curated.tsv"

# basic clinical fields to keep (shared across the curated cohorts)
CLINICAL_FIELDS = [
    "ID", "MRI_ID", "EEG_ID", "cohort", "site_id", "mri_machine",
    "relation_to_proposant", "FID", "population", "population_group",
    "pheno_ASD", "pheno_ID", "pheno_ADHD", "Sex", "age_yrs",
    "total_IQ", "verbal_IQ", "performance_IQ", "IQ_type_corrected",
    "SRS_rawscore", "SRS_tscore", "RBS-R_total", "ssp_total",
    "ADI_social_interaction", "ADI_communication", "ADI_rrb",
    "ADOS_sa_css", "ADOS_rrb_css", "ADOS_total_css",
    "vabsabcabc_standard", "vabsdscoresc_dss", "vabsdscoresd_dss", "vabsdscoress_dss",
]


def normalise_id(value):
    """Normalise an id so it matches across tables: strip whitespace, a trailing
    '.0' (float-inferred ids) and a leading 'sub-' (INOVAND fMRI ids)."""
    s = str(value).strip()
    if s in ("", "nan", "None", "NaN"):
        return None
    if s.endswith(".0"):
        s = s[:-2]
    if s.startswith("sub-"):
        s = s[len("sub-"):]
    if s.isdigit():                       # BIDS zero-padding: "0003" (func) == "3" (clinical)
        s = s.lstrip("0") or "0"          # ponytail: numeric-only, LEAP barcodes untouched
    return s


def load_clinical(path, cohort, timepoint):
    """Read one curated clinical file; keep shared fields; tag cohort+timepoint."""
    df = pd.read_csv(path, sep="\t", low_memory=False)
    df["ID"] = df["ID"].astype(str)
    df = df.drop_duplicates(subset="ID", keep="first")   # guard against dup rows
    keep = [c for c in CLINICAL_FIELDS if c in df.columns]
    df = df[keep].copy()
    df["cohort"] = cohort          # clean label
    df["timepoint"] = timepoint    # T1 / T2 / W1
    return df


def id_set(path, sep, id_cols, row_filter=None):
    """Normalised identifier set from a roster, optionally row-filtered first."""
    df = pd.read_csv(path, sep=sep, low_memory=False)
    if row_filter is not None:
        df = df[row_filter(df)]
    ids = set()
    for col in id_cols:
        if col in df.columns:
            ids.update(normalise_id(v) for v in df[col])
    ids.discard(None)
    return ids


def build_rosters():
    """Timepoint-aware anat/func sets + flat EEG/WGS sets."""
    fs = pd.read_csv(ANAT_ROSTER, sep="\t", low_memory=False)
    fs["MRI_ID"] = fs["MRI_ID"].astype(str)
    fs["cohort"] = fs["cohort"].astype(str)

    def fs_ids(mask):
        return {normalise_id(v) for v in fs.loc[mask, "MRI_ID"]} - {None}

    anat = {
        "LEAP:T1": fs_ids(fs["cohort"] == "LEAP_T1"),
        "LEAP:T2": fs_ids(fs["cohort"] == "LEAP_T2"),
        "INOVAND:W1": fs_ids(fs["cohort"].str.startswith("INOVAND")),
    }
    func = {
        "LEAP:T1": id_set(FUNC_ROSTER_V11, ",", ["ID"],
                          lambda d: (d["cohort"] == "LEAP") & (d["session"].astype(str) == "ses-01")),
        "LEAP:T2": id_set(FUNC_ROSTER_V11, ",", ["ID"],
                          lambda d: (d["cohort"] == "LEAP") & (d["session"].astype(str) == "ses-02")),
        "INOVAND:W1": id_set(FUNC_ROSTER_V11, ",", ["ID"],
                             lambda d: d["cohort"] == "INOVAND"),
    }
    eeg = id_set(EEG_ROSTER, ",", ["ID"])          # subject-level (no timepoint)
    wgs = id_set(WGS_ROSTER, "\t", ["ID"])         # per-person
    return anat, func, eeg, wgs


def main():
    base = pd.concat([
        load_clinical(LEAP_T1_CLINICAL, "LEAP", "T1"),
        load_clinical(LEAP_T2_CLINICAL, "LEAP", "T2"),
        load_clinical(INOVAND_CLINICAL, "INOVAND", "W1"),
    ], ignore_index=True)
    print(f"base rows (subject x timepoint): {len(base)}")
    print(f"  by cohort/timepoint: {base.groupby(['cohort','timepoint']).size().to_dict()}")
    print(f"  relations: {base['relation_to_proposant'].value_counts(dropna=False).to_dict()}")

    anat, func, eeg, wgs = build_rosters()

    def flags(row):
        keys = {normalise_id(row["ID"]), normalise_id(row["MRI_ID"])} - {None}
        bucket = f"{row['cohort']}:{row['timepoint']}"
        return pd.Series({
            "has_MRI_ANAT": bool(keys & anat.get(bucket, set())),
            "has_MRI_FUNC": bool(keys & func.get(bucket, set())),
            "has_EEG": bool(keys & eeg),
            "has_WGS": bool(keys & wgs),
        })

    f = base.apply(flags, axis=1)
    base = pd.concat([base, f], axis=1)
    base["has_MRI"] = base["has_MRI_ANAT"] | base["has_MRI_FUNC"]

    # order: identity/clinical first, flags last
    flag_cols = ["has_MRI", "has_MRI_ANAT", "has_MRI_FUNC", "has_EEG", "has_WGS"]
    front = ["ID", "timepoint", "cohort"]
    base = base[front + [c for c in base.columns if c not in front + flag_cols] + flag_cols]

    print("\navailability (TRUE counts):")
    for c in flag_cols:
        print(f"  {c:14s} {int(base[c].sum())}")

    base.to_csv(OUTPUT, sep="\t", index=False)
    print(f"\nWrote {len(base)} rows x {base.shape[1]} cols -> {OUTPUT}")


if __name__ == "__main__":
    main()
