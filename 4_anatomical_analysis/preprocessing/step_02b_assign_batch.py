"""
Step 2b — Build the ComBat batch column for the NDD-focused pipeline.

Merges with the curated clinical TSVs using cohort-specific key strategies,
then constructs the `batch` column.

    cohort               merge key (FS <-> CLN)            batch label
    ----------------------------------------------------------------------
    LEAP_W1/W2/W3        MRI_ID <-> MRI_ID                 mri_machine (e.g. KCL)
    INOVAND_T1/T2        MRI_ID <-> strip_sub(MRI_ID)       mri_machine (Ingenia/Intera/Ingenia 3T)
    INFOR                MRI_ID <-> strip_sub(MRI_ID)       mri_machine (Ingenia 3T)
    PIP                  MRI_ID <-> MRI_ID                  mri_machine (numeric site code)
    CANDY_PIP_W1/W2      MRI_ID <-> MRI_ID                  mri_machine (numeric site code)

Unmatched subjects fall back to a per-cohort label:
    INOVAND_T1, INOVAND_T2  ->  "INOVAND_unmatched"
    INFOR                   ->  "INFOR_unmatched"
    LEAP_W1/2/3             ->  "LEAP_unmatched"
    PIP / CANDY_PIP_W1/W2   ->  dropped (no mri_machine)
"""
import re
import pandas as pd

from config import load_clinical, BATCH_COL


def _strip_sub(s: object) -> str:
    """sub-0004 -> 4   |   124 -> 124   |   124.0 -> 124"""
    s = str(s)
    if s.endswith(".0"):
        s = s[:-2]
    m = re.match(r"^sub-?0*(\d+)$", s)
    if m:
        return m.group(1)
    return s.lstrip("0") or s


def _load_clinical() -> pd.DataFrame:
    cln = load_clinical(low_memory=False, dtype=str)[
        ["ID", "MRI_ID", "cohort", "mri_machine", "site_id", "age_yrs", "Sex"]
    ]
    cln["MRI_ID_n"] = cln["MRI_ID"].map(_strip_sub)
    cln["age_yrs"]  = pd.to_numeric(cln["age_yrs"], errors="coerce")
    cln["Sex"]      = pd.to_numeric(cln["Sex"],     errors="coerce")
    return cln


def _clinical_subset_for_cohort(cohort: str, cln: pd.DataFrame) -> pd.DataFrame:
    """Return clinical rows that should be used to enrich rows of this FS
    cohort. LEAP_W1 <-> LEAP_T1 etc. for per-wave age matching.
    """
    if cohort.startswith("LEAP_W"):
        wave = cohort.replace("LEAP_W", "LEAP_T")
        sub = cln[cln["cohort"] == wave]
    elif cohort.startswith("INOVAND"):
        sub = cln[cln["cohort"].str.startswith("INOVAND", na=False)]
    elif cohort == "INFOR":
        sub = cln[cln["cohort"] == "INFOR"]
    elif cohort == "PIP":
        sub = cln[cln["cohort"] == "PIP"]
    elif cohort.startswith("CANDY_PIP"):
        sub = cln[cln["cohort"] == cohort]
    else:
        return cln.iloc[0:0]

    # De-dupe per MRI_ID_n: prefer rows with non-null age_yrs and mri_machine
    # so the kept row carries the maximum useful info.
    sub = sub.assign(
        _has_age=sub["age_yrs"].notna().astype(int),
        _has_machine=sub["mri_machine"].notna().astype(int),
    )
    sub = sub.sort_values(["_has_machine", "_has_age"], ascending=False)
    sub = sub.drop_duplicates("MRI_ID_n", keep="first")
    return sub.drop(columns=["_has_age", "_has_machine"])


def assign_batch(df: pd.DataFrame) -> pd.DataFrame:
    cln = _load_clinical()

    df = df.copy()
    df["MRI_ID"] = df["MRI_ID"].astype(str).str.replace(r"\.0$", "", regex=True)
    # FS file already has an empty `mri_machine` column — drop it before merging
    if "mri_machine" in df.columns:
        df = df.drop(columns="mri_machine")

    df["mri_machine"] = pd.NA

    # The FS file's demographics are unreliable for INOVAND/INFOR (placeholder
    # values: age=9.0 / Sex=1.0 for everyone). Clinical is authoritative, so
    # we OVERRIDE — not just fill — age_yrs and Sex from clinical whenever
    # the merge succeeds.
    df["age_yrs"] = pd.to_numeric(df["age_yrs"], errors="coerce")
    df["Sex"]     = pd.to_numeric(df["Sex"],     errors="coerce")

    n_age_set = n_sex_set = 0
    for cohort, sub in df.groupby("cohort", sort=False):
        cln_sub = _clinical_subset_for_cohort(cohort, cln)
        if cln_sub.empty:
            print(f"[02b] {cohort:18s} no clinical rule — skip")
            continue
        cln_idx = cln_sub.set_index("MRI_ID_n")
        keys = sub["MRI_ID"]

        machines    = keys.map(cln_idx["mri_machine"])
        cln_age     = keys.map(cln_idx["age_yrs"])
        cln_sex     = keys.map(cln_idx["Sex"])

        # Override (not fill) when clinical has a value
        age_ovr = cln_age.notna()
        sex_ovr = cln_sex.notna()
        df.loc[sub.index, "mri_machine"] = machines.values
        df.loc[sub.index[age_ovr], "age_yrs"] = cln_age[age_ovr].values
        df.loc[sub.index[sex_ovr], "Sex"]     = cln_sex[sex_ovr].values
        n_age_set += int(age_ovr.sum())
        n_sex_set += int(sex_ovr.sum())

        print(f"[02b] {cohort:18s} matched {machines.notna().sum():4d}/{len(sub):4d}  "
              f"(age overridden from CLN={int(age_ovr.sum())}, "
              f"sex={int(sex_ovr.sum())})")

    # Recompute age^2 from the (possibly new) age_yrs
    df["age2"] = df["age_yrs"] ** 2
    print(f"[02b] demographics overridden from clinical: "
          f"age={n_age_set}, sex={n_sex_set}, age^2 recomputed")

    # Drop unmatched subjects: no clinical entry -> no scanner -> no batch.
    n_total = len(df)
    unmatched = df["mri_machine"].isna()
    if unmatched.any():
        per_cohort = df.loc[unmatched, "cohort"].value_counts().to_dict()
        print(f"[02b] dropping {int(unmatched.sum()):,} unmatched subjects  "
              f"(no clinical mri_machine):  {per_cohort}")
        df = df.loc[~unmatched].copy()

    df[BATCH_COL] = df["mri_machine"].astype("string")
    n_batch = df[BATCH_COL].nunique()
    print(f"[02b] kept {len(df):,}/{n_total:,}  ->  {n_batch} unique batches")
    return df


if __name__ == "__main__":
    from step_01_load    import load_dataframe
    from step_02_filter_qc import restrict_to_ndd, filter_qc
    df = restrict_to_ndd(load_dataframe())
    df = filter_qc(df)
    df = assign_batch(df)
    print("\nBatch label distribution:")
    print(df[BATCH_COL].value_counts())
