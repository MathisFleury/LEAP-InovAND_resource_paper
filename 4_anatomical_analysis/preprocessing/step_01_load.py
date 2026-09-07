"""
Step 1 — Load the FreeSurfer 8.1 QC+clinical TSV into a DataFrame.

The `_ML_corrected` input has the per-scan FS measures and QC labels, but
its clinical merge dropped some cohorts (e.g. INFOR: 41 QC-pass rows with
all clinical fields NaN). We rescue those rows by filling `age_yrs` /
`Sex` from the curated per-cohort clinical TSVs only where the ML file is
NaN — preserving the per-session ages the ML file got right for
longitudinal cohorts (LEAP_W2/W3, INOVAND, PIP).
"""
import re
import pandas as pd
from config import INPUT_TSV, load_clinical


def _normalize_mri_id(s) -> str | None:
    if pd.isna(s):
        return None
    s = str(s).strip()
    if s.lower().startswith("sub-"):
        s = s[4:]
    if s.endswith(".0") and s[:-2].isdigit():
        s = s[:-2]
    if s.isdigit() and len(s) > 1:
        s = s.lstrip("0") or "0"
    return s


def _normalize_cohort(c) -> str:
    if pd.isna(c):
        return ""
    return re.sub(r"_(W|T)\d+$", "", str(c).strip())


def _fill_clinical_from_curated(df: pd.DataFrame) -> pd.DataFrame:
    """Fill NaN age_yrs / Sex from the curated clinical TSVs on (cohort, MRI_ID)."""
    cln = load_clinical(low_memory=False,
                        usecols=["MRI_ID", "cohort", "Sex", "age_yrs", "population_group"])
    cln = cln.dropna(subset=["MRI_ID"]).copy()
    cln["_mri"] = cln["MRI_ID"].astype(str).map(_normalize_mri_id)
    cln["_coh"] = cln["cohort"].map(_normalize_cohort)
    cln = cln.drop_duplicates(["_coh", "_mri"], keep="first")

    df = df.copy()
    df["_mri"] = df["MRI_ID"].astype(str).map(_normalize_mri_id)
    df["_coh"] = df["cohort"].map(_normalize_cohort)
    joined = df.merge(
        cln[["_coh", "_mri", "Sex", "age_yrs", "population_group"]],
        on=["_coh", "_mri"], how="left", suffixes=("", "_cln"),
    )

    filled_by_col = {}
    for col in ("age_yrs", "Sex", "population_group"):
        cln_col = f"{col}_cln"
        if cln_col not in joined.columns:
            continue
        before_nan = joined[col].isna().sum()
        joined[col] = joined[col].fillna(joined[cln_col])
        after_nan = joined[col].isna().sum()
        filled_by_col[col] = before_nan - after_nan
        joined = joined.drop(columns=[cln_col])

    joined = joined.drop(columns=["_mri", "_coh"])
    parts = ", ".join(f"{k}=+{v}" for k, v in filled_by_col.items())
    print(f"[01] clinical fallback from per-cohort curated TSVs: {parts}")

    # Recompute age2 wherever age_yrs is present but age2 is missing.
    # Rows recovered by the fallback above have age_yrs but the ML file's
    # age2 column was NaN — leaving NaN in the design matrix makes ComBat
    # singular.
    if "age2" in joined.columns:
        age = pd.to_numeric(joined["age_yrs"], errors="coerce")
        mask = joined["age2"].isna() & age.notna()
        n_fixed = int(mask.sum())
        if n_fixed:
            joined.loc[mask, "age2"] = (age[mask] ** 2).values
            print(f"[01] recomputed age2 from age_yrs for {n_fixed} row(s)")
    return joined


def load_dataframe() -> pd.DataFrame:
    df = pd.read_csv(INPUT_TSV, sep="\t", low_memory=False)
    print(f"[01] loaded  {INPUT_TSV.name}  ->  {len(df):,} rows x {df.shape[1]:,} cols")
    df = _fill_clinical_from_curated(df)
    return df


if __name__ == "__main__":
    df = load_dataframe()
    print(df[["ID", "cohort", "QC_seg", "age_yrs", "Sex"]].head())
