"""
Step 1 — Load LEAP and INOVAND XCP-D subject tables and merge framewise
displacement QC.

Input:
    df_xcp_subjects_LEAP.csv     per-scan metadata + connectivity path (LEAP)
    df_xcp_subjects_INOVAND.csv  per-scan metadata + connectivity path (INOVAND)
    qc_xcp_d_LEAP.csv            (subject, session, run, mean_fd, std_fd)
    qc_xcp_d_INOVAND.csv         (subject, session, run, mean_fd, std_fd)

Output:
    Single DataFrame with one row per (ID, session, run) carrying the
    connectivity matrix path, demographics, cohort, machine info, and
    mean_fd. INOVAND's ID column is rebuilt from ID_genet so the final
    keying matches the structural pipeline + clinical tables.
"""
import pandas as pd

from config import (
    XCP_D_SUBJECTS_LEAP, XCP_D_SUBJECTS_INOVAND,
    QC_XCP_LEAP, QC_XCP_INOVAND,
    AGE_COL,
)


def _load_leap() -> pd.DataFrame:
    df = pd.read_csv(XCP_D_SUBJECTS_LEAP, low_memory=False)
    df["ID"] = df["ID"].astype(str)

    qc = pd.read_csv(QC_XCP_LEAP)
    qc["ID"] = qc["subject"].str[4:]              # sub-XYZ → XYZ
    qc = qc[["ID", "session", "run", "mean_fd"]]

    df = df.merge(qc, on=["ID", "session", "run"], how="inner")
    print(f"[01] LEAP    loaded  {len(df):,} rows after QC merge")
    return df


def _load_inovand() -> pd.DataFrame:
    df = pd.read_csv(XCP_D_SUBJECTS_INOVAND, low_memory=False)
    df["ID"] = df["ID"].astype(str)

    # ID_genet is the cross-cohort key used by the clinical tables; promote
    # it to `ID`. Fall back to the BIDS MRI_ID when ID_genet is missing.
    df = df.drop(columns=["ID"]).rename(columns={"ID_genet": "ID"})

    qc = pd.read_csv(QC_XCP_INOVAND)
    qc = qc.rename(columns={"subject": "MRI_ID"})
    qc = qc[["MRI_ID", "session", "run", "mean_fd"]]

    df = df.merge(qc, on=["MRI_ID", "session", "run"], how="inner")
    df["ID"] = df["ID"].fillna(df["MRI_ID"])
    print(f"[01] INOVAND loaded  {len(df):,} rows after QC merge")
    return df


def load_dataframe() -> pd.DataFrame:
    """Load LEAP + INOVAND subject tables, concatenate, drop NaN ages, then
    attach the curated clinical diagnosis (replaces reliance on the old
    `control_status` column)."""
    df_leap    = _load_leap()
    df_inovand = _load_inovand()
    df = pd.concat([df_leap, df_inovand], ignore_index=True, sort=False)
    n0 = len(df)
    df = df.dropna(subset=[AGE_COL]).copy()
    print(f"[01] concat  →  {len(df):,} rows  (dropped {n0 - len(df):,} with NaN {AGE_COL})")

    # Curated clinical merge (PopulationS1_curated + clinical_matched). Same
    # canonical-key approach as 4_anatomical_analysis.
    from clinical_merge import merge_curated_clinical
    df = merge_curated_clinical(df)
    return df


if __name__ == "__main__":
    df = load_dataframe()
    print(df[["ID", "cohort", "session", "run", "mean_fd", AGE_COL, "Sex"]].head())
