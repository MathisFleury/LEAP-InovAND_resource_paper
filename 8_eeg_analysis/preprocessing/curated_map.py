"""
Curated clinical demographics + k-means Cluster, single-sourced for both EEG
sections.

Demographics (age_yrs, Sex, PopulationS1, Population1) come from the BROAD curated
per-cohort clinical TSVs — not the small clustering table — so EEG subjects that
were never clustered (no IQ+SRS) are still kept. The k-means Cluster is left-joined
from the curated cluster table (NaN for NT / unclustered subjects).

INOVAND `Subject_ID` (sub-XXXX) maps to the clinical ID via the INOVAND_EEG_bids
crosswalk (ID lookup only). LEAP `ID` already equals the clinical ID.
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from config import (LEAP_CLINICAL, INOVAND_CLINICAL, CURATED_CLUSTERS,  # noqa: E402
                    INOVAND_BIDS_MAP)

_FROZEN_DEMO = ["PopulationS1", "Population1", "age_yrs", "Sex", "Cluster", "control_status"]


def load_curated_demographics() -> pd.DataFrame:
    """ID, age_yrs, Sex, PopulationS1, Population1 from both curated clinical TSVs."""
    frames = []
    for f in (LEAP_CLINICAL, INOVAND_CLINICAL):
        d = pd.read_table(f, low_memory=False)
        d = d[["ID", "age_yrs", "Sex", "population_group", "population"]].copy()
        frames.append(d)
    demo = pd.concat(frames, ignore_index=True).rename(
        columns={"population_group": "PopulationS1", "population": "Population1"})
    demo["ID"] = demo["ID"].astype(str)
    demo = demo.dropna(subset=["PopulationS1"]).drop_duplicates("ID")
    return demo


def load_curated_clusters() -> pd.DataFrame:
    """ID, Cluster (curated k-means) — left-joined so NT/unclustered stay NaN."""
    c = pd.read_csv(CURATED_CLUSTERS, low_memory=False)[["ID", "Cluster"]].copy()
    c["ID"] = c["ID"].astype(str)
    return c


def _inovand_ids_to_clinical(df: pd.DataFrame) -> pd.DataFrame:
    """Resolve INOVAND Subject_ID (sub-XXXX) to clinical ID via the crosswalk."""
    df = df.copy()
    df["ID"] = df["ID"].astype(str)
    xwalk = (pd.read_csv(INOVAND_BIDS_MAP, low_memory=False)[["ID", "INOVAND_EEG_bids"]].dropna())
    sub_to_id = dict(zip(xwalk["INOVAND_EEG_bids"].astype(str), xwalk["ID"].astype(str)))
    ino = df["cohort"] == "INOVAND"
    df.loc[ino, "ID"] = df.loc[ino, "Subject_ID"].astype(str).map(sub_to_id)
    return df


def attach_curated_demographics(df: pd.DataFrame) -> pd.DataFrame:
    """df needs cohort + ID (LEAP clinical) + Subject_ID (INOVAND sub-XXXX)."""
    df = _inovand_ids_to_clinical(df)
    df = df.drop(columns=[c for c in _FROZEN_DEMO if c in df.columns])
    df = df.merge(load_curated_demographics(), on="ID", how="inner")
    return df.merge(load_curated_clusters(), on="ID", how="left")
