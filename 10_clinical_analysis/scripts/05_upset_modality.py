#!/usr/bin/env python3
# =============================================================================
# 05 - UpSet plot: data-modality availability by phenotype (curated cohort)
# =============================================================================
# Port of eeg_mri-pipeline/.../demographics_cohort/fig_upset_plot.py, on the
# CURATED cohort: WGS x MRI x EEG availability intersections, with stacked bars
# coloured by phenotype (NT / Relatives / IDD / Autism).
#
#   phenotype + MRI/EEG : curated roster (CLUSTERS_CURATED_FILE); MRI/EEG derived
#                         from MRI_ID / EEG_ID being populated.
#   WGS                 : joined by ID from the genetic roster's wgs_score.
#
# Output: figures/upset_modality_by_phenotype.pdf (+ .png)
# =============================================================================

import os
import sys
import pandas as pd
import matplotlib.pyplot as plt
from upsetplot import UpSet

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _config import CLUSTERS_FILE, FIGURES_DIR

PALETTE = {"Autism": "#8991FA", "Relatives": "#9AD5D3",
           "IDD": "#D8A4CB", "NT": "#C1C2BC"}
ORDER = ["NT", "Relatives", "IDD", "Autism"]


def load():
    """Entire cohort from the single roster where MRI / EEG / WGS availability
    are all recorded on the SAME ID (df_clusters_complete_kmeans.csv). Deriving
    MRI/EEG from the curated MRI_ID/EEG_ID and joining WGS from v4 loses the
    intersections, because curated LEAP IDs (numeric barcodes) don't match the
    v4 / roster C0733-* IDs. wgs_score here (n=2535) matches gnomAD-v4 (n=2531)."""
    df = pd.read_csv(CLUSTERS_FILE, low_memory=False).drop_duplicates("ID")
    df = df[~df.get("gene_first", False).fillna(False).astype(bool)]  # drop genetics-first ascertained
    df["MRI"] = df["MRI"].fillna(False).astype(bool)
    df["EEG"] = df["EEG"].fillna(False).astype(bool)
    df["WGS"] = df["wgs_score"].fillna(False).astype(bool)
    df["PopulationS1"] = df["PopulationS1"].replace({"TD": "NT", "ID": "IDD"})
    df = df[df["PopulationS1"].isin(ORDER)]
    print(f"Entire cohort: {len(df)}  (MRI={df['MRI'].sum()}, "
          f"EEG={df['EEG'].sum()}, WGS={df['WGS'].sum()})")
    return df


def main():
    df = load()
    df["PopulationS1"] = pd.Categorical(df["PopulationS1"], categories=ORDER, ordered=True)
    print(df.groupby("PopulationS1", observed=True)[["WGS", "MRI", "EEG"]].sum())

    # Boolean multi-index for UpSet (WGS x MRI x EEG).
    data = (df.set_index(df["WGS"] == 1)
              .set_index(df["MRI"] == 1, append=True)
              .set_index(df["EEG"] == 1, append=True))

    fig = plt.figure(figsize=(10, 15))
    upset = UpSet(data, intersection_plot_elements=0, element_size=50,
                  show_counts=True, show_percentages="{:.1%}",
                  orientation="horizontal", sort_by="cardinality")
    upset.add_stacked_bars(by="PopulationS1", colors=PALETTE,
                           title="Count by phenotype", elements=10)
    upset.plot(fig=fig)
    plt.suptitle(" ")
    plt.show()   # size/save interactively


if __name__ == "__main__":
    main()
