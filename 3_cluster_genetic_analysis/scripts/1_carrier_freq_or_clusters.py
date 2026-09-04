#!/usr/bin/env python3
# =============================================================================
# 1 - Carrier Frequencies & Odds Ratios by Clinical Cluster (gnomad v4,
#     current curated k-means clustering)
# =============================================================================
# Cluster-level (C1/C2/C3 + NT + IDD) analogue of
# 2_genetic_analysis/1_carrier_freq_or.py: constrained carrier freq + OR on
# the gnomad v4 matrix, for both constraint definitions (_contraint_ and
# _contraint_wo_PLI_), DEL+LoF and DEL+LoF+Miss, run for PAN and EUR-only
# ancestries. No returnable-variants (diag_genetic) bar.
#
# Moved here from 2_genetic_analysis (which is now population-level only) so
# cluster-level analysis lives in the "cluster_X_analysis" section, matching
# the 4_anatomical_analysis/5_cluster_anatomical_analysis and
# 6_functional_analysis/7_cluster_functional_analysis convention. Every other
# clustering-method variant here (legacy/gmm/manual/spark_lof, hg19/hg38-v2)
# is superseded and lives in legacy/3_cluster_genetic_analysis/.
#
# (PGS-by-cluster (top25/bottom25 OR forest) is only in the legacy hg19
# cluster script; it doesn't use the gene list, so it's unchanged here.)
# Outputs -> figures/tables with CLUSTER in the name.
# =============================================================================

import os
import sys
import pandas as pd
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _config import (
    FIGURES_DIR_CURRENT, TABLES_DIR_CURRENT, CARRIER_V4, CLUSTERS_FILE,
    PALETTE_CLUSTERS, ORDER_CLUSTERS,
    LABELS_CONSTRAINED, COLS_DELLOF_CONSTRAINED,
    LABELS_DELLOFMISS_CONSTRAINED, COLS_DELLOFMISS_CONSTRAINED,
    LABELS_DELLOF_CONSTRAINED_WOPLI, COLS_DELLOF_CONSTRAINED_WOPLI,
    LABELS_DELLOFMISS_CONSTRAINED_WOPLI, COLS_DELLOFMISS_CONSTRAINED_WOPLI,
)

# Reuse plotting/OR functions from 2_genetic_analysis's shared engine
from importlib.util import spec_from_file_location, module_from_spec
_script_dir = os.path.dirname(os.path.abspath(__file__))
_shared_engine = os.path.normpath(os.path.join(
    _script_dir, "..", "..", "2_genetic_analysis", "scripts", "_carrier_freq_or_shared.py"))
_spec = spec_from_file_location("script02", _shared_engine)
_mod = module_from_spec(_spec)
_spec.loader.exec_module(_mod)
plot_carrier_frequencies = _mod.plot_carrier_frequencies
compute_and_plot_odds_ratios = _mod.compute_and_plot_odds_ratios

os.makedirs(FIGURES_DIR_CURRENT, exist_ok=True)
os.makedirs(TABLES_DIR_CURRENT, exist_ok=True)

# Same ANALYSES table as 1_carrier_freq_or.py
ANALYSES = [
    ("dellof",     "constraint",       COLS_DELLOF_CONSTRAINED,            LABELS_CONSTRAINED,                    "dellof_constraint"),
    ("dellof",     "constraint_woPLI", COLS_DELLOF_CONSTRAINED_WOPLI,      LABELS_DELLOF_CONSTRAINED_WOPLI,       "dellof_constraint_woPLI"),
    ("dellofmiss", "constraint",       COLS_DELLOFMISS_CONSTRAINED,        LABELS_DELLOFMISS_CONSTRAINED,         "dellofmiss_constraint"),
    ("dellofmiss", "constraint_woPLI", COLS_DELLOFMISS_CONSTRAINED_WOPLI,  LABELS_DELLOFMISS_CONSTRAINED_WOPLI,   "dellofmiss_constraint_woPLI"),
]


def main():
    print(f"Loading gnomad v4 carrier matrix: {CARRIER_V4}")
    df = pd.read_table(CARRIER_V4, low_memory=False)
    df = df[df["PopulationS1"] != "other"]
    df = df.drop_duplicates("ID")
    # NB: v4 already uses IDD/NT spelling — no relabel applied.

    if "Cluster" in df.columns:
        df = df.drop(columns=["Cluster"])

    print(f"Loading clusters from: {CLUSTERS_FILE}")
    clusters = pd.read_csv(CLUSTERS_FILE, low_memory=False)
    merge_cols = ["ID", "Cluster"]
    if "Relation_to_proposant" in clusters.columns:
        merge_cols.append("Relation_to_proposant")
    df = df.merge(clusters[merge_cols], how="left", on="ID")

    # C1/C2/C3 for Autism participants, plus NT and IDD participants (matches the legacy hg19 cluster script)
    def assign_cluster_nt(row):
        if row["PopulationS1"] == "NT" and row.get("Relation_to_proposant") == "participant":
            return "NT"
        elif row["PopulationS1"] == "IDD" and row.get("Relation_to_proposant") == "participant":
            return "IDD"
        elif pd.isna(row.get("Cluster")):
            return np.nan
        elif row["Cluster"] in ("C1", "C2", "C3") and row["PopulationS1"] == "Autism":
            return row["Cluster"]
        return np.nan

    df["Cluster_with_NT"] = df.apply(assign_cluster_nt, axis=1)

    subsets = [("PAN", df), ("EUR", df[df["EUR_ancestry"] == True])]

    for anc, df_anc in subsets:
        for variant, defn, cols, labels, stem in ANALYSES:
            cols = [c for c in cols if c != "diag_genetic"]  # no returnable-variants bar
            missing = [c for c in cols if c not in df_anc.columns]
            if missing:
                print(f"  Skipping {anc}_CLUSTER_{stem}: missing columns {missing}")
                continue
            print(f"\n=== {anc} CLUSTER {variant} {defn} ===")
            plot_carrier_frequencies(
                df_anc, "Cluster_with_NT", cols,
                PALETTE_CLUSTERS, ORDER_CLUSTERS, labels,
                f"{anc}_CLUSTER_{stem}_freq.pdf", ylim=(0, 1), figures_dir=FIGURES_DIR_CURRENT,
            )
            compute_and_plot_odds_ratios(
                df_anc, "Cluster_with_NT", cols,
                PALETTE_CLUSTERS, ORDER_CLUSTERS, labels,
                f"{anc}_CLUSTER_{stem}_or.pdf", tables_dir=TABLES_DIR_CURRENT, figures_dir=FIGURES_DIR_CURRENT,
            )

    print("\n=== Done. v4 cluster figures/tables in:", FIGURES_DIR_CURRENT, "===")


if __name__ == "__main__":
    main()
