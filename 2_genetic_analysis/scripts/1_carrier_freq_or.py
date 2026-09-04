#!/usr/bin/env python3
# =============================================================================
# 1 - Carrier Frequencies & Odds Ratios (gnomad v4, PAN, two constraint defs)
# =============================================================================
# Runs the headline PAN PopulationS1 constrained freq + OR analysis on the
# gnomad v4 pre-annotated carrier matrix (CARRIER_V4), for BOTH constraint
# definitions requested:
#   - "_contraint_"        (standard, pLI-filtered)
#   - "_contraint_wo_PLI_" (pLI filter dropped)
# Each is run for DEL+LoF and for DEL+LoF+Miss. Outputs -> figures/tables.
#
# Plot/OR functions are reused from _carrier_freq_or_shared.py. Unlike the
# legacy hg19/hg38 scripts, the v4 file already
# uses IDD/NT spelling and a boolean diag_genetic column, so no relabel is done.
# =============================================================================

import os
import sys
import pandas as pd
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _config import (
    FIGURES_DIR_V4, TABLES_DIR_V4, CARRIER_V4,
    PALETTE_FREQ, ORDER_FREQ,
    LABELS_CONSTRAINED, COLS_DELLOF_CONSTRAINED,
    LABELS_DELLOFMISS_CONSTRAINED, COLS_DELLOFMISS_CONSTRAINED,
    LABELS_DELLOF_CONSTRAINED_WOPLI, COLS_DELLOF_CONSTRAINED_WOPLI,
    LABELS_DELLOFMISS_CONSTRAINED_WOPLI, COLS_DELLOFMISS_CONSTRAINED_WOPLI,
)

# Reuse plotting/OR functions from script 3
from importlib.util import spec_from_file_location, module_from_spec
_script_dir = os.path.dirname(os.path.abspath(__file__))
_spec = spec_from_file_location("script02", os.path.join(_script_dir, "_carrier_freq_or_shared.py"))
_mod = module_from_spec(_spec)
_spec.loader.exec_module(_mod)
plot_carrier_frequencies = _mod.plot_carrier_frequencies
compute_and_plot_odds_ratios = _mod.compute_and_plot_odds_ratios

os.makedirs(FIGURES_DIR_V4, exist_ok=True)
os.makedirs(TABLES_DIR_V4, exist_ok=True)

# (variant, constraint-definition) -> columns, labels, output stem
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
    # NB: v4 already uses IDD/NT and ships diag_genetic — no relabel applied.

    # Build PAN PopulationS1 frequency grouping (same logic as scripts 3/5)
    def assign_pop_freq(row):
        if row["PopulationS1"] == "Autism":
            return "Autism"
        elif row["PopulationS1"] == "NT":
            return "NT"
        elif row["PopulationS1"] == "IDD":
            return "IDD"
        elif row["Population_undiagnosed1"] == "Undiagnosed siblings":
            return "Undiagnosed siblings"
        elif row["Population_undiagnosed1"] == "Undiagnosed parents":
            return "Undiagnosed parents"
        return np.nan

    df["PopulationS1_freq"] = df.apply(assign_pop_freq, axis=1)

    # All ancestries (PAN) + EUR-only
    subsets = [("PAN", df), ("EUR", df[df["EUR_ancestry"] == True])]

    for anc, df_anc in subsets:
        for variant, defn, cols, labels, stem in ANALYSES:
            cols = [c for c in cols if c != "diag_genetic"]  # no returnable-variants bar
            missing = [c for c in cols if c not in df_anc.columns]
            if missing:
                print(f"  Skipping {anc}_{stem}: missing columns {missing}")
                continue
            print(f"\n=== {anc} {variant} {defn} ===")
            plot_carrier_frequencies(
                df_anc, "PopulationS1_freq", cols,
                PALETTE_FREQ, ORDER_FREQ, labels,
                f"{anc}_{stem}_freq.pdf", ylim=(0, 1), figures_dir=FIGURES_DIR_V4,
            )
            compute_and_plot_odds_ratios(
                df_anc, "PopulationS1_freq", cols,
                PALETTE_FREQ, ORDER_FREQ, labels,
                f"{anc}_{stem}_or.pdf", tables_dir=TABLES_DIR_V4, figures_dir=FIGURES_DIR_V4,
            )

    print("\n=== Done. v4 figures/tables in:", FIGURES_DIR_V4, "===")


if __name__ == "__main__":
    main()
