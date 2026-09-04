#!/usr/bin/env python3
"""
Assemble the gnomad v4 input table for the common(PGS) x rare(carrier) interaction
figure (3_interaction_common_rare.R).

v4 is the priority genetics regime (current priority regime): carriers come from the v4
matrix (CARRIER_V4, IDD/NT spelling) rather than the v2 columns baked into the
curated cluster df. Curated throughout — no frozen INDIVIDUALS_METRICS:
    v4 carrier matrix  x  curated k-means clusters (ID; brings Cluster +
      relation_to_proposant -> Relation_to_proposant)
      + ancestry-specific v4 PGS (barcode)

Two ancestry versions, per 15_v4_cluster_clusters.py's convention:
    PAN  pan-ancestry PGS, whole cohort
    EUR  European-only PGS, restricted to EUR_ancestry == True
Writes one CSV each. Renames EUR_ancestry -> EUR_ancestry_y for the R figure; PGS
trait names keep their "/" so R's read.csv make.names-mangles them to the "." form
the figure references.

Output: outputs/tables/df_v4_interaction_inputs_{PAN,EUR}.csv  (needs $IMG5 mounted).
Figures then land in outputs/figures/, flat and PAN_/EUR_-prefixed (3_interaction_common_rare.R).
"""
import os
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _config import CARRIER_V4, PGS_PAN_V4, PGS_EUR_V4, PGS_TRAITS, CLUSTERS_CURATED_KMEANS_FILE

_SECTION = Path(__file__).resolve().parent.parent
OUT_DIR = Path(os.environ.get("V4_INTERACTION_DIR", _SECTION / "outputs" / "tables"))
OUT_DIR.mkdir(parents=True, exist_ok=True)
PGS_FILES = {"PAN": PGS_PAN_V4, "EUR": PGS_EUR_V4}


def main() -> int:
    base = pd.read_table(CARRIER_V4, low_memory=False)
    base = base[base["PopulationS1"] != "other"].drop_duplicates("ID")

    clusters = pd.read_csv(CLUSTERS_CURATED_KMEANS_FILE, low_memory=False)
    clusters = clusters[["ID", "Cluster", "relation_to_proposant"]].rename(
        columns={"relation_to_proposant": "Relation_to_proposant"})
    base = base.merge(clusters, how="left", on="ID")

    for anc, pgs_file in PGS_FILES.items():
        df = base[base["EUR_ancestry"] == True] if anc == "EUR" else base  # noqa: E712
        praw = pd.read_table(pgs_file)
        traits = [t for t in PGS_TRAITS if t in praw.columns]
        pgs = praw[["sample_id"] + traits].rename(columns={"sample_id": "barcode"})
        df = df.merge(pgs, how="left", on="barcode").rename(columns={"EUR_ancestry": "EUR_ancestry_y"})
        out = OUT_DIR / f"df_v4_interaction_inputs_{anc}.csv"
        df.to_csv(out, index=False)
        print(f"[{anc}] wrote {out}  (n={len(df)})  clusters={df['Cluster'].value_counts(dropna=False).to_dict()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
