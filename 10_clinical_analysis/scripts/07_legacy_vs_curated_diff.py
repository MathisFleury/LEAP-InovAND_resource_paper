#!/usr/bin/env python3
# =============================================================================
# 06 - Legacy vs curated clustering-input diff
# =============================================================================
# Both use the SAME method (k-means on IQ x SRS); the only difference is the
# clinical source — frozen/paper (df_clusters_complete_kmeans.csv) vs the
# curated re-derivation (individuals_metrics_with_clusters_curated.csv). This
# isolates the DATA changes: cohort membership, IQ, SRS-2, cluster label,
# population label — per subject.
#
# Outputs (10_clinical_analysis/outputs/tables/):
#   legacy_vs_curated_diff.csv      one row per ID that changed (what & how)
#   legacy_vs_curated_summary.csv   counts per change type
# =============================================================================

import os
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _config import CLUSTERS_FILE, CLUSTERS_CURATED_FILE, TABLES_DIR  # noqa: E402

FIELDS = ["IQ", "SRS_tscore", "Cluster", "PopulationS1", "Population1"]
TOL = 1e-6  # numeric equality tolerance for IQ / SRS


def _load(path):
    d = pd.read_csv(path, low_memory=False).drop_duplicates("ID")
    d["ID"] = d["ID"].astype(str)
    return d.set_index("ID")


def _differs(a, b):
    if pd.isna(a) and pd.isna(b):
        return False
    if pd.isna(a) != pd.isna(b):
        return True
    if isinstance(a, (int, float, np.floating)) and isinstance(b, (int, float, np.floating)):
        return abs(float(a) - float(b)) > TOL
    return str(a) != str(b)


def main():
    leg, cur = _load(CLUSTERS_FILE), _load(CLUSTERS_CURATED_FILE)
    leg_ids, cur_ids = set(leg.index), set(cur.index)
    both = sorted(leg_ids & cur_ids)

    # --- cohort membership ---
    only_leg = sorted(leg_ids - cur_ids)
    only_cur = sorted(cur_ids - leg_ids)

    # --- per-field diffs on shared IDs ---
    rows = []
    for i in both:
        rl, rc = leg.loc[i], cur.loc[i]
        changed = {f: (rl.get(f), rc.get(f)) for f in FIELDS if _differs(rl.get(f), rc.get(f))}
        if changed:
            rec = {"ID": i, "changed_fields": ";".join(changed)}
            for f, (lv, cv) in changed.items():
                rec[f"{f}_legacy"], rec[f"{f}_curated"] = lv, cv
            rows.append(rec)
    diff = pd.DataFrame(rows)

    os.makedirs(TABLES_DIR, exist_ok=True)
    diff_path = os.path.join(TABLES_DIR, "legacy_vs_curated_diff.csv")
    diff.to_csv(diff_path, index=False)

    # --- summary ---
    def n_changed(f):
        return int(diff["changed_fields"].str.contains(rf"\b{f}\b").sum()) if not diff.empty else 0

    summ = pd.DataFrame([
        {"metric": "IDs in legacy", "n": len(leg_ids)},
        {"metric": "IDs in curated", "n": len(cur_ids)},
        {"metric": "IDs in both", "n": len(both)},
        {"metric": "only in legacy", "n": len(only_leg)},
        {"metric": "only in curated", "n": len(only_cur)},
        {"metric": "shared IDs with any change", "n": 0 if diff.empty else len(diff)},
        {"metric": "IQ changed", "n": n_changed("IQ")},
        {"metric": "SRS_tscore changed", "n": n_changed("SRS_tscore")},
        {"metric": "Cluster changed", "n": n_changed("Cluster")},
        {"metric": "PopulationS1 changed", "n": n_changed("PopulationS1")},
        {"metric": "Population1 changed", "n": n_changed("Population1")},
    ])
    summ_path = os.path.join(TABLES_DIR, "legacy_vs_curated_summary.csv")
    summ.to_csv(summ_path, index=False)

    # cluster transition breakdown (among clustered-in-both)
    if not diff.empty and "Cluster_legacy" in diff.columns:
        cl = diff.dropna(subset=["Cluster_legacy", "Cluster_curated"])
        trans = (cl.groupby(["Cluster_legacy", "Cluster_curated"]).size()
                 .reset_index(name="n").sort_values("n", ascending=False))
        trans.to_csv(os.path.join(TABLES_DIR, "legacy_vs_curated_cluster_transitions.csv"), index=False)

    print(summ.to_string(index=False))
    print(f"\nSaved: {diff_path}\n       {summ_path}")
    if not diff.empty:
        med = {f: (diff[f"{f}_legacy"] - diff[f"{f}_curated"]).abs().median()
               for f in ("IQ", "SRS_tscore")
               if f"{f}_legacy" in diff.columns and pd.api.types.is_numeric_dtype(diff[f"{f}_legacy"])}
        print("median |legacy-curated| among changed:", {k: round(v, 2) for k, v in med.items()})


if __name__ == "__main__":
    main()
