#!/usr/bin/env python3
"""
Most prominent networks / regions per cluster x direction (from the NBS edges).

Reads the per-cluster NBS-component edge lists (06_cluster_nbs.py) and, for each
cluster x direction (hypo = cluster<NT, hyper = cluster>NT), ranks:
  - cortical + subcortical NETWORKS by the number of component edges incident to
    their nodes (via the 4S156 `network_label`);
  - individual REGIONS the same way, split cortical vs subcortical.

Edge involvement = an edge contributes +1 to each of its two endpoints'
network / region. "prominence" = that summed incident-edge count (a node touched
by many component edges is a hub of the affected sub-network).

Outputs (outputs/tables/nbs/):
  nbs_prominent_networks.csv   cluster,direction,scope,network,n_edges,rank,pct
  nbs_prominent_regions.csv    cluster,direction,scope,region,network,n_edges,rank
"""
import os
from pathlib import Path
import pandas as pd

_SECTION = Path(__file__).resolve().parent.parent
NBS_DIR = _SECTION / "outputs" / "tables" / os.environ.get("NBS_SUBDIR", "nbs")
ATLAS_FILE = ("/Users/mfleury/POSTDOC/LIBRAIRY/eeg_mri-pipeline/ressources/atlases/"
              "SCHAEFER/atlas-4S156Parcels/atlas-4S156Parcels_dseg.tsv")

CLUSTERS = ["C1", "C2", "C3"]
DIRS = ["hypo", "hyper"]
TOP_REGIONS = 8   # per cluster x direction x scope


def main() -> int:
    atlas = pd.read_csv(ATLAS_FILE, sep="\t")
    net_of = dict(zip(atlas["label"], atlas["network_label"]))
    # Cortical/subcortical split from label_7network (robust to network_label
    # quirks, e.g. RH-MN is mislabelled "Limbic"): only true cortical Schaefer
    # parcels carry a "7Networks_*" label_7network; subcortical rows are "n/a".
    l7_of = dict(zip(atlas["label"], atlas["label_7network"].astype(str)))

    def scope_of_node(node):
        l7 = str(l7_of.get(node, ""))
        return "cortical" if l7.startswith("7Networks") else "subcortical"

    net_rows, reg_rows = [], []
    for C in CLUSTERS:
        fp = NBS_DIR / f"cluster_{C}_nbs_edges.csv"
        if not fp.exists():
            print(f"  missing {fp.name}; skipping {C}")
            continue
        d = pd.read_csv(fp)
        for direction in DIRS:
            sub = d[d.direction == direction]
            # explode each edge into its two endpoints
            eps = []
            for region in sub.region:
                clean = region[len("con_"):] if region.startswith("con_") else region
                parts = clean.split("/", 1)
                eps.append(parts[0])
                if len(parts) > 1:
                    eps.append(parts[1])
            if not eps:
                continue
            ep = pd.Series(eps, name="node")
            net = ep.map(net_of)
            valid = net.notna()
            ep, net = ep[valid], net[valid]
            scope = ep.map(scope_of_node)
            total = len(ep)

            # --- networks ranked by incident-edge count ---
            nc = (pd.DataFrame({"network": net, "scope": scope})
                  .value_counts().reset_index(name="n_edges"))
            nc = nc.sort_values("n_edges", ascending=False).reset_index(drop=True)
            for rank, r in enumerate(nc.itertuples(), 1):
                net_rows.append({"cluster": C, "direction": direction,
                                 "scope": r.scope, "network": r.network,
                                 "n_edges": int(r.n_edges), "rank": rank,
                                 "pct": round(100 * r.n_edges / total, 1)})

            # --- individual regions ranked, split cortical / subcortical ---
            rc = (pd.DataFrame({"region": ep.values, "network": net.values,
                                "scope": scope.values})
                  .value_counts().reset_index(name="n_edges"))
            for sc in ("cortical", "subcortical"):
                top = rc[rc.scope == sc].sort_values("n_edges", ascending=False).head(TOP_REGIONS)
                for rank, r in enumerate(top.itertuples(), 1):
                    reg_rows.append({"cluster": C, "direction": direction,
                                     "scope": sc, "region": r.region,
                                     "network": r.network, "n_edges": int(r.n_edges),
                                     "rank": rank})

    net_df = pd.DataFrame(net_rows)
    reg_df = pd.DataFrame(reg_rows)
    net_out = NBS_DIR / "nbs_prominent_networks.csv"
    reg_out = NBS_DIR / "nbs_prominent_regions.csv"
    net_df.to_csv(net_out, index=False)
    reg_df.to_csv(reg_out, index=False)
    print(f"Saved: {net_out}\nSaved: {reg_out}\n")

    # human-readable digest: top network (cortical + subcortical) per cell
    print("=== most prominent network per cluster x direction ===")
    for C in CLUSTERS:
        for direction in DIRS:
            cell = net_df[(net_df.cluster == C) & (net_df.direction == direction)]
            if cell.empty:
                print(f"  {C} {direction:5s}: (none)"); continue
            top_cort = cell[cell.scope == "cortical"].head(1)
            top_sub = cell[cell.scope == "subcortical"].head(1)
            ct = (f"{top_cort.network.iloc[0]} ({top_cort.n_edges.iloc[0]})"
                  if len(top_cort) else "-")
            sb = (f"{top_sub.network.iloc[0]} ({top_sub.n_edges.iloc[0]})"
                  if len(top_sub) else "-")
            print(f"  {C} {direction:5s}: cortical={ct:28s} subcortical={sb}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
