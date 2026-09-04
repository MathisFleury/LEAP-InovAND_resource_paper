"""Network-pair N_raw / N_norm summary of the significant NBS component edges,
reproducing the manuscript's original within/between-network reporting for the
no-GSR concat (primary) pipeline.

N_raw  = number of significant NBS edges between two networks.
N_norm = N_raw / number of possible connections between those networks
         (within: n*(n-1)/2; between: n_a*n_b), using the 4S atlas node counts.

# ponytail: reads the same atlas the analysis uses; no reimplementation of NBS.
"""
import os
import pandas as pd
from collections import Counter

SEC = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ATLAS = '/Users/mfleury/POSTDOC/LIBRAIRY/eeg_mri-pipeline/ressources/atlases/SCHAEFER/atlas-4S156Parcels/atlas-4S156Parcels_dseg.tsv'
V = os.environ.get("XCPD_VARIANT_TAG", "nogsr_concat")
NBS = f"{SEC}/outputs/tables/nbs_wholegroup_{V}/nbs_edges.csv"
OUT = f"{SEC}/outputs/{V}/tables/autism_vs_td_network_proportions.csv"


def main():
    atlas = pd.read_csv(ATLAS, sep='\t')
    lab2net = dict(zip(atlas['label'], atlas['network_label']))
    # node counts must reflect the ROIs actually RETAINED in the connectivity
    # (after coverage QC), not the full atlas — else N_norm denominators are wrong.
    feats = open(f"{SEC}/preprocessing/outputs/nogsr/features_used_nogsr.txt").read().split()
    rois = set()
    for f in feats:
        a, b = f.replace('con_', '').split('/')
        rois.update((a, b))
    npn = Counter(lab2net.get(r) for r in rois)
    print(f"ROIs retained in connectivity: {len(rois)} of {len(atlas)} atlas parcels")
    low = {k: v for k, v in npn.items() if v <= 2}
    if low:
        print(f"LOW-COVERAGE networks (<=2 parcels retained) — N_norm unreliable: {low}")
    nbs = pd.read_csv(NBS)

    rows = []
    for direction in ['hypo', 'hyper']:
        sub = nbs[nbs['direction'] == direction]
        raw = Counter()
        for reg in sub['region']:
            a, b = reg.replace('con_', '').split('/')
            na, nb = lab2net.get(a), lab2net.get(b)
            if na is None or nb is None:
                continue
            raw[tuple(sorted((na, nb)))] += 1
        for (na, nb), n in raw.items():
            poss = npn[na] * (npn[na] - 1) / 2 if na == nb else npn[na] * npn[nb]
            rows.append(dict(direction=direction, net_a=na, net_b=nb, within=na == nb,
                             N_raw=n, possible=int(poss), N_norm=round(n / poss, 4)))
    df = pd.DataFrame(rows).sort_values(['direction', 'N_norm'], ascending=[True, False])
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    df.to_csv(OUT, index=False)
    for d in ['hypo', 'hyper']:
        print(f"\n=== {d} — top network pairs by N_norm ===")
        print(df[df.direction == d].head(10).to_string(index=False))
    print(f"\nwrote -> {OUT}")
    print(f"atlas node counts: {dict(npn)}  (total {len(atlas)})")


if __name__ == "__main__":
    main()
