"""Age distribution across the curated k-means clusters (Reviewer #1 age follow-up).

Closes the gap that the interaction analysis (01-06) doesn't cover: does cluster
membership relate to age? Answer: yes, modestly (C1 oldest, C2 youngest) — reported
so the manuscript can state age is a cluster covariate rather than claim clusters are
age-neutral. Kruskal-Wallis + pairwise Mann-Whitney (Bonferroni).

# ponytail: one-shot descriptive script, no config layer — reads the label file directly.
"""
import pandas as pd, numpy as np
from scipy import stats
from itertools import combinations
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
LABELS = ROOT / "1_clustering/outputs/curated/tables/individuals_metrics_with_clusters_curated.csv"
OUT = Path(__file__).resolve().parents[1] / "outputs/tables"
OUT.mkdir(parents=True, exist_ok=True)

AUT = {"Autism without IDD", "Autism with IDD", "Autism to exclude"}


def summarise(df):
    rows = []
    for c in ["C1", "C2", "C3"]:
        a = df.loc[df["Cluster"] == c, "age"].dropna()
        rows.append(dict(cluster=c, n=len(a), mean=a.mean(), sd=a.std(),
                         median=a.median(), q1=a.quantile(.25), q3=a.quantile(.75),
                         min=a.min(), max=a.max()))
    return pd.DataFrame(rows)


def kw_and_pairwise(df, label):
    g = [df.loc[df["Cluster"] == c, "age"].dropna() for c in ["C1", "C2", "C3"]]
    g = [x for x in g if len(x)]
    H, p = stats.kruskal(*g)
    print(f"\n[{label}] Kruskal-Wallis: H={H:.3f}, p={p:.3g}")
    labs = ["C1", "C2", "C3"]
    out = [dict(comparison=label, test="kruskal", stat=H, p=p, p_bonf=p)]
    if label == "all_labelled":  # pairwise only for the full set (3 non-empty groups)
        for i, j in combinations(range(3), 2):
            U, pp = stats.mannwhitneyu(g[i], g[j], alternative="two-sided")
            print(f"  {labs[i]} vs {labs[j]}: p={pp:.3g}, p_bonf={min(pp*3,1):.3g}")
            out.append(dict(comparison=f"{labs[i]}_vs_{labs[j]}", test="mannwhitneyu",
                            stat=U, p=pp, p_bonf=min(pp * 3, 1)))
    return out


def main():
    df = pd.read_csv(LABELS)
    df = df[df["Cluster"].notna()].copy()
    df["age"] = pd.to_numeric(df["age_yrs"], errors="coerce")
    # The curated cohort keeps growing (was 1037 when this script was written,
    # per root CLAUDE.md); a hard equality check goes stale on every update, so
    # guard only against loading an empty/garbage file instead.
    assert len(df) > 900, f"suspiciously few labelled rows: {len(df)}"
    print(f"[info] {len(df)} labelled rows (was 1037 when this script was written)")

    tests = []
    print("== All labelled ==")
    s_all = summarise(df); print(s_all.round(2).to_string(index=False))
    s_all.to_csv(OUT / "age_by_cluster_all.csv", index=False)
    tests += kw_and_pairwise(df, "all_labelled")

    aut = df[df["population"].isin(AUT)]
    print("\n== Autistic only ==")
    s_aut = summarise(aut); print(s_aut.round(2).to_string(index=False))
    s_aut.to_csv(OUT / "age_by_cluster_autistic.csv", index=False)
    tests += kw_and_pairwise(aut, "autistic_only")

    pd.DataFrame(tests).to_csv(OUT / "age_by_cluster_tests.csv", index=False)
    print(f"\nwrote -> {OUT}/age_by_cluster_{{all,autistic,tests}}.csv")


if __name__ == "__main__":
    main()
