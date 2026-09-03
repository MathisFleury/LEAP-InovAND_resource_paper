"""
Compare the legacy fMRI preprocessing output (preprocessing_cohort.py) with
the refactored pipeline_ndd-style output in
`6_functional_analysis/preprocessing/outputs/df_conn_cohort_norm.csv`.

Reports:
  • sample overlap (legacy only / new only / both)
  • demographics of the common N (age, sex, cohort, group, mean_fd)
  • per-feature Pearson correlation between legacy and new
    (one value per `con_*` column, computed on common subjects)
  • per-subject Pearson correlation between legacy and new
    (one value per common subject, across 12 090 features)
  • Autism vs NT Welch t-test agreement: per-feature t-stats from both
    pipelines on the common-subject subset, plus FDR-significant overlap

Outputs land under:
    6_functional_analysis/preprocessing/outputs_validate_refactor/
(kept separate from preprocessing/outputs/, the real pipeline output this
script only reads from -- never write into that directory from here).
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind, pearsonr
from statsmodels.stats.multitest import multipletests

_SCRIPT_DIR = Path(__file__).resolve().parent
OUT_DIR = _SCRIPT_DIR / "outputs_validate_refactor"
OUT_DIR.mkdir(parents=True, exist_ok=True)

LEGACY_CSV = Path(
    "/Users/mfleury/POSTDOC/LIBRAIRY/eeg_mri-pipeline/results/dataframes/fmri/df_conn_cohort_norm.csv"
)
NEW_CSV = _SCRIPT_DIR / "outputs" / "df_conn_cohort_norm.csv"


def _load_with_dedup(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, low_memory=False)
    df["ID"] = df["ID"].astype(str)
    n0 = len(df)
    df = df.drop_duplicates(subset="ID", keep="first").reset_index(drop=True)
    if n0 != len(df):
        print(f"   dropped {n0 - len(df)} duplicate IDs in {path.name}")
    return df


def _con_columns(df: pd.DataFrame) -> list[str]:
    return [c for c in df.columns if c.startswith("con_")]


def sample_overlap(df_l: pd.DataFrame, df_n: pd.DataFrame) -> dict[str, int]:
    ids_l = set(df_l["ID"])
    ids_n = set(df_n["ID"])
    return {
        "legacy_N": len(ids_l),
        "new_N": len(ids_n),
        "common_N": len(ids_l & ids_n),
        "legacy_only_N": len(ids_l - ids_n),
        "new_only_N": len(ids_n - ids_l),
    }


def common_subjects_demographics(
    df_l: pd.DataFrame, df_n: pd.DataFrame, common_ids: list[str]
) -> pd.DataFrame:
    """Demographics drawn from the new dataframe (has clean group/cohort)."""
    sub_n = df_n.set_index("ID").loc[common_ids]
    rows = []
    rows.append({"feature": "N", "value": len(sub_n)})
    rows.append({"feature": "age_mean_yrs", "value": float(sub_n["age_yrs"].mean())})
    rows.append({"feature": "age_sd_yrs",   "value": float(sub_n["age_yrs"].std())})
    for col, label in [("Sex", "sex"), ("cohort", "cohort"), ("group", "group")]:
        if col in sub_n.columns:
            for k, v in sub_n[col].value_counts(dropna=False).to_dict().items():
                rows.append({"feature": f"{label}={k}", "value": int(v)})
    rows.append({"feature": "mean_fd_median_mm", "value": float(sub_n["mean_fd"].median())})
    return pd.DataFrame(rows)


def per_feature_correlation(
    df_l: pd.DataFrame, df_n: pd.DataFrame, common_ids: list[str], features: list[str]
) -> pd.DataFrame:
    L = df_l.set_index("ID").loc[common_ids, features].to_numpy(dtype=float)
    N = df_n.set_index("ID").loc[common_ids, features].to_numpy(dtype=float)

    # Column-wise Pearson via vectorised normalisation.
    Lc = L - L.mean(axis=0, keepdims=True)
    Nc = N - N.mean(axis=0, keepdims=True)
    num = (Lc * Nc).sum(axis=0)
    den = np.sqrt((Lc ** 2).sum(axis=0) * (Nc ** 2).sum(axis=0))
    r = np.divide(num, den, out=np.full_like(num, np.nan), where=den > 0)

    return pd.DataFrame({"feature": features, "r": r})


def per_subject_correlation(
    df_l: pd.DataFrame, df_n: pd.DataFrame, common_ids: list[str], features: list[str]
) -> pd.DataFrame:
    L = df_l.set_index("ID").loc[common_ids, features].to_numpy(dtype=float)
    N = df_n.set_index("ID").loc[common_ids, features].to_numpy(dtype=float)

    Lc = L - L.mean(axis=1, keepdims=True)
    Nc = N - N.mean(axis=1, keepdims=True)
    num = (Lc * Nc).sum(axis=1)
    den = np.sqrt((Lc ** 2).sum(axis=1) * (Nc ** 2).sum(axis=1))
    r = np.divide(num, den, out=np.full_like(num, np.nan), where=den > 0)

    return pd.DataFrame({"ID": common_ids, "r": r})


def autism_vs_nt_t_stats(
    df: pd.DataFrame, features: list[str], common_ids: list[str]
) -> pd.DataFrame:
    """Welch t-test (AutismS > 0, Control < 0) per feature on common-id sample.

    Uses the raw `PopulationS1` labels stored in both pipelines' CSVs:
    'AutismS' vs 'Control'. (The 'Autism' / 'NT' labels used downstream in
    01_autism_td_connectivity_analysis.py only appear after the kmeans
    cluster-table merge, which is irrelevant for a preprocessing audit.)
    """
    sub = df[df["ID"].isin(common_ids)].copy().set_index("ID")
    if "PopulationS1" not in sub.columns:
        return pd.DataFrame(columns=["feature", "t", "p"])
    sub = sub[sub["PopulationS1"].isin(["AutismS", "Control"])]
    aut = sub[sub["PopulationS1"] == "AutismS"][features].to_numpy(dtype=float)
    ctl = sub[sub["PopulationS1"] == "Control"][features].to_numpy(dtype=float)
    if aut.size == 0 or ctl.size == 0:
        return pd.DataFrame({"feature": features, "t": np.nan, "p": np.nan})
    t, p = ttest_ind(aut, ctl, equal_var=False, nan_policy="omit", axis=0)
    return pd.DataFrame({"feature": features, "t": np.asarray(t), "p": np.asarray(p)})


def make_summary_figure(
    feat_r: pd.DataFrame,
    subj_r: pd.DataFrame,
    t_legacy: pd.DataFrame,
    t_new: pd.DataFrame,
    out_path: Path,
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))

    ax = axes[0, 0]
    ax.hist(feat_r["r"].dropna(), bins=50, color="#5B9BD5", edgecolor="black", alpha=0.85)
    ax.axvline(feat_r["r"].median(), color="red", linestyle="--",
               label=f"median = {feat_r['r'].median():.3f}")
    ax.set_xlabel("Pearson r per feature  (legacy vs new, across common subjects)")
    ax.set_ylabel("# features")
    ax.set_title(f"Per-feature agreement  (n_features = {len(feat_r):,})")
    ax.legend()
    ax.grid(True, alpha=0.3)

    ax = axes[0, 1]
    ax.hist(subj_r["r"].dropna(), bins=50, color="#70AD47", edgecolor="black", alpha=0.85)
    ax.axvline(subj_r["r"].median(), color="red", linestyle="--",
               label=f"median = {subj_r['r'].median():.3f}")
    ax.set_xlabel("Pearson r per subject  (legacy vs new, across common features)")
    ax.set_ylabel("# subjects")
    ax.set_title(f"Per-subject agreement  (N = {len(subj_r):,})")
    ax.legend()
    ax.grid(True, alpha=0.3)

    ax = axes[1, 0]
    df_t = t_legacy.merge(t_new, on="feature", suffixes=("_legacy", "_new")).dropna()
    if not df_t.empty:
        r_t, _ = pearsonr(df_t["t_legacy"], df_t["t_new"])
        ax.scatter(df_t["t_legacy"], df_t["t_new"], s=4, alpha=0.35, color="#7030A0")
        lo = float(min(df_t["t_legacy"].min(), df_t["t_new"].min()))
        hi = float(max(df_t["t_legacy"].max(), df_t["t_new"].max()))
        ax.plot([lo, hi], [lo, hi], "r--", linewidth=1)
        ax.set_xlabel("legacy Welch t  (Autism vs NT)")
        ax.set_ylabel("new Welch t  (Autism vs NT)")
        ax.set_title(f"t-stat agreement   r = {r_t:.3f}")
        ax.grid(True, alpha=0.3)

    ax = axes[1, 1]
    if not df_t.empty:
        df_t["fdr_legacy"] = multipletests(df_t["p_legacy"], method="fdr_bh")[1]
        df_t["fdr_new"]    = multipletests(df_t["p_new"],    method="fdr_bh")[1]
        sig_legacy = set(df_t.loc[df_t["fdr_legacy"] < 0.05, "feature"])
        sig_new    = set(df_t.loc[df_t["fdr_new"]    < 0.05, "feature"])
        both       = sig_legacy & sig_new
        ax.bar(["legacy", "new", "both"],
               [len(sig_legacy), len(sig_new), len(both)],
               color=["#5B9BD5", "#ED7D31", "#70AD47"])
        for i, n in enumerate([len(sig_legacy), len(sig_new), len(both)]):
            ax.text(i, n, str(n), ha="center", va="bottom")
        ax.set_ylabel("# FDR-significant connections  (q < 0.05)")
        ax.set_title("Autism-vs-NT hits — overlap")
        ax.grid(True, axis="y", alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    print("=" * 70)
    print("LEGACY vs NEW fMRI PREPROCESSING — comparison")
    print("=" * 70)
    print(f"legacy: {LEGACY_CSV}")
    print(f"new   : {NEW_CSV}")

    df_l = _load_with_dedup(LEGACY_CSV)
    df_n = _load_with_dedup(NEW_CSV)

    overlap = sample_overlap(df_l, df_n)
    print("\n[overlap]")
    for k, v in overlap.items():
        print(f"  {k:<20s} {v:,}")

    common_ids = sorted(set(df_l["ID"]) & set(df_n["ID"]))
    features = sorted(set(_con_columns(df_l)) & set(_con_columns(df_n)))

    pd.DataFrame(list(overlap.items()), columns=["stat", "N"]).to_csv(
        OUT_DIR / "sample_overlap.csv", index=False
    )

    demo = common_subjects_demographics(df_l, df_n, common_ids)
    demo.to_csv(OUT_DIR / "common_subjects_demographics.csv", index=False)

    print(f"\n[per-feature correlation across {len(common_ids):,} common subjects]")
    feat_r = per_feature_correlation(df_l, df_n, common_ids, features)
    feat_r.to_csv(OUT_DIR / "per_feature_correlation.csv", index=False)
    print(f"  median r = {feat_r['r'].median():.3f}    mean r = {feat_r['r'].mean():.3f}")
    print(f"  |r|<0.2: {(feat_r['r'].abs() < 0.2).sum():,}    "
          f"|r|>0.8: {(feat_r['r'].abs() > 0.8).sum():,}")

    print(f"\n[per-subject correlation across {len(features):,} common features]")
    subj_r = per_subject_correlation(df_l, df_n, common_ids, features)
    subj_r.to_csv(OUT_DIR / "per_subject_correlation.csv", index=False)
    print(f"  median r = {subj_r['r'].median():.3f}    mean r = {subj_r['r'].mean():.3f}")
    print(f"  |r|<0.2: {(subj_r['r'].abs() < 0.2).sum():,}    "
          f"|r|>0.8: {(subj_r['r'].abs() > 0.8).sum():,}")

    print(f"\n[Welch t-test Autism vs NT — common-id subset]")
    t_l = autism_vs_nt_t_stats(df_l, features, common_ids)
    t_n = autism_vs_nt_t_stats(df_n, features, common_ids)
    df_t = t_l.merge(t_n, on="feature", suffixes=("_legacy", "_new")).dropna()
    df_t["fdr_legacy"] = multipletests(df_t["p_legacy"], method="fdr_bh")[1]
    df_t["fdr_new"]    = multipletests(df_t["p_new"],    method="fdr_bh")[1]
    df_t.to_csv(OUT_DIR / "autism_vs_nt_t_comparison.csv", index=False)
    r_t, _ = pearsonr(df_t["t_legacy"], df_t["t_new"])
    sig_legacy = set(df_t.loc[df_t["fdr_legacy"] < 0.05, "feature"])
    sig_new    = set(df_t.loc[df_t["fdr_new"]    < 0.05, "feature"])
    print(f"  r(t_legacy, t_new) = {r_t:.3f}")
    print(f"  FDR-sig legacy: {len(sig_legacy):4d}   "
          f"FDR-sig new: {len(sig_new):4d}   "
          f"overlap: {len(sig_legacy & sig_new):4d}")

    fig_path = OUT_DIR / "legacy_vs_new_summary.pdf"
    make_summary_figure(feat_r, subj_r, t_l, t_n, fig_path)
    print(f"\n  wrote  {fig_path.name}")

    print("\n✓ comparison complete  →  6_functional_analysis/preprocessing/outputs_validate_refactor/")


if __name__ == "__main__":
    main()
