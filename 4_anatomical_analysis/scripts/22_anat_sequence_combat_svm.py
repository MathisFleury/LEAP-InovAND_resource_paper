#!/usr/bin/env python3
"""
Acquisition-unit (scanner) harmonisation check for structural MRI -- SVM AUC
before/after ComBat, the anatomical companion to
6_functional_analysis/concat/scripts/13_sequence_combat_svm.py (same method,
applied to the FreeSurfer feature matrix instead of connectivity edges).

Method (matches 13_sequence_combat_svm.py exactly): one-vs-one linear-SVM
classification of scanner (the `batch` column, scanner-level per
4_anatomical_analysis/preprocessing) from the regional features, 5-fold CV
ROC AUC, for every scanner pair with >=MINN subjects, BEFORE ComBat (raw,
post feature-selection) and AFTER (harmonised). AUC ~0.5 after ComBat =>
acquisition-unit effect removed.

Both states are computed in-memory from the same preprocessing pipeline run
(step_01..step_03 = pre-ComBat features; + step_04 = post-ComBat), so pre and
post are guaranteed to be the exact same subjects/features/batch labels --
no join needed (unlike the functional script, which joins two separately-
saved files on `path`).

Output -> outputs/qc_sequence/anat_sequence_svm_auc_bymachine.{csv,pdf,png}
"""
import os
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.svm import SVC
from sklearn.decomposition import PCA
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score, StratifiedKFold

_SECTION = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_SECTION / "preprocessing"))

from step_01_load import load_dataframe
from step_02_filter_qc import (
    restrict_to_ndd, filter_qc, add_euler,
    drop_missing_covariates, drop_tiny_cohorts, keep_best_run,
)
from step_02b_assign_batch import assign_batch
from step_03_select_features import (
    select_features, restrict_to_numeric, drop_subjects_missing_features,
    drop_zero_variance_per_batch,
)
from step_04_combat import run_combat
from config import BATCH_COL

OUT = _SECTION / "outputs" / "qc_sequence"
OUT.mkdir(parents=True, exist_ok=True)
MINN = int(os.environ.get("MINN", "20"))     # min subjects/scanner to include
NPCA = 100
RNG = 42


def auc_pair(X, y):
    """5-fold CV ROC AUC, one-vs-one linear SVM on PCA-reduced features.
    n_components is capped below the smallest training-fold size (4/5 of n)."""
    ncomp = max(2, min(NPCA, int(0.7 * len(y)), X.shape[1]))
    clf = make_pipeline(StandardScaler(),
                        PCA(n_components=ncomp, random_state=RNG),
                        SVC(kernel="linear"))
    cv = StratifiedKFold(5, shuffle=True, random_state=RNG)
    return cross_val_score(clf, X, y, cv=cv, scoring="roc_auc").mean()


def main() -> int:
    # ---------- rebuild the pipeline up to (pre) and through (post) ComBat --
    df = load_dataframe()
    df = restrict_to_ndd(df)
    df = filter_qc(df)
    df = add_euler(df)
    df = assign_batch(df)
    df = keep_best_run(df)
    df = drop_missing_covariates(df)
    df = drop_tiny_cohorts(df)

    features = select_features(df)
    features = restrict_to_numeric(df, features)
    df = drop_subjects_missing_features(df, features)
    features = drop_zero_variance_per_batch(df, features)

    pre = df[features].to_numpy(float)                       # before ComBat
    df_post, _ = run_combat(df, features)
    post = df_post[features].to_numpy(float)                 # after ComBat
    batch = df[BATCH_COL].values

    scanners = [s for s, n in pd.Series(batch).value_counts().items() if n >= MINN]
    print(f"scanners with >= {MINN} subjects: {len(scanners)} -> {scanners}")
    print(f"anatomical features: {len(features)}")

    rows = []
    for a, b in combinations(scanners, 2):
        m = np.isin(batch, [a, b])
        y = (batch[m] == a).astype(int)
        if y.sum() < 5 or (len(y) - y.sum()) < 5:
            continue
        a_pre = auc_pair(np.nan_to_num(pre[m]), y)
        a_post = auc_pair(np.nan_to_num(post[m]), y)
        rows.append({"scanner_a": a, "scanner_b": b, "n_a": int(y.sum()),
                     "n_b": int(len(y) - y.sum()),
                     "auc_before": round(a_pre, 3), "auc_after": round(a_post, 3)})
        print(f"  {a:12s} vs {b:12s}: AUC before={a_pre:.3f}  after={a_post:.3f}")

    res = pd.DataFrame(rows)
    res.to_csv(OUT / "anat_sequence_svm_auc_bymachine.csv", index=False)
    print(f"\nMean one-vs-one AUC — before ComBat: {res.auc_before.mean():.3f} | "
          f"after ComBat: {res.auc_after.mean():.3f}  (0.5 = chance)")

    # Per-scanner mean AUC vs every other scanner -- matches the plot style of
    # 6_functional_analysis/concat/outputs/qc_sequence/sequence_svm_auc_bymachine.pdf
    # (a grouped bar chart per scanner, not a paired before/after scatter).
    all_scanners = sorted(set(res.scanner_a) | set(res.scanner_b))
    before_mean, after_mean = {}, {}
    for s in all_scanners:
        m = (res.scanner_a == s) | (res.scanner_b == s)
        before_mean[s] = res.loc[m, "auc_before"].mean()
        after_mean[s] = res.loc[m, "auc_after"].mean()
    order = sorted(all_scanners, key=lambda s: -before_mean[s])

    x = np.arange(len(order)); w = 0.38
    fig, ax = plt.subplots(figsize=(11, 5.5))
    ax.bar(x - w / 2, [before_mean[s] for s in order], width=w, color="#4A90E2", label="before ComBat")
    ax.bar(x + w / 2, [after_mean[s] for s in order], width=w, color="#7a0010", label="after ComBat")
    ax.axhline(0.5, ls="--", c="0.4", lw=1)
    ax.set_xticks(x); ax.set_xticklabels(order, rotation=45, ha="right")
    ax.set_ylabel("mean one-vs-one SVM ROC AUC\n(vs other scanners)")
    ax.set_title("sMRI — scanner classifiability before/after ComBat")
    ax.set_ylim(0, 1.05)
    ax.legend()
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT / "anat_sequence_svm_auc_bymachine.pdf", bbox_inches="tight")
    fig.savefig(OUT / "anat_sequence_svm_auc_bymachine.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {OUT}/anat_sequence_svm_auc_bymachine.csv + .pdf/.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
