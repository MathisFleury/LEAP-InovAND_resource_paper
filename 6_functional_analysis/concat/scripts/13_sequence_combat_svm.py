#!/usr/bin/env python3
"""
Acquisition-unit (scanner/sequence) harmonisation check — SVM AUC before/after
ComBat, the Supplementary-Fig-19b analogue extended from *site* to *scanner*.

Reviewer point: sequences differ between AND within sites; their influence
should be assessed and (if present) mitigated with ComBat. The ComBat batch here
is `mri_machine`, which is the finest acquisition-unit label we have — it
separates LEAP sites (each a distinct vendor/protocol) AND the multiple INOVAND
scanners/field strengths at a single site (Ingenia / Intera / Ingenia 3T /
1.5T variants). This script tests whether that acquisition-unit signal is
removed by ComBat.

Method (matches Supp 19b): one-vs-one linear-SVM classification of scanner from
the connectivity edges, 5-fold CV ROC AUC, for every scanner pair with >=MINN
subjects, BEFORE ComBat (raw Fisher-z) and AFTER (harmonised, as used in the
analysis). AUC ~0.5 after ComBat => acquisition-unit effect removed.

Both inputs are local; matched on the `path` column. Featured variant via
XCPD_VARIANT_TAG (default nogsr).
Outputs -> outputs/qc_sequence/: sequence_svm_auc.csv, sequence_svm_auc.{pdf,png}
"""
import os
from pathlib import Path
from itertools import combinations
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.svm import SVC
from sklearn.decomposition import PCA
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score, StratifiedKFold

_SECTION = Path(__file__).resolve().parent.parent
V = os.environ.get("XCPD_VARIANT_TAG", "nogsr")
PRE = _SECTION / "preprocessing" / "outputs" / V / f"df_conn_raw_{V}.csv"
POST = _SECTION / "preprocessing" / "outputs" / V / f"df_conn_cohort_norm_{V}.csv"
OUT = _SECTION / "outputs" / "qc_sequence"
OUT.mkdir(parents=True, exist_ok=True)
MINN = int(os.environ.get("MINN", "20"))     # min subjects/scanner to include
NPCA = 100
RNG = 42


def auc_pair(X, y):
    """5-fold CV ROC AUC, one-vs-one linear SVM on PCA-reduced edges.
    n_components is capped below the smallest training-fold size (4/5 of n)."""
    ncomp = max(2, min(NPCA, int(0.7 * len(y)), X.shape[1]))
    clf = make_pipeline(StandardScaler(),
                        PCA(n_components=ncomp, random_state=RNG),
                        SVC(kernel="linear"))
    cv = StratifiedKFold(5, shuffle=True, random_state=RNG)
    return cross_val_score(clf, X, y, cv=cv, scoring="roc_auc").mean()


def main() -> int:
    post = pd.read_csv(POST, low_memory=False)
    pre = pd.read_csv(PRE, low_memory=False)
    con = [c for c in post.columns if c.startswith("con_")]
    con_pre = [c for c in pre.columns if c.startswith("con_")]
    con = [c for c in con if c in con_pre]

    # attach scanner label to the raw (pre) rows via path; align to post subjects
    pre_m = pre.merge(post[["path", "mri_machine"]], on="path", how="inner")
    scanners = [s for s, n in post["mri_machine"].value_counts().items() if n >= MINN]
    print(f"[{V}] scanners with >= {MINN} subjects: {len(scanners)} -> {scanners}")
    print(f"connectivity edges: {len(con)}")

    rows = []
    for a, b in combinations(scanners, 2):
        # POST (harmonised)
        m = post["mri_machine"].isin([a, b])
        Xpost = post.loc[m, con].to_numpy(float); ypost = (post.loc[m, "mri_machine"] == a).astype(int).values
        # PRE (raw), same subjects via path
        mp = pre_m["mri_machine"].isin([a, b])
        Xpre = pre_m.loc[mp, con].to_numpy(float); ypre = (pre_m.loc[mp, "mri_machine"] == a).astype(int).values
        if ypost.sum() < 5 or (len(ypost) - ypost.sum()) < 5:
            continue
        a_pre = auc_pair(np.nan_to_num(Xpre), ypre)
        a_post = auc_pair(np.nan_to_num(Xpost), ypost)
        rows.append({"scanner_a": a, "scanner_b": b, "n_a": int(ypost.sum()),
                     "n_b": int(len(ypost) - ypost.sum()),
                     "auc_before": round(a_pre, 3), "auc_after": round(a_post, 3)})
        print(f"  {a:12s} vs {b:12s}: AUC before={a_pre:.3f}  after={a_post:.3f}")
    res = pd.DataFrame(rows)
    res.to_csv(OUT / "sequence_svm_auc.csv", index=False)
    print(f"\nMean one-vs-one AUC — before ComBat: {res.auc_before.mean():.3f} | "
          f"after ComBat: {res.auc_after.mean():.3f}  (0.5 = chance)")

    # plot: paired before/after
    fig, ax = plt.subplots(figsize=(6, 5))
    for _, r in res.iterrows():
        ax.plot([0, 1], [r.auc_before, r.auc_after], color="0.7", lw=0.8, zorder=1)
    ax.scatter(np.zeros(len(res)), res.auc_before, s=30, c="#4A90E2", zorder=2, label="before")
    ax.scatter(np.ones(len(res)), res.auc_after, s=30, c="#7a0010", zorder=2, label="after")
    ax.axhline(0.5, ls="--", c="0.4", lw=1)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["before ComBat", "after ComBat"])
    ax.set_ylabel("one-vs-one SVM ROC AUC (scanner classification)")
    ax.set_ylim(0.4, 1.02); ax.set_xlim(-0.3, 1.3)
    ax.set_title(f"Scanner/acquisition-unit classifiability ({V})\n"
                 f"mean AUC {res.auc_before.mean():.2f} -> {res.auc_after.mean():.2f}")
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    fig.tight_layout(); fig.savefig(OUT / "sequence_svm_auc.pdf", bbox_inches="tight")
    fig.savefig(OUT / "sequence_svm_auc.png", dpi=200, bbox_inches="tight"); plt.close(fig)
    print(f"Saved: {OUT}/sequence_svm_auc.csv + .pdf/.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
