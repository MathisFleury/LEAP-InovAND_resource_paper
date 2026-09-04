#!/usr/bin/env python3
"""
Step 1 of the XCP-D v0.11 post-processing: build connectivity from timeseries.

XCP-D 0.11 emits parcel mean timeseries (no connectivity matrix), so for each
subject-session-RUN we read the per-run 4S156 `stat-mean_timeseries.tsv`
(timepoints x 156 parcels), compute the Pearson correlation matrix, Fisher
z-transform it, and flatten the upper triangle into `con_<src>/<tgt>` features
(156*155/2 = 12090) — the same feature convention the downstream ComBat /
regression / analysis scripts expect.

Per-run (NOT concatenated): we deliberately skip XCP-D's --combine-runs output
and keep individual runs, each tagged with its own v0.11 mean_fd, so the
downstream lowest-FD dedup picks one best run per subject — mirroring
6_functional_analysis/non_concat. Proper cross-run concatenation is future work (TODO #1).

The canonical 156-parcel order is taken from the first timeseries found (the
atlas dseg lists 158 rows across sub-atlases, so we use the actual XCP-D output
columns and reindex every subject to that order).

Output: outputs/<variant>/df_conn_raw_<variant>.csv
Run for each pipeline:  XCPD_VARIANT=nogsr|gsr python step_01_build_connectivity.py
"""
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import COHORT_IDP, TS_GLOB, CONN_RAW_CSV, VARIANT, ATLAS_FILE, OUT_BASE  # noqa: E402

# Connectivity estimator:
#   emp : empirical sample correlation (np.corrcoef) + per-subject IQR rescale
#         -> arctanh(±3)  [matches the original preprocessing_cohort.py]
#   lw  : Ledoit-Wolf shrinkage covariance -> correlation -> Fisher-z arctanh
#         [matches "Pearson correlation with Ledoit-Wolf shrinkage estimator"]
CONN_METHOD = os.environ.get("CONN_METHOD", "emp").lower()


def _lw_corr(X: np.ndarray) -> np.ndarray:
    """Ledoit-Wolf shrinkage correlation matrix; NaN for bad (NaN/flat) parcels."""
    from sklearn.covariance import LedoitWolf
    n = X.shape[1]
    C = np.full((n, n), np.nan)
    good = np.isfinite(X).all(axis=0) & (np.nanstd(X, axis=0) > 0)
    if good.sum() >= 2:
        cov = LedoitWolf().fit(X[:, good]).covariance_
        d = np.sqrt(np.diag(cov))
        corr = cov / np.outer(d, d)
        idx = np.where(good)[0]
        C[np.ix_(idx, idx)] = corr
    return C


def _curation_spec():
    """From the 4S156 dseg, the parcels to collapse/drop before correlation.

    The emitted timeseries splits the thalamus into 14 nuclei (the coarse
    LH/RH_Thalamus in the dseg have no voxels), so we average the 7 LEFT nuclei
    into LH_Thalamus and the 7 RIGHT nuclei into RH_Thalamus → 2 columns (L and R
    kept separate). Cerebellum is reduced to Region1.
    """
    d = pd.read_csv(ATLAS_FILE, sep="\t")
    thal = d.loc[d["network_label"].eq("Thalamus"), "label"].tolist()
    thal = [t for t in thal if t not in ("LH_Thalamus", "RH_Thalamus")]  # 14 nuclei
    thal_lh = [t for t in thal if t.startswith("LH-")]
    thal_rh = [t for t in thal if t.startswith("RH-")]
    cereb_drop = [f"Cerebellar_Region{i}" for i in range(2, 11)]          # keep Region1
    return thal_lh, thal_rh, cereb_drop


THAL_LH, THAL_RH, CEREB_DROP = _curation_spec()


def curate(ts: pd.DataFrame) -> pd.DataFrame:
    """Average the 7 LH nuclei -> LH_Thalamus and 7 RH nuclei -> RH_Thalamus
    (2 columns), drop Cerebellar 2-10."""
    ts = ts.copy()
    ts["LH_Thalamus"] = ts[THAL_LH].mean(axis=1)
    ts["RH_Thalamus"] = ts[THAL_RH].mean(axis=1)
    return ts.drop(columns=THAL_LH + THAL_RH + CEREB_DROP)


def _read_mean_fd(ts_path: Path) -> float:
    """v0.11 mean_fd for this run, from the sibling linc_qc tsv."""
    qc = ts_path.with_name(
        ts_path.name.replace("_seg-4S156Parcels_stat-mean_timeseries", "_desc-linc_qc"))
    try:
        return float(pd.read_csv(qc, sep="\t")["mean_fd"].iloc[0])
    except Exception:
        return float("nan")


def main() -> int:
    # collect PER-RUN timeseries (NOT the --combine-runs concatenated file: a
    # subject with >1 run also emits a no-run combined file, which we skip;
    # single-run subjects emit only the no-run file, which we keep). Proper
    # cross-run concatenation is deferred — see TODO #1.
    found = []
    for cohort, idp in COHORT_IDP.items():
        if not idp.exists():
            print(f"  [{cohort}] IDP dir missing (variant not run yet?): {idp}")
            continue
        all_ts = sorted(idp.glob(f"sub-*/ses-*/func/{TS_GLOB}"))
        by_dir: dict[Path, list[Path]] = {}
        for p in all_ts:
            by_dir.setdefault(p.parent, []).append(p)
        files = []
        for funcdir, ps in by_dir.items():
            per_run = [p for p in ps if "_run-" in p.name]
            files += per_run if per_run else ps          # skip combined when runs exist
        print(f"  [{cohort}] {len(files)} per-run timeseries "
              f"({len(all_ts) - len(files)} combined files skipped)")
        found += [(cohort, p) for p in files]
    if not found:
        print("No timeseries found — nothing written."); return 1

    # canonical parcel order = columns of the first timeseries, AFTER curation
    # (thalamic nuclei collapsed to LH/RH_Thalamus, cerebellum -> Region1).
    labels = list(curate(pd.read_csv(found[0][1], sep="\t")).columns)
    N = len(labels)
    iu = np.triu_indices(N, k=1)
    # Name edges with the LATER-ordered node first, to match the legacy
    # (eeg_mri-pipeline) df_conn_cohort_norm orientation exactly: legacy used
    # the lower triangle (con_node[j]/node[i] for i<j). Connectivity is
    # symmetric so values are unchanged — only the label string aligns.
    features = [f"con_{labels[j]}/{labels[i]}" for i, j in zip(*iu)]
    print(f"[{VARIANT}] {N} parcels -> {len(features)} connectivity features")

    rows, feats = [], []
    for cohort, f in found:
        m = re.search(r"(sub-[^_/]+)_(ses-[^_]+)", f.name)
        sub, ses = (m.group(1), m.group(2)) if m else (f.parts[-4], "")
        rm = re.search(r"_(run-[^_]+)", f.name)
        run = rm.group(1) if rm else "run-1"
        try:
            ts = curate(pd.read_csv(f, sep="\t")).reindex(columns=labels)
        except Exception as e:
            print(f"    skip {f.name}: {e}"); continue
        with np.errstate(divide="ignore", invalid="ignore"):
            if CONN_METHOD == "lw":
                # Ledoit-Wolf shrinkage correlation -> standard Fisher-z.
                corr = _lw_corr(ts.values)
                z = np.arctanh(np.clip(corr, -0.999999, 0.999999))
            else:
                # empirical correlation + original IQR rescale -> arctanh(±3).
                r = np.corrcoef(ts.values, rowvar=False)   # full matrix incl. diagonal
                iqr = np.nanpercentile(r, 75) - np.nanpercentile(r, 25)
                rn = (r - np.nanmedian(r)) / iqr if iqr > 0 else (r - np.nanmedian(r))
                z = np.where(np.abs(rn) >= 0.999, np.sign(rn) * 3.0, np.arctanh(rn))
        edge = z[iu]
        feats.append(edge)
        rows.append({"ID": sub, "cohort": cohort, "session": ses, "run": run,
                     "mean_fd_v11": _read_mean_fd(f),
                     "n_nan_edges": int(np.isnan(edge).sum()), "path": str(f)})

    meta = pd.DataFrame(rows)
    conn = pd.DataFrame(np.vstack(feats), columns=features)
    out = pd.concat([meta.reset_index(drop=True), conn], axis=1)
    if CONN_METHOD == "lw":
        out_csv = OUT_BASE / f"{VARIANT}_lw" / f"df_conn_raw_{VARIANT}_lw.csv"
        out_csv.parent.mkdir(parents=True, exist_ok=True)
    else:
        out_csv = CONN_RAW_CSV
    out.to_csv(out_csv, index=False)
    print(f"[{CONN_METHOD}] Wrote {out_csv}")
    print(f"  rows={len(out)}  features={len(features)}  "
          f"cohorts={meta['cohort'].value_counts().to_dict()}")
    print(f"  subject-sessions with >=1 NaN edge: {int(conn.isna().any(axis=1).sum())}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
