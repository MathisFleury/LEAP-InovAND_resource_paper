"""
Step 3 — Load per-subject connectivity matrices, normalize, and flatten.

For each surviving row in `df`:
    1. Read the XCP-D-produced parcellated correlation matrix from
       `df.loc[i, 'path']`. Files are TSVs indexed by 'Node'.
    2. Skip subjects whose first row/column carries any NaN (corrupt XCP-D
       output) — these IDs are written to bad_subjects.txt.
    3. (Optional) Re-centre by median and rescale by inter-quartile range.
    4. (Optional) Apply Fisher z = arctanh(r), clipping |r| ≥ 0.999 to ±3.
    5. Detect ROI rows/columns that are NaN for any surviving subject and
       drop them from every subject's matrix (so every individual ends
       with the same ROI grid).
    6. Extract the lower triangle (excluding diagonal) and pivot to a
       wide table indexed by ID, columns named `con_{source}/{target}`.

Returns:
    df_conn         wide DataFrame (n_subjects × n_features) of connectivity
                    values, ID-indexed.
    df_aligned      subset of the input df, in the same order/index as df_conn
                    (subjects with corrupt matrices removed).
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from config import (
    APPLY_IQR_RESCALE, APPLY_FISHER_Z,
    OUT_BAD_SUBJECTS, OUT_NAN_MASK,
)


def _normalize_subject_matrix(mat: np.ndarray) -> np.ndarray:
    out = mat.copy()
    if APPLY_IQR_RESCALE:
        med = np.nanmedian(out)
        iqr = np.nanpercentile(out, 75) - np.nanpercentile(out, 25)
        if iqr > 0:
            out = (out - med) / iqr
    if APPLY_FISHER_Z:
        out = np.where(
            np.abs(out) >= 0.999,
            np.sign(out) * 3.0,
            np.arctanh(np.clip(out, -0.999, 0.999)),
        )
    return out


def _load_subject_matrix(path: str) -> tuple[np.ndarray | None, pd.DataFrame | None, str | None]:
    """Return (matrix, raw_df_for_index, reject_reason)."""
    df_sub = pd.read_table(path, index_col="Node")
    mat = df_sub.to_numpy(dtype=float)
    border_nan = (
        np.isnan(mat[1, :]).sum() + np.isnan(mat[:, 1]).sum()
    )
    if border_nan > 0:
        return None, None, "border_nan"
    if np.isnan(mat).all():
        return None, None, "all_nan"
    return _normalize_subject_matrix(mat), df_sub, None


def load_connectivity(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    matrices: list[np.ndarray] = []
    kept_idx: list[int] = []
    bad_ids: list[str] = []
    reference_df: pd.DataFrame | None = None

    for i, row in df.reset_index(drop=True).iterrows():
        mat, ref, reason = _load_subject_matrix(row["path"])
        if mat is None:
            bad_ids.append(str(row["ID"]))
            continue
        if reference_df is None:
            reference_df = ref
        matrices.append(mat)
        kept_idx.append(i)

    OUT_BAD_SUBJECTS.write_text("\n".join(bad_ids) + ("\n" if bad_ids else ""))
    print(f"[03] connectivity loaded  →  kept {len(matrices):,}, "
          f"rejected {len(bad_ids):,} (bad_subjects.txt written)")

    if not matrices:
        raise RuntimeError("no surviving connectivity matrices — check inputs")

    cube = np.array(matrices)                              # (n_sub, n_roi, n_roi)
    nan_rows = np.any(np.isnan(cube), axis=0)              # (n_roi, n_roi)
    drop_rows = nan_rows.any(axis=1)
    drop_cols = nan_rows.any(axis=0)
    np.save(OUT_NAN_MASK, nan_rows)

    # Drop NaN rows/cols (they are identical via symmetry, but be defensive).
    cube = cube[:, ~drop_rows, :][:, :, ~drop_cols]
    if reference_df is not None:
        clean_index   = reference_df.index[~drop_rows]
        clean_columns = reference_df.columns[~drop_cols]
    else:
        clean_index   = [f"ROI_{i}" for i in range(cube.shape[1])]
        clean_columns = [f"ROI_{i}" for i in range(cube.shape[2])]

    df_aligned = df.reset_index(drop=True).iloc[kept_idx].reset_index(drop=True)

    # Lower triangle (k=-1, excluding diagonal) flattened to long, then pivot.
    tril_mask = np.tril(np.ones(cube.shape[1:], dtype=bool), k=-1)
    feat_names = [
        f"con_{clean_index[r]}/{clean_columns[c]}"
        for r, c in zip(*np.where(tril_mask))
    ]
    n_feat = len(feat_names)
    flat = np.empty((cube.shape[0], n_feat), dtype=float)
    rr, cc = np.where(tril_mask)
    for k, (r, c) in enumerate(zip(rr, cc)):
        flat[:, k] = cube[:, r, c]

    df_conn = pd.DataFrame(flat, columns=feat_names)
    df_conn.insert(0, "ID", df_aligned["ID"].astype(str).values)
    df_conn = df_conn.drop_duplicates(subset="ID").set_index("ID")

    df_aligned = df_aligned[df_aligned["ID"].astype(str).isin(df_conn.index)].copy()
    df_aligned = df_aligned.drop_duplicates(subset="ID").set_index("ID").loc[df_conn.index]
    df_aligned = df_aligned.reset_index()

    print(f"[03] flattened  →  {df_conn.shape[0]:,} subjects × {df_conn.shape[1]:,} features")
    return df_conn, df_aligned


if __name__ == "__main__":
    from step_01_load import load_dataframe
    from step_02_filter_qc import filter_qc, restrict_to_groups, keep_best_session
    from step_02b_assign_batch import assign_batch

    df = assign_batch(keep_best_session(restrict_to_groups(filter_qc(load_dataframe()))))
    df_conn, df_aligned = load_connectivity(df)
    print(df_conn.head().iloc[:, :3])
