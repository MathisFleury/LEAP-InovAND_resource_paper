"""
Step 6 — Z-score each EEG feature across the FULL sample.

Reference = full sample (DX_COL=None): we no longer normalise to the NT /
control group. Mirrors the anat pipeline step_06_zscore with DX_COL=None.
NaNs are skipped per column (pandas mean/std default).
"""
import numpy as np


def zscore(df, features):
    df = df.copy()
    mu = df[features].mean()
    sd = df[features].std(ddof=1).replace(0, np.nan)
    df[features] = (df[features] - mu) / sd
    print(f"[06] z-scored {len(features)} feature(s) across N={len(df)} (full sample)")
    return df


if __name__ == "__main__":
    # ponytail self-check: z-scored column has mean 0, sd 1
    import pandas as pd
    demo = pd.DataFrame({"feat": np.arange(100, dtype=float)})
    out = zscore(demo, ["feat"])
    assert abs(out["feat"].mean()) < 1e-9
    assert abs(out["feat"].std(ddof=1) - 1.0) < 1e-9
    print("ok")
