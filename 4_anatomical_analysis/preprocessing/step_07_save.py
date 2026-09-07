"""
Step 7 — Save outputs.

Writes (into config.OUT_DIR, on $IMG5):
    freesurfer_zscore_qc_combat_regress.tsv   final z-scored TSV
    n_at_each_stage.tsv                       audit trail (N kept)
    features_used.txt                         feature list
    combat_model.pkl                          fitted ComBat model
    regression_betas.npy                      regression coefficients
"""
import pickle
import numpy as np
import pandas as pd

from config import OUT_DIR, ID_COLS, BATCH_COL, QC_COL, QC_TAG, EULER_COL, VARIANT_TAG


def save_zscored_tsv(df: pd.DataFrame, features: list[str]) -> None:
    meta_candidates = ID_COLS + [BATCH_COL, QC_COL, EULER_COL, "population_group"]
    meta, seen = [], set()
    for c in meta_candidates:                        # de-dupe, preserve order
        if c in df.columns and c not in seen:
            meta.append(c); seen.add(c)
    out = df[meta + features]
    fp   = OUT_DIR / f"freesurfer_zscore_{QC_TAG}_combat_regress{VARIANT_TAG}.tsv"
    out.to_csv(fp, sep="\t", index=False)
    print(f"[07] wrote  {fp}  ({len(out):,} rows x {out.shape[1]:,} cols)")


def save_audit(n_log: dict[str, int]) -> None:
    fp = OUT_DIR / f"n_at_each_stage_{QC_TAG}{VARIANT_TAG}.tsv"
    pd.DataFrame(list(n_log.items()), columns=["stage", "N"]).to_csv(
        fp, sep="\t", index=False
    )
    print(f"[07] wrote  {fp}")


def save_features(features: list[str]) -> None:
    fp = OUT_DIR / f"features_used_{QC_TAG}{VARIANT_TAG}.txt"
    fp.write_text("\n".join(features))
    print(f"[07] wrote  {fp}  ({len(features):,} features)")


def save_models(combat_model: dict, betas: np.ndarray) -> None:
    with open(OUT_DIR / f"combat_model_{QC_TAG}{VARIANT_TAG}.pkl", "wb") as f:
        pickle.dump(combat_model, f)
    np.save(OUT_DIR / f"regression_betas_{QC_TAG}{VARIANT_TAG}.npy", betas)
    print(f"[07] wrote  combat_model_{QC_TAG}{VARIANT_TAG}.pkl  +  "
          f"regression_betas_{QC_TAG}{VARIANT_TAG}.npy")


if __name__ == "__main__":
    print("Run via run_pipeline.py — depends on full pipeline output.")
