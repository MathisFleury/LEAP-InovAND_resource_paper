"""
NDD-focused anatomical z-scoring pipeline: INFOR + INOVAND + LEAP.

    load -> restrict to NDD -> QC filter -> euler -> assign batch (scanner)
         -> ComBat -> regress -> z-score -> save

Ported in-repo (2026-09) from eeg_mri-pipeline's pipeline_ndd, so the raw-data
build lives next to the analysis scripts that consume its output, matching
6_functional_analysis/*/preprocessing and 8_eeg_analysis/preprocessing.

Usage:
    /usr/local/bin/python3.11 4_anatomical_analysis/preprocessing/run_pipeline.py
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))          # this pipeline's own config + step modules

from step_01_load            import load_dataframe
from step_02_filter_qc       import (
    restrict_to_ndd, filter_qc, add_euler,
    drop_missing_covariates, drop_tiny_cohorts, keep_best_run,
)
from step_02b_assign_batch   import assign_batch
from step_03_select_features import (
    select_features, restrict_to_numeric, drop_subjects_missing_features,
    drop_zero_variance_per_batch,
)
from step_04_combat          import run_combat
from step_05_regress         import fit_and_residualise
from step_06_zscore          import zscore
from step_07_save            import (
    save_zscored_tsv, save_audit, save_features, save_models,
)


def main() -> None:
    n_log = {}

    # ---------- 1. load ------------------------------------------------
    df = load_dataframe()
    n_log["loaded"] = len(df)

    # ---------- 2. restrict + QC + euler -------------------------------
    df = restrict_to_ndd(df);          n_log["after_NDD_filter"]    = len(df)
    df = filter_qc(df);                n_log["after_QC"]            = len(df)
    df = add_euler(df)

    # ---------- 2b. assign scanner-level batch + fill demographics -----
    df = assign_batch(df)

    # ---------- 2c. dedup to best run, drop missing, drop tiny batches
    df = keep_best_run(df);            n_log["after_best_run"]          = len(df)
    df = drop_missing_covariates(df);  n_log["after_drop_missing"]      = len(df)
    df = drop_tiny_cohorts(df);        n_log["after_drop_tiny_batches"] = len(df)

    # ---------- 3. select FreeSurfer features --------------------------
    features = select_features(df)
    features = restrict_to_numeric(df, features)
    df = drop_subjects_missing_features(df, features)
    n_log["after_feature_complete"] = len(df)
    features = drop_zero_variance_per_batch(df, features)

    # ---------- 4. ComBat ----------------------------------------------
    df, combat_model = run_combat(df, features)

    # ---------- 5. regress nuisance covariates (per-feature loop) ------
    df, features, betas = fit_and_residualise(df, features)

    # ---------- 6. z-score ---------------------------------------------
    df = zscore(df, features);         n_log["zscored"]             = len(df)

    # ---------- 7. save ------------------------------------------------
    save_zscored_tsv(df, features)
    save_audit(n_log)
    save_features(features)
    save_models(combat_model, betas)

    print("\nNDD anatomical preprocessing pipeline complete")


if __name__ == "__main__":
    main()
