"""
LEAP + INOVAND fMRI connectivity preprocessing pipeline.

    load → QC filter (FD<1mm) → restrict ASD/Control → dedup → assign batch
         → load connectivity (Fisher z, NaN ROI handling)
         → ComBat (neuroHarmonize) → per-feature regression → z-score → save

Usage:
    cd 6_functional_analysis/non_concat/preprocessing
    /usr/local/bin/python3.11 run_pipeline.py

The output `outputs/df_conn_cohort_norm.csv` is a drop-in replacement for
the legacy file at `eeg_mri-pipeline/results/dataframes/fmri/df_conn_cohort_norm.csv`.
"""
from step_01_load             import load_dataframe
from step_02_filter_qc        import filter_qc, restrict_to_groups, keep_best_session
from step_02b_assign_batch    import assign_batch
from step_03_load_connectivity import load_connectivity
from step_04_combat           import run_combat
from step_05_regress          import fit_and_residualise
from step_06_zscore           import zscore
from step_07_save             import (
    save_zscored_csv, save_audit, save_features, save_models,
)
from config                   import APPLY_ZSCORE


def main() -> None:
    n_log = {}

    df = load_dataframe();          n_log["loaded"]              = len(df)
    df = filter_qc(df);              n_log["after_QC"]            = len(df)
    df = restrict_to_groups(df);     n_log["after_group_restrict"] = len(df)
    df = keep_best_session(df);      n_log["after_dedup"]         = len(df)
    df = assign_batch(df);           n_log["after_batch_assign"]  = len(df)

    df_conn, df_meta = load_connectivity(df)
    n_log["after_connectivity_load"] = len(df_conn)

    df_conn, combat_model = run_combat(df_conn, df_meta)
    df_conn, features, betas = fit_and_residualise(df_conn, df_meta)
    if APPLY_ZSCORE:
        df_conn = zscore(df_conn)
    n_log["final"] = len(df_conn)

    save_zscored_csv(df_conn, df_meta)
    save_audit(n_log)
    save_features(features)
    save_models(combat_model, betas)
    print("\n✓ fMRI preprocessing complete")


if __name__ == "__main__":
    main()
