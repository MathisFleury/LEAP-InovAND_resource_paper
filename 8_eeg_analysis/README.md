# 8. EEG Analysis

Resting-state EEG, whole-group autism vs NT: alpha peak frequency
(`01_eeg_alpha_peak_autism_td.py`) and multi-band relative power
(`02_eeg_power_bands_autism_td.py`).

```bash
cd 8_eeg_analysis/scripts
python3.11 run_eeg_analysis.py
```

## Data

- **Alpha peak** — rebuilt from curated raw features on `$IMG5`
  (`preprocessing/build_alpha_peak_corrected.py`, regressed + z-scored on
  the full sample), written to
  `preprocessing/outputs/Alpha_peak_combined_corrected.csv`. This is a
  *different, newer file that happens to share a filename* with the frozen
  `$LIB/eeg_mri-pipeline/results/dataset_paper/dataframes/Alpha_peak_combined_corrected.csv`
  — don't confuse the two.
- **Power spectrum** — no raw LEAP extraction exists on `$IMG5`, so this
  stays on the frozen
  `$LIB/eeg_mri-pipeline/results/dataset_paper/dataframes/Power_spectrum_combined_corrected.csv`,
  with curated demographics/cluster labels grafted on via
  `preprocessing/curated_map.py`.

`9_cluster_eeg_analysis` shares this section's `preprocessing/` outright
(no separate rebuild).

## Outputs

`outputs/tables/` — all CSVs (`stats/`, `fdr_corrected/`,
`covariate_effects/`, plus the top-level alpha-peak stats CSV).
`outputs/figures/` — all plots (the two alpha-peak violin PDFs, `plots/`
for the ~260 per-band/per-region power plots).
