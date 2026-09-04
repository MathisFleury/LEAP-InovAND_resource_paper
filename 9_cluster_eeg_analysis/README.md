# 9. Cluster EEG Analysis

Per-cluster analogue of `8_eeg_analysis`: alpha peak frequency and
multi-band relative power, autism clusters (C1/C2/C3) vs pooled NT.

```bash
cd 9_cluster_eeg_analysis/scripts
python3.11 run_cluster_eeg_analysis.py
```

## Data

Shares `8_eeg_analysis/preprocessing/` outright — same alpha-peak rebuild
and frozen power-spectrum table, no separate section-9 preprocessing. See
`8_eeg_analysis/README.md` for the data source details.

## Outputs

`outputs/tables/` — all CSVs (`stats/`, `fdr_corrected/`, plus the
top-level alpha-peak stats CSV). `outputs/figures/` — all plots
(`plots/` for per-band figures, the alpha-peak violin PDF, the power
synthesis heatmap).

An earlier, frozen paper-vintage run and a stale orphaned `figures/curated/`
subtree (no script writes there anymore) live in
`../legacy/9_cluster_eeg_analysis/` (local only, not part of this public
release).
