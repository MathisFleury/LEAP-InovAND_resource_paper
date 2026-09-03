#!/bin/bash
# Mirror every nogsr* sensitivity mode for the gsr (36P/6mm) variant.
# Chain per mode: [build raw] -> run_pipeline_v11 (ComBat/regress/z) -> analysis suite.
# Env-var semantics match run_pipeline_v11.py exactly. Fail-soft per step (logs, continues).
# ponytail: one driver over 9 modes beats 20+ hand-typed env-prefixed commands.
set -u
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SEC=$REPO/6-2_functional_analysis
PRE=$SEC/preprocessing
QC=$SEC/qc
ANA=$SEC/scripts
PY=/usr/local/bin/python3.11
export XCPD_VARIANT=gsr
cd "$PRE" || exit 1

echo "############ RAW BUILDS ############"
[ -f "$PRE/outputs/gsr/df_conn_raw_gsr.csv" ]        || $PY build_connectivity_v11.py
[ -f "$PRE/outputs/gsr_concat/df_conn_raw_gsr_concat.csv" ] || $PY build_connectivity_concat_v11.py
[ -f "$PRE/outputs/gsr_mixed/df_conn_raw_gsr_mixed.csv" ]   || $PY build_connectivity_mixed_v11.py

echo "############ PIPELINES (ComBat->regress->z) ############"
# base (gsr) already done
                      CONCAT=1                       $PY run_pipeline_v11.py
                      MIXED=1                        $PY run_pipeline_v11.py
                      REGRESS_FIRST=1                $PY run_pipeline_v11.py
                      COHORT_FILTER=INOVAND          $PY run_pipeline_v11.py
                      COHORT_FILTER=LEAP             $PY run_pipeline_v11.py
             CONCAT=1 COHORT_FILTER=INOVAND          $PY run_pipeline_v11.py
             CONCAT=1 COHORT_FILTER=LEAP             $PY run_pipeline_v11.py

echo "############ ANAT-QC-PASS (derive filtered norm from base gsr) ############"
cd "$QC" && $PY anat_qc_crossref.py
$PY - <<PYEOF
import pandas as pd, pathlib
PRE = pathlib.Path("$REPO/6-2_functional_analysis/preprocessing/outputs")
QC  = pathlib.Path("$REPO/6-2_functional_analysis/qc")
flags = pd.read_csv(QC / "anat_qc_flags_gsr.csv")
keep = set(flags.loc[flags["passed_anat_qc"] == True, "ID"].astype(str))
norm = pd.read_csv(PRE / "gsr" / "df_conn_cohort_norm_gsr.csv", low_memory=False)
norm["ID"] = norm["ID"].astype(str)
out = norm[norm["ID"].isin(keep)].reset_index(drop=True)
d = PRE / "gsr_anatQCpass"; d.mkdir(parents=True, exist_ok=True)
out.to_csv(d / "df_conn_cohort_norm_gsr_anatQCpass.csv", index=False)
print(f"[anatQCpass] kept {len(out)} / {len(norm)} subjects passing anat QC")
PYEOF

echo "############ ANALYSES (autism vs NT suite per mode) ############"
for sfx in _concat _mixed _regfirst _INOVAND _LEAP _concat_INOVAND _concat_LEAP _anatQCpass; do
  csv="$PRE/outputs/gsr${sfx}/df_conn_cohort_norm_gsr${sfx}.csv"
  out="$SEC/outputs/gsr${sfx}"
  if [ ! -f "$csv" ]; then echo "SKIP gsr${sfx} (no norm csv)"; continue; fi
  mkdir -p "$out"
  if AUTISM_TD_FMRI_CSV="$csv" AUTISM_TD_OUTPUT_DIR="$out" $PY "$ANA/run_functional_analysis.py" > "$out/_run.log" 2>&1; then
    echo "OK   gsr${sfx}"
  else
    echo "FAIL gsr${sfx}  (see $out/_run.log)"
  fi
done
echo "############ DONE ############"
