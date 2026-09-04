#!/usr/bin/env python3
"""
Euler-derived defects vs. ALL anatomical phenotypes (not just the 36
FDR-significant-for-group-difference ones) -- ggseg brain-map version of
R4.5's image-quality x phenotype effect-size request, per user request
(2026-08-31): a whole-brain figure instead of the 36-row table/forest plot.

For every one of the 232 features (68 thickness, 68 area, 68 grayvol,
28 subcortical volume), correlates defects = -mean_euler against the
phenotype, pooled (whole Autism+NT sample, matching 17_euler_confound_check
.py's "pooled" scope) -- reported as an EQUIVALENT COHEN'S D
(d = 2r / sqrt(1 - r^2), the standard r-to-d conversion), not as r itself,
per user request 2026-08-31: the reviewer specifically asked for an "effect
size", and d is what every other effect size in this manuscript is expressed
in, so reporting this confound check in the same units lets it be compared
directly against the reported group-difference effect sizes (|d| ~ 0.2-0.6)
rather than living on a different (correlation) scale.

Outputs ggseg-input CSVs in the same 5-column format
2_plot_anatomical_mri_brain_v2.R already reads (label, t_stat, p_val,
cohens_d, p_fdr) -- `cohens_d` holds this converted value (so ggseg colors by
it) and `p_fdr` holds the ORIGINAL group-difference BONFERRONI-corrected p
(`p_bonf`, not `p_fdr` -- confirmed via R4.4 of the reviewer response: the
main-text whole-cohort Results are Bonferroni-corrected; Supplementary Table
8's own FDR correction is a separate, more lenient criterion used only in
that table, and outlines a much larger, different set of regions -- 37 FDR
vs 10 Bonferroni across all four metrics. The column is still named `p_fdr`
here only because that's the field name the shared R plotting code expects;
its contents are the Bonferroni p) -- so the outline marks exactly the
regions the main-text Results paragraph reports as significant, not a
different, FDR-only set.

Output: outputs/figures/r_input_files_qc/euler_r_vs_pheno_<atlas>_<metric>.csv
"""
import sys
from pathlib import Path
from importlib import import_module

import numpy as np
import pandas as pd
from scipy.stats import pearsonr

sys.path.insert(0, str(Path(__file__).resolve().parent))
_m01 = import_module("1_anatomical_mri_autism_nt_v2")

_SECTION = Path(__file__).resolve().parent.parent
RINPUT = _SECTION / "outputs" / "figures" / "r_input_files"
OUT = _SECTION / "outputs" / "figures" / "r_input_files_qc"
OUT.mkdir(parents=True, exist_ok=True)

METRICS = [("thickness", "dk",    "t_stat_anat_dk_thickness_mri_autism_vs_control.csv"),
           ("area",      "dk",    "t_stat_anat_dk_area_mri_autism_vs_control.csv"),
           ("grayvol",   "dk",    "t_stat_anat_dk_grayvol_mri_autism_vs_control.csv"),
           ("volume",    "aseg",  "t_stat_anat_aseg_volume_mri_autism_vs_control.csv")]


def main():
    df = _m01.load_data()
    if "mean_euler" not in df.columns:
        raise SystemExit("mean_euler not in loaded data")
    df["defects"] = -df["mean_euler"]
    df = df[df["PopulationS1"].isin(["Autism", "NT"])]

    for metric, atlas, fname in METRICS:
        fp = RINPUT / fname
        if not fp.exists():
            print(f"  {metric}: {fname} missing -- skip"); continue
        group = pd.read_csv(fp)  # per-region GROUP-DIFFERENCE result: label,t_stat,p_val,cohens_d,p_fdr

        rows = []
        for _, g in group.iterrows():
            col = g["label"] if metric == "volume" else f"{g['label']}_{metric}"
            if col not in df.columns:
                continue
            x = df["defects"].values
            y = df[col].values
            ok = np.isfinite(x) & np.isfinite(y)
            if ok.sum() < 10:
                continue
            r, p = pearsonr(x[ok], y[ok])
            d_equiv = 2 * r / np.sqrt(1 - r ** 2)
            rows.append({"label": g["label"], "t_stat": np.nan, "p_val": round(p, 4),
                         "cohens_d": round(d_equiv, 4), "r": round(r, 4),
                         "p_fdr": g["p_bonf"], "n": int(ok.sum())})
        out = pd.DataFrame(rows)
        out.to_csv(OUT / f"euler_r_vs_pheno_{atlas}_{metric}.csv", index=False)
        n_group_sig = int((out["p_fdr"] < 0.05).sum())
        print(f"  {metric:9s} ({atlas}): {len(out)} features, "
              f"mean|d_equiv|={out['cohens_d'].abs().mean():.3f} "
              f"(mean|r|={out['r'].abs().mean():.3f})  "
              f"Bonferroni-significant group effect (outlined): {n_group_sig}")
    print(f"\nSaved: {OUT}/euler_r_vs_pheno_<atlas>_<metric>.csv")


if __name__ == "__main__":
    main()
