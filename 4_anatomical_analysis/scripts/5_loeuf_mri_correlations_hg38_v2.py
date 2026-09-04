#!/usr/bin/env python3
"""
LOEUF (hg38) × MRI correlation statistics — REVISION

Produces the stats CSV consumed by 6_loeuf_mri_brain_maps_hg38_v2.R.
Mirrors the hg38 path of ../../4_anatomical_analysis/scripts/11_loeuf_mri_interaction.py
but uses the curated MRI table (QC + ComBat + age/sex/eTIV-regression z-scored)
and the curated-pipeline loader/wave-selection logic from 1_anatomical_mri_autism_nt_v2.py.

Input MRI : the curated MRI table (resolved from 1_anatomical_mri_autism_nt_v2)
Input gen : df_carrier_genelist_DEL_LOF_MISS_withalphamissense_LOEUF_gnomadv4.tsv
Joining   : curated MRI → df_clusters_complete (canonical join keys) →
            ID → merge on carrier_annotations_hg38.tsv

Output:
  ../outputs/tables_genetics_hg38/loeuf_mri_interaction_statistics_hg38.csv
"""

from pathlib import Path
import re
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from scipy import stats as scipy_stats
from statsmodels.stats.multitest import multipletests

# Reuse loader / wave-selection / canonical-key logic from script 01_v2.
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module
_mod01 = import_module("1_anatomical_mri_autism_nt_v2")  # noqa: E402
import _config  # noqa: E402  shared paths + hyperparameters (see _config.py)

# =============================================================================
# PATHS / PARAMS (centralised in _config.py)
# =============================================================================
_SCRIPT_DIR  = Path(__file__).resolve().parent
_SECTION_DIR = _SCRIPT_DIR.parent                      # 4_anatomical_analysis/

CARRIER_HG38   = Path(_config.CARRIER_V4)       # gnomAD v4 carrier table
GENETIC_ROSTER = Path(_config.GENETIC_ROSTER)   # genetics<->MRI bridge (paper IDs)
TABLES_DIR     = Path(_config.OUTPUT_BASE) / "tables_genetics_hg38"
AUTISM_ONLY    = _config.AUTISM_ONLY            # restrict to autistic carriers

# =============================================================================
# Gene lists — all-carriers only.
# The constrained-only subset is intentionally dropped: running both subsets
# in the same FDR pool inflated the q-values.
# =============================================================================
GENE_LISTS_HG38_ALL = [
    ("dellof_syngo_best_score_hg38",      "SynGO delloF", "#4A90E2"),
    ("dellof_chromepitf_best_score_hg38", "CHROM delloF", "#8E44AD"),
]

# Subcortical ROIs: reuse the canonical set from 1_anatomical_mri_autism_nt_v2
# (ASEG_REGIONS). Labels in the CSV match those produced by 01 (e.g.
# "Left-Thalamus", "CC_Anterior") so the R script can reuse the existing
# rename_subcortical_labels() helper from 2_plot_anatomical_mri_brain_v2.R.
SUBCORTICAL_COLS = list(_mod01.ASEG_REGIONS)

# (new_col, lof_col, del_col) — min ignoring NaN
# v4 column names: lof/del per-gene-list LOEUF best scores.
COMBINE_HG38 = [
    ("dellof_syngo_best_score_hg38",
     "lof_syngo_genes_best_score",          "del_syngo_genes_best_score"),
    ("dellof_chromepitf_best_score_hg38",
     "lof_chromepitf_genes_best_score",     "del_chromepitf_genes_best_score"),
]


def neg_log10(series: pd.Series) -> pd.Series:
    s = series.dropna()
    return -np.log10(s[s > 0])


def load_curated_mri_with_canonical_id() -> pd.DataFrame:
    """Run the curated-pipeline wave-selection pipeline and attach canonical ID
    (the one used in df_clusters_complete and carrier_annotations_hg38)."""
    print(f"Loading curated MRI: {_mod01.MRI_FILE}")
    df = pd.read_csv(_mod01.MRI_FILE, sep="\t", low_memory=False)
    df["ID"] = df["ID"].astype(str)
    df = _mod01._rename_freesurfer_cols(df)

    df = _mod01.select_one_row_per_subject(df)
    df["_join_key"] = df.apply(_mod01._canonical_join_key_mri, axis=1)
    df = df.dropna(subset=["_join_key"])

    # Bridge through the genetic-aligned roster (IDs match the v4 carrier file),
    # NOT the curated clinical roster — see GENETIC_ROSTER note above.
    clu = pd.read_csv(GENETIC_ROSTER, low_memory=False)
    clu["_join_key"] = clu.apply(_mod01._canonical_join_key_clusters, axis=1)
    pheno = "population_group" if "population_group" in clu.columns else "PopulationS1"
    clu = (clu[["ID", "_join_key", "cohort", pheno]]
           .rename(columns={pheno: "PopulationS1"})
           .dropna(subset=["_join_key"]))
    clu["ID"] = clu["ID"].astype(str)

    df = df.drop(columns=[c for c in ["ID"] if c in df.columns])
    merged = clu.merge(df, on="_join_key", how="inner", suffixes=("", "_mri"))

    print(f"  curated MRI ready: {len(merged)} rows (canonical ID attached)")
    return merged


def load_carriers_hg38() -> pd.DataFrame:
    print(f"Loading carriers (hg38): {CARRIER_HG38}")
    g = pd.read_table(CARRIER_HG38, low_memory=False)
    g["ID"] = g["ID"].astype(str)
    # v4 already uses IDD/NT spelling (CLAUDE.md) — no relabel. Prefer the
    # `population_group` grouping when present, for consistency with the curated
    # MRI labels. (LOEUF correlations pool all carriers, so this does not change
    # the result — it only keeps the label convention consistent.)
    if "population_group" in g.columns:
        g["PopulationS1"] = g["population_group"]
    g = g[g["PopulationS1"] != "other"].drop_duplicates("ID")
    print(f"  carriers: {len(g)} rows")
    return g


def build_dataset() -> tuple[pd.DataFrame, list[str]]:
    mri  = load_curated_mri_with_canonical_id()
    gens = load_carriers_hg38()

    # Keep the v4 (genetics) PopulationS1 — it already folds the autism subtypes
    # (with/without IDD, "to exclude") into "Autism". Drop the MRI-side copies to
    # avoid a clash on merge.
    mri = mri.drop(columns=[c for c in ("cohort", "PopulationS1") if c in mri.columns])
    df = gens.merge(mri, on="ID", how="inner")
    print(f"  merged: {len(df)} individuals")

    # Restrict to autistic carriers only (see AUTISM_ONLY).
    if AUTISM_ONLY:
        before = len(df)
        df = df[df["PopulationS1"] == "Autism"].copy()
        print(f"  AUTISM-ONLY filter: {before} -> {len(df)} autistic individuals")


    # Combine lof + del → dellof (min ignoring NaN), same as script 11.
    for new_col, lof_col, del_col in COMBINE_HG38:
        lof_s = df[lof_col] if lof_col in df.columns else pd.Series(np.nan, index=df.index)
        del_s = df[del_col] if del_col in df.columns else pd.Series(np.nan, index=df.index)
        df[new_col] = np.fmin(lof_s.values, del_s.values)

    ct_cols = [c for c in df.columns
               if (c.endswith("_thickness") or c.endswith("_area")
                   or c.endswith("_grayvol"))
               and (c.startswith("lh_") or c.startswith("rh_"))
               and "MeanThickness" not in c and "WhiteSurfArea" not in c]
    sc_cols = [c for c in SUBCORTICAL_COLS if c in df.columns]
    missing_sc = [c for c in SUBCORTICAL_COLS if c not in df.columns]
    if missing_sc:
        print(f"  subcortical missing (skipped): {missing_sc}")
    print(f"  cortical columns: {len(ct_cols)} · subcortical columns: {len(sc_cols)}")
    return df, ct_cols, sc_cols


def run_correlations(df: pd.DataFrame, gene_lists, ct_cols, sc_cols, label: str) -> pd.DataFrame:
    print(f"\nCorrelations: {label}")
    rows = []
    for col, gl_label, _ in gene_lists:
        if col not in df.columns:
            print(f"  skip {col}: missing")
            continue
        x_full = neg_log10(df[col])
        n_carriers = len(x_full)
        for mri_col in list(ct_cols) + list(sc_cols):
            mri = df[mri_col].dropna()
            common = x_full.index.intersection(mri.index)
            if len(common) < 10:
                continue
            xz = scipy_stats.zscore(x_full.loc[common].values)
            yz = scipy_stats.zscore(mri.loc[common].values)
            r, p = pearsonr(xz, yz)
            m = re.match(r"^(lh|rh)_(.+)_(thickness|area|grayvol)$", mri_col)
            if m:
                hemi     = "left" if m.group(1) == "lh" else "right"
                region   = m.group(2)
                mri_type = m.group(3)
            else:
                # Subcortical: leave hemisphere/region empty — the R side
                # rebuilds the ggseg label from `mri_col` directly (lowercased)
                # via rename_subcortical_labels(), which handles midline
                # structures (Brain-Stem, CC_*, 3rd/4th-Ventricle) too.
                if mri_col.startswith("Left-"):
                    hemi = "left"
                elif mri_col.startswith("Right-"):
                    hemi = "right"
                else:
                    hemi = ""
                region   = mri_col
                mri_type = "subcortical"
            rows.append({
                "genetic_type":      "LOEUF",
                "genetic_feature":   col,
                "gene_list_label":   gl_label,
                "mri_col":           mri_col,
                "mri_type":          mri_type,
                "hemisphere":        hemi,
                "region":            region,
                "sample_size":       int(len(common)),
                "n_total_carriers":  int(n_carriers),
                "correlation":       float(r),
                "p_value":           float(p),
                "p_permutation":     np.nan,
            })
        print(f"  {gl_label}: {n_carriers} carriers")

    df_s = pd.DataFrame(rows)
    if df_s.empty:
        return df_s
    # FDR PER PANEL (one ggseg facet = one genetic_feature × one mri_type).
    # Matches the pooling used by 6_loeuf_mri_brain_maps_hg38_v2.R so the
    # `p_fdr` column in this CSV equals what's plotted.
    df_s["p_fdr"] = np.nan
    for _, sub in df_s.groupby(["genetic_feature", "mri_type"]):
        valid = sub["p_value"].notna()
        if valid.sum() == 0:
            continue
        _, pfdr, _, _ = multipletests(sub.loc[valid, "p_value"].values, method="fdr_bh")
        df_s.loc[sub.index[valid], "p_fdr"] = pfdr
    return df_s


# Short pathway tokens for the per-file names (genetic_feature -> token).
PATHWAY_SHORT = {
    "dellof_syngo_best_score_hg38":      "syngo",
    "dellof_chromepitf_best_score_hg38": "chromepitf",
}


def main() -> None:
    TABLES_DIR.mkdir(parents=True, exist_ok=True)

    df, ct_cols, sc_cols = build_dataset()
    stats = run_correlations(df, GENE_LISTS_HG38_ALL, ct_cols, sc_cols, "hg38 all carriers")
    if stats.empty:
        print("No stats produced — nothing written.")
        return

    # One r-input file per pathway × MRI feature. run_correlations already
    # computes FDR within each (genetic_feature, mri_type) group, so each file
    # carries its own independent FDR correction over its ROIs only.
    print("\nWriting per-pathway × per-feature r-input files:")
    for (feat, mtype), sub in stats.groupby(["genetic_feature", "mri_type"]):
        short = PATHWAY_SHORT.get(feat, feat)
        out = TABLES_DIR / f"loeuf_mri_{short}_{mtype}.csv"
        sub.to_csv(out, index=False)
        n_sig = int((sub["p_fdr"] < 0.05).sum())
        print(f"  {out.name}: {len(sub)} ROIs · FDR<0.05 {n_sig}")


if __name__ == "__main__":
    main()
