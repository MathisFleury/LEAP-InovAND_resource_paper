#!/usr/bin/env python3
# =============================================================================
# 01b - hg38 Carrier Annotation from Raw Data
# =============================================================================
# Builds carrier_annotations_hg38.tsv from:
#   - df_38_with_lof_carrier_status.tsv  (individual-level LoF carrier flags)
#   - df_variants_GRCh38.tsv             (variant-level, for metabolic LoF + lof_carrier)
#   - SV_DEL_annotation_12juin2025_ZM.tsv (DEL, same as hg19)
#
# The LoF carrier file already contains carrier flags for standard gene lists
# (hcndddom, hcnddxlinked, sparksfari1, eagle, syngo, chromepitf, constraint).
# This script adds:
#   - lof_carrier (overall unconstrained flag)
#   - Metabolic gene list annotations (LoF + DEL)
#   - Combined HCNDD dominant + X-linked (males only) flags
#   - DEL carrier flags + gene names + LOEUF scores
#   - dellof merged carrier flags
#
# Note: hg38 missense variant file not yet available. dellofmiss columns are
# not produced here. Script 04 gracefully skips missing columns.
#
# Output: tables_hg38/carrier_annotations_hg38.tsv
# =============================================================================

import os
import sys
import pandas as pd
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _config import (
    TABLES_DIR_HG38, GENELIST_PATH, INDIVIDUALS_METRICS, LOEUF_SCORES,
    DEL_FILE, MISSING_DENOVO_CNV, LOF_HG38_FILE, VARIANTS_HG38_FILE,
)

# Import helpers from script 01
from importlib.util import spec_from_file_location, module_from_spec
_script_dir = os.path.dirname(os.path.abspath(__file__))
_spec = spec_from_file_location("script01", os.path.join(_script_dir, "01_carrier_annotation.py"))
_mod = module_from_spec(_spec)
_spec.loader.exec_module(_mod)

load_gene_lists = _mod.load_gene_lists
load_loeuf_scores = _mod.load_loeuf_scores
annotate_carrier_flags = _mod.annotate_carrier_flags
add_xlinked_dominant_combined = _mod.add_xlinked_dominant_combined
extract_del_gene_names = _mod.extract_del_gene_names
assign_best_gene_scores = _mod.assign_best_gene_scores
combine_xlinked_scores = _mod.combine_xlinked_scores

os.makedirs(TABLES_DIR_HG38, exist_ok=True)


# ============================================================================
# hg38-specific helpers
# ============================================================================

def merge_carrier_flags_ext(df, prefix1, prefix2, out_prefix):
    """Merge carrier boolean flags (OR logic), including metabolic."""
    suffixes = [
        "_carrier",
        "_hcndddomv6_contraint_carrier", "_hcnddxlinked_contraint_carrier",
        "_hcndddom_xlinked_boyz_contraint_carrier",
        "_sparksfari1_contraint_carrier", "_eagle_contraint_carrier",
        "_syngo_contraint_carrier", "_chromepitf_contraint_carrier",
        "_syngo_chromepitf_contraint_carrier", "_contraint_carrier",
        "_metabolic_contraint_carrier", "_metabolic_carrier",
        "_hcndddomv6_carrier", "_hcnddxlinked_carrier",
        "_hcndddom_xlinked_boyz_carrier",
        "_sparksfari1_carrier", "_eagle_carrier",
        "_syngo_carrier", "_chromepitf_carrier",
        "_syngo_chromepitf_carrier",
    ]
    for s in suffixes:
        c1 = f"{prefix1}{s}"
        c2 = f"{prefix2}{s}"
        if c1 in df.columns and c2 in df.columns:
            df[f"{out_prefix}{s}"] = df[c1] | df[c2]
    return df


def extract_metabolic_del_genes(df_del, metabolic_genes, metabolic_constrained):
    """Extract metabolic gene names from DEL Gene_name, grouped per sample."""
    met_set = set(metabolic_genes)
    met_const_set = set(metabolic_constrained)

    df_temp = df_del[["SAMPLE", "Gene_name"]].copy()
    df_temp["del_metabolic_genes"] = df_del["Gene_name"].apply(
        lambda x: ";".join([g for g in str(x).split(";") if g in met_set]) or np.nan
    )
    df_temp["del_metabolic_constraint_genes"] = df_del["Gene_name"].apply(
        lambda x: ";".join([g for g in str(x).split(";") if g in met_const_set]) or np.nan
    )

    def _merge(series):
        non_na = series.dropna()
        if non_na.empty:
            return np.nan
        genes = []
        for item in non_na:
            genes.extend(item.split(";"))
        unique = sorted(set(g.strip() for g in genes if g.strip()))
        return ";".join(unique) if unique else np.nan

    merged = df_temp.groupby("SAMPLE", as_index=False).agg({
        "del_metabolic_genes": _merge,
        "del_metabolic_constraint_genes": _merge,
    })
    merged = merged.rename(columns={"SAMPLE": "barcode"})
    return merged


# ============================================================================
# Main pipeline
# ============================================================================

def main():
    # --- 1. Load base LoF carrier file ---
    print(f"Loading hg38 LoF carrier file: {LOF_HG38_FILE}")
    if not os.path.exists(LOF_HG38_FILE):
        print(f"Error: LoF file not found at {LOF_HG38_FILE}")
        print("Make sure the external drive is mounted.")
        return
    df = pd.read_table(LOF_HG38_FILE, low_memory=False)
    print(f"  {len(df)} individuals, {len(df.columns)} columns")

    # --- 2. Add FID from clinical data ---
    print(f"Loading FID from: {INDIVIDUALS_METRICS}")
    clinic = pd.read_table(INDIVIDUALS_METRICS)
    fid_map = clinic[["ID", "FID"]].drop_duplicates("ID")
    df = df.merge(fid_map, how="left", on="ID")

    # --- 3. Load gene lists (including metabolic) ---
    print(f"Loading gene lists: {GENELIST_PATH}")
    gene_lists, constrained_lists = load_gene_lists(GENELIST_PATH)

    # Add metabolic gene list (not in standard load_gene_lists)
    genetrek = pd.read_table(GENELIST_PATH)
    gene_lists["metabolic"] = genetrek[
        genetrek["Metabolic_process_GO_8152"] == True
    ]["symbol"].tolist()
    constrained_lists["metabolic"] = list(
        set(gene_lists["metabolic"]) & set(gene_lists["constraint"])
    )
    print(f"  Metabolic genes: {len(gene_lists['metabolic'])} total, "
          f"{len(constrained_lists['metabolic'])} constrained")

    # --- 4. Load LOEUF scores ---
    print(f"Loading LOEUF scores: {LOEUF_SCORES}")
    loeuf_dict = load_loeuf_scores(LOEUF_SCORES)

    # ================================================================
    # 5. LoF annotation from variants file
    # ================================================================
    print(f"\nLoading hg38 LoF variants: {VARIANTS_HG38_FILE}")
    if not os.path.exists(VARIANTS_HG38_FILE):
        print(f"Error: Variants file not found at {VARIANTS_HG38_FILE}")
        return
    variants = pd.read_table(VARIANTS_HG38_FILE, low_memory=False)
    print(f"  {len(variants)} variants loaded")

    # Filter: HC LoF, no flags, not in LCR
    # Note: file is pre-filtered to HC canonical LoF; VEP_LoF_flags uses "." for none
    lof_vars = variants[
        (variants["VEP_LoF"] == "HC")
        & ((variants["VEP_LoF_flags"].isna()) | (variants["VEP_LoF_flags"] == "."))
        & (variants["LCR"].astype(str) != "True")
    ].copy()
    print(f"  {len(lof_vars)} LoF variants after filtering")

    # 5a. lof_carrier flag (any LoF variant)
    df["lof_carrier"] = df["barcode"].isin(lof_vars["sample"].unique())
    print(f"  lof_carrier: {df['lof_carrier'].sum()} individuals")

    # 5b. Metabolic LoF carriers
    metabolic_set = set(gene_lists["metabolic"])
    metabolic_constraint_set = set(constrained_lists["metabolic"])

    lof_met_all = lof_vars[lof_vars["VEP_SYMBOL"].isin(metabolic_set)]
    lof_met_const = lof_vars[lof_vars["VEP_SYMBOL"].isin(metabolic_constraint_set)]

    df["lof_metabolic_carrier"] = df["barcode"].isin(lof_met_all["sample"].unique())
    df["lof_metabolic_contraint_carrier"] = df["barcode"].isin(lof_met_const["sample"].unique())
    print(f"  lof_metabolic_carrier: {df['lof_metabolic_carrier'].sum()}")
    print(f"  lof_metabolic_contraint_carrier: {df['lof_metabolic_contraint_carrier'].sum()}")

    # Metabolic LoF gene names
    if len(lof_met_all) > 0:
        met_genes = (
            lof_met_all.groupby("sample")["VEP_SYMBOL"]
            .apply(lambda x: ";".join(sorted(set(x))))
            .reset_index()
            .rename(columns={"sample": "barcode", "VEP_SYMBOL": "lof_metabolic_genes"})
        )
        df = df.merge(met_genes, how="left", on="barcode")
    else:
        df["lof_metabolic_genes"] = np.nan

    if len(lof_met_const) > 0:
        met_const_genes = (
            lof_met_const.groupby("sample")["VEP_SYMBOL"]
            .apply(lambda x: ";".join(sorted(set(x))))
            .reset_index()
            .rename(columns={"sample": "barcode", "VEP_SYMBOL": "lof_metabolic_constraint_genes"})
        )
        df = df.merge(met_const_genes, how="left", on="barcode")
    else:
        df["lof_metabolic_constraint_genes"] = np.nan

    # LOEUF scores for metabolic LoF gene columns
    met_gene_cols = [c for c in ["lof_metabolic_genes", "lof_metabolic_constraint_genes"]
                     if c in df.columns]
    df = assign_best_gene_scores(df, met_gene_cols, loeuf_dict)

    # 5c. Combined HCNDD dominant + X-linked recessive (males only) for LoF
    df["lof_hcndddom_xlinked_boyz_contraint_carrier"] = df.apply(
        lambda r: bool(r.get("lof_hcndddomv6_contraint_carrier", False))
        or (bool(r.get("lof_hcnddxlinked_contraint_carrier", False)) and r["Sex"] == 1),
        axis=1,
    )
    df["lof_hcndddom_xlinked_boyz_carrier"] = df.apply(
        lambda r: bool(r.get("lof_hcndddomv6_carrier", False))
        or (bool(r.get("lof_hcnddxlinked_carrier", False)) and r["Sex"] == 1),
        axis=1,
    )

    # ================================================================
    # 6. Process DELETIONS (same pipeline as hg19 script 01)
    # ================================================================
    print("\nProcessing deletions...")
    df_del = pd.read_table(DEL_FILE)
    df_del = df_del[df_del["coding"] != False]

    # Add missing de novo CNVs
    if os.path.exists(MISSING_DENOVO_CNV):
        temp = pd.read_csv(MISSING_DENOVO_CNV, sep=";")
        df_del = pd.concat([df_del, temp])
        df_del = df_del.drop_duplicates("LOCATION_SVTYPE_SAMPLE", keep="last")
    df_del = df_del.reset_index(drop=True)

    # Filter
    notes_to_exclude = [
        "flag_frequency", "leukemia", "non_coding",
        "non_coding, flag_frequency", "non_coding,flag_frequency",
    ]
    df_del = df_del[
        (~df_del["notes"].isin(notes_to_exclude))
        & (df_del["coding"] == True)
        & (df_del["outlier_DEL"] != True)
        & (df_del["wave_ind"] != True)
        & (df_del["B_loss_AFmax"] <= 0)
    ]
    df_del = df_del[
        (df_del["validation"] == True)
        | ((df_del["validation"].isna()) & (df_del["predicted_status"] == 1))
    ]
    print(f"  {len(df_del)} DEL after filtering")

    # Build DEL carrier lists (standard gene lists)
    del_carrier_constrained = {
        "hcndddom": df_del[df_del["constraint_hcndddom_v6_cds"] == True]["SAMPLE"].tolist(),
        "hcnddxlinked": df_del[df_del["constraint_hcnddxlinkedrec_v6_cds"] == True]["SAMPLE"].tolist(),
        "sparksfari1": df_del[df_del["constraint_spark_sfari1_cds"] == True]["SAMPLE"].tolist(),
        "eagle": df_del[df_del["constraint_eagle_ds_cds"] == True]["SAMPLE"].tolist(),
        "syngo": df_del[df_del["constraint_syngo_cds"] == True]["SAMPLE"].tolist(),
        "chromepitf": df_del[df_del["constraint_chromepitf_cds"] == True]["SAMPLE"].tolist(),
        "constraint": df_del[df_del["constraint_cds"] == True]["SAMPLE"].tolist(),
    }
    del_carrier_constrained["syngo_chromepitf"] = (
        del_carrier_constrained["syngo"] + del_carrier_constrained["chromepitf"]
    )

    del_carrier_all = {
        "hcndddom": df_del[df_del["hcndddom_v6_cds"] == True]["SAMPLE"].tolist(),
        "hcnddxlinked": df_del[df_del["hcnddxlinkedrec_v6_cds"] == True]["SAMPLE"].tolist(),
        "sparksfari1": df_del[df_del["spark_sfari1_cds"] == True]["SAMPLE"].tolist(),
        "eagle": df_del[df_del["eagle_ds_cds"] == True]["SAMPLE"].tolist(),
        "syngo": df_del[df_del["syngo_cds"] == True]["SAMPLE"].tolist(),
        "chromepitf": df_del[df_del["chromepitf_cds"] == True]["SAMPLE"].tolist(),
        "constraint": df_del["SAMPLE"].tolist(),
    }
    del_carrier_all["syngo_chromepitf"] = (
        del_carrier_all["syngo"] + del_carrier_all["chromepitf"]
    )

    # Add metabolic DEL carriers (check Gene_name against metabolic gene list)
    del_carrier_constrained["metabolic"] = df_del[
        df_del["Gene_name"].apply(
            lambda x: any(g in metabolic_constraint_set for g in str(x).split(";"))
        )
    ]["SAMPLE"].tolist()
    del_carrier_all["metabolic"] = df_del[
        df_del["Gene_name"].apply(
            lambda x: any(g in metabolic_set for g in str(x).split(";"))
        )
    ]["SAMPLE"].tolist()

    # Annotate standard DEL carrier flags
    df = annotate_carrier_flags(df, del_carrier_constrained, del_carrier_all, "del")

    # Add metabolic DEL carrier flags (not covered by annotate_carrier_flags)
    df["del_metabolic_contraint_carrier"] = df["barcode"].isin(del_carrier_constrained["metabolic"])
    df["del_metabolic_carrier"] = df["barcode"].isin(del_carrier_all["metabolic"])
    print(f"  del_metabolic_carrier: {df['del_metabolic_carrier'].sum()}")

    # Combined HCNDD dominant + X-linked (males only) for DEL
    df = add_xlinked_dominant_combined(
        df,
        del_carrier_constrained["hcndddom"], del_carrier_constrained["hcnddxlinked"],
        del_carrier_all["hcndddom"], del_carrier_all["hcnddxlinked"],
        "del",
    )

    # DEL gene names (standard lists)
    del_genes = extract_del_gene_names(
        df_del, gene_lists, constrained_lists, gene_lists["protein_coding"]
    )
    df = df.merge(
        del_genes.rename(columns={"del_SAMPLE": "barcode"}),
        how="left", on="barcode",
    )

    # Metabolic DEL gene names
    met_del_genes = extract_metabolic_del_genes(
        df_del, gene_lists["metabolic"], constrained_lists["metabolic"]
    )
    df = df.merge(met_del_genes, how="left", on="barcode")

    # LOEUF scores for all DEL gene columns
    del_gene_cols = [c for c in df.columns
                     if c.startswith("del_") and c.endswith("_genes")]
    df = assign_best_gene_scores(df, del_gene_cols, loeuf_dict)
    df = combine_xlinked_scores(df, "del")

    # ================================================================
    # 7. Merge DEL + LoF → dellof
    # ================================================================
    print("\nMerging DEL + LoF carrier flags...")
    df = merge_carrier_flags_ext(df, "del", "lof", "dellof")
    n_dellof = df["dellof_contraint_carrier"].sum() if "dellof_contraint_carrier" in df.columns else 0
    print(f"  dellof_contraint_carrier: {n_dellof}")

    # ================================================================
    # 8. Export
    # ================================================================
    output_path = os.path.join(TABLES_DIR_HG38, "carrier_annotations_hg38.tsv")
    df.to_csv(output_path, sep="\t", index=False)
    print(f"\nhg38 carrier annotations saved: {output_path}")
    print(f"Total individuals: {len(df)}")
    print(f"Columns: {len(df.columns)}")
    print("\nDone.")


if __name__ == "__main__":
    main()
