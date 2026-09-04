#!/usr/bin/env python3
# =============================================================================
# Carrier Annotation, hg19/gnomAD v2 (legacy -- superseded by gnomAD v4)
# =============================================================================
# Annotates individuals as carriers of rare variants (DEL, LoF, missense)
# across multiple gene lists (HCNDD, SPARK/SFARI, EAGLE, SynGO, ChromEpiTF).
# Assigns gnomAD LOEUF constraint scores to each carrier.
#
# Refactored from: script_zakaria/carrier_loeufscore.py
# Output: carrier_annotations.tsv (one row per individual, carrier flags + LOEUF scores)
# =============================================================================

import os
import sys
import pandas as pd
import numpy as np

_ORIGINAL_SCRIPTS_DIR = os.path.normpath(os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "..", "..", "..",
    "2_genetic_analysis", "scripts"))
sys.path.insert(0, _ORIGINAL_SCRIPTS_DIR)
from _config import (
    GENELIST_PATH, TABLES_DIR,
    INDIVIDUALS_METRICS, LOEUF_SCORES, DIAG_FILE, CARRIER_PREPROCESS,
    DEL_FILE, MISSING_DENOVO_CNV, SNV_FILE,
)

os.makedirs(TABLES_DIR, exist_ok=True)

# ============================================================================
# Helper functions
# ============================================================================

GENE_LIST_NAMES = [
    "hcndddom", "hcnddxlinked", "sparksfari1", "eagle",
    "syngo", "chromepitf", "syngo_chromepitf", "constraint",
]


def load_gene_lists(genelist_path):
    """Load gene lists from the GenesTrek/HGNC file and build constrained intersections."""
    genetrek = pd.read_table(genelist_path)

    dominant_inheritance = [
        "Autosomal Dominant",
        "Autosomal or pseudoautosomal dominant",
        "Autosomal or pseudoautosomal dominant - maternally imprinted",
        "Autosomal or pseudoautosomal dominant - paternally imprinted",
        "X-linked Dominant",
    ]
    xlinked_rec_inheritance = ["X-linked Recessive"]
    eagle_class_tokeep = ["Definitive", "Strong"]

    # Base gene lists
    lists = {}
    lists["hcndddom"] = genetrek[
        (genetrek["High Confidence NDD genes v6.3"] == True)
        & (genetrek["Inheritance_combined2"].isin(dominant_inheritance))
    ]["symbol"].tolist()

    lists["hcnddxlinked"] = genetrek[
        (genetrek["High Confidence NDD genes v6.3"] == True)
        & (genetrek["Inheritance_combined2"].isin(xlinked_rec_inheritance))
    ]["symbol"].tolist()

    lists["sparksfari1"] = genetrek[
        (genetrek["SFARI_1 genes (updated 2025-04-03)"] == True)
        | (genetrek["SPARK genes (updated 2025-01-07"] == True)
    ]["symbol"].tolist()

    lists["eagle"] = genetrek[
        genetrek["eagle_class"].isin(eagle_class_tokeep)
    ]["symbol"].tolist()

    lists["constraint"] = genetrek[
        genetrek["Constraint_combined"] == True
    ]["symbol"].tolist()

    lists["syngo"] = genetrek[
        genetrek["SynGO v1.2 (release 2023-12-01"] == True
    ]["symbol"].tolist()

    lists["chromepitf"] = genetrek[
        genetrek["Chrom_Epi_TF_combined"] == True
    ]["symbol"].tolist()

    lists["syngo_chromepitf"] = genetrek[
        (genetrek["Chrom_Epi_TF_combined"] == True)
        | (genetrek["SynGO v1.2 (release 2023-12-01"] == True)
    ]["symbol"].tolist()

    lists["protein_coding"] = genetrek[
        genetrek["locus_group"] == "protein-coding gene"
    ]["symbol"].tolist()

    # Constrained intersections
    constraint_set = set(lists["constraint"])
    constrained = {}
    for name in ["hcndddom", "hcnddxlinked", "sparksfari1", "eagle",
                  "syngo", "chromepitf", "syngo_chromepitf"]:
        constrained[name] = list(set(lists[name]) & constraint_set)

    return lists, constrained


def load_loeuf_scores(loeuf_path):
    """Load gnomAD LOEUF scores and return a gene→score dictionary."""
    or_tr = pd.read_table(loeuf_path, compression="gzip")
    or_tr = or_tr[or_tr["oe_lof_upper"].notna()]
    return dict(zip(or_tr["gene"], or_tr["oe_lof_upper"]))


def annotate_carrier_flags(clinic_df, carrier_lists_constrained, carrier_lists_all, prefix):
    """Add carrier boolean columns for each gene list (constrained + all).

    carrier_lists_constrained/all: dict mapping gene list name → list of sample barcodes.
    prefix: "del", "lof", or "miss".
    """
    for name in GENE_LIST_NAMES:
        if name == "constraint":
            clinic_df[f"{prefix}_contraint_carrier"] = clinic_df["barcode"].isin(
                carrier_lists_constrained.get(name, [])
            )
            clinic_df[f"{prefix}_carrier"] = clinic_df["barcode"].isin(
                carrier_lists_all.get(name, [])
            )
        elif name == "syngo_chromepitf":
            clinic_df[f"{prefix}_syngo_chromepitf_contraint_carrier"] = clinic_df["barcode"].isin(
                carrier_lists_constrained.get(name, [])
            )
            clinic_df[f"{prefix}_syngo_chromepitf_carrier"] = clinic_df["barcode"].isin(
                carrier_lists_all.get(name, [])
            )
        elif name in ("hcndddom", "hcnddxlinked"):
            suffix = "domv6" if name == "hcndddom" else "xlinked"
            clinic_df[f"{prefix}_hcndd{suffix}_contraint_carrier"] = clinic_df["barcode"].isin(
                carrier_lists_constrained.get(name, [])
            )
            clinic_df[f"{prefix}_hcndd{suffix}_carrier"] = clinic_df["barcode"].isin(
                carrier_lists_all.get(name, [])
            )
        else:
            clinic_df[f"{prefix}_{name}_contraint_carrier"] = clinic_df["barcode"].isin(
                carrier_lists_constrained.get(name, [])
            )
            clinic_df[f"{prefix}_{name}_carrier"] = clinic_df["barcode"].isin(
                carrier_lists_all.get(name, [])
            )

    return clinic_df


def add_xlinked_dominant_combined(clinic_df, carrier_list_dom_constrained,
                                   carrier_list_xlinked_constrained,
                                   carrier_list_dom_all, carrier_list_xlinked_all,
                                   prefix):
    """Combine HCNDD dominant + X-linked recessive (males only) carrier flags."""
    def _combine(row, dom_list, xlinked_list):
        if row["barcode"] in dom_list:
            return True
        elif row["barcode"] in xlinked_list and row["Sex"] == 1:
            return True
        return False

    clinic_df[f"{prefix}_hcndddom_xlinked_boyz_contraint_carrier"] = clinic_df.apply(
        lambda r: _combine(r, carrier_list_dom_constrained, carrier_list_xlinked_constrained), axis=1
    )
    clinic_df[f"{prefix}_hcndddom_xlinked_boyz_carrier"] = clinic_df.apply(
        lambda r: _combine(r, carrier_list_dom_all, carrier_list_xlinked_all), axis=1
    )
    return clinic_df


def extract_del_gene_names(df_del, gene_lists, constrained_lists, protein_coding_list):
    """Extract gene names from deletion Gene_name column for each gene list.
    Returns a merged dataframe (one row per SAMPLE) with gene name columns.
    """
    df_temp = df_del[["SAMPLE", "Gene_name"]].copy()

    all_lists = {
        "hcndddomv6_constraint_genes": constrained_lists["hcndddom"],
        "hcnddxlinkedrec_constraint_genes": constrained_lists["hcnddxlinked"],
        "sparksfari1_constraint_genes": constrained_lists["sparksfari1"],
        "eagles_ds_constraint_genes": constrained_lists["eagle"],
        "syngo_constraint_genes": constrained_lists["syngo"],
        "chromepitf_constraint_genes": constrained_lists["chromepitf"],
        "syngo_chromepitf_constraint_genes": constrained_lists["syngo_chromepitf"],
        "constraint_genes": gene_lists["constraint"],
        "hcndddomv6_genes": gene_lists["hcndddom"],
        "hcnddxlinkedrec_genes": gene_lists["hcnddxlinked"],
        "sparksfari1_genes": gene_lists["sparksfari1"],
        "eagles_ds_genes": gene_lists["eagle"],
        "syngo_genes": gene_lists["syngo"],
        "chromepitf_genes": gene_lists["chromepitf"],
        "syngo_chromepitf_genes": gene_lists["syngo_chromepitf"],
        "proteinconding_genes": protein_coding_list,
    }

    for col_name, gl in all_lists.items():
        gl_set = set(gl)
        df_temp[col_name] = df_del["Gene_name"].apply(
            lambda x: ";".join([g for g in str(x).split(";") if g in gl_set]) or np.nan
        )

    gene_columns = [c for c in df_temp.columns if c != "SAMPLE"]

    def merge_unique_genes(series):
        non_na = series.dropna()
        if non_na.empty:
            return np.nan
        genes = []
        for item in non_na:
            genes.extend(item.split(";"))
        unique_genes = sorted(set(g.strip() for g in genes if g.strip()))
        return ";".join(unique_genes)

    df_merged = df_temp.groupby("SAMPLE", as_index=False).agg(
        {col: merge_unique_genes for col in gene_columns}
    )
    df_merged.columns = ["del_" + c for c in df_merged.columns]
    return df_merged


def extract_snv_gene_names(df_snv, gene_lists, constrained_lists, protein_coding_list, sample_col="sample_id"):
    """Extract gene names from SNV symbol column for each gene list.
    Returns a merged dataframe (one row per sample) with gene name columns.
    """
    gene_col = "symbol"
    df_temp = df_snv[[sample_col, gene_col]].copy()

    all_lists = {
        "hcndddomv6_constraint_genes": constrained_lists["hcndddom"],
        "hcnddxlinkedrec_constraint_genes": constrained_lists["hcnddxlinked"],
        "sparksfari1_constraint_genes": constrained_lists["sparksfari1"],
        "eagles_ds_constraint_genes": constrained_lists["eagle"],
        "syngo_constraint_genes": constrained_lists["syngo"],
        "chromepitf_constraint_genes": constrained_lists["chromepitf"],
        "syngo_chromepitf_constraint_genes": constrained_lists["syngo_chromepitf"],
        "constraint_genes": gene_lists["constraint"],
        "hcndddomv6_genes": gene_lists["hcndddom"],
        "hcnddxlinkedrec_genes": gene_lists["hcnddxlinked"],
        "sparksfari1_genes": gene_lists["sparksfari1"],
        "eagles_ds_genes": gene_lists["eagle"],
        "syngo_genes": gene_lists["syngo"],
        "chromepitf_genes": gene_lists["chromepitf"],
        "syngo_chromepitf_genes": gene_lists["syngo_chromepitf"],
        "proteinconding_genes": protein_coding_list,
    }

    for col_name, gl in all_lists.items():
        gl_set = set(gl)
        df_temp[col_name] = df_snv[gene_col].apply(
            lambda x: x if x in gl_set else np.nan
        )

    gene_columns = [c for c in df_temp.columns if c != sample_col]

    def merge_unique_genes(series):
        non_na = series.dropna()
        if non_na.empty:
            return np.nan
        unique_genes = sorted(set(g.strip() for g in non_na if g.strip()))
        return ";".join(unique_genes)

    df_merged = df_temp.groupby(sample_col, as_index=False).agg(
        {col: merge_unique_genes for col in gene_columns}
    )
    return df_merged


def assign_best_gene_scores(df, columns, gene_scores, suffix="_best"):
    """For each gene-name column, find the gene with the lowest LOEUF score."""
    def _best(gene_str, scores_dict):
        if pd.isna(gene_str):
            return pd.Series([np.nan, np.nan])
        genes = str(gene_str).split(";")
        best_gene, best_score = np.nan, np.inf
        for g in genes:
            score = scores_dict.get(g, np.nan)
            if not pd.isna(score) and score < best_score:
                best_score = score
                best_gene = g
        if best_score == np.inf:
            return pd.Series([np.nan, np.nan])
        return pd.Series([best_gene, best_score])

    for col in columns:
        if col not in df.columns:
            continue
        new_cols = df[col].apply(lambda x: _best(x, gene_scores))
        df[col + suffix + "_gene"] = new_cols[0]
        df[col + suffix + "_score"] = new_cols[1]
    return df


def combine_xlinked_scores(df, prefix):
    """Combine dominant + X-linked recessive LOEUF scores (males only) for a variant type."""
    for constraint_suffix in ["_constraint_genes", "_genes"]:
        dom_score = f"{prefix}_hcndddomv6{constraint_suffix}_best_score"
        dom_gene = f"{prefix}_hcndddomv6{constraint_suffix}_best_gene"
        xlink_score = f"{prefix}_hcnddxlinkedrec{constraint_suffix}_best_score"
        xlink_gene = f"{prefix}_hcnddxlinkedrec{constraint_suffix}_best_gene"
        out_score = f"{prefix}_hcndddomv6_xlinked_boyz{constraint_suffix}_best_score"
        out_gene = f"{prefix}_hcndddomv6_xlinked_boyz{constraint_suffix}_best_gene"

        if not all(c in df.columns for c in [dom_score, xlink_score]):
            continue

        # X-linked score only applies to males (Sex==1)
        temp_xlinked = df.apply(
            lambda r: r[xlink_score] if r["Sex"] == 1 else -np.inf, axis=1
        )
        df[out_score] = df.apply(
            lambda r: max(r[dom_score] if not pd.isna(r[dom_score]) else -np.inf,
                          temp_xlinked[r.name]),
            axis=1
        )
        df[out_score] = df[out_score].replace(-np.inf, np.nan)

        df[out_gene] = df.apply(
            lambda r: r[dom_gene]
            if (r[dom_score] if not pd.isna(r[dom_score]) else -np.inf) >= temp_xlinked[r.name]
            else r[xlink_gene],
            axis=1
        )
    return df


def combine_variant_scores(df, prefix1, prefix2, out_prefix, mode="max"):
    """Combine LOEUF scores from two variant types.
    mode='max': keep the higher score (for dellof combination)
    mode='min': keep the lower score (for min_dellof/min_dellofmiss)
    """
    base_suffixes = [
        "hcndddomv6_constraint_genes", "hcnddxlinkedrec_constraint_genes",
        "hcndddomv6_xlinked_boyz_constraint_genes",
        "sparksfari1_constraint_genes", "eagles_ds_constraint_genes",
        "syngo_constraint_genes", "chromepitf_constraint_genes",
        "syngo_chromepitf_constraint_genes", "constraint_genes",
        "hcndddomv6_genes", "hcnddxlinkedrec_genes",
        "hcndddomv6_xlinked_boyz_genes",
        "sparksfari1_genes", "eagles_ds_genes",
        "syngo_genes", "chromepitf_genes",
        "syngo_chromepitf_genes", "proteinconding_genes",
    ]

    for suffix in base_suffixes:
        s1 = f"{prefix1}_{suffix}_best_score"
        g1 = f"{prefix1}_{suffix}_best_gene"
        s2 = f"{prefix2}_{suffix}_best_score"
        g2 = f"{prefix2}_{suffix}_best_gene"
        out_s = f"{out_prefix}_{suffix}_best_score"
        out_g = f"{out_prefix}_{suffix}_best_gene"

        if not all(c in df.columns for c in [s1, g1, s2, g2]):
            continue

        if mode == "max":
            cond = df[s1].fillna(-np.inf) >= df[s2].fillna(-np.inf)
        else:
            cond = df[s1].fillna(np.inf) <= df[s2].fillna(np.inf)

        df.loc[cond, out_s] = df[s1]
        df.loc[~cond, out_s] = df[s2]
        df.loc[cond, out_g] = df[g1]
        df.loc[~cond, out_g] = df[g2]
        df.loc[df[s1].isna() & df[s2].isna(), [out_s, out_g]] = np.nan

    return df


def merge_carrier_flags(df, prefix1, prefix2, out_prefix):
    """Merge carrier boolean flags from two variant types (OR logic)."""
    suffixes = [
        "_carrier", "_hcndddomv6_contraint_carrier", "_hcnddxlinked_contraint_carrier",
        "_hcndddom_xlinked_boyz_contraint_carrier",
        "_sparksfari1_contraint_carrier", "_eagle_contraint_carrier",
        "_syngo_contraint_carrier", "_chromepitf_contraint_carrier",
        "_syngo_chromepitf_contraint_carrier", "_contraint_carrier",
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


# ============================================================================
# Main pipeline
# ============================================================================

def main():
    # --- Check for pre-processed carrier file ---
    # If raw genetic files (SNV/CNV) are not available, use an existing
    # pre-processed carrier annotation file as input.
    if not os.path.exists(SNV_FILE) and os.path.exists(CARRIER_PREPROCESS):
        print(f"Raw genetic files not found. Using pre-processed carrier file:")
        print(f"  {CARRIER_PREPROCESS}")
        clinic_filt = pd.read_table(CARRIER_PREPROCESS, low_memory=False)
        output_path = os.path.join(TABLES_DIR, "carrier_annotations.tsv")
        clinic_filt.to_csv(output_path, sep="\t", index=False)
        print(f"Carrier annotations saved: {output_path}")
        print(f"Total individuals: {len(clinic_filt)}, Columns: {len(clinic_filt.columns)}")
        print("Done (pre-processed mode).")
        return

    # --- 1. Load clinical data ---
    print(f"Loading clinical data from: {INDIVIDUALS_METRICS}")
    clinic = pd.read_table(INDIVIDUALS_METRICS)
    clinic_filt = clinic[(clinic["wgs_score"] == True) & (clinic["gene_first"] == False)]
    clinic_filt = clinic_filt[
        ["ID", "Sex", "FID", "barcode", "cohort", "EUR_ancestry",
         "PopulationS1", "Population1", "Population_undiagnosed1"]
    ].copy()

    # --- 2. Load LOEUF scores ---
    print(f"Loading LOEUF scores from: {LOEUF_SCORES}")
    loeuf_dict = load_loeuf_scores(LOEUF_SCORES)

    # --- 3. Load gene lists ---
    print(f"Loading gene lists from: {GENELIST_PATH}")
    gene_lists, constrained_lists = load_gene_lists(GENELIST_PATH)

    # =====================================================================
    # 4. Process DELETIONS
    # =====================================================================
    print("Processing deletions...")
    df_del = pd.read_table(DEL_FILE)
    df_del = df_del[df_del["coding"] != False]

    # Add missing de novo CNVs
    if os.path.exists(MISSING_DENOVO_CNV):
        temp = pd.read_csv(MISSING_DENOVO_CNV, sep=";")
        df_del = pd.concat([df_del, temp])
        df_del = df_del.drop_duplicates("LOCATION_SVTYPE_SAMPLE", keep="last")
    df_del = df_del.reset_index(drop=True)

    # Filter: remove flag_frequency, non_coding, outliers, WAVE, keep absent from gnomAD
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

    # Keep validated or ML-predicted deletions
    df_del = df_del[
        (df_del["validation"] == True)
        | ((df_del["validation"].isna()) & (df_del["predicted_status"] == 1))
    ]

    # Build carrier lists
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

    # Annotate carrier flags
    clinic_filt = annotate_carrier_flags(clinic_filt, del_carrier_constrained, del_carrier_all, "del")
    clinic_filt = add_xlinked_dominant_combined(
        clinic_filt,
        del_carrier_constrained["hcndddom"], del_carrier_constrained["hcnddxlinked"],
        del_carrier_all["hcndddom"], del_carrier_all["hcnddxlinked"],
        "del",
    )

    # Extract gene names and assign LOEUF scores
    del_genes = extract_del_gene_names(df_del, gene_lists, constrained_lists, gene_lists["protein_coding"])
    clinic_filt = clinic_filt.merge(del_genes.rename(columns={"del_SAMPLE": "barcode"}), how="left", on="barcode")
    del_gene_cols = [c for c in clinic_filt.columns if c.startswith("del_") and c.endswith("_genes")]
    clinic_filt = assign_best_gene_scores(clinic_filt, del_gene_cols, loeuf_dict)
    clinic_filt = combine_xlinked_scores(clinic_filt, "del")

    # =====================================================================
    # 5. Process LoF SNVs
    # =====================================================================
    print("Processing LoF SNVs...")
    lof_df = pd.read_excel(SNV_FILE)
    lof_df_filt = lof_df[
        (lof_df["comment"] != "to_exclude")
        & (lof_df["visualized"] != "dup")
        & (lof_df["visualized"] != "mnv")
        & (lof_df["visualized"] != "to_check")
        & (lof_df["isLCR"] != True)
        & (lof_df["impact"] != "inframe_variant")
    ]
    lof_df_filt = lof_df_filt[
        ((lof_df_filt["visualized"].isna()) & (lof_df_filt["freddyscore"] == True))
        | (lof_df_filt["validation"] == True)
    ]
    lof_df_filt = lof_df_filt[
        (lof_df_filt["canonical"] == "YES")
        & (lof_df_filt["loftee_flags"].isna())
        & (lof_df_filt["loftee"] == "HC")
    ]

    sample = "sample_id"
    lof_carrier_constrained = {
        "hcndddom": lof_df_filt[(lof_df_filt["constraint"] == True) & (lof_df_filt["hcndddom_v6"] == True)][sample].tolist(),
        "hcnddxlinked": lof_df_filt[(lof_df_filt["constraint"] == True) & (lof_df_filt["hcnddxlinkedrec_v6"] == True)][sample].tolist(),
        "sparksfari1": lof_df_filt[(lof_df_filt["constraint"] == True) & (lof_df_filt["spark_sfari1"] == True)][sample].tolist(),
        "eagle": lof_df_filt[(lof_df_filt["constraint"] == True) & (lof_df_filt["eagle_ds"] == True)][sample].tolist(),
        "syngo": lof_df_filt[(lof_df_filt["constraint"] == True) & (lof_df_filt["syngo"] == True)][sample].tolist(),
        "chromepitf": lof_df_filt[(lof_df_filt["constraint"] == True) & (lof_df_filt["chromepitf"] == True)][sample].tolist(),
        "constraint": lof_df_filt[lof_df_filt["constraint"] == True][sample].tolist(),
    }
    lof_carrier_constrained["syngo_chromepitf"] = (
        lof_carrier_constrained["syngo"] + lof_carrier_constrained["chromepitf"]
    )

    lof_carrier_all = {
        "hcndddom": lof_df_filt[lof_df_filt["hcndddom_v6"] == True][sample].tolist(),
        "hcnddxlinked": lof_df_filt[lof_df_filt["hcnddxlinkedrec_v6"] == True][sample].tolist(),
        "sparksfari1": lof_df_filt[lof_df_filt["spark_sfari1"] == True][sample].tolist(),
        "eagle": lof_df_filt[lof_df_filt["eagle_ds"] == True][sample].tolist(),
        "syngo": lof_df_filt[lof_df_filt["syngo"] == True][sample].tolist(),
        "chromepitf": lof_df_filt[lof_df_filt["chromepitf"] == True][sample].tolist(),
        "constraint": lof_df_filt[sample].tolist(),
    }
    lof_carrier_all["syngo_chromepitf"] = (
        lof_carrier_all["syngo"] + lof_carrier_all["chromepitf"]
    )

    clinic_filt = annotate_carrier_flags(clinic_filt, lof_carrier_constrained, lof_carrier_all, "lof")
    clinic_filt = add_xlinked_dominant_combined(
        clinic_filt,
        lof_carrier_constrained["hcndddom"], lof_carrier_constrained["hcnddxlinked"],
        lof_carrier_all["hcndddom"], lof_carrier_all["hcnddxlinked"],
        "lof",
    )

    lof_genes = extract_snv_gene_names(lof_df_filt, gene_lists, constrained_lists, gene_lists["protein_coding"])
    lof_genes.columns = ["lof_" + c for c in lof_genes.columns]
    clinic_filt = clinic_filt.merge(
        lof_genes.rename(columns={"lof_sample_id": "barcode"}), how="left", on="barcode"
    )
    lof_gene_cols = [c for c in clinic_filt.columns if c.startswith("lof_") and c.endswith("_genes")]
    clinic_filt = assign_best_gene_scores(clinic_filt, lof_gene_cols, loeuf_dict)
    clinic_filt = combine_xlinked_scores(clinic_filt, "lof")

    # =====================================================================
    # 6. Process Missenses (AlphaMissense pathogenic)
    # =====================================================================
    print("Processing missense variants...")
    miss_df = pd.read_excel(SNV_FILE)
    miss_df_filt = miss_df[
        (miss_df["comment"] != "to_exclude")
        & (miss_df["visualized"] != "dup")
        & (miss_df["visualized"] != "mnv")
        & (miss_df["visualized"] != "to_check")
        & (miss_df["isLCR"] != True)
        & (miss_df["impact"] != "inframe_variant")
    ]
    miss_df_filt = miss_df_filt[
        (miss_df_filt["impact"] == "missense_variant")
        & (miss_df_filt["am_class"] == "pathogenic")
        & (miss_df["comment"] != "benign")
    ]
    miss_df_filt = miss_df_filt[
        ((miss_df_filt["visualized"].isna()) & (miss_df_filt["freddyscore"] == True))
        | (miss_df_filt["validation"] == True)
    ]

    miss_carrier_constrained = {
        "hcndddom": miss_df_filt[(miss_df_filt["constraint"] == True) & (miss_df_filt["hcndddom_v6"] == True)][sample].tolist(),
        "hcnddxlinked": miss_df_filt[(miss_df_filt["constraint"] == True) & (miss_df_filt["hcnddxlinkedrec_v6"] == True)][sample].tolist(),
        "sparksfari1": miss_df_filt[(miss_df_filt["constraint"] == True) & (miss_df_filt["spark_sfari1"] == True)][sample].tolist(),
        "eagle": miss_df_filt[(miss_df_filt["constraint"] == True) & (miss_df_filt["eagle_ds"] == True)][sample].tolist(),
        "syngo": miss_df_filt[(miss_df_filt["constraint"] == True) & (miss_df_filt["syngo"] == True)][sample].tolist(),
        "chromepitf": miss_df_filt[(miss_df_filt["constraint"] == True) & (miss_df_filt["chromepitf"] == True)][sample].tolist(),
        "constraint": miss_df_filt[miss_df_filt["constraint"] == True][sample].tolist(),
    }
    miss_carrier_constrained["syngo_chromepitf"] = (
        miss_carrier_constrained["syngo"] + miss_carrier_constrained["chromepitf"]
    )

    miss_carrier_all = {
        "hcndddom": miss_df_filt[miss_df_filt["hcndddom_v6"] == True][sample].tolist(),
        "hcnddxlinked": miss_df_filt[miss_df_filt["hcnddxlinkedrec_v6"] == True][sample].tolist(),
        "sparksfari1": miss_df_filt[miss_df_filt["spark_sfari1"] == True][sample].tolist(),
        "eagle": miss_df_filt[miss_df_filt["eagle_ds"] == True][sample].tolist(),
        "syngo": miss_df_filt[miss_df_filt["syngo"] == True][sample].tolist(),
        "chromepitf": miss_df_filt[miss_df_filt["chromepitf"] == True][sample].tolist(),
        "constraint": miss_df_filt[sample].tolist(),
    }
    miss_carrier_all["syngo_chromepitf"] = (
        miss_carrier_all["syngo"] + miss_carrier_all["chromepitf"]
    )

    clinic_filt = annotate_carrier_flags(clinic_filt, miss_carrier_constrained, miss_carrier_all, "miss")
    clinic_filt = add_xlinked_dominant_combined(
        clinic_filt,
        miss_carrier_constrained["hcndddom"], miss_carrier_constrained["hcnddxlinked"],
        miss_carrier_all["hcndddom"], miss_carrier_all["hcnddxlinked"],
        "miss",
    )

    miss_genes = extract_snv_gene_names(miss_df_filt, gene_lists, constrained_lists, gene_lists["protein_coding"])
    miss_genes.columns = ["miss_" + c for c in miss_genes.columns]
    clinic_filt = clinic_filt.merge(
        miss_genes.rename(columns={"miss_sample_id": "barcode"}), how="left", on="barcode"
    )
    miss_gene_cols = [c for c in clinic_filt.columns if c.startswith("miss_") and c.endswith("_genes")]
    clinic_filt = assign_best_gene_scores(clinic_filt, miss_gene_cols, loeuf_dict)
    clinic_filt = combine_xlinked_scores(clinic_filt, "miss")

    # =====================================================================
    # 7. Merge DEL + LoF → dellof
    # =====================================================================
    print("Merging DEL + LoF carrier flags...")
    clinic_filt = merge_carrier_flags(clinic_filt, "del", "lof", "dellof")
    clinic_filt = combine_variant_scores(clinic_filt, "lof", "del", "dellof", mode="max")

    # Also compute min scores for dellof
    clinic_filt = combine_variant_scores(clinic_filt, "lof", "del", "min_dellof", mode="min")

    # =====================================================================
    # 8. Merge DEL + LoF + Miss → dellofmiss
    # =====================================================================
    print("Merging DEL + LoF + Miss carrier flags...")
    clinic_filt = merge_carrier_flags(clinic_filt, "dellof", "miss", "dellofmiss")
    clinic_filt = combine_variant_scores(clinic_filt, "miss", "dellof", "min_dellofmiss", mode="min")

    # =====================================================================
    # 9. Export
    # =====================================================================
    output_path = os.path.join(TABLES_DIR, "carrier_annotations.tsv")
    clinic_filt.to_csv(output_path, sep="\t", index=False)
    print(f"\nCarrier annotations saved: {output_path}")
    print(f"Total individuals: {len(clinic_filt)}")
    print(f"Columns: {len(clinic_filt.columns)}")
    print("\nDone.")


if __name__ == "__main__":
    main()
