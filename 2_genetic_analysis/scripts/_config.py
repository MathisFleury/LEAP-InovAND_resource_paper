# =============================================================================
# Shared Configuration for Genetic Analysis Scripts
# =============================================================================
# Paths, palettes, label mappings used across all scripts in 2_genetic_analysis.
# =============================================================================

import os
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

# --- Paths ---
_script_dir = os.path.dirname(os.path.abspath(__file__))
_lib_dir = os.path.normpath(os.path.join(_script_dir, "..", "..", ".."))

OUTPUT_BASE = os.path.normpath(os.path.join(_script_dir, "..", "outputs"))
FIGURES_DIR = os.path.join(OUTPUT_BASE, "figures")
TABLES_DIR = os.path.join(OUTPUT_BASE, "tables")

# --- Input data files (hardcoded) ---
INDIVIDUALS_METRICS = os.path.join(
    _lib_dir, "imaging2genet", "0_input", "dataframes", "individuals_metrics.tsv"
)
GENELIST_PATH = os.path.join(
    _lib_dir, "eeg_mri-pipeline", "ressources", "genetics",
    "hgnc_complete_set_20250424_GeneList_Updated20260225.txt"
)
LOEUF_SCORES = os.path.join(
    _lib_dir, "LEAP-InovAND_resource", "script_zakaria",
    "gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz"
)
DEL_FILE = os.path.join(
    _lib_dir, "LEAP-InovAND_resource", "script_zakaria",
    "SV_DEL_annotation_12juin2025_ZM.tsv"
)
MISSING_DENOVO_CNV = os.path.join(
    _lib_dir, "LEAP-InovAND_resource", "script_zakaria", "missing_denovo_cnv.csv"
)
SNV_FILE = os.path.join(
    _lib_dir, "LEAP-InovAND_resource", "script_zakaria",
    "slivar_all_ensg_annotated1 3.xlsx"
)
DIAG_FILE = os.path.join(
    "/Volumes", "Imaging5", "EEG_MRI-MF", "ALL", "results", "tabular", "genetics",
    "diagnostic_clinic_JM", "diag_individuals_combined_260408.csv"
)
CARRIER_PREPROCESS = os.path.join(
    _lib_dir, "eeg_mri-pipeline", "results", "dataset_paper", "dataframes",
    "df_carrier_genelist_DEL_LOF_MISS_withalphamissense_sex_takenintoaccount.tsv"
)
CLUSTERS_FILE = os.path.join(
    _lib_dir, "eeg_mri-pipeline", "results", "dataset_paper", "dataframes",
    "df_clusters_complete_kmeans.csv"
)
PGS_FILE = os.path.join(
    _lib_dir, "eeg_mri-pipeline", "results", "dataset_paper", "dataframes",
    "InovAND-LEAP_pgs_SBayesRC_20250401.tsv"
)
CARRIER_HG38 = os.path.join(
    "/Volumes", "Imaging5", "EEG_MRI-MF", "ALL", "results", "tabular", "genetics",
    "df_carrier_genelist_DEL_LOF_MISS_withalphamissense_LOEUF_gnomadv2_grch38.tsv"
)
DUP_FILE = os.path.join(
    "/Volumes", "Imaging5", "EEG_MRI-MF", "ALL", "results", "tabular", "genetics",
    "dup_per_individual_paperset.tsv"
)
LOF_HG38_FILE = os.path.join(
    "/Volumes", "Imaging5", "EEG_MRI-MF", "ALL", "results", "tabular", "genetics",
    "df_38_with_lof_carrier_status.tsv"
)
VARIANTS_HG38_FILE = os.path.join(
    "/Volumes", "Imaging5", "EEG_MRI-MF", "ALL", "results", "tabular", "genetics",
    "df_variants_GRCh38.tsv"
)
FIGURES_DIR_HG38 = os.path.join(OUTPUT_BASE, "figures_hg38")
TABLES_DIR_HG38 = os.path.join(OUTPUT_BASE, "tables_hg38")

# --- Matplotlib defaults ---
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Helvetica", "Arial", "DejaVu Sans", "Liberation Sans"]
mpl.rcParams["pdf.fonttype"] = 42
sns.set_context("paper")

# --- Palettes and group orders ---
PALETTE_POPULATION1 = {
    "Relatives": "#9AD5D3",
    "IDD": "#D8A4CB",
    "NT": "#C1C2BC",
    "Autism with IDD": "#324095",
    "Autism without IDD": "#5CAEE1",
}
ORDER_POPULATION1 = ["Autism without IDD", "Autism with IDD", "IDD", "Relatives", "NT"]

PALETTE_POPULATIONS = {
    "Autism": "#8991FA",
    "Relatives": "#9AD5D3",
    "IDD": "#D8A4CB",
    "NT": "#C1C2BC",
}
ORDER_POPULATIONS = ["Autism", "IDD", "Relatives", "NT"]

PALETTE_FREQ = {
    "Autism": "#8991FA",
    "Undiagnosed siblings": "#B0D5D3",
    "Undiagnosed parents": "#9AD5D3",
    "Undiagnosed 2+": "#9AD5D3",
    "Relatives": "#9AD5D3",
    "IDD": "#D8A4CB",
    "NT": "#C1C2BC",
    "Autism with IDD": "#324095",
    "Autism without IDD": "#5CAEE1",
}
ORDER_FREQ = ["Autism", "IDD", "Undiagnosed siblings", "Undiagnosed parents", "NT"]
ORDER_OR = ["Autism", "IDD", "Undiagnosed siblings", "Undiagnosed parents"]

PALETTE_CLUSTERS = {
    "NT": "#C1C2BC",
    "C1": "#7A8B47",
    "C2": "#ff9fa0",
    "C3": "#e7ba52",
    "IDD": "#D8A4CB",
}
ORDER_CLUSTERS = ["NT", "C1", "C2", "C3", "IDD"]

PALETTE_ANCESTRY = {
    "AFR": "#DB5E57",
    "AMR": "#DBD057",
    "CSA": "#74DB57",
    "EAS": "#58DBAA",
    "EUR": "#579BDB",
    "MID": "#8557DB",
    "UNK": "#A9B3C7",
}
ORDER_ANCESTRY = ["AFR", "AMR", "CSA", "EAS", "EUR", "MID", "UNK"]

# --- Gene list label mappings ---
LABELS_CONSTRAINED = {
    "diag_genetic": "Returnable causative\nvariants",
    "dellof_contraint_carrier": "Constrained",
    "dellof_hcndddom_xlinked_boyz_contraint_carrier": "Constrained\nHCNDD dominant and\nrec X-linked",
    "dellof_sparksfari1_contraint_carrier": "Constrained\nSPARK SFARI1",
    "dellof_eagle_contraint_carrier": "Constrained\nEAGLE d&s",
    "dellof_syngo_contraint_carrier": "Constrained\nSynGO",
    "dellof_chromepitf_contraint_carrier": "Constrained\nChromEpiTF",
}

LABELS_ALL_GENES = {
    "dellof_contraint_carrier": "Constraint",
    "dellof_hcndddom_xlinked_boyz_carrier": "HCNDD dominant and\nrec X-linked",
    "dellof_sparksfari1_carrier": "SPARK_SFARI1",
    "dellof_eagle_carrier": "EAGLE d&s",
    "dellof_syngo_carrier": "SynGO",
    "dellof_chromepitf_carrier": "ChromEpiTF",
}

LABELS_DELLOFMISS_CONSTRAINED = {
    "diag_genetic": "Returnable causative\nvariants",
    "dellofmiss_contraint_carrier": "Constrained",
    "dellofmiss_hcndddom_xlinked_boyz_contraint_carrier": "Constrained\nHCNDD dominant and\nrec X-linked",
    "dellofmiss_sparksfari1_contraint_carrier": "Constrained\nSPARK SFARI1",
    "dellofmiss_eagle_contraint_carrier": "Constrained\nEAGLE d&s",
    "dellofmiss_syngo_contraint_carrier": "Constrained\nSynGO",
    "dellofmiss_chromepitf_contraint_carrier": "Constrained\nChromEpiTF",
}

LABELS_DELLOFMISS_ALL = {
    "dellofmiss_contraint_carrier": "Constraint",
    "dellofmiss_hcndddom_xlinked_boyz_carrier": "HCNDD dominant and\nrec X-linked",
    "dellofmiss_sparksfari1_carrier": "SPARK_SFARI1",
    "dellofmiss_eagle_carrier": "EAGLE d&s",
    "dellofmiss_syngo_carrier": "SynGO",
    "dellofmiss_chromepitf_carrier": "ChromEpiTF",
}

# --- Duplications (overall "Constrained" + 5 gene lists matching deletions) ---
# "dup_any_constraint_carrier" is a synthetic column built at runtime by
# OR-ing the 5 per-list dup constraint flags (the dup file has no overall flag).
# TODO: when DUP_FILE gains an overall "any constrained gene" carrier flag
# (e.g. dup_constraint_carrier, analogue of dellof_contraint_carrier), drop the
# synthetic name here and the OR computation in scripts 02/04, and update
# _DUP_ANY_CONSTRAINED below to reference that single column instead of the 5.
COLS_DUP_CONSTRAINED = [
    "dup_any_constraint_carrier",
    "dup_hcnddv7_constraint_carrier",
    "dup_sfari1_constraint_carrier",
    "dup_eagle_constraint_carrier",
    "dup_syngo_constraint_carrier",
    "dup_chromepitf_constraint_carrier",
]
LABELS_DUP_CONSTRAINED = {
    "dup_any_constraint_carrier": "Constrained",
    "dup_hcnddv7_constraint_carrier": "Constrained\nHCNDD dominant and\nrec X-linked",
    "dup_sfari1_constraint_carrier": "Constrained\nSPARK SFARI1",
    "dup_eagle_constraint_carrier": "Constrained\nEAGLE d&s",
    "dup_syngo_constraint_carrier": "Constrained\nSynGO",
    "dup_chromepitf_constraint_carrier": "Constrained\nChromEpiTF",
}

# --- Combined-variant panels (1 overall "Constrained" + 5 gene lists) ---
COMBO_PANEL_KEYS = ["constrained", "hcndd", "sparksfari1", "eagle", "syngo", "chromepitf"]
COMBO_PANEL_LABELS = {
    "constrained": "Constrained",
    "hcndd": "Constrained\nHCNDD dominant and\nrec X-linked",
    "sparksfari1": "Constrained\nSPARK SFARI1",
    "eagle": "Constrained\nEAGLE d&s",
    "syngo": "Constrained\nSynGO",
    "chromepitf": "Constrained\nChromEpiTF",
}
# Dup file has no overall "any constrained" flag; proxy by OR-ing the 5 per-list flags.
# TODO: replace this list with the single real column when DUP_FILE provides it.
_DUP_ANY_CONSTRAINED = [
    "dup_hcnddv7_constraint_carrier",
    "dup_sfari1_constraint_carrier",
    "dup_eagle_constraint_carrier",
    "dup_syngo_constraint_carrier",
    "dup_chromepitf_constraint_carrier",
]
COMBO_SOURCES_DELLOF = {
    "constrained": ["dellof_contraint_carrier"],
    "hcndd": ["dellof_hcndddom_xlinked_boyz_contraint_carrier"],
    "sparksfari1": ["dellof_sparksfari1_contraint_carrier"],
    "eagle": ["dellof_eagle_contraint_carrier"],
    "syngo": ["dellof_syngo_contraint_carrier"],
    "chromepitf": ["dellof_chromepitf_contraint_carrier"],
}
COMBO_SOURCES_DELLOFDUP = {
    "constrained": ["dellof_contraint_carrier"] + _DUP_ANY_CONSTRAINED,
    "hcndd": ["dellof_hcndddom_xlinked_boyz_contraint_carrier", "dup_hcnddv7_constraint_carrier"],
    "sparksfari1": ["dellof_sparksfari1_contraint_carrier", "dup_sfari1_constraint_carrier"],
    "eagle": ["dellof_eagle_contraint_carrier", "dup_eagle_constraint_carrier"],
    "syngo": ["dellof_syngo_contraint_carrier", "dup_syngo_constraint_carrier"],
    "chromepitf": ["dellof_chromepitf_contraint_carrier", "dup_chromepitf_constraint_carrier"],
}
COMBO_SOURCES_DELLOFDUPMISS = {
    "constrained": ["dellofmiss_contraint_carrier"] + _DUP_ANY_CONSTRAINED,
    "hcndd": ["dellofmiss_hcndddom_xlinked_boyz_contraint_carrier", "dup_hcnddv7_constraint_carrier"],
    "sparksfari1": ["dellofmiss_sparksfari1_contraint_carrier", "dup_sfari1_constraint_carrier"],
    "eagle": ["dellofmiss_eagle_contraint_carrier", "dup_eagle_constraint_carrier"],
    "syngo": ["dellofmiss_syngo_contraint_carrier", "dup_syngo_constraint_carrier"],
    "chromepitf": ["dellofmiss_chromepitf_contraint_carrier", "dup_chromepitf_constraint_carrier"],
}

# hg38 variant: the generated hg38 carrier file lacks dellofmiss_* columns,
# so build the DEL+LoF+Miss flag on the fly by OR-ing dellof_* + miss_*.
COMBO_SOURCES_DELLOFDUPMISS_HG38 = {
    "constrained": [
        "dellof_contraint_carrier",
        "miss_contraint_carrier",
    ] + _DUP_ANY_CONSTRAINED,
    "hcndd": [
        "dellof_hcndddom_xlinked_boyz_contraint_carrier",
        "miss_hcndddom_xlinked_boyz_contraint_carrier",
        "dup_hcnddv7_constraint_carrier",
    ],
    "sparksfari1": [
        "dellof_sparksfari1_contraint_carrier",
        "miss_sparksfari1_contraint_carrier",
        "dup_sfari1_constraint_carrier",
    ],
    "eagle": [
        "dellof_eagle_contraint_carrier",
        "miss_eagle_contraint_carrier",
        "dup_eagle_constraint_carrier",
    ],
    "syngo": [
        "dellof_syngo_contraint_carrier",
        "miss_syngo_contraint_carrier",
        "dup_syngo_constraint_carrier",
    ],
    "chromepitf": [
        "dellof_chromepitf_contraint_carrier",
        "miss_chromepitf_contraint_carrier",
        "dup_chromepitf_constraint_carrier",
    ],
}

LABELS_PGS = {
    "autism_grove2019": "Autism",
    "adhd_demontis2023": "ADHD",
    "int_savage2018": "Intelligence",
    "anxiety_purves2020": "Anxiety",
    "mdd_meng2024": "MDD",
    "ptsd_nievergelt2024": "PTSD",
    "ukbiobank_watanabe2019/felt_loved": "Felt loved as a child",
    "age_diagnosis_zhang2025/early": "Early autism diagnosis",
    "age_diagnosis_zhang2025/late": "Late autism diagnosis",
}

# --- Constrained gene list columns (DEL+LoF) ---
COLS_DELLOF_CONSTRAINED = [
    "diag_genetic",
    "dellof_contraint_carrier",
    "dellof_hcndddom_xlinked_boyz_contraint_carrier",
    "dellof_sparksfari1_contraint_carrier",
    "dellof_eagle_contraint_carrier",
    "dellof_syngo_contraint_carrier",
    "dellof_chromepitf_contraint_carrier",
]

COLS_DELLOF_ALL = [
    "dellof_contraint_carrier",
    "dellof_hcndddom_xlinked_boyz_carrier",
    "dellof_sparksfari1_carrier",
    "dellof_eagle_carrier",
    "dellof_syngo_carrier",
    "dellof_chromepitf_carrier",
]

COLS_DELLOFMISS_CONSTRAINED = [
    "diag_genetic",
    "dellofmiss_contraint_carrier",
    "dellofmiss_hcndddom_xlinked_boyz_contraint_carrier",
    "dellofmiss_sparksfari1_contraint_carrier",
    "dellofmiss_eagle_contraint_carrier",
    "dellofmiss_syngo_contraint_carrier",
    "dellofmiss_chromepitf_contraint_carrier",
]

COLS_DELLOFMISS_ALL = [
    "dellofmiss_contraint_carrier",
    "dellofmiss_hcndddom_xlinked_boyz_carrier",
    "dellofmiss_sparksfari1_carrier",
    "dellofmiss_eagle_carrier",
    "dellofmiss_syngo_carrier",
    "dellofmiss_chromepitf_carrier",
]

PGS_TRAITS = [
    "autism_grove2019",
    "age_diagnosis_zhang2025/early",
    "ptsd_nievergelt2024",
    "anxiety_purves2020",
    "mdd_meng2024",
    "adhd_demontis2023",
    "age_diagnosis_zhang2025/late",
    "int_savage2018",
    "ukbiobank_watanabe2019/felt_loved",
]

# --- Extended column lists (hg38: + metabolic + syngo_chromepitf) ---
COLS_DELLOF_CONSTRAINED_EXT = COLS_DELLOF_CONSTRAINED + [
    "dellof_metabolic_contraint_carrier",
]
COLS_DELLOF_ALL_EXT = COLS_DELLOF_ALL + [
    "dellof_metabolic_carrier",
]
COLS_DELLOFMISS_CONSTRAINED_EXT = COLS_DELLOFMISS_CONSTRAINED + [
    "dellofmiss_metabolic_contraint_carrier",
]
COLS_DELLOFMISS_ALL_EXT = COLS_DELLOFMISS_ALL + [
    "dellofmiss_metabolic_carrier",
]

# --- Miss-only column lists and labels (for hg38 where dellofmiss unavailable) ---
COLS_MISS_CONSTRAINED = [
    "diag_genetic",
    "miss_contraint_carrier",
    "miss_hcndddom_xlinked_boyz_contraint_carrier",
    "miss_sparksfari1_contraint_carrier",
    "miss_eagle_contraint_carrier",
    "miss_syngo_contraint_carrier",
    "miss_chromepitf_contraint_carrier",
]
COLS_MISS_ALL = [
    "miss_contraint_carrier",
    "miss_hcndddom_xlinked_boyz_carrier",
    "miss_sparksfari1_carrier",
    "miss_eagle_carrier",
    "miss_syngo_carrier",
    "miss_chromepitf_carrier",
]
LABELS_MISS_CONSTRAINED = {
    "diag_genetic": "Returnable causative\nvariants",
    "miss_contraint_carrier": "Constrained",
    "miss_hcndddom_xlinked_boyz_contraint_carrier": "Constrained\nHCNDD dominant and\nrec X-linked",
    "miss_sparksfari1_contraint_carrier": "Constrained\nSPARK SFARI1",
    "miss_eagle_contraint_carrier": "Constrained\nEAGLE d&s",
    "miss_syngo_contraint_carrier": "Constrained\nSynGO",
    "miss_chromepitf_contraint_carrier": "Constrained\nChromEpiTF",
}
LABELS_MISS_ALL = {
    "miss_contraint_carrier": "Constraint",
    "miss_hcndddom_xlinked_boyz_carrier": "HCNDD dominant and\nrec X-linked",
    "miss_sparksfari1_carrier": "SPARK_SFARI1",
    "miss_eagle_carrier": "EAGLE d&s",
    "miss_syngo_carrier": "SynGO",
    "miss_chromepitf_carrier": "ChromEpiTF",
}

# --- Extended label maps (hg38: + metabolic + syngo_chromepitf) ---
LABELS_CONSTRAINED_EXT = {
    **LABELS_CONSTRAINED,
    "dellof_metabolic_contraint_carrier": "Constrained\nMetabolic",
}
LABELS_ALL_GENES_EXT = {
    **LABELS_ALL_GENES,
    "dellof_metabolic_carrier": "Metabolic",
}
LABELS_DELLOFMISS_CONSTRAINED_EXT = {
    **LABELS_DELLOFMISS_CONSTRAINED,
    "dellofmiss_metabolic_contraint_carrier": "Constrained\nMetabolic",
}
LABELS_DELLOFMISS_ALL_EXT = {
    **LABELS_DELLOFMISS_ALL,
    "dellofmiss_metabolic_carrier": "Metabolic",
}
