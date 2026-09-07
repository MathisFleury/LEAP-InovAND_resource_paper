# =============================================================================
# Unified project configuration — single source of truth
# =============================================================================
# One config for the whole LEAP-INOVAND project. Every section's _config.py /
# preprocessing/config.py imports from here (mostly as a thin re-export stub),
# so shared paths, palettes, label maps, gene-column lists and the cross-modality
# analysis policy live in ONE place.
#
# Design notes
# ------------
# * Section-relative OUTPUT dirs (outputs/figures, outputs/tables, ...) are NOT
#   defined here — each section computes its own from its own directory (use the
#   `outputs()` helper below, or keep the two/three lines in the section stub).
#   A single flat FIGURES_DIR would otherwise collide across sections.
# * Covariate policies differ by modality, so they are namespaced
#   (FUNC_* vs EEG_*); a stub aliases the bare name it needs.
# * Curated is the priority regime (see CLAUDE.md): CLUSTERS_DEFAULT points at
#   the curated k-means labels; the frozen paper labels stay available as
#   CLUSTERS_KMEANS_FROZEN.
# =============================================================================

import os

# --- Anchors -----------------------------------------------------------------
PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))          # LEAP-InovAND_resource
LIB = os.path.normpath(os.path.join(PROJECT_DIR, ".."))           # $LIB (parent of repo)
IMG5 = "/Volumes/Imaging5/EEG_MRI-MF"                              # $IMG5 (mounted volume)

_EEGMRI = os.path.join(LIB, "eeg_mri-pipeline")
_DATASET_PAPER = os.path.join(_EEGMRI, "results", "dataset_paper", "dataframes")
_GENETICS_IMG5 = os.path.join(IMG5, "ALL", "results", "tabular", "genetics")


# --- Matplotlib defaults -----------------------------------------------------
def apply_mpl_defaults():
    """Project-wide matplotlib/seaborn defaults (Helvetica, editable PDF text)."""
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import seaborn as sns
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["font.sans-serif"] = ["Helvetica", "Arial", "DejaVu Sans", "Liberation Sans"]
    mpl.rcParams["pdf.fonttype"] = 42
    sns.set_context("paper")


# Applied at import for backward compatibility (the old _config.py did the same).
try:
    apply_mpl_defaults()
except Exception:
    pass  # headless / matplotlib missing — scripts that plot will import it themselves


# --- Output-dir helper -------------------------------------------------------
def outputs(section_dir, curated=False, name="outputs"):
    """Return (OUTPUT_BASE, FIGURES_DIR, TABLES_DIR) for a section.

    section_dir : the section root (…/2_genetic_analysis), typically
                  os.path.dirname(scripts_dir).
    curated     : if True, nest under outputs/curated/ (priority regime).
    """
    base = os.path.join(section_dir, name, *(("curated",) if curated else ()))
    return base, os.path.join(base, "figures"), os.path.join(base, "tables")


# =============================================================================
# CLUSTER LABELS
# =============================================================================
# Curated is the priority regime → CLUSTERS_DEFAULT = curated k-means.
CLUSTERS_KMEANS_FROZEN = os.path.join(_DATASET_PAPER, "df_clusters_complete_kmeans.csv")
CLUSTERS_WARD_FROZEN = os.path.join(_DATASET_PAPER, "df_clusters_complete.csv")
CLUSTERS_CURATED = os.path.join(
    PROJECT_DIR, "1_clustering", "outputs", "tables",
    "individuals_metrics_with_clusters_curated.csv",
)
CLUSTERS_DEFAULT = CLUSTERS_CURATED
# Legacy aliases used by existing scripts (kept so imports don't break).
CLUSTERS_FILE = CLUSTERS_KMEANS_FROZEN                 # genetic/clinical roster (frozen)
CLUSTERS_CURATED_FILE = CLUSTERS_CURATED               # used by 10_clinical
CLUSTERS_CURATED_KMEANS_FILE = CLUSTERS_CURATED        # used by 3_cluster_genetic
CLUSTERS_LEGACY_FILE = os.path.join(
    PROJECT_DIR, "1_clustering", "outputs", "tables", "cluster_assignments.csv",
)

# =============================================================================
# CLINICAL / SHARED DATA FILES
# =============================================================================
INDIVIDUALS_METRICS = os.path.join(LIB, "imaging2genet", "0_input", "dataframes",
                                   "individuals_metrics.tsv")

# Curated per-cohort clinical TSVs — broad demographics roster (age_yrs, Sex,
# population_group, population, EEG_ID). Far larger than the clustering table, so
# preferred for demographics; the cluster table only adds the k-means Cluster.
LEAP_CLINICAL_CURATED = os.path.join(
    IMG5, "LEAP", "_clinical_data", "curated", "LEAP_clinical_curated_t1.tsv")
INOVAND_CLINICAL_CURATED = os.path.join(
    IMG5, "INOVAND", "_clinical_data", "curated", "INOVAND_clinical_curated_merged.tsv")

# =============================================================================
# GENETICS DATA FILES
# =============================================================================
GENELIST_PATH = os.path.join(_EEGMRI, "ressources", "genetics",
                             "hgnc_complete_set_20250424_GeneList_Updated20260225.txt")
LOEUF_SCORES = os.path.join(PROJECT_DIR, "script_zakaria",
                            "gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz")
DEL_FILE = os.path.join(PROJECT_DIR, "script_zakaria", "SV_DEL_annotation_12juin2025_ZM.tsv")
MISSING_DENOVO_CNV = os.path.join(PROJECT_DIR, "script_zakaria", "missing_denovo_cnv.csv")
SNV_FILE = os.path.join(PROJECT_DIR, "script_zakaria", "slivar_all_ensg_annotated1 3.xlsx")
DIAG_FILE = os.path.join(_GENETICS_IMG5, "diagnostic_clinic_JM",
                         "diag_individuals_combined_260408.csv")
CARRIER_PREPROCESS = os.path.join(
    _DATASET_PAPER, "df_carrier_genelist_DEL_LOF_MISS_withalphamissense_sex_takenintoaccount.tsv")
PGS_FILE = os.path.join(_DATASET_PAPER, "InovAND-LEAP_pgs_SBayesRC_20250401.tsv")
CARRIER_HG38 = os.path.join(
    _GENETICS_IMG5, "df_carrier_genelist_DEL_LOF_MISS_withalphamissense_LOEUF_gnomadv2_grch38.tsv")
DUP_FILE = os.path.join(_GENETICS_IMG5, "dup_per_individual_paperset_V2.tsv")
LOF_HG38_FILE = os.path.join(_GENETICS_IMG5, "df_38_with_lof_carrier_status.tsv")
VARIANTS_HG38_FILE = os.path.join(_GENETICS_IMG5, "df_variants_GRCh38.tsv")
# gnomad v4 pre-annotated carrier matrix (already IDD/NT spelling — do NOT relabel).
CARRIER_V4 = os.path.join(
    _GENETICS_IMG5, "df_carrier_genelist_DEL_LOF_MISS_withalphamissense_LOEUF_gnomadv4.tsv")
PGS_PAN_V4 = os.path.join(_GENETICS_IMG5, "InovAND-LEAP.pgs.panancestry.v4.tsv")
PGS_EUR_V4 = os.path.join(_GENETICS_IMG5, "InovAND-LEAP.pgs.europeans.v4.tsv")

# =============================================================================
# ANATOMICAL (MRI) DATA FILES
# =============================================================================
# CANONICAL anat regime = curated FreeSurfer (QC + ComBat + age/sex/eTIV/euler
# regression, z-scored), + CURATED clinical. Use MRI_CURATED for new anatomical
# analyses; take clinical features (SRS/IQ/RBS-R/VABS) from the CURATED cluster
# table (CLUSTERS_CURATED), NOT from any embedded clinical columns.
#
# Column naming in this table is raw FreeSurfer 8.1 stats — load it via 05's
# `load_curated_mri_with_canonical_id()` (4_anatomical_analysis/),
# which renames to tidy lh_/rh_*_{thickness,area,grayvol} + subcortical and
# attaches a canonical ID.
# Switch qc1<->qc12 by changing the filename here. As of 2026-07 the regression
# is fit on the FULL sample (not NT) — see the pipeline_ndd FIT_REFERENCE_COL note
# and 5_cluster_anatomical_analysis report §2.5. Override with MRI_NDD_FILE to
# point analyses at a variant (e.g. the `_noeuler` sensitivity table).
# Built in-repo by 4_anatomical_analysis/preprocessing/run_pipeline.py (2026-09
# port) into its own outputs/ dir, matching 6_functional_analysis/ and
# 8_eeg_analysis/'s preprocessing/outputs/ pattern.
MRI_CURATED = os.environ.get("MRI_NDD_FILE") or os.path.join(
    PROJECT_DIR, "4_anatomical_analysis", "preprocessing", "outputs",
    "freesurfer_zscore_qc1_combat_regress.tsv")
MRI_FILE = MRI_CURATED  # legacy alias (curated-pipeline scripts import MRI_FILE)

# FROZEN / paper ComBat table (imaging2genet). Region columns use the tidy
# lh_/rh_*_{thickness,area,volume} names AND embed frozen clinical columns.
# Paper reproduction only — it is NOT the canonical anat and its embedded
# clinical is frozen (use curated instead). Cortical volume suffix is _volume.
MRI_ANAT_COMBAT_FROZEN = os.path.join(LIB, "imaging2genet", "0_input", "dataframes",
                                      "MRI_ANAT_INOVAND_LEAP_COMBAT.tsv")

# ggseg aseg (subcortical) labels are CamelCase and mostly match FreeSurfer
# names; the only remaps needed are Thalamus -> *-Thalamus-Proper and
# 3rd/4th-Ventricle -> x3rd/x4th-ventricle. Do NOT lowercase aseg labels.
ASEG_LABEL_FIXES = {
    "Left-Thalamus": "Left-Thalamus-Proper", "Right-Thalamus": "Right-Thalamus-Proper",
    "3rd-Ventricle": "x3rd-ventricle", "4th-Ventricle": "x4th-ventricle",
}

# =============================================================================
# PALETTES & GROUP ORDERS
# =============================================================================
PALETTE_CLUSTERS = {"NT": "#C1C2BC", "C1": "#7A8B47", "C2": "#ff9fa0",
                    "C3": "#e7ba52", "IDD": "#D8A4CB"}
ORDER_CLUSTERS = ["NT", "C1", "C2", "C3", "IDD"]

PALETTE_POPULATION1 = {"Relatives": "#9AD5D3", "IDD": "#D8A4CB", "NT": "#C1C2BC",
                       "Autism with IDD": "#324095", "Autism without IDD": "#5CAEE1"}
ORDER_POPULATION1 = ["Autism without IDD", "Autism with IDD", "IDD", "Relatives", "NT"]

PALETTE_POPULATIONS = {"Autism": "#8991FA", "Relatives": "#9AD5D3",
                       "IDD": "#D8A4CB", "NT": "#C1C2BC"}
ORDER_POPULATIONS = ["Autism", "IDD", "Relatives", "NT"]

PALETTE_FREQ = {"Autism": "#8991FA", "Undiagnosed siblings": "#B0D5D3",
                "Undiagnosed parents": "#9AD5D3", "Undiagnosed 2+": "#9AD5D3",
                "Relatives": "#9AD5D3", "IDD": "#D8A4CB", "NT": "#C1C2BC",
                "Autism with IDD": "#324095", "Autism without IDD": "#5CAEE1"}
ORDER_FREQ = ["Autism", "IDD", "Undiagnosed siblings", "Undiagnosed parents", "NT"]
ORDER_OR = ["Autism", "IDD", "Undiagnosed siblings", "Undiagnosed parents"]

PALETTE_ANCESTRY = {"AFR": "#DB5E57", "AMR": "#DBD057", "CSA": "#74DB57",
                    "EAS": "#58DBAA", "EUR": "#579BDB", "MID": "#8557DB", "UNK": "#A9B3C7"}
ORDER_ANCESTRY = ["AFR", "AMR", "CSA", "EAS", "EUR", "MID", "UNK"]

# Group colours used by clinical/EEG violin scripts.
COLORS = {"Autism": "#8991FA", "Relatives": "#9AD5D3", "IDD": "#D8A4CB", "TD": "#C1C2BC",
          "Autism with IDD": "#324095", "Autism without IDD": "#5CAEE1", "NDD": "#8991FA"}
COLOR_AUTISM_MEDIAN = "#8991FA"
COLOR_NT_MEDIAN = "#C1C2BC"

# =============================================================================
# GENE-LIST COLUMN & LABEL MAPPINGS (genetic + cluster-genetic sections)
# =============================================================================
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

COLS_DUP_CONSTRAINED = [
    "dup_any_constraint_carrier", "dup_hcnddv7_constraint_carrier",
    "dup_sfari1_constraint_carrier", "dup_eagle_constraint_carrier",
    "dup_syngo_constraint_carrier", "dup_chromepitf_constraint_carrier",
]
LABELS_DUP_CONSTRAINED = {
    "dup_any_constraint_carrier": "Constrained",
    "dup_hcnddv7_constraint_carrier": "Constrained\nHCNDD dominant and\nrec X-linked",
    "dup_sfari1_constraint_carrier": "Constrained\nSPARK SFARI1",
    "dup_eagle_constraint_carrier": "Constrained\nEAGLE d&s",
    "dup_syngo_constraint_carrier": "Constrained\nSynGO",
    "dup_chromepitf_constraint_carrier": "Constrained\nChromEpiTF",
}
COLS_DUP_ALL = [
    "dup_any_constraint_carrier", "dup_hcnddv7_carrier", "dup_sfari1_carrier",
    "dup_eagle_carrier", "dup_syngo_carrier", "dup_chromepitf_carrier",
]
LABELS_DUP_ALL = {
    "dup_any_constraint_carrier": "Constraint",
    "dup_hcnddv7_carrier": "HCNDD dominant and\nrec X-linked",
    "dup_sfari1_carrier": "SPARK_SFARI1",
    "dup_eagle_carrier": "EAGLE d&s",
    "dup_syngo_carrier": "SynGO",
    "dup_chromepitf_carrier": "ChromEpiTF",
}

COMBO_PANEL_KEYS = ["constrained", "hcndd", "sparksfari1", "eagle", "syngo", "chromepitf"]
COMBO_PANEL_LABELS = {
    "constrained": "Constrained",
    "hcndd": "Constrained\nHCNDD dominant and\nrec X-linked",
    "sparksfari1": "Constrained\nSPARK SFARI1",
    "eagle": "Constrained\nEAGLE d&s",
    "syngo": "Constrained\nSynGO",
    "chromepitf": "Constrained\nChromEpiTF",
}
_DUP_ANY_CONSTRAINED = [
    "dup_hcnddv7_constraint_carrier", "dup_sfari1_constraint_carrier",
    "dup_eagle_constraint_carrier", "dup_syngo_constraint_carrier",
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
COMBO_SOURCES_DELLOFDUPMISS_HG38 = {
    "constrained": ["dellof_contraint_carrier", "miss_contraint_carrier"] + _DUP_ANY_CONSTRAINED,
    "hcndd": ["dellof_hcndddom_xlinked_boyz_contraint_carrier",
              "miss_hcndddom_xlinked_boyz_contraint_carrier", "dup_hcnddv7_constraint_carrier"],
    "sparksfari1": ["dellof_sparksfari1_contraint_carrier",
                    "miss_sparksfari1_contraint_carrier", "dup_sfari1_constraint_carrier"],
    "eagle": ["dellof_eagle_contraint_carrier", "miss_eagle_contraint_carrier",
              "dup_eagle_constraint_carrier"],
    "syngo": ["dellof_syngo_contraint_carrier", "miss_syngo_contraint_carrier",
              "dup_syngo_constraint_carrier"],
    "chromepitf": ["dellof_chromepitf_contraint_carrier", "miss_chromepitf_contraint_carrier",
                   "dup_chromepitf_constraint_carrier"],
}

LABELS_PGS = {
    "autism_grove2019": "Autism", "adhd_demontis2023": "ADHD",
    "int_savage2018": "Intelligence", "anxiety_purves2020": "Anxiety",
    "mdd_meng2024": "MDD", "ptsd_nievergelt2024": "PTSD",
    "ukbiobank_watanabe2019/felt_loved": "Felt loved as a child",
    "age_diagnosis_zhang2025/early": "Early autism diagnosis",
    "age_diagnosis_zhang2025/late": "Late autism diagnosis",
}
PGS_TRAITS = [
    "autism_grove2019", "age_diagnosis_zhang2025/early", "ptsd_nievergelt2024",
    "anxiety_purves2020", "mdd_meng2024", "adhd_demontis2023",
    "age_diagnosis_zhang2025/late", "int_savage2018", "ukbiobank_watanabe2019/felt_loved",
]

COLS_DELLOF_CONSTRAINED = [
    "diag_genetic", "dellof_contraint_carrier",
    "dellof_hcndddom_xlinked_boyz_contraint_carrier", "dellof_sparksfari1_contraint_carrier",
    "dellof_eagle_contraint_carrier", "dellof_syngo_contraint_carrier",
    "dellof_chromepitf_contraint_carrier",
]
COLS_DELLOF_ALL = [
    "dellof_contraint_carrier", "dellof_hcndddom_xlinked_boyz_carrier",
    "dellof_sparksfari1_carrier", "dellof_eagle_carrier", "dellof_syngo_carrier",
    "dellof_chromepitf_carrier",
]
COLS_DELLOFMISS_CONSTRAINED = [
    "diag_genetic", "dellofmiss_contraint_carrier",
    "dellofmiss_hcndddom_xlinked_boyz_contraint_carrier", "dellofmiss_sparksfari1_contraint_carrier",
    "dellofmiss_eagle_contraint_carrier", "dellofmiss_syngo_contraint_carrier",
    "dellofmiss_chromepitf_contraint_carrier",
]
COLS_DELLOFMISS_ALL = [
    "dellofmiss_contraint_carrier", "dellofmiss_hcndddom_xlinked_boyz_carrier",
    "dellofmiss_sparksfari1_carrier", "dellofmiss_eagle_carrier",
    "dellofmiss_syngo_carrier", "dellofmiss_chromepitf_carrier",
]
COLS_DELLOF_CONSTRAINED_WOPLI = [
    "diag_genetic", "dellof_contraint_wo_PLI_carrier",
    "dellof_hcndddom_xlinked_boyz_contraint_wo_PLI_carrier",
    "dellof_sparksfari1_contraint_wo_PLI_carrier", "dellof_eagle_contraint_wo_PLI_carrier",
    "dellof_syngo_contraint_wo_PLI_carrier", "dellof_chromepitf_contraint_wo_PLI_carrier",
]
LABELS_DELLOF_CONSTRAINED_WOPLI = {
    "diag_genetic": "Returnable causative\nvariants",
    "dellof_contraint_wo_PLI_carrier": "Constrained",
    "dellof_hcndddom_xlinked_boyz_contraint_wo_PLI_carrier": "Constrained\nHCNDD dominant and\nrec X-linked",
    "dellof_sparksfari1_contraint_wo_PLI_carrier": "Constrained\nSPARK SFARI1",
    "dellof_eagle_contraint_wo_PLI_carrier": "Constrained\nEAGLE d&s",
    "dellof_syngo_contraint_wo_PLI_carrier": "Constrained\nSynGO",
    "dellof_chromepitf_contraint_wo_PLI_carrier": "Constrained\nChromEpiTF",
}
COLS_DELLOFMISS_CONSTRAINED_WOPLI = [
    "diag_genetic", "dellofmiss_contraint_wo_PLI_carrier",
    "dellofmiss_hcndddom_xlinked_boyz_contraint_wo_PLI_carrier",
    "dellofmiss_sparksfari1_contraint_wo_PLI_carrier", "dellofmiss_eagle_contraint_wo_PLI_carrier",
    "dellofmiss_syngo_contraint_wo_PLI_carrier", "dellofmiss_chromepitf_contraint_wo_PLI_carrier",
]
LABELS_DELLOFMISS_CONSTRAINED_WOPLI = {
    "diag_genetic": "Returnable causative\nvariants",
    "dellofmiss_contraint_wo_PLI_carrier": "Constrained",
    "dellofmiss_hcndddom_xlinked_boyz_contraint_wo_PLI_carrier": "Constrained\nHCNDD dominant and\nrec X-linked",
    "dellofmiss_sparksfari1_contraint_wo_PLI_carrier": "Constrained\nSPARK SFARI1",
    "dellofmiss_eagle_contraint_wo_PLI_carrier": "Constrained\nEAGLE d&s",
    "dellofmiss_syngo_contraint_wo_PLI_carrier": "Constrained\nSynGO",
    "dellofmiss_chromepitf_contraint_wo_PLI_carrier": "Constrained\nChromEpiTF",
}

# Extended (hg38: + metabolic)
COLS_DELLOF_CONSTRAINED_EXT = COLS_DELLOF_CONSTRAINED + ["dellof_metabolic_contraint_carrier"]
COLS_DELLOF_ALL_EXT = COLS_DELLOF_ALL + ["dellof_metabolic_carrier"]
COLS_DELLOFMISS_CONSTRAINED_EXT = COLS_DELLOFMISS_CONSTRAINED + ["dellofmiss_metabolic_contraint_carrier"]
COLS_DELLOFMISS_ALL_EXT = COLS_DELLOFMISS_ALL + ["dellofmiss_metabolic_carrier"]
LABELS_CONSTRAINED_EXT = {**LABELS_CONSTRAINED, "dellof_metabolic_contraint_carrier": "Constrained\nMetabolic"}
LABELS_ALL_GENES_EXT = {**LABELS_ALL_GENES, "dellof_metabolic_carrier": "Metabolic"}
LABELS_DELLOFMISS_CONSTRAINED_EXT = {**LABELS_DELLOFMISS_CONSTRAINED, "dellofmiss_metabolic_contraint_carrier": "Constrained\nMetabolic"}
LABELS_DELLOFMISS_ALL_EXT = {**LABELS_DELLOFMISS_ALL, "dellofmiss_metabolic_carrier": "Metabolic"}

# Miss-only (hg38 where dellofmiss unavailable)
COLS_MISS_CONSTRAINED = [
    "diag_genetic", "miss_contraint_carrier", "miss_hcndddom_xlinked_boyz_contraint_carrier",
    "miss_sparksfari1_contraint_carrier", "miss_eagle_contraint_carrier",
    "miss_syngo_contraint_carrier", "miss_chromepitf_contraint_carrier",
]
COLS_MISS_ALL = [
    "miss_contraint_carrier", "miss_hcndddom_xlinked_boyz_carrier", "miss_sparksfari1_carrier",
    "miss_eagle_carrier", "miss_syngo_carrier", "miss_chromepitf_carrier",
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

# =============================================================================
# CROSS-MODALITY ANALYSIS POLICY  (namespaced — covariates differ by modality)
# =============================================================================
# Reference sample for nuisance regression + z-score = FULL sample (DX_COL=None),
# matching the anat pipeline. We no longer z-score against the NT/control group.
DX_COL = None
CONTROL_VALUE = None

# EEG (8_/9_): regress age/age²/sex, then full-sample z-score.
EEG_REGRESS_COVARS = ["age_yrs", "age_yrs_sq", "Sex"]

# Functional (6_/6-2): regress age/age²/sex/mean_fd; preserve age/age²/sex in ComBat.
# `cohort` is NOT included: the ComBat batch (scanner) is perfectly nested within
# cohort (each scanner belongs to one cohort), so scanner-level ComBat already
# removes cohort effects, and cohort is collinear with the batch design — keeping
# it as a preserved covariate is rank-deficient/no-op, and regressing it is
# redundant after harmonisation.
FUNC_REGRESS_COVARS = ["age_yrs", "age_yrs_sq", "Sex", "mean_fd"]
FUNC_COMBAT_COVARS = ["age_yrs", "age_yrs_sq", "Sex"]

# Anatomical curated-pipeline LOEUF analyses.
AUTISM_ONLY = False
N_PERM = 1000
SEED = 43
