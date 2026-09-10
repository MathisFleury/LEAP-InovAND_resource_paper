# =============================================================================
# LEAP biological sample (biomaterials) collection table.
#
# Source: /Volumes/Imaging5/EEG_MRI-MF/LEAP/_clinical_data/raw/BiosamplesTable_forThomas.docx
# (transcribed verbatim below; no InovAND biosample table exists at this path,
# so this is LEAP-only). Same gt-table style as the MRI acquisition tables
# (gt_table_mri_LEAP.R). Location column dropped -- not needed for public release.
# =============================================================================

library(dplyr)
library(gt)

OUT_DIR <- "/Users/mfleury/POSTDOC/LIBRAIRY/LEAP-InovAND_resource/_resources"

biosamples <- tibble::tribble(
  ~Timepoint,                    ~Sample_Type,                        ~Kit_Volume,                                                                                                                          ~N_Samples, ~N_Participants,                                    ~Purpose,               ~Storage,
  "LEAP Timepoints 1 and 2", "Whole Blood (EDTA)",                "1 x 10ml EDTA tubes",                                                                                                                "2,112",    "560 + 64 Twins/Siblings + 408 Mothers + 253 Fathers", "Genetics",             "-80°C",
  "LEAP Timepoints 1 and 2", "Whole Blood (Tempus)",              "2 x 3ml TEMPUS tubes",                                                                                                                "1,202",    "561 + 63 Twins/Siblings",                             "RNA, Inflammation",   "-80°C",
  "LEAP Timepoints 1 and 2", "Blood Plasma",                      "1 x 10ml EDTA tubes to generate ~4 x 1.8ml screw-cap microcentrifuge tubes",                                                        "2,249",    "522 + 58 Twins/Siblings",                             "Serotonin, Inflammation", "-80°C",
  "LEAP Timepoints 1 and 2", "Extracted DNA",                     "8x96-well plates (neat DNA, blood); 11x96-well plates (200µl neat DNA, saliva); 11x96-well plates (remaining neat DNA, saliva); all not normalised", "-", "-",                                                   "Extracted DNA",       "-80°C",
  "LEAP Timepoints 1 and 2", "Saliva",                            "1 x 2.5ml OG-500 ORAgene saliva kits",                                                                                                "1,605",    "544 + 76 Twins/Siblings + 408 Mothers + 315 Fathers", "Genetics, Epigenetics", "Room temperature",
  "LEAP Timepoints 1 and 2", "Urine",                             "2 x 10ml urine tube",                                                                                                                 "3,665",    "618 + 408 Twins/Siblings + 432 Mothers + 321 Fathers","p-Cresol",             "-80°C",
  "LEAP Timepoints 1 and 2", "Hairs",                             "50 mL centrifuge tubes, MEF media + 50 mg/ml Gentamycin + 15 mM HEPES",                                                              ">2,539",   "311 + 102 Mothers + 35 Fathers",                      "iPSC",                 "-",
  "LEAP Timepoint 3",        "Whole Blood (EDTA)*",               "1 x 5ml minimum in 6ml EDTA tubes",                                                                                                   "204",      "130 + 22 Fathers + 52 Mothers",                       "Genetics",             "-80°C",
  "LEAP Timepoint 3",        "Whole Blood (EDTA)",                "1 x 1ml (approx.) in 6ml EDTA tubes",                                                                                                 "298",      "298",                                                  "Serotonin",            "-80°C",
  "LEAP Timepoint 3",        "Whole Blood (EDTA)",                "1 x 5ml minimum in 6ml EDTA tubes",                                                                                                   "250",      "250",                                                  "Biobank",              "-80°C",
  "LEAP Timepoint 3",        "Whole Blood (Tempus)*",             "2 x 3ml TEMPUS tubes",                                                                                                                "264",      "150",                                                  "RNA, Inflammation",   "-80°C",
  "LEAP Timepoint 3",        "Whole Blood (without anticoagulant)","1 x 8.5ml BD SST II Advance tubes",                                                                                                  "304",      "304",                                                  "Inflammation",         "-80°C",
  "LEAP Timepoint 3",        "Blood Plasma*",                     "~4 x 0.5ml into 2ml screw-cap microcentrifuge tubes",                                                                                 "371",      "103",                                                  "Inflammation",         "-80°C",
  "LEAP Timepoint 3",        "Blood Serum",                       "~5 x 0.5ml into 2ml screw-cap microcentrifuge tubes",                                                                                 "1,361",    "304",                                                  "Inflammation",         "-80°C",
  "LEAP Timepoint 3",        "Saliva**",                          "1 x 2.5ml OG-500 ORAgene saliva kits",                                                                                                "229",      "87 + 64 Fathers + 70 Mothers",                        "Genetics",             "Room temperature",
  "LEAP Timepoint 3",        "Stool",                             "OM-200 GeneGut kit",                                                                                                                  "239",      "239",                                                  "Microbiome",           "-80°C",
)

gt_tbl <- biosamples %>%
  gt(groupname_col = "Timepoint") %>%
  tab_header(title = md("**Supplementary Table. LEAP Biological Sample (Biomaterials) Collection**")) %>%
  cols_label(
    Sample_Type = "Sample Type", Kit_Volume = "Kit / Volume", N_Samples = "N Samples",
    N_Participants = "N Unique Participants", Purpose = "Purpose", Storage = "Storage Condition"
  ) %>%
  tab_source_note(source_note = "*Only if not already collected/processed at LEAP Timepoints 1/2. **Only if whole blood for genetics not available/processed.") %>%
  tab_options(table.width = pct(100), table.font.size = px(10))

gt::gtsave(gt_tbl, file.path(OUT_DIR, "Supp_Table_LEAP_Biosamples.pdf"))
cat("Saved Supp_Table_LEAP_Biosamples.pdf\n")
