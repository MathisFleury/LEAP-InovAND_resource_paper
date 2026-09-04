#!/usr/bin/env Rscript
# =============================================================================
# gt tables: autism vs NT MRI stats (REVISED anat + CURATED clinical)
# =============================================================================
# Publication-ready gt tables (t, p, q, Cohen's d; Left/Right side by side) for
#   - cortical thickness   - surface area   - subcortical volume
# Recomputed on the revised regressed FreeSurfer + curated clinical (autism vs
# NT), i.e. the per-region stats written by 1_anatomical_mri_autism_nt_v2.py
# (run with CURATED_CLINICAL=1). Adapted from
# eeg_mri-pipeline/.../tables_creation/gt_table_anat_creation.R.
#
# Input : ../outputs/figures/r_input_files/t_stat_anat_*.csv
#           (label, t_stat, p_val, cohens_d, p_fdr, p_bonf)
# Output: ../outputs/figures/table_{cortical_thickness,surface_area,subcortical_volume}.{pdf,html}
# =============================================================================
suppressMessages({ library(dplyr); library(tidyr); library(gt); library(stringr) })

args <- commandArgs(trailingOnly = FALSE)
sp <- sub("--file=", "", args[grep("--file=", args)])
script_dir <- if (length(sp) == 0) getwd() else dirname(normalizePath(sp))
base_dir <- file.path(dirname(script_dir), Sys.getenv("LOEUF_OUT_DIR", "outputs"), "figures")
r_input_dir <- file.path(base_dir, "r_input_files")

separate_hemispheres <- function(df) df %>%
  mutate(hemisphere = case_when(str_detect(label, "^lh_|^left-|^Left-") ~ "Left",
                                str_detect(label, "^rh_|^right-|^Right-") ~ "Right",
                                TRUE ~ "Bilateral"),
         roi_name = label %>% str_remove("^(lh_|rh_|left-|right-|Left-|Right-)") %>%
           str_to_title() %>% str_replace_all("_", " ")) %>%
  mutate(roi_name = case_when(
    str_detect(roi_name, "Thalamus Proper") ~ "Thalamus",
    str_detect(roi_name, "Ventral Dc") ~ "Ventral DC",
    str_detect(roi_name, "Cerebellum White Matter") ~ "Cerebellum White Matter",
    str_detect(roi_name, "Cerebellum Cortex") ~ "Cerebellum Cortex",
    str_detect(roi_name, "Lateral Ventricle") ~ "Lateral Ventricle",
    str_detect(roi_name, "X3rd") ~ "3rd Ventricle", str_detect(roi_name, "X4th") ~ "4th Ventricle",
    str_detect(roi_name, "Brain Stem") ~ "Brain Stem",
    str_detect(roi_name, "^Cc ") ~ str_replace(roi_name, "^Cc", "CC"),
    TRUE ~ roi_name))

fmt_p <- function(p) vapply(as.numeric(p), function(x)
  if (is.na(x)) NA_character_ else if (x < 1e-4) formatC(x, format = "e", digits = 2)
  else as.character(round(x, 5)), character(1))

make_table <- function(df, measure) {
  d <- separate_hemispheres(df) %>% filter(hemisphere != "Bilateral") %>%
    mutate(fdr = p_bonf, significant = !is.na(p_bonf) & p_bonf < 0.05)
  side <- function(dd, sfx) dd %>%
    transmute(roi_name, !!paste0("t", sfx) := round(t_stat, 2),
              !!paste0("p", sfx) := fmt_p(p_val), !!paste0("q", sfx) := fmt_p(fdr),
              !!paste0("d", sfx) := round(cohens_d, 3), !!paste0("sig", sfx) := significant)
  tab <- full_join(side(filter(d, hemisphere == "Left"), "_L"),
                   side(filter(d, hemisphere == "Right"), "_R"), by = "roi_name") %>%
    arrange(roi_name) %>% mutate(sep = "") %>%
    select(roi_name, t_L, p_L, q_L, d_L, sig_L, sep, t_R, p_R, q_R, d_R, sig_R)
  g <- tab %>% gt() %>%
    tab_header(title = md(sprintf("**%s — Autism vs NT** (revised anat, curated clinical)", measure))) %>%
    tab_spanner("Left Hemisphere", c(t_L, p_L, q_L, d_L)) %>%
    tab_spanner("Right Hemisphere", c(t_R, p_R, q_R, d_R)) %>%
    cols_label(roi_name = "ROI", t_L = "t", p_L = "p", q_L = "p_bonf", d_L = "d", sep = "",
               t_R = "t", p_R = "p", q_R = "p_bonf", d_R = "d") %>%
    tab_style(cell_borders(sides = "left", weight = px(2)), cells_body(columns = sep)) %>%
    tab_style(cell_text(weight = "bold"), cells_column_labels()) %>%
    tab_style(cell_text(weight = "bold"), cells_column_spanners()) %>%
    tab_style(cell_text(style = "italic"), cells_body(columns = roi_name)) %>%
    tab_style(cell_fill(color = "lightyellow"),
              cells_body(rows = tab$sig_L == TRUE, columns = c(t_L, p_L, q_L, d_L))) %>%
    tab_style(cell_fill(color = "lightyellow"),
              cells_body(rows = tab$sig_R == TRUE, columns = c(t_R, p_R, q_R, d_R))) %>%
    tab_options(table.font.size = px(9), column_labels.font.size = px(9),
                data_row.padding = px(2)) %>%
    tab_source_note("Bonferroni corrected; highlighted cells = p_bonf < 0.05. Positive d = Autism > NT.") %>%
    cols_hide(c(sig_L, sig_R))
  g
}

SPECS <- list(
  list(file = "t_stat_anat_dk_thickness_mri_autism_vs_control.csv", meas = "Cortical thickness", out = "table_cortical_thickness"),
  list(file = "t_stat_anat_dk_area_mri_autism_vs_control.csv",      meas = "Surface area",       out = "table_surface_area"),
  list(file = "t_stat_anat_aseg_volume_mri_autism_vs_control.csv",  meas = "Subcortical volume", out = "table_subcortical_volume"))

for (s in SPECS) {
  fp <- file.path(r_input_dir, s$file)
  if (!file.exists(fp)) { warning(paste("missing", fp)); next }
  g <- make_table(read.csv(fp, stringsAsFactors = FALSE), s$meas)
  # HTML never needs Chrome; PDF/PNG go through webshot2 -> chromote -> headless
  # Chrome, which can fail to launch (debug port). Write HTML always, attempt the
  # rasterised formats best-effort so a Chrome hiccup doesn't lose the table.
  html_out <- file.path(base_dir, paste0(s$out, ".html"))
  gtsave(g, html_out); cat("Saved:", html_out, "\n")
  for (ext in c("pdf", "png")) {
    out <- file.path(base_dir, paste0(s$out, ".", ext))
    ok <- tryCatch({ gtsave(g, out); TRUE }, error = function(e) {
      cat("WARN: could not render", ext, "(", conditionMessage(e), ")\n"); FALSE })
    if (ok) cat("Saved:", out, "\n")
  }
}
