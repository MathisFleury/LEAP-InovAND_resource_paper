#!/usr/bin/env Rscript
# =============================================================================
# gt tables: per-cluster MRI stats (Cluster vs NT), REVISED anat
# =============================================================================
# Same design as 4_anatomical_analysis/curated/scripts/16_gt_tables_anat_v2.R
# (t, p, p_fdr, Cohen's d; Left/Right side by side), one table set per
# cluster (C1/C2/C3) for
#   - cortical thickness   - surface area   - subcortical volume
# Reads the per-cluster stats written by 01_generate_cluster_mri_inputs_v2.py
# (cluster_vs_nt()). Uses FDR (not Bonferroni) for significance.
#
# Input : ../outputs/figures/r_input_files/t_stat_cluster_<C>_<metric>_<atlas>_mri_cluster_vs_td.csv
#           (label, t_stat, p_val, cohens_d, p_fdr)
# Output: ../outputs/figures/table_cluster_<C>_{cortical_thickness,surface_area,subcortical_volume}.{pdf,html,png}
# =============================================================================
suppressMessages({ library(dplyr); library(tidyr); library(gt); library(stringr) })

args <- commandArgs(trailingOnly = FALSE)
sp <- sub("--file=", "", args[grep("--file=", args)])
script_dir <- if (length(sp) == 0) getwd() else dirname(normalizePath(sp))
base_dir <- file.path(dirname(script_dir), Sys.getenv("CLUSTER_OUT_DIR", "outputs"), "figures")
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

make_table <- function(df, measure, cluster) {
  d <- separate_hemispheres(df) %>% filter(hemisphere != "Bilateral") %>%
    mutate(fdr = p_fdr, significant = !is.na(p_fdr) & p_fdr < 0.05)
  side <- function(dd, sfx) dd %>%
    transmute(roi_name, !!paste0("t", sfx) := round(t_stat, 2),
              !!paste0("p", sfx) := fmt_p(p_val), !!paste0("q", sfx) := fmt_p(fdr),
              !!paste0("d", sfx) := round(cohens_d, 3), !!paste0("sig", sfx) := significant)
  tab <- full_join(side(filter(d, hemisphere == "Left"), "_L"),
                   side(filter(d, hemisphere == "Right"), "_R"), by = "roi_name") %>%
    arrange(roi_name) %>% mutate(sep = "") %>%
    select(roi_name, t_L, p_L, q_L, d_L, sig_L, sep, t_R, p_R, q_R, d_R, sig_R)
  g <- tab %>% gt() %>%
    tab_header(title = md(sprintf("**%s — %s vs NT** (revised anat)", measure, cluster))) %>%
    tab_spanner("Left Hemisphere", c(t_L, p_L, q_L, d_L)) %>%
    tab_spanner("Right Hemisphere", c(t_R, p_R, q_R, d_R)) %>%
    cols_label(roi_name = "ROI", t_L = "t", p_L = "p", q_L = "p_fdr", d_L = "d", sep = "",
               t_R = "t", p_R = "p", q_R = "p_fdr", d_R = "d") %>%
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
    tab_source_note("FDR (Benjamini-Hochberg) corrected; highlighted cells = p_fdr < 0.05. Positive d = cluster > NT.") %>%
    cols_hide(c(sig_L, sig_R))
  g
}

MEASURES <- list(
  list(metric = "thickness", atlas = "dk",   meas = "Cortical thickness", out = "cortical_thickness"),
  list(metric = "area",      atlas = "dk",   meas = "Surface area",       out = "surface_area"),
  list(metric = "volume",    atlas = "aseg", meas = "Subcortical volume", out = "subcortical_volume"))

for (cluster in c("C1", "C2", "C3")) {
  for (m in MEASURES) {
    fp <- file.path(r_input_dir, sprintf("t_stat_cluster_%s_%s_%s_mri_cluster_vs_td.csv",
                                          cluster, m$metric, m$atlas))
    if (!file.exists(fp)) { warning(paste("missing", fp)); next }
    g <- make_table(read.csv(fp, stringsAsFactors = FALSE), m$meas, cluster)
    out_stub <- sprintf("table_cluster_%s_%s", cluster, m$out)
    # HTML never needs Chrome; PDF/PNG go through webshot2 -> chromote -> headless
    # Chrome, which can fail to launch (debug port). Write HTML always, attempt
    # the rasterised formats best-effort so a Chrome hiccup doesn't lose the table.
    html_out <- file.path(base_dir, paste0(out_stub, ".html"))
    gtsave(g, html_out); cat("Saved:", html_out, "\n")
    for (ext in c("pdf", "png")) {
      out <- file.path(base_dir, paste0(out_stub, ".", ext))
      ok <- tryCatch({ gtsave(g, out); TRUE }, error = function(e) {
        cat("WARN: could not render", ext, "(", conditionMessage(e), ")\n"); FALSE })
      if (ok) cat("Saved:", out, "\n")
    }
  }
}
