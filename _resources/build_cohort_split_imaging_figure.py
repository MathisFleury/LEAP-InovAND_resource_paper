#!/usr/bin/env python3.11
"""
OPTIONAL supplementary figure — primary imaging results (Autism vs NT),
computed SEPARATELY for INOVAND and LEAP (not pooled), as a robustness check.

Standalone, on-demand -- not wired into any run_all.sh / numbered orchestrator.

Structural: reuses 4_anatomical_analysis/scripts/1_anatomical_mri_autism_nt.py's
own load_data()/_run_metric() (imported directly, no file changes to that
script), filtered to one cohort at a time, then renders each cohort's 4-panel
brain map via the existing shared ggseg engine
(_combined_brain_figure_shared.R, parameterised via its own ANAT_* env vars).

Functional: reuses 6_functional_analysis/concat/scripts/03_nbs_test.py
unchanged, run twice via subprocess against cohort-pre-filtered connectivity
CSVs (that script is already env-var parameterised: AUTISM_TD_FMRI_CSV,
AUTISM_TD_OUTPUT_DIR). Full N_PERM=5000, not shortened.

NOTE on scope (deliberate, see report): this does NOT replicate the full
Figure 6b/3 cortical+subcortical ggseg brain-map chain (08->10->09R->11R) --
that chain has its own fragile env-var contract (a prior bug is documented
inline in 6_functional_analysis/concat/run_all.sh) and is a lot of surface
to duplicate correctly for an optional/exploratory figure. Instead this uses
03_nbs_test.py's own built-in plots (null-distribution + network-block-count
matrix per direction), which are self-contained and already verified correct.

Working files land under _resources/cohort_split_tmp/ (gitignored via the
repo's blanket outputs/ rule? NO -- this is under _resources/, not an
outputs/ dir, so it is NOT auto-ignored; kept out of git explicitly below).

Output: _resources/Figure_Imaging_Cohort_Split.pdf
"""
import importlib.util
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import fitz  # PyMuPDF

ROOT = Path(__file__).resolve().parent.parent
PY = "/usr/local/bin/python3.11"
TMP = Path(__file__).resolve().parent / "cohort_split_tmp"
OUT_PDF = Path(__file__).resolve().parent / "Figure_Imaging_Cohort_Split.pdf"

# --------------------------------------------------------------------- structural
ANAT_SCRIPTS = ROOT / "4_anatomical_analysis" / "scripts"
ANAT_SHARED_R = ANAT_SCRIPTS / "_combined_brain_figure_shared.R"

# --------------------------------------------------------------------- functional
FUNC_SCRIPTS = ROOT / "6_functional_analysis" / "concat" / "scripts"
FUNC_CONN_CSV = (ROOT / "6_functional_analysis" / "concat" / "preprocessing" /
                  "outputs" / "nogsr_concat" / "df_conn_cohort_norm_nogsr_concat.csv")
N_PERM = os.environ.get("N_PERM", "5000")


def _load_anat_module():
    spec = importlib.util.spec_from_file_location(
        "anat_autism_nt", ANAT_SCRIPTS / "1_anatomical_mri_autism_nt.py")
    m = importlib.util.module_from_spec(spec)
    sys.path.insert(0, str(ANAT_SCRIPTS))  # so its own `import _config` resolves
    spec.loader.exec_module(m)
    return m


def build_structural(cohort_label, cohort_prefixes, out_dir):
    """cohort_prefixes: tuple of `cohort` column prefixes to keep (e.g. ('LEAP',))."""
    m = _load_anat_module()
    df = m.load_data()
    mask = df["cohort"].astype(str).str.startswith(tuple(cohort_prefixes))
    sub = df[mask].copy()
    n_a = int((sub["PopulationS1"] == "Autism").sum())
    n_nt = int((sub["PopulationS1"] == "NT").sum())
    print(f"[structural/{cohort_label}] N = {len(sub)}  (Autism={n_a}, NT={n_nt})")

    r_input_dir = out_dir / "r_input_files"
    r_input_dir.mkdir(parents=True, exist_ok=True)
    for metric, atlas, fname in [
        ("thickness", "dk", "t_stat_anat_dk_thickness_mri_autism_vs_control.csv"),
        ("area", "dk", "t_stat_anat_dk_area_mri_autism_vs_control.csv"),
        ("grayvol", "dk", "t_stat_anat_dk_grayvol_mri_autism_vs_control.csv"),
        ("volume", "aseg", "t_stat_anat_aseg_volume_mri_autism_vs_control.csv"),
    ]:
        res = m._run_metric(sub, metric, atlas)
        if res.empty:
            print(f"  [skip] {metric}: no testable region")
            continue
        res.to_csv(r_input_dir / fname, index=False)
        n_fdr = int((res["p_fdr"] < 0.05).sum())
        print(f"  {metric:10s}: {len(res):3d} regions, FDR<0.05 = {n_fdr}")

    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, ANAT_INPUT_DIR=str(r_input_dir), ANAT_OUTPUT_DIR=str(fig_dir),
               ANAT_COMPARISON=f"Autism vs NT -- {cohort_label} only")
    subprocess.run(["Rscript", str(ANAT_SHARED_R)], cwd=str(ANAT_SCRIPTS), env=env, check=True)
    pdf = fig_dir / "combined_brain_maps_fdr_cohensd.pdf"
    if not pdf.exists():
        raise SystemExit(f"expected {pdf} was not produced")
    return pdf, n_a, n_nt


def build_functional(cohort_label, cohort_value, out_dir):
    out_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(FUNC_CONN_CSV, low_memory=False)
    sub = df[df["cohort"] == cohort_value].copy()
    n_a = int((sub["population_group"] == "Autism").sum())
    n_nt = int((sub["population_group"].replace("TD", "NT") == "NT").sum())
    print(f"[functional/{cohort_label}] N = {len(sub)}  (Autism={n_a}, NT={n_nt})")

    conn_csv = out_dir / f"df_conn_{cohort_label}.csv"
    sub.to_csv(conn_csv, index=False)

    nbs_dir = out_dir / "nbs"
    nbs_dir.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, AUTISM_TD_FMRI_CSV=str(conn_csv), AUTISM_TD_OUTPUT_DIR=str(nbs_dir),
               N_PERM=N_PERM)
    subprocess.run([PY, "03_nbs_test.py"], cwd=str(FUNC_SCRIPTS), env=env, check=True)

    summary = {}
    for d in ("hyper", "hypo"):
        csv = nbs_dir / "tables" / f"nbs_components_{d}.csv"
        if csv.exists():
            r = pd.read_csv(csv)
            sig = r[r.p_fwer < 0.05]
            summary[d] = (len(sig), int(r.n_edges.max()) if len(r) else 0,
                          float(r.p_fwer.min()) if len(r) else float("nan"))
    return nbs_dir / "figures", summary, n_a, n_nt


def _page_to_pixmap_rect(doc_out_page, src_pdf_path, rect):
    src = fitz.open(src_pdf_path)
    doc_out_page.show_pdf_page(rect, src, 0)
    src.close()


def main():
    if TMP.exists():
        shutil.rmtree(TMP)
    TMP.mkdir(parents=True)

    cohorts = [
        ("INOVAND", ("INOVAND", "INFOR"), "INOVAND"),
        ("LEAP", ("LEAP",), "LEAP"),
    ]

    struct_results = {}
    func_results = {}
    for label, anat_prefixes, func_value in cohorts:
        sdir = TMP / f"structural_{label}"
        struct_results[label] = build_structural(label, anat_prefixes, sdir)
        fdir = TMP / f"functional_{label}"
        func_results[label] = build_functional(label, func_value, fdir)

    # ---------------- report ----------------
    print("\n" + "=" * 70)
    print("SUMMARY")
    for label, _, _ in cohorts:
        _, n_a_s, n_nt_s = struct_results[label]
        _, summary, n_a_f, n_nt_f = func_results[label]
        print(f"\n{label}: structural N=Autism {n_a_s}/NT {n_nt_s}; "
              f"functional N=Autism {n_a_f}/NT {n_nt_f}")
        for d, (nsig, maxedges, minp) in summary.items():
            print(f"  functional {d}: {nsig} sig component(s) (FWER<0.05), "
                  f"largest {maxedges} edges, min p_fwer={minp:.4f}")

    # ---------------- compose ----------------
    # Label band is tall enough for the letter and the descriptive text on
    # separate lines -- side-by-side placement collided with the next
    # column's label for long strings (caught by rendering and inspecting
    # the composite before finalising).
    CM = 28.3465
    MARGIN = 0.7 * CM
    GAP = 0.5 * CM
    LABEL_BAND = 1.1 * CM
    col_w = 8.6 * CM
    struct_h = col_w * 1.85   # tall 4-panel stacked brain maps
    func_h = 6.0 * CM         # null-distribution pair, shorter

    page_w = 2 * col_w + GAP + 2 * MARGIN
    page_h = MARGIN + LABEL_BAND + struct_h + GAP + LABEL_BAND + func_h + MARGIN

    doc = fitz.open()
    page = doc.new_page(width=page_w, height=page_h)
    labels = ["a", "b"]
    y = MARGIN
    x_positions = [MARGIN, MARGIN + col_w + GAP]

    # row 1: structural, INOVAND | LEAP
    for i, (label, _, _) in enumerate(cohorts):
        pdf_path, n_a, n_nt = struct_results[label][0], struct_results[label][1], struct_results[label][2]
        page.insert_text((x_positions[i], y + 16), labels[i],
                          fontsize=18, fontname="Helvetica-Bold")
        page.insert_text((x_positions[i], y + LABEL_BAND - 6),
                          f"{label} structural (Autism n={n_a}, NT n={n_nt})",
                          fontsize=10, fontname="Helvetica")
        _page_to_pixmap_rect(page, pdf_path,
                              fitz.Rect(x_positions[i], y + LABEL_BAND, x_positions[i] + col_w,
                                        y + LABEL_BAND + struct_h))
    y += LABEL_BAND + struct_h + GAP

    # row 2: functional, INOVAND | LEAP -- null-distribution plots (hyper+hypo side by side)
    labels2 = ["c", "d"]
    for i, (label, _, _) in enumerate(cohorts):
        fig_dir, summary, n_a, n_nt = func_results[label]
        page.insert_text((x_positions[i], y + 16), labels2[i],
                          fontsize=18, fontname="Helvetica-Bold")
        page.insert_text((x_positions[i], y + LABEL_BAND - 6),
                          f"{label} functional NBS (Autism n={n_a}, NT n={n_nt})",
                          fontsize=10, fontname="Helvetica")
        sub_w = col_w / 2
        for j, d in enumerate(("hyper", "hypo")):
            p = fig_dir / f"nbs_null_distribution_{d}.pdf"
            if p.exists():
                _page_to_pixmap_rect(page, p,
                                      fitz.Rect(x_positions[i] + j * sub_w, y + LABEL_BAND,
                                                x_positions[i] + (j + 1) * sub_w, y + LABEL_BAND + func_h))
    doc.save(OUT_PDF)
    print(f"\nSaved: {OUT_PDF}  ({page_w / CM:.1f} x {page_h / CM:.1f} cm)")


if __name__ == "__main__":
    main()
