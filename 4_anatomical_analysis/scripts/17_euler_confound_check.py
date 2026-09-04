#!/usr/bin/env python3
"""
Image-quality (Euler number) confound check for the anatomical Autism-vs-NT
group differences — REVISION pipeline only.

Structural analogue of 6-2's 12_motion_confound_check.py. Reviewer ask: report
effect sizes of the association between image quality (Euler / surface defects)
and the imaging phenotypes that show group differences, to disambiguate whether
quality could be driving the Autism-vs-NT effects.

Quality convention: `defects = -mean_euler` (higher = worse reconstruction, so
the sign logic matches the FD analysis). Autism typically has more negative
Euler (more defects).

Produces (outputs/tables/qc_euler/ for the CSVs, outputs/figures/qc_euler/ for the plot):
  (1) Euler group difference (Autism vs NT): Welch t, Cohen's d.
  (2) defects <-> phenotype association for each FDR-significant region
      (thickness / area / subcortical volume): Pearson r pooled + within-group,
      sign compared to the group-effect direction (sign of Cohen's d). Same sign
      => candidate confound; opposite => quality not driving it.
  (3) Euler raincloud/violin by population group (standard naming; "Autism to
      exclude" folded into "Autism without IDD").

Phenotypes are already Euler-regressed in this pipeline, so near-zero
associations are the expected/reassuring outcome; this quantifies it explicitly.
"""
import os
import sys
from importlib import import_module
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind, pearsonr

sys.path.insert(0, str(Path(__file__).resolve().parent))
_m01 = import_module("1_anatomical_mri_autism_nt")
import _config  # noqa: E402

_SECTION = Path(__file__).resolve().parent.parent
RINPUT = _SECTION / "outputs" / "tables" / "r_input_files"
TABLES_OUT = _SECTION / "outputs" / "tables" / "qc_euler"
OUT = _SECTION / "outputs" / "figures" / "qc_euler"
TABLES_OUT.mkdir(parents=True, exist_ok=True)
OUT.mkdir(parents=True, exist_ok=True)
CURATED = Path("/Users/mfleury/POSTDOC/LIBRAIRY/LEAP-InovAND_resource/"
               "1_clustering/outputs/tables/individuals_metrics_with_clusters_curated.csv")

POP_REMAP = {"Autism to exclude": "Autism without IDD", "Autism": "Autism without IDD"}
GROUPS = ["NT", "Relatives", "IDD", "Autism without IDD", "Autism with IDD"]
# significant-region tables (metric, atlas file token, column suffix)
METRICS = [("thickness", "t_stat_anat_dk_thickness_mri_autism_vs_control.csv", "_thickness"),
           ("area",      "t_stat_anat_dk_area_mri_autism_vs_control.csv",      "_area"),
           ("volume",    "t_stat_anat_aseg_volume_mri_autism_vs_control.csv",  "")]


def cohens_d(a, b):
    na, nb = len(a), len(b)
    sd = np.sqrt(((na-1)*np.var(a, ddof=1) + (nb-1)*np.var(b, ddof=1)) / (na+nb-2))
    return (np.mean(a) - np.mean(b)) / sd


def main() -> int:
    df = _m01.load_data()                       # Autism/NT + tidy regional cols + mean_euler
    if "mean_euler" not in df.columns:
        print("ERROR: mean_euler not in loaded data"); return 1
    df["defects"] = -df["mean_euler"]           # higher = worse
    a = df["PopulationS1"] == "Autism"; n = df["PopulationS1"] == "NT"

    # (1) Euler group difference
    ea, en = df.loc[a, "mean_euler"].dropna().values, df.loc[n, "mean_euler"].dropna().values
    t, p = ttest_ind(ea, en, equal_var=False); d = cohens_d(ea, en)
    print(f"[Euler group diff] Autism {ea.mean():.1f} vs NT {en.mean():.1f} | "
          f"Welch t={t:.2f} p={p:.2e} Cohen's d={d:+.3f} (negative = autism worse)")
    pd.DataFrame([{"measure": "mean_euler", "autism_mean": ea.mean(), "nt_mean": en.mean(),
                   "welch_t": t, "p": p, "cohens_d": d, "n_autism": len(ea), "n_nt": len(en)}]
                 ).to_csv(TABLES_OUT / "euler_effect_sizes.csv", index=False)

    # (2) defects <-> phenotype for FDR-significant regions
    rows = []
    for metric, fname, suf in METRICS:
        fp = RINPUT / fname
        if not fp.exists():
            print(f"  {metric}: {fname} missing"); continue
        st = pd.read_csv(fp)
        sig = st[st["p_fdr"] < 0.05]
        for _, r in sig.iterrows():
            col = f"{r['label']}{suf}"
            if col not in df.columns:
                continue
            gdir = np.sign(r["cohens_d"])       # group-effect direction
            for scope, m in [("pooled", a | n), ("within_NT", n), ("within_Autism", a)]:
                x = df.loc[m, "defects"].values; y = df.loc[m, col].values
                ok = np.isfinite(x) & np.isfinite(y)
                if ok.sum() < 10:
                    continue
                rr, pr = pearsonr(x[ok], y[ok])
                rows.append({"metric": metric, "region": r["label"], "group_cohens_d": round(r["cohens_d"], 3),
                             "scope": scope, "r_defects_vs_pheno": round(rr, 3), "p": round(pr, 4), "n": int(ok.sum()),
                             "same_dir_as_group": bool(np.sign(rr) == gdir)})
    assoc = pd.DataFrame(rows)
    assoc.to_csv(TABLES_OUT / "euler_phenotype_association.csv", index=False)
    # summary per metric (within-Autism scope = the interpretable one)
    print("\n[defects <-> phenotype] FDR-significant regions (within-Autism scope):")
    for metric, _, _ in METRICS:
        w = assoc[(assoc.metric == metric) & (assoc.scope == "within_Autism")]
        if len(w) == 0:
            continue
        same = int(w.same_dir_as_group.sum())
        print(f"  {metric:9s}: {len(w)} sig regions | mean|r|={w.r_defects_vs_pheno.abs().mean():.3f} "
              f"max|r|={w.r_defects_vs_pheno.abs().max():.3f} | same-direction-as-group: {same}/{len(w)}")

    # (3) Euler violin by population
    make_violin()
    print(f"\nSaved: {TABLES_OUT}/euler_effect_sizes.csv, euler_phenotype_association.csv; {OUT}/euler_violin.pdf")
    return 0


def make_violin():
    mri = pd.read_csv(_config.MRI_FILE, sep="\t", low_memory=False)
    mri["_jk"] = mri.apply(_m01._canonical_join_key_mri, axis=1)
    cur = pd.read_csv(CURATED, low_memory=False)
    cur["_jk"] = cur.apply(_m01._canonical_join_key_clusters, axis=1)
    pop = dict(zip(cur["_jk"], cur["population"]))
    mri["pop"] = mri["_jk"].map(pop)
    mri = mri.dropna(subset=["pop", "mean_euler"])
    mri["pop"] = mri["pop"].replace(POP_REMAP)
    order = [g for g in GROUPS if (mri["pop"] == g).any()]

    fig, (ax, axall) = plt.subplots(2, 1, figsize=(13, 6.0), height_ratios=[len(order), 1.4])
    def draw(axis, groups, data_by):
        for i, g in enumerate(groups):
            v = data_by(g)
            if len(v) == 0: continue
            axis.violinplot(v, positions=[i], vert=False, widths=0.9, showextrema=False)
            for b in axis.collections[-1:]:
                b.set_facecolor("0.85"); b.set_edgecolor("none"); b.set_alpha(0.9)
            axis.scatter(v, np.full(len(v), i) + np.random.uniform(-0.12, 0.12, len(v)),
                         s=4, c="0.25", alpha=0.5, zorder=3, edgecolors="none")
            q1, med, q3 = np.percentile(v, [25, 50, 75])
            axis.plot([q1, q3], [i, i], color="#7a0010", lw=4, zorder=4, solid_capstyle="round")
            axis.scatter([med], [i], s=90, color="#7a0010", zorder=5)
    draw(ax, order, lambda g: mri.loc[mri["pop"] == g, "mean_euler"].values)
    ax.set_yticks(range(len(order))); ax.set_yticklabels(order); ax.invert_yaxis()
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    draw(axall, ["All"], lambda g: mri["mean_euler"].values)
    axall.set_yticks([0]); axall.set_yticklabels(["All"])
    for s in ("top", "right"): axall.spines[s].set_visible(False)
    lo = np.percentile(mri["mean_euler"], 1)
    ax.set_xlim(lo, 5); axall.set_xlim(lo, 5)
    axall.set_xlabel("mean Euler number (higher = better reconstruction)")
    fig.tight_layout()
    fig.savefig(OUT / "euler_violin.pdf", bbox_inches="tight")
    fig.savefig(OUT / "euler_violin.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    raise SystemExit(main())
