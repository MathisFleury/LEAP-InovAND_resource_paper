#!/usr/bin/env python3
"""
Motion (mean FD) confound check for the fMRI NBS group differences (6-2).

Reviewer ask: report effect sizes of the association between image-quality
(mean FD) and the connectivity phenotypes that show group differences, to
disambiguate whether quality could be driving the Autism-vs-NT effects.

Produces, for the featured variant (default nogsr_concat):
  (1) FD group difference (Autism vs NT): Welch t, Cohen's d.
  (2) FD <-> subnetwork-connectivity association for the significant NBS
      hyper / hypo components: Pearson r (+95% CI) computed pooled (group
      partialled out), within-NT, and within-Autism. The sign is compared with
      the group-effect direction (hyper = Autism>NT positive; hypo = Autism<NT
      negative). If FD tracks connectivity OPPOSITE to the group effect, motion
      cannot be driving it; if SAME direction, motion is a candidate confound.
  (3) FD raincloud/violin by population group (+ an "All" row), mirroring
      fd_violinplot.png.

Outputs -> outputs/qc_motion/:
  fd_effect_sizes.csv, fd_connectivity_association.csv, fd_violin.{pdf,png}
"""
import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind, pearsonr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _sensitivity_utils import filter_connectivity_cols  # noqa: E402

_SECTION = Path(__file__).resolve().parent.parent
VARIANT = os.environ.get("XCPD_VARIANT_TAG", "nogsr_concat")
CONN = _SECTION / "preprocessing" / "outputs" / VARIANT / f"df_conn_cohort_norm_{VARIANT}.csv"
NBS_DIR = _SECTION / "outputs" / VARIANT / "nbs_double"
CURATED = Path("/Users/mfleury/POSTDOC/LIBRAIRY/LEAP-InovAND_resource/"
               "1_clustering/outputs/curated/tables/individuals_metrics_with_clusters_curated.csv")
OUT = _SECTION / "outputs" / "qc_motion"
OUT.mkdir(parents=True, exist_ok=True)

# Standard population naming. "Autism to exclude" (and bare "Autism", subtype
# unspecified) are folded into "Autism without IDD".
POP_REMAP = {"Autism to exclude": "Autism without IDD", "Autism": "Autism without IDD"}
GROUPS = [("NT", "NT"), ("Relatives", "Relatives"), ("IDD", "IDD"),
          ("Autism without IDD", "Autism without IDD"), ("Autism with IDD", "Autism with IDD")]


def cohens_d(a, b):
    na, nb = len(a), len(b)
    sd = np.sqrt(((na - 1) * np.var(a, ddof=1) + (nb - 1) * np.var(b, ddof=1)) / (na + nb - 2))
    return (np.mean(a) - np.mean(b)) / sd


def r_ci(r, n):
    if n < 4:
        return (np.nan, np.nan)
    z = np.arctanh(r); se = 1 / np.sqrt(n - 3)
    return tuple(np.tanh([z - 1.96 * se, z + 1.96 * se]))


def subnetwork_mean(df, conn_cols, edge_file):
    """Per-subject mean connectivity over the NBS-component edges in edge_file."""
    if not edge_file.exists():
        return None
    e = pd.read_csv(edge_file)
    cols = [f"con_{s}/{t}" for s, t in zip(e["source"], e["target"])]
    cols = [c for c in cols if c in conn_cols]
    if not cols:
        return None
    return df[cols].mean(axis=1)


def main() -> int:
    df = pd.read_csv(CONN, low_memory=False)
    df["ID"] = df["ID"].astype(str)
    conn_cols = filter_connectivity_cols(df.columns.tolist())
    df["dx"] = df["population_group"].replace("TD", "NT")

    # finer population label by ID (for the violin)
    cur = pd.read_csv(CURATED, low_memory=False)
    cur["ID"] = cur["ID"].astype(str)
    fine = dict(zip(cur["ID"], cur["population"]))
    df["pop_fine"] = df["ID"].map(fine).fillna(df["population_group"]).replace(POP_REMAP)

    # ---------- (1) FD group difference ----------
    fd_a = df.loc[df.dx == "Autism", "mean_fd"].dropna().values
    fd_n = df.loc[df.dx == "NT", "mean_fd"].dropna().values
    t, p = ttest_ind(fd_a, fd_n, equal_var=False)
    d = cohens_d(fd_a, fd_n)
    rows = [{"measure": "mean_FD", "contrast": "Autism vs NT",
             "autism_mean": np.mean(fd_a), "nt_mean": np.mean(fd_n),
             "welch_t": t, "p": p, "cohens_d": d, "n_autism": len(fd_a), "n_nt": len(fd_n)}]
    pd.DataFrame(rows).to_csv(OUT / "fd_effect_sizes.csv", index=False)
    print(f"[FD group diff] Autism {np.mean(fd_a):.3f} vs NT {np.mean(fd_n):.3f} mm | "
          f"Welch t={t:.2f} p={p:.2e} Cohen's d={d:+.3f}")

    # ---------- (2) FD <-> subnetwork connectivity ----------
    an = df[df.dx.isin(["Autism", "NT"])].copy()
    assoc = []
    for direction, grp_sign in [("hyper", +1), ("hypo", -1)]:
        ef = NBS_DIR / f"autism_vs_td_{direction}connectivity_full_data.csv"
        sm = subnetwork_mean(an, conn_cols, ef)
        if sm is None:
            print(f"  {direction}: no edges/file — skipping"); continue
        an[f"sn_{direction}"] = sm
        fd = an["mean_fd"]
        for scope, mask in [("pooled", an.dx.isin(["Autism", "NT"])),
                            ("within_NT", an.dx == "NT"), ("within_Autism", an.dx == "Autism")]:
            x = fd[mask].values; y = an.loc[mask, f"sn_{direction}"].values
            ok = np.isfinite(x) & np.isfinite(y)
            r, pr = pearsonr(x[ok], y[ok]); lo, hi = r_ci(r, ok.sum())
            group_dir = "Autism>NT (+)" if grp_sign > 0 else "Autism<NT (-)"
            # same-direction confound if FD-conn r has same sign as the group effect
            verdict = ("SAME direction as group effect -> candidate confound"
                       if np.sign(r) == grp_sign else
                       "OPPOSITE to group effect -> motion not driving it")
            assoc.append({"subnetwork": direction, "group_direction": group_dir, "scope": scope,
                          "r_FD_vs_conn": round(r, 3), "ci_lo": round(lo, 3), "ci_hi": round(hi, 3),
                          "p": pr, "n": int(ok.sum()), "interpretation": verdict})
            print(f"  {direction:5s} [{scope:13s}] r(FD,conn)={r:+.3f} [{lo:+.2f},{hi:+.2f}] "
                  f"p={pr:.3f}  | group={group_dir} -> {verdict}")
    pd.DataFrame(assoc).to_csv(OUT / "fd_connectivity_association.csv", index=False)

    # ---------- (3) FD violin (reference style) ----------
    make_violin(df)
    print(f"\nSaved: {OUT}/fd_effect_sizes.csv, fd_connectivity_association.csv, fd_violin.pdf")
    return 0


def make_violin(df):
    order = [g for g, _ in GROUPS if (df["pop_fine"] == g).any()]
    labels = dict(GROUPS)
    fig, (ax, axall) = plt.subplots(2, 1, figsize=(13, 6.5), height_ratios=[len(order), 1.4],
                                    sharex=True)
    def draw(axis, groups, data_by):
        for i, g in enumerate(groups):
            v = data_by(g)
            if len(v) == 0:
                continue
            axis.violinplot(v, positions=[i], vert=False, widths=0.9, showextrema=False)
            for b in axis.collections[-1:]:
                b.set_facecolor("0.85"); b.set_edgecolor("none"); b.set_alpha(0.9)
            axis.scatter(v, np.full(len(v), i) + np.random.uniform(-0.12, 0.12, len(v)),
                         s=4, c="0.25", alpha=0.5, zorder=3, edgecolors="none")
            q1, med, q3 = np.percentile(v, [25, 50, 75])
            axis.plot([q1, q3], [i, i], color="#7a0010", lw=4, zorder=4, solid_capstyle="round")
            axis.scatter([med], [i], s=90, color="#7a0010", zorder=5)
    draw(ax, order, lambda g: df.loc[df.pop_fine == g, "mean_fd"].dropna().values)
    ax.set_yticks(range(len(order))); ax.set_yticklabels([labels[g] for g in order])
    ax.invert_yaxis(); ax.set_xlim(0, 1.0)
    for s in ("top", "right"): ax.spines[s].set_visible(False)
    draw(axall, ["All"], lambda g: df["mean_fd"].dropna().values)
    axall.set_yticks([0]); axall.set_yticklabels(["All"]); axall.set_xlim(0, 1.0)
    for s in ("top", "right"): axall.spines[s].set_visible(False)
    axall.set_xlabel("framewise displacement — FD (mm)")
    fig.tight_layout()
    fig.savefig(OUT / "fd_violin.pdf", bbox_inches="tight")
    fig.savefig(OUT / "fd_violin.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    raise SystemExit(main())
