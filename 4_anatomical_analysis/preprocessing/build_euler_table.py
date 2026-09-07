"""
One-time extraction: walk each NDD cohort's recon-all logs and pull the
orig.nofix Euler numbers (lheno, rheno).

Output -> OUT_DIR/euler_numbers_recon_all.tsv  with columns:
    cohort_family   MRI_ID_n   ses_id   lh_euler_nofix   rh_euler_nofix   mean_euler_nofix

Run once before the pipeline so step_02b can merge it in (already-present
output is reused — no need to re-run unless new scans were added).

    /usr/local/bin/python3.11 build_euler_table.py
"""
import re
import pandas as pd
from pathlib import Path

from config import OUT_DIR, IMG5_ROOT

COHORT_ROOTS = {
    "INOVAND":   IMG5_ROOT / "INOVAND" / "_nipoppy" / "derivatives" / "freesurfer" / "8.1.0" / "output",
    "LEAP":      IMG5_ROOT / "LEAP" / "_nipoppy" / "derivatives" / "freesurfer" / "8.1.0" / "output",
    "INFOR":     IMG5_ROOT / "INFOR" / "_nipoppy" / "derivatives" / "freesurfer" / "8.1.0" / "output",
    "PIP":       IMG5_ROOT / "PIP" / "_nipoppy" / "derivatives" / "freesurfer" / "8.1.0" / "output_new",
    "CANDY_PIP": IMG5_ROOT / "CANDY_PIP" / "_nipoppy" / "derivatives" / "freesurfer" / "8.1.0" / "output_new",
}

EULER_PAT = re.compile(r"orig\.nofix\s+lheno\s*=\s*(-?\d+)\s*,\s*rheno\s*=\s*(-?\d+)")
SUB_PAT   = re.compile(r"sub-([^_]+)_ses-(\d+)")


def _strip_zeros(s: str) -> str:
    """0001 -> 1 ;  100693509718 -> 100693509718 ;  0 -> 0"""
    return s.lstrip("0") or s


def _parse_log(log_path: Path) -> tuple[int, int] | None:
    """Return (lheno, rheno) from a recon-all.log, or None if not present."""
    try:
        with open(log_path, "r", errors="ignore") as f:
            for line in f:
                m = EULER_PAT.search(line)
                if m:
                    return int(m.group(1)), int(m.group(2))
    except OSError:
        return None
    return None


def _walk_cohort(cohort_family: str, root: Path) -> list[dict]:
    rows = []
    if not root.exists():
        print(f"[euler] {cohort_family}: root does not exist  ({root})")
        return rows
    subj_dirs = [d for d in root.iterdir() if d.is_dir()]
    print(f"[euler] {cohort_family}: walking {len(subj_dirs):,} subject_ses dirs")
    n_ok = n_skip = 0
    for d in subj_dirs:
        m = SUB_PAT.match(d.name)
        if not m:
            n_skip += 1; continue
        log = d / "scripts" / "recon-all.log"
        if not log.exists():
            n_skip += 1; continue
        e = _parse_log(log)
        if e is None:
            n_skip += 1; continue
        lh, rh = e
        rows.append({
            "cohort_family":    cohort_family,
            "MRI_ID_n":         _strip_zeros(m.group(1)),
            "ses_id":           _strip_zeros(m.group(2)),
            "lh_euler_nofix":   lh,
            "rh_euler_nofix":   rh,
            "mean_euler_nofix": (lh + rh) / 2,
        })
        n_ok += 1
    print(f"[euler]   parsed {n_ok:,}   skipped {n_skip:,}")
    return rows


def main() -> None:
    all_rows = []
    for fam, root in COHORT_ROOTS.items():
        all_rows.extend(_walk_cohort(fam, root))
    df = pd.DataFrame(all_rows)
    out = OUT_DIR / "euler_numbers_recon_all.tsv"
    df.to_csv(out, sep="\t", index=False)
    print(f"\n[euler] wrote {out}  ({len(df):,} rows)")
    print(df.groupby("cohort_family")["mean_euler_nofix"].describe().round(1))


if __name__ == "__main__":
    main()
