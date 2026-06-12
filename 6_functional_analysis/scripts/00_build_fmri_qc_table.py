#!/usr/bin/env python3
"""
Build a per-subject fMRI QC table from XCP-D outputs.

For every subject row in df_conn_cohort_norm.csv we recover the matching
*_desc-linc_qc.tsv (XCP-D) and *_desc-denoised_bold.json (for TR) and compute:

    minutes_quality_data = num_retained_volumes * TR / 60

Outputs:
    outputs/tables/fmri_qc_per_subject.csv

This table feeds the two reviewer-requested sensitivity analyses:
  - 6-min quality-data cutoff (Autism-vs-NT and cluster-vs-NT)
  - mean_fd + minutes_quality_data as covariates
"""

import json
import re
import sys
from pathlib import Path

import pandas as pd

_SCRIPT_DIR = Path(__file__).parent
_SECTION_DIR = _SCRIPT_DIR.parent
OUTPUT_DIR = _SECTION_DIR / 'outputs' / 'tables'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

EEG_MRI_RESULTS = Path('/Users/mfleury/POSTDOC/LIBRAIRY/eeg_mri-pipeline/results')
FMRI_CONN_FILE = _SECTION_DIR / 'preprocessing' / 'outputs' / 'df_conn_cohort_norm.csv'

LINC_QC_COLS = [
    'mean_fd',
    'mean_fd_post_censoring',
    'num_dummy_volumes',
    'num_censored_volumes',
    'num_retained_volumes',
]


def derive_linc_qc_path(conn_path: str) -> Path:
    """The connectivity .tsv and the linc_qc .tsv share the same prefix."""
    return Path(
        conn_path.replace(
            '_seg-4S156Parcels_stat-pearsoncorrelation_relmat.tsv',
            '_desc-linc_qc.tsv',
        )
    )


def derive_bold_json_path(conn_path: str) -> Path:
    return Path(
        conn_path.replace(
            '_seg-4S156Parcels_stat-pearsoncorrelation_relmat.tsv',
            '_res-2_desc-denoised_bold.json',
        )
    )


_FALLBACK_RE = re.compile(r'_run-\d+_space-')


def find_fallback_linc_qc(conn_path: str) -> Path | None:
    """
    When the run referenced in conn_path no longer exists on disk (e.g. the
    primary connectivity was generated from a previous xcp_d run that has
    since been re-numbered), pick any linc_qc.tsv from the same func/ folder.
    Use the alphabetically-first match so the choice is deterministic.
    """
    parent = derive_linc_qc_path(conn_path).parent
    if not parent.is_dir():
        return None
    candidates = sorted(parent.glob('*_desc-linc_qc.tsv'))
    return candidates[0] if candidates else None


def find_fallback_bold_json(conn_path: str) -> Path | None:
    parent = derive_bold_json_path(conn_path).parent
    if not parent.is_dir():
        return None
    candidates = sorted(parent.glob('*_res-2_desc-denoised_bold.json'))
    return candidates[0] if candidates else None


def parse_linc_qc(qc_path: Path) -> dict:
    if not qc_path.is_file():
        return {c: None for c in LINC_QC_COLS}
    df = pd.read_csv(qc_path, sep='\t')
    if len(df) == 0:
        return {c: None for c in LINC_QC_COLS}
    row = df.iloc[0]
    return {c: row[c] if c in df.columns else None for c in LINC_QC_COLS}


def parse_repetition_time(json_path: Path) -> float | None:
    if not json_path.is_file():
        return None
    try:
        with open(json_path) as f:
            meta = json.load(f)
        return float(meta.get('RepetitionTime')) if meta.get('RepetitionTime') is not None else None
    except Exception:
        return None


def main() -> int:
    if not FMRI_CONN_FILE.exists():
        raise FileNotFoundError(f'Connectivity dataframe not found: {FMRI_CONN_FILE}')

    df = pd.read_csv(FMRI_CONN_FILE, low_memory=False, usecols=['ID', 'path', 'session', 'run', 'cohort'])
    df['ID'] = df['ID'].astype(str)
    print(f'Loaded {len(df)} subject rows from df_conn_cohort_norm.csv')

    records = []
    n_qc_exact = n_qc_fallback = n_qc_missing = 0
    n_tr_exact = n_tr_fallback = n_tr_missing = 0

    for _, row in df.iterrows():
        conn_path = row['path']

        # --- linc_qc.tsv ------------------------------------------------------
        qc_path = derive_linc_qc_path(conn_path)
        qc_source = 'exact'
        if not qc_path.is_file():
            fb = find_fallback_linc_qc(conn_path)
            if fb is not None:
                qc_path = fb
                qc_source = 'fallback_same_session'
            else:
                qc_source = 'missing'

        qc_values = parse_linc_qc(qc_path) if qc_source != 'missing' else {c: None for c in LINC_QC_COLS}

        if qc_source == 'exact':
            n_qc_exact += 1
        elif qc_source == 'fallback_same_session':
            n_qc_fallback += 1
        else:
            n_qc_missing += 1

        # --- RepetitionTime ---------------------------------------------------
        json_path = derive_bold_json_path(conn_path)
        tr_source = 'exact'
        tr_sec = parse_repetition_time(json_path)
        if tr_sec is None:
            fb = find_fallback_bold_json(conn_path)
            if fb is not None:
                tr_sec = parse_repetition_time(fb)
                tr_source = 'fallback_same_session'
            else:
                tr_source = 'missing'

        if tr_source == 'exact' and tr_sec is not None:
            n_tr_exact += 1
        elif tr_source == 'fallback_same_session' and tr_sec is not None:
            n_tr_fallback += 1
        else:
            n_tr_missing += 1

        # --- compute minutes of quality data ----------------------------------
        n_retained = qc_values.get('num_retained_volumes')
        minutes_quality = None
        if n_retained is not None and tr_sec is not None:
            try:
                minutes_quality = float(n_retained) * float(tr_sec) / 60.0
            except (TypeError, ValueError):
                minutes_quality = None

        records.append({
            'ID': row['ID'],
            'cohort': row['cohort'],
            'session': row['session'],
            'run': row['run'],
            'path_conn': conn_path,
            'path_linc_qc': str(qc_path) if qc_source != 'missing' else '',
            'linc_qc_source': qc_source,
            'TR_sec': tr_sec,
            'TR_source': tr_source,
            'num_dummy_volumes': qc_values.get('num_dummy_volumes'),
            'num_censored_volumes': qc_values.get('num_censored_volumes'),
            'num_retained_volumes': qc_values.get('num_retained_volumes'),
            'mean_fd': qc_values.get('mean_fd'),
            'mean_fd_post_censoring': qc_values.get('mean_fd_post_censoring'),
            'minutes_quality_data': minutes_quality,
        })

    out_df = pd.DataFrame.from_records(records)

    print('\nProvenance summary:')
    print(f'  linc_qc.tsv  : {n_qc_exact} exact / {n_qc_fallback} fallback / {n_qc_missing} missing')
    print(f'  bold.json TR : {n_tr_exact} exact / {n_tr_fallback} fallback / {n_tr_missing} missing')

    n_minutes = out_df['minutes_quality_data'].notna().sum()
    print(f'\nminutes_quality_data computed for {n_minutes}/{len(out_df)} subjects')
    if n_minutes:
        s = out_df['minutes_quality_data'].dropna()
        print(f'  min / median / max  : {s.min():.2f} / {s.median():.2f} / {s.max():.2f} min')
        print(f'  fraction <  3 min   : {(s < 3).mean():.1%}')
        print(f'  fraction <  6 min   : {(s < 6).mean():.1%}')
        print(f'  fraction >= 6 min   : {(s >= 6).mean():.1%}')

    n_fd = out_df['mean_fd'].notna().sum()
    if n_fd:
        s_fd = out_df['mean_fd'].dropna()
        print(f'\nmean_fd computed for {n_fd}/{len(out_df)} subjects')
        print(f'  fraction <  1 mm    : {(s_fd < 1).mean():.1%}')
        print(f'  fraction <  0.5 mm  : {(s_fd < 0.5).mean():.1%}')
        print(f'  fraction <  0.3 mm  : {(s_fd < 0.3).mean():.1%}')

    out_path = OUTPUT_DIR / 'fmri_qc_per_subject.csv'
    out_df.to_csv(out_path, index=False)
    print(f'\nSaved QC table: {out_path}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
