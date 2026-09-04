# Functional connectivity methods — section 6 vs 7 (and 6 vs "6-block")

> **Note (2026-09):** this note was written against the `non_concat/` (legacy
> XCP-D 0.8) connectivity table below, which was the shared input at the time.
> The manuscript-reported whole-group and per-cluster analyses now read
> `concat/preprocessing/outputs/nogsr_concat/df_conn_cohort_norm_nogsr_concat.csv`
> instead (see `README.md`) — the *inference logic* described here (edge /
> network-block / NBS, whole-group vs per-cluster) still applies unchanged,
> only the input table and cohort size differ from what's written below.

All functional-connectivity analyses share **one input and one control group**;
they differ only in (a) *who is on the "case" side* of the contrast and (b) *what
inference is applied*. This note pins down those differences.

## Shared foundation (as originally run, non_concat pipeline)

- **Connectivity**: `6_functional_analysis/non_concat/preprocessing/outputs/df_conn_cohort_norm.csv`
  — xcp_d / 4S156Parcels, **156 nodes → 10,731 edges** per subject, ComBat-
  harmonised, curated cohort.
- **Diagnosis**: `population_group` from the curated clinical merge (TD → NADT relabelled **NT**). Cluster labels (k-means) attached from `df_clusters_complete_kmeans.csv`.
- **Control group (NT)**: the **same** pooled NT set (~296) is the reference for
  *every* contrast below.
- **Per-edge statistic**: **Welch t** (unequal variance), autism-side − NT.
- **Direction split**: **hypo** = case < NT (t<0); **hyper** = case > NT (t>0).

So nothing about the data or the per-edge test changes between 6, "6-block", and 7.
Two things change: the **case group** and the **inference/thresholding**.

---

## Section 6 — whole-group Autism vs NT

**Case side = ALL autistic participants pooled** (~330) vs NT (~296). One contrast.
Question: *does autism as a whole differ from NT in connectivity?*

Section 6 applies **three levels of inference** to that single contrast
(`run_functional_analysis.py` runs them in order):

1. **Edge-level** (`01_…`) — Welch t per edge, **BH-FDR q<0.05**.
2. **Network-block** (`06_network_block_permutation.py`) — this is the
   **"6-block"** analysis: average each subject's edges within every network-pair
   block (11 networks → 66 blocks), Welch t per block, **label-permutation p**
   (5000 perms) + BH-FDR across the 66 blocks. Motivation: edge-level FDR is
   underpowered for distributed effects and ignores network dependence; the
   block test gives each **network pair** a real p-value.
3. **Component-level** (`07_nbs_test.py`, `08_…`) — **NBS** (Zalesky 2010):
   primary edge threshold p<0.01, connected components, **FWER by permutation**.

Plus two **sensitivity re-runs** of the edge-level analysis (same contrast):
- `01b` — restrict to ≥6 min post-censoring data & mean FD < 1 mm.
- `01c` — OLS `connectivity ~ group + mean_fd + minutes_quality_data`
  (motion + data-quantity covariates on the group term).

**In short:** "6" = the edge-level autism-vs-NT map; **"6-block"** = the same
contrast tested at the **network-block** level (66 network pairs, permutation).
Same data, same groups — a *coarser, dependence-aware* unit of inference.

---

## Section 7 — Autism *cluster* vs NT (stratified)

**Case side = each IQ×SRS cluster separately**: C1 (n≈60), C2 (≈110), C3 (≈78),
each vs the **same** NT (~296). Three contrasts instead of one.
Question: *does each clinical subgroup differ from NT, and do the subgroups
differ from each other?*

Section 7 mirrors section 6's machinery per cluster:
- edge-level Welch t + FDR (`01_…`),
- network-matrix summary (`04_…`, the cluster analogue of 6's block view),
- **NBS per cluster** (`06_cluster_nbs.py`) — this feeds **figure_6b**,
- prominent networks/regions digest (`09_…`).

`figure_6b` = the per-cluster NBS/edge maps (cortical Schaefer surface +
subcortical insets), one column per cluster, hypo/hyper blocks.

---

## One-line distinctions

| | case group | inference unit | # contrasts |
|---|---|---|---|
| **6** (edge) | all autism | single edge (FDR) | 1 |
| **6-block** (network-block) | all autism | network **pair** (11×11 → 66 blocks, permutation) | 1 |
| **7** (cluster) | each cluster C1/C2/C3 | edge / component, **per cluster** | 3 |

- **6 → 6-block**: same autism-vs-NT contrast, moved from *edge* to *network-block*
  inference (more power for distributed effects, respects network dependence).
- **6 → 7**: same edge-level machinery, but the autism side is **split by
  cluster** — turning one whole-group map into three subgroup maps (figure_6b).

## Thresholding note (curated cohort)

On the curated data, **BH-FDR wipes out almost everything** (per-cluster: 0/0/1
edges). So figure_6b is shown at either **NBS** (largest component, primary
p<0.01; no component reaches FWER<0.05) or **p<0.05 uncorrected** — stated in the
caption. The network-block ("6-block") test is the better-powered whole-group
counterpart when edge-FDR is empty.
