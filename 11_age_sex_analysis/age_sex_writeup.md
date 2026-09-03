# Effects of Age and Sex on Anatomical and Functional Group-Level Results

## Methods

### Samples
Age and sex effects were tested in the autistic versus neurotypical (NT)
contrast for both imaging modalities. The anatomical sample comprised 1262
participants (897 autism, 365 NT) with quality-controlled, ComBat-harmonised
FreeSurfer estimates; the functional sample comprised 720 participants (409
autism, 311 NT) with quality-controlled, ComBat-harmonised functional
connectivity (4S156Parcels atlas). Age ranged from 1.7–55.0 years
(anatomical) and 2.3–54.4 years (functional).

### Anatomical analysis
For each cortical region (thickness, surface area, cortical volume; Desikan–
Killiany atlas) and subcortical volume (ASEG atlas), we fitted two ordinary
least-squares models on the autism+NT sample:

- *age model*: `region ~ diagnosis * age + sex + cohort`
- *sex model*: `region ~ diagnosis * sex + age + cohort`

with age mean-centred and cohort (INOVAND/LEAP) included to absorb residual
site variance. The single interaction coefficient (diagnosis × age or
diagnosis × sex) was retained per region. p-values were corrected for
multiple comparisons using the Benjamini–Hochberg false discovery rate (FDR)
**within each metric** (232 regions total: 68 thickness, 68 surface area, 68
cortical volume, 28 subcortical volume), consistent with the main
autism-versus-NT anatomical analysis. A region was considered to show a
robust interaction if it survived FDR < 0.05 within its metric. Standardised
effect sizes were obtained from the coefficient t-statistic
(d = 2t / √df_resid).

### Functional analysis
Because edge-level mass-univariate testing across 10,731 connections is
strongly underpowered for an interaction effect, the functional interactions
were assessed with the Network-Based Statistic (NBS; Zalesky et al., 2010),
which provides component-level, family-wise-error (FWER)-controlled
inference exploiting network topology. For each edge we computed the
diagnosis × age (and diagnosis × sex) interaction t-statistic from the model
`edge ~ diagnosis + age + sex + cohort + diagnosis:modulator`. By the
Frisch–Waugh–Lovell theorem this equals the t-statistic of the edge
residualised on the reduced (no-interaction) design regressed on the
similarly residualised interaction regressor, allowing all edges to be
evaluated simultaneously. Supra-threshold edges (primary edge-forming
threshold p < 0.01) were grouped into connected components in the 156-node
graph, separately for the positive (autism > NT) and negative (autism < NT)
directions, with component size defined as the number of edges. The null
distribution of the largest component size was built from 5000 permutations
of the orthogonalised interaction regressor (Smith et al., 2007), which is
exchangeable under the null of no partial association. The FWER p-value of
each observed component was the proportion of permutations whose largest
component equalled or exceeded it. Analyses used Python 3.11 (statsmodels,
scipy) and R 4.x (ggseg) for visualisation.

## Results

### Age by diagnosis
At the anatomical level, three cortical regions showed a diagnosis × age
interaction surviving FDR correction within their metric: cortical thickness
of the right (t = −3.69, p_FDR = 0.016, d = −0.21) and left (t = −3.38,
p_FDR = 0.025, d = −0.19) rostral middle frontal gyrus, and surface area of
the left medial orbitofrontal cortex (t = −3.44, p_FDR = 0.041, d = −0.19).
In each case the negative interaction reflects a shallower age-related
increase (autism flat versus a positive NT slope) in the harmonised metric.
No subcortical volume showed a significant interaction (minimum p_FDR = 0.84),
and effect sizes were small throughout (|d| ≤ 0.21).

However, these interactions were not robust to the marked age imbalance
between groups. Below age 10 the sample contained 420 autistic but only 51 NT
participants (under age 6: 223 vs 6), so the NT age slope is poorly supported
at the young end. In a sensitivity analysis restricting to age ≥ 10 — where
both groups are adequately represented — the interaction disappeared in both
regions (left: t = −0.16, p = 0.88; right: t = −1.39, p = 0.17), and it was
already abolished or attenuated below significance when excluding only the
under-6s (age ≥ 6 left: t = −1.24, p = 0.22; right: t = −2.09, p = 0.037).
The full-sample effect is therefore driven by extrapolation of the NT
trajectory through an age range with very few NT observations and should not
be interpreted as a robust developmental difference.

At the functional level, the NBS revealed a candidate network component in
the autism > NT direction at the primary threshold of p = 0.01, but it did
not reach significance (largest component = 31 edges, p_FWER = 0.52); the
autism < NT component was likewise non-significant (70 edges, p_FWER = 0.22).
No individual edge survived FDR correction (minimum FDR-corrected p = 0.59;
maximum |d| = 0.30). The double network matrices of these components (Figure,
functional) show connections distributed across visual, control and attention
networks and, proportionally, between the thalamus and somatomotor, attention
and limbic systems — but in this sample the functional age-by-diagnosis effect
remained below the significance threshold.

### Sex by diagnosis
No anatomical region showed a diagnosis × sex interaction surviving FDR
correction in any metric (minimum p_FDR = 0.14, subcortical volume; 9/232
regions at uncorrected p < 0.05; maximum |d| = 0.19). At the functional
level, the largest autism > NT component (146 edges) reached an uncorrected
component-level p_FWER = 0.019, but this does not survive correction for the
family of interaction tests performed (≈ 6 comparisons; corrected p ≈ 0.11);
the autism < NT component was non-significant (41 edges, p_FWER = 0.44). No
edge survived FDR correction (minimum FDR-corrected p = 0.29; maximum
|d| = 0.32).

## Discussion

We found marginal evidence of an age-by-diagnosis interaction in prefrontal
cortical structure (rostral middle frontal thickness bilaterally, left medial
orbitofrontal surface area) that survived within-metric FDR correction in the
full sample. This effect did not, however, survive our robustness checks: it
was abolished when the analysis was restricted to the age range in which both
groups are adequately represented (age ≥ 10), and was already attenuated below
significance when only the under-6s were excluded. This sensitivity reflects
the strong age imbalance between groups — the NT group is sparse below age 10
(6 NT versus 223 autistic participants under age 6) — so that the full-sample
interaction is driven largely by extrapolation of the NT trajectory through a
poorly supported age range rather than by a genuine difference in maturational
slope. We therefore did not replicate previously reported developmental effects
in autism, and conclude that there is no robust general developmental trajectory
distinguishing autistic from NT development after early childhood. This may
reflect heterogeneous developmental trajectories within the autism group. At
the functional level the age-by-diagnosis interaction likewise did not reach
component-level significance and no edge survived FDR correction, consistent
across modalities with the absence of a robust developmental effect. Any
sub-threshold trends — the prefrontal anatomical interaction and the functional
age-by-diagnosis component — should be interpreted with caution.

We did not observe a significant sex-by-diagnosis interaction in either
modality, suggesting limited evidence of sexually distinct patterns of altered
structure or functional connectivity previously reported in autism. Adequate
recruitment of women and girls with autism remains a challenge for the field —
the present sample is approximately 2.5:1 male-to-female — and these results
require replication in samples with a more balanced male-to-female ratio.
