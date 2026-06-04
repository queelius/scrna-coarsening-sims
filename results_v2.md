# Validation v2: Poisson-Beta DGP at realistic scRNA-seq scale

**Setup change vs v1**: count generator upgraded from negative binomial to
Poisson-Beta two-state (the canonical bursty-transcription model used by
SymSim, Splat, and the methodological literature since Kim & Marioni 2013).
Scale upgraded from 800 cells x 5 genes to 2000 cells x 50 genes drawn from
a log-normal mu distribution (range 0.09 to 13.1, median 1.95).

## Headline result: cell-total consistency is exact regardless of DGP

| `alpha_1` | observed mean total | naive predicted total | diff |
|---|---|---|---|
| 0.0 | 27.27 | 27.27 | +1.7e-04 |
| 0.3 | 46.23 | 46.23 | -2.0e-04 |
| 0.6 | 70.72 | 70.72 | -1.2e-04 |
| 1.0 | 100.05 | 100.05 | +3.1e-05 |

The naive ZINB's predicted observed mean equals the empirical observed mean
to ~1e-4 across all 50 genes summed. The misspecification result holds
exactly when the count distribution is more realistic (and bursty) than the
NB the ZINB assumes.

## Per-gene `mu` relative bias (across 30 reps, 50 genes)

Quantiles (5 / 25 / 50 / 75 / 95):

| `alpha_1` | naive median (range) | oracle median (range) |
|---|---|---|
| 0.0 (C2 holds, but ZINB misspec for PB) | 0.226 (0.107 to 0.510) | 0.023 (-0.009 to 0.064) |
| 0.3 (mild violation)                    | 0.473 (0.335 to 0.657) | 0.020 (-0.010 to 0.041) |
| 0.6 (moderate)                          | 0.639 (0.475 to 0.863) | 0.023 (-0.005 to 0.040) |
| 1.0 (severe)                            | 0.781 (0.498 to 1.069) | 0.018 ( 0.000 to 0.031) |

**Naive bias decomposes into two sources**:
1. *Intrinsic distributional misspec* (ZINB fitted to Poisson-Beta data):
   ~23% median bias at `alpha_1 = 0`. This exists even when C2 holds.
2. *Dropout-mechanism misspec* (C2 violation): the additional bias as
   `alpha_1` grows. At severe violation (`alpha_1 = 1.0`), median bias
   reaches 78% (additional 55 percentage points beyond the intrinsic).

**The dropout contribution is the identifiable one.** Spike-ins recover
oracle mu within 2-3% median bias at every violation level. They cannot
fix the intrinsic distributional misspec (that requires switching the
estimator class), but they can fix the dropout part with high precision.

## Bias by mu quintile (the differential-expression-testing implication)

Mean naive relative bias by gene-mu quintile (Q1 = lowest expression):

| `alpha_1` | Q1 (low mu) | Q2 | Q3 | Q4 | Q5 (high mu) | spread |
|---|---|---|---|---|---|---|
| 0.0 | +40% | +33% | +25% | +18% | +12% | 28 pp |
| 0.3 | +63% | +55% | +49% | +41% | +36% | 27 pp |
| 0.6 | +83% | +72% | +65% | +57% | +50% | 33 pp |
| 1.0 | +104% | +89% | +77% | +64% | +52% | 52 pp |

**Replicated and intensified from v1.** Low-expression genes are more
biased than high-expression genes, and the gradient steepens with
violation severity. Under severe violation, the lowest-quintile genes are
biased by more than 100% (more than double the true rate) while the
highest-quintile genes are biased by ~50%.

**Implication**: a naive ZINB analysis distorts low-vs-high gene
comparisons in a structured way. Differential-expression tests applied
to naive estimates would systematically over-call low-expression genes
as differentially expressed (since their naive rates are inflated more).
This is the empirical hook the methods paper can lead with.

## Oracle cell-total mu recovery

| `alpha_1` | true sum mu | oracle sum mu | rel error |
|---|---|---|---|
| 0.0 | 149.19 | 155.06 | +3.9% |
| 0.3 | 149.19 | 152.73 | +2.4% |
| 0.6 | 149.19 | 153.30 | +2.8% |
| 1.0 | 149.19 | 151.78 | +1.7% |

Oracle recovers true sum within 4% even at severe violation, with 30
replicates and 200 spike-in cells. (Slight upward bias in oracle reflects
the residual ZINB-vs-PB distributional misspec; the dropout part is
recovered cleanly.)

## What v2 adds beyond v1

1. **Robustness to count distribution**: cell-total consistency holds for
   bursty Poisson-Beta data, not just NB. The misspecification result is
   not an artifact of the v1 NB-on-NB matched DGP/estimator pair.
2. **Realistic scale**: 50 genes from log-normal mu distribution (0.09 to
   13.1) rather than 5 hand-picked genes. The bias pattern holds across
   the realistic mu range.
3. **Bias decomposition**: distributional-misspec contribution is separable
   from dropout-misspec contribution. Spike-ins fix the latter, not the
   former. This nuances the paper's claim: spike-ins restore identifiability
   of the dropout function, but the choice of estimator class still matters.
4. **Strengthened DE-testing implication**: the low-mu-bias gradient is
   stronger at scale (52pp spread at severe violation vs 32pp in v1). The
   empirical claim about false-positive DE on low-expression genes is more
   defensible.
5. **Computational feasibility**: 30 reps x 4 alpha levels at 2000x50 scale
   ran in ~80 seconds total. Scaling to a real-data benchmark on tabula
   muris or similar is well within reach.

## What v2 does NOT yet add

- Different dropout functional form (only logistic-in-log-y tested; the
  alternative exponential-decay form is implemented in `apply_dropout`
  but not yet swept).
- Comparison against published scRNA-seq simulators (SymSim, Splat).
- Real-data application: the test on a public scRNA-seq dataset where
  spike-ins are actually present.

These are the remaining items if Jingyi Li's collaboration wants to
push toward a Genome Biology submission with empirical bench. For a
methods-only submission to Biostatistics or AOAS, v2 evidence is
sufficient.

## Files

- `sim.R`: now contains both v1 (`sim_scrna`) and v2 (`sim_scrna_pb`,
  `sim_pb_gene`, `realistic_mu`, `apply_dropout`) generators.
- `run_pb.R`: v2 driver
- `results_pb.rds`: full v2 per-replicate results
