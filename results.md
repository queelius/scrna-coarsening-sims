# Validation: scRNA-seq zero inflation as a masked-cause coarsening problem

**Date**: 2026-05-07
**Hypothesis tested**: a naive zero-inflated negative-binomial (ZINB) estimator
fit to scRNA-seq counts under C2 violation (dropout depends on the latent count
Y) converges to a pseudo-true value where:
1. per-gene mean expression `mu_j` is biased,
2. cell-total expression `sum_j (1 - pi_j) * mu_j` is unbiased for the observed
   mean,
3. spike-in controls (true Y known on a subset of cells) restore identifiability.

This is the scRNA-seq analogue of the misspecification result in `mdrelax`
for masked-cause series-system likelihood under Relaxed C3.

## Setup

- 5 genes with true means `mu = (0.5, 1, 2, 5, 10)`, NB dispersion `phi = 1`
- 800 cells, 80 replicates per condition
- Dropout function: `pi(y) = plogis(2 - alpha_1 * log(y + 1))`
- Severity grid: `alpha_1 ∈ {0, 0.3, 0.6, 1.0}` (0 = C2 holds)
- Oracle uses 10% of cells (n=80) with known Y (spike-in analogue) to estimate
  `(alpha_0, alpha_1)` via logistic regression on dropped indicator, then
  refits per-gene mu with the dropout function fixed.

## Per-gene `mu` bias (mean across 80 reps)

| `alpha_1` | true mu | naive bias | naive % bias | oracle bias |
|---|---|---|---|---|
| 0   (C2 holds)        | 0.5 / 1 / 2 / 5 / 10 | +0.06 / -0.01 / -0.07 / +0.11 / -0.05 | ~0% across the board | similar |
| 0.3 (mild violation)  | 0.5 / 1 / 2 / 5 / 10 | +0.24 / +0.17 / +0.41 / +1.07 / +2.09 | +49% / +17% / +20% / +21% / +21% | ≤ 0.04 |
| 0.6 (moderate)        | 0.5 / 1 / 2 / 5 / 10 | +0.24 / +0.41 / +0.76 / +1.95 / +3.54 | +47% / +41% / +38% / +39% / +35% | ≤ 0.11 |
| 1.0 (severe)          | 0.5 / 1 / 2 / 5 / 10 | +0.36 / +0.64 / +1.18 / +2.44 / +3.88 | +71% / +64% / +59% / +49% / +39% | ≤ 0.13 |

**Observation 1**: naive mu is consistently biased upward, with absolute bias
scaling roughly linearly with `alpha_1`. The naive estimator interprets the
dropout-induced upward shift of the surviving non-zeros as evidence that the
true expression rate is higher than it is.

**Observation 2**: the oracle, using only 10% of cells as spike-ins to fit the
dropout function, recovers per-gene mu to within Monte Carlo error at all
violation levels. This validates the singleton-candidate-set / spike-in
identification argument.

## Cell-total expression: naive predicted vs. observed

| `alpha_1` | observed mean total | naive predicted total | difference |
|---|---|---|---|
| 0   | 2.209  | 2.209  | < 1e-4 |
| 0.3 | 4.105  | 4.105  | < 1e-4 |
| 0.6 | 6.924  | 6.924  | < 1e-4 |
| 1.0 | 11.069 | 11.069 | < 1e-4 |

**Observation 3**: cell-total naive prediction matches observed cell-total mean
to numerical precision at every violation level. This is the **cell-total
consistency** that mdrelax's misspecification theorem predicts: the naive
estimator absorbs the dropout asymmetry into per-component rates while still
producing the right marginal observed-data prediction.

## Cell-total true `mu` recovery (oracle)

| `alpha_1` | true sum `mu` | oracle sum `mu` | error |
|---|---|---|---|
| 0   | 18.5 | 18.64 | +0.14 |
| 0.3 | 18.5 | 18.41 | -0.09 |
| 0.6 | 18.5 | 18.35 | -0.15 |
| 1.0 | 18.5 | 18.68 | +0.18 |

**Observation 4**: oracle recovers true cell-total `mu` within ±1% at all
violation levels. (Naive cannot, because the naive total prediction tracks the
observed mean, not the true mu sum. The gap between observed and true is the
dropped expression mass.)

## Subfinding: relative-bias pattern under stronger violation

Under `alpha_1 = 1.0` (severe), the relative bias decreases with true
expression: 71% for `mu = 0.5`, falling monotonically to 39% for `mu = 10`.
Under `alpha_1 = 0.6` the pattern is flatter (47% to 35%). Under `alpha_1 =
0.3` it inverts at the low end (49% for `mu = 0.5`, 17% for `mu = 1`, then
~20% above).

**Implication for the paper**: lowly-expressed genes are differentially mis-
categorized under the naive ZINB, which has direct consequences for
differential-expression testing. A naive analysis would systematically
overestimate the rate of low-expression genes by more than high-expression
ones, plausibly causing false-positive DE calls on the low end. This is a
**testable empirical claim** to add to the paper.

## What this validates

- **Cell-total consistency** (cross-referencer agent's prediction): exact, at
  numerical precision.
- **Per-gene rate bias**: confirmed, with an interpretable scaling pattern.
- **Spike-in identifiability**: confirmed (oracle recovers truth from 10%
  spike-ins).

## What this does NOT yet validate

- The translation of `mdrelax` Theorem 8 (rank-of-candidate-set-matrix
  identifiability) to the dense scRNA-seq case where the within-cell mechanism
  is what matters. The simulation here uses the simplest C2-violation model
  (logistic in log-y); more complex within-cell dependence patterns (gene-gene
  correlation in dropout, cell-quality effects) would need further study.
- Whether the relative-bias pattern (lower-mu → higher relative bias)
  persists under different dropout functions or holds for more realistic
  scRNA-seq DGPs (e.g. SymSim, Splat).

## Files

- `sim.R`: data generation, naive ZINB MLE, oracle estimation
- `run.R`: driver across `alpha_1` grid
- `results.rds`: full per-replicate results

## Next step

If proceeding to a paper: replicate this study using a more realistic
generative model (e.g., the Poisson-Beta two-state expression model used in
Jiang et al. 2022) and reach out to Jingyi Jessica Li at UCLA. She's the
last author of the anchor paper, actively works on this controversy, and
would be the natural collaborator. The methodological contribution
(coarsening framework + identifiability theorem port) is on the user's side;
the empirical motivation (real scRNA-seq datasets) is on hers.
