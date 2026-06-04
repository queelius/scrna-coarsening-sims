# Validation v7: ERCC-vs-endogenous capture-efficiency sensitivity

**The DECENT paper notes**: ERCC molecules and endogenous transcripts have
different capture efficiencies (poly-A length, sequence properties). If we
fit dropout from ERCC and apply to genes, we are using pi_ERCC instead of
the true pi_endo. **How much does this gap cost?**

This is the most important caveat for any ERCC-based methodology. v7
quantifies it directly.

## Setup

- 2000 cells, 50 genes (log-normal mu), Poisson-Beta DGP
- Endogenous dropout: pi_endo(y) = plogis(1.5 - 0.6 * log(y+1))
- ERCC dropout:       pi_ERCC(y) = plogis(1.5 + shift - 0.6 * log(y+1))
- Positive shift = ERCC has more dropout than endogenous (the realistic
  direction per DECENT)
- Sweep shift in {-1, -0.5, 0, 0.5, 1, 1.5, 2}
- 92 ERCC at log-uniform concentrations from 0.5 to 50
- 25 reps per condition

## Per-gene relative bias (median across 50 genes)

| shift | naive | oracle ERCC | oracle TRUE-endo (reference) |
|---|---|---|---|
| -1.0 | +60.5% | **-39.0%** | +2.0% |
| -0.5 | +64.3% | -23.3% | +2.3% |
| **0.0** | +61.8% | **+2.0%** | +2.1% |
| +0.5 | +62.6% | +34.4% | +2.0% |
| +1.0 | +63.2% | **+57.5%** | +2.2% |
| +1.5 | +66.9% | +73.2% | +2.4% |
| +2.0 | +63.0% | **+80.7%** | +2.2% |

## Headline findings

1. **At matched dropout (shift = 0)**: oracle ERCC matches oracle true-endo
   within MC error. This is the v3d result restated.

2. **The "useful zone" for ERCC-based correction is |shift| <= 0.5**:
   - At shift = +-0.5, oracle bias is 23-34%, still substantially better
     than naive's 62%
   - At shift = +-1.0, oracle bias is comparable to naive
   - At |shift| > 1, **oracle is worse than naive**

3. **The bias direction is informative**:
   - shift > 0 (ERCC has more dropout): oracle OVER-estimates mu (positive
     bias) because it thinks dropout is more severe than reality
   - shift < 0 (ERCC has less dropout): oracle UNDER-estimates mu

4. **ERCC dropout-fit recovery is precise**: alpha_0 estimates match the
   ERCC truth (alpha_0_endo + shift) at every shift level to within 0.01.
   The fit is doing exactly what it should: recovering pi_ERCC. The bias
   is from applying pi_ERCC to genes whose true dropout is pi_endo.

5. **Mu-stratified bias is dramatic at high shift**:

| shift | oracle ERCC bias Q1 (low mu) | Q5 (high mu) | spread |
|---|---|---|---|
| 0.0 | +0.2% | +3.0% | 3pp |
| +1.0 | +106% | +33% | 73pp |
| +2.0 | +221% | +47% | 174pp |

   Low-expression genes are catastrophically affected when the ERCC-endo
   gap is large.

## How realistic are these shift values?

From the literature on ERCC capture efficiency in scRNA-seq:
- DECENT paper: endogenous "much more efficiently captured than ERCC"
- Typical capture-efficiency ratio: 2-3x for endogenous vs ERCC
- In logit terms, this corresponds to shift in [log(2), log(3)] = [0.7, 1.1]
- **Real-world ERCC-vs-endogenous shifts are right at the "useful zone" boundary**

This means the methodology, applied naively to real ERCC data, would
likely give modest improvement over naive ZINB at best. The headline
finding for the paper is: **ERCC-based correction is useful only when
ERCC and endogenous capture efficiencies are well-matched, which they
typically are not in practice.**

## What this means for the paper

The paper now has three options for addressing this caveat:

### Option A: acknowledge as fundamental limit (most honest)

Frame the methodology as a **theoretical framework** that recovers the
true dropout function from spike-ins. Empirically, this is useful only
when the spike-in capture matches the endogenous capture. **All ERCC-
based methods (TASC, DECENT, ours) face this same fundamental limit.**
The framework's contribution is identifiability theorems and form-
robustness, not solving the ERCC-endo transferability problem.

This would be a clear, defensible, and honest paper. It gives
practitioners a sensitivity table they can use to decide whether to
apply spike-in correction or not.

### Option B: propose a calibration correction

Estimate the ERCC-endo shift from data using housekeeping genes (genes
with known stable expression across conditions, e.g., GAPDH, ACTB). Fit
the shift such that the oracle estimate of housekeeping-gene mu matches
the population mean.

This is more work (~1 week) but gives the paper a stronger story:
"existing methods have this problem, our framework also identifies it
and proposes a calibration".

### Option C: position as upstream-of-correction tool

Frame the methodology not as a corrector, but as a **diagnostic** for
detecting when spike-in correction is appropriate. The framework
predicts cell-total consistency and identifies ERCC-vs-endogenous shifts
through the gap between predicted and observed marginals.

This is more abstract but harder to scoop and more conceptually clean.

## Recommendation

**Option A (acknowledge as fundamental limit) for the first paper, with
Option B (calibration correction) as natural follow-up work.**

Reasons:
- Option A is honest and complete with the current simulation portfolio.
- Option B is a substantial extension that could be its own paper.
- The first paper's contribution is the framework + identifiability +
  form/heterogeneity-robustness, all of which stand independent of the
  ERCC-endo issue.

The paper's claim becomes:

> Within the regime where spike-in capture matches endogenous capture
> (shift in [-0.5, +0.5]), the methodology recovers per-gene mu within
> 2% bias under matched conditions and within 23-34% bias under modest
> mismatch. Beyond this regime (|shift| >= 1), spike-in correction is
> comparable to or worse than no correction. Realistic ERCC-vs-endo
> shifts in current scRNA-seq protocols are at the boundary of this
> useful zone, suggesting that future protocol improvements (better-
> matched spike-ins) or methodological extensions (calibration via
> housekeeping genes) would substantially improve practical utility.

This is a well-bounded, honest, and useful methodological contribution.

## What this changes about real-data application

The Tabula Muris real-data application becomes more nuanced. We can no
longer claim "oracle gives unbiased per-gene mu". Instead, we can:

1. Apply oracle and compute the implied (alpha_0_ERCC, alpha_1) values.
2. Compute housekeeping-gene mu under both naive and oracle.
3. Note that the oracle-vs-naive difference for housekeeping genes
   reflects the ERCC-endo capture gap.
4. Use this to back out an estimated shift and report a corrected
   sensitivity range.

This is more honest and more interesting than just running the
methodology and reporting numbers.

## Synthesis: validation status after v7

| Property | Status |
|---|---|
| Cell-total consistency | exact, structural |
| Per-gene oracle recovery, matched dropout | ~2% bias |
| Per-gene oracle recovery, mismatched dropout form | 3-6% bias |
| **Per-gene oracle recovery, ERCC-endo capture gap** | **23-34% at \|shift\|=0.5; equal to naive at \|shift\|=1; worse at \|shift\|>1** |
| Spike-in fraction requirements | as low as 0.5% |
| ERCC realistic design (matched capture) | identical to perfect spike |
| Power gain on DE testing | 10pp at moderate signals |
| Cell-level capture variation (within-cell) | oracle robust to sigma_q = 1.0 |
| Cell-type heterogeneity, population mu | oracle ~2% bias |
| Per-cluster mu recovery | oracle ~2% bias per cluster |
| Spike-fit attenuation (alpha_1=0 boundary) | exists; LRT regularization fails; BIC/shrinkage TBD |
| **Methodology has zero intrinsic bias when dropout known** | **confirmed (v5)** |
| **Useful zone for spike-in correction** | **\|shift\| <= 0.5 (matched-capture regime)** |

## Files

- `run_v7_ercc_endo.R`: driver
- `results_v7_ercc_endo.rds`: full per-rep results
- `log_v7_ercc_endo.txt`: stdout log
