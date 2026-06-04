# Real-data application: feasibility investigation

**Goal**: determine whether real-data application of the methodology is
feasible without external collaboration, what dataset to use, and what
existing methods to position against.

---

## 1. Accessible ERCC-equipped public scRNA-seq datasets

| Dataset | Cells | ERCCs | Access | Tissues / cell types | Notes |
|---|---|---|---|---|---|
| **Tabula Muris** (FACS Smart-seq2) | 53,760 | 92 | [figshare DOI 10.6084/m9.figshare.5715040](https://doi.org/10.6084/m9.figshare.5715040), Bioconductor `TabulaMurisData::TabulaMurisSmartSeq2()` | 20 mouse organs, all annotated cell types | Largest, best-curated; ERCC counts in `altExp(sce, "ERCC")` |
| Lun 416B cell line | ~384 | 92 | Bioconductor `scRNAseq::LunSpikeInData()` | Mouse 416B cells, oncogene perturbation | Small, two-condition DE benchmark |
| Tung iPSC | 864 | 92 | Bioconductor `scRNAseq::TungData()` | Human iPSCs, three individuals | Fluidigm C1 |
| Kolodziejczyk ESC | 704 | 92 | Bioconductor `scRNAseq::KolodziejczykESCData()` | Mouse ESCs, three culture conditions | Smart-seq2 |
| Buettner ESC (cell-cycle) | 288 | 92 | Bioconductor `scRNAseq::BuettnerESCData()` | Mouse ESCs, three cell-cycle phases | Fluidigm C1 |

**Verdict**: real-data application is fully feasible without collaboration.
**Tabula Muris is the natural primary dataset** (size, curation, multiple
tissues, well-annotated cell types). For paper presentation, picking 1-2
tissues (e.g., spleen, marrow) and showing the methodology end-to-end is
sufficient.

ERCC concentration values (the c_k vector required by
`fit_dropout_ercc_logistic`) are bundled in the `ERCC` annotation file
shipped by Thermo Fisher (publication 4455352, ERCC92_data table) and
mirrored in many R packages (e.g., `absSimSeq::ERCC92_data`).

---

## 2. Existing methods that overlap with the framework

The methodology is **not novel as an ERCC-based dropout-correction
estimator**. Two existing methods solve closely related problems:

### TASC (Jia et al. 2017, Nucleic Acids Research)

- Empirical Bayes hierarchical model
- Uses ERCC spike-ins to estimate cell-specific dropout rates and
  amplification bias
- Borrows information across cells via a shared distribution
- Hierarchical mixture model for biological variance + DE detection
- [Paper](https://academic.oup.com/nar/article/45/19/10978/4210934),
  [GitHub](https://github.com/scrna-seq/TASC)

### DECENT (Ye et al. 2019, Bioinformatics)

- Beta-binomial capture model (gene- and cell-specific)
- Differential expression specific
- ERCC-calibrated, but works without spike-ins too
- Explicitly acknowledges the limitation that **ERCC vs endogenous
  transcripts have different capture efficiencies** (poly-A length, GC
  content, sequence properties)
- [Paper](https://academic.oup.com/bioinformatics/article/35/24/5155/5514046)

### What this means for the paper

The paper's contribution **is not "build an ERCC-based dropout corrector"**.
Instead, it is:

1. **Theoretical unification**: framing scRNA-seq dropout as a coarsening
   problem and importing identifiability theorems from masked-data
   reliability statistics. TASC and DECENT solve the inference problem
   without naming it as masked-data; the unification provides a precise
   identifiability characterization (Theorem 8 ports cleanly).
2. **Form-robustness theorems**: simulation evidence (v3a, v3b) showing
   spike-in correction works regardless of dropout function family
   (logistic vs exponential), with a 1-3pp cost from form mismatch.
   Neither TASC nor DECENT investigates this systematically.
3. **Cell-total consistency** (v1, v2, v3a): a structural property that
   holds for any zero-inflated estimator under any dropout function.
   This explains why naive ZINB still predicts observed counts well even
   when per-gene mu is biased: the saturated marginal mean is preserved.
4. **Cell-type heterogeneity robustness** (v4-celltypes, v5-cluster):
   oracle recovers population marginal AND per-cluster mu within ~2%
   bias. This is a clean robustness result that follows from ZINB's
   mean-matching property.
5. **Spike-in fit attenuation diagnosis** (v5-diagnostic) and proposed
   regularization (v6-regularization, in progress): identifies a
   practical pitfall of MLE-based spike-in fitting when alpha_1 is near
   the boundary, and shows that LRT-based selection between constant
   and Y-dependent dropout fixes it.

**Honest paper claim**: "We provide a theoretical framework that explains
when ERCC-based dropout correction works and when it fails, importing
identifiability machinery from a different statistical literature
(masked-data reliability inference). Existing tools (TASC, DECENT) are
recovered as special cases of the framework. The framework predicts and
the simulations confirm cell-total consistency, form-robustness, and
heterogeneity-robustness, none of which are systematically investigated
in prior work."

This is a less-bold but more-defensible positioning. It is also a more
sustainable contribution: a framework paper is harder to scoop than a
methods paper.

---

## 3. ERCC-vs-endogenous capture-efficiency difference (the critical caveat)

**The DECENT paper explicitly flags**: "ERCC molecules ... are captured
and amplified less efficiently than endogenous transcripts ... models
estimated using spike-ins may not be entirely appropriate for endogenous
transcripts".

In our framework's vocabulary: the dropout function for ERCC,
pi_ERCC(y), is not the same as the dropout function for endogenous
transcripts, pi_endo(y). Spike-in fitting recovers pi_ERCC, but applying
it to endogenous genes is approximate.

This is a substantive caveat for the paper. Three response strategies:

1. **Acknowledge and bound**: report the magnitude of expected
   ERCC-vs-endogenous difference from existing literature (typically
   ~2-3x in capture efficiency) and provide a sensitivity analysis: how
   much does the methodology degrade when the true endogenous dropout
   differs from the spike-in-fit dropout?
2. **Bias-correct**: introduce an adjustment factor c (say from
   literature) such that pi_endo(y) = pi_ERCC(c * y). Estimate c from
   data using genes with known calibration (housekeeping genes), or
   take c from literature.
3. **Ignore in the methodology paper, address in the application**:
   document the limitation and demonstrate empirically that even with
   this caveat, the corrected estimates are useful (e.g., cell-type
   markers are sharper than under naive).

Strategy 1 is cleanest and adds a v7 simulation: misspecified-target
oracle. v7 would test what happens when the dropout function differs
between training (spike-ins) and test (endogenous genes).

---

## 4. Recommended next steps (no collaboration required)

**Phase A (cheap, do next)**: v7 misspecified-target sensitivity analysis.
Generate data with different dropout functions for spike-ins vs genes
(varying ratios of capture efficiency), apply the spike-in oracle to
genes, quantify bias as a function of the spike-in-vs-endogenous gap.
This produces the practitioner-facing "if ERCC capture is X% lower than
endogenous, your bias is Y%" table the paper needs.

**Phase B (engineering, ~1 day)**: Tabula Muris real-data application.
Download the spleen subset, fit the methodology (naive ZINB, plain
spike-in oracle, regularized spike-in oracle), compare to:
- Cell-type marker gene expression (do oracle estimates sharpen
  marker-gene contrast?)
- Published TASC and DECENT results (run their R packages on the same
  data and compare estimates)
- Known biological controls (housekeeping genes should have low
  variance)

**Phase C (paper-writing, ~1-2 weeks)**: draft `papers/scrna-coarsening/`
based on:
- The bridge analysis (`scrna-masked-bridge-report.md`)
- The full simulation portfolio (v1-v7)
- The Tabula Muris real-data application
- Direct comparison with TASC and DECENT

The Jingyi Li outreach can happen at any of these stages, but is not
required for any of them. The paper can stand on simulation + Tabula
Muris + comparison with TASC/DECENT alone.

---

## 5. Updated prior-art table

| Reference | What it does | Differentiator from this work |
|---|---|---|
| Jia et al. 2017 (TASC) | ERCC-based empirical-Bayes dropout estimation | Hierarchical estimator, not a framework. Doesn't formalize identifiability or form-robustness. |
| Ye et al. 2019 (DECENT) | Beta-binomial capture for DE | DE-specific; explicitly notes ERCC-vs-endogenous limitation; doesn't analyze cell-type heterogeneity robustness. |
| Risso et al. 2018 (ZINB-WaVE) | ZINB factor model with covariates | Allows logit-linear dropout in covariates; doesn't use ERCC directly; doesn't analyze identifiability. |
| Finak et al. 2015 (MAST) | Hurdle model for DE | Two-part model; treats dropout as observed; doesn't use ERCC. |
| Pierson & Yau 2015 (ZIFA) | Factor analysis with zero inflation | Dropout is function of expression magnitude; no ERCC integration; doesn't analyze identifiability. |
| **This work** | Theoretical unification via masked-data coarsening | Identifiability theorems, form-robustness, cell-type-heterogeneity robustness, regularization diagnosis. Explains when existing methods work and provides a framework for new methods. |

---

## 6. Conclusion of investigation

Real-data application is feasible. The methodology is novel as a **theoretical framework** but not as a **dropout-correction tool**. The
paper should reposition as a framework paper that subsumes TASC and
DECENT as special cases and adds: identifiability characterization,
form-robustness, heterogeneity-robustness, and regularization. The
critical caveat (ERCC-vs-endogenous capture difference) needs explicit
sensitivity analysis (v7) before the paper is submission-ready.

The Jingyi Li outreach is no longer essential. The paper can be drafted
and submitted solo with Tabula Muris as the empirical anchor.
