# V8: Real-data application on Tabula Muris spleen.
# Pipeline:
#   1. Load Spleen-counts.csv: 23,341 genes + 92 ERCCs x 1,718 cells
#   2. Match cell barcodes to annotations (cell types)
#   3. Filter cells (QC: >= 500 genes detected, >= 50,000 reads)
#   4. Filter genes (expressed in >= 10 cells); keep top variable genes
#   5. Load ERCC concentrations from Thermo Fisher table
#   6. Fit ERCC dropout function via fit_dropout_ercc_logistic
#   7. Run naive ZINB and oracle on selected genes
#   8. Estimate ERCC-vs-endogenous shift via housekeeping genes
#   9. Report findings

suppressPackageStartupMessages({
    library(data.table)
})
source("sim.R")

set.seed(20260511)

# -----------------------------------------------------------------
# 1. Load count matrix
# -----------------------------------------------------------------
cat("Loading Spleen-counts.csv...\n")
t0 <- Sys.time()
counts <- fread("tabula_muris/FACS/Spleen-counts.csv", data.table = FALSE)
gene_names <- counts[[1]]
counts <- as.matrix(counts[, -1, drop = FALSE])
rownames(counts) <- gene_names
cat("  ", nrow(counts), "rows (genes+ERCC) x", ncol(counts), "cells\n")
cat("  ", round(as.numeric(Sys.time() - t0), 1), "s\n")

# Split into genes and ERCC
ercc_mask <- grepl("^ERCC-", gene_names)
cat("  ERCC rows:", sum(ercc_mask), "\n")
X_ercc <- t(counts[ercc_mask, , drop = FALSE])  # cells x ERCC
X_genes <- t(counts[!ercc_mask, , drop = FALSE])  # cells x genes
cat("  X_genes shape:", dim(X_genes)[1], "x", dim(X_genes)[2], "\n")
cat("  X_ercc shape :", dim(X_ercc)[1], "x", dim(X_ercc)[2], "\n")

# -----------------------------------------------------------------
# 2. Cell QC
# -----------------------------------------------------------------
cell_total_reads <- rowSums(X_genes)
cell_genes_detected <- rowSums(X_genes > 0)
cat("\nCell QC summary:\n")
cat("  median total reads:", median(cell_total_reads), "\n")
cat("  median genes detected:", median(cell_genes_detected), "\n")
qc_pass <- cell_total_reads >= 50000 & cell_genes_detected >= 500
cat("  cells passing QC (>=50k reads, >=500 genes):", sum(qc_pass), "/", length(qc_pass), "\n")

X_genes <- X_genes[qc_pass, , drop = FALSE]
X_ercc  <- X_ercc[qc_pass, , drop = FALSE]
n_cells <- nrow(X_genes)

# -----------------------------------------------------------------
# 3. Gene filtering: keep genes expressed in >= 10 cells
# -----------------------------------------------------------------
gene_n_cells <- colSums(X_genes > 0)
keep_genes <- gene_n_cells >= 10
cat("\nGenes expressed in >=10 cells:", sum(keep_genes), "/", length(keep_genes), "\n")
X_genes <- X_genes[, keep_genes, drop = FALSE]

# Compute mean expression per gene
gene_means_obs <- colMeans(X_genes)
cat("  gene mean summary:", round(quantile(gene_means_obs, c(0.05, 0.5, 0.95)), 3), "\n")

# -----------------------------------------------------------------
# 4. ERCC concentrations
# -----------------------------------------------------------------
ercc_table <- fread("/tmp/ercc.tsv", data.table = FALSE)
# Match ERCC IDs in our matrix to the table
ercc_ids_in_data <- colnames(X_ercc)
ercc_match <- match(ercc_ids_in_data, ercc_table$`ERCC ID`)
cat("\nERCC ID match: ", sum(!is.na(ercc_match)), "/", length(ercc_match), "matched to concentration table\n")
# Keep only matched ERCCs
keep_ercc <- !is.na(ercc_match)
X_ercc <- X_ercc[, keep_ercc, drop = FALSE]
ercc_match <- ercc_match[keep_ercc]
ercc_conc_mix1 <- ercc_table$`concentration in Mix 1 (attomoles/ul)`[ercc_match]
cat("  ERCC concentration range (attomoles/ul): [", min(ercc_conc_mix1), ",", max(ercc_conc_mix1), "]\n")

# Convert ERCC concentration to expected count per cell. The exact conversion
# depends on lysis volume and dilution, but for fitting we only need RELATIVE
# concentrations. The ERCC mix follows ratios spanning ~6 orders of magnitude.
# For dropout fitting via marginal likelihood, we need expected count per cell.
# Approximate: rescale so that the highest-conc ERCC has expected count
# matching its observed non-zero mean.
ercc_obs_mean <- colMeans(X_ercc)
cat("  ERCC observed mean range:", round(quantile(ercc_obs_mean, c(0.05, 0.5, 0.95)), 2), "\n")
# Scale concentrations so that they correspond to expected pre-dropout counts
# Use the top-5 ERCCs as anchors: their observed mean should approximate
# their expected count * (1 - dropout_rate). Without a calibration, use a
# linear scaling: c_expected = ercc_conc_mix1 * scale_factor where
# scale_factor matches the median (above-zero) observed-vs-conc ratio.
nonzero_ercc <- ercc_obs_mean > 0.5
if (sum(nonzero_ercc) > 5) {
    scale_factor <- median(ercc_obs_mean[nonzero_ercc] / ercc_conc_mix1[nonzero_ercc])
} else {
    scale_factor <- 1
}
ercc_expected <- ercc_conc_mix1 * scale_factor
cat("  ERCC scale factor (mean-ratio):", signif(scale_factor, 3), "\n")
cat("  ERCC expected count range:", round(quantile(ercc_expected, c(0.05, 0.5, 0.95)), 3), "\n")

# -----------------------------------------------------------------
# 5. Fit ERCC dropout function via marginal likelihood
# -----------------------------------------------------------------
cat("\nFitting ERCC dropout function...\n")
t0 <- Sys.time()
drop_fit_ercc <- fit_dropout_ercc_logistic(ercc_expected, X_ercc)
cat("  fit time:", round(as.numeric(Sys.time() - t0), 1), "s\n")
cat("  ERCC-fit alpha_0:", round(drop_fit_ercc$params$alpha_0, 3),
    " alpha_1:", round(drop_fit_ercc$params$alpha_1, 3), "\n")
cat("  pi(y=0)   =", round(plogis(drop_fit_ercc$params$alpha_0), 3), "\n")
cat("  pi(y=10)  =", round(plogis(drop_fit_ercc$params$alpha_0 -
                                  drop_fit_ercc$params$alpha_1 * log(11)), 3), "\n")
cat("  pi(y=100) =", round(plogis(drop_fit_ercc$params$alpha_0 -
                                  drop_fit_ercc$params$alpha_1 * log(101)), 3), "\n")

# -----------------------------------------------------------------
# 6. Pick gene subset for fitting (top variable + housekeeping)
# -----------------------------------------------------------------
gene_var_obs <- apply(X_genes, 2, var)
gene_cv2 <- gene_var_obs / (colMeans(X_genes)^2 + 1e-3)
# Top 200 by mean expression (avoid noise-dominated low-expression)
top_by_mean <- order(gene_means_obs, decreasing = TRUE)[1:200]

# Mouse housekeeping genes commonly used for normalization
hk_genes <- c("Gapdh", "Actb", "B2m", "Hprt", "Tbp", "Ywhaz", "Hmbs", "Tfrc",
              "Sdha", "Ubc", "Pgk1", "Rpl13a", "Hprt1")
hk_present <- intersect(hk_genes, colnames(X_genes))
cat("\nHousekeeping genes present in dataset:", paste(hk_present, collapse = ", "), "\n")
hk_idx <- which(colnames(X_genes) %in% hk_present)

# Combined gene set
gene_subset_idx <- unique(c(top_by_mean, hk_idx))
cat("Gene subset for fitting:", length(gene_subset_idx), "genes\n")
X_sub <- X_genes[, gene_subset_idx, drop = FALSE]
gene_subset_names <- colnames(X_genes)[gene_subset_idx]

# -----------------------------------------------------------------
# 7. Naive ZINB and ERCC-oracle fits
# -----------------------------------------------------------------
cat("\nFitting naive ZINB (", ncol(X_sub), "genes)...\n")
t0 <- Sys.time()
naive <- fit_naive_zinb(X_sub)
cat("  ", round(as.numeric(Sys.time() - t0), 1), "s,",
    sum(is.na(naive[, "mu_hat"])), "failures\n")

cat("Fitting ERCC oracle (", ncol(X_sub), "genes)...\n")
t0 <- Sys.time()
oracle <- fit_oracle_zinb_v3(X_sub, drop_fit_ercc)
cat("  ", round(as.numeric(Sys.time() - t0), 1), "s,",
    sum(is.na(oracle[, "mu_hat"])), "failures\n")

# -----------------------------------------------------------------
# 8. Compare estimates
# -----------------------------------------------------------------
naive_mu <- naive[, "mu_hat"]
oracle_mu <- oracle[, "mu_hat"]
naive_pi  <- naive[, "pi_hat"]
obs_mean <- colMeans(X_sub)

# Cell-total consistency check
naive_pred_mean <- (1 - naive_pi) * naive_mu
cat("\n=== Cell-total consistency ===\n")
cat("  obs mean range :", round(quantile(obs_mean, c(0.05, 0.5, 0.95)), 3), "\n")
cat("  naive (1-pi)*mu:", round(quantile(naive_pred_mean, c(0.05, 0.5, 0.95)), 3), "\n")
cat("  median |diff|  :", signif(median(abs(naive_pred_mean - obs_mean), na.rm = TRUE), 3), "\n")

# Naive vs oracle ratios
ratio <- oracle_mu / naive_mu
cat("\n=== Oracle / Naive mu ratios ===\n")
cat("  ratio quantiles (5/25/50/75/95):",
    round(quantile(ratio, c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE), 2), "\n")

# Housekeeping gene results
cat("\n=== Housekeeping genes ===\n")
hk_in_subset <- which(gene_subset_names %in% hk_present)
cat(sprintf("  %-10s  %10s  %10s  %10s  %10s\n",
            "gene", "obs_mean", "naive_mu", "oracle_mu", "oracle/naive"))
for (i in hk_in_subset) {
    cat(sprintf("  %-10s  %10.2f  %10.2f  %10.2f  %10.3f\n",
                gene_subset_names[i],
                obs_mean[i], naive_mu[i], oracle_mu[i],
                oracle_mu[i] / naive_mu[i]))
}

# -----------------------------------------------------------------
# 9. Estimate ERCC-vs-endogenous shift heuristic
# -----------------------------------------------------------------
# Logic: if oracle is unbiased, oracle_mu / observed_mean should match
# 1/(1 - pi_endo(mu)). If pi_ERCC and pi_endo differ by a shift, the
# observed oracle/observed ratio will be off.
# A simple heuristic: compute the shift such that the oracle estimate of
# the highest-expressed (most identifiable) genes matches their observed
# mean as if dropout were minimal.
hk_obs <- obs_mean[hk_in_subset]
hk_naive <- naive_mu[hk_in_subset]
hk_oracle <- oracle_mu[hk_in_subset]
hk_naive_pi <- naive_pi[hk_in_subset]
cat("\n=== ERCC-vs-endogenous shift estimate (housekeeping-based) ===\n")
cat("  housekeeping naive pi:", round(hk_naive_pi, 3), "\n")
cat("  Implied dropout from naive (1-(1-pi)) =", round(hk_naive_pi, 3), "\n")
# The naive (1-pi)*mu = obs_mean by construction. So pi tells us the dropout share.
# Compare to ERCC-fit dropout at observed mean expressions.
hk_ercc_dropout <- plogis(drop_fit_ercc$params$alpha_0 -
                           drop_fit_ercc$params$alpha_1 * log(hk_obs + 1))
cat("  ERCC-fit dropout at obs_mean :", round(hk_ercc_dropout, 3), "\n")
cat("  difference (ERCC - naive)    :", round(hk_ercc_dropout - hk_naive_pi, 3), "\n")
cat("  median ERCC-naive dropout difference:",
    round(median(hk_ercc_dropout - hk_naive_pi, na.rm = TRUE), 3), "\n")

# Save results
saveRDS(list(
    drop_fit_ercc = drop_fit_ercc,
    naive = naive, oracle = oracle,
    obs_mean = obs_mean,
    gene_subset_names = gene_subset_names,
    hk_genes_present = hk_present,
    n_cells = n_cells,
    n_ercc = ncol(X_ercc),
    n_genes_total = ncol(X_genes),
    ercc_expected_conc = ercc_expected
), "results_v8_tabula_muris.rds")
cat("\nDone. Saved results_v8_tabula_muris.rds\n")
