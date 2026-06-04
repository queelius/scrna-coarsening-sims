# V9: Compare framework's oracle vs TASC and DECENT on Tabula Muris spleen.
# Runs after BiocManager::install("DECENT", ...) succeeds.

suppressPackageStartupMessages({
    library(data.table)
})
source("sim.R")

set.seed(20260512)

# -----------------------------------------------------------------
# Load Tabula Muris spleen (same setup as v8)
# -----------------------------------------------------------------

cat("Loading data...\n")
counts <- fread("tabula_muris/FACS/Spleen-counts.csv", data.table = FALSE)
gene_names <- counts[[1]]
counts <- as.matrix(counts[, -1, drop = FALSE])
rownames(counts) <- gene_names
ercc_mask <- grepl("^ERCC-", gene_names)
X_ercc_full <- t(counts[ercc_mask, , drop = FALSE])
X_genes_full <- t(counts[!ercc_mask, , drop = FALSE])

# QC
qc_pass <- rowSums(X_genes_full) >= 50000 & rowSums(X_genes_full > 0) >= 500
X_ercc <- X_ercc_full[qc_pass, , drop = FALSE]
X_genes <- X_genes_full[qc_pass, , drop = FALSE]
n_cells <- nrow(X_genes)

# Filter genes
keep_genes <- colSums(X_genes > 0) >= 10
X_genes <- X_genes[, keep_genes, drop = FALSE]

# ERCC concentrations
ercc_table <- fread("/tmp/ercc.tsv", data.table = FALSE)
ercc_match <- match(colnames(X_ercc), ercc_table$`ERCC ID`)
keep_ercc <- !is.na(ercc_match)
X_ercc <- X_ercc[, keep_ercc, drop = FALSE]
ercc_conc_mix1 <- ercc_table$`concentration in Mix 1 (attomoles/ul)`[ercc_match[keep_ercc]]
nonzero_ercc <- colMeans(X_ercc) > 0.5
scale_factor <- median(colMeans(X_ercc)[nonzero_ercc] / ercc_conc_mix1[nonzero_ercc])
ercc_expected <- ercc_conc_mix1 * scale_factor

# Pick gene subset (same as v8 for comparability)
gene_means_obs <- colMeans(X_genes)
top_by_mean <- order(gene_means_obs, decreasing = TRUE)[1:200]
hk_genes <- c("Gapdh", "Actb", "B2m", "Hprt", "Tbp", "Ywhaz", "Hmbs", "Tfrc",
              "Sdha", "Ubc", "Pgk1", "Rpl13a")
hk_idx <- which(colnames(X_genes) %in% hk_genes)
gene_subset_idx <- unique(c(top_by_mean, hk_idx))
X_sub <- X_genes[, gene_subset_idx, drop = FALSE]
gene_subset_names <- colnames(X_genes)[gene_subset_idx]

cat("Gene subset:", length(gene_subset_idx), "genes,", n_cells, "cells\n")

# -----------------------------------------------------------------
# Our framework's oracle (from v8)
# -----------------------------------------------------------------

cat("\n=== Our framework: ERCC-fit oracle ===\n")
t0 <- Sys.time()
drop_fit_ercc <- fit_dropout_ercc_logistic(ercc_expected, X_ercc)
cat("Dropout fit:", round(as.numeric(Sys.time() - t0), 1), "s,",
    "alpha_0 =", round(drop_fit_ercc$params$alpha_0, 3),
    "alpha_1 =", round(drop_fit_ercc$params$alpha_1, 3), "\n")

t0 <- Sys.time()
naive <- fit_naive_zinb(X_sub)
oracle <- fit_oracle_zinb_v3(X_sub, drop_fit_ercc)
cat("Naive + oracle fit:", round(as.numeric(Sys.time() - t0), 1), "s\n")

# -----------------------------------------------------------------
# DECENT
# -----------------------------------------------------------------

decent_mu <- rep(NA_real_, ncol(X_sub))
if (requireNamespace("DECENT", quietly = TRUE)) {
    cat("\n=== DECENT comparison ===\n")
    library(DECENT)
    # DECENT expects: a count matrix + ERCC counts (optional).
    # The decent() function fits per-gene parameters with capture-aware NB.
    # Reference: Ye, Speed, Salim 2019.
    # API may differ; check vignette before final run.
    t0 <- Sys.time()
    tryCatch({
        # Need to pass full count matrix not subset to use DECENT properly,
        # but for runtime keep to subset.
        # decent_result <- decent(data.obs = t(X_sub), X = NULL, use.spikes = FALSE)
        # Returning: per-gene mean estimates accessible via summary()
        cat("(DECENT API integration TBD; need to read package vignette)\n")
    }, error = function(e) {
        cat("DECENT failed:", conditionMessage(e), "\n")
    })
    cat("DECENT total:", round(as.numeric(Sys.time() - t0), 1), "s\n")
} else {
    cat("\nDECENT not installed; skipping comparison\n")
}

# -----------------------------------------------------------------
# TASC
# -----------------------------------------------------------------

tasc_mu <- rep(NA_real_, ncol(X_sub))
if (requireNamespace("TASC", quietly = TRUE)) {
    cat("\n=== TASC comparison ===\n")
    library(TASC)
    # TASC expects: cell-by-gene matrix + ERCC matrix
    # Reference: Jia et al. 2017
    # The package fits per-cell dropout via empirical Bayes
    cat("(TASC API integration TBD; need to read package vignette)\n")
} else {
    cat("\nTASC not installed; skipping comparison\n")
}

# -----------------------------------------------------------------
# Summary table
# -----------------------------------------------------------------

obs_mean <- colMeans(X_sub)
hk_in_subset <- which(gene_subset_names %in% hk_genes)

cat("\n=== Housekeeping comparison ===\n")
cat(sprintf("  %-10s  %10s  %10s  %10s  %10s  %10s\n",
            "gene", "obs_mean", "naive_mu", "ours", "DECENT", "TASC"))
for (i in hk_in_subset) {
    cat(sprintf("  %-10s  %10.2f  %10.2f  %10.2f  %10s  %10s\n",
                gene_subset_names[i], obs_mean[i],
                naive[i, "mu_hat"], oracle[i, "mu_hat"],
                ifelse(is.na(decent_mu[i]), "n/a", sprintf("%.2f", decent_mu[i])),
                ifelse(is.na(tasc_mu[i]), "n/a", sprintf("%.2f", tasc_mu[i]))))
}

saveRDS(list(naive = naive, oracle = oracle, decent_mu = decent_mu, tasc_mu = tasc_mu,
             obs_mean = obs_mean, gene_subset_names = gene_subset_names,
             hk_in_subset = hk_in_subset),
        "results_v9_method_comparison.rds")
cat("\nDone. Saved results_v9_method_comparison.rds\n")
