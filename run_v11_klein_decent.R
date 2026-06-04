# V11: DECENT vs spike-fit oracle on Klein 2015 inDrops K562 dataset.
# Addresses methodology-auditor finding M-A4 (UMI-protocol head-to-head
# with DECENT).
#
# Dataset: GSE65525 (Klein et al. 2015 Cell), inDrops UMI scRNA-seq.
# We use the K562_cells file (240 cells, 23k+ genes + 102 ERCC),
# which contains both endogenous transcripts and ERCC spike-ins in the
# same UMI-counting protocol that DECENT is designed for.
#
# Pipeline:
#   1. Load K562_cells.csv (gene-by-cell, includes ERCC rows)
#   2. Cell QC (>= 200 genes detected, >= 1000 reads -- UMI scale)
#   3. Gene filter (expressed in >= 10 cells); subset to 2000 top by mean
#   4. Fit ERCC dropout via fit_dropout_ercc_logistic (our spike-fit oracle)
#   5. Run DECENT::fitNoDE on the same subset
#   6. Compare per-gene mu estimates: oracle vs DECENT
#   7. Save side-by-side results

suppressPackageStartupMessages({
    library(data.table)
})
source("sim.R")

# DECENT may be slow or memory-heavy; load only when needed.
require_decent <- function() {
    if (!requireNamespace("DECENT", quietly = TRUE)) {
        stop("DECENT not installed")
    }
    invisible(NULL)
}

set.seed(20260523)

# ----------------------------------------------------------------------
# 1. Load Klein 2015 K562_cells.csv
# ----------------------------------------------------------------------
cat("Loading Klein 2015 K562_cells.csv...\n")
t0 <- Sys.time()
counts <- fread("klein2015/GSM1599500_K562_cells.csv", data.table = FALSE,
                header = FALSE)
gene_names <- counts[[1]]
counts <- as.matrix(counts[, -1, drop = FALSE])
rownames(counts) <- gene_names
mode(counts) <- "integer"
cat("  ", nrow(counts), "rows x", ncol(counts), "cells in",
    round(as.numeric(Sys.time() - t0), 1), "s\n")
cat("  count range:", range(counts), "\n")

# Drop any rows with empty/NA gene names
keep_named <- !is.na(gene_names) & nchar(gene_names) > 0
counts <- counts[keep_named, , drop = FALSE]
gene_names <- gene_names[keep_named]

ercc_mask <- grepl("^ERCC-", gene_names)
cat("  ERCC rows:", sum(ercc_mask), "\n")
X_ercc  <- t(counts[ercc_mask,  , drop = FALSE])  # cells x ERCC
X_genes <- t(counts[!ercc_mask, , drop = FALSE])  # cells x genes
cat("  X_genes:", dim(X_genes)[1], "x", dim(X_genes)[2], "\n")
cat("  X_ercc :", dim(X_ercc)[1],  "x", dim(X_ercc)[2],  "\n")

# ----------------------------------------------------------------------
# 2. Cell QC (UMI scale)
# ----------------------------------------------------------------------
cell_total <- rowSums(X_genes)
cell_genes_det <- rowSums(X_genes > 0)
cat("\nCell QC:\n")
cat("  median total UMIs:", median(cell_total), "\n")
cat("  median genes detected:", median(cell_genes_det), "\n")
qc_pass <- cell_total >= 1000 & cell_genes_det >= 200
cat("  cells passing QC:", sum(qc_pass), "/", length(qc_pass), "\n")
X_genes <- X_genes[qc_pass, , drop = FALSE]
X_ercc  <- X_ercc[qc_pass,  , drop = FALSE]
n_cells <- nrow(X_genes)

# ----------------------------------------------------------------------
# 3. Gene filter and subset (top 2000 by mean for tractability)
# ----------------------------------------------------------------------
gene_n_cells <- colSums(X_genes > 0)
keep_genes <- gene_n_cells >= 10
cat("\nGenes expressed in >=10 cells:", sum(keep_genes), "/", length(keep_genes), "\n")
X_genes <- X_genes[, keep_genes, drop = FALSE]

gene_means_obs <- colMeans(X_genes)
# Pick top 2000 by mean for tractable DECENT runtime
n_target <- min(2000, ncol(X_genes))
top_idx <- order(gene_means_obs, decreasing = TRUE)[1:n_target]
X_sub <- X_genes[, top_idx, drop = FALSE]
gene_sub_names <- colnames(X_genes)[top_idx]
cat("Subsetting to top", n_target, "genes by mean\n")
cat("  subset mean range:", round(range(colMeans(X_sub)), 3), "\n")

# ----------------------------------------------------------------------
# 4. ERCC dropout fit (our oracle)
# ----------------------------------------------------------------------
# For Klein 2015, ERCC concentrations are nominal Mix 1 attomoles
# from the standard Thermo Fisher table. Use the same /tmp/ercc.tsv.
ercc_table <- fread("/tmp/ercc.tsv", data.table = FALSE)
ercc_match <- match(colnames(X_ercc), ercc_table$`ERCC ID`)
keep_ercc <- !is.na(ercc_match)
X_ercc_sub <- X_ercc[, keep_ercc, drop = FALSE]
ercc_match <- ercc_match[keep_ercc]
ercc_conc <- ercc_table$`concentration in Mix 1 (attomoles/ul)`[ercc_match]
cat("\nERCC matched:", ncol(X_ercc_sub), "\n")
cat("  ERCC concentration range:", signif(min(ercc_conc), 2), "to",
    signif(max(ercc_conc), 2), "\n")

ercc_obs_mean <- colMeans(X_ercc_sub)
nonzero <- ercc_obs_mean > 0.5
scale_factor <- if (sum(nonzero) > 5) {
    median(ercc_obs_mean[nonzero] / ercc_conc[nonzero])
} else 1
ercc_expected <- ercc_conc * scale_factor
cat("  scale factor:", signif(scale_factor, 3), "\n")
cat("  ERCC expected count range:", round(quantile(ercc_expected, c(0.05, 0.5, 0.95)), 3), "\n")

cat("\nFitting ERCC dropout (spike-fit oracle)...\n")
t0 <- Sys.time()
drop_fit <- fit_dropout_ercc_logistic(ercc_expected, X_ercc_sub)
cat("  fit time:", round(as.numeric(Sys.time() - t0), 1), "s\n")
cat("  alpha_0:", round(drop_fit$params$alpha_0, 3),
    "  alpha_1:", round(drop_fit$params$alpha_1, 3), "\n")
cat("  pi(y=0)=",  round(plogis(drop_fit$params$alpha_0), 3),
    " pi(y=10)=", round(plogis(drop_fit$params$alpha_0 -
                                drop_fit$params$alpha_1 * log(11)), 3),
    " pi(y=100)=", round(plogis(drop_fit$params$alpha_0 -
                                drop_fit$params$alpha_1 * log(101)), 3), "\n")

cat("\nRunning spike-fit oracle on", ncol(X_sub), "genes...\n")
t0 <- Sys.time()
oracle <- fit_oracle_zinb_v3(X_sub, drop_fit)
cat("  ", round(as.numeric(Sys.time() - t0), 1), "s,",
    sum(is.na(oracle[, "mu_hat"])), "failures\n")

cat("\nFitting naive ZINB on", ncol(X_sub), "genes...\n")
t0 <- Sys.time()
naive <- fit_naive_zinb(X_sub)
cat("  ", round(as.numeric(Sys.time() - t0), 1), "s,",
    sum(is.na(naive[, "mu_hat"])), "failures\n")

# ----------------------------------------------------------------------
# 5. DECENT
# ----------------------------------------------------------------------
cat("\nLoading DECENT...\n")
require_decent()
library(DECENT)

# DECENT expects: data (genes x cells), spikes (genes x cells), cell.type
# K562 is a single cell line so all cells share one cell.type label.
# DECENT::decent() fits the full no-DE model when CE.range is provided.
# Use fitNoDE for unconditional dropout-corrected fit.

data_mat <- t(X_sub)            # genes x cells
spikes_mat <- t(X_ercc_sub)     # ERCC genes x cells
cat("  DECENT input: data =", dim(data_mat)[1], "x", dim(data_mat)[2],
    "  spikes =", dim(spikes_mat)[1], "x", dim(spikes_mat)[2], "\n")

# Reduce gene set further if needed for runtime: subset to 500 for first attempt
n_decent <- min(500, nrow(data_mat))
if (nrow(data_mat) > n_decent) {
    cat("  Reducing to top", n_decent, "genes for DECENT runtime budget\n")
    # match the subset to oracle for fair comparison
    oracle_means_for_subset <- colMeans(X_sub)
    decent_pick <- order(oracle_means_for_subset, decreasing = TRUE)[1:n_decent]
    data_mat <- data_mat[decent_pick, , drop = FALSE]
    gene_sub_decent <- gene_sub_names[decent_pick]
} else {
    gene_sub_decent <- gene_sub_names
}
cat("  DECENT genes:", nrow(data_mat), "\n")

cat("\nFitting DECENT fitNoDE...\n")
t0 <- Sys.time()
decent_fit <- tryCatch({
    DECENT::fitNoDE(
        data.obs    = data_mat,
        spikes      = spikes_mat,
        spike.conc  = ercc_expected,
        use.spikes  = TRUE,
        CE.range    = c(0.02, 0.1),
        tau.init    = c(-5, 0),
        tau.global  = TRUE,
        tau.est     = "endo",
        normalize   = "ML",
        GQ.approx   = TRUE,
        maxit       = 20,
        parallel    = FALSE
    )
}, error = function(e) {
    cat("DECENT failed:", conditionMessage(e), "\n")
    NULL
})
elapsed <- round(as.numeric(Sys.time() - t0), 1)
cat("  DECENT runtime:", elapsed, "s\n")

# ----------------------------------------------------------------------
# 6. Save results
# ----------------------------------------------------------------------
results <- list(
    dataset = "Klein2015 GSE65525 K562_cells",
    n_cells = n_cells,
    n_ercc  = ncol(X_ercc_sub),
    n_genes_subset = ncol(X_sub),
    drop_fit = drop_fit,
    naive = naive,
    oracle = oracle,
    decent = decent_fit,
    decent_elapsed_s = elapsed,
    gene_sub_names = gene_sub_names,
    gene_sub_decent = gene_sub_decent,
    obs_mean = colMeans(X_sub),
    ercc_expected = ercc_expected
)
saveRDS(results, "results_v11_klein_decent.rds")
cat("\nSaved results_v11_klein_decent.rds\n")

# ----------------------------------------------------------------------
# 7. Quick summary
# ----------------------------------------------------------------------
if (!is.null(decent_fit)) {
    cat("\n=== DECENT vs spike-fit oracle summary ===\n")
    # DECENT estimates: decent_fit$em.fit$par.DE has the gene-level estimates
    # For fitNoDE, look at par.noDE structure. Print head.
    cat("  decent_fit class:", class(decent_fit), "\n")
    cat("  names(decent_fit):", paste(names(decent_fit), collapse = ", "), "\n")
} else {
    cat("\n=== DECENT did not return a fit (see error above) ===\n")
}
cat("\nDone.\n")
