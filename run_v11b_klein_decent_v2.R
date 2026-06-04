# V11b: Retry DECENT on Klein 2015 with corrected spike-in scaling.
#
# Diagnosis from v11: DECENT produced NA log-likelihoods on first
# iteration. Examining the DECENT example data showed sp.true on the
# order of 100s (with sp.obs row sums in the 1000s, implying capture
# efficiency around 5 percent). Our prior scaling used median(obs/conc)
# yielding ercc_expected medians around 0.004, far too small for DECENT.
#
# This retry estimates the per-ERCC true count via the inverse-CE
# correction: sp.true_k = sp.obs.rowmean_k / CE, with CE = 0.05 (the
# midpoint of DECENT's documented CE.range). We then pass sp.true on
# the same scale DECENT expects.

suppressPackageStartupMessages({
    library(data.table)
    library(DECENT)
})
source("sim.R")

set.seed(20260523)

# Reload Klein K562 data
cat("Loading Klein 2015 K562_cells.csv...\n")
counts <- fread("klein2015/GSM1599500_K562_cells.csv", data.table = FALSE,
                header = FALSE)
gene_names <- counts[[1]]
counts <- as.matrix(counts[, -1, drop = FALSE])
rownames(counts) <- gene_names
mode(counts) <- "integer"

keep_named <- !is.na(gene_names) & nchar(gene_names) > 0
counts <- counts[keep_named, , drop = FALSE]
gene_names <- gene_names[keep_named]

ercc_mask <- grepl("^ERCC-", gene_names)
X_ercc  <- t(counts[ercc_mask,  , drop = FALSE])
X_genes <- t(counts[!ercc_mask, , drop = FALSE])

cell_total <- rowSums(X_genes)
cell_genes_det <- rowSums(X_genes > 0)
qc_pass <- cell_total >= 1000 & cell_genes_det >= 200
X_genes <- X_genes[qc_pass, , drop = FALSE]
X_ercc  <- X_ercc[qc_pass,  , drop = FALSE]

gene_n_cells <- colSums(X_genes > 0)
keep_genes <- gene_n_cells >= 10
X_genes <- X_genes[, keep_genes, drop = FALSE]

gene_means_obs <- colMeans(X_genes)
n_target <- 500
top_idx <- order(gene_means_obs, decreasing = TRUE)[1:n_target]
X_sub <- X_genes[, top_idx, drop = FALSE]

cat("X_sub:", dim(X_sub), "  X_ercc:", dim(X_ercc), "\n")

# DECENT-style spike-true: inverse-CE scaling
# Use only the ERCC species with at least some observed expression.
# DECENT's example has CE around 5 percent (1/16 ratio approx).
# We estimate per-ERCC sp.true from observed row means assuming CE = 0.05.
CE_assumed <- 0.05
ercc_obs_mean <- colMeans(X_ercc)
# Drop ERCCs that are completely zero (DECENT cannot handle).
nonzero_ercc <- ercc_obs_mean > 0.5
X_ercc_keep <- X_ercc[, nonzero_ercc, drop = FALSE]
sp_true <- colMeans(X_ercc_keep) / CE_assumed
cat("\nERCC with non-trivial signal:", ncol(X_ercc_keep), "\n")
cat("sp.true range:", round(range(sp_true), 1), "\n")
cat("sp.true (first 5):", round(sp_true[1:5], 1), "\n")

# DECENT inputs (genes x cells for both)
data_mat <- t(X_sub)
spikes_mat <- t(X_ercc_keep)
cat("\nDECENT inputs: data", dim(data_mat), "  spikes", dim(spikes_mat), "\n")

cat("\nFitting DECENT fitNoDE with CE-corrected sp.true...\n")
t0 <- Sys.time()
decent_fit <- tryCatch({
    DECENT::fitNoDE(
        data.obs    = data_mat,
        spikes      = spikes_mat,
        spike.conc  = sp_true,
        use.spikes  = TRUE,
        CE.range    = c(0.02, 0.10),
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
cat("  runtime:", elapsed, "s\n")

if (!is.null(decent_fit)) {
    cat("\n=== DECENT SUCCEEDED ===\n")
    cat("  class:", class(decent_fit), "\n")
    cat("  names:", paste(names(decent_fit), collapse = ", "), "\n")
} else {
    cat("\n=== DECENT FAILED on v11b too ===\n")
}

saveRDS(list(decent_fit = decent_fit, sp_true = sp_true,
             elapsed_s = elapsed), "results_v11b_klein_decent.rds")
cat("\nSaved.\n")
