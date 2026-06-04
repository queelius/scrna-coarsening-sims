# V9: Direct comparison with DECENT on Tabula Muris.
# DECENT was designed for UMI scRNA-seq (low counts). Tabula Muris is
# Smart-seq2 with high counts. To make a fair comparison, we restrict
# to genes with observed mean below DECENT's tested regime (~mean < 200)
# and downsample cells.

suppressPackageStartupMessages({
    library(DECENT); library(data.table)
})
source("sim.R")
set.seed(20260512)

# -----------------------------------------------------------------
# Load and preprocess Tabula Muris spleen (same as v8)
# -----------------------------------------------------------------
counts <- fread("tabula_muris/FACS/Spleen-counts.csv", data.table = FALSE)
gn <- counts[[1]]; counts <- as.matrix(counts[, -1]); rownames(counts) <- gn
X_ercc_full <- t(counts[grepl("^ERCC-", gn), , drop = FALSE])
X_genes_full <- t(counts[!grepl("^ERCC-", gn), , drop = FALSE])
qc <- rowSums(X_genes_full) >= 50000 & rowSums(X_genes_full > 0) >= 500
X_ercc <- X_ercc_full[qc, ]; X_genes <- X_genes_full[qc, ]
X_genes <- X_genes[, colSums(X_genes > 0) >= 10]

ercc_table <- fread("/tmp/ercc.tsv", data.table = FALSE)
em <- match(colnames(X_ercc), ercc_table$`ERCC ID`)
X_ercc <- X_ercc[, !is.na(em)]
ercc_conc <- ercc_table$`concentration in Mix 1 (attomoles/ul)`[em[!is.na(em)]]

# Pick genes in DECENT's tested regime (mean expression < 200) plus housekeeping
gene_means <- colMeans(X_genes)
mid_idx <- which(gene_means >= 1 & gene_means <= 200)
hk <- c("Gapdh","Actb","B2m","Hprt","Tbp","Ywhaz","Hmbs","Tfrc","Sdha","Ubc","Pgk1","Rpl13a")
hk_idx <- which(colnames(X_genes) %in% hk)
hk_in_range <- intersect(hk_idx, mid_idx)
# Keep all hk in range plus a sample of other mid-range genes
set.seed(20260512)
other_mid <- setdiff(mid_idx, hk_idx)
other_sample <- sample(other_mid, min(60, length(other_mid)))
gene_idx <- unique(c(hk_in_range, other_sample))
X_sub <- X_genes[, gene_idx]
gene_names_sub <- colnames(X_genes)[gene_idx]
cat("Selected", length(gene_idx), "genes (",
    length(hk_in_range), "housekeeping in range,",
    length(other_sample), "other mid-range)\n")

# Cells: subsample to keep DECENT runtime manageable
cell_idx <- sample.int(nrow(X_sub), 500)
X_sub_t <- t(X_sub[cell_idx, ])  # genes x cells
X_ercc_t <- t(X_ercc[cell_idx, ])

cat("Subset for comparison:", nrow(X_sub_t), "genes x", ncol(X_sub_t), "cells\n")
cat("Gene mean range:", round(quantile(rowMeans(X_sub_t), c(0.05, 0.5, 0.95))), "\n")

# -----------------------------------------------------------------
# Our framework: ERCC oracle
# -----------------------------------------------------------------
cat("\n=== Our framework ===\n")
t0 <- Sys.time()
drop_fit <- fit_dropout_ercc_logistic(ercc_conc, X_ercc[cell_idx, ])
cat("Dropout fit:", round(as.numeric(Sys.time() - t0), 1), "s; ",
    "alpha_0 =", round(drop_fit$params$alpha_0, 3),
    "alpha_1 =", round(drop_fit$params$alpha_1, 3), "\n")

t0 <- Sys.time()
naive <- fit_naive_zinb(t(X_sub_t))  # back to cells x genes for our fitter
oracle <- fit_oracle_zinb_v3(t(X_sub_t), drop_fit)
cat("Naive + oracle fit:", round(as.numeric(Sys.time() - t0), 1), "s\n")

# -----------------------------------------------------------------
# DECENT: same data
# -----------------------------------------------------------------
cat("\n=== DECENT ===\n")
t0 <- Sys.time()
decent_fit <- tryCatch(
    fitNoDE(data.obs = X_sub_t, spikes = X_ercc_t, spike.conc = ercc_conc,
            use.spikes = TRUE, CE.range = c(0.02, 0.1),
            tau.global = TRUE, tau.init = c(-5, 0), tau.est = "endo",
            normalize = "ML", GQ.approx = TRUE, maxit = 10, parallel = FALSE),
    error = function(e) { cat("DECENT failed:", conditionMessage(e), "\n"); NULL })
cat("DECENT fit:", round(as.numeric(Sys.time() - t0), 1), "s\n")

# -----------------------------------------------------------------
# Comparison
# -----------------------------------------------------------------
obs_mean <- rowMeans(X_sub_t)
naive_mu <- naive[, "mu_hat"]
oracle_mu <- oracle[, "mu_hat"]
decent_mu <- if (!is.null(decent_fit)) decent_fit$est.mu else rep(NA_real_, nrow(X_sub_t))

# Save full comparison
df <- data.frame(
    gene = gene_names_sub,
    obs_mean = obs_mean,
    naive_mu = naive_mu,
    oracle_mu = oracle_mu,
    decent_mu = decent_mu,
    is_hk = gene_names_sub %in% hk)

# Print housekeeping comparison
cat("\n=== Housekeeping genes: side-by-side ===\n")
cat(sprintf("  %-10s  %10s  %10s  %10s  %10s\n",
            "gene", "obs_mean", "naive_mu", "oracle_mu", "DECENT_mu"))
hk_rows <- df[df$is_hk, ]
hk_rows <- hk_rows[order(hk_rows$obs_mean, decreasing = TRUE), ]
for (i in seq_len(nrow(hk_rows))) {
    cat(sprintf("  %-10s  %10.2f  %10.2f  %10.2f  %10s\n",
                hk_rows$gene[i], hk_rows$obs_mean[i],
                hk_rows$naive_mu[i], hk_rows$oracle_mu[i],
                ifelse(is.na(hk_rows$decent_mu[i]), "NA",
                       sprintf("%.2f", hk_rows$decent_mu[i]))))
}

# Cross-method correlations
cat("\n=== Cross-method correlations (all genes) ===\n")
cat(sprintf("  cor(naive, oracle) : %.3f\n",
            cor(naive_mu, oracle_mu, use = "pair")))
cat(sprintf("  cor(naive, DECENT) : %.3f\n",
            cor(naive_mu, decent_mu, use = "pair")))
cat(sprintf("  cor(oracle, DECENT): %.3f\n",
            cor(oracle_mu, decent_mu, use = "pair")))

saveRDS(list(df = df, drop_fit = drop_fit, decent_fit = decent_fit),
        "results_v9_decent_compare.rds")
cat("\nDone. Saved results_v9_decent_compare.rds\n")
