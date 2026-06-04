# V10: Cross-tissue Tabula Muris analysis (Marrow + Liver) mirroring
# the v8 spleen pipeline. Addresses methodology-auditor finding M-A3
# (cross-tissue generalization) by running the same naive-vs-oracle
# comparison on two additional FACS Smart-seq2 tissues.
#
# Pipeline (per tissue):
#   1. Load tissue-counts.csv from local FACS/ unzip
#   2. Cell QC: >= 500 genes detected, >= 50,000 reads
#   3. Gene filter: expressed in >= 10 cells
#   4. ERCC concentration matching (Thermo Fisher cms_095046.txt)
#   5. Fit ERCC dropout function (logistic-in-log-y)
#   6. Naive ZINB + ERCC oracle on top 200 by mean + housekeeping
#   7. Cell-total consistency check
#   8. Oracle/naive ratio quantiles
#   9. Housekeeping-gene comparison table
#  10. Save to results_v10_<tissue>.rds

suppressPackageStartupMessages({
    library(data.table)
})
source("sim.R")

run_one_tissue <- function(tissue, counts_path,
                           ercc_table_path = "/tmp/ercc.tsv",
                           seed = 20260523) {
    cat("\n=================================================\n")
    cat("  Tabula Muris", tissue, "\n")
    cat("=================================================\n")
    set.seed(seed)

    # 1. Load
    cat("Loading", counts_path, "...\n")
    t0 <- Sys.time()
    counts <- fread(counts_path, data.table = FALSE)
    gene_names <- counts[[1]]
    counts <- as.matrix(counts[, -1, drop = FALSE])
    rownames(counts) <- gene_names
    cat("  ", nrow(counts), "rows x", ncol(counts), "cells in",
        round(as.numeric(Sys.time() - t0), 1), "s\n")

    ercc_mask <- grepl("^ERCC-", gene_names)
    cat("  ERCC rows:", sum(ercc_mask), "\n")
    X_ercc  <- t(counts[ercc_mask,  , drop = FALSE])
    X_genes <- t(counts[!ercc_mask, , drop = FALSE])

    # 2. Cell QC
    cell_total_reads     <- rowSums(X_genes)
    cell_genes_detected  <- rowSums(X_genes > 0)
    qc_pass <- cell_total_reads >= 50000 & cell_genes_detected >= 500
    cat("  cells passing QC:", sum(qc_pass), "/", length(qc_pass), "\n")
    X_genes <- X_genes[qc_pass, , drop = FALSE]
    X_ercc  <- X_ercc[qc_pass,  , drop = FALSE]
    n_cells <- nrow(X_genes)

    # 3. Gene filter
    gene_n_cells <- colSums(X_genes > 0)
    keep_genes   <- gene_n_cells >= 10
    cat("  genes expressed in >=10 cells:", sum(keep_genes), "\n")
    X_genes <- X_genes[, keep_genes, drop = FALSE]

    # 4. ERCC concentrations
    ercc_table <- fread(ercc_table_path, data.table = FALSE)
    ercc_match <- match(colnames(X_ercc), ercc_table$`ERCC ID`)
    keep_ercc  <- !is.na(ercc_match)
    X_ercc     <- X_ercc[, keep_ercc, drop = FALSE]
    ercc_match <- ercc_match[keep_ercc]
    ercc_conc  <- ercc_table$`concentration in Mix 1 (attomoles/ul)`[ercc_match]
    cat("  ERCC matched:", ncol(X_ercc), " conc range [",
        signif(min(ercc_conc), 2), ",", signif(max(ercc_conc), 2), "]\n")

    # 5. ERCC dropout fit
    ercc_obs_mean <- colMeans(X_ercc)
    nonzero <- ercc_obs_mean > 0.5
    scale_factor <- if (sum(nonzero) > 5) {
        median(ercc_obs_mean[nonzero] / ercc_conc[nonzero])
    } else 1
    ercc_expected <- ercc_conc * scale_factor
    cat("\nFitting ERCC dropout function...\n")
    t0 <- Sys.time()
    drop_fit <- fit_dropout_ercc_logistic(ercc_expected, X_ercc)
    cat("  fit time:", round(as.numeric(Sys.time() - t0), 1), "s\n")
    cat("  alpha_0:", round(drop_fit$params$alpha_0, 3),
        "  alpha_1:", round(drop_fit$params$alpha_1, 3), "\n")
    cat("  pi(y=0)=",   round(plogis(drop_fit$params$alpha_0), 3),
        " pi(y=10)=",   round(plogis(drop_fit$params$alpha_0 -
                                     drop_fit$params$alpha_1 * log(11)), 3),
        " pi(y=100)=",  round(plogis(drop_fit$params$alpha_0 -
                                     drop_fit$params$alpha_1 * log(101)), 3), "\n")

    # 6. Pick gene subset
    gene_means_obs <- colMeans(X_genes)
    top_by_mean    <- order(gene_means_obs, decreasing = TRUE)[1:200]
    hk_genes <- c("Gapdh", "Actb", "B2m", "Hprt", "Tbp", "Ywhaz", "Hmbs",
                  "Tfrc", "Sdha", "Ubc", "Pgk1", "Rpl13a", "Hprt1")
    hk_present <- intersect(hk_genes, colnames(X_genes))
    cat("\nHousekeeping genes present:", paste(hk_present, collapse = ", "), "\n")
    hk_idx <- which(colnames(X_genes) %in% hk_present)
    gene_subset_idx <- unique(c(top_by_mean, hk_idx))
    cat("Gene subset:", length(gene_subset_idx), "genes\n")
    X_sub <- X_genes[, gene_subset_idx, drop = FALSE]
    gene_subset_names <- colnames(X_genes)[gene_subset_idx]

    # 7. Fits
    cat("\nFitting naive ZINB...\n")
    t0 <- Sys.time()
    naive <- fit_naive_zinb(X_sub)
    cat("  ", round(as.numeric(Sys.time() - t0), 1), "s,",
        sum(is.na(naive[, "mu_hat"])), "failures\n")
    cat("Fitting ERCC oracle...\n")
    t0 <- Sys.time()
    oracle <- fit_oracle_zinb_v3(X_sub, drop_fit)
    cat("  ", round(as.numeric(Sys.time() - t0), 1), "s,",
        sum(is.na(oracle[, "mu_hat"])), "failures\n")

    naive_mu  <- naive[, "mu_hat"]
    oracle_mu <- oracle[, "mu_hat"]
    naive_pi  <- naive[, "pi_hat"]
    obs_mean  <- colMeans(X_sub)

    # 8. Cell-total consistency
    naive_pred_mean <- (1 - naive_pi) * naive_mu
    diff_abs <- abs(naive_pred_mean - obs_mean)
    cat("\n=== Cell-total consistency ===\n")
    cat("  obs mean (q5/50/95):",
        round(quantile(obs_mean, c(0.05, 0.5, 0.95)), 3), "\n")
    cat("  naive (1-pi)*mu    :",
        round(quantile(naive_pred_mean, c(0.05, 0.5, 0.95)), 3), "\n")
    cat("  median |diff|      :", signif(median(diff_abs, na.rm = TRUE), 3), "\n")
    cat("  max |diff|         :", signif(max(diff_abs, na.rm = TRUE), 3), "\n")
    cat("  n genes            :", length(naive_mu), "\n")

    # 9. Oracle/naive ratios
    ratio <- oracle_mu / naive_mu
    cat("\n=== Oracle / Naive ratios ===\n")
    cat("  q5/25/50/75/95:",
        round(quantile(ratio, c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE), 2), "\n")

    # 10. Housekeeping table
    cat("\n=== Housekeeping genes ===\n")
    hk_in_subset <- which(gene_subset_names %in% hk_present)
    hk_table <- data.frame(
        gene      = gene_subset_names[hk_in_subset],
        obs_mean  = obs_mean[hk_in_subset],
        naive_mu  = naive_mu[hk_in_subset],
        oracle_mu = oracle_mu[hk_in_subset],
        ratio     = oracle_mu[hk_in_subset] / naive_mu[hk_in_subset]
    )
    print(hk_table, row.names = FALSE)

    # Bias bound diagnostic: implied logit shift
    # If oracle is unbiased, oracle/observed should be 1/(1-pi_endo).
    # The naive pi (from ZINB fit) tells us empirical dropout share at obs_mean.
    # ERCC-fit pi at obs_mean tells us what ERCC predicts at the same mean.
    # Shift = logit(ERCC pi) - logit(naive pi) at obs_mean.
    cat("\n=== Implied ERCC-endogenous shift (housekeeping) ===\n")
    hk_obs <- obs_mean[hk_in_subset]
    hk_naive_pi <- naive_pi[hk_in_subset]
    hk_ercc_pi  <- plogis(drop_fit$params$alpha_0 -
                          drop_fit$params$alpha_1 * log(hk_obs + 1))
    safe_logit <- function(p) {
        p_c <- pmin(pmax(p, 0.001), 0.999)
        log(p_c / (1 - p_c))
    }
    implied_s <- safe_logit(hk_ercc_pi) - safe_logit(hk_naive_pi)
    cat("  ERCC pi  :", round(hk_ercc_pi, 3), "\n")
    cat("  naive pi :", round(hk_naive_pi, 3), "\n")
    cat("  shift s  :", round(implied_s, 3), "\n")
    cat("  median |s| (housekeeping):",
        round(median(abs(implied_s), na.rm = TRUE), 3), "\n")

    saveRDS(list(
        tissue            = tissue,
        drop_fit          = drop_fit,
        naive             = naive,
        oracle            = oracle,
        obs_mean          = obs_mean,
        gene_subset_names = gene_subset_names,
        hk_genes_present  = hk_present,
        hk_table          = hk_table,
        cell_total_median_abs_diff = median(diff_abs, na.rm = TRUE),
        cell_total_max_abs_diff    = max(diff_abs, na.rm = TRUE),
        ratio_quantiles  = round(quantile(ratio, c(0.05, 0.25, 0.5, 0.75, 0.95),
                                          na.rm = TRUE), 3),
        implied_shift_median = round(median(abs(implied_s), na.rm = TRUE), 3),
        n_cells           = n_cells,
        n_ercc            = ncol(X_ercc),
        n_genes_total     = ncol(X_genes)
    ), sprintf("results_v10_%s.rds", tolower(tissue)))
    cat("\nSaved results_v10_", tolower(tissue), ".rds\n", sep = "")
}

# Run both tissues
run_one_tissue("Marrow", "tabula_muris/FACS/Marrow-counts.csv")
run_one_tissue("Liver",  "tabula_muris/FACS/Liver-counts.csv")
cat("\n=== ALL TISSUES DONE ===\n")
