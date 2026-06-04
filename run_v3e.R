# V3e: differential-expression testing impact.
# Two-condition simulation. Question: does naive ZINB produce biased
# log-fold-change estimates? Does oracle correction restore unbiased
# logFC?
#
# Design:
#   1000 cells per condition (A, B), 50 genes total
#   40 genes are non-DE: mu_A_j = mu_B_j
#   10 genes are DE: mu_B_j = 2 * mu_A_j  (true logFC = log(2) = 0.693)
#   Both conditions experience same logistic dropout

source("sim.R")
set.seed(20260509)

n_cells_per_cond <- 1000
n_genes <- 50
n_de <- 10
de_factor <- 2
true_logfc_de <- log(de_factor)

mu_A <- realistic_mu(n_genes, log_mean = 0.5, log_sd = 1.3, seed = 1)
mu_B <- mu_A
de_idx <- sample.int(n_genes, n_de)
mu_B[de_idx] <- mu_B[de_idx] * de_factor
true_logfc <- log(mu_B / mu_A)
cat("True DE genes:", sort(de_idx), "\n")
cat("True logFC: nonDE all 0;", n_de, "DE at log(", de_factor, ") =",
    round(true_logfc_de, 3), "\n\n")

alpha_0 <- 1.5
alpha_1_grid <- c(0.0, 0.6, 1.0)
n_reps <- 25
spike_in_frac <- 0.10

one_rep_v3e <- function(alpha_1) {
    # Generate condition A
    sim_A <- sim_pb_with_dropout(n_cells_per_cond, mu_A, pi_logistic,
                                 params = list(alpha_0 = alpha_0,
                                               alpha_1 = alpha_1))
    # Generate condition B
    sim_B <- sim_pb_with_dropout(n_cells_per_cond, mu_B, pi_logistic,
                                 params = list(alpha_0 = alpha_0,
                                               alpha_1 = alpha_1))

    # Naive fit per condition
    naive_A <- fit_naive_zinb(sim_A$X)
    naive_B <- fit_naive_zinb(sim_B$X)
    naive_logfc <- log(naive_B[, "mu_hat"] / naive_A[, "mu_hat"])

    # Oracle fit: pool spike-ins from both conditions for dropout fit
    n_spike <- ceiling(spike_in_frac * n_cells_per_cond)
    sa_idx <- sample.int(n_cells_per_cond, n_spike)
    sb_idx <- sample.int(n_cells_per_cond, n_spike)
    Y_pool <- rbind(sim_A$Y[sa_idx, ], sim_B$Y[sb_idx, ])
    X_pool <- rbind(sim_A$X[sa_idx, ], sim_B$X[sb_idx, ])

    drop_fit <- fit_dropout_logistic(Y_pool, X_pool)

    if (is.null(drop_fit)) {
        oracle_logfc <- rep(NA_real_, n_genes)
    } else {
        oracle_A <- fit_oracle_zinb_v3(sim_A$X, drop_fit)
        oracle_B <- fit_oracle_zinb_v3(sim_B$X, drop_fit)
        oracle_logfc <- log(oracle_B[, "mu_hat"] / oracle_A[, "mu_hat"])
    }

    list(alpha_1 = alpha_1,
         naive_logfc = naive_logfc,
         oracle_logfc = oracle_logfc,
         naive_mu_A = naive_A[, "mu_hat"],
         naive_mu_B = naive_B[, "mu_hat"])
}

results <- list()
for (a1 in alpha_1_grid) {
    cat(sprintf("---- alpha_1 = %g ----\n", a1))
    t0 <- Sys.time()
    reps <- replicate(n_reps, one_rep_v3e(a1), simplify = FALSE)
    cat("  ", n_reps, "reps in", round(as.numeric(Sys.time() - t0), 1), "s\n")

    naive_logfc_mat  <- t(sapply(reps, `[[`, "naive_logfc"))
    oracle_logfc_mat <- t(sapply(reps, `[[`, "oracle_logfc"))

    naive_mean  <- colMeans(naive_logfc_mat,  na.rm = TRUE)
    oracle_mean <- colMeans(oracle_logfc_mat, na.rm = TRUE)
    naive_sd  <- apply(naive_logfc_mat,  2, sd, na.rm = TRUE)
    oracle_sd <- apply(oracle_logfc_mat, 2, sd, na.rm = TRUE)

    # Bias and power summaries
    is_de <- seq_len(n_genes) %in% de_idx
    cat("  ---- non-DE genes (true logFC = 0) ----\n")
    cat("  naive  logFC: mean =", round(mean(naive_mean[!is_de]),  3),
        " (|max| =", round(max(abs(naive_mean[!is_de])), 3),
        "); SD across reps median =", round(median(naive_sd[!is_de]),  3), "\n")
    cat("  oracle logFC: mean =", round(mean(oracle_mean[!is_de]), 3),
        " (|max| =", round(max(abs(oracle_mean[!is_de])), 3),
        "); SD across reps median =", round(median(oracle_sd[!is_de]), 3), "\n")
    cat("  ---- DE genes (true logFC =", round(true_logfc_de, 3), ") ----\n")
    cat("  naive  logFC: mean =", round(mean(naive_mean[is_de]),  3),
        " bias =", round(mean(naive_mean[is_de]) - true_logfc_de, 3),
        "; SD across reps median =", round(median(naive_sd[is_de]),  3), "\n")
    cat("  oracle logFC: mean =", round(mean(oracle_mean[is_de]), 3),
        " bias =", round(mean(oracle_mean[is_de]) - true_logfc_de, 3),
        "; SD across reps median =", round(median(oracle_sd[is_de]), 3), "\n")

    # Discovery rate at z > 2 threshold (rough Wald-style call)
    naive_z  <- naive_mean  / pmax(naive_sd / sqrt(n_reps),  1e-8)
    oracle_z <- oracle_mean / pmax(oracle_sd / sqrt(n_reps), 1e-8)
    naive_called  <- abs(naive_z)  > 2
    oracle_called <- abs(oracle_z) > 2
    cat("  discovery rates (|z| > 2):\n")
    cat("    naive:  type-I = ", round(mean(naive_called[!is_de]),  3),
        "  power = ", round(mean(naive_called[is_de]),  3), "\n")
    cat("    oracle: type-I = ", round(mean(oracle_called[!is_de]), 3),
        "  power = ", round(mean(oracle_called[is_de]), 3), "\n")

    results[[as.character(a1)]] <- list(
        alpha_1 = a1, is_de = is_de, true_logfc = true_logfc,
        naive_logfc_mean = naive_mean,
        oracle_logfc_mean = oracle_mean,
        naive_logfc_sd = naive_sd,
        oracle_logfc_sd = oracle_sd
    )
}

saveRDS(results, "results_v3e.rds")
cat("\nDone. Saved results_v3e.rds\n")
