# V4-monotone: bias non-monotonicity at intermediate alpha_1.
# Question: is the v3e finding (bias goes -0.013, +0.015, -0.069 across
# alpha_1 = 0, 0.6, 1.0) MC noise or a real non-monotone curve?
# Run finer grid with more reps to resolve.

source("sim.R")
set.seed(20260509)

n_cells_per_cond <- 1000
n_genes <- 50
n_de <- 10
spike_in_frac <- 0.10
de_factor <- 2; true_logfc_de <- log(de_factor)

mu_A <- realistic_mu(n_genes, log_mean = 0.5, log_sd = 1.3, seed = 1)
mu_B <- mu_A
de_idx <- sample.int(n_genes, n_de)
mu_B[de_idx] <- mu_B[de_idx] * de_factor

alpha_1_grid <- c(0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.3)
n_reps <- 80
alpha_0 <- 1.5

one_rep <- function(alpha_1) {
    sim_A <- sim_pb_with_dropout(n_cells_per_cond, mu_A, pi_logistic,
                                 list(alpha_0 = alpha_0, alpha_1 = alpha_1))
    sim_B <- sim_pb_with_dropout(n_cells_per_cond, mu_B, pi_logistic,
                                 list(alpha_0 = alpha_0, alpha_1 = alpha_1))
    naive_A <- fit_naive_zinb(sim_A$X)
    naive_B <- fit_naive_zinb(sim_B$X)
    naive_logfc <- log(naive_B[, "mu_hat"] / naive_A[, "mu_hat"])
    n_spike <- ceiling(spike_in_frac * n_cells_per_cond)
    Y_pool <- rbind(sim_A$Y[sample.int(n_cells_per_cond, n_spike), ],
                    sim_B$Y[sample.int(n_cells_per_cond, n_spike), ])
    X_pool <- rbind(sim_A$X[sample.int(n_cells_per_cond, n_spike), ],
                    sim_B$X[sample.int(n_cells_per_cond, n_spike), ])
    drop_fit <- fit_dropout_logistic(Y_pool, X_pool)
    if (is.null(drop_fit)) {
        oracle_logfc <- rep(NA_real_, n_genes)
    } else {
        oracle_A <- fit_oracle_zinb_v3(sim_A$X, drop_fit)
        oracle_B <- fit_oracle_zinb_v3(sim_B$X, drop_fit)
        oracle_logfc <- log(oracle_B[, "mu_hat"] / oracle_A[, "mu_hat"])
    }
    list(naive_logfc = naive_logfc, oracle_logfc = oracle_logfc)
}

results <- list()
for (a1 in alpha_1_grid) {
    cat(sprintf("---- alpha_1 = %g ----\n", a1))
    t0 <- Sys.time()
    reps <- replicate(n_reps, one_rep(a1), simplify = FALSE)
    cat("  ", n_reps, "reps in",
        round(as.numeric(Sys.time() - t0, units = "secs"), 1), "s\n")
    n_log <- t(sapply(reps, `[[`, "naive_logfc"))
    o_log <- t(sapply(reps, `[[`, "oracle_logfc"))
    nm <- colMeans(n_log, na.rm = TRUE); ns <- apply(n_log, 2, sd, na.rm = TRUE)
    om <- colMeans(o_log, na.rm = TRUE); os <- apply(o_log, 2, sd, na.rm = TRUE)
    is_de <- seq_len(n_genes) %in% de_idx
    de_naive_bias <- mean(nm[is_de]) - true_logfc_de
    de_oracle_bias <- mean(om[is_de]) - true_logfc_de
    # Per-gene SE of bias estimate (across n_reps)
    se_naive_bias <- sqrt(mean((ns[is_de])^2) / n_reps / sum(is_de))
    cat(sprintf("  DE bias: naive = %+.4f (SE %.4f), oracle = %+.4f\n",
                de_naive_bias, se_naive_bias, de_oracle_bias))
    cat(sprintf("  Non-DE: naive mean abs = %.3f max = %.3f, oracle mean abs = %.3f max = %.3f\n",
                mean(abs(nm[!is_de])), max(abs(nm[!is_de])),
                mean(abs(om[!is_de])), max(abs(om[!is_de]))))
    results[[sprintf("a%.2f", a1)]] <- list(
        alpha_1 = a1, naive_logfc = n_log, oracle_logfc = o_log,
        is_de = is_de, naive_mean = nm, oracle_mean = om,
        naive_sd = ns, oracle_sd = os,
        de_naive_bias = de_naive_bias, de_oracle_bias = de_oracle_bias,
        se_naive_bias = se_naive_bias
    )
}
saveRDS(results, "results_v4_monotone.rds")
cat("Done.\n")
