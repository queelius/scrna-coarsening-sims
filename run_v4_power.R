# V4-power: power at smaller DE signals (1.3x, 1.5x).
# Question: where does naive ZINB lose more power than oracle?

source("sim.R")
set.seed(20260509)

n_cells_per_cond <- 1000
n_genes <- 50
n_de <- 10
spike_in_frac <- 0.10

mu_A <- realistic_mu(n_genes, log_mean = 0.5, log_sd = 1.3, seed = 1)
de_idx <- sample.int(n_genes, n_de)

de_factor_grid <- c(1.3, 1.5, 2.0)
alpha_1_grid <- c(0.0, 0.6, 1.0)
n_reps <- 40
alpha_0 <- 1.5

one_rep <- function(de_factor, alpha_1) {
    mu_B <- mu_A
    mu_B[de_idx] <- mu_B[de_idx] * de_factor
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
for (df in de_factor_grid) {
    for (a1 in alpha_1_grid) {
        true_logfc_de <- log(df)
        cat(sprintf("---- de_factor = %g, alpha_1 = %g ----\n", df, a1))
        t0 <- Sys.time()
        reps <- replicate(n_reps, one_rep(df, a1), simplify = FALSE)
        cat("  ", n_reps, "reps in",
            round(as.numeric(Sys.time() - t0, units = "secs"), 1), "s\n")
        n_log <- t(sapply(reps, `[[`, "naive_logfc"))
        o_log <- t(sapply(reps, `[[`, "oracle_logfc"))
        nm <- colMeans(n_log, na.rm = TRUE); ns <- apply(n_log, 2, sd, na.rm = TRUE)
        om <- colMeans(o_log, na.rm = TRUE); os <- apply(o_log, 2, sd, na.rm = TRUE)
        is_de <- seq_len(n_genes) %in% de_idx
        n_z <- nm / pmax(ns / sqrt(n_reps), 1e-8)
        o_z <- om / pmax(os / sqrt(n_reps), 1e-8)
        cat(sprintf("  DE genes  : naive bias=%.3f sd=%.3f  oracle bias=%.3f sd=%.3f\n",
                    mean(nm[is_de]) - true_logfc_de, median(ns[is_de]),
                    mean(om[is_de]) - true_logfc_de, median(os[is_de])))
        cat(sprintf("  power |z|>2 : naive=%.2f  oracle=%.2f\n",
                    mean(abs(n_z[is_de]) > 2), mean(abs(o_z[is_de]) > 2)))
        cat(sprintf("  type-I |z|>2: naive=%.3f  oracle=%.3f\n",
                    mean(abs(n_z[!is_de]) > 2), mean(abs(o_z[!is_de]) > 2)))
        results[[sprintf("df%.2f_a%.2f", df, a1)]] <- list(
            de_factor = df, alpha_1 = a1, true_logfc_de = true_logfc_de,
            naive_logfc = n_log, oracle_logfc = o_log,
            is_de = is_de, naive_mean = nm, oracle_mean = om,
            naive_sd = ns, oracle_sd = os
        )
    }
}
saveRDS(results, "results_v4_power.rds")
cat("Done.\n")
