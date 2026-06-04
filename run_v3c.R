# V3c: spike-in fraction sweep.
# Where does the oracle break? Sweep spike-in fraction from 0.5% to 20%
# at fixed truth (logistic, alpha_1 = 0.6) and matched-form oracle.

source("sim.R")
set.seed(20260509)

n_cells <- 2000
n_genes <- 50
mu_vec  <- realistic_mu(n_genes, log_mean = 0.5, log_sd = 1.3, seed = 1)

frac_grid <- c(0.005, 0.01, 0.02, 0.05, 0.10, 0.20)
n_reps_at <- function(frac) if (frac < 0.02) 60 else if (frac < 0.05) 40 else 25

alpha_0 <- 1.5; alpha_1 <- 0.6
truth_dropout <- list(form = "logistic", pi_fn = pi_logistic,
                      params = list(alpha_0 = alpha_0, alpha_1 = alpha_1))

one_rep_v3c <- function(frac) {
    sim <- sim_pb_with_dropout(n_cells, mu_vec, pi_logistic,
                               params = list(alpha_0 = alpha_0, alpha_1 = alpha_1))
    n_spike <- max(1, ceiling(frac * n_cells))
    spike_idx <- sample.int(n_cells, n_spike)
    Y_spike <- sim$Y[spike_idx, , drop = FALSE]
    X_spike <- sim$X[spike_idx, , drop = FALSE]

    drop_fit <- fit_dropout_logistic(Y_spike, X_spike)
    if (is.null(drop_fit)) {
        oracle_mu <- rep(NA_real_, n_genes)
    } else {
        oracle_mu <- fit_oracle_zinb_v3(sim$X, drop_fit)[, "mu_hat"]
    }
    list(frac = frac, n_spike = n_spike,
         oracle_mu = oracle_mu,
         drop_fit = if (is.null(drop_fit)) c(NA, NA) else c(drop_fit$params$alpha_0,
                                                              drop_fit$params$alpha_1))
}

results <- list()
for (frac in frac_grid) {
    n_reps <- n_reps_at(frac)
    cat(sprintf("---- frac = %g (n_spike = %d, n_reps = %d) ----\n",
                frac, max(1, ceiling(frac * n_cells)), n_reps))
    t0 <- Sys.time()
    reps <- replicate(n_reps, one_rep_v3c(frac), simplify = FALSE)
    cat("  ", n_reps, "reps in", round(as.numeric(Sys.time() - t0), 1), "s\n")

    om_mu_mat <- t(sapply(reps, `[[`, "oracle_mu"))
    drop_pars_mat <- t(sapply(reps, `[[`, "drop_fit"))

    om_relbias_mean <- (colMeans(om_mu_mat, na.rm = TRUE) - mu_vec) / mu_vec
    om_relbias_sd <- apply(om_mu_mat, 2,
                           function(v) sd(v / mean(v, na.rm = TRUE), na.rm = TRUE))
    failure_rate <- mean(is.na(om_mu_mat))

    cat("  median relbias:", round(median(om_relbias_mean, na.rm = TRUE), 4), "\n")
    cat("  IQR  relbias :",
        round(quantile(om_relbias_mean, 0.25, na.rm = TRUE), 4), "to",
        round(quantile(om_relbias_mean, 0.75, na.rm = TRUE), 4), "\n")
    cat("  rel-bias SD across reps (median):",
        round(median(om_relbias_sd, na.rm = TRUE), 4), "\n")
    cat("  oracle failure rate:", round(failure_rate, 3), "\n")
    cat("  dropout-fit alpha_0 mean:", round(mean(drop_pars_mat[, 1], na.rm = TRUE), 3),
        " (truth", alpha_0, ")\n")
    cat("  dropout-fit alpha_1 mean:", round(mean(drop_pars_mat[, 2], na.rm = TRUE), 3),
        " (truth", alpha_1, ")\n")

    results[[as.character(frac)]] <- list(
        frac = frac, n_reps = n_reps,
        oracle_mu = om_mu_mat,
        drop_pars = drop_pars_mat,
        om_relbias_mean = om_relbias_mean,
        om_relbias_sd = om_relbias_sd,
        failure_rate = failure_rate
    )
}

saveRDS(results, "results_v3c.rds")
cat("\nDone. Saved results_v3c.rds\n")
