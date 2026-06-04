# V2 driver: Poisson-Beta two-state DGP at realistic scRNA-seq scale.
# Tests whether the misspecification result (cell-total consistency,
# per-gene mu bias) is robust to a richer (and more realistic) DGP.

source("sim.R")

set.seed(20260508)

# Realistic scale: 2000 cells, 50 genes from log-normal mu distribution
n_cells <- 2000
n_genes <- 50
mu_vec  <- realistic_mu(n_genes, log_mean = 0.5, log_sd = 1.3, seed = 1)

cat("True mu summary (50 genes from log-normal):\n")
print(summary(mu_vec))
cat("Quartiles:", round(quantile(mu_vec, c(0.25, 0.5, 0.75)), 3), "\n\n")

n_reps   <- 30
alpha_0  <- 1.5
alpha_1_grid <- c(0.0, 0.3, 0.6, 1.0)

one_replicate_pb <- function(n_cells, mu_vec, alpha_0, alpha_1,
                             a = 0.5, burst_factor = 5,
                             spike_in_frac = 0.10) {
    sim <- sim_scrna_pb(n_cells, mu_vec, a, burst_factor,
                        alpha_0, alpha_1, dropout_type = "logistic")
    naive <- fit_naive_zinb(sim$X)

    n_spike <- ceiling(spike_in_frac * n_cells)
    spike_idx <- sample.int(n_cells, n_spike)
    Y_spike <- sim$Y[spike_idx, , drop = FALSE]
    X_spike <- sim$X[spike_idx, , drop = FALSE]
    drop_pars <- fit_oracle_dropout(Y_spike, X_spike)
    if (any(is.na(drop_pars))) {
        oracle <- matrix(NA_real_, length(mu_vec), 2,
                         dimnames = list(NULL, c("mu_hat", "phi_hat")))
    } else {
        oracle <- fit_oracle_zinb(sim$X, drop_pars[1], drop_pars[2])
    }

    list(
        truth_mu = mu_vec,
        naive_mu = naive[, "mu_hat"],
        naive_pi = naive[, "pi_hat"],
        oracle_mu = oracle[, "mu_hat"],
        observed_mean = colMeans(sim$X),
        naive_predicted_mean = (1 - naive[, "pi_hat"]) * naive[, "mu_hat"],
        drop_pars = drop_pars
    )
}

results <- list()
for (a1 in alpha_1_grid) {
    cat("---- alpha_1 =", a1, "----\n")
    t0 <- Sys.time()
    reps <- replicate(n_reps,
                      one_replicate_pb(n_cells, mu_vec,
                                       alpha_0 = alpha_0, alpha_1 = a1),
                      simplify = FALSE)
    cat("  ", n_reps, "reps in", round(as.numeric(Sys.time() - t0), 1), "s\n")

    naive_mu_mat   <- t(sapply(reps, `[[`, "naive_mu"))
    oracle_mu_mat  <- t(sapply(reps, `[[`, "oracle_mu"))
    obs_mean_mat   <- t(sapply(reps, `[[`, "observed_mean"))
    naive_pred_mat <- t(sapply(reps, `[[`, "naive_predicted_mean"))

    naive_mu_mean  <- colMeans(naive_mu_mat,  na.rm = TRUE)
    oracle_mu_mean <- colMeans(oracle_mu_mat, na.rm = TRUE)
    obs_mean       <- colMeans(obs_mean_mat)
    naive_pred_mean <- colMeans(naive_pred_mat, na.rm = TRUE)

    naive_mu_bias  <- naive_mu_mean - mu_vec
    oracle_mu_bias <- oracle_mu_mean - mu_vec
    naive_mu_relbias  <- naive_mu_bias / mu_vec
    oracle_mu_relbias <- oracle_mu_bias / mu_vec
    obs_pred_diff <- naive_pred_mean - obs_mean

    quartiles <- function(v) round(quantile(v, c(0.05, 0.25, 0.5, 0.75, 0.95),
                                            na.rm = TRUE), 3)
    cat("  naive  mu rel. bias quantiles (5/25/50/75/95):", quartiles(naive_mu_relbias), "\n")
    cat("  oracle mu rel. bias quantiles (5/25/50/75/95):", quartiles(oracle_mu_relbias), "\n")
    cat("  cell-total observed:", round(sum(obs_mean), 3), "\n")
    cat("  cell-total naive E[X]:", round(sum(naive_pred_mean, na.rm = TRUE), 3),
        " diff =", signif(sum(obs_pred_diff, na.rm = TRUE), 3), "\n")
    cat("  cell-total true mu:", round(sum(mu_vec), 3),
        "  oracle mu sum:", round(sum(oracle_mu_mean, na.rm = TRUE), 3), "\n")

    # Stratify rel-bias by mu quintile
    mu_q <- cut(mu_vec, breaks = quantile(mu_vec, seq(0, 1, 0.2)),
                include.lowest = TRUE, labels = paste0("Q", 1:5))
    cat("  naive rel. bias by mu quintile:\n")
    by_quint <- tapply(naive_mu_relbias, mu_q, mean, na.rm = TRUE)
    print(round(by_quint, 3))

    results[[as.character(a1)]] <- list(
        alpha_1 = a1,
        mu_vec = mu_vec,
        naive_mu_mean = naive_mu_mean,
        oracle_mu_mean = oracle_mu_mean,
        obs_mean = obs_mean,
        naive_pred_mean = naive_pred_mean,
        naive_mu_relbias = naive_mu_relbias,
        oracle_mu_relbias = oracle_mu_relbias,
        cell_total_obs = sum(obs_mean),
        cell_total_naive = sum(naive_pred_mean, na.rm = TRUE),
        cell_total_oracle = sum(oracle_mu_mean, na.rm = TRUE),
        cell_total_truth = sum(mu_vec),
        by_quint = by_quint
    )
}

saveRDS(results, "results_pb.rds")
cat("\nDone. Saved results_pb.rds\n")
