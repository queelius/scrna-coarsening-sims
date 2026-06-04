# Driver: run the misspecification study over multiple replicates and
# multiple C2 violation severities.

source("sim.R")

set.seed(20260507)

mu_vec   <- c(0.5, 1, 2, 5, 10)   # 5 genes spanning low / mid / high expression
phi      <- 1
n_cells  <- 800
n_reps   <- 80
alpha_0  <- 2.0                   # baseline dropout intercept
alpha_1_grid <- c(0.0, 0.3, 0.6, 1.0)  # severity of C2 violation

results <- list()
for (a1 in alpha_1_grid) {
    cat("\nC2 violation severity alpha_1 =", a1, "\n")
    reps <- replicate(n_reps,
                      one_replicate(n_cells, mu_vec, phi,
                                    alpha_0 = alpha_0, alpha_1 = a1,
                                    spike_in_frac = 0.10),
                      simplify = FALSE)

    naive_mu_mat   <- t(sapply(reps, `[[`, "naive_mu"))
    oracle_mu_mat  <- t(sapply(reps, `[[`, "oracle_mu"))
    obs_mean_mat   <- t(sapply(reps, `[[`, "observed_mean"))
    naive_pred_mat <- t(sapply(reps, `[[`, "naive_predicted_mean"))

    summary_tbl <- data.frame(
        gene             = seq_along(mu_vec),
        true_mu          = mu_vec,
        naive_mu_mean    = colMeans(naive_mu_mat,  na.rm = TRUE),
        naive_mu_bias    = colMeans(naive_mu_mat,  na.rm = TRUE) - mu_vec,
        oracle_mu_mean   = colMeans(oracle_mu_mat, na.rm = TRUE),
        oracle_mu_bias   = colMeans(oracle_mu_mat, na.rm = TRUE) - mu_vec,
        obs_mean         = colMeans(obs_mean_mat),
        naive_pred_mean  = colMeans(naive_pred_mat, na.rm = TRUE),
        obs_pred_diff    = colMeans(naive_pred_mat, na.rm = TRUE) - colMeans(obs_mean_mat)
    )
    cat("  per-gene summary (means across", n_reps, "reps):\n")
    print(summary_tbl, digits = 3, row.names = FALSE)

    cell_total_obs    <- sum(colMeans(obs_mean_mat))
    cell_total_naive  <- sum(colMeans(naive_pred_mat, na.rm = TRUE))
    cell_total_oracle <- sum(colMeans(oracle_mu_mat, na.rm = TRUE))
    cat(sprintf("  cell-total observed mean    : %.3f\n", cell_total_obs))
    cat(sprintf("  cell-total naive E[X] sum   : %.3f   (target: observed mean)\n", cell_total_naive))
    cat(sprintf("  cell-total true mu sum      : %.3f   (target: oracle mu sum)\n", sum(mu_vec)))
    cat(sprintf("  cell-total oracle mu sum    : %.3f\n", cell_total_oracle))

    results[[as.character(a1)]] <- list(
        alpha_1 = a1,
        per_gene = summary_tbl,
        cell_total_obs = cell_total_obs,
        cell_total_naive = cell_total_naive,
        cell_total_oracle = cell_total_oracle,
        cell_total_truth_mu = sum(mu_vec)
    )
}

saveRDS(results, "results.rds")
cat("\nDone. Saved results.rds\n")
