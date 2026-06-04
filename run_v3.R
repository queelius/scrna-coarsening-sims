# V3 driver: dropout-form robustness (v3a) and oracle misspecification (v3b).
# Question 1 (v3a): does cell-total consistency hold under exponential dropout?
# Question 2 (v3b): how bad is per-gene mu recovery when oracle assumes the
#   wrong dropout family?

source("sim.R")

set.seed(20260509)

n_cells <- 2000
n_genes <- 50
mu_vec  <- realistic_mu(n_genes, log_mean = 0.5, log_sd = 1.3, seed = 1)

n_reps <- 25
spike_in_frac <- 0.10

# Severity grids per form. Picked so that the dropout rate at median mu (~2)
# spans roughly 0.50 to 0.90 in both forms.
sev_logistic <- c(0.0, 0.3, 0.6, 1.0)   # alpha_1 in pi_logistic
sev_exp      <- c(0.0, 0.05, 0.15, 0.30)  # alpha_1 in pi_exponential
truth_intercept <- list(logistic = 1.5, exponential = 0.7)  # alpha_0 / pi0

one_replicate_v3 <- function(truth_form, alpha_1) {
    if (truth_form == "logistic") {
        params_truth <- list(alpha_0 = truth_intercept$logistic, alpha_1 = alpha_1)
        sim <- sim_pb_with_dropout(n_cells, mu_vec, pi_logistic, params_truth)
    } else {
        params_truth <- list(pi0 = truth_intercept$exponential, alpha_1 = alpha_1)
        sim <- sim_pb_with_dropout(n_cells, mu_vec, pi_exponential, params_truth)
    }

    naive <- fit_naive_zinb(sim$X)

    n_spike <- ceiling(spike_in_frac * n_cells)
    spike_idx <- sample.int(n_cells, n_spike)
    Y_spike <- sim$Y[spike_idx, , drop = FALSE]
    X_spike <- sim$X[spike_idx, , drop = FALSE]

    if (truth_form == "logistic") {
        drop_matched    <- fit_dropout_logistic(Y_spike, X_spike)
        drop_mismatched <- fit_dropout_exponential(Y_spike, X_spike)
    } else {
        drop_matched    <- fit_dropout_exponential(Y_spike, X_spike)
        drop_mismatched <- fit_dropout_logistic(Y_spike, X_spike)
    }

    safe_oracle <- function(fit) {
        if (is.null(fit)) {
            matrix(NA_real_, n_genes, 2,
                   dimnames = list(NULL, c("mu_hat", "phi_hat")))
        } else {
            fit_oracle_zinb_v3(sim$X, fit)
        }
    }

    oracle_matched    <- safe_oracle(drop_matched)
    oracle_mismatched <- safe_oracle(drop_mismatched)

    list(
        truth_form = truth_form, alpha_1 = alpha_1,
        truth_mu = mu_vec,
        naive_mu = naive[, "mu_hat"],
        naive_pi = naive[, "pi_hat"],
        oracle_matched_mu    = oracle_matched[, "mu_hat"],
        oracle_mismatched_mu = oracle_mismatched[, "mu_hat"],
        observed_mean        = colMeans(sim$X),
        naive_predicted_mean = (1 - naive[, "pi_hat"]) * naive[, "mu_hat"]
    )
}

# Build full config grid
configs <- rbind(
    data.frame(truth_form = "logistic",    alpha_1 = sev_logistic, stringsAsFactors = FALSE),
    data.frame(truth_form = "exponential", alpha_1 = sev_exp,      stringsAsFactors = FALSE)
)

results <- list()
for (i in seq_len(nrow(configs))) {
    tf <- configs$truth_form[i]; a1 <- configs$alpha_1[i]
    cat(sprintf("---- truth = %-11s  alpha_1 = %-6g ----\n", tf, a1))
    t0 <- Sys.time()
    reps <- replicate(n_reps, one_replicate_v3(tf, a1), simplify = FALSE)
    cat("  ", n_reps, "reps in", round(as.numeric(Sys.time() - t0), 1), "s\n")

    naive_mu_mat <- t(sapply(reps, `[[`, "naive_mu"))
    om_mu_mat    <- t(sapply(reps, `[[`, "oracle_matched_mu"))
    omm_mu_mat   <- t(sapply(reps, `[[`, "oracle_mismatched_mu"))
    obs_mean_mat <- t(sapply(reps, `[[`, "observed_mean"))
    naive_pred_mat <- t(sapply(reps, `[[`, "naive_predicted_mean"))

    naive_mu_mean <- colMeans(naive_mu_mat,  na.rm = TRUE)
    om_mu_mean    <- colMeans(om_mu_mat,     na.rm = TRUE)
    omm_mu_mean   <- colMeans(omm_mu_mat,    na.rm = TRUE)
    obs_mean      <- colMeans(obs_mean_mat)
    naive_pred    <- colMeans(naive_pred_mat, na.rm = TRUE)

    naive_relbias <- (naive_mu_mean - mu_vec) / mu_vec
    om_relbias    <- (om_mu_mean    - mu_vec) / mu_vec
    omm_relbias   <- (omm_mu_mean   - mu_vec) / mu_vec

    qf <- function(v) round(quantile(v, c(0.05, 0.5, 0.95), na.rm = TRUE), 3)
    cat("  naive       relbias (5/50/95):", qf(naive_relbias), "\n")
    cat("  oracle MATCH    (5/50/95):", qf(om_relbias), "\n")
    cat("  oracle MISMATCH (5/50/95):", qf(omm_relbias), "\n")
    cat(sprintf("  cell-total: obs = %.3f, naive pred = %.3f, diff = %.2e\n",
                sum(obs_mean), sum(naive_pred), sum(naive_pred) - sum(obs_mean)))
    cat(sprintf("  cell-total mu sums: truth %.2f, oracle MATCH %.2f, oracle MISMATCH %.2f\n",
                sum(mu_vec), sum(om_mu_mean), sum(omm_mu_mean)))

    results[[paste(tf, a1, sep = "_")]] <- list(
        truth_form = tf, alpha_1 = a1, mu_vec = mu_vec,
        naive_mu_mean = naive_mu_mean,
        oracle_matched_mu_mean = om_mu_mean,
        oracle_mismatched_mu_mean = omm_mu_mean,
        obs_mean = obs_mean,
        naive_pred = naive_pred,
        naive_relbias = naive_relbias,
        oracle_matched_relbias = om_relbias,
        oracle_mismatched_relbias = omm_relbias,
        cell_total_obs = sum(obs_mean),
        cell_total_naive = sum(naive_pred),
        cell_total_oracle_matched = sum(om_mu_mean),
        cell_total_oracle_mismatched = sum(omm_mu_mean),
        cell_total_truth = sum(mu_vec)
    )
}

saveRDS(results, "results_v3.rds")
cat("\nDone. Saved results_v3.rds\n")
