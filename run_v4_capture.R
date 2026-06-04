# V4-capture: cell-level capture-efficiency variation.
# DGP modification: each cell has a baseline-dropout shift delta_i ~ N(0, sigma_q),
# so effective dropout is plogis(alpha_0 + delta_i - alpha_1*log(y+1)).
# This induces gene-gene correlation in dropout within cells.
# Question: does the standard oracle (assuming homogeneous dropout) still
# work, or does cell-level variation break it?

source("sim.R")
set.seed(20260509)

n_cells  <- 2000
n_genes  <- 50
mu_vec   <- realistic_mu(n_genes, log_mean = 0.5, log_sd = 1.3, seed = 1)
spike_in_frac <- 0.10

alpha_0      <- 1.5
alpha_1      <- 0.6
sigma_q_grid <- c(0.0, 0.25, 0.5, 1.0)   # 0 = homogeneous; 1 = heavy variation
n_reps       <- 30

#' Apply logistic dropout with cell-level intercept shift delta_i.
apply_dropout_celleffect <- function(Y, alpha_0, alpha_1, delta) {
    eta_base <- alpha_0 - alpha_1 * log(Y + 1)
    eta <- sweep(eta_base, 1, delta, "+")
    pi_drop <- pmin(pmax(plogis(eta), 0), 1)
    keep <- matrix(rbinom(length(Y), 1, 1 - pi_drop), nrow(Y), ncol(Y))
    Y * keep
}

sim_with_capture <- function(n_cells, mu_vec, alpha_0, alpha_1, sigma_q,
                             a = 0.5, burst_factor = 5) {
    p_genes <- length(mu_vec)
    Y <- matrix(0, n_cells, p_genes)
    for (j in seq_len(p_genes)) {
        Y[, j] <- sim_pb_gene(n_cells, mu_vec[j], a, burst_factor)
    }
    delta <- rnorm(n_cells, 0, sigma_q)
    X <- apply_dropout_celleffect(Y, alpha_0, alpha_1, delta)
    list(Y = Y, X = X, delta = delta, mu_vec = mu_vec)
}

one_rep <- function(sigma_q) {
    sim <- sim_with_capture(n_cells, mu_vec, alpha_0, alpha_1, sigma_q)
    naive <- fit_naive_zinb(sim$X)
    n_spike <- ceiling(spike_in_frac * n_cells)
    spike_idx <- sample.int(n_cells, n_spike)
    Y_spike <- sim$Y[spike_idx, , drop = FALSE]
    X_spike <- sim$X[spike_idx, , drop = FALSE]
    drop_fit <- fit_dropout_logistic(Y_spike, X_spike)
    if (is.null(drop_fit)) {
        oracle_mu <- rep(NA_real_, n_genes)
    } else {
        oracle_mu <- fit_oracle_zinb_v3(sim$X, drop_fit)[, "mu_hat"]
    }
    list(sigma_q = sigma_q,
         naive_mu = naive[, "mu_hat"], oracle_mu = oracle_mu,
         observed_mean = colMeans(sim$X),
         drop_fit_pars = if (is.null(drop_fit)) c(NA, NA) else
             c(drop_fit$params$alpha_0, drop_fit$params$alpha_1))
}

results <- list()
for (sq in sigma_q_grid) {
    cat(sprintf("---- sigma_q = %g ----\n", sq))
    t0 <- Sys.time()
    reps <- replicate(n_reps, one_rep(sq), simplify = FALSE)
    cat("  ", n_reps, "reps in",
        round(as.numeric(Sys.time() - t0, units = "secs"), 1), "s\n")
    n_mu <- t(sapply(reps, `[[`, "naive_mu"))
    o_mu <- t(sapply(reps, `[[`, "oracle_mu"))
    drop_pars <- t(sapply(reps, `[[`, "drop_fit_pars"))

    naive_rb  <- (colMeans(n_mu, na.rm = TRUE) - mu_vec) / mu_vec
    oracle_rb <- (colMeans(o_mu, na.rm = TRUE) - mu_vec) / mu_vec

    qf <- function(v) round(quantile(v, c(0.05, 0.5, 0.95), na.rm = TRUE), 3)
    cat("  naive  rb (5/50/95):", qf(naive_rb), "\n")
    cat("  oracle rb (5/50/95):", qf(oracle_rb), "\n")
    cat(sprintf("  drop fit alpha_0 mean = %.3f (truth %.2f), alpha_1 mean = %.3f (truth %.2f)\n",
                mean(drop_pars[, 1], na.rm = TRUE), alpha_0,
                mean(drop_pars[, 2], na.rm = TRUE), alpha_1))
    results[[sprintf("sq%.2f", sq)]] <- list(
        sigma_q = sq, naive_mu = n_mu, oracle_mu = o_mu,
        naive_relbias = naive_rb, oracle_relbias = oracle_rb,
        drop_pars = drop_pars
    )
}
saveRDS(results, "results_v4_capture.rds")
cat("Done.\n")
