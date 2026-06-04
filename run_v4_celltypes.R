# V4-celltypes: cell-type heterogeneity.
# DGP modification: 50/50 mix of two cell types A and B with different mu
# profiles. Naive ZINB and oracle both fit single per-gene mu (the wrong
# model). Question: does dropout correction still help when the data
# violate homogeneity?

source("sim.R")
set.seed(20260509)

n_cells  <- 2000
n_genes  <- 50
spike_in_frac <- 0.10
alpha_0  <- 1.5
alpha_1_grid <- c(0.0, 0.6, 1.0)
n_reps   <- 30

# Two cell types with different mu profiles.
mu_A <- realistic_mu(n_genes, log_mean = 0,   log_sd = 1.0, seed = 11)
mu_B <- realistic_mu(n_genes, log_mean = 1.5, log_sd = 1.0, seed = 22)
mu_pop <- 0.5 * mu_A + 0.5 * mu_B  # population marginal mean

sim_two_celltypes <- function(n_cells, mu_A, mu_B, alpha_0, alpha_1,
                              a = 0.5, burst_factor = 5) {
    n_A <- n_cells %/% 2
    n_B <- n_cells - n_A
    p_genes <- length(mu_A)
    Y <- matrix(0, n_cells, p_genes)
    for (j in seq_len(p_genes)) {
        Y[1:n_A, j] <- sim_pb_gene(n_A, mu_A[j], a, burst_factor)
        Y[(n_A + 1):n_cells, j] <- sim_pb_gene(n_B, mu_B[j], a, burst_factor)
    }
    cell_type <- c(rep("A", n_A), rep("B", n_B))
    X <- apply_dropout_fn(Y, pi_logistic,
                          list(alpha_0 = alpha_0, alpha_1 = alpha_1))
    list(Y = Y, X = X, cell_type = cell_type)
}

one_rep <- function(alpha_1) {
    sim <- sim_two_celltypes(n_cells, mu_A, mu_B, alpha_0, alpha_1)
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
    list(naive_mu = naive[, "mu_hat"], oracle_mu = oracle_mu,
         observed_mean = colMeans(sim$X))
}

results <- list()
for (a1 in alpha_1_grid) {
    cat(sprintf("---- alpha_1 = %g ----\n", a1))
    t0 <- Sys.time()
    reps <- replicate(n_reps, one_rep(a1), simplify = FALSE)
    cat("  ", n_reps, "reps in",
        round(as.numeric(Sys.time() - t0, units = "secs"), 1), "s\n")

    n_mu <- t(sapply(reps, `[[`, "naive_mu"))
    o_mu <- t(sapply(reps, `[[`, "oracle_mu"))

    # Compare to TWO targets: population marginal mu_pop and a "wrong-but-natural" target
    naive_rb_pop <- (colMeans(n_mu, na.rm = TRUE) - mu_pop) / mu_pop
    oracle_rb_pop <- (colMeans(o_mu, na.rm = TRUE) - mu_pop) / mu_pop

    qf <- function(v) round(quantile(v, c(0.05, 0.5, 0.95), na.rm = TRUE), 3)
    cat("  vs population mean (mu_pop):\n")
    cat("    naive  rb (5/50/95):", qf(naive_rb_pop), "\n")
    cat("    oracle rb (5/50/95):", qf(oracle_rb_pop), "\n")
    # Also report the absolute-bias improvement
    abs_imp <- median(abs(naive_rb_pop) - abs(oracle_rb_pop))
    cat(sprintf("  median |naive_rb| - |oracle_rb| = %+.3f (positive = oracle helps)\n",
                abs_imp))

    results[[sprintf("a%.2f", a1)]] <- list(
        alpha_1 = a1, naive_mu = n_mu, oracle_mu = o_mu,
        naive_rb_pop = naive_rb_pop, oracle_rb_pop = oracle_rb_pop,
        mu_pop = mu_pop, mu_A = mu_A, mu_B = mu_B
    )
}
saveRDS(results, "results_v4_celltypes.rds")
cat("Done.\n")
