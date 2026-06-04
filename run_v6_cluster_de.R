# V6-cluster-de: per-cluster DE benchmark with explicit non-DE blocks.
# Setup: 35 non-DE genes (true logFC = 0 between cluster A and B),
# plus 15 DE genes split as 5 each at de_factor = 1.5, 2, 3.
# Cleanly assesses type-I and stratified power.

source("sim.R")
set.seed(20260511)

n_cells_per_cluster <- 1000
n_total <- 2 * n_cells_per_cluster
spike_in_frac <- 0.10
alpha_0  <- 1.5
alpha_1_grid <- c(0.0, 0.3, 0.6, 1.0)
n_reps   <- 30

# Gene structure: 35 non-DE + 5*1.5x + 5*2x + 5*3x = 50 total
n_nonde <- 35
de_factors <- c(rep(1.5, 5), rep(2, 5), rep(3, 5))
n_de <- length(de_factors)
n_genes <- n_nonde + n_de

# Common mu profile, then DE genes get the factor applied to one cluster
mu_base <- realistic_mu(n_genes, log_mean = 0.5, log_sd = 1.0, seed = 33)
mu_A <- mu_base
mu_B <- mu_base
de_idx <- (n_nonde + 1):n_genes
mu_B[de_idx] <- mu_B[de_idx] * de_factors
true_cluster_logfc <- log(mu_B / mu_A)

cat("Gene structure:\n")
cat("  non-DE genes:", n_nonde, "\n")
cat("  DE genes (5 each at 1.5x, 2x, 3x):", n_de, "\n")
cat("  mu_base summary:", round(summary(mu_base), 2), "\n\n")

sim_two_celltypes <- function(n_per, mu_A, mu_B, alpha_0, alpha_1,
                              a = 0.5, burst_factor = 5) {
    p_genes <- length(mu_A)
    Y <- matrix(0, 2 * n_per, p_genes)
    for (j in seq_len(p_genes)) {
        Y[1:n_per, j] <- sim_pb_gene(n_per, mu_A[j], a, burst_factor)
        Y[(n_per + 1):(2 * n_per), j] <- sim_pb_gene(n_per, mu_B[j], a, burst_factor)
    }
    X <- apply_dropout_fn(Y, pi_logistic,
                          list(alpha_0 = alpha_0, alpha_1 = alpha_1))
    list(Y = Y, X = X)
}

one_rep <- function(alpha_1) {
    sim <- sim_two_celltypes(n_cells_per_cluster, mu_A, mu_B, alpha_0, alpha_1)
    A_mask <- 1:n_cells_per_cluster
    B_mask <- (n_cells_per_cluster + 1):n_total
    X_A <- sim$X[A_mask, ]; X_B <- sim$X[B_mask, ]
    Y_A <- sim$Y[A_mask, ]; Y_B <- sim$Y[B_mask, ]

    naive_A <- fit_naive_zinb(X_A)
    naive_B <- fit_naive_zinb(X_B)
    naive_logfc <- log(naive_B[, "mu_hat"] / naive_A[, "mu_hat"])

    n_spike_per <- ceiling(spike_in_frac * n_cells_per_cluster)
    sa <- sample.int(n_cells_per_cluster, n_spike_per)
    sb <- sample.int(n_cells_per_cluster, n_spike_per)
    Y_pool <- rbind(Y_A[sa, ], Y_B[sb, ])
    X_pool <- rbind(X_A[sa, ], X_B[sb, ])
    drop_fit <- fit_dropout_logistic(Y_pool, X_pool)
    drop_reg <- fit_dropout_logistic_regularized(Y_pool, X_pool)

    safe_logfc <- function(fit) {
        if (is.null(fit)) return(rep(NA_real_, n_genes))
        oA <- fit_oracle_zinb_v3(X_A, fit)
        oB <- fit_oracle_zinb_v3(X_B, fit)
        log(oB[, "mu_hat"] / oA[, "mu_hat"])
    }
    oracle_logfc <- safe_logfc(drop_fit)
    reg_oracle_logfc <- safe_logfc(drop_reg)

    list(naive_logfc = naive_logfc, oracle_logfc = oracle_logfc,
         reg_oracle_logfc = reg_oracle_logfc)
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
    r_log <- t(sapply(reps, `[[`, "reg_oracle_logfc"))

    nm <- colMeans(n_log, na.rm = TRUE); ns <- apply(n_log, 2, sd, na.rm = TRUE)
    om <- colMeans(o_log, na.rm = TRUE); os <- apply(o_log, 2, sd, na.rm = TRUE)
    rm <- colMeans(r_log, na.rm = TRUE); rs <- apply(r_log, 2, sd, na.rm = TRUE)

    is_nonde <- seq_len(n_genes) <= n_nonde
    de_15 <- (n_nonde + 1):(n_nonde + 5)
    de_20 <- (n_nonde + 6):(n_nonde + 10)
    de_30 <- (n_nonde + 11):(n_nonde + 15)

    # logFC bias on non-DE (should be 0)
    cat(sprintf("  non-DE logFC mean abs: naive = %.3f, oracle = %.3f, REG = %.3f\n",
                mean(abs(nm[is_nonde])),
                mean(abs(om[is_nonde])),
                mean(abs(rm[is_nonde]))))
    # logFC bias on DE strata
    for (group in list(list(idx = de_15, name = "1.5x", true = log(1.5)),
                       list(idx = de_20, name = "2x",   true = log(2)),
                       list(idx = de_30, name = "3x",   true = log(3)))) {
        cat(sprintf("  DE %s : naive bias = %+.3f, oracle = %+.3f, REG = %+.3f\n",
                    group$name,
                    mean(nm[group$idx]) - group$true,
                    mean(om[group$idx]) - group$true,
                    mean(rm[group$idx]) - group$true))
    }

    # DE call rates at |z| > 2
    n_z <- nm / pmax(ns / sqrt(n_reps), 1e-8)
    o_z <- om / pmax(os / sqrt(n_reps), 1e-8)
    r_z <- rm / pmax(rs / sqrt(n_reps), 1e-8)
    cat(sprintf("  type-I |z|>2 (over %d non-DE): naive = %.3f, oracle = %.3f, REG = %.3f\n",
                n_nonde,
                mean(abs(n_z[is_nonde]) > 2),
                mean(abs(o_z[is_nonde]) > 2),
                mean(abs(r_z[is_nonde]) > 2)))
    for (group in list(list(idx = de_15, name = "1.5x"),
                       list(idx = de_20, name = "2x"),
                       list(idx = de_30, name = "3x"))) {
        cat(sprintf("  power |z|>2 at %s: naive = %.2f, oracle = %.2f, REG = %.2f\n",
                    group$name,
                    mean(abs(n_z[group$idx]) > 2),
                    mean(abs(o_z[group$idx]) > 2),
                    mean(abs(r_z[group$idx]) > 2)))
    }

    results[[sprintf("a%.2f", a1)]] <- list(
        alpha_1 = a1,
        naive_logfc = n_log, oracle_logfc = o_log, reg_logfc = r_log,
        n_mean = nm, o_mean = om, r_mean = rm,
        n_sd = ns, o_sd = os, r_sd = rs,
        is_nonde = is_nonde, de_15 = de_15, de_20 = de_20, de_30 = de_30,
        true_cluster_logfc = true_cluster_logfc
    )
}
saveRDS(results, "results_v6_cluster_de.rds")
cat("Done.\n")
