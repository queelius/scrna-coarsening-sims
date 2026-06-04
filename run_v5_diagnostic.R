# V5-diagnostic: decompose the oracle bias at alpha_1 = 0.
# Three estimators run on the same data:
#   1. naive ZINB (constant per-gene pi)
#   2. spike-in oracle (fits dropout from spike-ins, uses it)
#   3. true-dropout oracle (uses the true (alpha_0, alpha_1), no fit noise)
# Compare DE bias trajectories. If true-dropout oracle has the same
# alpha_1 = 0 bias as spike-in oracle, the artifact is ZINB-vs-PB
# distributional misspec. If true-dropout oracle is unbiased, it is
# spike-in fit noise.

source("sim.R")
set.seed(20260511)

n_cells_per_cond <- 1000
n_genes <- 50
n_de <- 10
de_factor <- 2; true_logfc_de <- log(de_factor)
spike_in_frac <- 0.10
alpha_0 <- 1.5
alpha_1_grid <- c(0.0, 0.3, 0.6, 1.0)
n_reps <- 80

mu_A <- realistic_mu(n_genes, log_mean = 0.5, log_sd = 1.3, seed = 1)
mu_B <- mu_A
de_idx <- sample.int(n_genes, n_de)
mu_B[de_idx] <- mu_B[de_idx] * de_factor

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
    drop_fit_spike <- fit_dropout_logistic(Y_pool, X_pool)

    drop_fit_true <- list(form = "logistic", pi_fn = pi_logistic,
                          params = list(alpha_0 = alpha_0, alpha_1 = alpha_1))

    safe_logfc <- function(fit) {
        if (is.null(fit)) return(rep(NA_real_, n_genes))
        oA <- fit_oracle_zinb_v3(sim_A$X, fit)
        oB <- fit_oracle_zinb_v3(sim_B$X, fit)
        log(oB[, "mu_hat"] / oA[, "mu_hat"])
    }
    spike_logfc <- safe_logfc(drop_fit_spike)
    true_logfc  <- safe_logfc(drop_fit_true)

    list(alpha_1 = alpha_1,
         naive_logfc = naive_logfc,
         spike_oracle_logfc = spike_logfc,
         true_oracle_logfc = true_logfc,
         spike_drop_pars = if (is.null(drop_fit_spike)) c(NA, NA) else
             c(drop_fit_spike$params$alpha_0, drop_fit_spike$params$alpha_1))
}

results <- list()
for (a1 in alpha_1_grid) {
    cat(sprintf("---- alpha_1 = %g ----\n", a1))
    t0 <- Sys.time()
    reps <- replicate(n_reps, one_rep(a1), simplify = FALSE)
    cat("  ", n_reps, "reps in",
        round(as.numeric(Sys.time() - t0, units = "secs"), 1), "s\n")

    n_log <- t(sapply(reps, `[[`, "naive_logfc"))
    s_log <- t(sapply(reps, `[[`, "spike_oracle_logfc"))
    t_log <- t(sapply(reps, `[[`, "true_oracle_logfc"))
    drop_pars <- t(sapply(reps, `[[`, "spike_drop_pars"))

    is_de <- seq_len(n_genes) %in% de_idx
    n_mean <- colMeans(n_log, na.rm = TRUE)
    s_mean <- colMeans(s_log, na.rm = TRUE)
    t_mean <- colMeans(t_log, na.rm = TRUE)

    de_naive_bias  <- mean(n_mean[is_de]) - true_logfc_de
    de_spike_bias  <- mean(s_mean[is_de]) - true_logfc_de
    de_true_bias   <- mean(t_mean[is_de]) - true_logfc_de
    nd_naive_max  <- max(abs(n_mean[!is_de]))
    nd_spike_max  <- max(abs(s_mean[!is_de]))
    nd_true_max   <- max(abs(t_mean[!is_de]))

    cat(sprintf("  DE bias: naive = %+.4f, spike-oracle = %+.4f, TRUE-oracle = %+.4f\n",
                de_naive_bias, de_spike_bias, de_true_bias))
    cat(sprintf("  Non-DE max abs: naive = %.3f, spike-oracle = %.3f, TRUE-oracle = %.3f\n",
                nd_naive_max, nd_spike_max, nd_true_max))
    cat(sprintf("  spike-fit alpha_0 mean = %.3f (truth %.2f); alpha_1 mean = %.3f (truth %.2f)\n",
                mean(drop_pars[, 1], na.rm = TRUE), alpha_0,
                mean(drop_pars[, 2], na.rm = TRUE), a1))
    # Per-gene-type breakdown for DE: how does the bias depend on mu?
    cat("  per-DE-gene bias (true-oracle) by true mu:\n")
    de_mus <- mu_A[de_idx]; de_bias_t <- t_mean[is_de] - true_logfc_de
    o_idx <- order(de_mus)
    for (k in o_idx) {
        cat(sprintf("    mu_A = %5.2f, bias = %+.4f\n", de_mus[k], de_bias_t[k]))
    }

    results[[sprintf("a%.2f", a1)]] <- list(
        alpha_1 = a1,
        naive_logfc = n_log, spike_oracle_logfc = s_log,
        true_oracle_logfc = t_log,
        is_de = is_de,
        de_naive_bias = de_naive_bias,
        de_spike_bias = de_spike_bias,
        de_true_bias  = de_true_bias,
        spike_drop_pars = drop_pars
    )
}
saveRDS(results, "results_v5_diagnostic.rds")
cat("Done.\n")
