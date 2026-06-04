# V6-regularization: address the alpha_1 = 0 spike-fit attenuation by
# adding LRT-based selection between constant-dropout and Y-dependent
# (logistic) dropout. Compare four estimators:
#   1. naive ZINB
#   2. plain spike-in oracle (current)
#   3. regularized spike-in oracle (LRT-selected)
#   4. true-dropout oracle (no fit; reference)

source("sim.R")
set.seed(20260511)

n_cells_per_cond <- 1000
n_genes <- 50
n_de <- 10
de_factor <- 2; true_logfc_de <- log(de_factor)
spike_in_frac <- 0.10
alpha_0 <- 1.5
alpha_1_grid <- c(0.0, 0.3, 0.6, 1.0)
n_reps <- 60

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

    drop_plain <- fit_dropout_logistic(Y_pool, X_pool)
    drop_reg   <- fit_dropout_logistic_regularized(Y_pool, X_pool)
    drop_true  <- list(form = "logistic", pi_fn = pi_logistic,
                       params = list(alpha_0 = alpha_0, alpha_1 = alpha_1))

    safe_logfc <- function(fit) {
        if (is.null(fit)) return(rep(NA_real_, n_genes))
        oA <- fit_oracle_zinb_v3(sim_A$X, fit)
        oB <- fit_oracle_zinb_v3(sim_B$X, fit)
        log(oB[, "mu_hat"] / oA[, "mu_hat"])
    }

    list(alpha_1 = alpha_1,
         naive_logfc = naive_logfc,
         plain_oracle_logfc = safe_logfc(drop_plain),
         reg_oracle_logfc   = safe_logfc(drop_reg),
         true_oracle_logfc  = safe_logfc(drop_true),
         plain_pars = if (is.null(drop_plain)) c(NA, NA) else
             c(drop_plain$params$alpha_0, drop_plain$params$alpha_1),
         reg_pars = if (is.null(drop_reg)) c(NA, NA) else
             c(drop_reg$params$alpha_0, drop_reg$params$alpha_1),
         reg_selection = if (is.null(drop_reg)) NA else drop_reg$selection,
         reg_LRT_p = if (is.null(drop_reg)) NA else drop_reg$LRT_p)
}

results <- list()
for (a1 in alpha_1_grid) {
    cat(sprintf("---- alpha_1 = %g ----\n", a1))
    t0 <- Sys.time()
    reps <- replicate(n_reps, one_rep(a1), simplify = FALSE)
    cat("  ", n_reps, "reps in",
        round(as.numeric(Sys.time() - t0, units = "secs"), 1), "s\n")

    n_log <- t(sapply(reps, `[[`, "naive_logfc"))
    p_log <- t(sapply(reps, `[[`, "plain_oracle_logfc"))
    r_log <- t(sapply(reps, `[[`, "reg_oracle_logfc"))
    t_log <- t(sapply(reps, `[[`, "true_oracle_logfc"))
    plain_pars <- t(sapply(reps, `[[`, "plain_pars"))
    reg_pars   <- t(sapply(reps, `[[`, "reg_pars"))
    reg_sel    <- sapply(reps, `[[`, "reg_selection")

    is_de <- seq_len(n_genes) %in% de_idx
    n_mean <- colMeans(n_log, na.rm = TRUE)
    p_mean <- colMeans(p_log, na.rm = TRUE)
    r_mean <- colMeans(r_log, na.rm = TRUE)
    t_mean <- colMeans(t_log, na.rm = TRUE)

    de_naive_bias <- mean(n_mean[is_de]) - true_logfc_de
    de_plain_bias <- mean(p_mean[is_de]) - true_logfc_de
    de_reg_bias   <- mean(r_mean[is_de]) - true_logfc_de
    de_true_bias  <- mean(t_mean[is_de]) - true_logfc_de

    cat(sprintf("  DE bias: naive = %+.4f, PLAIN = %+.4f, REG = %+.4f, TRUE = %+.4f\n",
                de_naive_bias, de_plain_bias, de_reg_bias, de_true_bias))
    cat(sprintf("  Selection rates over %d reps:\n", n_reps))
    sel_table <- table(reg_sel)
    print(sel_table)
    cat(sprintf("  PLAIN spike-fit alpha_0 = %.3f, alpha_1 = %.3f (truth %.2f, %.2f)\n",
                mean(plain_pars[, 1], na.rm = TRUE),
                mean(plain_pars[, 2], na.rm = TRUE), alpha_0, a1))
    cat(sprintf("  REG   spike-fit alpha_0 = %.3f, alpha_1 = %.3f (truth %.2f, %.2f)\n",
                mean(reg_pars[, 1], na.rm = TRUE),
                mean(reg_pars[, 2], na.rm = TRUE), alpha_0, a1))

    results[[sprintf("a%.2f", a1)]] <- list(
        alpha_1 = a1,
        naive_logfc = n_log, plain_logfc = p_log,
        reg_logfc = r_log, true_logfc = t_log,
        is_de = is_de,
        de_naive_bias = de_naive_bias,
        de_plain_bias = de_plain_bias,
        de_reg_bias = de_reg_bias,
        de_true_bias = de_true_bias,
        plain_pars = plain_pars, reg_pars = reg_pars,
        reg_sel = reg_sel
    )
}
saveRDS(results, "results_v6_regularization.rds")
cat("Done.\n")
