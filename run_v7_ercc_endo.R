# V7: ERCC-vs-endogenous capture-efficiency sensitivity.
# The DECENT paper notes that ERCC molecules and endogenous transcripts
# have different capture efficiencies. If we fit dropout from ERCC and
# apply to genes, we're using pi_ERCC instead of the true pi_endo.
# Question: how much bias does this gap induce?
#
# Parametrization: pi_endo(y) = plogis(alpha_0 - alpha_1*log(y+1))
#                  pi_ERCC(y) = plogis(alpha_0 + shift - alpha_1*log(y+1))
# Positive shift = ERCC has more dropout than endogenous (the realistic
# direction per DECENT). Sweep shift values.

source("sim.R")
set.seed(20260511)

n_cells <- 2000
n_genes <- 50
mu_vec  <- realistic_mu(n_genes, log_mean = 0.5, log_sd = 1.3, seed = 1)

alpha_0_endo <- 1.5
alpha_1      <- 0.6
shift_grid   <- c(-1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0)
n_reps       <- 25

# ERCC design (realistic)
n_ercc <- 92
ercc_concentrations <- exp(seq(log(0.5), log(50), length.out = n_ercc))

one_rep <- function(shift) {
    # Genes: Poisson-Beta with endogenous dropout
    sim_gene <- sim_pb_with_dropout(
        n_cells, mu_vec, pi_logistic,
        params = list(alpha_0 = alpha_0_endo, alpha_1 = alpha_1))

    # ERCC: Poisson true counts, ERCC-specific dropout (shifted alpha_0)
    Y_ercc <- matrix(0, n_cells, n_ercc)
    for (k in seq_len(n_ercc)) {
        Y_ercc[, k] <- rpois(n_cells, ercc_concentrations[k])
    }
    X_ercc <- apply_dropout_fn(
        Y_ercc, pi_logistic,
        list(alpha_0 = alpha_0_endo + shift, alpha_1 = alpha_1))

    # Oracle using ERCC-derived dropout
    drop_ercc <- fit_dropout_ercc_logistic(ercc_concentrations, X_ercc)
    # Oracle using TRUE endogenous dropout (reference)
    drop_true_endo <- list(form = "logistic", pi_fn = pi_logistic,
                           params = list(alpha_0 = alpha_0_endo,
                                         alpha_1 = alpha_1))

    naive <- fit_naive_zinb(sim_gene$X)

    safe_oracle <- function(fit) {
        if (is.null(fit)) return(rep(NA_real_, n_genes))
        fit_oracle_zinb_v3(sim_gene$X, fit)[, "mu_hat"]
    }
    oracle_ercc <- safe_oracle(drop_ercc)
    oracle_true <- safe_oracle(drop_true_endo)

    list(shift = shift,
         naive_mu = naive[, "mu_hat"],
         oracle_ercc_mu = oracle_ercc,
         oracle_true_mu = oracle_true,
         observed_mean = colMeans(sim_gene$X),
         drop_ercc_pars = if (is.null(drop_ercc)) c(NA, NA) else
             c(drop_ercc$params$alpha_0, drop_ercc$params$alpha_1))
}

results <- list()
for (sh in shift_grid) {
    cat(sprintf("---- shift = %+g (alpha_0_ERCC = %g) ----\n", sh, alpha_0_endo + sh))
    t0 <- Sys.time()
    reps <- replicate(n_reps, one_rep(sh), simplify = FALSE)
    cat("  ", n_reps, "reps in",
        round(as.numeric(Sys.time() - t0, units = "secs"), 1), "s\n")

    n_mu_mat   <- t(sapply(reps, `[[`, "naive_mu"))
    oe_mu_mat  <- t(sapply(reps, `[[`, "oracle_ercc_mu"))
    ot_mu_mat  <- t(sapply(reps, `[[`, "oracle_true_mu"))
    drop_pars  <- t(sapply(reps, `[[`, "drop_ercc_pars"))

    naive_rb     <- (colMeans(n_mu_mat,  na.rm = TRUE) - mu_vec) / mu_vec
    oracle_e_rb  <- (colMeans(oe_mu_mat, na.rm = TRUE) - mu_vec) / mu_vec
    oracle_t_rb  <- (colMeans(ot_mu_mat, na.rm = TRUE) - mu_vec) / mu_vec

    qf <- function(v) round(quantile(v, c(0.05, 0.5, 0.95), na.rm = TRUE), 3)
    cat("  naive          relbias (5/50/95):", qf(naive_rb), "\n")
    cat("  oracle ERCC    relbias (5/50/95):", qf(oracle_e_rb), "\n")
    cat("  oracle TRUE-endo (5/50/95):", qf(oracle_t_rb), "\n")
    cat(sprintf("  ERCC-fit alpha_0 mean = %.3f (truth at shift: %.2f),  alpha_1 mean = %.3f (truth %.2f)\n",
                mean(drop_pars[, 1], na.rm = TRUE), alpha_0_endo + sh,
                mean(drop_pars[, 2], na.rm = TRUE), alpha_1))

    # Per-mu-quintile naive vs oracle ERCC
    mu_q <- cut(mu_vec, breaks = quantile(mu_vec, seq(0, 1, 0.2)),
                include.lowest = TRUE, labels = paste0("Q", 1:5))
    cat("  oracle ERCC rel-bias by mu quintile:\n")
    by_quint <- tapply(oracle_e_rb, mu_q, function(v) round(mean(v, na.rm = TRUE), 3))
    print(by_quint)

    results[[sprintf("shift_%+g", sh)]] <- list(
        shift = sh,
        naive_relbias = naive_rb,
        oracle_ercc_relbias = oracle_e_rb,
        oracle_true_relbias = oracle_t_rb,
        drop_ercc_pars = drop_pars,
        by_quint = by_quint
    )
}
saveRDS(results, "results_v7_ercc_endo.rds")
cat("\nDone. Saved results_v7_ercc_endo.rds\n")
