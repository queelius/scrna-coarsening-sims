# V3d: realistic ERCC spike-in design.
# 92 ERCC transcripts at known concentrations spiked into all cells.
# Compare oracle quality vs the v3 "perfect spike-in" baseline.

source("sim.R")
set.seed(20260509)

n_cells <- 2000
n_genes <- 50
mu_vec  <- realistic_mu(n_genes, log_mean = 0.5, log_sd = 1.3, seed = 1)

# ERCC design: 92 transcripts, log-uniform concentrations from 0.5 to 50.
# Each ERCC has Poisson-distributed true count Y_ik ~ Poisson(c_k) per cell.
n_ercc <- 92
ercc_concentrations <- exp(seq(log(0.5), log(50), length.out = n_ercc))

alpha_0 <- 1.5
alpha_1_grid <- c(0.0, 0.3, 0.6, 1.0)
n_reps <- 25

# Generate gene data, generate ERCC data with same dropout function, fit
# dropout from ERCC, run oracle on genes.
one_rep_v3d <- function(alpha_1) {
    # Real genes: Poisson-Beta with logistic dropout
    sim_gene <- sim_pb_with_dropout(n_cells, mu_vec, pi_logistic,
                                    params = list(alpha_0 = alpha_0,
                                                  alpha_1 = alpha_1))
    # ERCC: Poisson true counts, same dropout
    Y_ercc <- matrix(0, n_cells, n_ercc)
    for (k in seq_len(n_ercc)) {
        Y_ercc[, k] <- rpois(n_cells, ercc_concentrations[k])
    }
    X_ercc <- apply_dropout_fn(Y_ercc, pi_logistic,
                               list(alpha_0 = alpha_0, alpha_1 = alpha_1))

    # ERCC-realistic dropout fit (uses only c_k, marginalizes over Y)
    drop_ercc <- fit_dropout_ercc_logistic(ercc_concentrations, X_ercc)

    # Compare against v3 perfect-spike-in baseline
    # Use 200-cell perfect spike-in (= 10% as before) on same gene data
    spike_idx <- sample.int(n_cells, 200)
    drop_perfect <- fit_dropout_logistic(sim_gene$Y[spike_idx, , drop = FALSE],
                                         sim_gene$X[spike_idx, , drop = FALSE])

    naive <- fit_naive_zinb(sim_gene$X)

    safe_oracle <- function(fit) {
        if (is.null(fit)) {
            return(rep(NA_real_, n_genes))
        }
        fit_oracle_zinb_v3(sim_gene$X, fit)[, "mu_hat"]
    }
    oracle_ercc <- safe_oracle(drop_ercc)
    oracle_perfect <- safe_oracle(drop_perfect)

    list(alpha_1 = alpha_1,
         naive_mu = naive[, "mu_hat"],
         oracle_ercc_mu = oracle_ercc,
         oracle_perfect_mu = oracle_perfect,
         observed_mean = colMeans(sim_gene$X),
         drop_ercc_pars = if (is.null(drop_ercc)) c(NA, NA) else
             c(drop_ercc$params$alpha_0, drop_ercc$params$alpha_1),
         drop_perfect_pars = if (is.null(drop_perfect)) c(NA, NA) else
             c(drop_perfect$params$alpha_0, drop_perfect$params$alpha_1))
}

results <- list()
for (a1 in alpha_1_grid) {
    cat(sprintf("---- alpha_1 = %g ----\n", a1))
    t0 <- Sys.time()
    reps <- replicate(n_reps, one_rep_v3d(a1), simplify = FALSE)
    cat("  ", n_reps, "reps in", round(as.numeric(Sys.time() - t0), 1), "s\n")

    naive_mu_mat   <- t(sapply(reps, `[[`, "naive_mu"))
    ercc_mu_mat    <- t(sapply(reps, `[[`, "oracle_ercc_mu"))
    perfect_mu_mat <- t(sapply(reps, `[[`, "oracle_perfect_mu"))
    drop_ercc_mat <- t(sapply(reps, `[[`, "drop_ercc_pars"))
    drop_perf_mat <- t(sapply(reps, `[[`, "drop_perfect_pars"))

    naive_rb   <- (colMeans(naive_mu_mat,   na.rm = TRUE) - mu_vec) / mu_vec
    ercc_rb    <- (colMeans(ercc_mu_mat,    na.rm = TRUE) - mu_vec) / mu_vec
    perfect_rb <- (colMeans(perfect_mu_mat, na.rm = TRUE) - mu_vec) / mu_vec

    qf <- function(v) round(quantile(v, c(0.05, 0.5, 0.95), na.rm = TRUE), 3)
    cat("  naive          relbias (5/50/95):", qf(naive_rb), "\n")
    cat("  oracle PERFECT relbias (5/50/95):", qf(perfect_rb), "\n")
    cat("  oracle ERCC    relbias (5/50/95):", qf(ercc_rb), "\n")
    cat("  ERCC-fit alpha_0 mean:", round(mean(drop_ercc_mat[, 1], na.rm = TRUE), 3),
        " alpha_1:", round(mean(drop_ercc_mat[, 2], na.rm = TRUE), 3),
        " (truth", alpha_0, ",", a1, ")\n")
    cat("  PERFECT-fit alpha_0 mean:", round(mean(drop_perf_mat[, 1], na.rm = TRUE), 3),
        " alpha_1:", round(mean(drop_perf_mat[, 2], na.rm = TRUE), 3), "\n")

    results[[as.character(a1)]] <- list(
        alpha_1 = a1,
        naive_relbias = naive_rb,
        oracle_perfect_relbias = perfect_rb,
        oracle_ercc_relbias = ercc_rb,
        drop_ercc_pars = drop_ercc_mat,
        drop_perfect_pars = drop_perf_mat
    )
}

saveRDS(results, "results_v3d.rds")
cat("\nDone. Saved results_v3d.rds\n")
