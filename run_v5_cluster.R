# V5-cluster: per-cluster mu recovery.
# Same two-cell-type DGP as v4-celltypes, but fit ZINB per cluster using
# TRUE cell-type labels. Question: does oracle preserve cell-type-specific
# mu estimates and per-cluster DE detection?
#
# Within each cluster, the data is homogeneous, so we should expect oracle
# to do as well as in the IID case (~2% bias) and naive to do its usual ~50%.
# The novel test is per-cluster DE: are cluster-specific differences
# (mu_A_j vs mu_B_j) recovered correctly by oracle but distorted by naive?

source("sim.R")
set.seed(20260511)

n_cells <- 2000
n_genes <- 50
spike_in_frac <- 0.10
alpha_0 <- 1.5
alpha_1_grid <- c(0.0, 0.6, 1.0)
n_reps <- 30

mu_A <- realistic_mu(n_genes, log_mean = 0,   log_sd = 1.0, seed = 11)
mu_B <- realistic_mu(n_genes, log_mean = 1.5, log_sd = 1.0, seed = 22)
true_cluster_logfc <- log(mu_B / mu_A)

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
    A_mask <- sim$cell_type == "A"
    B_mask <- sim$cell_type == "B"
    X_A <- sim$X[A_mask, ]; X_B <- sim$X[B_mask, ]
    Y_A <- sim$Y[A_mask, ]; Y_B <- sim$Y[B_mask, ]

    naive_A <- fit_naive_zinb(X_A)
    naive_B <- fit_naive_zinb(X_B)
    naive_logfc <- log(naive_B[, "mu_hat"] / naive_A[, "mu_hat"])

    # Pool spike-ins from both clusters for dropout fit
    n_spike_per <- ceiling(spike_in_frac * sum(A_mask))
    sa <- sample.int(nrow(X_A), n_spike_per)
    sb <- sample.int(nrow(X_B), n_spike_per)
    Y_pool <- rbind(Y_A[sa, ], Y_B[sb, ])
    X_pool <- rbind(X_A[sa, ], X_B[sb, ])
    drop_fit <- fit_dropout_logistic(Y_pool, X_pool)

    if (is.null(drop_fit)) {
        oracle_A_mu <- rep(NA_real_, n_genes)
        oracle_B_mu <- rep(NA_real_, n_genes)
    } else {
        oracle_A <- fit_oracle_zinb_v3(X_A, drop_fit)
        oracle_B <- fit_oracle_zinb_v3(X_B, drop_fit)
        oracle_A_mu <- oracle_A[, "mu_hat"]
        oracle_B_mu <- oracle_B[, "mu_hat"]
    }
    oracle_logfc <- log(oracle_B_mu / oracle_A_mu)

    list(naive_A_mu = naive_A[, "mu_hat"], naive_B_mu = naive_B[, "mu_hat"],
         oracle_A_mu = oracle_A_mu, oracle_B_mu = oracle_B_mu,
         naive_logfc = naive_logfc, oracle_logfc = oracle_logfc)
}

results <- list()
for (a1 in alpha_1_grid) {
    cat(sprintf("---- alpha_1 = %g ----\n", a1))
    t0 <- Sys.time()
    reps <- replicate(n_reps, one_rep(a1), simplify = FALSE)
    cat("  ", n_reps, "reps in",
        round(as.numeric(Sys.time() - t0, units = "secs"), 1), "s\n")

    nA <- t(sapply(reps, `[[`, "naive_A_mu"))
    nB <- t(sapply(reps, `[[`, "naive_B_mu"))
    oA <- t(sapply(reps, `[[`, "oracle_A_mu"))
    oB <- t(sapply(reps, `[[`, "oracle_B_mu"))
    n_lfc <- t(sapply(reps, `[[`, "naive_logfc"))
    o_lfc <- t(sapply(reps, `[[`, "oracle_logfc"))

    nA_rb <- (colMeans(nA, na.rm = TRUE) - mu_A) / mu_A
    nB_rb <- (colMeans(nB, na.rm = TRUE) - mu_B) / mu_B
    oA_rb <- (colMeans(oA, na.rm = TRUE) - mu_A) / mu_A
    oB_rb <- (colMeans(oB, na.rm = TRUE) - mu_B) / mu_B

    n_lfc_bias <- colMeans(n_lfc, na.rm = TRUE) - true_cluster_logfc
    o_lfc_bias <- colMeans(o_lfc, na.rm = TRUE) - true_cluster_logfc

    qf <- function(v) round(quantile(v, c(0.05, 0.5, 0.95), na.rm = TRUE), 3)
    cat("  cluster A mu rel-bias (5/50/95):\n")
    cat("    naive :", qf(nA_rb), "\n")
    cat("    oracle:", qf(oA_rb), "\n")
    cat("  cluster B mu rel-bias (5/50/95):\n")
    cat("    naive :", qf(nB_rb), "\n")
    cat("    oracle:", qf(oB_rb), "\n")
    cat("  per-gene logFC bias (cluster-A vs cluster-B):\n")
    cat("    naive  abs-bias (5/50/95):", qf(abs(n_lfc_bias)), "\n")
    cat("    oracle abs-bias (5/50/95):", qf(abs(o_lfc_bias)), "\n")
    cat("    median |naive_bias| - |oracle_bias| =",
        round(median(abs(n_lfc_bias)) - median(abs(o_lfc_bias)), 4),
        "(positive = oracle helps)\n")

    # DE call rate at |z| > 2
    n_lfc_sd <- apply(n_lfc, 2, sd, na.rm = TRUE)
    o_lfc_sd <- apply(o_lfc, 2, sd, na.rm = TRUE)
    n_z <- colMeans(n_lfc, na.rm = TRUE) / pmax(n_lfc_sd / sqrt(n_reps), 1e-8)
    o_z <- colMeans(o_lfc, na.rm = TRUE) / pmax(o_lfc_sd / sqrt(n_reps), 1e-8)
    is_real_de <- abs(true_cluster_logfc) > log(1.5)
    cat(sprintf("  DE call rate at |z|>2 (true |logFC|>log(1.5) in %d genes):\n",
                sum(is_real_de)))
    cat(sprintf("    naive  power = %.2f, type-I = %.3f\n",
                mean(abs(n_z[is_real_de]) > 2),
                mean(abs(n_z[!is_real_de]) > 2)))
    cat(sprintf("    oracle power = %.2f, type-I = %.3f\n",
                mean(abs(o_z[is_real_de]) > 2),
                mean(abs(o_z[!is_real_de]) > 2)))

    results[[sprintf("a%.2f", a1)]] <- list(
        alpha_1 = a1,
        naive_A_mu = nA, naive_B_mu = nB,
        oracle_A_mu = oA, oracle_B_mu = oB,
        naive_A_relbias = nA_rb, naive_B_relbias = nB_rb,
        oracle_A_relbias = oA_rb, oracle_B_relbias = oB_rb,
        naive_logfc_bias = n_lfc_bias, oracle_logfc_bias = o_lfc_bias,
        true_cluster_logfc = true_cluster_logfc
    )
}
saveRDS(results, "results_v5_cluster.rds")
cat("Done.\n")
