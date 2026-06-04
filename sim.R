# scRNA-seq vs masked-cause reliability: misspecification study
# Mirrors mdrelax simulation_study.R idiom.
#
# Validates the cross-referencer agent's prediction: a naive ZINB fit
# under C2 violation (dropout depends on latent count Y) recovers
# observed-mean (cell-total expression) consistently but produces
# biased per-gene mu_j estimates.

suppressPackageStartupMessages(library(MASS))

# -----------------------------------------------------------------
# Data generation: NB counts + Y-dependent dropout
# -----------------------------------------------------------------

#' Simulate scRNA-seq-like data with controllable C2 violation severity.
#'
#' True counts: Y_ij ~ NB(mu_j, phi).
#' Observed: X_ij = 0 with prob pi(Y_ij), else X_ij = Y_ij.
#' pi(y) = plogis(alpha_0 - alpha_1 * log(y + 1)).
#'
#' alpha_1 = 0 -> constant dropout (C2 holds, naive ZINB unbiased).
#' alpha_1 > 0 -> dropout depends on Y (C2 violated, the case of interest).
sim_scrna <- function(n_cells, mu_vec, phi = 1, alpha_0 = 2, alpha_1 = 0.6) {
    p_genes <- length(mu_vec)
    Y <- matrix(0, n_cells, p_genes)
    for (j in seq_len(p_genes)) {
        Y[, j] <- rnbinom(n_cells, mu = mu_vec[j], size = phi)
    }
    pi_drop <- plogis(alpha_0 - alpha_1 * log(Y + 1))
    keep <- matrix(rbinom(length(Y), 1, 1 - pi_drop), n_cells, p_genes)
    X <- Y * keep
    list(Y = Y, X = X,
         truth = list(mu = mu_vec, phi = phi,
                      alpha_0 = alpha_0, alpha_1 = alpha_1))
}

# -----------------------------------------------------------------
# Richer DGP: Poisson-Beta two-state expression model (bursty)
# -----------------------------------------------------------------

#' Sample from Poisson-Beta two-state expression model.
#'
#' Each cell has an active-promoter probability drawn from Beta(a_j, b_j),
#' then count drawn from Poisson(lambda_j * p_active). Mean of marginal
#' Y is mu_j = lambda_j * a_j / (a_j + b_j).
#'
#' Parametrization: pick mu_j and a fixed burst-shape ratio. Set
#' lambda_j = mu_j * burst_factor, derive b_j from the mean equation.
#'
#' This is the standard scRNA-seq generative model (Kim & Marioni 2013,
#' SymSim, Splat use variants).
sim_pb_gene <- function(n_cells, mu, a = 0.5, burst_factor = 5) {
    lambda <- mu * burst_factor
    b <- a * (lambda / mu - 1)  # = a * (burst_factor - 1)
    p_active <- rbeta(n_cells, a, b)
    rpois(n_cells, lambda * p_active)
}

#' Generate gene means from a realistic log-normal distribution.
#' Real scRNA-seq mu's typically span 0.01 to ~100 with heavy right tail.
realistic_mu <- function(n_genes, log_mean = 0, log_sd = 1.5, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    mu <- exp(rnorm(n_genes, mean = log_mean, sd = log_sd))
    pmax(mu, 0.05)  # floor to avoid degenerate genes
}

#' Apply a dropout model to a true count matrix Y.
#' dropout_type: "logistic" (pi = plogis(alpha_0 - alpha_1 * log(y+1)))
#'           or "exponential" (pi = (1 - exp(-baseline)) * exp(-alpha_1 * y))
apply_dropout <- function(Y, alpha_0, alpha_1, dropout_type = "logistic") {
    if (dropout_type == "logistic") {
        pi_drop <- plogis(alpha_0 - alpha_1 * log(Y + 1))
    } else if (dropout_type == "exponential") {
        # alpha_0 sets baseline-zero probability when y is large; alpha_1 is decay
        pi_floor <- 1 - exp(-exp(alpha_0) * 1)  # high baseline
        pi_drop <- pi_floor * exp(-alpha_1 * Y)
    } else {
        stop("unknown dropout_type")
    }
    keep <- matrix(rbinom(length(Y), 1, 1 - pi_drop),
                   nrow(Y), ncol(Y))
    Y * keep
}

#' Full simulation with Poisson-Beta DGP and chosen dropout type.
sim_scrna_pb <- function(n_cells, mu_vec, a = 0.5, burst_factor = 5,
                         alpha_0 = 2, alpha_1 = 0.6,
                         dropout_type = "logistic") {
    p_genes <- length(mu_vec)
    Y <- matrix(0, n_cells, p_genes)
    for (j in seq_len(p_genes)) {
        Y[, j] <- sim_pb_gene(n_cells, mu_vec[j], a = a, burst_factor = burst_factor)
    }
    X <- apply_dropout(Y, alpha_0, alpha_1, dropout_type)
    list(Y = Y, X = X,
         truth = list(mu = mu_vec, a = a, burst_factor = burst_factor,
                      alpha_0 = alpha_0, alpha_1 = alpha_1,
                      dropout_type = dropout_type))
}

# -----------------------------------------------------------------
# v3: explicit dropout-function abstraction for form-robustness tests
# -----------------------------------------------------------------

#' Logistic-in-log-y dropout: pi(y) = plogis(alpha_0 - alpha_1 * log(y+1))
pi_logistic <- function(y, alpha_0, alpha_1) {
    plogis(alpha_0 - alpha_1 * log(y + 1))
}

#' Exponential-in-y dropout: pi(y) = pi0 * exp(-alpha_1 * y)
#' pi0 is the dropout rate at y = 0; rate alpha_1 controls decay speed.
pi_exponential <- function(y, pi0, alpha_1) {
    pi0 * exp(-alpha_1 * y)
}

#' Apply a generic dropout function to a count matrix.
apply_dropout_fn <- function(Y, pi_fn, params) {
    pi_drop <- do.call(pi_fn, c(list(Y), params))
    pi_drop <- pmin(pmax(pi_drop, 0), 1)
    keep <- matrix(rbinom(length(Y), 1, 1 - pi_drop), nrow(Y), ncol(Y))
    Y * keep
}

#' Full PB simulation parametrized by dropout function and params.
sim_pb_with_dropout <- function(n_cells, mu_vec, pi_fn, params,
                                a = 0.5, burst_factor = 5) {
    p_genes <- length(mu_vec)
    Y <- matrix(0, n_cells, p_genes)
    for (j in seq_len(p_genes)) {
        Y[, j] <- sim_pb_gene(n_cells, mu_vec[j], a = a, burst_factor = burst_factor)
    }
    X <- apply_dropout_fn(Y, pi_fn, params)
    list(Y = Y, X = X, mu_vec = mu_vec, pi_fn = pi_fn, params = params)
}

#' Fit logistic dropout from spike-in (Y, X) pairs via direct MLE.
fit_dropout_logistic <- function(Y_spike, X_spike) {
    y <- as.vector(Y_spike); x <- as.vector(X_spike)
    eligible <- y > 0
    if (sum(eligible) < 20) return(NULL)
    dropped <- as.integer(x[eligible] == 0)
    y_e <- y[eligible]
    nll <- function(par) {
        pi_d <- plogis(par[1] - par[2] * log(y_e + 1))
        pi_d <- pmin(pmax(pi_d, 1e-10), 1 - 1e-10)
        -sum(dbinom(dropped, 1, pi_d, log = TRUE))
    }
    opt <- tryCatch(
        optim(c(0, 0.3), nll, method = "Nelder-Mead",
              control = list(maxit = 2000)),
        error = function(e) NULL
    )
    if (is.null(opt) || opt$convergence != 0) return(NULL)
    list(form = "logistic",
         pi_fn = pi_logistic,
         params = list(alpha_0 = opt$par[1], alpha_1 = opt$par[2]))
}

#' Fit exponential dropout from spike-in (Y, X) pairs via direct MLE.
fit_dropout_exponential <- function(Y_spike, X_spike) {
    y <- as.vector(Y_spike); x <- as.vector(X_spike)
    eligible <- y > 0
    if (sum(eligible) < 20) return(NULL)
    dropped <- as.integer(x[eligible] == 0)
    y_e <- y[eligible]
    nll <- function(par) {
        pi0 <- plogis(par[1]); a1 <- exp(par[2])
        pi_d <- pi0 * exp(-a1 * y_e)
        pi_d <- pmin(pmax(pi_d, 1e-10), 1 - 1e-10)
        -sum(dbinom(dropped, 1, pi_d, log = TRUE))
    }
    opt <- tryCatch(
        optim(c(qlogis(0.5), log(0.1)), nll, method = "Nelder-Mead",
              control = list(maxit = 2000)),
        error = function(e) NULL
    )
    if (is.null(opt) || opt$convergence != 0) return(NULL)
    list(form = "exponential",
         pi_fn = pi_exponential,
         params = list(pi0 = plogis(opt$par[1]), alpha_1 = exp(opt$par[2])))
}

#' Generic oracle ZINB NLL using a dropout fit object (form + params).
nll_oracle_gene_v3 <- function(par, x, dropout_fit, y_max = 200) {
    mu <- exp(par[1]); phi <- exp(par[2])
    y_grid <- 0:y_max
    p_y <- dnbinom(y_grid, mu = mu, size = phi)
    pi_drop <- do.call(dropout_fit$pi_fn, c(list(y_grid), dropout_fit$params))
    pi_drop <- pmin(pmax(pi_drop, 0), 1)
    p_x0 <- p_y[1] + sum(p_y[-1] * pi_drop[-1])
    is_zero <- x == 0
    ll_zero <- log(p_x0 + 1e-300) * sum(is_zero)
    pos_x <- x[!is_zero]
    pos_x_clip <- pmin(pos_x, y_max)
    ll_pos <- sum(log(p_y[pos_x_clip + 1] *
                          (1 - pi_drop[pos_x_clip + 1]) + 1e-300))
    -(ll_zero + ll_pos)
}

#' Regularized dropout fit: select between constant-dropout (alpha_1 = 0)
#' and Y-dependent (logistic) via likelihood-ratio test against chi-sq(1)
#' at p < 0.05. Addresses the alpha_1 = 0 spike-fit attenuation found in v5.
fit_dropout_logistic_regularized <- function(Y_spike, X_spike, alpha = 0.05) {
    y <- as.vector(Y_spike); x <- as.vector(X_spike)
    eligible <- y > 0
    if (sum(eligible) < 20) return(NULL)
    dropped <- as.integer(x[eligible] == 0)
    y_e <- y[eligible]

    # Constant-dropout fit
    nll_c <- function(par) {
        pi_d <- plogis(par[1])
        pi_d <- pmin(pmax(pi_d, 1e-10), 1 - 1e-10)
        -sum(dbinom(dropped, 1, pi_d, log = TRUE))
    }
    opt_c <- tryCatch(
        optim(0, nll_c, method = "Brent", lower = -10, upper = 10),
        error = function(e) NULL)

    # Logistic Y-dependent fit
    nll_l <- function(par) {
        pi_d <- plogis(par[1] - par[2] * log(y_e + 1))
        pi_d <- pmin(pmax(pi_d, 1e-10), 1 - 1e-10)
        -sum(dbinom(dropped, 1, pi_d, log = TRUE))
    }
    opt_l <- tryCatch(
        optim(c(0, 0.3), nll_l, method = "Nelder-Mead",
              control = list(maxit = 2000)),
        error = function(e) NULL)

    if (is.null(opt_l) || (!is.null(opt_l) && opt_l$convergence != 0)) {
        return(NULL)
    }
    if (is.null(opt_c)) {
        return(list(form = "logistic", pi_fn = pi_logistic,
                    params = list(alpha_0 = opt_l$par[1], alpha_1 = opt_l$par[2]),
                    selection = "logistic-only", LRT_p = NA))
    }

    LRT_stat <- 2 * (opt_c$value - opt_l$value)  # values are NLL, so positive when full fits better
    LRT_p <- pchisq(LRT_stat, df = 1, lower.tail = FALSE)

    if (LRT_p > alpha) {
        list(form = "logistic", pi_fn = pi_logistic,
             params = list(alpha_0 = opt_c$par, alpha_1 = 0),
             selection = "constant", LRT_p = LRT_p)
    } else {
        list(form = "logistic", pi_fn = pi_logistic,
             params = list(alpha_0 = opt_l$par[1], alpha_1 = opt_l$par[2]),
             selection = "logistic", LRT_p = LRT_p)
    }
}

#' Fit logistic dropout from ERCC-style spike-in data.
#' ERCCs have known concentration c_k but unobserved realized Y_ik.
#' Y_ik ~ Poisson(c_k). Marginal likelihood handles the latent Y.
#' This is the "realistic" version of fit_dropout_logistic.
fit_dropout_ercc_logistic <- function(C, X_ercc, y_max = 100) {
    n_cells <- nrow(X_ercc); K <- ncol(X_ercc)
    nll <- function(par) {
        a0 <- par[1]; a1 <- par[2]
        ll <- 0
        for (k in seq_len(K)) {
            c_k <- C[k]
            y_grid <- 0:y_max
            p_y <- dpois(y_grid, c_k)
            pi_drop <- plogis(a0 - a1 * log(y_grid + 1))
            pi_drop <- pmin(pmax(pi_drop, 0), 1)
            p_x0 <- p_y[1] + sum(p_y[-1] * pi_drop[-1])
            x_k <- X_ercc[, k]
            zero_count <- sum(x_k == 0)
            pos_x <- x_k[x_k > 0]
            pos_x_clip <- pmin(pos_x, y_max)
            ll <- ll + zero_count * log(p_x0 + 1e-300) +
                  sum(log(p_y[pos_x_clip + 1] *
                              (1 - pi_drop[pos_x_clip + 1]) + 1e-300))
        }
        -ll
    }
    opt <- tryCatch(
        optim(c(0, 0.3), nll, method = "Nelder-Mead",
              control = list(maxit = 2000)),
        error = function(e) NULL
    )
    if (is.null(opt) || opt$convergence != 0) return(NULL)
    list(form = "logistic",
         pi_fn = pi_logistic,
         params = list(alpha_0 = opt$par[1], alpha_1 = opt$par[2]))
}

#' Generic oracle ZINB fit using a dropout fit object.
#' y_max defaults to data-adaptive: max(y_max_floor, y_max_mult * max(x)),
#' capped at y_max_cap. Required for high-expression genes (e.g. real
#' scRNA-seq) where the latent Y grid must extend beyond the bulk of the
#' NB(mu, phi) mass to compute marginal P(X=0) correctly.
fit_oracle_zinb_v3 <- function(X, dropout_fit, y_max_floor = 500,
                               y_max_mult = 5, y_max_cap = 200000) {
    p_genes <- ncol(X)
    out <- matrix(NA_real_, p_genes, 2,
                  dimnames = list(NULL, c("mu_hat", "phi_hat")))
    for (j in seq_len(p_genes)) {
        x <- X[, j]
        x_max <- max(x)
        y_max <- min(max(y_max_floor, ceiling(y_max_mult * x_max)), y_max_cap)
        m_init <- max(mean(x[x > 0]), 0.5)
        start <- c(log(m_init), log(1))
        opt <- tryCatch(
            optim(start, nll_oracle_gene_v3, x = x,
                  dropout_fit = dropout_fit, y_max = y_max,
                  method = "Nelder-Mead",
                  control = list(maxit = 2000)),
            error = function(e) NULL
        )
        if (!is.null(opt) && opt$convergence == 0) {
            out[j, ] <- c(exp(opt$par[1]), exp(opt$par[2]))
        }
    }
    out
}

# -----------------------------------------------------------------
# Naive ZINB MLE (assumes per-gene constant dropout pi_j; misspecified)
# -----------------------------------------------------------------

#' Negative log-likelihood for naive ZINB on a single gene.
#' Parameters: log(mu), log(phi), logit(pi).
nll_zinb_gene <- function(par, x) {
    mu <- exp(par[1]); phi <- exp(par[2]); pi <- plogis(par[3])
    is_zero <- x == 0
    p_nb_zero <- dnbinom(0, mu = mu, size = phi)
    ll_zero_per <- log(pi + (1 - pi) * p_nb_zero + 1e-300)
    ll_pos  <- log(1 - pi + 1e-300) + dnbinom(x[!is_zero],
                                              mu = mu, size = phi, log = TRUE)
    -(ll_zero_per * sum(is_zero) + sum(ll_pos))
}

fit_naive_zinb <- function(X) {
    p_genes <- ncol(X)
    out <- matrix(NA_real_, p_genes, 3,
                  dimnames = list(NULL, c("mu_hat", "phi_hat", "pi_hat")))
    for (j in seq_len(p_genes)) {
        x <- X[, j]
        m_init <- max(mean(x[x > 0]), 0.5)
        start <- c(log(m_init), log(1), qlogis(min(max(mean(x == 0), 0.05), 0.95)))
        opt <- tryCatch(
            optim(start, nll_zinb_gene, x = x,
                  method = "Nelder-Mead",
                  control = list(maxit = 2000)),
            error = function(e) NULL
        )
        if (!is.null(opt) && opt$convergence == 0) {
            out[j, ] <- c(exp(opt$par[1]), exp(opt$par[2]), plogis(opt$par[3]))
        }
    }
    out
}

# -----------------------------------------------------------------
# Oracle: a subset of cells have spike-ins (true Y known).
# Use them to fit the dropout function pi(y) = plogis(a0 - a1 log(y+1)),
# then refit per-gene mu given the known dropout function.
# -----------------------------------------------------------------

fit_oracle_dropout <- function(Y_spike, X_spike) {
    y_flat <- as.vector(Y_spike); x_flat <- as.vector(X_spike)
    dropped <- as.integer(x_flat == 0 & y_flat > 0)
    eligible <- y_flat > 0
    if (sum(eligible) < 20 || sum(dropped[eligible]) < 5) return(c(NA, NA))
    fit <- tryCatch(
        glm(dropped[eligible] ~ log(y_flat[eligible] + 1),
            family = binomial()),
        error = function(e) NULL
    )
    if (is.null(fit)) c(NA, NA) else c(coef(fit)[1], -coef(fit)[2])
}

# Negative log-lik for ONE gene, with dropout function known (alpha_0, alpha_1).
# Marginalizes Y to get observed-X likelihood.
nll_oracle_gene <- function(par, x, alpha_0, alpha_1, y_max = 200) {
    mu <- exp(par[1]); phi <- exp(par[2])
    y_grid <- 0:y_max
    p_y <- dnbinom(y_grid, mu = mu, size = phi)
    pi_drop <- plogis(alpha_0 - alpha_1 * log(y_grid + 1))
    # P(X = 0) = sum_y p_y * pi_y + p_y(0) * (1 - pi_0) ... no wait.
    # X = 0 iff (Y = 0) OR (Y > 0 and dropped). All Y -> X=0 has prob pi_drop(y).
    # Wait: if Y = 0, X = 0 deterministically (Y * keep = 0). Dropout draw is irrelevant.
    # If Y > 0, X = 0 iff dropped, X = Y iff kept.
    # So P(X = 0) = p_y(0) + sum_{y>0} p_y * pi_drop(y).
    # P(X = x) for x > 0 = p_y(x) * (1 - pi_drop(x)).
    p_x0 <- p_y[1] + sum(p_y[-1] * pi_drop[-1])
    is_zero <- x == 0
    ll_zero <- log(p_x0 + 1e-300) * sum(is_zero)
    pos_x <- x[!is_zero]
    pos_x_clip <- pmin(pos_x, y_max)
    ll_pos <- sum(log(p_y[pos_x_clip + 1] * (1 - pi_drop[pos_x_clip + 1]) + 1e-300))
    -(ll_zero + ll_pos)
}

fit_oracle_zinb <- function(X, alpha_0, alpha_1) {
    p_genes <- ncol(X)
    out <- matrix(NA_real_, p_genes, 2, dimnames = list(NULL, c("mu_hat", "phi_hat")))
    for (j in seq_len(p_genes)) {
        x <- X[, j]
        m_init <- max(mean(x[x > 0]), 0.5)
        start <- c(log(m_init), log(1))
        opt <- tryCatch(
            optim(start, nll_oracle_gene, x = x,
                  alpha_0 = alpha_0, alpha_1 = alpha_1,
                  method = "Nelder-Mead",
                  control = list(maxit = 2000)),
            error = function(e) NULL
        )
        if (!is.null(opt) && opt$convergence == 0) {
            out[j, ] <- c(exp(opt$par[1]), exp(opt$par[2]))
        }
    }
    out
}

# -----------------------------------------------------------------
# One replicate of the study
# -----------------------------------------------------------------

one_replicate <- function(n_cells, mu_vec, phi = 1,
                          alpha_0 = 2, alpha_1 = 0.6,
                          spike_in_frac = 0.10) {
    sim <- sim_scrna(n_cells, mu_vec, phi, alpha_0, alpha_1)
    naive <- fit_naive_zinb(sim$X)

    n_spike <- ceiling(spike_in_frac * n_cells)
    spike_idx <- sample.int(n_cells, n_spike)
    Y_spike <- sim$Y[spike_idx, , drop = FALSE]
    X_spike <- sim$X[spike_idx, , drop = FALSE]
    drop_pars <- fit_oracle_dropout(Y_spike, X_spike)
    if (any(is.na(drop_pars))) {
        oracle <- matrix(NA_real_, length(mu_vec), 2,
                         dimnames = list(NULL, c("mu_hat", "phi_hat")))
    } else {
        oracle <- fit_oracle_zinb(sim$X, drop_pars[1], drop_pars[2])
    }

    list(
        truth_mu = mu_vec,
        naive_mu = naive[, "mu_hat"],
        naive_pi = naive[, "pi_hat"],
        oracle_mu = oracle[, "mu_hat"],
        observed_mean = colMeans(sim$X),
        truth_observed_mean = colMeans(sim$X),
        naive_predicted_mean = (1 - naive[, "pi_hat"]) * naive[, "mu_hat"],
        drop_pars = drop_pars
    )
}
