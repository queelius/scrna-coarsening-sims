# V11c: Compare DECENT vs spike-fit oracle per-gene mu estimates on
# Klein 2015 K562 inDrops data. Combines the v11 (oracle, naive) and
# v11b (DECENT) results.

suppressPackageStartupMessages({
    library(data.table)
})
source("sim.R")

v11  <- readRDS("results_v11_klein_decent.rds")
v11b <- readRDS("results_v11b_klein_decent.rds")

decent <- v11b$decent_fit
oracle <- v11$oracle
naive  <- v11$naive
gene_names_oracle <- v11$gene_sub_names
obs_mean_oracle  <- v11$obs_mean

# DECENT was run on top 500 genes by mean; same ordering as our top
# 500 from v11 (which itself was top 2000 then top 500 for DECENT).
n_decent <- length(decent$est.mu)
cat("DECENT est.mu length:", n_decent, "\n")
cat("DECENT est.pi0 length:", length(decent$est.pi0), "\n")
cat("Oracle est length:", nrow(oracle), "\n")

# Match: DECENT used top 500 by mean. v11 oracle used top 2000 by mean.
# Identify the common 500 by ordering on obs_mean within oracle.
top500_oracle <- order(obs_mean_oracle, decreasing = TRUE)[1:n_decent]
oracle_500 <- oracle[top500_oracle, , drop = FALSE]
naive_500  <- naive[top500_oracle, , drop = FALSE]
obs_500    <- obs_mean_oracle[top500_oracle]
gene_500   <- gene_names_oracle[top500_oracle]

oracle_mu <- oracle_500[, "mu_hat"]
naive_mu  <- naive_500[, "mu_hat"]
decent_mu <- decent$est.mu

cat("\nDECENT est.mu summary:\n")
print(summary(decent_mu))
cat("\nOracle mu summary:\n")
print(summary(oracle_mu))
cat("\nNaive mu summary:\n")
print(summary(naive_mu))
cat("\nObs mean summary:\n")
print(summary(obs_500))

# Per-gene comparison: ratios
oracle_naive_ratio <- oracle_mu / naive_mu
decent_naive_ratio <- decent_mu / naive_mu
decent_oracle_ratio <- decent_mu / oracle_mu

cat("\n=== oracle / naive ratio quantiles (top 500) ===\n")
print(round(quantile(oracle_naive_ratio,
                     c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE), 3))
cat("\n=== DECENT / naive ratio quantiles ===\n")
print(round(quantile(decent_naive_ratio,
                     c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE), 3))
cat("\n=== DECENT / oracle ratio quantiles ===\n")
print(round(quantile(decent_oracle_ratio,
                     c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE), 3))

# Correlation of estimates
cat("\n=== Pearson correlation of mu estimates (top 500) ===\n")
cat("  oracle vs naive  :", round(cor(oracle_mu, naive_mu, use = "complete.obs"), 4), "\n")
cat("  DECENT vs naive  :", round(cor(decent_mu, naive_mu, use = "complete.obs"), 4), "\n")
cat("  DECENT vs oracle :", round(cor(decent_mu, oracle_mu, use = "complete.obs"), 4), "\n")
cat("  DECENT vs obs    :", round(cor(decent_mu, obs_500, use = "complete.obs"), 4), "\n")
cat("  oracle vs obs    :", round(cor(oracle_mu, obs_500, use = "complete.obs"), 4), "\n")

cat("\n=== log-scale Pearson correlation ===\n")
lo <- function(x) log10(pmax(x, 0.01))
cat("  log(oracle) vs log(naive) :", round(cor(lo(oracle_mu), lo(naive_mu)), 4), "\n")
cat("  log(DECENT) vs log(naive) :", round(cor(lo(decent_mu), lo(naive_mu)), 4), "\n")
cat("  log(DECENT) vs log(oracle):", round(cor(lo(decent_mu), lo(oracle_mu)), 4), "\n")
cat("  log(DECENT) vs log(obs)   :", round(cor(lo(decent_mu), lo(obs_500)), 4), "\n")
cat("  log(oracle) vs log(obs)   :", round(cor(lo(oracle_mu), lo(obs_500)), 4), "\n")

# Save merged for inclusion in paper
saveRDS(list(
    obs_mean = obs_500,
    naive_mu = naive_mu,
    oracle_mu = oracle_mu,
    decent_mu = decent_mu,
    gene_names = gene_500,
    decent_CE = decent$CE,
    decent_tau1 = decent$tau1,
    decent_tau0 = decent$tau0,
    oracle_naive_quantiles = round(quantile(oracle_naive_ratio,
                                            c(0.05, 0.25, 0.5, 0.75, 0.95),
                                            na.rm = TRUE), 3),
    decent_naive_quantiles = round(quantile(decent_naive_ratio,
                                            c(0.05, 0.25, 0.5, 0.75, 0.95),
                                            na.rm = TRUE), 3),
    decent_oracle_quantiles = round(quantile(decent_oracle_ratio,
                                             c(0.05, 0.25, 0.5, 0.75, 0.95),
                                             na.rm = TRUE), 3),
    decent_class = class(decent),
    decent_CE_value = decent$CE
), "results_v11c_klein_compare.rds")
cat("\nSaved results_v11c_klein_compare.rds\n")
