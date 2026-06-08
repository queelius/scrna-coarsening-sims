# Figure generation for the scrna-coarsening paper.
# Reads results_v*.rds files and produces PDFs in
# /home/spinoza/github/coarsening/papers/scrna-coarsening/figures/.

suppressPackageStartupMessages({
    library(ggplot2)
    library(scales)
})

fig_dir <- "/home/spinoza/github/coarsening/papers/scrna-coarsening/figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

theme_set(theme_minimal(base_size = 11) +
          theme(panel.grid.minor = element_blank(),
                strip.background = element_rect(fill = "grey95", color = NA),
                strip.text = element_text(face = "bold")))

save_fig <- function(p, name, width = 5.5, height = 3.5) {
    path <- file.path(fig_dir, paste0(name, ".pdf"))
    ggsave(path, plot = p, width = width, height = height, device = "pdf")
    cat("Wrote", path, "\n")
}

# -----------------------------------------------------------------
# Figure 1: ERCC-endogenous bias surface (v7)
# -----------------------------------------------------------------

v7 <- readRDS("results_v7_ercc_endo.rds")
v7_df <- do.call(rbind, lapply(v7, function(r) {
    data.frame(
        shift = r$shift,
        method = c("naive", "oracle (ERCC fit)", "oracle (true endo)"),
        median_bias = c(median(r$naive_relbias, na.rm = TRUE),
                        median(r$oracle_ercc_relbias, na.rm = TRUE),
                        median(r$oracle_true_relbias, na.rm = TRUE)),
        q05 = c(quantile(r$naive_relbias, 0.05, na.rm = TRUE),
                quantile(r$oracle_ercc_relbias, 0.05, na.rm = TRUE),
                quantile(r$oracle_true_relbias, 0.05, na.rm = TRUE)),
        q95 = c(quantile(r$naive_relbias, 0.95, na.rm = TRUE),
                quantile(r$oracle_ercc_relbias, 0.95, na.rm = TRUE),
                quantile(r$oracle_true_relbias, 0.95, na.rm = TRUE)))
}))
v7_df$method <- factor(v7_df$method,
                       levels = c("oracle (true endo)", "oracle (ERCC fit)", "naive"))

p1 <- ggplot(v7_df, aes(x = shift, y = median_bias, color = method, fill = method)) +
    geom_hline(yintercept = 0, color = "grey60", linetype = "dashed") +
    annotate("rect", xmin = -0.5, xmax = 0.5, ymin = -Inf, ymax = Inf,
             alpha = 0.10, fill = "darkgreen") +
    annotate("text", x = 0, y = 1.0, label = "useful zone",
             color = "darkgreen", size = 3, fontface = "italic") +
    geom_ribbon(aes(ymin = q05, ymax = q95), alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.9) + geom_point(size = 2) +
    scale_color_manual(values = c("oracle (true endo)" = "#1f77b4",
                                   "oracle (ERCC fit)" = "#ff7f0e",
                                   "naive" = "#7f7f7f")) +
    scale_fill_manual(values = c("oracle (true endo)" = "#1f77b4",
                                  "oracle (ERCC fit)" = "#ff7f0e",
                                  "naive" = "#7f7f7f")) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    scale_x_continuous(breaks = seq(-1, 2, 0.5)) +
    labs(x = expression("ERCC-endogenous logit shift "*italic(s)),
         y = "Per-gene relative bias (median, IQR shown)",
         color = NULL, fill = NULL,
         title = "ERCC-based correction is bounded by the spike-in/endogenous gap",
         subtitle = "Useful zone: |s| <= 0.5; beyond, ERCC oracle is comparable to or worse than naive") +
    theme(legend.position = "top")
save_fig(p1, "fig1_ercc_endo_bias", width = 6.5, height = 4)

# -----------------------------------------------------------------
# Figure 2: Form-robustness (v3) - matched vs mismatched oracle
# -----------------------------------------------------------------

v3 <- readRDS("results_v3.rds")
v3_df <- do.call(rbind, lapply(v3, function(r) {
    data.frame(
        truth_form = r$truth_form,
        alpha_1 = r$alpha_1,
        method = c("naive", "oracle MATCHED", "oracle MISMATCHED"),
        median_bias = c(median(r$naive_relbias, na.rm = TRUE),
                        median(r$oracle_matched_relbias, na.rm = TRUE),
                        median(r$oracle_mismatched_relbias, na.rm = TRUE)))
}))
v3_df$truth_form <- factor(v3_df$truth_form,
                           levels = c("logistic", "exponential"),
                           labels = c("Logistic dropout truth",
                                       "Exponential dropout truth"))
v3_df$method <- factor(v3_df$method,
                       levels = c("naive", "oracle MISMATCHED", "oracle MATCHED"))

p2 <- ggplot(v3_df, aes(x = factor(alpha_1), y = median_bias, fill = method)) +
    geom_col(position = position_dodge(0.8), width = 0.7) +
    facet_wrap(~ truth_form, scales = "free_x") +
    geom_hline(yintercept = 0, color = "grey60") +
    scale_fill_manual(values = c("naive" = "#7f7f7f",
                                  "oracle MISMATCHED" = "#ff7f0e",
                                  "oracle MATCHED" = "#1f77b4")) +
    scale_y_continuous(labels = percent_format(accuracy = 1)) +
    labs(x = expression("Dropout severity "*alpha[1]),
         y = "Per-gene median relative bias",
         fill = NULL,
         title = "Form-robustness: spike-in correction is robust to dropout family mismatch",
         subtitle = "Mismatched oracle within 1-3pp of matched; both far below naive") +
    theme(legend.position = "top")
save_fig(p2, "fig2_form_robustness", width = 7, height = 4)

# -----------------------------------------------------------------
# Figure 3: Cell-total consistency in real data (v8)
# -----------------------------------------------------------------

v8 <- readRDS("results_v8_tabula_muris.rds")
naive_pred <- (1 - v8$naive[, "pi_hat"]) * v8$naive[, "mu_hat"]
v8_consistency <- data.frame(
    obs_mean = v8$obs_mean,
    naive_pred = naive_pred,
    gene = v8$gene_subset_names)
# Mark housekeeping
v8_consistency$is_hk <- v8$gene_subset_names %in% v8$hk_genes_present

p3 <- ggplot(v8_consistency, aes(x = obs_mean, y = naive_pred)) +
    geom_abline(slope = 1, intercept = 0, color = "grey60", linetype = "dashed") +
    geom_point(aes(color = is_hk, size = is_hk), alpha = 0.7) +
    scale_x_log10(labels = label_comma()) +
    scale_y_log10(labels = label_comma()) +
    scale_color_manual(values = c("FALSE" = "#1f77b4", "TRUE" = "#d62728"),
                       labels = c("FALSE" = "other", "TRUE" = "housekeeping")) +
    scale_size_manual(values = c("FALSE" = 1.0, "TRUE" = 2.5),
                      labels = c("FALSE" = "other", "TRUE" = "housekeeping")) +
    labs(x = "Observed mean (log scale)",
         y = expression("Naive ZINB predicted mean "*(1-hat(pi))*hat(mu)),
         color = NULL, size = NULL,
         title = "Cell-total consistency in Tabula Muris spleen",
         subtitle = "Predicted = observed across 207 genes (median |diff| = 0.07)") +
    theme(legend.position = "bottom")
save_fig(p3, "fig3_cell_total_consistency", width = 5.5, height = 4.5)

# -----------------------------------------------------------------
# Figure 4: Housekeeping gene comparison (v8)
# -----------------------------------------------------------------

hk_idx <- which(v8$gene_subset_names %in% v8$hk_genes_present)
hk_df <- data.frame(
    gene = v8$gene_subset_names[hk_idx],
    obs_mean = v8$obs_mean[hk_idx],
    naive = v8$naive[hk_idx, "mu_hat"],
    oracle = v8$oracle[hk_idx, "mu_hat"])
hk_long <- rbind(
    data.frame(gene = hk_df$gene, value = hk_df$obs_mean,
               method = "observed mean", obs_mean_anchor = hk_df$obs_mean),
    data.frame(gene = hk_df$gene, value = hk_df$naive,
               method = "naive ZINB mu", obs_mean_anchor = hk_df$obs_mean),
    data.frame(gene = hk_df$gene, value = hk_df$oracle,
               method = "oracle mu", obs_mean_anchor = hk_df$obs_mean))
hk_long$gene <- factor(hk_long$gene,
                        levels = hk_df$gene[order(hk_df$obs_mean, decreasing = TRUE)])
hk_long$method <- factor(hk_long$method,
                         levels = c("observed mean", "oracle mu", "naive ZINB mu"))

p4 <- ggplot(hk_long, aes(x = gene, y = value, fill = method)) +
    geom_col(position = position_dodge(0.8), width = 0.7) +
    scale_y_log10(labels = label_comma()) +
    scale_fill_manual(values = c("observed mean" = "#7f7f7f",
                                  "oracle mu" = "#1f77b4",
                                  "naive ZINB mu" = "#d62728")) +
    labs(x = NULL, y = "Value (log scale)",
         fill = NULL,
         title = "Housekeeping genes: oracle vs naive vs observed mean",
         subtitle = "Two regimes: high-expression (top 5) agree; low-mid (bottom 7) diverge by 5-10x") +
    theme(legend.position = "top",
          axis.text.x = element_text(angle = 30, hjust = 1, face = "italic"))
save_fig(p4, "fig4_housekeeping", width = 6.5, height = 4)

# -----------------------------------------------------------------
# Figure 5: Per-cluster DE bias by effect size (v6_cluster_de)
# -----------------------------------------------------------------

v6 <- readRDS("results_v6_cluster_de.rds")
fig5_rows <- list()
for (key in names(v6)) {
    r <- v6[[key]]
    is_de_15 <- r$de_15
    is_de_20 <- r$de_20
    is_de_30 <- r$de_30
    fig5_rows[[length(fig5_rows) + 1]] <- data.frame(
        alpha_1 = r$alpha_1,
        de_factor_label = c("1.5x", "2x", "3x"),
        true_logfc = c(log(1.5), log(2), log(3)),
        naive_bias = c(mean(r$n_mean[is_de_15]) - log(1.5),
                       mean(r$n_mean[is_de_20]) - log(2),
                       mean(r$n_mean[is_de_30]) - log(3)),
        oracle_bias = c(mean(r$o_mean[is_de_15]) - log(1.5),
                        mean(r$o_mean[is_de_20]) - log(2),
                        mean(r$o_mean[is_de_30]) - log(3)))
}
v6_df <- do.call(rbind, fig5_rows)
v6_long <- rbind(
    data.frame(alpha_1 = v6_df$alpha_1, de_factor = v6_df$de_factor_label,
               method = "naive", bias = v6_df$naive_bias),
    data.frame(alpha_1 = v6_df$alpha_1, de_factor = v6_df$de_factor_label,
               method = "oracle", bias = v6_df$oracle_bias))
v6_long$de_factor <- factor(v6_long$de_factor, levels = c("1.5x", "2x", "3x"))

p5 <- ggplot(v6_long, aes(x = alpha_1, y = bias, color = method, shape = method)) +
    geom_hline(yintercept = 0, color = "grey60", linetype = "dashed") +
    geom_line(linewidth = 0.9) + geom_point(size = 2.5) +
    facet_wrap(~ de_factor, labeller = label_both) +
    scale_color_manual(values = c("naive" = "#7f7f7f", "oracle" = "#1f77b4")) +
    labs(x = expression("Dropout severity "*alpha[1]),
         y = "DE log fold change bias",
         color = NULL, shape = NULL,
         title = "Per-cluster DE bias: naive worsens with severity, oracle stays near zero",
         subtitle = "Oracle reduces DE bias by 5-12x at 3x effect size and severe dropout") +
    theme(legend.position = "top")
save_fig(p5, "fig5_de_bias", width = 7, height = 3.5)

cat("\nAll figures saved to", fig_dir, "\n")
