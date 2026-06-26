###############################################################################
# 03_ptau_protein_age_trend_AD.R
# Age-trend tracking analysis: AD-spectrum (EMCI + LMCI + AD) only
# - Same method as 02, but restricted to cognitively impaired subjects
# - Which proteins' age-trajectory best matches PTAU in the AD continuum?
###############################################################################

# ===========================================================================
# 0. Setup
# ===========================================================================
library(readxl)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)

set.seed(42)

# ---- Configurable parameters ----
LOESS_SPAN      <- 0.75
MIN_VALID_PAIRS <- 30
TOP_N_PLOT      <- 10
TOP_N_EXPORT    <- 20

# Output directory
out_dir <- "output/ptau_protein_age_trend_AD"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("LOESS span: %.2f | Min valid pairs: %d | Top N plot: %d\n",
            LOESS_SPAN, MIN_VALID_PAIRS, TOP_N_PLOT))
cat("Population: AD spectrum only (EMCI + LMCI + AD)\n")

# ===========================================================================
# 1. Load data
# ===========================================================================
cat("\n========== Loading data ==========\n")
df_raw <- read_excel("master_data.xlsx", sheet = "Sheet1")
cat(sprintf("Raw data: %d rows x %d columns\n", nrow(df_raw), ncol(df_raw)))

dict <- read.csv("protein_raw_data/protein dict.csv", stringsAsFactors = FALSE)
cat(sprintf("Protein dictionary: %d entries\n", nrow(dict)))

# ===========================================================================
# 2. Filter: baseline → non-missing PTAU & AGE → AD spectrum DX
# ===========================================================================
cat("\n========== Filtering samples ==========\n")

cat("VISCODE2 distribution:\n")
print(table(df_raw$VISCODE2, useNA = "always"))

# Baseline only
df_bl <- df_raw[df_raw$VISCODE2 == "bl", ]
cat(sprintf("Baseline samples: %d rows\n", nrow(df_bl)))

# Non-missing PTAU and AGE
df_analysis <- df_bl[!is.na(df_bl$PTAU) & !is.na(df_bl$AGE), ]
cat(sprintf("Baseline with non-missing PTAU & AGE: %d samples\n", nrow(df_analysis)))

# Report full DX distribution before filtering
cat("\nDX distribution before AD-spectrum filter:\n")
print(table(df_analysis$DX, useNA = "always"))

# Filter to EMCI + LMCI + AD only
ad_spectrum <- c("EMCI", "LMCI", "AD")
df_analysis <- df_analysis[df_analysis$DX %in% ad_spectrum, ]
cat(sprintf("\nAfter filtering to EMCI + LMCI + AD: %d samples\n", nrow(df_analysis)))

cat("\nDX distribution (analysis set):\n")
print(table(df_analysis$DX, useNA = "always"))

# Report AGE distribution
cat(sprintf("\nAGE range: %.1f - %.1f | Median: %.1f\n",
            min(df_analysis$AGE), max(df_analysis$AGE), median(df_analysis$AGE)))

if (nrow(df_analysis) < MIN_VALID_PAIRS) {
  stop(sprintf("Insufficient samples: %d < %d (MIN_VALID_PAIRS). Cannot proceed.",
               nrow(df_analysis), MIN_VALID_PAIRS))
}

# ===========================================================================
# 3. Identify protein columns
# ===========================================================================
cat("\n========== Identifying protein columns ==========\n")
protein_cols <- grep("^X[0-9]+\\.[0-9]+$", colnames(df_analysis), value = TRUE)
cat(sprintf("SomaScan protein columns identified: %d\n", length(protein_cols)))

# ===========================================================================
# 4. Convert protein columns to numeric
# ===========================================================================
cat("\n========== Converting protein columns to numeric ==========\n")

convert_protein <- function(x) {
  x[x == "NA" | x == ""] <- NA
  as.numeric(x)
}

na_string_counts <- sapply(df_analysis[protein_cols[1:min(10, length(protein_cols))]],
                           function(x) sum(x == "NA", na.rm = TRUE))
cat("Example 'NA' string counts in first 10 protein columns:\n")
print(na_string_counts)

df_analysis[protein_cols] <- lapply(df_analysis[protein_cols], convert_protein)

cat("\nAfter conversion - sample protein column types:\n")
print(sapply(df_analysis[protein_cols[1:min(5, length(protein_cols))]], class))

# ===========================================================================
# 5. Filter proteins by missingness
# ===========================================================================
cat("\n========== Filtering proteins by missingness (>20%) ==========\n")

missing_rate <- sapply(df_analysis[protein_cols], function(x) mean(is.na(x)) * 100)

proteins_high_missing <- names(missing_rate[missing_rate > 20])
proteins_kept         <- names(missing_rate[missing_rate <= 20])

cat(sprintf("Proteins with >20%% missing: %d (removed)\n", length(proteins_high_missing)))
cat(sprintf("Proteins with <=20%% missing: %d (kept)\n", length(proteins_kept)))

cat("\nTop 10 proteins with highest missing rate:\n")
print(head(sort(missing_rate, decreasing = TRUE), 10))

if (length(proteins_high_missing) > 0) {
  removed_df <- data.frame(
    protein_id      = proteins_high_missing,
    missing_rate_pct = round(missing_rate[proteins_high_missing], 2),
    stringsAsFactors = FALSE
  )
  write.csv(removed_df,
            file.path(out_dir, "removed_proteins_high_missing.csv"),
            row.names = FALSE)
}

# ===========================================================================
# 6. Log2 transformation
# ===========================================================================
cat("\n========== Log2 transformation ==========\n")

cat(sprintf("PTAU range before log2: %.2f - %.2f\n",
            min(df_analysis$PTAU), max(df_analysis$PTAU)))
df_analysis$PTAU_log2 <- log2(df_analysis$PTAU)
cat(sprintf("PTAU_log2 range: %.2f - %.2f\n",
            min(df_analysis$PTAU_log2), max(df_analysis$PTAU_log2)))

cat("Log2 transforming protein values...\n")
protein_mat   <- as.matrix(df_analysis[proteins_kept])
protein_log2  <- log2(protein_mat)

neg_count <- sum(protein_mat <= 0, na.rm = TRUE)
cat(sprintf("Protein values <= 0: %d\n", neg_count))
if (neg_count > 0) {
  cat("  Replacing -Inf/NaN with NA\n")
  protein_log2[is.infinite(protein_log2) | is.nan(protein_log2)] <- NA
}

cat(sprintf("Log2 protein matrix: %d x %d\n", nrow(protein_log2), ncol(protein_log2)))

# ===========================================================================
# 7. PTAU LOESS trend (fitted once)
# ===========================================================================
cat("\n========== Fitting PTAU LOESS trend ==========\n")

age_vals   <- df_analysis$AGE
ptau_log2  <- df_analysis$PTAU_log2

ptau_valid  <- !is.na(age_vals) & !is.na(ptau_log2)
cat(sprintf("Valid (AGE, PTAU) pairs: %d / %d\n", sum(ptau_valid), nrow(df_analysis)))

ptau_loess <- tryCatch(
  loess(ptau_log2[ptau_valid] ~ age_vals[ptau_valid], span = LOESS_SPAN),
  error = function(e) NULL
)

if (is.null(ptau_loess)) {
  stop("PTAU LOESS fitting failed. Cannot proceed with analysis.")
}
cat("PTAU LOESS fit: success\n")

ptau_fitted <- rep(NA_real_, nrow(df_analysis))
ptau_fitted[ptau_valid] <- fitted(ptau_loess)

ptau_fitted_z <- (ptau_fitted - mean(ptau_fitted, na.rm = TRUE)) /
                 sd(ptau_fitted, na.rm = TRUE)

ptau_age_cor <- cor.test(ptau_fitted[ptau_valid], age_vals[ptau_valid],
                         method = "spearman", exact = FALSE)
cat(sprintf("PTAU LOESS trend vs AGE: Spearman rho = %.4f, p = %.2e\n",
            ptau_age_cor$estimate, ptau_age_cor$p.value))
if (ptau_age_cor$estimate > 0) {
  cat("  -> PTAU INCREASES with age\n")
} else {
  cat("  -> PTAU DECREASES with age\n")
}

# ===========================================================================
# 8. Per-protein LOESS and trend correlation
# ===========================================================================
cat("\n========== Running per-protein LOESS trend correlation ==========\n")

n_proteins <- ncol(protein_log2)
protein_names <- colnames(protein_log2)

results <- data.frame(
  protein_id      = protein_names,
  n_valid         = NA_integer_,
  loess_corr      = NA_real_,
  loess_corr_abs  = NA_real_,
  loess_mae_z     = NA_real_,
  loess_converged = FALSE,
  effective_span  = NA_real_,
  stringsAsFactors = FALSE
)

cat(sprintf("Testing %d proteins...\n", n_proteins))
pb <- txtProgressBar(min = 0, max = n_proteins, style = 3)

for (i in seq_len(n_proteins)) {
  setTxtProgressBar(pb, i)

  prot_log2  <- protein_log2[, i]
  valid_idx  <- !is.na(age_vals) & !is.na(prot_log2)
  n_valid    <- sum(valid_idx)
  results$n_valid[i] <- n_valid

  if (n_valid >= MIN_VALID_PAIRS) {
    eff_span <- max(LOESS_SPAN, 15 / n_valid)
    results$effective_span[i] <- eff_span

    loess_fit <- tryCatch(
      loess(prot_log2[valid_idx] ~ age_vals[valid_idx], span = eff_span),
      error = function(e) NULL
    )

    if (!is.null(loess_fit)) {
      results$loess_converged[i] <- TRUE

      prot_fitted <- rep(NA_real_, length(prot_log2))
      prot_fitted[valid_idx] <- fitted(loess_fit)

      shared_valid <- which(valid_idx & !is.na(ptau_fitted))
      if (length(shared_valid) >= 10) {
        cor_test <- tryCatch(
          cor.test(ptau_fitted[shared_valid], prot_fitted[shared_valid],
                   method = "spearman", exact = FALSE),
          error = function(e) NULL
        )
        if (!is.null(cor_test)) {
          results$loess_corr[i] <- cor_test$estimate
          results$loess_corr_abs[i] <- abs(cor_test$estimate)
        }
      }

      if (!is.na(results$loess_corr[i])) {
        prot_fitted_z <- (prot_fitted - mean(prot_fitted, na.rm = TRUE)) /
                         sd(prot_fitted, na.rm = TRUE)
        shared_z <- which(shared_valid %in% which(!is.na(prot_fitted_z)))
        if (length(shared_z) >= 10) {
          results$loess_mae_z[i] <- mean(abs(
            ptau_fitted_z[shared_z] - prot_fitted_z[shared_z]
          ), na.rm = TRUE)
        }
      }
    }
  }
}
close(pb)
cat("\n")

# ===========================================================================
# 9. Post-process results
# ===========================================================================
cat("\n========== Post-processing results ==========\n")

results_valid  <- results[!is.na(results$loess_corr), ]
results_failed <- results[is.na(results$loess_corr), ]

cat(sprintf("Proteins successfully ranked: %d / %d\n",
            nrow(results_valid), n_proteins))
cat(sprintf("Proteins with failed LOESS/correlation: %d\n", nrow(results_failed)))

if (nrow(results_failed) > 0) {
  n_low_n      <- sum(results_failed$n_valid < MIN_VALID_PAIRS, na.rm = TRUE)
  n_no_loess   <- sum(results_failed$n_valid >= MIN_VALID_PAIRS &
                      !results_failed$loess_converged, na.rm = TRUE)
  cat(sprintf("  Insufficient valid pairs (< %d): %d\n", MIN_VALID_PAIRS, n_low_n))
  cat(sprintf("  LOESS failed to converge: %d\n", n_no_loess))
}

results_valid <- results_valid[order(results_valid$loess_corr, decreasing = TRUE), ]
results_valid$rank <- seq_len(nrow(results_valid))

cat(sprintf("\nloess_corr distribution:\n"))
cat(sprintf("  Min:    %.4f\n", min(results_valid$loess_corr)))
cat(sprintf("  Median: %.4f\n", median(results_valid$loess_corr)))
cat(sprintf("  Max:    %.4f\n", max(results_valid$loess_corr)))

# ===========================================================================
# 10. Annotate with gene symbols
# ===========================================================================
cat("\n========== Annotating with gene symbols ==========\n")

annot_lookup <- dict[, c("Analytes", "EntrezGeneSymbol", "TargetFullName", "Target", "UniProt")]
colnames(annot_lookup)[1] <- "protein_id"

results_valid <- merge(results_valid, annot_lookup, by = "protein_id",
                        all.x = TRUE, sort = FALSE)
results_valid <- results_valid[order(results_valid$loess_corr, decreasing = TRUE), ]
results_valid$rank <- seq_len(nrow(results_valid))

n_annotated <- sum(!is.na(results_valid$EntrezGeneSymbol) &
                   results_valid$EntrezGeneSymbol != "")
cat(sprintf("Proteins with gene symbol: %d / %d\n", n_annotated, nrow(results_valid)))

# ===========================================================================
# 11. Helper: gene label with duplicate handling
# ===========================================================================
get_gene_label <- function(df) {
  label <- ifelse(is.na(df$EntrezGeneSymbol) | df$EntrezGeneSymbol == "",
                  df$protein_id, df$EntrezGeneSymbol)
  dup_genes <- label[duplicated(label) | duplicated(label, fromLast = TRUE)]
  for (g in unique(dup_genes)) {
    idx <- which(label == g)
    label[idx] <- paste0(g, " (", df$protein_id[idx], ")")
  }
  return(label)
}

# ===========================================================================
# 12. Print top results
# ===========================================================================
cat("\n========== Top 20 Proteins (best tracking PTAU age-trend in AD spectrum) ==========\n")
top20 <- head(results_valid, TOP_N_EXPORT)
top20$gene_label <- get_gene_label(top20)
print(data.frame(
  Rank        = top20$rank,
  Gene        = top20$EntrezGeneSymbol,
  Target      = top20$Target,
  Trend_Corr  = round(top20$loess_corr, 4),
  MAE_z       = round(top20$loess_mae_z, 4),
  N           = top20$n_valid,
  Eff_Span    = round(top20$effective_span, 2)
))

cat("\n========== Bottom 10 Proteins (best INVERSE tracking) ==========\n")
bottom10 <- tail(results_valid, 10)
bottom10 <- bottom10[order(bottom10$loess_corr), ]
bottom10$gene_label <- get_gene_label(bottom10)
print(data.frame(
  Rank        = bottom10$rank,
  Gene        = bottom10$EntrezGeneSymbol,
  Target      = bottom10$Target,
  Trend_Corr  = round(bottom10$loess_corr, 4),
  MAE_z       = round(bottom10$loess_mae_z, 4),
  N           = bottom10$n_valid
))

# ===========================================================================
# 13. Save CSVs
# ===========================================================================
cat("\n========== Saving CSVs ==========\n")

write.csv(results_valid,
          file.path(out_dir, "all_proteins_age_trend_ranking_AD.csv"),
          row.names = FALSE)
cat(sprintf("Saved: all_proteins_age_trend_ranking_AD.csv (%d rows)\n", nrow(results_valid)))

write.csv(head(results_valid, TOP_N_EXPORT),
          file.path(out_dir, "top20_best_tracking_proteins_AD.csv"),
          row.names = FALSE)
cat(sprintf("Saved: top20_best_tracking_proteins_AD.csv\n"))

if (nrow(results_failed) > 0) {
  results_failed$failure_reason <- ifelse(
    results_failed$n_valid < MIN_VALID_PAIRS,
    sprintf("Insufficient valid pairs (%d < %d)", results_failed$n_valid, MIN_VALID_PAIRS),
    ifelse(!results_failed$loess_converged,
           "LOESS failed to converge",
           "Spearman correlation returned NA")
  )
  write.csv(results_failed,
            file.path(out_dir, "removed_proteins_loess_failure_AD.csv"),
            row.names = FALSE)
  cat(sprintf("Saved: removed_proteins_loess_failure_AD.csv (%d rows)\n", nrow(results_failed)))
}

# ===========================================================================
# 14. Combined plot: PTAU + top N proteins (z-score normalized)
# ===========================================================================
cat("\n========== Generating combined trend plot ==========\n")

top_plot <- head(results_valid, TOP_N_PLOT)
top_plot$gene_label <- get_gene_label(top_plot)

ptau_raw_z <- (ptau_log2 - mean(ptau_log2, na.rm = TRUE)) /
              sd(ptau_log2, na.rm = TRUE)

scatter_data <- data.frame(
  AGE      = age_vals,
  value_z  = ptau_raw_z,
  series   = "PTAU (raw)"
)

line_data_list <- list()
line_data_list[[1]] <- data.frame(
  AGE      = age_vals,
  fitted_z = ptau_fitted_z,
  series   = "PTAU",
  rank     = 0,
  corr_val = NA_real_
)

for (k in seq_len(nrow(top_plot))) {
  prot_id    <- top_plot$protein_id[k]
  prot_idx   <- match(prot_id, colnames(protein_log2))
  prot_log2v <- protein_log2[, prot_idx]
  valid_idx  <- !is.na(age_vals) & !is.na(prot_log2v)

  if (sum(valid_idx) < MIN_VALID_PAIRS) next

  eff_span <- max(LOESS_SPAN, 15 / sum(valid_idx))
  prot_fit <- tryCatch(
    loess(prot_log2v[valid_idx] ~ age_vals[valid_idx], span = eff_span),
    error = function(e) NULL
  )
  if (is.null(prot_fit)) next

  prot_fitted <- rep(NA_real_, length(prot_log2v))
  prot_fitted[valid_idx] <- fitted(prot_fit)
  prot_fitted_z <- (prot_fitted - mean(prot_fitted, na.rm = TRUE)) /
                   sd(prot_fitted, na.rm = TRUE)

  gene_lbl <- top_plot$gene_label[k]
  corr_val <- top_plot$loess_corr[k]
  line_data_list[[k + 1]] <- data.frame(
    AGE      = age_vals,
    fitted_z = prot_fitted_z,
    series   = sprintf("%d. %s", k, gene_lbl),
    rank     = k,
    corr_val = corr_val
  )
}
plot_lines <- do.call(rbind, line_data_list)
plot_lines$series <- factor(plot_lines$series, levels = unique(plot_lines$series))

n_prot_lines <- length(unique(plot_lines$series)) - 1
if (n_prot_lines > 0) {
  prot_colors <- hcl.colors(n_prot_lines, "Dynamic")
  names(prot_colors) <- setdiff(unique(plot_lines$series), "PTAU")
  all_colors <- c("PTAU" = "black", prot_colors)
} else {
  all_colors <- c("PTAU" = "black")
}

legend_labels <- setNames(as.character(plot_lines$series), plot_lines$series)
for (s in unique(plot_lines$series)) {
  if (s != "PTAU") {
    rv <- plot_lines$corr_val[plot_lines$series == s][1]
    if (!is.na(rv)) {
      legend_labels[s] <- sprintf("%s  (ρ=%.3f)", s, rv)
    }
  }
}

combined_p <- ggplot() +
  geom_point(data = scatter_data,
             aes(x = AGE, y = value_z),
             alpha = 0.12, size = 0.6, color = "grey50") +
  geom_line(data = plot_lines,
            aes(x = AGE, y = fitted_z, color = series, linewidth = series)) +
  scale_color_manual(values = all_colors, labels = legend_labels) +
  scale_linewidth_manual(values = setNames(
    c(1.4, rep(0.8, n_prot_lines)),
    c("PTAU", setdiff(unique(plot_lines$series), "PTAU"))
  ), guide = "none") +
  labs(
    title = "Proteins Whose Age-Trajectory Best Tracks PTAU (AD Spectrum)",
    subtitle = sprintf(
      "LOESS span = %.2f | Z-scored | %d subjects (EMCI+LMCI+AD) | Top %d of %d proteins",
      LOESS_SPAN, nrow(df_analysis), TOP_N_PLOT, nrow(results_valid)
    ),
    x = "Age (years)",
    y = "Z-score of LOESS Fitted Values",
    color = "Series (ranked by\ntrend correlation)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.text        = element_text(size = 8),
    legend.title       = element_text(size = 9),
    plot.title         = element_text(face = "bold"),
    plot.subtitle      = element_text(size = 7.5, color = "grey40"),
    legend.position    = "right",
    legend.key.height  = unit(0.35, "cm"),
    legend.key.width   = unit(1.2, "cm")
  ) +
  guides(color = guide_legend(ncol = 1, byrow = TRUE))

ggsave(file.path(out_dir, "combined_ptau_top10_age_trend_AD.pdf"),
       combined_p, width = 10, height = 7, device = "pdf")
ggsave(file.path(out_dir, "combined_ptau_top10_age_trend_AD.png"),
       combined_p, width = 10, height = 7, dpi = 300, device = "png")
cat("Saved: combined_ptau_top10_age_trend_AD.pdf/png\n")

# ===========================================================================
# 15. Individual trend plots for top N proteins
# ===========================================================================
cat("\n========== Generating individual trend plots ==========\n")

indiv_plots <- list()

for (k in seq_len(min(TOP_N_PLOT, nrow(top_plot)))) {
  prot_id    <- top_plot$protein_id[k]
  prot_idx   <- match(prot_id, colnames(protein_log2))
  prot_log2v <- protein_log2[, prot_idx]
  valid_idx  <- !is.na(age_vals) & !is.na(prot_log2v)

  if (sum(valid_idx) < MIN_VALID_PAIRS) next

  eff_span <- max(LOESS_SPAN, 15 / sum(valid_idx))
  prot_fit <- tryCatch(
    loess(prot_log2v[valid_idx] ~ age_vals[valid_idx], span = eff_span),
    error = function(e) NULL
  )
  if (is.null(prot_fit)) next

  prot_fitted <- rep(NA_real_, length(prot_log2v))
  prot_fitted[valid_idx] <- fitted(prot_fit)

  gene_lbl <- top_plot$gene_label[k]
  corr_val <- top_plot$loess_corr[k]
  mae_val  <- top_plot$loess_mae_z[k]

  plot_df <- data.frame(
    AGE           = age_vals,
    prot_log2     = prot_log2v,
    prot_fitted   = prot_fitted
  )

  p <- ggplot() +
    geom_point(data = plot_df,
               aes(x = AGE, y = prot_log2),
               alpha = 0.3, size = 0.8, color = "#2c7bb6") +
    geom_line(data = plot_df,
              aes(x = AGE, y = prot_fitted),
              color = "#2c7bb6", linewidth = 1.2) +
    labs(
      title   = sprintf("Rank #%d: %s", k, gene_lbl),
      subtitle = sprintf(
        "Trend corr ρ = %.4f | MAE_z = %.4f | N = %d | %s | AD spectrum",
        corr_val, mae_val, top_plot$n_valid[k], top_plot$Target[k]
      ),
      x = "Age (years)",
      y = sprintf("log2(%s)", gene_lbl)
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title    = element_text(face = "bold"),
      plot.subtitle = element_text(size = 8, color = "grey40")
    )

  indiv_plots[[k]] <- p

  safe_name <- gsub("[^A-Za-z0-9_]", "_", gene_lbl)
  ggsave(file.path(out_dir, sprintf("indiv_trend_%02d_%s_AD.pdf", k, safe_name)),
         p, width = 7, height = 5, device = "pdf")
  ggsave(file.path(out_dir, sprintf("indiv_trend_%02d_%s_AD.png", k, safe_name)),
         p, width = 7, height = 5, dpi = 300, device = "png")
}

if (length(indiv_plots) > 0) {
  n_indiv <- length(indiv_plots)
  ncol <- if (n_indiv <= 5) 1 else 2
  multi_p <- wrap_plots(indiv_plots, ncol = ncol) +
    plot_annotation(
      title   = "Top Proteins Tracking PTAU Age-Trend — AD Spectrum (EMCI+LMCI+AD)",
      subtitle = sprintf("log2(protein) vs AGE with LOESS trend | %d subjects",
                         nrow(df_analysis)),
      theme = theme(plot.title = element_text(face = "bold", size = 14))
    )

  multi_w <- ifelse(ncol == 1, 8, 14)
  multi_h <- ceiling(n_indiv / ncol) * 5

  ggsave(file.path(out_dir, "indiv_trends_multi_panel_AD.pdf"),
         multi_p, width = multi_w, height = multi_h, device = "pdf",
         limitsize = FALSE)
  ggsave(file.path(out_dir, "indiv_trends_multi_panel_AD.png"),
         multi_p, width = multi_w, height = multi_h, dpi = 300, device = "png",
         limitsize = FALSE)
  cat(sprintf("Saved: %d individual trend plots + multi-panel figure\n", n_indiv))
}

# ===========================================================================
# 16. Analysis summary
# ===========================================================================
cat("\n")
cat("========== Analysis Summary ==========\n")
cat(sprintf("Population: AD spectrum (EMCI + LMCI + AD)\n"))
cat(sprintf("Input:  %d baseline samples (EMCI=%d, LMCI=%d, AD=%d)\n",
            nrow(df_analysis),
            sum(df_analysis$DX == "EMCI"),
            sum(df_analysis$DX == "LMCI"),
            sum(df_analysis$DX == "AD")))
cat(sprintf("AGE range: %.1f - %.1f years (median %.1f)\n",
            min(df_analysis$AGE), max(df_analysis$AGE), median(df_analysis$AGE)))
cat(sprintf("Proteins identified: %d\n", length(protein_cols)))
cat(sprintf("  Removed (>20%% NA): %d\n", length(proteins_high_missing)))
cat(sprintf("  Kept for analysis: %d\n", length(proteins_kept)))
cat(sprintf("  Successfully ranked: %d\n", nrow(results_valid)))
cat(sprintf("  LOESS/correlation failed: %d\n", nrow(results_failed)))
cat(sprintf("\nPTAU vs AGE trend:\n"))
cat(sprintf("  Spearman rho = %.4f (p = %.2e)\n",
            ptau_age_cor$estimate, ptau_age_cor$p.value))
cat(sprintf("  Direction: %s\n",
            ifelse(ptau_age_cor$estimate > 0, "INCREASING with age", "DECREASING with age")))
cat(sprintf("\nTop 5 best-tracking proteins:\n"))
for (i in 1:min(5, nrow(results_valid))) {
  g <- results_valid$EntrezGeneSymbol[i]
  t <- results_valid$Target[i]
  r <- results_valid$loess_corr[i]
  m <- results_valid$loess_mae_z[i]
  cat(sprintf("  #%d: %s (%s)  rho=%.4f  MAE_z=%.4f\n", i, g, t, r, m))
}
cat(sprintf("\nOutput directory: %s\n", normalizePath(out_dir)))
cat("\n========== Done ==========\n")
