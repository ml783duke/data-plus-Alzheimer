###############################################################################
# 04_lifestyle_ptau.R
# Phase 1: Which lifestyle factors are directly associated with PTAU?
# - For each lifestyle variable: log2(PTAU) ~ lifestyle + AGE + SEX
# - FDR correction across all tested lifestyle variables
# - Forest plot + summary table
###############################################################################

# ===========================================================================
# 0. Setup
# ===========================================================================
library(readxl)
library(dplyr)
library(ggplot2)

set.seed(42)

out_dir <- "output/lifestyle_ptau"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Lifestyle variables to exclude (too sparse or all-missing)
EXCLUDE_VARS <- c(
  "NPIK1","NPIK2","NPIK3","NPIK4","NPIK5","NPIK6","NPIK7","NPIK8",
  "NPIK9A","NPIK9B","NPIK9C",    # NPI sub-items (<10% coverage)
  "QID46_9",                      # Social activities (2.5%)
  "QID44_1_3","QID44_2_3",       # Loneliness (0%)
  "PTWORK",                       # Occupation (0%)
  "MH14AALCH"                     # Drinks/day among abusers only (2.8%)
)

cat("========== Phase 1: Lifestyle → PTAU ==========\n")
cat(sprintf("Excluding %d sparse variables\n", length(EXCLUDE_VARS)))

# ===========================================================================
# 1. Load data
# ===========================================================================
cat("\n========== Loading data ==========\n")
master <- read_excel("master_data.xlsx", sheet = "Sheet1")
cat(sprintf("Master: %d rows x %d columns\n", nrow(master), ncol(master)))

# ===========================================================================
# 2. Filter to analysis set: baseline + PTAU + AGE + SEX
# ===========================================================================
cat("\n========== Filtering samples ==========\n")

# Baseline only
df_bl <- master[master$VISCODE2 == "bl", ]
cat(sprintf("Baseline: %d rows\n", nrow(df_bl)))

# Non-missing PTAU, AGE, and PTGENDER
df_analysis <- df_bl[!is.na(df_bl$PTAU) & !is.na(df_bl$AGE) & !is.na(df_bl$PTGENDER), ]
cat(sprintf("Baseline + PTAU + AGE + SEX: %d samples\n", nrow(df_analysis)))

# Log2 transform PTAU
df_analysis$PTAU_log2 <- log2(df_analysis$PTAU)
df_analysis$SEX <- ifelse(df_analysis$PTGENDER == "Male", 1, 0)

# ===========================================================================
# 3. Identify lifestyle variables
# ===========================================================================
cat("\n========== Identifying lifestyle variables ==========\n")

mmse_pos <- which(colnames(df_analysis) == "MMSE")
prot_start <- grep("^X[0-9]+[.][0-9]+$", colnames(df_analysis))[1]
lifestyle_cols <- colnames(df_analysis)[(mmse_pos + 1):(prot_start - 1)]
cat(sprintf("Lifestyle columns: %d (cols %d-%d)\n",
            length(lifestyle_cols), mmse_pos + 1, prot_start - 1))

# Filter to testable variables
test_vars <- setdiff(lifestyle_cols, EXCLUDE_VARS)
cat(sprintf("After excluding sparse: %d variables to test\n", length(test_vars)))

# Report coverage for each
cat("\nCoverage in analysis set:\n")
for (v in test_vars) {
  n <- sum(!is.na(df_analysis[[v]]))
  n_unique <- length(unique(na.omit(df_analysis[[v]])))
  cat(sprintf("  %-20s: %d / %d (%.1f%%)  unique=%d\n",
              v, n, nrow(df_analysis), 100*n/nrow(df_analysis), n_unique))
}

# ===========================================================================
# 4. For each lifestyle variable: log2(PTAU) ~ lifestyle + AGE + SEX
# ===========================================================================
cat("\n========== Running regressions ==========\n")

results <- data.frame(
  variable      = test_vars,
  n             = NA_integer_,
  n_unique      = NA_integer_,
  beta          = NA_real_,
  se            = NA_real_,
  t_value       = NA_real_,
  p_value       = NA_real_,
  r_squared     = NA_real_,
  adj_r_squared = NA_real_,
  stringsAsFactors = FALSE
)

for (i in seq_along(test_vars)) {
  v <- test_vars[i]
  df_sub <- df_analysis[!is.na(df_analysis[[v]]), ]
  results$n[i] <- nrow(df_sub)
  results$n_unique[i] <- length(unique(df_sub[[v]]))

  if (nrow(df_sub) >= 30 && results$n_unique[i] >= 2) {
    # Standardize lifestyle variable for comparable effect sizes
    x_val <- df_sub[[v]]
    x_sd <- sd(x_val, na.rm = TRUE)
    if (x_sd > 0) {
      x_scaled <- (x_val - mean(x_val, na.rm = TRUE)) / x_sd
    } else {
      x_scaled <- x_val
    }

    fit <- tryCatch(
      lm(PTAU_log2 ~ x_scaled + AGE + SEX, data = df_sub),
      error = function(e) NULL
    )

    if (!is.null(fit)) {
      s <- summary(fit)
      coefs <- s$coefficients
      if ("x_scaled" %in% rownames(coefs)) {
        results$beta[i]    <- coefs["x_scaled", "Estimate"]
        results$se[i]      <- coefs["x_scaled", "Std. Error"]
        results$t_value[i] <- coefs["x_scaled", "t value"]
        results$p_value[i] <- coefs["x_scaled", "Pr(>|t|)"]
      }
      results$r_squared[i]     <- s$r.squared
      results$adj_r_squared[i] <- s$adj.r.squared
    }
  }
}

# Remove untestable
results <- results[!is.na(results$p_value), ]
cat(sprintf("Successfully tested: %d / %d lifestyle variables\n",
            nrow(results), length(test_vars)))

# ===========================================================================
# 5. FDR correction
# ===========================================================================
cat("\n========== FDR correction ==========\n")

results$fdr <- p.adjust(results$p_value, method = "BH")
results$significant_fdr05 <- results$fdr < 0.05
results$significant_fdr10 <- results$fdr < 0.10

n_sig_05 <- sum(results$significant_fdr05)
n_sig_10 <- sum(results$significant_fdr10)
cat(sprintf("Significant at FDR < 0.05: %d\n", n_sig_05))
cat(sprintf("Significant at FDR < 0.10: %d\n", n_sig_10))

# Sort by p-value
results <- results[order(results$p_value), ]

# ===========================================================================
# 6. Add variable labels
# ===========================================================================
var_labels <- c(
  "DHA"           = "DHA (n-3 PUFA)",
  "EPA"           = "EPA (n-3 PUFA)",
  "HCys"          = "Homocysteine",
  "BSH_FA20.5_w3" = "EPA oxylipin (HPH)",
  "BSH_FA22.6_w3" = "DHA oxylipin (HPH)",
  "BSL_FA20.5_w3" = "EPA oxylipin (LPH)",
  "BSL_FA22.6_w3" = "DHA oxylipin (LPH)",
  "X424"          = "Palmitate (16:0)",
  "X439"          = "Stearate (18:0)",
  "X100008930"    = "Oleate (18:1)",
  "X2050"         = "Linoleate (18:2)",
  "X100000665"    = "Arachidonate (20:4)",
  "X100001181"    = "DHA (22:6) metabolite",
  "X848"          = "Cotinine",
  "X100002717"    = "Hydroxycotinine 1",
  "X100002719"    = "Hydroxycotinine 2",
  "X100004494"    = "Nicotine metabolite",
  "NPIK"          = "NPI Sleep (NPIK)",
  "NPIKTOT"       = "NPI Total Score",
  "NPIQ_K"        = "NPIQ Sleep (NPIQ-K)",
  "NPIQ_KSEV"     = "NPIQ Sleep Severity",
  "MH14ALCH"      = "Alcohol Abuse History",
  "MH16SMOK"      = "Smoking Status",
  "MH16ASMOK"     = "Smoking Packs/Day",
  "MH16BSMOK"     = "Smoking Years",
  "MH16CSMOK"     = "Smoking Pack-Years",
  "PTWORKHS"      = "Work Hours/Week"
)
results$label <- var_labels[results$variable]
results$label[is.na(results$label)] <- results$variable[is.na(results$label)]

# Categorize
results$category <- dplyr::case_when(
  results$variable %in% c("DHA","EPA","HCys","BSH_FA20.5_w3","BSH_FA22.6_w3",
                           "BSL_FA20.5_w3","BSL_FA22.6_w3",
                           "X424","X439","X100008930","X2050","X100000665","X100001181")
    ~ "Nutrition",
  results$variable %in% c("X848","X100002717","X100002719","X100004494",
                           "MH16SMOK","MH16ASMOK","MH16BSMOK","MH16CSMOK")
    ~ "Smoking",
  results$variable %in% c("MH14ALCH") ~ "Alcohol",
  results$variable %in% c("NPIK","NPIKTOT","NPIK1","NPIK2","NPIK3","NPIK4",
                           "NPIK5","NPIK6","NPIK7","NPIK8","NPIK9A","NPIK9B",
                           "NPIK9C","NPIQ_K","NPIQ_KSEV")
    ~ "Sleep",
  results$variable %in% c("PTWORKHS") ~ "Physical Activity",
  TRUE ~ "Social"
)

# ===========================================================================
# 7. Print results
# ===========================================================================
cat("\n========== Lifestyle → PTAU Results (FDR-corrected) ==========\n")
print(data.frame(
  Variable    = results$label,
  Category    = results$category,
  N           = results$n,
  Beta        = round(results$beta, 4),
  SE          = round(results$se, 4),
  P           = formatC(results$p_value, format = "e", digits = 2),
  FDR         = formatC(results$fdr, format = "e", digits = 2),
  Sig05       = results$significant_fdr05,
  stringsAsFactors = FALSE
), row.names = FALSE)

# ===========================================================================
# 8. Save CSV
# ===========================================================================
write.csv(results, file.path(out_dir, "lifestyle_ptau_associations.csv"), row.names = FALSE)
cat(sprintf("\nSaved: lifestyle_ptau_associations.csv\n"))

# ===========================================================================
# 9. Forest plot
# ===========================================================================
cat("\n========== Generating forest plot ==========\n")

plot_data <- results
plot_data <- plot_data[order(plot_data$beta), ]
plot_data$var_label <- factor(plot_data$label, levels = plot_data$label)

# Color by FDR significance
plot_data$sig_level <- "FDR ≥ 0.10"
plot_data$sig_level[plot_data$significant_fdr10] <- "FDR < 0.10"
plot_data$sig_level[plot_data$significant_fdr05] <- "FDR < 0.05"
plot_data$sig_level <- factor(plot_data$sig_level,
                               levels = c("FDR ≥ 0.10", "FDR < 0.10", "FDR < 0.05"))

# Compute CI
plot_data$ci_low  <- plot_data$beta - 1.96 * plot_data$se
plot_data$ci_high <- plot_data$beta + 1.96 * plot_data$se

p_forest <- ggplot(plot_data, aes(x = beta, y = var_label, color = sig_level)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", alpha = 0.6) +
  geom_point(aes(size = n), alpha = 0.8) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.2, alpha = 0.8) +
  scale_color_manual(
    values = c("FDR ≥ 0.10" = "grey60", "FDR < 0.10" = "#d95f02", "FDR < 0.05" = "#1b9e77"),
    name = "Significance"
  ) +
  scale_size_continuous(name = "N", range = c(1.5, 4)) +
  labs(
    title = "Lifestyle Factors Associated with log2(PTAU)",
    subtitle = sprintf(
      "Adjusted for AGE + SEX | %d subjects | Beta = per 1-SD change in lifestyle variable | %d variables tested",
      nrow(df_analysis), nrow(results)
    ),
    x = expression("Standardized β (change in log"[2]*" PTAU per 1-SD)"),
    y = ""
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(size = 7.5, color = "grey40"),
    legend.position = "bottom",
    axis.text.y  = element_text(size = 8)
  )

ggsave(file.path(out_dir, "lifestyle_ptau_forest.pdf"), p_forest,
       width = 9, height = 7, device = "pdf")
ggsave(file.path(out_dir, "lifestyle_ptau_forest.png"), p_forest,
       width = 9, height = 7, dpi = 300, device = "png")
cat("Saved: lifestyle_ptau_forest.pdf/png\n")

# ===========================================================================
# 10. Summary
# ===========================================================================
cat("\n")
cat("========== Analysis Summary ==========\n")
cat(sprintf("Subjects: %d baseline + PTAU + AGE + SEX\n", nrow(df_analysis)))
cat(sprintf("Lifestyle variables tested: %d\n", nrow(results)))
cat(sprintf("  Excluded (too sparse): %d\n", length(EXCLUDE_VARS)))
cat(sprintf("FDR < 0.05: %d\n", n_sig_05))
cat(sprintf("FDR < 0.10: %d\n", n_sig_10))

if (n_sig_05 > 0) {
  sig_vars <- results$label[results$significant_fdr05]
  cat(sprintf("\nSignificant lifestyle factors (FDR < 0.05):\n"))
  for (i in which(results$significant_fdr05)) {
    cat(sprintf("  %s | β=%.4f | P=%.2e | FDR=%.2e\n",
                results$label[i], results$beta[i],
                results$p_value[i], results$fdr[i]))
  }
} else if (n_sig_10 > 0) {
  sig_vars <- results$label[results$significant_fdr10]
  cat(sprintf("\nSuggestive lifestyle factors (FDR < 0.10):\n"))
  for (i in which(results$significant_fdr10)) {
    cat(sprintf("  %s | β=%.4f | P=%.2e | FDR=%.2e\n",
                results$label[i], results$beta[i],
                results$p_value[i], results$fdr[i]))
  }
} else {
  cat("\nNo lifestyle variable reached FDR < 0.10 significance.\n")
  cat("Top associations (by raw p-value):\n")
  top5 <- head(results, 5)
  for (i in 1:nrow(top5)) {
    cat(sprintf("  %s | β=%.4f | P=%.2e\n",
                top5$label[i], top5$beta[i], top5$p_value[i]))
  }
}

cat(sprintf("\nOutput directory: %s\n", normalizePath(out_dir)))
cat("\n========== Done ==========\n")
