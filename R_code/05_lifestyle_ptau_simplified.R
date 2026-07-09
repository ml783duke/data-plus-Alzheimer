###############################################################################
# 05_lifestyle_ptau_simplified.R
# Simplified lifestyle-PTAU analysis: 7 key variables
# - Univariate: log2(PTAU) ~ lifestyle + AGE + SEX (each separately, FDR)
# - Multivariable: log2(PTAU) ~ all 7 + AGE + SEX
###############################################################################

# ===========================================================================
# 0. Setup
# ===========================================================================
library(readxl)
library(dplyr)
library(ggplot2)

set.seed(42)

out_dir <- "output/lifestyle_ptau_simplified"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# 7 variables selected
LIFESTYLE_VARS <- c("DHA", "EPA", "HCys", "NPIK", "NPIKTOT", "MH14ALCH", "MH16SMOK")

# Labels
var_labels <- c(
  "DHA"       = "DHA (n-3 PUFA)",
  "EPA"       = "EPA (n-3 PUFA)",
  "HCys"      = "Homocysteine",
  "NPIK"      = "NPI Sleep Domain",
  "NPIKTOT"   = "NPI Total Score",
  "MH14ALCH"  = "Alcohol Abuse History",
  "MH16SMOK"  = "Smoking Status"
)

# Categories for coloring
var_cat <- c(
  "DHA"       = "Nutrition",
  "EPA"       = "Nutrition",
  "HCys"      = "Nutrition",
  "NPIK"      = "Neuropsychiatric",
  "NPIKTOT"   = "Neuropsychiatric",
  "MH14ALCH"  = "Substance Use",
  "MH16SMOK"  = "Substance Use"
)

cat("========== Simplified Lifestyle → PTAU ==========\n")
cat(sprintf("Testing %d variables\n", length(LIFESTYLE_VARS)))

# ===========================================================================
# 1. Load and filter data
# ===========================================================================
cat("\n========== Loading data ==========\n")
master <- read_excel("master_data.xlsx", sheet = "Sheet1")

df_bl <- master[master$VISCODE2 == "bl", ]
df_analysis <- df_bl[!is.na(df_bl$PTAU) & !is.na(df_bl$AGE) & !is.na(df_bl$PTGENDER), ]

df_analysis$PTAU_log2 <- log2(df_analysis$PTAU)
df_analysis$SEX <- ifelse(df_analysis$PTGENDER == "Male", 1L, 0L)

cat(sprintf("Analysis sample: %d subjects\n", nrow(df_analysis)))

# Report coverage
cat("\nVariable coverage:\n")
for (v in LIFESTYLE_VARS) {
  n <- sum(!is.na(df_analysis[[v]]))
  n_unique <- length(unique(na.omit(df_analysis[[v]])))
  cat(sprintf("  %-12s: %d / %d (%.1f%%)  unique=%d\n",
              v, n, nrow(df_analysis), 100*n/nrow(df_analysis), n_unique))
}

# ===========================================================================
# 2. Univariate models (each separately, adjusted for AGE + SEX)
# ===========================================================================
cat("\n========== Univariate models ==========\n")

uni_results <- data.frame(
  variable  = LIFESTYLE_VARS,
  n         = NA_integer_,
  beta      = NA_real_,
  se        = NA_real_,
  ci_low    = NA_real_,
  ci_high   = NA_real_,
  t_value   = NA_real_,
  p_value   = NA_real_,
  model_r2  = NA_real_,
  stringsAsFactors = FALSE
)

for (i in seq_along(LIFESTYLE_VARS)) {
  v <- LIFESTYLE_VARS[i]
  df_sub <- df_analysis[!is.na(df_analysis[[v]]), ]
  uni_results$n[i] <- nrow(df_sub)

  # Standardize lifestyle variable
  x_val <- df_sub[[v]]
  x_scaled <- as.vector(scale(x_val))

  fit <- lm(PTAU_log2 ~ x_scaled + AGE + SEX, data = df_sub)
  s <- summary(fit)

  uni_results$beta[i]      <- s$coefficients["x_scaled", "Estimate"]
  uni_results$se[i]        <- s$coefficients["x_scaled", "Std. Error"]
  uni_results$t_value[i]   <- s$coefficients["x_scaled", "t value"]
  uni_results$p_value[i]   <- s$coefficients["x_scaled", "Pr(>|t|)"]
  uni_results$model_r2[i]  <- s$r.squared
  uni_results$ci_low[i]    <- uni_results$beta[i] - 1.96 * uni_results$se[i]
  uni_results$ci_high[i]   <- uni_results$beta[i] + 1.96 * uni_results$se[i]
}

# FDR correction
uni_results$fdr <- p.adjust(uni_results$p_value, method = "BH")
uni_results$significant <- uni_results$fdr < 0.05

# Sort by p-value
uni_results <- uni_results[order(uni_results$p_value), ]

# Add labels
uni_results$label    <- var_labels[uni_results$variable]
uni_results$category <- var_cat[uni_results$variable]

cat(sprintf("FDR < 0.05: %d / %d\n", sum(uni_results$significant), nrow(uni_results)))

# ===========================================================================
# 3. Multivariable model
# ===========================================================================
cat("\n========== Multivariable model ==========\n")

# Build formula
var_list <- paste(LIFESTYLE_VARS, collapse = " + ")
formula_multi <- as.formula(paste("PTAU_log2 ~", var_list, "+ AGE + SEX"))

# Standardize all lifestyle vars for comparable effect sizes
df_multi <- df_analysis[, c("PTAU_log2", LIFESTYLE_VARS, "AGE", "SEX")]
df_multi_complete <- na.omit(df_multi)

for (v in LIFESTYLE_VARS) {
  df_multi_complete[[v]] <- as.vector(scale(df_multi_complete[[v]]))
}

fit_multi <- lm(formula_multi, data = df_multi_complete)

cat(sprintf("Complete cases: %d\n", nrow(df_multi_complete)))

# ===========================================================================
# 4. Print results
# ===========================================================================
cat("\n========== Univariate Results (FDR corrected) ==========\n")
print(data.frame(
  Variable    = uni_results$label,
  N           = uni_results$n,
  Beta        = round(uni_results$beta, 4),
  CI          = sprintf("[%.4f, %.4f]", uni_results$ci_low, uni_results$ci_high),
  P           = formatC(uni_results$p_value, format = "e", digits = 2),
  FDR         = formatC(uni_results$fdr, format = "e", digits = 2),
  Sig05       = uni_results$significant,
  stringsAsFactors = FALSE
), row.names = FALSE)

cat("\n========== Multivariable Model ==========\n")
s_multi <- summary(fit_multi)
cat(sprintf("Complete cases: %d\n", nrow(df_multi_complete)))
cat(sprintf("Model R² = %.4f, Adjusted R² = %.4f\n",
            s_multi$r.squared, s_multi$adj.r.squared))
cat(sprintf("Model F(%d, %d) = %.2f, p = %.2e\n\n",
            s_multi$fstatistic[2], s_multi$fstatistic[3],
            s_multi$fstatistic[1],
            pf(s_multi$fstatistic[1], s_multi$fstatistic[2],
               s_multi$fstatistic[3], lower.tail = FALSE)))

cat("Coefficients:\n")
multi_coefs <- as.data.frame(s_multi$coefficients[, c("Estimate", "Std. Error", "Pr(>|t|)")])
multi_coefs$Variable <- rownames(multi_coefs)
multi_coefs <- multi_coefs[multi_coefs$Variable != "(Intercept)", ]
print(multi_coefs, row.names = FALSE)

# ===========================================================================
# 5. Save tables
# ===========================================================================
cat("\n========== Saving tables ==========\n")

write.csv(uni_results,
          file.path(out_dir, "univariate_lifestyle_ptau.csv"),
          row.names = FALSE)
cat("Saved: univariate_lifestyle_ptau.csv\n")

# Multivariable coefficients with CI
multi_export <- data.frame(
  Variable    = rownames(s_multi$coefficients),
  Beta        = s_multi$coefficients[, "Estimate"],
  SE          = s_multi$coefficients[, "Std. Error"],
  t_value     = s_multi$coefficients[, "t value"],
  P_value     = s_multi$coefficients[, "Pr(>|t|)"],
  row.names   = NULL
)
multi_export$CI_low  <- multi_export$Beta - 1.96 * multi_export$SE
multi_export$CI_high <- multi_export$Beta + 1.96 * multi_export$SE

write.csv(multi_export,
          file.path(out_dir, "multivariable_lifestyle_ptau.csv"),
          row.names = FALSE)
cat("Saved: multivariable_lifestyle_ptau.csv\n")

# ===========================================================================
# 6. Forest plot (univariate results)
# ===========================================================================
cat("\n========== Generating figures ==========\n")

# Sort by beta for forest plot
plot_uni <- uni_results
plot_uni <- plot_uni[order(plot_uni$beta), ]
plot_uni$y_label <- factor(plot_uni$label, levels = plot_uni$label)

# Significance annotation
plot_uni$sig_anno <- ""
plot_uni$sig_anno[plot_uni$significant] <- "*"

p_forest <- ggplot(plot_uni, aes(x = beta, y = y_label, color = category)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", alpha = 0.6) +
  geom_point(aes(size = n), alpha = 0.9) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.15, alpha = 0.9) +
  geom_text(aes(x = ifelse(beta > 0, ci_high + 0.015, ci_low - 0.015),
                label = sig_anno),
            size = 6, color = "black", vjust = 0.5) +
  scale_color_manual(
    values = c("Nutrition" = "#2c7bb6", "Neuropsychiatric" = "#d7191c",
               "Substance Use" = "#fdae61"),
    name = ""
  ) +
  scale_size_continuous(name = "N", range = c(2.5, 5)) +
  labs(
    title = "Lifestyle Factors Associated with PTAU",
    subtitle = sprintf(
      "Univariate models: log2(PTAU) ~ variable + AGE + SEX | %d subjects | * FDR < 0.05",
      nrow(df_analysis)
    ),
    x = expression("Standardized β (change in log"[2]*" PTAU per 1 SD)"),
    y = ""
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(size = 9, color = "grey40"),
    legend.position  = "bottom",
    axis.text.y      = element_text(size = 10),
    panel.grid.major.y = element_blank()
  )

ggsave(file.path(out_dir, "lifestyle_ptau_univariate_forest.png"),
       p_forest, width = 9, height = 4.5, dpi = 300, device = "png")
cat("Saved: lifestyle_ptau_univariate_forest.png\n")

# ===========================================================================
# 7. Multivariable coefficient plot
# ===========================================================================

# Extract lifestyle + covariate coefficients from multivariable model
multi_plot_data <- multi_export[multi_export$Variable != "(Intercept)", ]
multi_plot_data$label <- ifelse(
  multi_plot_data$Variable %in% names(var_labels),
  var_labels[multi_plot_data$Variable],
  multi_plot_data$Variable
)
multi_plot_data$is_lifestyle <- multi_plot_data$Variable %in% LIFESTYLE_VARS
multi_plot_data <- multi_plot_data[order(multi_plot_data$Beta), ]
multi_plot_data$y_label <- factor(multi_plot_data$label, levels = multi_plot_data$label)

# P-value annotation
multi_plot_data$p_label <- sprintf("P = %.3f", multi_plot_data$P_value)
multi_plot_data$p_label[multi_plot_data$P_value < 0.001] <- "P < 0.001"

p_multi <- ggplot(multi_plot_data, aes(x = Beta, y = y_label,
                                        color = is_lifestyle)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", alpha = 0.6) +
  geom_point(size = 3, alpha = 0.9) +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.2, alpha = 0.9) +
  geom_text(aes(x = max(CI_high) * 0.85, label = p_label),
            size = 3.2, hjust = 1, color = "grey30") +
  scale_color_manual(values = c("TRUE" = "#2c7bb6", "FALSE" = "grey50"),
                      labels = c("TRUE" = "Lifestyle", "FALSE" = "Covariate"),
                      name = "") +
  labs(
    title = "Multivariable Model: Lifestyle Factors + PTAU",
    subtitle = sprintf(
      "log2(PTAU) ~ DHA + EPA + HCys + NPIK + NPIKTOT + MH14ALCH + MH16SMOK + AGE + SEX\n%d complete cases | Model R² = %.3f | Adjusted R² = %.3f",
      nrow(df_multi_complete), s_multi$r.squared, s_multi$adj.r.squared
    ),
    x = expression("Standardized β (change in log"[2]*" PTAU per 1 SD)"),
    y = ""
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 8, color = "grey40"),
    legend.position  = "bottom",
    axis.text.y      = element_text(size = 10),
    panel.grid.major.y = element_blank()
  )

ggsave(file.path(out_dir, "lifestyle_ptau_multivariable.png"),
       p_multi, width = 9, height = 5, dpi = 300, device = "png")
cat("Saved: lifestyle_ptau_multivariable.png\n")

# ===========================================================================
# 8. Combined table (univariate + multivariable side by side)
# ===========================================================================
cat("\n========== Combined Table ==========\n")

# Build comparison table
comp_table <- data.frame(
  Variable        = var_labels[LIFESTYLE_VARS],
  Category        = var_cat[LIFESTYLE_VARS],
  stringsAsFactors = FALSE
)
rownames(comp_table) <- LIFESTYLE_VARS

# Add univariate results
comp_table$Uni_Beta <- NA_real_
comp_table$Uni_P    <- NA_real_
comp_table$Uni_FDR  <- NA_real_
for (v in LIFESTYLE_VARS) {
  r <- uni_results[uni_results$variable == v, ]
  if (nrow(r) == 1) {
    comp_table[v, "Uni_Beta"] <- r$beta
    comp_table[v, "Uni_P"]    <- r$p_value
    comp_table[v, "Uni_FDR"]  <- r$fdr
  }
}

# Add multivariable results
comp_table$Multi_Beta <- NA_real_
comp_table$Multi_P    <- NA_real_
for (v in LIFESTYLE_VARS) {
  if (v %in% multi_export$Variable) {
    r <- multi_export[multi_export$Variable == v, ]
    comp_table[v, "Multi_Beta"] <- r$Beta
    comp_table[v, "Multi_P"]    <- r$P_value
  }
}

# Round for display
comp_table_display <- comp_table
comp_table_display$Uni_Beta  <- round(comp_table$Uni_Beta, 4)
comp_table_display$Uni_P     <- formatC(comp_table$Uni_P, format = "e", digits = 2)
comp_table_display$Uni_FDR   <- formatC(comp_table$Uni_FDR, format = "e", digits = 2)
comp_table_display$Multi_Beta <- round(comp_table$Multi_Beta, 4)
comp_table_display$Multi_P   <- formatC(comp_table$Multi_P, format = "e", digits = 2)

# Add significance markers
comp_table_display$Uni_Sig <- ""
comp_table_display$Uni_Sig[comp_table$Uni_FDR < 0.05] <- "*"

cat("\nCombined Uni- vs Multi-variable Results:\n")
print(comp_table_display[, c("Variable", "Category", "Uni_Beta", "Uni_P", "Uni_FDR",
                              "Uni_Sig", "Multi_Beta", "Multi_P")], row.names = FALSE)

write.csv(comp_table_display,
          file.path(out_dir, "lifestyle_ptau_combined_table.csv"),
          row.names = FALSE)
cat("\nSaved: lifestyle_ptau_combined_table.csv\n")

# ===========================================================================
# 9. Summary
# ===========================================================================
cat("\n")
cat("========== Analysis Summary ==========\n")
cat(sprintf("Subjects: %d (baseline + PTAU + AGE + SEX)\n", nrow(df_analysis)))
cat(sprintf("Variables tested: %d\n", length(LIFESTYLE_VARS)))
cat(sprintf("Univariate: %d / %d significant (FDR < 0.05)\n",
            sum(uni_results$significant), nrow(uni_results)))
cat(sprintf("Multivariable: %d complete cases, R² = %.4f, Adj R² = %.4f\n",
            nrow(df_multi_complete), s_multi$r.squared, s_multi$adj.r.squared))

if (sum(uni_results$significant) > 0) {
  cat("\nSignificant lifestyle factors (FDR < 0.05):\n")
  for (i in which(uni_results$significant)) {
    cat(sprintf("  %s: β = %.4f, P = %.2e, FDR = %.2e\n",
                uni_results$label[i], uni_results$beta[i],
                uni_results$p_value[i], uni_results$fdr[i]))
  }
}

cat(sprintf("\nOutput directory: %s\n", normalizePath(out_dir)))
cat("\n========== Done ==========\n")
