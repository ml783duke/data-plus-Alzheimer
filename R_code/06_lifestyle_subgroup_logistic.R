###############################################################################
# 06_lifestyle_subgroup_logistic.R
# Lifestyle → Subgroup: Logistic regression analysis
# - Univariate: Subgroup ~ lifestyle_var (each separately, FDR corrected)
# - Multivariable: Subgroup ~ all 7 lifestyle variables
# - Variable types handled correctly:
#   * Continuous (DHA, EPA, HCys, NPIKTOT): scaled → OR per 1-SD
#   * Binary    (NPIK, MH14ALCH, MH16SMOK): 0/1 → OR for presence vs absence
###############################################################################

# ===========================================================================
# 0. Setup
# ===========================================================================
library(readxl)
library(dplyr)
library(ggplot2)

set.seed(42)

out_dir <- "output/lifestyle_subgroup_logistic"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Variable definitions with types
LIFESTYLE_VARS <- c("DHA", "EPA", "HCys", "NPIK", "NPIKTOT", "MH14ALCH", "MH16SMOK")

# Variable type: "continuous" (scale → per 1-SD) or "binary" (0/1 → per unit)
VAR_TYPE <- c(
  "DHA"       = "continuous",
  "EPA"       = "continuous",
  "HCys"      = "continuous",
  "NPIK"      = "binary",
  "NPIKTOT"   = "continuous",   # ordinal count 0-8, treated as continuous per clinical convention
  "MH14ALCH"  = "binary",
  "MH16SMOK"  = "binary"
)

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

# OR interpretation per variable type
or_unit <- c(
  "DHA"       = "per 1-SD",
  "EPA"       = "per 1-SD",
  "HCys"      = "per 1-SD",
  "NPIK"      = "1 vs 0",
  "NPIKTOT"   = "per 1-SD",
  "MH14ALCH"  = "1 vs 0",
  "MH16SMOK"  = "1 vs 0"
)

cat("========== Lifestyle → Subgroup (Logistic) ==========\n")
cat(sprintf("Testing %d variables (%d continuous, %d binary)\n",
            length(LIFESTYLE_VARS),
            sum(VAR_TYPE == "continuous"),
            sum(VAR_TYPE == "binary")))

# ===========================================================================
# 1. Load and filter data
# ===========================================================================
cat("\n========== Loading data ==========\n")
master <- read_excel("master_data.xlsx", sheet = "Sheet1")

# Baseline only + non-missing Subgroup
df_bl <- master[master$VISCODE2 == "bl", ]
df_analysis <- df_bl[!is.na(df_bl$Subgroup), ]

# Ensure Subgroup is binary (0/1)
cat(sprintf("Subgroup distribution:\n"))
cat(sprintf("  0: %d\n", sum(df_analysis$Subgroup == 0)))
cat(sprintf("  1: %d\n", sum(df_analysis$Subgroup == 1)))
cat(sprintf("  Total: %d\n", nrow(df_analysis)))

# Report coverage
cat("\nVariable coverage (among those with Subgroup data):\n")
for (v in LIFESTYLE_VARS) {
  n <- sum(!is.na(df_analysis[[v]]))
  n_unique <- length(unique(na.omit(df_analysis[[v]])))
  cat(sprintf("  %-12s [%s]: %d / %d (%.1f%%)  unique=%d\n",
              v, VAR_TYPE[v], n, nrow(df_analysis), 100*n/nrow(df_analysis), n_unique))
}

# ===========================================================================
# 2. Univariate logistic models (each separately)
# ===========================================================================
cat("\n========== Univariate logistic models ==========\n")

uni_results <- data.frame(
  variable    = LIFESTYLE_VARS,
  var_type    = VAR_TYPE[LIFESTYLE_VARS],
  n           = NA_integer_,
  n_cases     = NA_integer_,
  odds_ratio  = NA_real_,
  or_ci_low   = NA_real_,
  or_ci_high  = NA_real_,
  beta        = NA_real_,
  se          = NA_real_,
  z_value     = NA_real_,
  p_value     = NA_real_,
  aic         = NA_real_,
  stringsAsFactors = FALSE
)

for (i in seq_along(LIFESTYLE_VARS)) {
  v <- LIFESTYLE_VARS[i]
  vtype <- VAR_TYPE[v]
  df_sub <- df_analysis[!is.na(df_analysis[[v]]), ]
  uni_results$n[i] <- nrow(df_sub)
  uni_results$n_cases[i] <- sum(df_sub$Subgroup == 1)

  x_val <- df_sub[[v]]

  if (vtype == "continuous") {
    # Scale: OR per 1-SD increase
    x_sd <- sd(x_val, na.rm = TRUE)
    if (x_sd > 0) {
      x_model <- as.vector((x_val - mean(x_val, na.rm = TRUE)) / x_sd)
    } else {
      x_model <- x_val
    }
  } else {
    # Binary: keep 0/1, OR for presence (1) vs absence (0)
    x_model <- x_val
  }

  fit <- tryCatch(
    glm(Subgroup ~ x_model, data = df_sub, family = binomial),
    error = function(e) NULL
  )

  if (!is.null(fit)) {
    s <- summary(fit)
    coefs <- s$coefficients
    if ("x_model" %in% rownames(coefs)) {
      uni_results$beta[i]    <- coefs["x_model", "Estimate"]
      uni_results$se[i]      <- coefs["x_model", "Std. Error"]
      uni_results$z_value[i] <- coefs["x_model", "z value"]
      uni_results$p_value[i] <- coefs["x_model", "Pr(>|z|)"]
      uni_results$odds_ratio[i] <- exp(uni_results$beta[i])
      uni_results$or_ci_low[i]  <- exp(uni_results$beta[i] - 1.96 * uni_results$se[i])
      uni_results$or_ci_high[i] <- exp(uni_results$beta[i] + 1.96 * uni_results$se[i])
    }
    uni_results$aic[i] <- fit$aic
  }
}

# FDR correction
uni_results <- uni_results[!is.na(uni_results$p_value), ]
uni_results$fdr <- p.adjust(uni_results$p_value, method = "BH")
uni_results$significant <- uni_results$fdr < 0.05

# Sort by p-value
uni_results <- uni_results[order(uni_results$p_value), ]

# Add labels
uni_results$label    <- var_labels[uni_results$variable]
uni_results$category <- var_cat[uni_results$variable]
uni_results$or_unit  <- or_unit[uni_results$variable]

cat(sprintf("Successfully tested: %d / %d variables\n", nrow(uni_results), length(LIFESTYLE_VARS)))
cat(sprintf("FDR < 0.05: %d / %d\n", sum(uni_results$significant), nrow(uni_results)))

# ===========================================================================
# 3. Multivariable logistic model
# ===========================================================================
cat("\n========== Multivariable logistic model ==========\n")

# Build formula
var_list <- paste(LIFESTYLE_VARS, collapse = " + ")
formula_multi <- as.formula(paste("Subgroup ~", var_list))

# Prepare data: scale continuous, keep binary as-is
df_multi <- df_analysis[, c("Subgroup", LIFESTYLE_VARS)]
df_multi_complete <- na.omit(df_multi)

# Scale continuous variables only; keep binary as 0/1
for (v in LIFESTYLE_VARS) {
  if (VAR_TYPE[v] == "continuous") {
    df_multi_complete[[v]] <- as.vector(scale(df_multi_complete[[v]]))
  }
  # binary variables stay as 0/1
}

fit_multi <- glm(formula_multi, data = df_multi_complete, family = binomial)

cat(sprintf("Complete cases: %d\n", nrow(df_multi_complete)))
cat(sprintf("Cases (Subgroup=1): %d\n", sum(df_multi_complete$Subgroup == 1)))

# ===========================================================================
# 4. Print results
# ===========================================================================
cat("\n========== Univariate Results (FDR corrected) ==========\n")
print(data.frame(
  Variable    = uni_results$label,
  Type        = uni_results$var_type,
  N           = uni_results$n,
  Cases       = uni_results$n_cases,
  OR          = sprintf("%.3f", uni_results$odds_ratio),
  CI95        = sprintf("[%.3f, %.3f]", uni_results$or_ci_low, uni_results$or_ci_high),
  Unit        = uni_results$or_unit,
  P           = formatC(uni_results$p_value, format = "e", digits = 2),
  FDR         = formatC(uni_results$fdr, format = "e", digits = 2),
  Sig05       = uni_results$significant,
  stringsAsFactors = FALSE
), row.names = FALSE)

cat("\n========== Multivariable Model ==========\n")
s_multi <- summary(fit_multi)
cat(sprintf("Complete cases: %d (%d subgroup=1)\n",
            nrow(df_multi_complete), sum(df_multi_complete$Subgroup == 1)))
cat(sprintf("Null deviance: %.2f  |  Residual deviance: %.2f\n",
            s_multi$null.deviance, s_multi$deviance))
cat(sprintf("AIC: %.2f\n", fit_multi$aic))

cat("\nCoefficients:\n")
multi_coefs <- as.data.frame(s_multi$coefficients[, c("Estimate", "Std. Error", "z value", "Pr(>|z|)")])
multi_coefs$Variable <- rownames(multi_coefs)
multi_coefs <- multi_coefs[multi_coefs$Variable != "(Intercept)", ]
multi_coefs$OR <- exp(multi_coefs$Estimate)
multi_coefs$CI_low  <- exp(multi_coefs$Estimate - 1.96 * multi_coefs[["Std. Error"]])
multi_coefs$CI_high <- exp(multi_coefs$Estimate + 1.96 * multi_coefs[["Std. Error"]])
multi_coefs$Unit <- ifelse(multi_coefs$Variable %in% names(or_unit),
                            or_unit[multi_coefs$Variable], "")
print(multi_coefs, row.names = FALSE)

# ===========================================================================
# 5. Save tables
# ===========================================================================
cat("\n========== Saving tables ==========\n")

write.csv(uni_results,
          file.path(out_dir, "univariate_lifestyle_subgroup.csv"),
          row.names = FALSE)
cat("Saved: univariate_lifestyle_subgroup.csv\n")

# Multivariable coefficients with OR and CI
multi_export <- data.frame(
  Variable   = rownames(s_multi$coefficients),
  Beta       = s_multi$coefficients[, "Estimate"],
  SE         = s_multi$coefficients[, "Std. Error"],
  z_value    = s_multi$coefficients[, "z value"],
  P_value    = s_multi$coefficients[, "Pr(>|z|)"],
  row.names  = NULL
)
multi_export$OR      <- exp(multi_export$Beta)
multi_export$CI_low  <- exp(multi_export$Beta - 1.96 * multi_export$SE)
multi_export$CI_high <- exp(multi_export$Beta + 1.96 * multi_export$SE)

write.csv(multi_export,
          file.path(out_dir, "multivariable_lifestyle_subgroup.csv"),
          row.names = FALSE)
cat("Saved: multivariable_lifestyle_subgroup.csv\n")

# ===========================================================================
# 6. Forest plot (univariate ORs)
# ===========================================================================
cat("\n========== Generating figures ==========\n")

plot_uni <- uni_results
plot_uni <- plot_uni[order(plot_uni$odds_ratio), ]
plot_uni$y_label <- factor(plot_uni$label, levels = plot_uni$label)

# Significance annotation
plot_uni$sig_anno <- ""
plot_uni$sig_anno[plot_uni$significant] <- "*"

# Add variable type indicator to label
plot_uni$label_with_unit <- paste0(plot_uni$label, "  [", plot_uni$or_unit, "]")
plot_uni$y_label_full <- factor(plot_uni$label_with_unit,
                                 levels = plot_uni$label_with_unit)

# Separate shape for continuous vs binary
plot_uni$point_shape <- ifelse(plot_uni$var_type == "continuous", 16, 17)

p_forest <- ggplot(plot_uni, aes(x = odds_ratio, y = y_label_full, color = category)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", alpha = 0.6) +
  geom_point(aes(size = n, shape = var_type), alpha = 0.9) +
  geom_errorbarh(aes(xmin = or_ci_low, xmax = or_ci_high), height = 0.15, alpha = 0.9) +
  scale_shape_manual(
    values = c("continuous" = 16, "binary" = 17),
    labels = c("continuous" = "Continuous (per 1-SD)", "binary" = "Binary (1 vs 0)"),
    name = "Variable Type"
  ) +
  scale_x_log10(
    breaks = c(0.1, 0.25, 0.5, 1, 2, 4, 8),
    labels = c("0.1", "0.25", "0.5", "1", "2", "4", "8")
  ) +
  scale_color_manual(
    values = c("Nutrition" = "#2c7bb6", "Neuropsychiatric" = "#d7191c",
               "Substance Use" = "#fdae61"),
    name = ""
  ) +
  scale_size_continuous(name = "N", range = c(2.5, 5)) +
  labs(
    title = "Lifestyle Factors Associated with Subgroup Membership",
    subtitle = sprintf(
      "%d subjects (%d cases) | * FDR < 0.05 | Scale: continuous vars per 1-SD; binary vars 1 vs 0",
      nrow(df_analysis), sum(df_analysis$Subgroup == 1)
    ),
    x = "Odds Ratio (log scale)",
    y = ""
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 14),
    plot.subtitle    = element_text(size = 8.5, color = "grey40"),
    legend.position  = "bottom",
    axis.text.y      = element_text(size = 9),
    panel.grid.major.y = element_blank()
  )

ggsave(file.path(out_dir, "lifestyle_subgroup_univariate_forest.png"),
       p_forest, width = 10, height = 5, dpi = 300, device = "png")
cat("Saved: lifestyle_subgroup_univariate_forest.png\n")

# ===========================================================================
# 7. Multivariable coefficient plot (OR)
# ===========================================================================

multi_plot_data <- multi_export[multi_export$Variable != "(Intercept)", ]
multi_plot_data$label <- ifelse(
  multi_plot_data$Variable %in% names(var_labels),
  var_labels[multi_plot_data$Variable],
  multi_plot_data$Variable
)
multi_plot_data$is_lifestyle <- multi_plot_data$Variable %in% LIFESTYLE_VARS
multi_plot_data$or_unit <- ifelse(
  multi_plot_data$Variable %in% names(or_unit),
  or_unit[multi_plot_data$Variable],
  ""
)
# Add unit annotation
multi_plot_data$label_with_unit <- ifelse(
  multi_plot_data$or_unit != "",
  paste0(multi_plot_data$label, "  [", multi_plot_data$or_unit, "]"),
  multi_plot_data$label
)

multi_plot_data <- multi_plot_data[order(multi_plot_data$OR), ]
multi_plot_data$y_label <- factor(multi_plot_data$label_with_unit,
                                   levels = multi_plot_data$label_with_unit)

# P-value annotation
multi_plot_data$p_label <- sprintf("P = %.3f", multi_plot_data$P_value)
multi_plot_data$p_label[multi_plot_data$P_value < 0.001] <- "P < 0.001"

# Point shape by variable type
multi_plot_data$pt_type <- ifelse(
  multi_plot_data$Variable %in% names(VAR_TYPE),
  VAR_TYPE[multi_plot_data$Variable],
  "continuous"
)

p_multi <- ggplot(multi_plot_data, aes(x = OR, y = y_label, color = is_lifestyle)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", alpha = 0.6) +
  geom_point(aes(shape = pt_type), size = 3, alpha = 0.9) +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.2, alpha = 0.9) +
  geom_text(aes(x = max(pmin(CI_high, 20)) * 0.55,
                label = p_label),
            size = 3.2, hjust = 1, color = "grey30") +
  scale_x_log10() +
  scale_shape_manual(
    values = c("continuous" = 16, "binary" = 17),
    labels = c("continuous" = "Continuous (per 1-SD)", "binary" = "Binary (1 vs 0)"),
    name = "Variable Type"
  ) +
  scale_color_manual(values = c("TRUE" = "#2c7bb6", "FALSE" = "grey50"),
                      labels = c("TRUE" = "Lifestyle", "FALSE" = "Covariate"),
                      name = "") +
  labs(
    title = "Multivariable Model: Lifestyle Factors → Subgroup",
    subtitle = sprintf(
      "Subgroup ~ DHA + EPA + HCys + NPIK + NPIKTOT + MH14ALCH + MH16SMOK\n%d complete cases (%d cases) | AIC = %.1f",
      nrow(df_multi_complete), sum(df_multi_complete$Subgroup == 1), fit_multi$aic
    ),
    x = "Odds Ratio (log scale)",
    y = ""
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.subtitle    = element_text(size = 8, color = "grey40"),
    legend.position  = "bottom",
    axis.text.y      = element_text(size = 9),
    panel.grid.major.y = element_blank()
  )

ggsave(file.path(out_dir, "lifestyle_subgroup_multivariable.png"),
       p_multi, width = 10, height = 5, dpi = 300, device = "png")
cat("Saved: lifestyle_subgroup_multivariable.png\n")

# ===========================================================================
# 8. Combined table (univariate + multivariable)
# ===========================================================================
cat("\n========== Combined Table ==========\n")

# Build comparison table
comp_table <- data.frame(
  Variable   = var_labels[LIFESTYLE_VARS],
  Category   = var_cat[LIFESTYLE_VARS],
  Type       = VAR_TYPE[LIFESTYLE_VARS],
  Unit       = or_unit[LIFESTYLE_VARS],
  stringsAsFactors = FALSE
)
rownames(comp_table) <- LIFESTYLE_VARS

# Add univariate results
comp_table$Uni_OR   <- NA_real_
comp_table$Uni_CI   <- NA_character_
comp_table$Uni_P    <- NA_real_
comp_table$Uni_FDR  <- NA_real_
comp_table$Uni_N    <- NA_integer_

for (v in LIFESTYLE_VARS) {
  r <- uni_results[uni_results$variable == v, ]
  if (nrow(r) == 1) {
    comp_table[v, "Uni_OR"]   <- r$odds_ratio
    comp_table[v, "Uni_CI"]   <- sprintf("[%.3f, %.3f]", r$or_ci_low, r$or_ci_high)
    comp_table[v, "Uni_P"]    <- r$p_value
    comp_table[v, "Uni_FDR"]  <- r$fdr
    comp_table[v, "Uni_N"]    <- r$n
  }
}

# Add multivariable results
comp_table$Multi_OR  <- NA_real_
comp_table$Multi_CI  <- NA_character_
comp_table$Multi_P   <- NA_real_
for (v in LIFESTYLE_VARS) {
  if (v %in% multi_export$Variable) {
    r <- multi_export[multi_export$Variable == v, ]
    comp_table[v, "Multi_OR"]  <- r$OR
    comp_table[v, "Multi_CI"]  <- sprintf("[%.3f, %.3f]", r$CI_low, r$CI_high)
    comp_table[v, "Multi_P"]   <- r$P_value
  }
}

# Round for display
comp_table_display <- comp_table
comp_table_display$Uni_OR    <- sprintf("%.3f", comp_table$Uni_OR)
comp_table_display$Uni_P     <- formatC(comp_table$Uni_P, format = "e", digits = 2)
comp_table_display$Uni_FDR   <- formatC(comp_table$Uni_FDR, format = "e", digits = 2)
comp_table_display$Multi_OR  <- sprintf("%.3f", comp_table$Multi_OR)
comp_table_display$Multi_P   <- formatC(comp_table$Multi_P, format = "e", digits = 2)

# Add FDR significance marker
comp_table_display$Uni_Sig <- ""
comp_table_display$Uni_Sig[comp_table$Uni_FDR < 0.05] <- "*"

cat("\nCombined Uni- vs Multi-variable Results:\n")
print(comp_table_display[, c("Variable", "Type", "Unit", "Uni_N", "Uni_OR", "Uni_CI",
                              "Uni_P", "Uni_FDR", "Uni_Sig", "Multi_OR", "Multi_CI",
                              "Multi_P")], row.names = FALSE)

write.csv(comp_table_display,
          file.path(out_dir, "lifestyle_subgroup_combined_table.csv"),
          row.names = FALSE)
cat("\nSaved: lifestyle_subgroup_combined_table.csv\n")

# ===========================================================================
# 9. Subgroup distribution by lifestyle variables
# ===========================================================================
cat("\n========== Distribution plots ==========\n")

# Subgroup distribution overview
p_dist <- ggplot(df_analysis, aes(x = factor(Subgroup), fill = factor(Subgroup))) +
  geom_bar(width = 0.6, alpha = 0.85) +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5, size = 5) +
  scale_fill_manual(values = c("0" = "#2c7bb6", "1" = "#d7191c"),
                     labels = c("0" = "Subgroup 0", "1" = "Subgroup 1"),
                     name = "") +
  labs(
    title = "Subgroup Distribution in Analysis Sample",
    subtitle = sprintf("Baseline subjects with subtype data: %d total (%d Subgroup=1)",
                       nrow(df_analysis), sum(df_analysis$Subgroup == 1)),
    x = "Subgroup",
    y = "Count"
  ) +
  theme_bw(base_size = 13) +
  theme(legend.position = "none")

ggsave(file.path(out_dir, "subgroup_distribution.png"),
       p_dist, width = 5, height = 4.5, dpi = 300, device = "png")
cat("Saved: subgroup_distribution.png\n")

# Boxplot: continuous variables by Subgroup
for (v in LIFESTYLE_VARS[VAR_TYPE[LIFESTYLE_VARS] == "continuous"]) {
  df_v <- df_analysis[!is.na(df_analysis[[v]]), ]
  df_v$Subgroup_f <- factor(df_v$Subgroup)

  p_box <- ggplot(df_v, aes(x = Subgroup_f, y = .data[[v]], fill = Subgroup_f)) +
    geom_boxplot(alpha = 0.7, width = 0.5, outlier.size = 1.5) +
    geom_jitter(width = 0.15, alpha = 0.4, size = 1.5) +
    scale_fill_manual(values = c("0" = "#2c7bb6", "1" = "#d7191c"), guide = "none") +
    labs(
      title = paste(var_labels[v], "by Subgroup"),
      subtitle = sprintf("n=%d (0=%d, 1=%d) | Continuous, OR per 1-SD",
                         nrow(df_v), sum(df_v$Subgroup == 0), sum(df_v$Subgroup == 1)),
      x = "Subgroup",
      y = var_labels[v]
    ) +
    theme_bw(base_size = 12)

  safe_v <- gsub("[^A-Za-z0-9_]", "_", v)
  ggsave(file.path(out_dir, sprintf("boxplot_%s_by_subgroup.png", safe_v)),
         p_box, width = 5.5, height = 4, dpi = 300, device = "png")
}
cat("Saved: boxplot_*_by_subgroup.png for continuous variables\n")

# Bar plot: binary variables by Subgroup
for (v in LIFESTYLE_VARS[VAR_TYPE[LIFESTYLE_VARS] == "binary"]) {
  df_v <- df_analysis[!is.na(df_analysis[[v]]), ]
  df_v$Subgroup_f <- factor(df_v$Subgroup)
  df_v[[v]] <- factor(df_v[[v]])

  # Compute percentages
  p_bar <- ggplot(df_v, aes(x = .data[[v]], fill = Subgroup_f)) +
    geom_bar(position = "fill", width = 0.5, alpha = 0.85) +
    scale_fill_manual(values = c("0" = "#2c7bb6", "1" = "#d7191c"),
                       labels = c("0" = "Subgroup 0", "1" = "Subgroup 1"),
                       name = "") +
    scale_y_continuous(labels = scales::percent) +
    labs(
      title = paste(var_labels[v], "by Subgroup"),
      subtitle = sprintf("n=%d (0=%d, 1=%d) | Binary, OR = 1 vs 0",
                         nrow(df_v), sum(df_v$Subgroup == 0), sum(df_v$Subgroup == 1)),
      x = var_labels[v],
      y = "Proportion"
    ) +
    theme_bw(base_size = 12)

  safe_v <- gsub("[^A-Za-z0-9_]", "_", v)
  ggsave(file.path(out_dir, sprintf("barplot_%s_by_subgroup.png", safe_v)),
         p_bar, width = 5, height = 4, dpi = 300, device = "png")
}
cat("Saved: barplot_*_by_subgroup.png for binary variables\n")

# ===========================================================================
# 10. Summary
# ===========================================================================
cat("\n")
cat("========== Analysis Summary ==========\n")
cat(sprintf("Total subjects with Subgroup data: %d (Subgroup=1: %d)\n",
            nrow(df_analysis), sum(df_analysis$Subgroup == 1)))
cat(sprintf("Variables tested: %d (%d continuous, %d binary)\n",
            nrow(uni_results),
            sum(uni_results$var_type == "continuous"),
            sum(uni_results$var_type == "binary")))
cat(sprintf("Univariate: %d / %d significant (FDR < 0.05)\n",
            sum(uni_results$significant), nrow(uni_results)))
cat(sprintf("Multivariable: %d complete cases (%d cases), AIC = %.1f\n",
            nrow(df_multi_complete), sum(df_multi_complete$Subgroup == 1),
            fit_multi$aic))

if (sum(uni_results$significant) > 0) {
  cat("\nSignificant lifestyle factors (FDR < 0.05):\n")
  for (i in which(uni_results$significant)) {
    cat(sprintf("  %s: OR = %.3f [%s], P = %.2e, FDR = %.2e\n",
                uni_results$label[i], uni_results$odds_ratio[i],
                uni_results$or_unit[i],
                uni_results$p_value[i], uni_results$fdr[i]))
  }
} else {
  cat("\nNo lifestyle variable reached FDR < 0.05 significance.\n")
  cat("Top associations (by raw p-value):\n")
  top5 <- head(uni_results, min(5, nrow(uni_results)))
  for (i in 1:nrow(top5)) {
    cat(sprintf("  %s: OR = %.3f [%.3f, %.3f] (%s), P = %.2e, FDR = %.2e\n",
                top5$label[i], top5$odds_ratio[i],
                top5$or_ci_low[i], top5$or_ci_high[i],
                top5$or_unit[i],
                top5$p_value[i], top5$fdr[i]))
  }
}

cat(sprintf("\nOutput directory: %s\n", normalizePath(out_dir)))
cat("\n========== Done ==========\n")
