###############################################################################
# 11_lifestyle_7var_ptau217.R
# Lifestyle → pTau217 analysis using lifestyle_7var_ptau.csv
# Same logic as 05_lifestyle_ptau_simplified.R but:
#  - Input: lifestyle_7var_ptau.csv (longitudinal, merged with master for SEX)
#  - Outcome: pTau217 (both plasma & CSF)
#  - Univariate: log2(pTau217) ~ lifestyle + AGE + SEX (FDR corrected)
#  - Multivariable: log2(pTau217) ~ all 7 + AGE + SEX
###############################################################################

library(readxl)
library(dplyr)
library(ggplot2)
library(ggtext)

set.seed(42)

out_dir <- "output/lifestyle_7var_ptau217"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

LIFESTYLE_VARS <- c("DHA", "EPA", "HCys", "NPIK", "NPIKTOT", "MH14ALCH", "MH16SMOK")

VAR_TYPE <- c("DHA"="continuous", "EPA"="continuous", "HCys"="continuous",
              "NPIK"="binary", "NPIKTOT"="continuous",
              "MH14ALCH"="binary", "MH16SMOK"="binary")

var_labels <- c(
  "DHA"="DHA (n-3 PUFA)", "EPA"="EPA (n-3 PUFA)", "HCys"="Homocysteine",
  "NPIK"="NPI Sleep Domain", "NPIKTOT"="NPI Total Score",
  "MH14ALCH"="Alcohol Abuse History", "MH16SMOK"="Smoking Status"
)
var_cat <- c(
  "DHA (n-3 PUFA)"="Nutrition", "EPA (n-3 PUFA)"="Nutrition", "Homocysteine"="Nutrition",
  "NPI Sleep Domain"="Neuropsychiatric", "NPI Total Score"="Neuropsychiatric",
  "Alcohol Abuse History"="Substance Use", "Smoking Status"="Substance Use"
)

cat("========== Lifestyle → pTau217 (Plasma & CSF) ==========\n")

# ---- 1. Load and prepare data ----
cat("\n--- Loading data ---\n")
csv <- read.csv("lifestyle_7var_ptau.csv", stringsAsFactors = FALSE)
master <- read_excel("master_data.xlsx", sheet = "Sheet1")

# Merge with master to get SEX
sex_lookup <- master[!duplicated(master$RID), c("RID", "PTGENDER")]
csv <- merge(csv, sex_lookup, by = "RID", all.x = TRUE)
csv$SEX <- ifelse(csv$PTGENDER == "Male", 1L, 0L)
cat(sprintf("CSV: %d rows, SEX merged: %d/%d\n", nrow(csv), sum(!is.na(csv$SEX)), nrow(csv)))

# ---- 2. Analysis function ----
run_ptau_analysis <- function(df, ptau_col, age_col, label, n_row) {
  cat(sprintf("\n========== %s ==========\n", label))

  # Filter to complete pTau + AGE (SEX optional)
  df_a <- df[!is.na(df[[ptau_col]]) & !is.na(df[[age_col]]), ]
  #CSF ALREADY LOGGED BUT NOT PLASMA
  if (ptau_col == "ptau217_csf") {
    df_a$PTAU_log2 <- df_a[[ptau_col]]
  } else {
    df_a$PTAU_log2 <- log2(df_a[[ptau_col]])
  }
  df_a$AGE_val <- df_a[[age_col]]
  # Only use SEX if >50% of subjects have it
  sex_coverage <- sum(!is.na(df_a$SEX)) / nrow(df_a)
  has_sex <- sex_coverage > 0.5 && length(unique(na.omit(df_a$SEX))) >= 2
  cat(sprintf("Complete cases (pTau+AGE): %d | with SEX: %d | SEX in model: %s\n",
              nrow(df_a), sum(!is.na(df_a$SEX)), ifelse(has_sex, "yes", "no")))

  # ---- Univariate ----
  uni_results <- data.frame(
    variable = LIFESTYLE_VARS, n = NA_integer_, n_cases = NA_integer_,
    beta = NA_real_, se = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
    t_value = NA_real_, p_value = NA_real_, model_r2 = NA_real_
  )

  for (i in seq_along(LIFESTYLE_VARS)) {
    v <- LIFESTYLE_VARS[i]
    df_sub <- df_a[!is.na(df_a[[v]]), ]
    uni_results$n[i] <- nrow(df_sub)

    x_val <- df_sub[[v]]
    if (VAR_TYPE[v] == "continuous") {
      x_sd <- sd(x_val, na.rm = TRUE)
      x_model <- if (x_sd > 0) as.vector(scale(x_val)) else x_val
    } else {
      x_model <- x_val
    }

    if (has_sex) {
      fit <- lm(PTAU_log2 ~ x_model + AGE_val + SEX, data = df_sub)
    } else {
      fit <- lm(PTAU_log2 ~ x_model + AGE_val, data = df_sub)
    }
    s <- summary(fit)
    uni_results$beta[i]   <- s$coefficients["x_model", "Estimate"]
    uni_results$se[i]     <- s$coefficients["x_model", "Std. Error"]
    uni_results$t_value[i]<- s$coefficients["x_model", "t value"]
    uni_results$p_value[i]<- s$coefficients["x_model", "Pr(>|t|)"]
    uni_results$model_r2[i]<- s$r.squared
    uni_results$ci_low[i] <- uni_results$beta[i] - 1.96 * uni_results$se[i]
    uni_results$ci_high[i]<- uni_results$beta[i] + 1.96 * uni_results$se[i]
  }

  uni_results$fdr <- p.adjust(uni_results$p_value, method = "BH")
  uni_results$significant <- uni_results$fdr < 0.05
  uni_results <- uni_results[order(uni_results$p_value), ]
  uni_results$label <- var_labels[uni_results$variable]
  uni_results$category <- var_cat[uni_results$variable]

  cat(sprintf("FDR < 0.05: %d / %d\n", sum(uni_results$significant), nrow(uni_results)))

  # ---- Multivariable ----
  var_list <- paste(LIFESTYLE_VARS, collapse = " + ")
  if (has_sex) {
    formula_multi <- as.formula(paste("PTAU_log2 ~", var_list, "+ AGE_val + SEX"))
    df_multi <- df_a[, c("PTAU_log2", LIFESTYLE_VARS, "AGE_val", "SEX")]
  } else {
    formula_multi <- as.formula(paste("PTAU_log2 ~", var_list, "+ AGE_val"))
    df_multi <- df_a[, c("PTAU_log2", LIFESTYLE_VARS, "AGE_val")]
  }
  df_multi_complete <- na.omit(df_multi)
  for (v in LIFESTYLE_VARS) {
    if (VAR_TYPE[v] == "continuous") {
      df_multi_complete[[v]] <- as.vector(scale(df_multi_complete[[v]]))
    }
  }

  fit_multi <- lm(formula_multi, data = df_multi_complete)
  s_multi <- summary(fit_multi)
  cat(sprintf("Complete cases: %d | R²=%.4f Adj R²=%.4f | P=%.2e\n",
              nrow(df_multi_complete), s_multi$r.squared,
              s_multi$adj.r.squared,
              pf(s_multi$fstatistic[1], s_multi$fstatistic[2],
                 s_multi$fstatistic[3], lower.tail = FALSE)))
  cat(sprintf("AGE P=%.4f", s_multi$coefficients["AGE_val", "Pr(>|t|)"]))
  if (has_sex) cat(sprintf(" | SEX P=%.4f", s_multi$coefficients["SEX", "Pr(>|t|)"]))
  cat("\n")

  # ---- Forest plot ----
  plot_uni <- uni_results
  plot_uni <- plot_uni[order(plot_uni$beta), ]
  plot_uni$y_label <- factor(plot_uni$label, levels = plot_uni$label)
  plot_uni$sig_anno <- ifelse(plot_uni$significant, "*", "")

  p_forest <- ggplot(plot_uni, aes(x = beta, y = y_label, color = category)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", alpha = 0.6) +
    geom_point(aes(size = n), alpha = 0.9) +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.15, alpha = 0.9) +
    scale_color_manual(values = c("Nutrition"="#2c7bb6","Neuropsychiatric"="#d7191c",
                                   "Substance Use"="#fdae61"), name = "") +
    scale_size_continuous(name = "N", range = c(2.5, 5)) +
    labs(
      title = paste("Lifestyle Factors Associated with", label),
      subtitle = sprintf("Univariate: log2(pTau) ~ variable + AGE%s | %d subjects | * FDR<0.05",
                         ifelse(has_sex, " + SEX", ""), nrow(df_a)),
      x = "Standardized Effect (per 1-SD or 1 vs 0)", y = ""
    ) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(size = 8.5, color = "grey40"),
          legend.position = "bottom", panel.grid.major.y = element_blank())

  ggsave(file.path(out_dir, sprintf("forest_%s.png", gsub(" ", "_", tolower(label)))),
         p_forest, width = 9, height = 4.5, dpi = 300)
  cat(sprintf("Saved: forest_%s\n", tolower(label)))

  # ---- Multivariable coefficient plot ----
  multi_extract <- as.data.frame(s_multi$coefficients)
  multi_extract$Variable <- rownames(multi_extract)
  colnames(multi_extract)[1:4] <- c("Beta", "SE", "t", "P")

  # Keep lifestyle + AGE + SEX
  multi_plot <- multi_extract[!grepl("Intercept", multi_extract$Variable), ]
  multi_plot$label <- ifelse(multi_plot$Variable %in% names(var_labels),
                              var_labels[multi_plot$Variable], multi_plot$Variable)
  multi_plot$is_lifestyle <- multi_plot$Variable %in% LIFESTYLE_VARS
  multi_plot$ci_low  <- multi_plot$Beta - 1.96 * multi_plot$SE
  multi_plot$ci_high <- multi_plot$Beta + 1.96 * multi_plot$SE
  multi_plot <- multi_plot[order(multi_plot$Beta), ]
  multi_plot$y_label <- factor(multi_plot$label, levels = multi_plot$label)
  multi_plot$p_label <- sprintf("P=%.3f", multi_plot$P)
  multi_plot$p_label[multi_plot$P < 0.001] <- "P<0.001"
  multi_plot$lsgroup <- case_when(
    multi_plot$label %in% c("Homocysteine", "EPA (n-3 PUFA)", "DHA (n-3 PUFA)") ~ "Nutrition",
    multi_plot$label %in% c("NPI Total Score", "NPI Sleep Domain") ~ "Neuropsychiatric",
    multi_plot$label %in% c("Alcohol Abuse History", "Smoking Status") ~ "Substance Use",
    TRUE ~ "Covariate"
    
  )
  multi_plot$label_color <- case_when(
    multi_plot$lsgroup == "Nutrition" ~ "#D4A017",
    multi_plot$lsgroup == "Neuropsychiatric" ~ "#C44E52",
    multi_plot$lsgroup == "Substance Use" ~ "#D2691E",
    TRUE ~ "black"
  )
  multi_plot$label_html <- paste0(
    "<span style='color:", 
    multi_plot$label_color,
    "'>",
    multi_plot$label,
    "</span>"
  )


  p_multi <- ggplot(multi_plot, aes(x = Beta, y = y_label, color = lsgroup)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", alpha = 0.6) +
    geom_point(size = 3, alpha = 0.9) +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.2, alpha = 0.9) +
    geom_text(
      aes(x = Beta, y = y_label, label = p_label),
      nudge_y = 0.30,
      size = 3.2,
      color = "grey30"
    ) +
    scale_color_manual(
      values = c(
        "Nutrition" = "#D4A017",
        "Neuropsychiatric" = "#C44E52",
        "Substance Use" = "#D2691E"
      ),
      na.value = "grey30",
      name = "Category"
    ) +
    scale_y_discrete(
      labels = multi_plot$label_html
    ) +
    labs(
      title = paste("Lifestyle Factors Are Not Directly Associated with", label),
      subtitle = sprintf("Multivariate: log(pTau) ~ 7 lifestyle + AGE%s | %d complete cases",
                         ifelse(has_sex, " + SEX", ""),
                         nrow(df_multi_complete), s_multi$r.squared, s_multi$adj.r.squared),
      x = "Standardized Effect", y = ""
    ) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(size = 8, color = "grey40"),
          legend.position = "none") +
    theme(
      axis.text.y = ggtext::element_markdown(size = 11, face = "bold")
    )

    print(p_multi)
  ggsave(file.path(out_dir, sprintf("multivariable_%s.png", gsub(" ", "_", tolower(label)))),
         p_multi, width = 9, height = 5, dpi = 300)
  cat(sprintf("Saved: multivariable_%s\n", tolower(label)))

  # ---- Combined table ----
  multi_extract <- as.data.frame(s_multi$coefficients)
  multi_extract$Variable <- rownames(multi_extract)
  colnames(multi_extract)[1:4] <- c("Beta", "SE", "t", "P")

  comp_table <- data.frame(Variable = var_labels[LIFESTYLE_VARS], Category = var_cat[LIFESTYLE_VARS])
  for (v in LIFESTYLE_VARS) {
    r <- uni_results[uni_results$variable == v, ]
    comp_table[v, "Uni_Effect"] <- sprintf("%.4f", r$beta)
    comp_table[v, "Uni_P"] <- formatC(r$p_value, format="e", digits=2)
    comp_table[v, "Uni_FDR"] <- formatC(r$fdr, format="e", digits=2)
    comp_table[v, "Uni_N"] <- r$n
    mr <- multi_extract[multi_extract$Variable == v, ]
    if (nrow(mr) == 1) {
      comp_table[v, "Multi_Effect"] <- sprintf("%.4f", mr$Beta)
      comp_table[v, "Multi_P"] <- formatC(mr$P, format="e", digits=2)
    } else {
      comp_table[v, "Multi_Effect"] <- "-"
      comp_table[v, "Multi_P"] <- "-"
    }
  }
  comp_table$Uni_Sig <- ifelse(uni_results$significant[match(LIFESTYLE_VARS, uni_results$variable)], "*", "")
  print(comp_table)
  write.csv(comp_table, file.path(out_dir, sprintf("combined_table_%s.csv", gsub(" ", "_", tolower(label)))), row.names = FALSE)

  return(list(uni = uni_results, multi = s_multi, n = nrow(df_a)))
}

# ---- 3. Run: CSF pTau217 ----
res_csf <- run_ptau_analysis(csv, "ptau217_csf", "age_at_csf_ptau",
                              "CSF pTau217", nrow(csv))

# ---- 4. Run: Plasma pTau217 ----
res_plasma <- run_ptau_analysis(csv, "ptau217_plasma", "age_at_plasma_ptau",
                                 "Plasma pTau217", nrow(csv))

# ---- 5. Comparison forest plot (CSF vs Plasma side-by-side) ----
cat("\n========== CSF vs Plasma Comparison ==========\n")

compare_uni <- merge(
  res_csf$uni[, c("variable", "beta", "p_value", "fdr")],
  res_plasma$uni[, c("variable", "beta", "p_value", "fdr")],
  by = "variable", suffixes = c("_CSF", "_Plasma")
)
compare_uni$label <- var_labels[compare_uni$variable]
compare_uni$category <- var_cat[compare_uni$variable]
compare_uni$fdr_min <- pmin(compare_uni$fdr_CSF, compare_uni$fdr_Plasma)
compare_uni$significant <- compare_uni$fdr_min < 0.05

# Long format for side-by-side
comp_long <- rbind(
  data.frame(compare_uni, Assay = "CSF pTau217", beta = compare_uni$beta_CSF,
             p_value = compare_uni$p_value_CSF),
  data.frame(compare_uni, Assay = "Plasma pTau217", beta = compare_uni$beta_Plasma,
             p_value = compare_uni$p_value_Plasma)
)
comp_long$label <- compare_uni$label[match(comp_long$variable, compare_uni$variable)]
comp_long$category <- compare_uni$category[match(comp_long$variable, compare_uni$variable)]

# Order by CSF beta
comp_long$label <- factor(comp_long$label, levels = compare_uni$label[order(compare_uni$beta_CSF)])

p_comp <- ggplot(comp_long, aes(x = beta, y = label, color = Assay, shape = Assay)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 3, alpha = 0.85, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c("CSF pTau217" = "#2c7bb6", "Plasma pTau217" = "#d7191c")) +
  labs(
    title = "Lifestyle → pTau217: CSF vs Plasma Comparison",
    subtitle = sprintf("CSF: %d subjects | Plasma: %d subjects | Univariate: log2(pTau) ~ variable + AGE",
                       res_csf$n, res_plasma$n),
    x = "Effect Size", y = ""
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 8, color = "grey40"),
        legend.position = "bottom")

ggsave(file.path(out_dir, "comparison_csf_vs_plasma.png"),
       p_comp, width = 9, height = 4.5, dpi = 300)
cat("Saved: comparison_csf_vs_plasma.png\n")

# ---- 6. Summary ----
cat("\n")
cat("========== Summary ==========\n")
cat(sprintf("Input: %d rows (longitudinal)\n", nrow(csv)))
cat(sprintf("CSF pTau217: %d subjects with complete data\n", res_csf$n))
cat(sprintf("Plasma pTau217: %d subjects with complete data\n", res_plasma$n))
cat(sprintf("CSF FDR<0.05: %d/%d\n", sum(res_csf$uni$significant), nrow(res_csf$uni)))
cat(sprintf("Plasma FDR<0.05: %d/%d\n", sum(res_plasma$uni$significant), nrow(res_plasma$uni)))
cat(sprintf("\nOutput: %s\n", normalizePath(out_dir)))
cat("========== Done ==========\n")
