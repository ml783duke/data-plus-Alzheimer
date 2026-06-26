###############################################################################
# 04_ptau_protein_age_adjusted.R
# Age-adjusted PTAU-protein association
# - Remove AGE effect from both log2(PTAU) and log2(protein) via linear residuals
# - Spearman correlation on residuals → "age-independent PTAU-protein association"
# - FDR correction, volcano plot, comparison with original Spearman results
# - Overlap analysis with 60 curated PTAU-associated proteins
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

out_dir <- "output/ptau_protein_age_adjusted"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ===========================================================================
# 1. Load data
# ===========================================================================
cat("\n========== Loading data ==========\n")
df_raw <- read_excel("master_data.xlsx", sheet = "Sheet1")
cat(sprintf("Raw data: %d rows x %d columns\n", nrow(df_raw), ncol(df_raw)))

dict <- read.csv("protein_raw_data/protein dict.csv", stringsAsFactors = FALSE)
cat(sprintf("Protein dictionary: %d entries\n", nrow(dict)))

# Load curated 60 PTAU-associated proteins
curated_pos <- read_excel("output/ptau_protein_correlation/PTAU_positive_negative_proteins.xlsx",
                           sheet = "Positive Proteins", skip = 1)
curated_neg <- read_excel("output/ptau_protein_correlation/PTAU_positive_negative_proteins.xlsx",
                           sheet = "Negative Proteins", skip = 1)
curated_all <- rbind(curated_pos, curated_neg)
colnames(curated_all)[2] <- "protein_id"
curated_ids <- curated_all[["protein_id"]]
cat(sprintf("Curated PTAU-associated proteins: %d\n", length(curated_ids)))

# ===========================================================================
# 2. Filter to baseline, non-missing PTAU & AGE
# ===========================================================================
cat("\n========== Filtering samples ==========\n")
df_bl <- df_raw[df_raw$VISCODE2 == "bl", ]
cat(sprintf("Baseline samples: %d\n", nrow(df_bl)))

df_analysis <- df_bl[!is.na(df_bl$PTAU) & !is.na(df_bl$AGE), ]
cat(sprintf("Baseline with non-missing PTAU & AGE: %d\n", nrow(df_analysis)))

cat(sprintf("AGE range: %.1f - %.1f (median %.1f)\n",
            min(df_analysis$AGE), max(df_analysis$AGE), median(df_analysis$AGE)))

# ===========================================================================
# 3. Protein columns & numeric conversion
# ===========================================================================
cat("\n========== Processing protein columns ==========\n")
protein_cols <- grep("^X[0-9]+\\.[0-9]+$", colnames(df_analysis), value = TRUE)
cat(sprintf("SomaScan protein columns: %d\n", length(protein_cols)))

convert_protein <- function(x) {
  x[x == "NA" | x == ""] <- NA
  as.numeric(x)
}
df_analysis[protein_cols] <- lapply(df_analysis[protein_cols], convert_protein)

# ===========================================================================
# 4. Filter proteins by missingness
# ===========================================================================
cat("\n========== Filtering proteins by missingness (>20%) ==========\n")
missing_rate <- sapply(df_analysis[protein_cols], function(x) mean(is.na(x)) * 100)
proteins_kept <- names(missing_rate[missing_rate <= 20])
cat(sprintf("Proteins kept: %d / %d\n", length(proteins_kept), length(protein_cols)))

# ===========================================================================
# 5. Log2 transform
# ===========================================================================
cat("\n========== Log2 transformation ==========\n")
df_analysis$PTAU_log2 <- log2(df_analysis$PTAU)
protein_mat  <- as.matrix(df_analysis[proteins_kept])
protein_log2 <- log2(protein_mat)
protein_log2[is.infinite(protein_log2) | is.nan(protein_log2)] <- NA
cat(sprintf("Log2 protein matrix: %d x %d\n", nrow(protein_log2), ncol(protein_log2)))

# ===========================================================================
# 6. Age-adjusted residual correlation (core analysis)
# ===========================================================================
cat("\n========== Age-Adjusted Residual Correlation ==========\n")

age_vals  <- df_analysis$AGE
ptau_log2 <- df_analysis$PTAU_log2

# ---- Step 1: Remove age effect from PTAU ----
# Use lm on complete cases
ptau_lm <- lm(ptau_log2 ~ age_vals)
ptau_residuals <- residuals(ptau_lm)

cat(sprintf("PTAU ~ AGE: R² = %.4f, p = %.2e\n",
            summary(ptau_lm)$r.squared,
            summary(ptau_lm)$coefficients[2, 4]))

# ---- Step 2: For each protein, remove age effect, then correlate residuals ----
n_proteins <- ncol(protein_log2)
results <- data.frame(
  protein_id          = colnames(protein_log2),
  n_valid             = NA_integer_,
  raw_spearman_rho    = NA_real_,    # PTAU vs protein (no age adjustment)
  residual_spearman_rho = NA_real_,  # PTAU residuals vs protein residuals
  p_value             = NA_real_,
  age_effect_r2       = NA_real_,    # R² of protein ~ AGE (how much age explains)
  stringsAsFactors = FALSE
)

cat(sprintf("Testing %d proteins...\n", n_proteins))
pb <- txtProgressBar(min = 0, max = n_proteins, style = 3)

for (i in seq_len(n_proteins)) {
  setTxtProgressBar(pb, i)

  prot_log2 <- protein_log2[, i]

  # Valid pairs (non-missing AGE, PTAU, protein)
  valid_idx <- !is.na(prot_log2)
  prot_valid <- prot_log2[valid_idx]
  age_valid  <- age_vals[valid_idx]
  ptau_valid <- ptau_log2[valid_idx]
  n_valid    <- length(prot_valid)
  results$n_valid[i] <- n_valid

  if (n_valid >= 30) {
    # Raw Spearman (no age adjustment)
    raw_cor <- tryCatch(
      cor.test(prot_valid, ptau_valid, method = "spearman", exact = FALSE),
      error = function(e) NULL
    )
    if (!is.null(raw_cor)) {
      results$raw_spearman_rho[i] <- raw_cor$estimate
    }

    # Remove age effect from protein
    prot_lm <- tryCatch(
      lm(prot_valid ~ age_valid),
      error = function(e) NULL
    )
    if (!is.null(prot_lm)) {
      results$age_effect_r2[i] <- summary(prot_lm)$r.squared
      prot_residuals <- residuals(prot_lm)

      # Correlate residuals (PTAU residuals vs protein residuals)
      res_cor <- tryCatch(
        cor.test(ptau_residuals[valid_idx], prot_residuals,
                 method = "spearman", exact = FALSE),
        error = function(e) NULL
      )
      if (!is.null(res_cor)) {
        results$residual_spearman_rho[i] <- res_cor$estimate
        results$p_value[i]              <- res_cor$p.value
      }
    }
  }
}
close(pb)
cat("\n")

# ===========================================================================
# 7. FDR correction & post-processing
# ===========================================================================
cat("\n========== FDR Correction ==========\n")

results <- results[!is.na(results$p_value), ]
results$fdr <- p.adjust(results$p_value, method = "BH")
results$significant_fdr05 <- results$fdr < 0.05
results$significant_fdr01 <- results$fdr < 0.01

# Significance category for plot coloring
results$sig_cat <- "Not significant"
results$sig_cat[results$significant_fdr01] <- "FDR < 0.01"
results$sig_cat[results$significant_fdr05 & !results$significant_fdr01] <- "FDR < 0.05"
results$sig_cat <- factor(results$sig_cat,
                           levels = c("FDR < 0.01", "FDR < 0.05", "Not significant"))

# Direction
results$direction <- ifelse(results$residual_spearman_rho > 0, "Positive", "Negative")

# Compute "age-dependence shift": how much did the correlation change after age adjustment?
results$rho_shift <- results$residual_spearman_rho - results$raw_spearman_rho

n_sig_05 <- sum(results$significant_fdr05)
n_sig_01 <- sum(results$significant_fdr01)
cat(sprintf("Proteins tested: %d\n", nrow(results)))
cat(sprintf("FDR < 0.05: %d\n", n_sig_05))
cat(sprintf("FDR < 0.01: %d\n", n_sig_01))
if (n_sig_05 > 0) {
  cat(sprintf("  Positive: %d\n", sum(results$significant_fdr05 & results$direction == "Positive")))
  cat(sprintf("  Negative: %d\n", sum(results$significant_fdr05 & results$direction == "Negative")))
}

cat(sprintf("\nRaw Spearman rho range:     %.4f to %.4f\n",
            min(results$raw_spearman_rho, na.rm = TRUE),
            max(results$raw_spearman_rho, na.rm = TRUE)))
cat(sprintf("Residual Spearman rho range: %.4f to %.4f\n",
            min(results$residual_spearman_rho, na.rm = TRUE),
            max(results$residual_spearman_rho, na.rm = TRUE)))
cat(sprintf("Rho shift range:             %.4f to %.4f\n",
            min(results$rho_shift, na.rm = TRUE),
            max(results$rho_shift, na.rm = TRUE)))
cat(sprintf("Median age R² of proteins:   %.4f\n",
            median(results$age_effect_r2, na.rm = TRUE)))

# ===========================================================================
# 8. Annotate with gene symbols
# ===========================================================================
cat("\n========== Annotating with gene symbols ==========\n")
annot_lookup <- dict[, c("Analytes", "EntrezGeneSymbol", "TargetFullName", "Target", "UniProt")]
colnames(annot_lookup)[1] <- "protein_id"

results <- merge(results, annot_lookup, by = "protein_id", all.x = TRUE, sort = FALSE)
results <- results[order(results$fdr), ]

n_annotated <- sum(!is.na(results$EntrezGeneSymbol) & results$EntrezGeneSymbol != "")
cat(sprintf("Proteins with gene symbol: %d / %d\n", n_annotated, nrow(results)))

# ===========================================================================
# 9. Print top results
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

cat("\n========== Top 20 Age-Adjusted PTAU-Associated Proteins ==========\n")
top20 <- head(results, 20)
top20$gene_label <- get_gene_label(top20)
print(data.frame(
  Rank        = 1:20,
  Gene        = top20$EntrezGeneSymbol,
  Target      = top20$Target,
  Resid_Rho   = round(top20$residual_spearman_rho, 4),
  Raw_Rho     = round(top20$raw_spearman_rho, 4),
  Rho_Shift   = round(top20$rho_shift, 4),
  FDR         = formatC(top20$fdr, format = "e", digits = 2),
  N           = top20$n_valid
))

# ===========================================================================
# 10. Overlap with 60 curated PTAU-associated proteins
# ===========================================================================
cat("\n========== Overlap with Curated 60 PTAU Proteins ==========\n")

# Among the curated proteins, what are their ranks in the age-adjusted analysis?
curated_in_results <- results[results$protein_id %in% curated_ids, ]
curated_in_results <- curated_in_results[order(curated_in_results$fdr), ]
cat(sprintf("Curated proteins found in age-adjusted results: %d / %d\n",
            nrow(curated_in_results), length(curated_ids)))

# Show curated proteins sorted by age-adjusted FDR
cat("\nCurated 60 proteins ranked by age-adjusted FDR:\n")
for (i in seq_len(min(30, nrow(curated_in_results)))) {
  r <- curated_in_results[i, ]
  cat(sprintf("  FDR #%d: %s (%s) | Resid_ρ=%.4f | Raw_ρ=%.4f | FDR=%.2e | Excel: %s\n",
              which(results$protein_id == r$protein_id),
              r$EntrezGeneSymbol, r$Target,
              r$residual_spearman_rho, r$raw_spearman_rho, r$fdr,
              curated_all[["Direction"]][curated_all[["protein_id"]] == r$protein_id][1]))
}

# Overlap: curated proteins that are FDR < 0.05 in age-adjusted analysis
curated_sig <- curated_in_results[curated_in_results$significant_fdr05, ]
cat(sprintf("\nCurated proteins with FDR < 0.05 (age-adjusted): %d\n", nrow(curated_sig)))

# Top 30 age-adjusted vs curated
top30_adj <- head(results, 30)
overlap_adj_curated <- intersect(top30_adj$protein_id, curated_ids)
cat(sprintf("Top 30 age-adjusted ∩ Curated 60: %d proteins\n", length(overlap_adj_curated)))
if (length(overlap_adj_curated) > 0) {
  for (pid in overlap_adj_curated) {
    r <- results[results$protein_id == pid, ]
    cat(sprintf("  %s (%s) | Resid_ρ=%.4f | FDR=%.2e\n",
                r$EntrezGeneSymbol, r$Target, r$residual_spearman_rho, r$fdr))
  }
}

# Save curated overlap table
write.csv(curated_in_results,
          file.path(out_dir, "curated_60_proteins_age_adjusted.csv"),
          row.names = FALSE)

# ===========================================================================
# 11. Save CSVs
# ===========================================================================
cat("\n========== Saving CSVs ==========\n")

write.csv(results, file.path(out_dir, "all_proteins_age_adjusted.csv"), row.names = FALSE)
cat(sprintf("Saved: all_proteins_age_adjusted.csv (%d rows)\n", nrow(results)))

sig05 <- results[results$significant_fdr05, ]
write.csv(sig05, file.path(out_dir, "significant_proteins_fdr05.csv"), row.names = FALSE)
cat(sprintf("Saved: significant_proteins_fdr05.csv (%d rows)\n", nrow(sig05)))

# ===========================================================================
# 12. Visualization
# ===========================================================================
cat("\n========== Generating plots ==========\n")

# ---- 12a. Volcano Plot (age-adjusted) ----
plot_data <- results
plot_data$neg_log10_fdr <- -log10(plot_data$fdr)
fdr_cap <- 1e-16
plot_data$neg_log10_fdr_capped <- pmin(plot_data$neg_log10_fdr, -log10(fdr_cap))


# Label top 20 significant with gene symbols
plot_data$label <- ""
top_label <- head(plot_data[plot_data$significant_fdr05, ], 20)
top_label$gene_label <- get_gene_label(top_label)
plot_data$label[match(top_label$protein_id, plot_data$protein_id)] <- top_label$gene_label

# Mark curated proteins
plot_data$is_curated <- plot_data$protein_id %in% curated_ids

volcano_p <- ggplot(plot_data, aes(x = residual_spearman_rho, y = neg_log10_fdr_capped)) +
  geom_point(data = subset(plot_data, !is_curated),
             aes(color = sig_cat), alpha = 0.5, size = 0.7) +
  geom_point(data = subset(plot_data, is_curated),
             aes(fill = sig_cat), color = "black", stroke = 0.3,
             shape = 21, size = 1.5, alpha = 0.9) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", alpha = 0.6) +
  geom_hline(yintercept = -log10(0.01), linetype = "dotted", color = "grey40", alpha = 0.4) +
  geom_text_repel(
    data = subset(plot_data, label != ""),
    aes(label = label),
    size = 2.5, max.overlaps = 25, box.padding = 0.4,
    point.padding = 0.2, segment.size = 0.2, color = "black", fontface = "italic"
  ) +
  scale_color_manual(values = c("FDR < 0.01" = "#d7191c", "FDR < 0.05" = "#fdae61",
                                 "Not significant" = "grey70")) +
  scale_fill_manual(values = c("FDR < 0.01" = "#d7191c", "FDR < 0.05" = "#fdae61",
                                "Not significant" = "grey70"), guide = "none") +
  labs(
    title = "Age-Adjusted PTAU-Protein Association",
    subtitle = sprintf("Spearman ρ of residuals (age removed) | %d proteins tested | %d with FDR < 0.05 | ■ = curated proteins",
                       nrow(results), n_sig_05),
    x = "Residual Spearman ρ (age-adjusted)",
    y = expression(-log[10](FDR)),
    color = "Significance"
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), plot.subtitle = element_text(size = 8))

ggsave(file.path(out_dir, "volcano_age_adjusted.pdf"), volcano_p, width = 10, height = 7.5, device = "pdf")
ggsave(file.path(out_dir, "volcano_age_adjusted.png"), volcano_p, width = 10, height = 7.5, dpi = 300, device = "png")
cat("Saved: volcano_age_adjusted.pdf/png\n")

# ---- 12b. Rho Shift Plot: Raw rho vs Adjusted rho ----
rho_shift_p <- ggplot(results, aes(x = raw_spearman_rho, y = residual_spearman_rho)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50", alpha = 0.5) +
  geom_point(aes(color = sig_cat), alpha = 0.4, size = 0.6) +
  # Highlight curated proteins
  geom_point(data = subset(results, protein_id %in% curated_ids),
             color = "black", fill = "#ffff33", shape = 21, stroke = 0.5,
             size = 1.8, alpha = 0.85) +
  scale_color_manual(values = c("FDR < 0.01" = "#d7191c", "FDR < 0.05" = "#fdae61",
                                 "Not significant" = "grey70")) +
  labs(
    title = "Effect of Age Adjustment on PTAU-Protein Correlation",
    subtitle = sprintf("Diagonal = no age effect | ▲ = curated proteins (%d) | %d tested",
                       length(curated_ids), nrow(results)),
    x = "Raw Spearman ρ (log2(PTAU) vs log2(Protein))",
    y = "Age-Adjusted Spearman ρ (residuals after removing AGE)",
    color = "Significance\n(age-adjusted)"
  ) +
  coord_fixed(xlim = c(-1, 1), ylim = c(-1, 1)) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), plot.subtitle = element_text(size = 7.5))

ggsave(file.path(out_dir, "rho_shift_raw_vs_adjusted.pdf"), rho_shift_p, width = 8.5, height = 7.5, device = "pdf")
ggsave(file.path(out_dir, "rho_shift_raw_vs_adjusted.png"), rho_shift_p, width = 8.5, height = 7.5, dpi = 300, device = "png")
cat("Saved: rho_shift_raw_vs_adjusted.pdf/png\n")

# ---- 12c. Top 30 bar plot (age-adjusted, positive & negative separately) ----
top30_positive <- head(results[results$direction == "Positive", ], 30)
top30_negative <- head(results[results$direction == "Negative", ], 30)
top30_positive <- top30_positive[order(top30_positive$residual_spearman_rho), ]
top30_negative <- top30_negative[order(top30_negative$residual_spearman_rho, decreasing = TRUE), ]

top30_pos_neg <- rbind(top30_positive, top30_negative)
top30_pos_neg$gene_label <- get_gene_label(top30_pos_neg)
top30_pos_neg$gene_label <- factor(top30_pos_neg$gene_label,
                                    levels = top30_pos_neg$gene_label)

# Color by curated status
top30_pos_neg$is_curated <- top30_pos_neg$protein_id %in% curated_ids
top30_pos_neg$bar_color <- ifelse(top30_pos_neg$is_curated, "#d7191c", "#4575b4")

bar_p <- ggplot(top30_pos_neg, aes(x = gene_label, y = residual_spearman_rho)) +
  geom_bar(stat = "identity", aes(fill = bar_color), width = 0.7) +
  geom_text(aes(
    label = paste0(round(residual_spearman_rho, 3)),
    hjust = ifelse(residual_spearman_rho > 0, -0.15, 1.15)
  ), size = 2.2, color = "grey30") +
  scale_fill_identity() +
  coord_flip() +
  labs(
    title = "Top 30 Age-Adjusted PTAU-Associated Proteins",
    subtitle = sprintf("Red = curated protein | FDR < 0.05 | %d baseline subjects | AGE effect removed",
                       nrow(df_analysis)),
    x = "", y = "Age-Adjusted Spearman ρ"
  ) +
  theme_bw(base_size = 10) +
  theme(plot.title = element_text(face = "bold"), plot.subtitle = element_text(size = 7.5))

ggsave(file.path(out_dir, "top30_barplot_age_adjusted.pdf"), bar_p, width = 10, height = 8.5, device = "pdf")
ggsave(file.path(out_dir, "top30_barplot_age_adjusted.png"), bar_p, width = 10, height = 8.5, dpi = 300, device = "png")
cat("Saved: top30_barplot_age_adjusted.pdf/png\n")

# ---- 12d. Curated proteins: raw vs adjusted correlation shift ----
curated_plot_data <- results[results$protein_id %in% curated_ids, ]
curated_plot_data$gene_label <- get_gene_label(curated_plot_data)
curated_plot_data <- curated_plot_data[order(curated_plot_data$residual_spearman_rho), ]
curated_plot_data$gene_label <- factor(curated_plot_data$gene_label,
                                        levels = curated_plot_data$gene_label)

curated_shift_p <- ggplot(curated_plot_data, aes(y = gene_label)) +
  geom_segment(aes(x = raw_spearman_rho, xend = residual_spearman_rho,
                   yend = gene_label),
               arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
               color = "grey50", linewidth = 0.6) +
  geom_point(aes(x = raw_spearman_rho), color = "#fc8d59", size = 2.5) +
  geom_point(aes(x = residual_spearman_rho), color = "#2c7bb6", size = 2.5) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey60") +
  labs(
    title = "Effect of Age Adjustment on Curated PTAU-Associated Proteins",
    subtitle = sprintf("Orange = Raw ρ | Blue = Age-adjusted ρ | %d proteins", nrow(curated_plot_data)),
    x = "Spearman ρ", y = ""
  ) +
  theme_bw(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 7.5),
    axis.text.y = element_text(face = "italic", size = 6.5)
  )

ggsave(file.path(out_dir, "curated_proteins_rho_shift.pdf"), curated_shift_p,
       width = 10, height = 14, device = "pdf", limitsize = FALSE)
ggsave(file.path(out_dir, "curated_proteins_rho_shift.png"), curated_shift_p,
       width = 10, height = 14, dpi = 300, device = "png", limitsize = FALSE)
cat("Saved: curated_proteins_rho_shift.pdf/png\n")

# ---- 12e. Combined figure ----
combined <- (volcano_p | (rho_shift_p + theme(legend.position = "none"))) /
  (bar_p + theme(legend.position = "none")) +
  plot_layout(heights = c(1, 1.2)) +
  plot_annotation(
    title = "Age-Adjusted PTAU-Protein Association Analysis",
    subtitle = sprintf("%d baseline subjects | AGE effect removed via linear regression residuals",
                       nrow(df_analysis)),
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

ggsave(file.path(out_dir, "combined_age_adjusted_analysis.pdf"), combined,
       width = 18, height = 16, device = "pdf", limitsize = FALSE)
ggsave(file.path(out_dir, "combined_age_adjusted_analysis.png"), combined,
       width = 18, height = 16, dpi = 300, device = "png", limitsize = FALSE)
cat("Saved: combined_age_adjusted_analysis.pdf/png\n")

# ===========================================================================
# 13. Summary
# ===========================================================================
cat("\n")
cat("========== Age-Adjusted Analysis Summary ==========\n")
cat(sprintf("Input: %d baseline samples with non-missing PTAU & AGE\n", nrow(df_analysis)))
cat(sprintf("  AGE range: %.1f - %.1f years\n", min(df_analysis$AGE), max(df_analysis$AGE)))
cat(sprintf("  PTAU ~ AGE: R² = %.4f\n", summary(ptau_lm)$r.squared))
cat(sprintf("Proteins tested: %d\n", nrow(results)))
cat(sprintf("FDR < 0.05: %d (Pos: %d, Neg: %d)\n",
            n_sig_05,
            sum(results$significant_fdr05 & results$direction == "Positive"),
            sum(results$significant_fdr05 & results$direction == "Negative")))
cat(sprintf("FDR < 0.01: %d\n", n_sig_01))

cat(sprintf("\nOverlap with Curated 60 Proteins:\n"))
cat(sprintf("  Curated proteins with FDR < 0.05 (age-adjusted): %d / %d\n",
            nrow(curated_sig), length(curated_ids)))
cat(sprintf("  Top 30 age-adjusted ∩ Curated 60: %d\n", length(overlap_adj_curated)))

cat(sprintf("\nTop 5 age-adjusted (by FDR):\n"))
for (i in 1:min(5, nrow(results))) {
  cat(sprintf("  #%d: %s (%s) | Resid_ρ=%.4f | Raw_ρ=%.4f | FDR=%.2e\n",
              i, results$EntrezGeneSymbol[i], results$Target[i],
              results$residual_spearman_rho[i], results$raw_spearman_rho[i],
              results$fdr[i]))
}

cat(sprintf("\nOutput directory: %s\n", normalizePath(out_dir)))
cat("\n========== Done ==========\n")
