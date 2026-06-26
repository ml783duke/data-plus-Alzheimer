###############################################################################
# 01_ptau_protein_correlation.R
# Spearman correlation analysis between SomaScan proteins and PTAU
# - Baseline (bl) samples only
# - Log2 transformation of PTAU and protein RFU values
# - FDR correction (Benjamini-Hochberg)
# - Volcano plot & top-hits visualization
###############################################################################

# ===========================================================================
# 0. Setup: libraries
# ===========================================================================
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(patchwork)

set.seed(42)

# Output directory
out_dir <- "output/ptau_protein_correlation"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ===========================================================================
# 1. Load data
# ===========================================================================
cat("\n========== Loading data ==========\n")
df_raw <- read_excel("master_data.xlsx", sheet = "Sheet1")
cat(sprintf("Raw data: %d rows x %d columns\n", nrow(df_raw), ncol(df_raw)))

# Load protein dictionary for annotation
dict <- read.csv("protein_raw_data/protein dict.csv", stringsAsFactors = FALSE)
cat(sprintf("Protein dictionary: %d entries\n", nrow(dict)))

# ===========================================================================
# 2. Select baseline samples only
# ===========================================================================
cat("\n========== Filtering to baseline ==========\n")

# Report visit distribution
cat("VISCODE2 distribution:\n")
print(table(df_raw$VISCODE2, useNA = "always"))

# Baseline only
df_bl <- df_raw[df_raw$VISCODE2 == "bl", ]
cat(sprintf("Baseline samples: %d rows\n", nrow(df_bl)))

# How many have PTAU?
cat(sprintf("Baseline with non-missing PTAU: %d\n", sum(!is.na(df_bl$PTAU))))

# ===========================================================================
# 3. Identify columns
# ===========================================================================
# Protein columns (SomaScan, starting with X followed by digits and dot)
protein_cols <- grep("^X[0-9]+\\.[0-9]+$", colnames(df_bl), value = TRUE)
cat(sprintf("SomaScan protein columns identified: %d\n", length(protein_cols)))

# Metadata columns to keep for reference
meta_cols <- c("RID", "PTAU", "TAU", "ABETA42", "ABETA40", "ptau181/ab42",
               "DX", "AGE", "PTGENDER", "PTEDUCAT", "APOE4", "FDG", "CDRSB", "MMSE")

# ===========================================================================
# 4. Convert protein columns to numeric (they are character with "NA" strings)
# ===========================================================================
cat("\n========== Converting protein columns to numeric ==========\n")

convert_protein <- function(x) {
  x[x == "NA" | x == ""] <- NA
  as.numeric(x)
}

# Count how many values are "NA" strings before conversion
na_string_counts <- sapply(df_bl[protein_cols[1:10]], function(x) sum(x == "NA", na.rm = TRUE))
cat("Example 'NA' string counts in first 10 protein columns:\n")
print(na_string_counts)

# Convert all protein columns
df_bl[protein_cols] <- lapply(df_bl[protein_cols], convert_protein)

# Verify conversion
cat("\nAfter conversion - sample protein column types:\n")
print(sapply(df_bl[protein_cols[1:5]], class))

# ===========================================================================
# 5. Filter proteins by missingness
# ===========================================================================
cat("\n========== Filtering proteins by missingness (>20%) ==========\n")

# Calculate missing rate per protein (among baseline samples with non-missing PTAU)
df_analysis <- df_bl[!is.na(df_bl$PTAU), ]
cat(sprintf("Analysis set (non-missing PTAU): %d samples\n", nrow(df_analysis)))

missing_rate <- sapply(df_analysis[protein_cols], function(x) mean(is.na(x)) * 100)

# Proteins to remove
proteins_high_missing <- names(missing_rate[missing_rate > 20])
proteins_kept <- names(missing_rate[missing_rate <= 20])

cat(sprintf("Proteins with >20%% missing: %d (will be removed)\n", length(proteins_high_missing)))
cat(sprintf("Proteins with <=20%% missing: %d (kept for analysis)\n", length(proteins_kept)))

# Show top missing proteins
cat("\nTop 10 proteins with highest missing rate:\n")
print(head(sort(missing_rate, decreasing = TRUE), 10))

# Save list of removed proteins for reference
removed_df <- data.frame(
  protein_id = proteins_high_missing,
  missing_rate_pct = round(missing_rate[proteins_high_missing], 2),
  stringsAsFactors = FALSE
)

# ===========================================================================
# 6. Log2 transformation
# ===========================================================================
cat("\n========== Log2 transformation ==========\n")

# PTAU: add 1 pseudo-count to avoid log(0), but PTAU is >0 so log2 directly
cat(sprintf("PTAU range before log2: %.2f - %.2f\n", min(df_analysis$PTAU), max(df_analysis$PTAU)))
df_analysis$PTAU_log2 <- log2(df_analysis$PTAU)
cat(sprintf("PTAU_log2 range: %.2f - %.2f\n", min(df_analysis$PTAU_log2), max(df_analysis$PTAU_log2)))

# Proteins: log2 transform (values are RFU, always >0 after removing NAs)
cat("Log2 transforming protein values...\n")
protein_mat <- as.matrix(df_analysis[proteins_kept])
protein_log2 <- log2(protein_mat)

# Check for any issues
neg_count <- sum(protein_mat <= 0, na.rm = TRUE)
cat(sprintf("Protein values <= 0: %d (will produce -Inf/NaN in log2)\n", neg_count))
if (neg_count > 0) {
  cat("  Replacing -Inf/NaN with NA\n")
  protein_log2[is.infinite(protein_log2) | is.nan(protein_log2)] <- NA
}

cat(sprintf("Log2 protein matrix: %d x %d\n", nrow(protein_log2), ncol(protein_log2)))

# ===========================================================================
# 7. Spearman correlation: each protein vs PTAU_log2
# ===========================================================================
cat("\n========== Running Spearman correlations ==========\n")

ptau_vals <- df_analysis$PTAU_log2
n_proteins <- ncol(protein_log2)

# Pre-allocate results
results <- data.frame(
  protein_id   = colnames(protein_log2),
  n             = NA_integer_,
  spearman_rho = NA_real_,
  p_value      = NA_real_,
  stringsAsFactors = FALSE
)

for (i in seq_len(n_proteins)) {
  prot_vals <- protein_log2[, i]
  valid_idx  <- !is.na(prot_vals)
  n_valid    <- sum(valid_idx)
  results$n[i] <- n_valid

  if (n_valid >= 10) {  # require at least 10 pairs
    test <- tryCatch(
      cor.test(prot_vals[valid_idx], ptau_vals[valid_idx],
               method = "spearman", exact = FALSE),
      error = function(e) NULL
    )
    if (!is.null(test)) {
      results$spearman_rho[i] <- test$estimate
      results$p_value[i]      <- test$p.value
    }
  }
}

# Remove proteins that couldn't be tested
results <- results[!is.na(results$p_value), ]
cat(sprintf("Proteins successfully tested: %d / %d\n", nrow(results), n_proteins))

# ===========================================================================
# 8. FDR correction (Benjamini-Hochberg)
# ===========================================================================
cat("\n========== FDR correction ==========\n")

results$fdr <- p.adjust(results$p_value, method = "BH")

# Significance at two thresholds
results$significant_fdr05 <- results$fdr < 0.05
results$significant_fdr01 <- results$fdr < 0.01

# Direction
results$direction <- ifelse(results$spearman_rho > 0, "Positive", "Negative")

n_sig_05 <- sum(results$significant_fdr05)
n_sig_01 <- sum(results$significant_fdr01)
cat(sprintf("Significant at FDR < 0.05: %d proteins\n", n_sig_05))
cat(sprintf("Significant at FDR < 0.01: %d proteins\n", n_sig_01))

if (n_sig_05 > 0) {
  cat(sprintf("  Positive correlations: %d\n", sum(results$significant_fdr05 & results$direction == "Positive")))
  cat(sprintf("  Negative correlations: %d\n", sum(results$significant_fdr05 & results$direction == "Negative")))
}

# ===========================================================================
# 9. Annotate results with gene names from protein dictionary
# ===========================================================================
cat("\n========== Annotating with gene symbols ==========\n")

# Create lookup: SeqId -> GeneSymbol, TargetFullName
# Use Analytes column ("X10000.28" format) matching protein column names
annot_lookup <- dict[, c("Analytes", "EntrezGeneSymbol", "TargetFullName", "Target", "UniProt")]
colnames(annot_lookup)[1] <- "protein_id"

results <- merge(results, annot_lookup, by = "protein_id", all.x = TRUE, sort = FALSE)

# Sort by FDR
results <- results[order(results$fdr), ]

# ===========================================================================
# 10. Print top results
# ===========================================================================
cat("\n========== Top 20 proteins (by FDR) ==========\n")
top20 <- head(results, 20)
print(data.frame(
  Gene        = top20$EntrezGeneSymbol,
  Protein     = top20$Target,
  Rho         = round(top20$spearman_rho, 3),
  P           = formatC(top20$p_value, format = "e", digits = 2),
  FDR         = formatC(top20$fdr, format = "e", digits = 2),
  N           = top20$n,
  Direction   = top20$direction
))

# ===========================================================================
# 11. Summary report
# ===========================================================================
cat("\n")
cat("========== Analysis Summary ==========\n")
cat(sprintf("Input: %d baseline samples, %d with non-missing PTAU\n",
            nrow(df_bl), nrow(df_analysis)))
cat(sprintf("Proteins tested: %d\n", nrow(results)))
cat(sprintf("Proteins removed (>20%% NA): %d\n", length(proteins_high_missing)))
cat(sprintf("FDR < 0.05: %d significant proteins\n", n_sig_05))
cat(sprintf("FDR < 0.01: %d significant proteins\n", n_sig_01))
cat(sprintf("Rho range: %.3f to %.3f\n", min(results$spearman_rho), max(results$spearman_rho)))

# ===========================================================================
# 12. Visualization
# ===========================================================================
cat("\n========== Generating plots ==========\n")

# --- Helper: get display label (gene symbol preferred, protein_id fallback) ---
get_gene_label <- function(df) {
  label <- ifelse(is.na(df$EntrezGeneSymbol) | df$EntrezGeneSymbol == "",
                  df$protein_id, df$EntrezGeneSymbol)
  # If same gene appears multiple times (different SeqIds), add protein_id prefix
  dup_genes <- label[duplicated(label) | duplicated(label, fromLast = TRUE)]
  for (g in unique(dup_genes)) {
    idx <- which(label == g)
    label[idx] <- paste0(g, " (", df$protein_id[idx], ")")
  }
  return(label)
}

# --- Volcano Plot ---
plot_data <- results
plot_data$neg_log10_fdr <- -log10(plot_data$fdr)

# Cap very small FDR for plotting
fdr_cap <- 1e-16
plot_data$neg_log10_fdr_capped <- pmin(plot_data$neg_log10_fdr, -log10(fdr_cap))

# Significance category for coloring
plot_data$sig_cat <- "Not significant"
plot_data$sig_cat[plot_data$significant_fdr01] <- "FDR < 0.01"
plot_data$sig_cat[plot_data$significant_fdr05 & !plot_data$significant_fdr01] <- "FDR < 0.05"
plot_data$sig_cat <- factor(plot_data$sig_cat, levels = c("FDR < 0.01", "FDR < 0.05", "Not significant"))

# Label top 15 significant proteins (by FDR) with gene symbols
plot_data$label <- ""
top_label <- head(plot_data[plot_data$significant_fdr05, ], 15)
top_label$gene_label <- get_gene_label(top_label)
plot_data$label[match(top_label$protein_id, plot_data$protein_id)] <- top_label$gene_label

volcano_p <- ggplot(plot_data, aes(x = spearman_rho, y = neg_log10_fdr_capped,
                                    color = sig_cat, label = label)) +
  geom_point(alpha = 0.6, size = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", alpha = 0.6) +
  geom_hline(yintercept = -log10(0.01), linetype = "dotted", color = "grey40", alpha = 0.4) +
  geom_text_repel(
    size = 2.8,
    max.overlaps = 20,
    box.padding = 0.4,
    point.padding = 0.2,
    segment.size = 0.2,
    color = "black",
    fontface = "italic"
  ) +
  scale_color_manual(
    values = c("FDR < 0.01" = "#E41A1C",
               "FDR < 0.05" = "#377EB8",
               "Not significant" = "grey70")
  ) +
  labs(
    title = "Proteins Associated with PTAU",
    subtitle = sprintf("Spearman correlation | %d baseline samples | %d proteins tested\nFDR < 0.05: %d sig | FDR < 0.01: %d sig",
                       nrow(df_analysis), nrow(results), n_sig_05, n_sig_01),
    x = expression("Spearman" ~ rho),
    y = expression(-log[10]("FDR")),
    color = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = c(0.85, 0.85),
    legend.background = element_rect(fill = "white", color = "grey80", size = 0.3),
    legend.margin = margin(4, 6, 4, 4),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 8, color = "grey40")
  )

ggsave(file.path(out_dir, "volcano_ptau_proteins.pdf"), volcano_p,
       width = 9, height = 7, device = "pdf")
ggsave(file.path(out_dir, "volcano_ptau_proteins.png"), volcano_p,
       width = 9, height = 7, dpi = 300)

# --- Top 30 Positive (rho > 0, sorted by FDR) ---
pos30 <- head(results[results$spearman_rho > 0, ], 30)
pos30$gene_label <- get_gene_label(pos30)
pos30$gene_label <- factor(pos30$gene_label, levels = rev(pos30$gene_label))
pos30$neg_log10_fdr <- -log10(pos30$fdr)

bar_pos <- ggplot(pos30, aes(x = neg_log10_fdr, y = gene_label, fill = neg_log10_fdr)) +
  geom_col() +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey40", alpha = 0.6) +
  geom_vline(xintercept = -log10(0.01), linetype = "dotted", color = "grey40", alpha = 0.4) +
  scale_fill_gradient(low = "#FDD49E", high = "#B2182B", name = expression(-log[10](FDR))) +
  labs(
    title = sprintf("Top 30 Positively Correlated (n=%d)", nrow(pos30)),
    subtitle = "Higher protein = Higher PTAU | Dashed: FDR=0.05, Dotted: FDR=0.01",
    x = expression(-log[10]("FDR")),
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 8, color = "grey40"),
    axis.text.y = element_text(face = "italic", size = 9)
  )

ggsave(file.path(out_dir, "top30_positive_barplot.pdf"), bar_pos,
       width = 10, height = 9, device = "pdf")
ggsave(file.path(out_dir, "top30_positive_barplot.png"), bar_pos,
       width = 10, height = 9, dpi = 300)

# --- Top 30 Negative (rho < 0, sorted by FDR) ---
neg30 <- head(results[results$spearman_rho < 0, ], 30)
neg30$gene_label <- get_gene_label(neg30)
neg30$gene_label <- factor(neg30$gene_label, levels = rev(neg30$gene_label))
neg30$neg_log10_fdr <- -log10(neg30$fdr)

bar_neg <- ggplot(neg30, aes(x = neg_log10_fdr, y = gene_label, fill = neg_log10_fdr)) +
  geom_col() +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey40", alpha = 0.6) +
  geom_vline(xintercept = -log10(0.01), linetype = "dotted", color = "grey40", alpha = 0.4) +
  scale_fill_gradient(low = "#D1E5F0", high = "#2166AC", name = expression(-log[10](FDR))) +
  labs(
    title = sprintf("Top 30 Negatively Correlated (n=%d)", nrow(neg30)),
    subtitle = "Higher protein = Lower PTAU | Dashed: FDR=0.05, Dotted: FDR=0.01",
    x = expression(-log[10]("FDR")),
    y = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 8, color = "grey40"),
    axis.text.y = element_text(face = "italic", size = 9)
  )

ggsave(file.path(out_dir, "top30_negative_barplot.pdf"), bar_neg,
       width = 10, height = 9, device = "pdf")
ggsave(file.path(out_dir, "top30_negative_barplot.png"), bar_neg,
       width = 10, height = 9, dpi = 300)

# --- Combined figure: Volcano + Positive + Negative ---
combined <- (volcano_p | (bar_pos / bar_neg)) +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(
    title = "SomaScan Proteomics vs log2(PTAU) — Spearman Correlation Analysis",
    caption = sprintf("Generated: %s | %d baseline samples | %d proteins tested",
                      Sys.Date(), nrow(df_analysis), nrow(results))
  )
ggsave(file.path(out_dir, "combined_ptau_protein_analysis.pdf"), combined,
       width = 18, height = 14, device = "pdf")
ggsave(file.path(out_dir, "combined_ptau_protein_analysis.png"), combined,
       width = 18, height = 14, dpi = 300)

cat(sprintf("Plots saved to %s/\n", out_dir))

# ===========================================================================
# 13. Save results table
# ===========================================================================
cat("\n========== Saving results ==========\n")

# Full results (sorted by FDR)
write.csv(results, file.path(out_dir, "all_protein_ptau_correlations.csv"),
          row.names = FALSE)

# Significant only (FDR < 0.05)
sig_results <- results[results$significant_fdr05, ]
write.csv(sig_results, file.path(out_dir, "significant_proteins_fdr05.csv"),
          row.names = FALSE)

# Top 30 positive
write.csv(pos30[, c("protein_id", "EntrezGeneSymbol", "TargetFullName",
                     "spearman_rho", "p_value", "fdr", "n")],
          file.path(out_dir, "top30_positive_correlations.csv"),
          row.names = FALSE)

# Top 30 negative
write.csv(neg30[, c("protein_id", "EntrezGeneSymbol", "TargetFullName",
                     "spearman_rho", "p_value", "fdr", "n")],
          file.path(out_dir, "top30_negative_correlations.csv"),
          row.names = FALSE)

# Removed proteins
write.csv(removed_df, file.path(out_dir, "removed_proteins_high_missing.csv"),
          row.names = FALSE)

cat("Results saved:\n")
cat("  - all_protein_ptau_correlations.csv (full results)\n")
cat("  - significant_proteins_fdr05.csv (FDR<0.05 hits)\n")
cat("  - top30_positive_correlations.csv (top 30 positive)\n")
cat("  - top30_negative_correlations.csv (top 30 negative)\n")
cat("  - removed_proteins_high_missing.csv (filtered proteins)\n")

cat("\n========== Done ==========\n")
