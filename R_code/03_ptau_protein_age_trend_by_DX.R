###############################################################################
# 03_ptau_protein_age_trend_by_DX.R
# Age-trend tracking analysis STRATIFIED by disease stage (DX)
# - EMCI, LMCI, AD only (exclude CN)
# - Within each DX group: LOESS trend correlation (protein vs PTAU across AGE)
# - Compare top tracking proteins across disease stages
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
LOESS_SPAN       <- 0.75
MIN_VALID_PAIRS  <- 20     # lower threshold for smaller DX subgroups
TOP_N_PLOT       <- 10     # proteins per group in individual plots
TOP_N_COMPARE    <- 20     # top N per group for cross-group comparison

# DX groups to analyze
DX_GROUPS <- c("EMCI", "LMCI", "AD")

# Output directory
out_dir <- "output/ptau_protein_age_trend_by_DX"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("DX groups: %s\n", paste(DX_GROUPS, collapse = ", ")))
cat(sprintf("LOESS span: %.2f | Min valid pairs: %d\n", LOESS_SPAN, MIN_VALID_PAIRS))

# ===========================================================================
# 1. Load data
# ===========================================================================
cat("\n========== Loading data ==========\n")
df_raw <- read_excel("master_data.xlsx", sheet = "Sheet1")
cat(sprintf("Raw data: %d rows x %d columns\n", nrow(df_raw), ncol(df_raw)))

dict <- read.csv("protein_raw_data/protein dict.csv", stringsAsFactors = FALSE)
cat(sprintf("Protein dictionary: %d entries\n", nrow(dict)))

# ===========================================================================
# 2. Filter to baseline, non-missing PTAU & AGE, target DX groups
# ===========================================================================
cat("\n========== Filtering samples ==========\n")

df_bl <- df_raw[df_raw$VISCODE2 == "bl", ]
cat(sprintf("Baseline samples: %d\n", nrow(df_bl)))

# Full DX distribution
cat("DX distribution (baseline):\n")
print(table(df_bl$DX, useNA = "always"))

# Filter to non-missing PTAU & AGE
df_analysis <- df_bl[!is.na(df_bl$PTAU) & !is.na(df_bl$AGE), ]
cat(sprintf("Baseline with non-missing PTAU & AGE: %d\n", nrow(df_analysis)))

# Filter to target DX groups
df_target <- df_analysis[df_analysis$DX %in% DX_GROUPS, ]
cat(sprintf("Filtered to EMCI + LMCI + AD: %d samples\n", nrow(df_target)))

# Per-group counts
for (g in DX_GROUPS) {
  n_g <- sum(df_target$DX == g)
  cat(sprintf("  %s: %d\n", g, n_g))
}

# ===========================================================================
# 3. Identify and convert protein columns (once, globally)
# ===========================================================================
cat("\n========== Processing protein columns ==========\n")
protein_cols <- grep("^X[0-9]+\\.[0-9]+$", colnames(df_target), value = TRUE)
cat(sprintf("SomaScan protein columns: %d\n", length(protein_cols)))

convert_protein <- function(x) {
  x[x == "NA" | x == ""] <- NA
  as.numeric(x)
}

df_target[protein_cols] <- lapply(df_target[protein_cols], convert_protein)

# ===========================================================================
# 4. Protein dictionary lookup
# ===========================================================================
annot_lookup <- dict[, c("Analytes", "EntrezGeneSymbol", "TargetFullName", "Target", "UniProt")]
colnames(annot_lookup)[1] <- "protein_id"

# ===========================================================================
# 5. Helper functions
# ===========================================================================

#' Run age-trend analysis for one DX group
run_dx_analysis <- function(df_group, dx_label) {
  cat(sprintf("\n========== DX = %s (%d samples) ==========\n", dx_label, nrow(df_group)))

  # --- 5a. Filter proteins by missingness within this group ---
  missing_rate <- sapply(df_group[protein_cols], function(x) mean(is.na(x)) * 100)
  proteins_kept    <- names(missing_rate[missing_rate <= 20])
  proteins_removed <- names(missing_rate[missing_rate > 20])
  cat(sprintf("  Proteins kept (<=20%% NA): %d | Removed: %d\n",
              length(proteins_kept), length(proteins_removed)))

  # --- 5b. Log2 transform ---
  df_group$PTAU_log2 <- log2(df_group$PTAU)
  protein_mat  <- as.matrix(df_group[proteins_kept])
  protein_log2 <- log2(protein_mat)
  protein_log2[is.infinite(protein_log2) | is.nan(protein_log2)] <- NA

  age_vals  <- df_group$AGE
  ptau_log2 <- df_group$PTAU_log2

  # --- 5c. PTAU LOESS ---
  ptau_valid <- !is.na(age_vals) & !is.na(ptau_log2)
  cat(sprintf("  Valid (AGE, PTAU) pairs: %d\n", sum(ptau_valid)))

  ptau_loess <- tryCatch(
    loess(ptau_log2[ptau_valid] ~ age_vals[ptau_valid], span = LOESS_SPAN),
    error = function(e) NULL
  )
  if (is.null(ptau_loess)) {
    cat("  ERROR: PTAU LOESS failed!\n")
    return(NULL)
  }

  ptau_fitted <- rep(NA_real_, nrow(df_group))
  ptau_fitted[ptau_valid] <- fitted(ptau_loess)

  # PTAU trend direction
  ptau_age_cor <- cor.test(ptau_fitted[ptau_valid], age_vals[ptau_valid],
                           method = "spearman", exact = FALSE)
  cat(sprintf("  PTAU vs AGE: Spearman rho = %.4f (p = %.2e)\n",
              ptau_age_cor$estimate, ptau_age_cor$p.value))

  # --- 5d. Per-protein LOESS + trend correlation ---
  n_proteins <- ncol(protein_log2)
  results <- data.frame(
    protein_id      = colnames(protein_log2),
    n_valid         = NA_integer_,
    loess_corr      = NA_real_,
    loess_mae_z     = NA_real_,
    loess_converged = FALSE,
    stringsAsFactors = FALSE
  )

  for (i in seq_len(n_proteins)) {
    prot_log2 <- protein_log2[, i]
    valid_idx <- !is.na(age_vals) & !is.na(prot_log2)
    n_valid   <- sum(valid_idx)
    results$n_valid[i] <- n_valid

    if (n_valid >= MIN_VALID_PAIRS) {
      eff_span <- max(LOESS_SPAN, 15 / n_valid)

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
          }
        }

        # Z-score MAE
        if (!is.na(results$loess_corr[i])) {
          prot_fitted_z <- (prot_fitted - mean(prot_fitted, na.rm = TRUE)) /
                           sd(prot_fitted, na.rm = TRUE)
          ptau_fitted_z <- (ptau_fitted - mean(ptau_fitted, na.rm = TRUE)) /
                           sd(ptau_fitted, na.rm = TRUE)
          shared_z <- shared_valid[!is.na(prot_fitted_z[shared_valid])]
          if (length(shared_z) >= 10) {
            results$loess_mae_z[i] <- mean(abs(
              ptau_fitted_z[shared_z] - prot_fitted_z[shared_z]
            ), na.rm = TRUE)
          }
        }
      }
    }
  }

  # Filter successful and sort
  results_valid <- results[!is.na(results$loess_corr), ]
  results_valid <- results_valid[order(results_valid$loess_corr, decreasing = TRUE), ]
  results_valid$rank <- seq_len(nrow(results_valid))

  # Annotate
  results_valid <- merge(results_valid, annot_lookup, by = "protein_id",
                          all.x = TRUE, sort = FALSE)
  results_valid <- results_valid[order(results_valid$loess_corr, decreasing = TRUE), ]
  results_valid$rank <- seq_len(nrow(results_valid))

  cat(sprintf("  Proteins successfully ranked: %d / %d\n",
              nrow(results_valid), n_proteins))
  cat(sprintf("  loess_corr range: %.4f to %.4f (median %.4f)\n",
              min(results_valid$loess_corr), max(results_valid$loess_corr),
              median(results_valid$loess_corr)))

  # Return list with all needed objects
  list(
    dx_label        = dx_label,
    n_samples       = nrow(df_group),
    ptau_age_rho    = ptau_age_cor$estimate,
    ptau_age_p      = ptau_age_cor$p.value,
    results         = results_valid,
    age_vals        = age_vals,
    ptau_log2       = ptau_log2,
    ptau_fitted     = ptau_fitted,
    protein_log2    = protein_log2,
    proteins_kept   = proteins_kept
  )
}

# ===========================================================================
# 6. Run analysis for each DX group
# ===========================================================================
dx_results <- list()
for (g in DX_GROUPS) {
  df_g <- df_target[df_target$DX == g, ]
  res  <- run_dx_analysis(df_g, g)
  dx_results[[g]] <- res
}

# ===========================================================================
# 7. Save per-group CSVs
# ===========================================================================
cat("\n========== Saving per-group CSVs ==========\n")
for (g in DX_GROUPS) {
  res <- dx_results[[g]]
  if (!is.null(res)) {
    fname <- sprintf("age_trend_ranking_%s.csv", g)
    write.csv(res$results, file.path(out_dir, fname), row.names = FALSE)
    cat(sprintf("Saved: %s (%d proteins)\n", fname, nrow(res$results)))
  }
}

# ===========================================================================
# 8. Cross-group comparison
# ===========================================================================
cat("\n========== Cross-Group Comparison ==========\n")

# Build a comparison table: for each protein, show its rank in each group
all_protein_ids <- unique(unlist(lapply(dx_results, function(r) {
  if (!is.null(r)) r$results$protein_id else NULL
})))
cat(sprintf("Total unique proteins across all groups: %d\n", length(all_protein_ids)))

comparison <- data.frame(
  protein_id = all_protein_ids,
  stringsAsFactors = FALSE
)

for (g in DX_GROUPS) {
  res <- dx_results[[g]]
  if (!is.null(res)) {
    comparison[[paste0("rank_", g)]] <- res$results$rank[match(
      all_protein_ids, res$results$protein_id
    )]
    comparison[[paste0("loess_corr_", g)]] <- res$results$loess_corr[match(
      all_protein_ids, res$results$protein_id
    )]
  }
}

# Compute intersection: proteins in top TOP_N_COMPARE of ALL groups
top_sets <- list()
for (g in DX_GROUPS) {
  res <- dx_results[[g]]
  if (!is.null(res)) {
    top_sets[[g]] <- head(res$results$protein_id, TOP_N_COMPARE)
    cat(sprintf("Top %d proteins in %s: %d\n", TOP_N_COMPARE, g, length(top_sets[[g]])))
  }
}

# Intersection across all 3 groups
top_intersect <- Reduce(intersect, top_sets)
cat(sprintf("\nProteins in Top %d of ALL 3 groups (EMCI ∩ LMCI ∩ AD): %d\n",
            TOP_N_COMPARE, length(top_intersect)))

# Pairwise overlaps
for (i in 1:(length(DX_GROUPS)-1)) {
  for (j in (i+1):length(DX_GROUPS)) {
    g1 <- DX_GROUPS[i]; g2 <- DX_GROUPS[j]
    overlap <- intersect(top_sets[[g1]], top_sets[[g2]])
    cat(sprintf("  %s ∩ %s: %d\n", g1, g2, length(overlap)))
  }
}

# Annotate and print the conserved proteins (in intersection)
if (length(top_intersect) > 0) {
  # Get annotations from any group's results
  annot_sub <- annot_lookup[annot_lookup$protein_id %in% top_intersect, ]

  cat("\n--- Conserved Top-Tracking Proteins (in Top 20 of ALL 3 groups) ---\n")
  for (pid in top_intersect) {
    gene <- annot_sub$EntrezGeneSymbol[annot_sub$protein_id == pid]
    target <- annot_sub$Target[annot_sub$protein_id == pid]
    ranks <- sapply(DX_GROUPS, function(g) {
      idx <- which(dx_results[[g]]$results$protein_id == pid)
      if (length(idx) > 0) dx_results[[g]]$results$rank[idx] else NA
    })
    corrs <- sapply(DX_GROUPS, function(g) {
      idx <- which(dx_results[[g]]$results$protein_id == pid)
      if (length(idx) > 0) round(dx_results[[g]]$results$loess_corr[idx], 4) else NA
    })
    cat(sprintf("  %s (%s): Ranks %s | Corrs %s\n",
                pid, gene,
                paste(sprintf("%s=#%d", DX_GROUPS, ranks), collapse = ", "),
                paste(sprintf("%s=%.4f", DX_GROUPS, corrs), collapse = ", ")))
  }
}

# Save comparison table
write.csv(comparison, file.path(out_dir, "cross_group_comparison.csv"), row.names = FALSE)
cat(sprintf("\nSaved: cross_group_comparison.csv (%d proteins)\n", nrow(comparison)))

# ===========================================================================
# 9. Visualization: Per-group combined plots
# ===========================================================================
cat("\n========== Generating per-group trend plots ==========\n")

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

for (g in DX_GROUPS) {
  res <- dx_results[[g]]
  if (is.null(res)) next

  top10 <- head(res$results, TOP_N_PLOT)
  top10$gene_label <- get_gene_label(top10)

  age_vals  <- res$age_vals
  ptau_log2 <- res$ptau_log2
  ptau_fitted <- res$ptau_fitted

  # PTAU z-scores
  ptau_raw_z    <- (ptau_log2 - mean(ptau_log2, na.rm = TRUE)) / sd(ptau_log2, na.rm = TRUE)
  ptau_fitted_z <- (ptau_fitted - mean(ptau_fitted, na.rm = TRUE)) / sd(ptau_fitted, na.rm = TRUE)

  # Build line data
  line_list <- list()
  line_list[[1]] <- data.frame(
    AGE      = age_vals,
    fitted_z = ptau_fitted_z,
    series   = "PTAU",
    stringsAsFactors = FALSE
  )

  for (k in seq_len(nrow(top10))) {
    pid <- top10$protein_id[k]
    pidx <- match(pid, colnames(res$protein_log2))
    if (is.na(pidx)) next
    prot_log2v <- res$protein_log2[, pidx]
    valid_idx <- !is.na(age_vals) & !is.na(prot_log2v)
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

    gene_lbl <- top10$gene_label[k]
    line_list[[k + 1]] <- data.frame(
      AGE      = age_vals,
      fitted_z = prot_fitted_z,
      series   = sprintf("%d. %s", k, gene_lbl),
      stringsAsFactors = FALSE
    )
  }

  plot_lines <- do.call(rbind, line_list)
  plot_lines$series <- factor(plot_lines$series, levels = unique(plot_lines$series))

  # Colors
  n_prot_lines <- length(unique(plot_lines$series)) - 1
  if (n_prot_lines > 0) {
    prot_colors <- hcl.colors(n_prot_lines, "Dynamic")
    names(prot_colors) <- setdiff(unique(plot_lines$series), "PTAU")
    all_colors <- c("PTAU" = "black", prot_colors)
  } else {
    all_colors <- c("PTAU" = "black")
  }

  p <- ggplot() +
    geom_point(data = data.frame(AGE = age_vals, value_z = ptau_raw_z),
               aes(x = AGE, y = value_z),
               alpha = 0.15, size = 0.6, color = "grey50") +
    geom_line(data = plot_lines,
              aes(x = AGE, y = fitted_z, color = series), linewidth = 1.0) +
    scale_color_manual(values = all_colors) +
    labs(
      title = sprintf("Proteins Tracking PTAU Age-Trend — %s", g),
      subtitle = sprintf("LOESS span = %.2f | Z-scored | %d subjects | Top %d proteins",
                         LOESS_SPAN, res$n_samples, TOP_N_PLOT),
      x = "Age (years)",
      y = "Z-score of LOESS Fitted Values",
      color = "Series"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.text  = element_text(size = 7.5),
      plot.title   = element_text(face = "bold"),
      plot.subtitle = element_text(size = 7.5, color = "grey40")
    )

  ggsave(file.path(out_dir, sprintf("trend_plot_%s.pdf", g)),
         p, width = 10, height = 7, device = "pdf")
  ggsave(file.path(out_dir, sprintf("trend_plot_%s.png", g)),
         p, width = 10, height = 7, dpi = 300, device = "png")
  cat(sprintf("Saved: trend_plot_%s.pdf/png\n", g))
}

# ===========================================================================
# 10. Cross-group comparison visualization
# ===========================================================================
cat("\n========== Generating cross-group comparison plot ==========\n")

# Plot A: PTAU trend in each DX group (z-scored, overlaid)
ptau_trend_list <- list()
for (g in DX_GROUPS) {
  res <- dx_results[[g]]
  if (is.null(res)) next
  age_vals <- res$age_vals
  ptau_fitted <- res$ptau_fitted
  ptau_fitted_z <- (ptau_fitted - mean(ptau_fitted, na.rm = TRUE)) / sd(ptau_fitted, na.rm = TRUE)
  ptau_trend_list[[g]] <- data.frame(
    AGE      = age_vals,
    fitted_z = ptau_fitted_z,
    DX       = g,
    stringsAsFactors = FALSE
  )
}
ptau_trends <- do.call(rbind, ptau_trend_list)

# Plot B: Multi-panel comparison — top proteins in intersection
# Show the conserved proteins' LOESS curves across all 3 groups

# Use top 5 from intersection (or fewer if less overlap)
n_conserved <- min(5, length(top_intersect))
if (n_conserved > 0) {
  conserved_pids <- top_intersect[1:n_conserved]

  conserved_plot_list <- list()
  for (pid in conserved_pids) {
    gene <- annot_lookup$EntrezGeneSymbol[annot_lookup$protein_id == pid]
    if (is.na(gene) || gene == "") gene <- pid

    # Collect fitted values across groups
    cdata_list <- list()
    for (g in DX_GROUPS) {
      res <- dx_results[[g]]
      if (is.null(res)) next
      pidx <- match(pid, colnames(res$protein_log2))
      if (is.na(pidx)) next

      prot_log2v <- res$protein_log2[, pidx]
      age_vals <- res$age_vals
      valid_idx <- !is.na(age_vals) & !is.na(prot_log2v)
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

      cdata_list[[g]] <- data.frame(
        AGE      = age_vals,
        fitted_z = prot_fitted_z,
        Series   = g,
        stringsAsFactors = FALSE
      )
    }

    if (length(cdata_list) > 0) {
      cdata <- do.call(rbind, cdata_list)

      # Also add PTAU trend for each group as dashed line
      ptau_overlay <- ptau_trends
      ptau_overlay$Series <- ptau_overlay$DX

      p <- ggplot() +
        geom_line(data = cdata, aes(x = AGE, y = fitted_z, color = Series),
                  linewidth = 1.0) +
        scale_color_manual(values = c("EMCI" = "#2c7bb6", "LMCI" = "#fdae61", "AD" = "#d7191c")) +
        labs(
          title = sprintf("%s (%s)", gene, pid),
          subtitle = "LOESS fitted trend across disease stages",
          x = "Age (years)",
          y = "Z-score"
        ) +
        theme_bw(base_size = 10) +
        theme(plot.title = element_text(face = "bold"))

      conserved_plot_list[[pid]] <- p
    }
  }

  if (length(conserved_plot_list) > 0) {
    ncol_cp <- min(3, length(conserved_plot_list))
    cp_multi <- wrap_plots(conserved_plot_list, ncol = ncol_cp) +
      plot_annotation(
        title = "Conserved Top-Tracking Proteins Across Disease Stages",
        subtitle = sprintf("Proteins in Top %d of EMCI, LMCI, and AD",
                           TOP_N_COMPARE),
        theme = theme(plot.title = element_text(face = "bold", size = 13))
      )

    ggsave(file.path(out_dir, "conserved_proteins_cross_DX.pdf"),
           cp_multi, width = 14, height = 8, device = "pdf", limitsize = FALSE)
    ggsave(file.path(out_dir, "conserved_proteins_cross_DX.png"),
           cp_multi, width = 14, height = 8, dpi = 300, device = "png", limitsize = FALSE)
    cat("Saved: conserved_proteins_cross_DX.pdf/png\n")
  }
}

# Plot C: Heatmap-style — top protein ranks across groups
# Show top 30 from each group in a matrix

top_n_heatmap <- 30
heatmap_pids <- unique(unlist(lapply(top_sets, function(s) s[1:min(top_n_heatmap, length(s))])))
# For each of these, get rank in each group
heatmap_data <- data.frame(
  protein_id = heatmap_pids,
  stringsAsFactors = FALSE
)
for (g in DX_GROUPS) {
  res <- dx_results[[g]]
  if (!is.null(res)) {
    heatmap_data[[paste0("rank_", g)]] <- sapply(heatmap_pids, function(pid) {
      idx <- which(res$results$protein_id == pid)
      if (length(idx) > 0) res$results$rank[idx] else NA
    })
  }
}

# Add gene symbols
heatmap_data <- merge(heatmap_data, annot_lookup[, c("protein_id", "EntrezGeneSymbol")],
                       by = "protein_id", all.x = TRUE)
heatmap_data$label <- ifelse(is.na(heatmap_data$EntrezGeneSymbol) | heatmap_data$EntrezGeneSymbol == "",
                              heatmap_data$protein_id, heatmap_data$EntrezGeneSymbol)

# Compute average rank across groups for sorting
heatmap_data$avg_rank <- rowMeans(
  heatmap_data[, paste0("rank_", DX_GROUPS)], na.rm = TRUE
)
heatmap_data <- heatmap_data[order(heatmap_data$avg_rank), ]
heatmap_data <- head(heatmap_data, 40)  # top 40 by average rank

# Reshape for ggplot
heatmap_long <- data.frame()
for (g in DX_GROUPS) {
  heatmap_long <- rbind(heatmap_long, data.frame(
    protein_id = heatmap_data$protein_id,
    label      = heatmap_data$label,
    DX         = g,
    rank       = heatmap_data[[paste0("rank_", g)]],
    stringsAsFactors = FALSE
  ))
}
heatmap_long$label <- factor(heatmap_long$label,
                              levels = rev(heatmap_data$label))
# Cap rank for display
heatmap_long$rank_capped <- pmin(heatmap_long$rank, 200)

hm_p <- ggplot(heatmap_long, aes(x = DX, y = label, fill = rank_capped)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(is.na(rank), "-", as.character(rank))),
            size = 2.5, color = "grey20") +
  scale_fill_gradientn(
    colors = c("#d7191c", "#fdae61", "#ffffbf", "#abd9e9", "#2c7bb6"),
    values = scales::rescale(c(1, 5, 20, 50, 200)),
    na.value = "grey90",
    name = "Rank\n(capped at 200)"
  ) +
  labs(
    title = "Top Protein Ranks Across Disease Stages",
    subtitle = sprintf("Top %d proteins per group (sorted by average rank across groups)", top_n_heatmap),
    x = "",
    y = ""
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.y = element_text(face = "italic", size = 7),
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(out_dir, "rank_heatmap_cross_DX.pdf"),
       hm_p, width = 8, height = 10, device = "pdf")
ggsave(file.path(out_dir, "rank_heatmap_cross_DX.png"),
       hm_p, width = 8, height = 10, dpi = 300, device = "png")
cat("Saved: rank_heatmap_cross_DX.pdf/png\n")

# ===========================================================================
# 11. PTAU trend comparison plot
# ===========================================================================
ptau_trend_p <- ggplot(ptau_trends, aes(x = AGE, y = fitted_z, color = DX)) +
  geom_line(linewidth = 1.2) +
  scale_color_manual(values = c("EMCI" = "#2c7bb6", "LMCI" = "#fdae61", "AD" = "#d7191c")) +
  labs(
    title = "PTAU LOESS Trend Across Age — By Disease Stage",
    subtitle = sprintf("Z-scored fitted values | EMCI: %d, LMCI: %d, AD: %d subjects",
                       dx_results[["EMCI"]]$n_samples,
                       dx_results[["LMCI"]]$n_samples,
                       dx_results[["AD"]]$n_samples),
    x = "Age (years)",
    y = "Z-score of LOESS-Fitted log2(PTAU)",
    color = "Disease Stage"
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(out_dir, "ptau_trend_by_DX.pdf"),
       ptau_trend_p, width = 8, height = 5.5, device = "pdf")
ggsave(file.path(out_dir, "ptau_trend_by_DX.png"),
       ptau_trend_p, width = 8, height = 5.5, dpi = 300, device = "png")
cat("Saved: ptau_trend_by_DX.pdf/png\n")

# ===========================================================================
# 12. Summary report
# ===========================================================================
cat("\n")
cat("========== Cross-Group Analysis Summary ==========\n")
for (g in DX_GROUPS) {
  res <- dx_results[[g]]
  if (!is.null(res)) {
    cat(sprintf("\n%s (n=%d):\n", g, res$n_samples))
    cat(sprintf("  PTAU vs AGE: rho = %.4f (p = %.2e)\n", res$ptau_age_rho, res$ptau_age_p))
    cat(sprintf("  Proteins ranked: %d\n", nrow(res$results)))
    top5 <- head(res$results, 5)
    for (i in 1:nrow(top5)) {
      gs <- top5$EntrezGeneSymbol[i]
      if (is.na(gs) || gs == "") gs <- top5$protein_id[i]
      cat(sprintf("    #%d: %s (ρ=%.4f)\n", i, gs, top5$loess_corr[i]))
    }
  }
}

cat(sprintf("\nConserved (in Top %d of all 3 groups): %d proteins\n",
            TOP_N_COMPARE, length(top_intersect)))
cat(sprintf("Pairwise overlaps: EMCI∩LMCI=%d, EMCI∩AD=%d, LMCI∩AD=%d\n",
            length(intersect(top_sets[["EMCI"]], top_sets[["LMCI"]])),
            length(intersect(top_sets[["EMCI"]], top_sets[["AD"]])),
            length(intersect(top_sets[["LMCI"]], top_sets[["AD"]]))))

cat(sprintf("\nOutput directory: %s\n", normalizePath(out_dir)))
cat("\n========== Done ==========\n")
