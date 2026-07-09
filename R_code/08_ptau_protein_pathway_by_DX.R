###############################################################################
# 08_ptau_protein_pathway_by_DX.R
# PTAU-Protein Correlation + OmniPath Pathway Analysis STRATIFIED by DX
# - Within each DX group: Spearman correlation → FDR → OmniPath MAPT/hub analysis
# - Compare: which MAPT interactors are significant only in certain stages?
###############################################################################

library(readxl)
library(dplyr)
library(ggplot2)
library(igraph)
library(ggraph)

set.seed(42)

out_dir <- "output/ptau_protein_pathway_by_DX"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

DX_GROUPS <- c("CN", "EMCI", "LMCI", "AD")

cat("========== PTAU-Protein Pathway Analysis by DX ==========\n")

# ===========================================================================
# 1. Load and prepare data
# ===========================================================================
cat("\n--- Loading data ---\n")
master <- read_excel("master_data.xlsx", sheet = "Sheet1")
dict <- read.csv("protein_raw_data/protein dict.csv", stringsAsFactors = FALSE)

df_bl <- master[master$VISCODE2 == "bl", ]
df_bl <- df_bl[!is.na(df_bl$PTAU) & !is.na(df_bl$DX), ]

protein_cols <- grep("^X[0-9]+\\.[0-9]+$", colnames(df_bl), value = TRUE)
cat(sprintf("Baseline with PTAU+DX: %d subjects, %d protein columns\n",
            nrow(df_bl), length(protein_cols)))

# Convert proteins
convert_protein <- function(x) { x[x == "NA" | x == ""] <- NA; as.numeric(x) }
df_bl[protein_cols] <- lapply(df_bl[protein_cols], convert_protein)

# Download OmniPath once
cat("Downloading OmniPath interactions...\n")
omni_url <- "https://omnipathdb.org/interactions?format=tsv&fields=sources,references&genesymbols=1"
tmp <- tempfile(fileext = ".tsv")
download.file(omni_url, tmp, method = "auto")
omni <- read.delim(tmp, stringsAsFactors = FALSE)

# ===========================================================================
# 2. Per-DX analysis function
# ===========================================================================
run_dx_pathway <- function(df_dx, dx_label) {
  cat(sprintf("\n========== DX = %s (%d subjects) ==========\n", dx_label, nrow(df_dx)))

  # --- Protein filtering ---
  missing_rate <- sapply(df_dx[protein_cols], function(x) mean(is.na(x)) * 100)
  proteins_kept <- names(missing_rate[missing_rate <= 20])
  cat(sprintf("  Proteins kept: %d / %d\n", length(proteins_kept), length(protein_cols)))

  # --- Log2 transform ---
  df_dx$PTAU_log2 <- log2(df_dx$PTAU)
  protein_mat <- as.matrix(df_dx[proteins_kept])
  protein_log2 <- log2(protein_mat)
  protein_log2[is.infinite(protein_log2) | is.nan(protein_log2)] <- NA

  # --- Spearman correlation ---
  ptau_vals <- df_dx$PTAU_log2
  n_prot <- ncol(protein_log2)
  results <- data.frame(protein_id = colnames(protein_log2), n = NA_integer_,
                        spearman_rho = NA_real_, p_value = NA_real_,
                        stringsAsFactors = FALSE)
  for (i in seq_len(n_prot)) {
    pv <- protein_log2[, i]; vi <- !is.na(pv); nv <- sum(vi)
    results$n[i] <- nv
    if (nv >= 10) {
      tst <- tryCatch(cor.test(pv[vi], ptau_vals[vi], method = "spearman", exact = FALSE),
                      error = function(e) NULL)
      if (!is.null(tst)) { results$spearman_rho[i] <- tst$estimate; results$p_value[i] <- tst$p.value }
    }
  }
  results <- results[!is.na(results$p_value), ]
  results$fdr <- p.adjust(results$p_value, method = "BH")
  results$significant_fdr05 <- results$fdr < 0.05

  # Merge gene symbols
  results <- merge(results,
                   dict[, c("Analytes", "EntrezGeneSymbol", "Target", "TargetFullName")],
                   by.x = "protein_id", by.y = "Analytes", all.x = TRUE)

  n_sig <- sum(results$significant_fdr05)
  cat(sprintf("  FDR<0.05: %d / %d\n", n_sig, nrow(results)))

  if (n_sig < 10) {
    cat("  WARNING: Too few significant proteins, skipping OmniPath analysis\n")
    return(list(results = results, n_sig = n_sig, mapt_genes = character(0),
                hub_genes = character(0)))
  }

  # --- Gene-level summary ---
  sig <- results[results$significant_fdr05, ]
  sig_genes <- unique(sig$EntrezGeneSymbol[!is.na(sig$EntrezGeneSymbol) & sig$EntrezGeneSymbol != ""])
  cat(sprintf("  Sig genes: %d\n", length(sig_genes)))

  gene_rho <- sig %>%
    filter(!is.na(EntrezGeneSymbol) & EntrezGeneSymbol != "") %>%
    group_by(EntrezGeneSymbol) %>%
    summarise(rho_median = median(spearman_rho, na.rm = TRUE), .groups = "drop")

  # --- MAPT interactors in this DX ---
  mapt_dx <- omni[
    (omni$source_genesymbol == "MAPT" | omni$target_genesymbol == "MAPT") &
    (omni$source_genesymbol %in% sig_genes | omni$target_genesymbol %in% sig_genes),
  ]
  mapt_genes_dx <- setdiff(unique(c(mapt_dx$source_genesymbol, mapt_dx$target_genesymbol)), "MAPT")

  # Annotate with rho
  if (length(mapt_genes_dx) > 0) {
    mapt_rho_dx <- data.frame(gene = mapt_genes_dx, stringsAsFactors = FALSE)
    mapt_rho_dx$rho <- gene_rho$rho_median[match(mapt_rho_dx$gene, gene_rho$EntrezGeneSymbol)]
    mapt_rho_dx$DX <- dx_label
  } else {
    mapt_rho_dx <- data.frame()
  }

  # --- Hub proteins ---
  inter_sig <- omni[
    omni$source_genesymbol %in% sig_genes & omni$target_genesymbol %in% sig_genes,
  ]
  if (nrow(inter_sig) > 0) {
    hc <- table(c(inter_sig$source_genesymbol, inter_sig$target_genesymbol))
    hub_dx <- data.frame(gene = names(hc), n_int = as.integer(hc), stringsAsFactors = FALSE)
    hub_dx$rho <- gene_rho$rho_median[match(hub_dx$gene, gene_rho$EntrezGeneSymbol)]
    hub_dx <- hub_dx[order(hub_dx$n_int, decreasing = TRUE), ]
  } else {
    hub_dx <- data.frame()
  }

  cat(sprintf("  MAPT interactors: %d | Sig-sig interactions: %d\n",
              length(mapt_genes_dx), nrow(inter_sig)))

  list(
    results    = results,
    n_sig      = n_sig,
    n_sig_genes = length(sig_genes),
    mapt_genes = mapt_genes_dx,
    mapt_rho   = mapt_rho_dx,
    hub_df     = hub_dx,
    gene_rho   = gene_rho
  )
}

# ===========================================================================
# 3. Run per DX group
# ===========================================================================
dx_pathway <- list()
for (g in DX_GROUPS) {
  df_g <- df_bl[df_bl$DX == g, ]
  if (nrow(df_g) < 30) { cat(sprintf("DX=%s: insufficient N (%d), skipping\n", g, nrow(df_g))); next }
  dx_pathway[[g]] <- run_dx_pathway(df_g, g)
}

# ===========================================================================
# 4. Cross-DX comparison of MAPT interactors
# ===========================================================================
cat("\n========== Cross-DX MAPT Interactor Comparison ==========\n")

# Collect all MAPT genes across groups
all_mapt <- unique(unlist(lapply(dx_pathway, function(x) x$mapt_genes)))
cat(sprintf("Total unique MAPT interactors across DX: %d\n", length(all_mapt)))

# Build comparison matrix
mapt_compare <- data.frame(gene = all_mapt, stringsAsFactors = FALSE)
for (g in DX_GROUPS) {
  res <- dx_pathway[[g]]
  if (is.null(res)) next
  mapt_compare[[paste0("sig_", g)]] <- mapt_compare$gene %in% res$mapt_genes
  mapt_compare[[paste0("rho_", g)]] <- NA_real_
  for (i in seq_len(nrow(mapt_compare))) {
    gn <- mapt_compare$gene[i]
    if (gn %in% res$mapt_genes) {
      mapt_compare[[paste0("rho_", g)]][i] <- res$mapt_rho$rho[match(gn, res$mapt_rho$gene)]
    }
  }
}
mapt_compare$n_groups <- rowSums(mapt_compare[, paste0("sig_", DX_GROUPS)], na.rm = TRUE)
mapt_compare <- mapt_compare[order(mapt_compare$n_groups, decreasing = TRUE), ]

write.csv(mapt_compare, file.path(out_dir, "mapt_interactors_by_DX.csv"), row.names = FALSE)

# Print summary
cat(sprintf("\nMAPT interactors found in:\n"))
for (k in 4:1) {
  cat(sprintf("  %d groups: %d genes\n", k, sum(mapt_compare$n_groups == k)))
  if (k >= 2 && sum(mapt_compare$n_groups == k) > 0) {
    cat(sprintf("    %s\n", paste(mapt_compare$gene[mapt_compare$n_groups == k], collapse = ", ")))
  }
}

# ===========================================================================
# 5. Visualizations
# ===========================================================================
cat("\n--- Generating figures ---\n")

# Figure 1: MAPT interactor rho heatmap across DX
if (nrow(mapt_compare) > 0) {
  # Select genes found in >=1 group with rho data
  plot_genes <- mapt_compare$gene[mapt_compare$n_groups >= 1]
  if (length(plot_genes) > 40) plot_genes <- mapt_compare$gene[mapt_compare$n_groups >= 2]

  hm_data <- data.frame()
  for (g in DX_GROUPS) {
    hm_data <- rbind(hm_data, data.frame(
      gene = plot_genes,
      DX   = g,
      rho  = mapt_compare[[paste0("rho_", g)]][match(plot_genes, mapt_compare$gene)],
      sig  = mapt_compare[[paste0("sig_", g)]][match(plot_genes, mapt_compare$gene)],
      stringsAsFactors = FALSE
    ))
  }
  # Sort genes by average rho
  gene_order <- hm_data %>% group_by(gene) %>%
    summarise(avg_rho = mean(rho, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(avg_rho))
  hm_data$gene <- factor(hm_data$gene, levels = gene_order$gene)

  p1 <- ggplot(hm_data, aes(x = DX, y = gene, fill = rho)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = ifelse(sig, sprintf("%.2f", rho), "·")), size = 3) +
    scale_fill_gradient2(low = "#2c7bb6", mid = "white", high = "#d7191c",
                          midpoint = 0, na.value = "grey90", name = "PTAU rho") +
    labs(
      title = "MAPT Interactors: PTAU Correlation by Disease Stage",
      subtitle = sprintf("%d MAPT-interacting genes | · = not FDR-significant in that group",
                         length(plot_genes)),
      x = "", y = ""
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          axis.text.y = element_text(face = "italic", size = 8))

  ggsave(file.path(out_dir, "mapt_heatmap_by_DX.png"),
         p1, width = 8, height = max(6, length(plot_genes) * 0.3), dpi = 200, device = "png")
  cat("Saved: mapt_heatmap_by_DX.png\n")
}

# Figure 2: Top 10 hubs per DX group comparison
hub_all <- data.frame()
for (g in DX_GROUPS) {
  res <- dx_pathway[[g]]
  if (is.null(res) || nrow(res$hub_df) == 0) next
  top10 <- head(res$hub_df, 10)
  top10$DX <- g
  top10$rank <- 1:nrow(top10)
  hub_all <- rbind(hub_all, top10)
}

if (nrow(hub_all) > 0) {
  hub_all$label <- paste0(hub_all$DX, " #", hub_all$rank, ": ", hub_all$gene)
  hub_all <- hub_all[order(hub_all$n_int), ]
  hub_all$label <- factor(hub_all$label, levels = hub_all$label)

  p2 <- ggplot(hub_all, aes(x = n_int, y = label, fill = rho)) +
    geom_col(alpha = 0.85) +
    scale_fill_gradient2(low = "#2c7bb6", mid = "white", high = "#d7191c",
                          midpoint = 0, name = "PTAU rho") +
    facet_wrap(~DX, scales = "free_y", ncol = 1) +
    labs(
      title = "Top 10 Hub Proteins per Disease Stage",
      subtitle = "Most connected proteins among PTAU-significant set",
      x = "Number of Interactions", y = ""
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

  ggsave(file.path(out_dir, "hub_comparison_by_DX.png"),
         p2, width = 10, height = 12, dpi = 200, device = "png")
  cat("Saved: hub_comparison_by_DX.png\n")
}

# Figure 3: Hub overlap network (genes that are hub in multiple DX groups)
hub_genes_all <- lapply(DX_GROUPS, function(g) {
  res <- dx_pathway[[g]]
  if (is.null(res) || nrow(res$hub_df) == 0) return(character(0))
  head(res$hub_df$gene, 15)
})
names(hub_genes_all) <- DX_GROUPS

# Count how many DX groups each hub gene appears in
all_hub_genes <- unique(unlist(hub_genes_all))
hub_overlap <- data.frame(gene = all_hub_genes, stringsAsFactors = FALSE)
for (g in DX_GROUPS) {
  hub_overlap[[g]] <- hub_overlap$gene %in% hub_genes_all[[g]]
}
hub_overlap$n_dx <- rowSums(hub_overlap[, DX_GROUPS])

# Get avg rho across groups
hub_overlap$avg_rho <- sapply(hub_overlap$gene, function(gn) {
  rho_vals <- sapply(DX_GROUPS, function(g) {
    res <- dx_pathway[[g]]
    if (is.null(res)) return(NA)
    res$hub_df$rho[match(gn, res$hub_df$gene)]
  })
  mean(rho_vals, na.rm = TRUE)
})

hub_overlap <- hub_overlap[order(hub_overlap$n_dx, decreasing = TRUE), ]

if (sum(hub_overlap$n_dx >= 2) > 0) {
  hub_label <- hub_overlap$gene
  hub_label[hub_overlap$n_dx < 2] <- ""

  p3 <- ggplot(hub_overlap, aes(x = CN, y = EMCI, color = avg_rho, size = n_dx)) +
    geom_jitter(width = 0.1, height = 0.1, alpha = 0.8) +
    geom_text(aes(label = ifelse(n_dx >= 2, gene, "")), vjust = -1, size = 3, color = "black") +
    scale_color_gradient2(low = "#2c7bb6", mid = "white", high = "#d7191c",
                           midpoint = 0, name = "Avg PTAU rho") +
    scale_size_continuous(range = c(2, 6), name = "DX groups") +
    labs(
      title = "Hub Protein Overlap: CN vs EMCI",
      subtitle = "Each point = a top-15 hub gene | Size = # DX groups",
      x = "Hub in CN", y = "Hub in EMCI"
    ) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

  ggsave(file.path(out_dir, "hub_overlap_scatter.png"),
         p3, width = 8, height = 7, dpi = 200, device = "png")
  cat("Saved: hub_overlap_scatter.png\n")
}

# ===========================================================================
# 6. Summary table
# ===========================================================================
cat("\n========== Per-DX Summary ==========\n")
summary_tbl <- data.frame(DX = DX_GROUPS, stringsAsFactors = FALSE)
summary_tbl$N <- sapply(DX_GROUPS, function(g) sum(df_bl$DX == g))
summary_tbl$Proteins_tested <- sapply(DX_GROUPS, function(g) {
  if (is.null(dx_pathway[[g]])) return(NA); nrow(dx_pathway[[g]]$results)
})
summary_tbl$FDR_sig <- sapply(DX_GROUPS, function(g) {
  if (is.null(dx_pathway[[g]])) return(NA); dx_pathway[[g]]$n_sig
})
summary_tbl$Sig_genes <- sapply(DX_GROUPS, function(g) {
  if (is.null(dx_pathway[[g]])) return(NA); dx_pathway[[g]]$n_sig_genes
})
summary_tbl$MAPT_interactors <- sapply(DX_GROUPS, function(g) {
  if (is.null(dx_pathway[[g]])) return(NA); length(dx_pathway[[g]]$mapt_genes)
})

print(summary_tbl)
write.csv(summary_tbl, file.path(out_dir, "per_dx_summary.csv"), row.names = FALSE)

cat(sprintf("\nOutput directory: %s\n", normalizePath(out_dir)))
cat("\n========== Done ==========\n")
