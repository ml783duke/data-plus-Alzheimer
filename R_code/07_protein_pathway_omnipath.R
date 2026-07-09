###############################################################################
# 07_protein_pathway_omnipath.R
# Pathway & Network Analysis of PTAU-Significant Proteins using OmniPath
# - Direct download from omnipathdb.org (OmnipathR package has API issues)
# - Query protein-protein interactions among FDR-significant proteins
# - MAPT-related interactions + hub proteins + network visualization
###############################################################################

library(dplyr)
library(ggplot2)

set.seed(42)

out_dir <- "output/protein_pathway_omnipath"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("========== PTAU-Significant Protein Pathway Analysis (OmniPath) ==========\n")

# ===========================================================================
# 1. Load FDR-significant proteins and map to gene symbols
# ===========================================================================
cat("\n--- Loading significant proteins ---\n")

all_cor <- read.csv("output/ptau_protein_correlation/all_protein_ptau_correlations.csv",
                    stringsAsFactors = FALSE)
cat(sprintf("All proteins: %d\n", nrow(all_cor)))

# Remove potentially broken gene columns, merge fresh
all_cor <- all_cor[, !colnames(all_cor) %in% c("EntrezGeneSymbol", "TargetFullName", "Target", "UniProt")]

dict <- read.csv("protein_raw_data/protein dict.csv", stringsAsFactors = FALSE)
all_cor <- merge(all_cor, dict[, c("Analytes", "EntrezGeneSymbol", "Target", "TargetFullName", "UniProt")],
                 by.x = "protein_id", by.y = "Analytes", all.x = TRUE)

sig <- all_cor[all_cor$significant_fdr05 == TRUE, ]
cat(sprintf("FDR<0.05 proteins: %d\n", nrow(sig)))

sig_genes <- unique(sig$EntrezGeneSymbol[!is.na(sig$EntrezGeneSymbol) & sig$EntrezGeneSymbol != ""])
cat(sprintf("Unique gene symbols in sig set: %d\n", length(sig_genes)))

# Build gene-level rho lookup (median rho per gene, for multi-protein genes)
gene_rho <- sig %>%
  filter(!is.na(EntrezGeneSymbol) & EntrezGeneSymbol != "") %>%
  group_by(EntrezGeneSymbol) %>%
  summarise(
    rho_median = median(spearman_rho, na.rm = TRUE),
    rho_mean   = mean(spearman_rho, na.rm = TRUE),
    n_proteins = n(),
    .groups = "drop"
  )
cat(sprintf("Gene-level rho table: %d genes\n", nrow(gene_rho)))

# ===========================================================================
# 2. Download OmniPath interactions directly
# ===========================================================================
cat("\n--- Downloading OmniPath interactions ---\n")

omni_url <- "https://omnipathdb.org/interactions?format=tsv&fields=sources,references&genesymbols=1"
tmp <- tempfile(fileext = ".tsv")
download.file(omni_url, tmp, method = "auto")
omni <- read.delim(tmp, stringsAsFactors = FALSE)
cat(sprintf("OmniPath interactions: %d\n", nrow(omni)))

# Use source/target as gene symbols (they are HGNC symbols in OmniPath)
# Map our sig genes against them
sig_set <- sig_genes

# Filter: >=1 gene in sig set (using gene symbol columns)
inter_in_sig <- omni[omni$source_genesymbol %in% sig_set | omni$target_genesymbol %in% sig_set, ]
cat(sprintf("Interactions with >=1 sig gene: %d\n", nrow(inter_in_sig)))

# Both sides in sig set
inter_both_sig <- omni[omni$source_genesymbol %in% sig_set & omni$target_genesymbol %in% sig_set, ]
cat(sprintf("Interactions between 2 sig genes: %d\n", nrow(inter_both_sig)))

# MAPT interactions
inter_mapt <- omni[omni$source_genesymbol == "MAPT" | omni$target_genesymbol == "MAPT", ]
cat(sprintf("Total MAPT interactions in OmniPath: %d\n", nrow(inter_mapt)))

# MAPT with sig gene
mapt_sig <- inter_mapt[inter_mapt$source_genesymbol %in% sig_set | inter_mapt$target_genesymbol %in% sig_set, ]
cat(sprintf("MAPT interactions with sig gene: %d\n", nrow(mapt_sig)))

if (nrow(mapt_sig) == 0 && nrow(inter_both_sig) == 0) {
  cat("\nWARNING: No interactions found between sig proteins and MAPT.\n")
  cat("This is expected: OmniPath focuses on signaling interactions,\n")
  cat("and many PTAU-correlated proteins may not directly interact with tau.\n")
  cat("Trying broader search...\n")

  # Fallback: genes near MAPT in interaction space (2-hop)
}

# ===========================================================================
# 3. Build summary tables
# ===========================================================================
cat("\n--- Building summary tables ---\n")

# Table 1: MAPT-interacting sig proteins
if (nrow(mapt_sig) > 0) {
  mapt_partners <- setdiff(unique(c(mapt_sig$source_genesymbol, mapt_sig$target_genesymbol)), "MAPT")

  mapt_tbl <- data.frame(
    gene = mapt_partners,
    stringsAsFactors = FALSE
  )
  mapt_tbl$rho_median <- gene_rho$rho_median[match(mapt_tbl$gene, gene_rho$EntrezGeneSymbol)]
  mapt_tbl$rho_mean <- gene_rho$rho_mean[match(mapt_tbl$gene, gene_rho$EntrezGeneSymbol)]
  mapt_tbl$direction <- ifelse(mapt_tbl$rho_median > 0, "Positive", "Negative")
  mapt_tbl <- mapt_tbl[order(abs(mapt_tbl$rho_median), decreasing = TRUE), ]

  write.csv(mapt_tbl, file.path(out_dir, "mapt_interactors_in_sig.csv"), row.names = FALSE)
  cat(sprintf("MAPT interactors table: %d genes\n", nrow(mapt_tbl)))
  print(mapt_tbl, row.names = FALSE)
}

# Table 2: Hub proteins in the sig-sig network
if (nrow(inter_both_sig) > 0) {
  hub_counts <- table(c(inter_both_sig$source_genesymbol, inter_both_sig$target_genesymbol))
  hub_df <- data.frame(
    gene = names(hub_counts),
    n_interactions = as.integer(hub_counts),
    stringsAsFactors = FALSE
  )
  hub_df$rho_median <- gene_rho$rho_median[match(hub_df$gene, gene_rho$EntrezGeneSymbol)]
  hub_df$rho_mean <- gene_rho$rho_mean[match(hub_df$gene, gene_rho$EntrezGeneSymbol)]
  hub_df <- hub_df[order(hub_df$n_interactions, decreasing = TRUE), ]
  write.csv(hub_df, file.path(out_dir, "hub_proteins_in_network.csv"), row.names = FALSE)
  cat(sprintf("\nHub proteins: %d (from %d sig-sig interactions)\n",
              nrow(hub_df), nrow(inter_both_sig)))
} else {
  hub_df <- data.frame()
}

# ===========================================================================
# 4. Network statistics summary
# ===========================================================================
cat("\n--- Network Statistics ---\n")

# Total sig genes with any OmniPath interaction
sig_with_inter <- unique(c(inter_in_sig$source_genesymbol[inter_in_sig$source_genesymbol %in% sig_set],
                           inter_in_sig$target_genesymbol[inter_in_sig$target_genesymbol %in% sig_set]))
cat(sprintf("Sig genes with >=1 OmniPath interaction: %d / %d (%.1f%%)\n",
            length(sig_with_inter), length(sig_set),
            100 * length(sig_with_inter) / length(sig_set)))

cat(sprintf("Sig genes interacting with MAPT: %d\n",
            length(setdiff(unique(c(mapt_sig$source_genesymbol, mapt_sig$target_genesymbol)), "MAPT"))))

# Top 10 genes by most interactions within sig-sig network
if (nrow(hub_df) > 0 && nrow(hub_df) >= 10) {
  cat("\nTop 10 hub genes (most sig-sig interactions):\n")
  for (i in 1:min(10, nrow(hub_df))) {
    cat(sprintf("  %d. %s: %d interactions (rho=%.4f)\n",
                i, hub_df$gene[i], hub_df$n_interactions[i], hub_df$rho_median[i]))
  }
}

# ===========================================================================
# 5. Visualizations
# ===========================================================================
cat("\n--- Generating figures ---\n")

# Figure 1: MAPT interactors colored by rho
if (nrow(mapt_sig) > 0 && nrow(mapt_tbl) > 0) {
  plot_mapt <- mapt_tbl
  plot_mapt <- plot_mapt[order(plot_mapt$rho_median), ]
  plot_mapt$gene <- factor(plot_mapt$gene, levels = plot_mapt$gene)

  p1 <- ggplot(plot_mapt, aes(x = rho_median, y = gene, fill = direction)) +
    geom_col(alpha = 0.85) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    scale_fill_manual(values = c("Positive" = "#d7191c", "Negative" = "#2c7bb6")) +
    labs(
      title = "MAPT (Tau) Interacting Proteins in PTAU-Significant Set",
      subtitle = sprintf("%d genes interact with MAPT (OmniPath) | FDR<0.05 | Colored by PTAU correlation",
                         nrow(plot_mapt)),
      x = "Median Spearman rho with log2(PTAU)",
      y = ""
    ) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))

  ggsave(file.path(out_dir, "mapt_interactors_rho.png"),
         p1, width = 10, height = max(3, nrow(plot_mapt) * 0.35), dpi = 300, device = "png")
  cat("Saved: mapt_interactors_rho.png\n")
}

# Figure 2: Top hub proteins
if (nrow(hub_df) > 0 && nrow(hub_df) >= 5) {
  top_hubs <- head(hub_df, 20)
  top_hubs <- top_hubs[order(top_hubs$n_interactions), ]
  top_hubs$gene <- factor(top_hubs$gene, levels = top_hubs$gene)

  p2 <- ggplot(top_hubs, aes(x = n_interactions, y = gene)) +
    geom_col(fill = "#2c7bb6", alpha = 0.85) +
    labs(
      title = "Top Hub Proteins in PTAU-Significant Network",
      subtitle = sprintf("Most connected proteins among %d FDR<0.05 proteins (OmniPath PPI)",
                         nrow(sig)),
      x = "Number of Interactions (within sig set)",
      y = ""
    ) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))

  ggsave(file.path(out_dir, "hub_proteins_network.png"),
         p2, width = 9, height = 5, dpi = 300, device = "png")
  cat("Saved: hub_proteins_network.png\n")
}

# Figure 3: Connectivity vs. Correlation scatter
if (nrow(hub_df) > 0) {
  hub_df$is_hub <- hub_df$n_interactions >= quantile(hub_df$n_interactions, 0.9, na.rm = TRUE)
  hub_df$label <- ifelse(hub_df$is_hub | abs(hub_df$rho_median) > 0.7,
                         hub_df$gene, "")

  p3 <- ggplot(hub_df, aes(x = rho_median, y = n_interactions)) +
    geom_point(aes(color = is_hub, size = n_interactions), alpha = 0.7) +
    geom_text(aes(label = label), vjust = -0.8, size = 3) +
    scale_color_manual(values = c("TRUE" = "#d7191c", "FALSE" = "grey60"), guide = "none") +
    scale_size_continuous(range = c(1.5, 5), guide = "none") +
    labs(
      title = "Protein Connectivity vs. PTAU Correlation",
      subtitle = sprintf("%d genes with interactions in sig set | Red = top 10%% hubs",
                         nrow(hub_df)),
      x = "Median Spearman rho with PTAU",
      y = "Number of PPI Interactions (degree)"
    ) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))

  ggsave(file.path(out_dir, "degree_vs_rho.png"),
         p3, width = 9, height = 6, dpi = 300, device = "png")
  cat("Saved: degree_vs_rho.png\n")
}

# Figure 4: Rho distribution of sig genes vs. genes with MAPT interactions
rho_dist <- data.frame(
  rho   = gene_rho$rho_median,
  group = "All Sig",
  stringsAsFactors = FALSE
)
if (nrow(mapt_sig) > 0) {
  mapt_genes_sig <- setdiff(unique(c(mapt_sig$source_genesymbol, mapt_sig$target_genesymbol)), "MAPT")
  rho_dist <- rbind(rho_dist, data.frame(
    rho   = gene_rho$rho_median[gene_rho$EntrezGeneSymbol %in% mapt_genes_sig],
    group = "MAPT Interactors",
    stringsAsFactors = FALSE
  ))
}

p4 <- ggplot(rho_dist, aes(x = rho, fill = group)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("All Sig" = "grey60", "MAPT Interactors" = "#d7191c")) +
  labs(
    title = "PTAU Correlation Distribution: MAPT Interactors vs. All",
    subtitle = "Are MAPT-interacting proteins enriched for stronger PTAU correlations?",
    x = "Median Spearman rho with PTAU",
    y = "Density",
    fill = ""
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

ggsave(file.path(out_dir, "rho_distribution_mapt.png"),
       p4, width = 8, height = 5, dpi = 300, device = "png")
cat("Saved: rho_distribution_mapt.png\n")

# Figure 5: Top 25 PTAU correlates with MAPT interaction status
top25 <- head(sig[order(sig$spearman_rho, decreasing = TRUE), ], 25)
top25$gene_label <- ifelse(
  is.na(top25$EntrezGeneSymbol) | top25$EntrezGeneSymbol == "",
  top25$protein_id, top25$EntrezGeneSymbol
)
top25$mapt_interactor <- top25$EntrezGeneSymbol %in%
  c(mapt_sig$source_genesymbol, mapt_sig$target_genesymbol)
top25 <- top25[!duplicated(top25$gene_label), ]
top25 <- head(top25, 20)
top25 <- top25[order(top25$spearman_rho), ]
top25$gene_label <- factor(top25$gene_label, levels = top25$gene_label)

p5 <- ggplot(top25, aes(x = spearman_rho, y = gene_label, fill = mapt_interactor)) +
  geom_col(alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(
    values = c("TRUE" = "#d7191c", "FALSE" = "grey70"),
    labels = c("TRUE" = "MAPT Interactor (OmniPath)", "FALSE" = "No known MAPT interaction"),
    name = ""
  ) +
  labs(
    title = "Top 20 PTAU Correlates: MAPT Interaction Status",
    subtitle = "OmniPath curated interactions (omnipathdb.org)",
    x = "Spearman rho with log2(PTAU)",
    y = ""
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.y = element_text(face = "italic"),
    legend.position = "bottom"
  )

ggsave(file.path(out_dir, "top20_mapt_overlap.png"),
       p5, width = 9, height = 5.5, dpi = 300, device = "png")
cat("Saved: top20_mapt_overlap.png\n")

# ===========================================================================
# 6. Export network for Cytoscape (if needed)
# ===========================================================================
if (nrow(inter_both_sig) > 0) {
  # Edge list for Cytoscape
  edges <- inter_both_sig[, c("source_genesymbol", "target_genesymbol", "is_directed",
                               "is_stimulation", "is_inhibition", "references")]
  colnames(edges)[1:2] <- c("source", "target")
  write.csv(edges, file.path(out_dir, "sig_sig_edges_cytoscape.csv"), row.names = FALSE)

  # Node attributes
  nodes <- gene_rho[gene_rho$EntrezGeneSymbol %in%
                    unique(c(edges$source, edges$target)), ]
  colnames(nodes)[1] <- "gene"
  write.csv(nodes, file.path(out_dir, "sig_sig_nodes_cytoscape.csv"), row.names = FALSE)
  cat("\nSaved: Cytoscape-compatible edge/node tables\n")
}

# ===========================================================================
# 7. Summary
# ===========================================================================
cat("\n")
cat("========== Analysis Summary ==========\n")
cat(sprintf("Input proteins: %d FDR-significant (%d unique genes)\n",
            nrow(sig), length(sig_genes)))
cat(sprintf("OmniPath total interactions: %d\n", nrow(omni)))
cat(sprintf("Sig genes with any OmniPath connection: %d\n", length(sig_with_inter)))
cat(sprintf("Sig-sig interactions: %d\n", nrow(inter_both_sig)))
cat(sprintf("MAPT interactions in OmniPath: %d total, %d with sig genes\n",
            nrow(inter_mapt), nrow(mapt_sig)))

cat(sprintf("\nOutput directory: %s\n", normalizePath(out_dir)))
cat("\n========== Done ==========\n")
