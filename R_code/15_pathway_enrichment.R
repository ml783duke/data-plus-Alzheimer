###############################################################################
# 15_pathway_enrichment.R
# Data-driven pathway enrichment (GO/KEGG/Reactome) of pTau-significant proteins
# - Full population (all subjects)
# - ORA (Over-Representation Analysis) + GSEA-style visualization
###############################################################################

library(readxl); library(dplyr); library(ggplot2)
library(clusterProfiler); library(org.Hs.eg.db); library(enrichplot)
library(ReactomePA); library(patchwork)

set.seed(42)
out_dir <- "output/pathway_enrichment_full"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("========== Pathway Enrichment: Full Population ==========\n")

# ---- 1. Load & correlate ----
cat("\n--- Loading ---\n")
df <- read.csv("lifestyle_7var_ptau_proteomics.csv")
df <- df[!is.na(df$ptau217_csf) & !is.na(df$dx_entry), ]
protein_cols <- colnames(df)[21:ncol(df)]

# Convert to numeric
for (pc in protein_cols) {
  df[[pc]] <- suppressWarnings(as.numeric(df[[pc]]))
}

# Correlation
df$PTAU_log2 <- log2(df$ptau217_csf)
ptau_vals <- df$PTAU_log2

missing_rate <- sapply(df[protein_cols], function(x) mean(is.na(x))*100)
proteins_kept <- names(missing_rate[missing_rate <= 20])
pl2 <- log2(as.matrix(df[proteins_kept]))
pl2[is.infinite(pl2) | is.nan(pl2)] <- NA

results <- data.frame(protein_id = colnames(pl2), n = NA_integer_,
                      spearman_rho = NA_real_, p_value = NA_real_)
for (i in seq_len(ncol(pl2))) {
  pv <- pl2[, i]; vi <- !is.na(pv); nv <- sum(vi); results$n[i] <- nv
  if (nv >= 10) {
    tst <- tryCatch(cor.test(pv[vi], ptau_vals[vi], method = "spearman", exact = FALSE),
                    error = function(e) NULL)
    if (!is.null(tst)) { results$spearman_rho[i] <- tst$estimate; results$p_value[i] <- tst$p.value }
  }
}
results <- results[!is.na(results$p_value), ]
results$fdr <- p.adjust(results$p_value, method = "BH")
results$GeneSymbol <- gsub("\\.[0-9]+\\.[0-9]+$", "", results$protein_id)

sig <- results[results$fdr < 0.05, ]
sig_genes <- unique(sig$GeneSymbol[!is.na(sig$GeneSymbol) & sig$GeneSymbol != ""])
all_genes <- unique(results$GeneSymbol[!is.na(results$GeneSymbol) & results$GeneSymbol != ""])
cat(sprintf("FDR<0.05: %d proteins -> %d genes (from %d total)\n",
            nrow(sig), length(sig_genes), length(all_genes)))

# ---- 2. Convert to Entrez IDs ----
cat("\n--- Converting to Entrez ---\n")
sig_entrez <- bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
all_entrez <- bitr(all_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
cat(sprintf("Sig Entrez: %d / All Entrez: %d\n", nrow(sig_entrez), nrow(all_entrez)))

# ---- 3. GO BP enrichment ----
cat("\n--- GO BP ---\n")
ego <- enrichGO(gene = sig_entrez$ENTREZID, universe = all_entrez$ENTREZID,
                OrgDb = org.Hs.eg.db, ont = "BP", pvalueCutoff = 0.05,
                qvalueCutoff = 0.2, readable = TRUE)
cat(sprintf("GO BP terms: %d\n", nrow(ego@result)))

if (nrow(ego@result) > 0) {
  # Dotplot
  p_go <- dotplot(ego, showCategory = 20, font.size = 10) +
    labs(title = "GO Biological Process Enrichment", subtitle = "pTau-significant proteins (ORA)")
  ggsave(file.path(out_dir, "go_bp_dotplot.png"), p_go, width = 11, height = 8, dpi = 200)

  # Cnetplot for top 5 terms
  tryCatch({
    p_cnet <- cnetplot(ego, showCategory = 5, circular = FALSE, colorEdge = TRUE) +
      labs(title = "GO BP: Gene-Concept Network (Top 5)")
    ggsave(file.path(out_dir, "go_bp_cnet.png"), p_cnet, width = 14, height = 10, dpi = 200)
  }, error = function(e) message("cnetplot failed: ", e$message))

  write.csv(ego@result, file.path(out_dir, "go_bp_enrichment.csv"), row.names = FALSE)
}

# ---- 4. Reactome enrichment ----
cat("\n--- Reactome ---\n")
reac <- enrichPathway(gene = sig_entrez$ENTREZID, universe = all_entrez$ENTREZID,
                      organism = "human", pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2, readable = TRUE)
cat(sprintf("Reactome terms: %d\n", nrow(reac@result)))

if (nrow(reac@result) > 0) {
  p_reac <- dotplot(reac, showCategory = 20, font.size = 10) +
    labs(title = "Reactome Pathway Enrichment", subtitle = "pTau-significant proteins (ORA)")
  ggsave(file.path(out_dir, "reactome_dotplot.png"), p_reac, width = 11, height = 8, dpi = 200)

  write.csv(reac@result, file.path(out_dir, "reactome_enrichment.csv"), row.names = FALSE)
}

# ---- 5. KEGG enrichment ----
cat("\n--- KEGG ---\n")
kegg <- enrichKEGG(gene = sig_entrez$ENTREZID, universe = all_entrez$ENTREZID,
                   organism = "hsa", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
# Make readable
if (nrow(kegg@result) > 0) {
  kegg <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
} else {
  kegg@result <- data.frame()
}
cat(sprintf("KEGG terms: %d\n", nrow(kegg@result)))

if (nrow(kegg@result) > 0) {
  p_kegg <- dotplot(kegg, showCategory = 20, font.size = 10) +
    labs(title = "KEGG Pathway Enrichment", subtitle = "pTau-significant proteins (ORA)")
  ggsave(file.path(out_dir, "kegg_dotplot.png"), p_kegg, width = 11, height = 8, dpi = 200)

  write.csv(kegg@result, file.path(out_dir, "kegg_enrichment.csv"), row.names = FALSE)
}

# ---- 6. GSEA (ranked by rho, no FDR cut) ----
cat("\n--- GSEA (ranked list) ---\n")

# Build ranked gene list
gene_rho <- results %>%
  filter(!is.na(GeneSymbol) & GeneSymbol != "") %>%
  group_by(GeneSymbol) %>%
  summarise(rho = median(spearman_rho, na.rm = TRUE), .groups = "drop")
ranked <- gene_rho$rho; names(ranked) <- gene_rho$GeneSymbol
ranked <- sort(ranked, decreasing = TRUE)

# Convert to Entrez for GSEA
rank_entrez_df <- bitr(names(ranked), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
rank_entrez <- ranked[rank_entrez_df$SYMBOL]
names(rank_entrez) <- rank_entrez_df$ENTREZID[match(names(rank_entrez), rank_entrez_df$SYMBOL)]
rank_entrez <- sort(rank_entrez, decreasing = TRUE)

gsea_go <- gseGO(geneList = rank_entrez, OrgDb = org.Hs.eg.db, ont = "BP",
                 pvalueCutoff = 0.05, eps = 0)
cat(sprintf("GSEA GO terms: %d\n", nrow(gsea_go@result)))

if (nrow(gsea_go@result) > 0) {
  # Ridgeplot
  p_ridge <- ridgeplot(gsea_go, showCategory = 15) +
    labs(title = "GSEA GO BP: pTau Correlation (Ranked)", subtitle = "Activated (red) vs Suppressed (blue)")
  ggsave(file.path(out_dir, "gsea_ridgeplot.png"), p_ridge, width = 11, height = 8, dpi = 200)

  # GSEA plot (top 6)
  p_gsea <- gseaplot2(gsea_go, geneSetID = 1:6, pvalue_table = TRUE)
  ggsave(file.path(out_dir, "gsea_running_score.png"), p_gsea, width = 12, height = 10, dpi = 200)

  write.csv(gsea_go@result, file.path(out_dir, "gsea_go_bp.csv"), row.names = FALSE)
}

# ---- 7. Summary table of top pathways ----
cat("\n--- Summary ---\n")
top_paths <- data.frame()
for (src in c("GO_BP", "Reactome", "KEGG")) {
  obj <- if (src == "GO_BP") ego else if (src == "Reactome") reac else kegg
  if (nrow(obj@result) > 0) {
    top <- head(obj@result[, c("ID", "Description", "p.adjust", "Count")], 5)
    top$Source <- src
    top_paths <- rbind(top_paths, top)
  }
}
if (nrow(top_paths) > 0) {
  top_paths <- top_paths[order(top_paths$p.adjust), ]
  write.csv(top_paths, file.path(out_dir, "top_pathways_summary.csv"), row.names = FALSE)
  cat("\nTop enriched pathways across all databases:\n")
  print(top_paths, row.names = FALSE)
}

cat(sprintf("\nOutput: %s\n", normalizePath(out_dir)))
cat("========== Done ==========\n")
