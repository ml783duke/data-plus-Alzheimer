###############################################################################
# 16_pathway_enrichment_by_DX.R
# Data-driven pathway enrichment (GO/Reactome) stratified by DX
# - Per DX: ORA of pTau-significant proteins
# - Compare enriched pathways across disease stages
###############################################################################

library(dplyr); library(ggplot2); library(clusterProfiler)
library(org.Hs.eg.db); library(enrichplot); library(ReactomePA); library(patchwork)

set.seed(42)
out_dir <- "output/pathway_enrichment_by_DX"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("========== Pathway Enrichment by DX ==========\n")

# ---- 1. Load & prepare ----
df <- read.csv("lifestyle_7var_ptau_proteomics.csv")
df <- df[!is.na(df$ptau217_csf) & !is.na(df$dx_entry), ]
protein_cols <- colnames(df)[21:ncol(df)]
for (pc in protein_cols) df[[pc]] <- suppressWarnings(as.numeric(df[[pc]]))

DX_GROUPS <- c("CN", "EMCI", "LMCI", "AD")
all_enrich <- list()

# ---- 2. Per-DX enrichment ----
for (g in DX_GROUPS) {
  cat(sprintf("\n========== %s ==========\n", g))
  df_g <- df[df$dx_entry == g, ]
  if (nrow(df_g) < 50) { cat(sprintf("  N=%d, skipping\n", nrow(df_g))); next }

  cat(sprintf("  N=%d subjects\n", nrow(df_g)))

  # Correlation
  df_g$PTAU_log2 <- log2(df_g$ptau217_csf)
  mr <- sapply(df_g[protein_cols], function(x) mean(is.na(x))*100)
  pk <- names(mr[mr <= 30])  # relaxed for small groups
  pl2 <- log2(as.matrix(df_g[pk])); pl2[is.infinite(pl2)|is.nan(pl2)] <- NA

  results <- data.frame(protein_id=colnames(pl2), n=NA_integer_,
                        spearman_rho=NA_real_, p_value=NA_real_)
  for (i in seq_len(ncol(pl2))) {
    pv <- pl2[,i]; vi <- !is.na(pv); nv <- sum(vi); results$n[i] <- nv
    if (nv >= 10) {
      tst <- tryCatch(cor.test(pv[vi], df_g$PTAU_log2[vi], method="spearman", exact=FALSE),
                      error=function(e) NULL)
      if (!is.null(tst)) { results$spearman_rho[i] <- tst$estimate; results$p_value[i] <- tst$p.value }
    }
  }
  results <- results[!is.na(results$p_value), ]
  results$fdr <- p.adjust(results$p_value, method="BH")
  results$gene <- gsub("\\.[0-9]+\\.[0-9]+$", "", results$protein_id)

  sig <- results[results$fdr < 0.05, ]
  sig_genes <- unique(sig$gene[!is.na(sig$gene) & sig$gene != ""])
  all_genes <- unique(results$gene[!is.na(results$gene) & results$gene != ""])
  cat(sprintf("  Sig: %d proteins -> %d genes (from %d)\n", nrow(sig), length(sig_genes), length(all_genes)))

  # Convert to Entrez
  sig_entrez <- tryCatch(bitr(sig_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db),
                         error=function(e) NULL)
  all_entrez <- tryCatch(bitr(all_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db),
                         error=function(e) NULL)
  if (is.null(sig_entrez) || nrow(sig_entrez) < 5 || is.null(all_entrez)) next

  # GO BP
  ego <- enrichGO(gene=sig_entrez$ENTREZID, universe=all_entrez$ENTREZID,
                  OrgDb=org.Hs.eg.db, ont="BP", pvalueCutoff=0.05, qvalueCutoff=0.2, readable=TRUE)

  # Reactome
  reac <- enrichPathway(gene=sig_entrez$ENTREZID, universe=all_entrez$ENTREZID,
                        organism="human", pvalueCutoff=0.05, qvalueCutoff=0.2, readable=TRUE)

  # Store top terms
  for (src_name in c("GO_BP")) {
    obj <- ego
    if (nrow(obj@result) > 0) {
      top <- head(obj@result[order(obj@result$p.adjust), ], 10)
      top$DX <- g; top$Source <- src_name
      all_enrich[[paste(g, src_name)]] <- top
    }
  }
  for (src_name in c("Reactome")) {
    obj <- reac
    if (nrow(obj@result) > 0) {
      top <- head(obj@result[order(obj@result$p.adjust), ], 10)
      top$DX <- g; top$Source <- src_name
      all_enrich[[paste(g, src_name)]] <- top
    }
  }

  cat(sprintf("  GO: %d | Reactome: %d terms\n",
              nrow(ego@result), nrow(reac@result)))
}

# ---- 3. Combine & compare ----
cat("\n========== Cross-DX Comparison ==========\n")

if (length(all_enrich) > 0) {
  enrich_all <- do.call(rbind, all_enrich)
  write.csv(enrich_all, file.path(out_dir, "enrichment_by_DX.csv"), row.names = FALSE)

  # Select top 5 per DX from GO BP
  go_dx <- enrich_all[enrich_all$Source == "GO_BP", ]
  go_dx <- go_dx %>%
    group_by(DX) %>%
    slice_min(p.adjust, n = 5) %>%
    ungroup()
  go_dx$label <- paste0(go_dx$DX, ": ", substr(go_dx$Description, 1, 55))
  go_dx <- go_dx[order(go_dx$p.adjust, decreasing = TRUE), ]
  go_dx$label <- factor(go_dx$label, levels = go_dx$label)

  p_go_dx <- ggplot(go_dx, aes(x = -log10(p.adjust), y = label, size = Count, color = DX)) +
    geom_point(alpha = 0.85) +
    scale_color_manual(values=c("CN"="#2c7bb6","EMCI"="#fdae61","LMCI"="#1b9e77","AD"="#d7191c")) +
    scale_size_continuous(range=c(2,7)) +
    labs(title="GO BP Enrichment by Disease Stage",
         subtitle="Top 5 pathways per DX | ORA of pTau-significant proteins",
         x="-log10(adjusted P)", y="") +
    theme_bw(11) + theme(plot.title=element_text(face="bold"))
  ggsave(file.path(out_dir, "go_bp_by_DX.png"), p_go_dx, width=11, height=8, dpi=200)

  # Reactome by DX
  rx_dx <- enrich_all[enrich_all$Source == "Reactome", ]
  rx_dx <- rx_dx %>%
    group_by(DX) %>%
    slice_min(p.adjust, n = 5) %>%
    ungroup()
  rx_dx$label <- paste0(rx_dx$DX, ": ", substr(rx_dx$Description, 1, 55))
  rx_dx <- rx_dx[order(rx_dx$p.adjust, decreasing = TRUE), ]
  rx_dx$label <- factor(rx_dx$label, levels = rx_dx$label)

  p_rx_dx <- ggplot(rx_dx, aes(x = -log10(p.adjust), y = label, size = Count, color = DX)) +
    geom_point(alpha = 0.85) +
    scale_color_manual(values=c("CN"="#2c7bb6","EMCI"="#fdae61","LMCI"="#1b9e77","AD"="#d7191c")) +
    scale_size_continuous(range=c(2,7)) +
    labs(title="Reactome Pathway Enrichment by Disease Stage",
         subtitle="Top 5 pathways per DX | ORA of pTau-significant proteins",
         x="-log10(adjusted P)", y="") +
    theme_bw(11) + theme(plot.title=element_text(face="bold"))
  ggsave(file.path(out_dir, "reactome_by_DX.png"), p_rx_dx, width=11, height=8, dpi=200)

  # ---- Heatmap: top terms shared/specific ----
  # Collect top pathways shared across >=2 DX groups
  term_counts <- table(enrich_all$Description)
  shared_terms <- names(term_counts[term_counts >= 2])

  if (length(shared_terms) > 0) {
    hm_data <- enrich_all[enrich_all$Description %in% shared_terms, ]
    hm_data <- hm_data %>%
      group_by(Description) %>%
      slice_min(p.adjust, n = 1) %>%
      ungroup()
    # Limit to top 30
    hm_data <- hm_data[order(hm_data$p.adjust), ]
    hm_data <- head(hm_data, 30)
    hm_data$Description <- factor(hm_data$Description, levels = rev(unique(hm_data$Description)))

    # Expand to all DX
    hm_expanded <- expand.grid(Description = levels(hm_data$Description), DX = DX_GROUPS)
    hm_expanded <- merge(hm_expanded, enrich_all[, c("Description","DX","p.adjust")],
                         by = c("Description","DX"), all.x = TRUE)

    p_hm <- ggplot(hm_expanded, aes(x = DX, y = Description, fill = -log10(p.adjust))) +
      geom_tile(color="white", linewidth=0.5) +
      scale_fill_gradientn(colors=c("grey90","#fdae61","#d7191c"), na.value="grey90",
                           name="-log10(P)") +
      labs(title="Pathway Enrichment Heatmap: Shared Across DX Groups",
           subtitle="GO BP + Reactome | Grey = not significant", x="", y="") +
      theme_minimal(10) + theme(plot.title=element_text(face="bold"))
    ggsave(file.path(out_dir, "pathway_heatmap_by_DX.png"), p_hm, width=8, height=7, dpi=200)
  }
}

# ---- 4. Summary ----
cat(sprintf("\nOutput: %s\n", normalizePath(out_dir)))
cat("========== Done ==========\n")
