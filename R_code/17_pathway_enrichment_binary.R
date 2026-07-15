###############################################################################
# 17_pathway_enrichment_binary.R
# Data-driven pathway enrichment: CN vs non-CN & AD vs non-AD
###############################################################################

library(dplyr); library(ggplot2); library(clusterProfiler)
library(org.Hs.eg.db); library(enrichplot); library(ReactomePA); library(patchwork)

set.seed(42)
out_dir <- "output/pathway_enrichment_binary"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("========== Pathway Enrichment: Binary Comparisons ==========\n")

# ---- Load ----
df <- read.csv("lifestyle_7var_ptau_proteomics.csv")
df <- df[!is.na(df$ptau217_csf) & !is.na(df$dx_entry), ]
protein_cols <- colnames(df)[21:ncol(df)]
for (pc in protein_cols) df[[pc]] <- suppressWarnings(as.numeric(df[[pc]]))

# Define comparisons
comparisons <- list(
  "CN_vs_nonCN" = list(group1 = "CN", group2 = c("EMCI","LMCI","AD"),
                         labels = c("CN","non-CN")),
  "AD_vs_nonAD" = list(group1 = "AD", group2 = c("CN","EMCI","LMCI"),
                         labels = c("AD","non-AD"))
)

all_enrich <- list()

for (comp_name in names(comparisons)) {
  comp <- comparisons[[comp_name]]
  cat(sprintf("\n========== %s ==========\n", comp_name))

  for (grp_idx in 1:2) {
    if (grp_idx == 1) {
      grp_label <- comp$labels[1]
      df_g <- df[df$dx_entry == comp$group1, ]
    } else {
      grp_label <- comp$labels[2]
      df_g <- df[df$dx_entry %in% comp$group2, ]
    }

    cat(sprintf("\n--- %s (N=%d) ---\n", grp_label, nrow(df_g)))
    if (nrow(df_g) < 50) next

    df_g$PTAU_log2 <- log2(df_g$ptau217_csf)
    mr <- sapply(df_g[protein_cols], function(x) mean(is.na(x))*100)
    pk <- names(mr[mr <= 30])
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
    cat(sprintf("Sig: %d proteins -> %d genes\n", nrow(sig), length(sig_genes)))

    sig_entrez <- tryCatch(bitr(sig_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db), error=function(e) NULL)
    all_entrez <- tryCatch(bitr(all_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db), error=function(e) NULL)
    if (is.null(sig_entrez) || nrow(sig_entrez) < 5) next

    # GO BP
    ego <- enrichGO(gene=sig_entrez$ENTREZID, universe=all_entrez$ENTREZID,
                    OrgDb=org.Hs.eg.db, ont="BP", pvalueCutoff=0.05, qvalueCutoff=0.2, readable=TRUE)
    if (nrow(ego@result) > 0) {
      top <- head(ego@result[order(ego@result$p.adjust), ], 15)
      top$Group <- grp_label; top$Comparison <- comp_name; top$Source <- "GO_BP"
      all_enrich[[paste(comp_name, grp_label, "GO")]] <- top
    }

    # Reactome
    reac <- enrichPathway(gene=sig_entrez$ENTREZID, universe=all_entrez$ENTREZID,
                          organism="human", pvalueCutoff=0.05, qvalueCutoff=0.2, readable=TRUE)
    if (nrow(reac@result) > 0) {
      top <- head(reac@result[order(reac@result$p.adjust), ], 15)
      top$Group <- grp_label; top$Comparison <- comp_name; top$Source <- "Reactome"
      all_enrich[[paste(comp_name, grp_label, "RX")]] <- top
    }

    cat(sprintf("GO: %d | Reactome: %d\n", nrow(ego@result), nrow(reac@result)))
  }
}

# ---- Combine & compare ----
cat("\n========== Generating Comparison Figures ==========\n")

enrich_all <- do.call(rbind, all_enrich)
write.csv(enrich_all, file.path(out_dir, "enrichment_binary_comparison.csv"), row.names = FALSE)

for (comp_name in names(comparisons)) {
  comp_data <- enrich_all[enrich_all$Comparison == comp_name, ]

  for (src in c("GO_BP", "Reactome")) {
    src_data <- comp_data[comp_data$Source == src, ]
    if (nrow(src_data) == 0) next

    # Top 10 per group
    src_plot <- src_data %>%
      group_by(Group) %>%
      slice_min(p.adjust, n = 10) %>%
      ungroup()

    src_plot$label <- substr(src_plot$Description, 1, 55)
    src_plot$label <- factor(src_plot$label, levels = rev(unique(src_plot$label[order(src_plot$p.adjust, decreasing=TRUE)])))

    p <- ggplot(src_plot, aes(x = -log10(p.adjust), y = label, size = Count, color = Group)) +
      geom_point(alpha = 0.85) +
      scale_color_manual(values = c("CN"="#2c7bb6","non-CN"="#d7191c",
                                     "AD"="#d7191c","non-AD"="#2c7bb6")) +
      scale_size_continuous(range = c(2, 7)) +
      facet_wrap(~ Group, scales = "free_y", ncol = 1) +
      labs(title = sprintf("%s Enrichment: %s", src, comp_name),
           subtitle = "Top 10 pathways per group | ORA of pTau-significant proteins",
           x = "-log10(adjusted P)", y = "") +
      theme_bw(11) + theme(plot.title = element_text(face = "bold"))

    ggsave(file.path(out_dir, sprintf("%s_%s.png", comp_name, src)),
           p, width = 11, height = 8, dpi = 200)
    cat(sprintf("Saved: %s_%s.png\n", comp_name, src))
  }

  # Shared vs unique top pathways
  groups <- unique(comp_data$Group)
  if (length(groups) == 2) {
    terms1 <- comp_data$Description[comp_data$Group == groups[1]]
    terms2 <- comp_data$Description[comp_data$Group == groups[2]]
    shared <- intersect(terms1, terms2)
    only1  <- setdiff(terms1, terms2)
    only2  <- setdiff(terms2, terms1)
    cat(sprintf("\n%s: shared=%d, %s-only=%d, %s-only=%d\n",
                comp_name, length(shared), groups[1], length(only1), groups[2], length(only2)))
  }
}

cat(sprintf("\nOutput: %s\n", normalizePath(out_dir)))
cat("========== Done ==========\n")
