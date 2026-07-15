###############################################################################
# 10_tau_mechanism_analysis.R
# From PPI Network → Tau Regulation Mechanism Model
# Steps:
#  1. Community detection (Louvain) on pTau-significant protein network
#  2. Pathway enrichment per module (Reactome + GO + KEGG)
#  3. Functional protein classification (kinases, phosphatases, etc.)
#  4. Tau regulation map construction
#  5. Statistical enrichment of tau-related gene sets
#  6. Biological interpretation & candidate prioritization
###############################################################################

library(readxl)
library(dplyr)
library(igraph)
library(ggraph)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(patchwork)

set.seed(42)

out_dir <- "output/tau_mechanism_analysis"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("========== From Network to Mechanism: Tau Regulation Model ==========\n")

# ===========================================================================
# 0. Load data & prepare
# ===========================================================================
cat("\n--- 0. Loading data ---\n")

master <- read.csv("lifestyle_7var_ptau_proteomics.csv")

df_bl <- master[!is.na(master$ptau217_csf) & !is.na(master$dx_entry), ]
protein_cols <- colnames(df_bl)[21:ncol(df_bl)]

convert_protein <- function(x) { x[x == "NA" | x == ""] <- NA; as.numeric(x) }
df_bl[protein_cols] <- lapply(df_bl[protein_cols], convert_protein)

# Run correlation
df_bl$PTAU_log2 <- log2(df_bl$ptau217_csf)
ptau_vals <- df_bl$PTAU_log2

missing_rate <- sapply(df_bl[protein_cols], function(x) mean(is.na(x)) * 100)
proteins_kept <- names(missing_rate[missing_rate <= 20])
protein_log2 <- log2(as.matrix(df_bl[proteins_kept]))
protein_log2[is.infinite(protein_log2) | is.nan(protein_log2)] <- NA

results <- data.frame(protein_id = colnames(protein_log2), n = NA_integer_,
                      spearman_rho = NA_real_, p_value = NA_real_)
for (i in seq_len(ncol(protein_log2))) {
  pv <- protein_log2[, i]; vi <- !is.na(pv); nv <- sum(vi)
  results$n[i] <- nv
  if (nv >= 10) {
    tst <- tryCatch(cor.test(pv[vi], ptau_vals[vi], method = "spearman", exact = FALSE), error = function(e) NULL)
    if (!is.null(tst)) { results$spearman_rho[i] <- tst$estimate; results$p_value[i] <- tst$p.value }
  }
}
results <- results[!is.na(results$p_value), ]
results$fdr <- p.adjust(results$p_value, method = "BH")
results$EntrezGeneSymbol <- gsub("\\.[0-9]+\\.[0-9]+$", "", results$protein_id)

sig <- results[results$fdr < 0.05, ]
sig_genes_all <- unique(sig$EntrezGeneSymbol[!is.na(sig$EntrezGeneSymbol) & sig$EntrezGeneSymbol != ""])
cat(sprintf("FDR<0.05: %d proteins → %d genes\n", nrow(sig), length(sig_genes_all)))

# Gene-level summary
gene_rho <- sig %>%
  filter(!is.na(EntrezGeneSymbol) & EntrezGeneSymbol != "") %>%
  group_by(EntrezGeneSymbol) %>%
  summarise(
    rho_median     = median(spearman_rho, na.rm = TRUE),
    rho_mean       = mean(spearman_rho, na.rm = TRUE),
    n_proteins     = n(),
    best_protein   = protein_id[which.max(abs(spearman_rho))],
    .groups = "drop"
  )

# Download OmniPath
cat("Downloading OmniPath...\n")
omni_url <- "https://omnipathdb.org/interactions?format=tsv&fields=sources,references&genesymbols=1"
tmp <- tempfile(fileext = ".tsv")
download.file(omni_url, tmp, method = "auto")
omni <- read.delim(tmp, stringsAsFactors = FALSE)

# Sig-sig interactions
inter_sig <- omni[
  omni$source_genesymbol %in% sig_genes_all &
  omni$target_genesymbol %in% sig_genes_all,
]
cat(sprintf("Sig-sig interactions: %d\n", nrow(inter_sig)))

# Build igraph
edges <- inter_sig[, c("source_genesymbol", "target_genesymbol")]
colnames(edges) <- c("from", "to")
g <- graph_from_data_frame(edges, directed = FALSE)

# Add node attributes
V(g)$rho        <- gene_rho$rho_median[match(V(g)$name, gene_rho$EntrezGeneSymbol)]
V(g)$rho_abs    <- abs(V(g)$rho)
V(g)$degree     <- degree(g)
V(g)$betweenness <- betweenness(g, normalized = TRUE)

cat(sprintf("Network: %d nodes, %d edges\n", vcount(g), ecount(g)))

# ===========================================================================
# Step 1: Louvain community detection
# ===========================================================================
cat("\n========== Step 1: Community Detection ==========\n")

set.seed(42)
communities <- cluster_louvain(g)

V(g)$module <- communities$membership
n_modules <- max(communities$membership)
cat(sprintf("Modules detected: %d\n", n_modules))

# Module statistics
module_stats <- data.frame()
for (m in 1:n_modules) {
  nodes_m <- V(g)$name[V(g)$module == m]
  n_m <- length(nodes_m)
  subg <- induced_subgraph(g, nodes_m)
  avg_deg <- if (n_m > 1) mean(degree(subg)) else 0
  top_hubs <- head(sort(degree(subg), decreasing = TRUE), 3)
  hub_names <- paste(names(top_hubs), collapse = ", ")
  avg_rho <- mean(V(g)$rho[V(g)$module == m], na.rm = TRUE)

  module_stats <- rbind(module_stats, data.frame(
    Module = m, Size = n_m, Avg_Degree = round(avg_deg, 1),
    Avg_Rho = round(avg_rho, 4),
    Top_Hubs = hub_names,
    stringsAsFactors = FALSE
  ))
}
module_stats <- module_stats[order(module_stats$Size, decreasing = TRUE), ]
cat("\nModule summary:\n")
print(module_stats)

# Visualize: community network
# Focus on modules with >=10 members for clarity
main_modules <- module_stats$Module[module_stats$Size >= 10]
g_main <- induced_subgraph(g, V(g)$name[V(g)$module %in% main_modules])

# For large networks, show top-connected nodes
if (vcount(g_main) > 200) {
  top_nodes <- V(g_main)$name[order(V(g_main)$degree, decreasing = TRUE)][1:200]
  g_main <- induced_subgraph(g_main, top_nodes)
}

set.seed(123)
# Community colors
n_mod_colors <- length(unique(V(g_main)$module))
mod_colors <- hcl.colors(n_mod_colors, "Dynamic")
names(mod_colors) <- sort(unique(V(g_main)$module))

p_comm <- ggraph(g_main, layout = "fr") +
  geom_edge_link(color = "grey90", alpha = 0.15) +
  geom_node_point(aes(size = degree, fill = factor(module)), shape = 21,
                  color = "grey30", stroke = 0.2, alpha = 0.85) +
  geom_node_text(aes(
    label = ifelse(degree > quantile(degree, 0.95) | abs(rho) > 0.7, name, ""),
    size = degree), repel = TRUE, max.overlaps = 50, box.padding = 0.3) +
  scale_fill_manual(values = mod_colors, name = "Module") +
  scale_size_continuous(range = c(1, 6), guide = "none") +
  labs(
    title = "pTau-Associated Protein Network: Louvain Communities",
    subtitle = sprintf("%d proteins, %d modules | Colors = functional communities | Labeled = top hubs",
                       vcount(g_main), length(unique(V(g_main)$module))),
    caption = "OmniPath PPI | Node size = degree"
  ) +
  theme_void() +
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 9, color = "grey40"))

ggsave(file.path(out_dir, "community_network.png"),
       p_comm, width = 16, height = 14, dpi = 200)
cat("Saved: community_network.png\n")

# ===========================================================================
# Step 2: Pathway enrichment per module
# ===========================================================================
cat("\n========== Step 2: Pathway Enrichment ==========\n")

# Tau-related pathway filter list
tau_pathway_keywords <- c(
  "tau", "MAPK", "PI3K", "AKT", "GSK", "CDK5", "phosphorylation",
  "dephosphorylation", "synaptic", "synapse", "neuroinflammation",
  "microglia", "cytoskeleton", "axon", "Alzheimer", "neurodegener",
  "amyloid", "APP", "calcium", "CAMK", "PKC", "SRC", "FYN",
  "phosphatase", "kinase", "ubiquitin", "autophagy", "apoptosis",
  "Wnt", "Notch", "mTOR"
)

all_enrich <- list()
module_enrich_summary <- data.frame()

for (m in sort(unique(module_stats$Module))) {
  nodes_m <- V(g)$name[V(g)$module == m]
  if (length(nodes_m) < 5) next
  cat(sprintf("\n--- Module %d (%d genes) ---\n", m, length(nodes_m)))

  # Convert to EntrezID for clusterProfiler
  entrez_ids <- tryCatch(
    bitr(nodes_m, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db),
    error = function(e) NULL
  )
  if (is.null(entrez_ids) || nrow(entrez_ids) < 3) {
    cat("  Too few genes with Entrez mapping, skipping.\n")
    next
  }

  # GO BP enrichment
  ego <- tryCatch(
    enrichGO(gene = entrez_ids$ENTREZID, OrgDb = org.Hs.eg.db,
             ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
             readable = TRUE),
    error = function(e) NULL
  )

  # Reactome enrichment
  reac <- tryCatch(
    enrichPathway(gene = entrez_ids$ENTREZID, organism = "human",
                  pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE),
    error = function(e) NULL
  )

  # Collect top enriched terms
  for (src in c("GO_BP", "Reactome")) {
    enrich_obj <- if (src == "GO_BP") ego else reac
    if (is.null(enrich_obj) || nrow(enrich_obj@result) == 0) next

    top <- head(enrich_obj@result, 10)
    top$Module <- m
    top$Source <- src
    top$Module_Size <- length(nodes_m)
    all_enrich[[paste(m, src, sep = "_")]] <- top

    # Mark tau-relevant terms
    top$tau_related <- grepl(paste(tau_pathway_keywords, collapse = "|"),
                             top$Description, ignore.case = TRUE)

    module_enrich_summary <- rbind(module_enrich_summary,
      top[top$tau_related, c("Module", "Source", "Description", "p.adjust", "Count")])
  }
}

# Combine enrichment results
if (length(all_enrich) > 0) {
  enrich_all <- do.call(rbind, all_enrich)
  enrich_all$tau_related <- grepl(paste(tau_pathway_keywords, collapse = "|"),
                                  enrich_all$Description, ignore.case = TRUE)

  # Dotplot: top tau-related pathways per module
  enrich_tau <- enrich_all[enrich_all$tau_related, ]
  if (nrow(enrich_tau) > 0) {
    # Select top 5 per module
    enrich_plot <- enrich_tau %>%
      group_by(Module) %>%
      slice_min(p.adjust, n = 5) %>%
      ungroup()

    enrich_plot$label <- sprintf("M%d: %s", enrich_plot$Module,
                                  substr(enrich_plot$Description, 1, 60))
    enrich_plot <- enrich_plot[order(enrich_plot$p.adjust, decreasing = TRUE), ]
    enrich_plot$label <- factor(enrich_plot$label, levels = enrich_plot$label)

    p_enrich <- ggplot(enrich_plot, aes(x = -log10(p.adjust), y = label,
                                         size = Count, color = Source)) +
      geom_point(alpha = 0.85) +
      scale_color_manual(values = c("GO_BP" = "#2c7bb6", "Reactome" = "#d7191c")) +
      scale_size_continuous(range = c(2, 7)) +
      facet_wrap(~ factor(Module), scales = "free_y", ncol = 1) +
      labs(
        title = "Tau-Relevant Pathway Enrichment by Network Module",
        subtitle = "GO Biological Process + Reactome | Top 5 terms per module",
        x = "-log10(adjusted P)", y = ""
      ) +
      theme_bw(base_size = 11) +
      theme(plot.title = element_text(face = "bold"))

    ggsave(file.path(out_dir, "module_pathway_enrichment.png"),
           p_enrich, width = 12, height = max(8, nrow(enrich_plot) * 0.3), dpi = 200)
    cat("Saved: module_pathway_enrichment.png\n")
  }

  write.csv(enrich_all, file.path(out_dir, "module_enrichment_all.csv"), row.names = FALSE)
}

# ===========================================================================
# Step 3: Functional protein classification
# ===========================================================================
cat("\n========== Step 3: Functional Classification ==========\n")

# Curated gene lists for tau-relevant functional categories
tau_categories <- list(
  "Tau Kinases" = c("GSK3B", "GSK3A", "CDK5", "MAPT", "MAPK1", "MAPK3", "MAPK8",
                    "MAPK9", "MAPK10", "MAPK11", "MAPK12", "MAPK13", "MAPK14",
                    "DYRK1A", "CSNK1D", "CSNK1E", "CSNK1G1", "CSNK2A1", "CSNK2A2",
                    "TTBK1", "TTBK2", "MARK1", "MARK2", "MARK3", "MARK4",
                    "PRKAA1", "PRKAA2", "NUAK1", "NUAK2", "PHKG1", "PHKG2",
                    "PRKCA", "PRKCB", "PRKCG", "PRKCD", "PRKCE", "PRKCZ", "PRKCI",
                    "PRKACA", "PRKACB", "PRKG1", "CAMK2A", "CAMK2B", "CAMK2D",
                    "PINK1", "LRRK2", "SRPK1", "SRPK2", "EIF2AK2", "PKN1",
                    "SGK1", "TTK", "CHEK2", "RPS6KA3", "RPS6KB1", "PDPK1"),
  "Tau Phosphatases" = c("PPP1CA", "PPP1CB", "PPP1CC", "PPP2CA", "PPP2CB",
                         "PPP2R1A", "PPP2R5A", "PPP2R5B", "PPP2R5C", "PPP2R5D",
                         "PPP2R5E", "PPP3CA", "PPP3CB", "PPP3CC", "PPP3R1",
                         "PPP5C", "PTEN", "DUSP1", "DUSP6", "CDKN3"),
  "Tau Folding/Stability" = c("PIN1", "HSP90AA1", "HSP90AB1", "HSPA1A", "HSPA8",
                              "BAG2", "BAG3", "CHIP", "STUB1", "FKBP5", "FKBP4",
                              "PPID", "PPIA", "PPIB", "PPIC", "PPIF",
                              "DNAJA1", "DNAJA2", "DNAJB1", "DNAJB6",
                              "HSPB1", "HSPB8", "CRYAB"),
  "Synaptic Proteins" = c("DLG4", "DLG2", "DLG1", "DLG3", "SYN1", "SYN2", "SYN3",
                          "SYP", "SYNGR1", "SYNGR3", "SNAP25", "SNAP91",
                          "STX1A", "STXBP1", "VAMP2", "SYT1", "SYT4", "SYT5",
                          "GRIN1", "GRIN2A", "GRIN2B", "GRIA1", "GRIA2", "GRIA3", "GRIA4",
                          "GABRA1", "GABRB2", "GABRG2", "SHANK2", "SHANK3",
                          "HOMER1", "HOMER2", "HOMER3", "ARC", "BDNF", "NGF",
                          "NTRK2", "NGFR", "GAP43", "NRGN", "STMN1", "STMN2",
                          "CAMK2A", "CAMK2B", "PRKCG", "YWHAH", "YWHAE", "YWHAG",
                          "YWHAB", "SNCA", "SNCB", "PSD95", "DLGAP1"),
  "Cytoskeleton / Axonal Transport" = c("MAP2", "MAP1B", "MAP1A", "TUBB", "TUBB2A",
                                        "TUBB3", "TUBB4A", "TUBA1A", "TUBA1B",
                                        "TUBA4A", "ACTB", "ACTG1", "ACTN1", "ACTN2",
                                        "ACTN4", "SPTAN1", "SPTBN1", "SPTBN2",
                                        "DCTN1", "DCTN2", "DYNC1H1", "DYNC1I2",
                                        "DYNC1LI1", "DYNC1LI2", "DYNLL1", "DYNLL2",
                                        "KIF1A", "KIF1B", "KIF2A", "KIF3A", "KIF5A",
                                        "KIF5B", "KIF5C", "KLC1", "KLC2", "KLC3",
                                        "KLC4", "DNM1", "DNM2", "DNM3", "CLTC",
                                        "CLTA", "AP2B1", "AP2A1", "AP2A2", "AP2M1"),
  "Microglia / Immune" = c("CD33", "TREM2", "TYROBP", "CSF1R", "CX3CR1",
                           "ITGAM", "ITGB2", "CD68", "CD14", "TLR2", "TLR4",
                           "MYD88", "NFKB1", "NFKB2", "RELA", "RELB", "NFKBIA",
                           "TNF", "IL1B", "IL6", "IL10", "TGFB1", "TGFB2",
                           "C1QA", "C1QB", "C1QC", "C3", "C4A", "C4B",
                           "SRC", "SYK", "LYN", "HCK", "FGR", "BTK"),
  "Astrocyte Activation" = c("GFAP", "S100B", "ALDH1L1", "AQP4", "GJA1", "GJB6",
                             "VIM", "NES", "SOX9", "STAT3", "NOTCH1", "NOTCH2",
                             "HES1", "HES5", "NFIA", "NFIB"),
  "Growth Factor Signaling" = c("AKT1", "AKT2", "AKT3", "MTOR", "RPS6KB1",
                                "RPS6KA1", "RPS6KA2", "RPS6KA3", "RPS6KA4",
                                "PIK3CA", "PIK3CB", "PIK3R1", "PIK3R2",
                                "IGF1R", "INSR", "IRS1", "IRS2", "EGFR",
                                "ERBB2", "ERBB3", "ERBB4", "FGFR1", "FGFR2",
                                "FGFR3", "NTRK1", "NTRK2", "NTRK3",
                                "PDGFRA", "PDGFRB", "VEGFA", "FLT1", "KDR",
                                "GRB2", "SOS1", "SOS2", "HRAS", "KRAS", "NRAS",
                                "BRAF", "RAF1", "MAP2K1", "MAP2K2"),
  "Metabolism" = c("HK1", "HK2", "PKM", "PFKM", "PFKL", "PFKP", "LDHA", "LDHB",
                   "PDHA1", "PDHB", "CS", "IDH1", "IDH2", "IDH3A", "OGDH",
                   "SDHA", "SDHB", "FH", "MDH1", "MDH2", "ACO1", "ACO2",
                   "MT-ND1", "MT-CO1", "MT-CO2", "MT-ATP6", "SLC2A1", "SLC2A3",
                   "SLC2A4", "GAPDH", "PGK1", "PGAM1", "ENO1", "ENO2",
                   "ALDOA", "ALDOC", "TPI1", "APOE", "CLU", "LPL", "LDLR",
                   "ABCA1", "ABCG1", "SOAT1", "CYP46A1")
)

# Classify each significant gene
sig_gene_info <- gene_rho %>%
  mutate(
    Category1 = NA_character_,
    Category2 = NA_character_
  )

for (i in seq_len(nrow(sig_gene_info))) {
  gn <- sig_gene_info$EntrezGeneSymbol[i]
  cats <- names(tau_categories)[sapply(tau_categories, function(gs) gn %in% gs)]
  if (length(cats) >= 1) {
    sig_gene_info$Category1[i] <- cats[1]
    if (length(cats) >= 2) sig_gene_info$Category2[i] <- cats[2]
  }
}
sig_gene_info$Category1[is.na(sig_gene_info$Category1)] <- "Others"

# Category summary
cat_summary <- sig_gene_info %>%
  group_by(Category1) %>%
  summarise(
    N = n(),
    Avg_Rho = round(mean(rho_median, na.rm = TRUE), 4),
    Top_Genes = paste(head(EntrezGeneSymbol[order(abs(rho_median), decreasing = TRUE)], 5), collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(desc(N))

cat("\nFunctional category summary:\n")
print(cat_summary)

write.csv(cat_summary, file.path(out_dir, "functional_category_summary.csv"), row.names = FALSE)
write.csv(sig_gene_info, file.path(out_dir, "gene_functional_classification.csv"), row.names = FALSE)

# Category barplot
cat_plot <- cat_summary[cat_summary$Category1 != "Others", ]
cat_plot$Category1 <- factor(cat_plot$Category1, levels = cat_plot$Category1[order(cat_plot$N)])

p_cat <- ggplot(cat_plot, aes(x = N, y = Category1, fill = Avg_Rho)) +
  geom_col(alpha = 0.85) +
  geom_text(aes(label = N), hjust = -0.2, size = 3.5) +
  scale_fill_gradient2(low = "#2c7bb6", mid = "white", high = "#d7191c",
                        midpoint = 0, name = "Mean PTAU rho") +
  labs(
    title = "Functional Classification of pTau-Associated Proteins",
    subtitle = sprintf("%d genes classified into tau-relevant functional categories",
                       sum(cat_summary$N)),
    x = "Number of Genes", y = ""
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(out_dir, "functional_categories.png"),
       p_cat, width = 10, height = 5, dpi = 200)
cat("Saved: functional_categories.png\n")

# ===========================================================================
# Step 4: Tau regulation map — MAPT-centered signaling
# ===========================================================================
cat("\n========== Step 4: Tau Regulation Signaling Map ==========\n")

# MAPT interactors in our sig set
mapt_edges <- omni[
  (omni$source_genesymbol == "MAPT" | omni$target_genesymbol == "MAPT") &
  (omni$source_genesymbol %in% sig_genes_all | omni$target_genesymbol %in% sig_genes_all),
]

mapt_genes <- unique(c(mapt_edges$source_genesymbol, mapt_edges$target_genesymbol))
mapt_genes <- setdiff(mapt_genes, "MAPT")
cat(sprintf("MAPT interactors in sig set: %d\n", length(mapt_genes)))

# Extract key tau-regulating genes from multiple axes
tau_axes <- list(
  "GSK3B axis"    = c("MAPT", "GSK3B", "GSK3A", "AKT1", "AKT2", "PTEN", "PDPK1",
                      "PIK3CA", "PIK3R1", "PRKCA", "PRKCB", "PRKCG", "PRKCZ",
                      "PPP1CA", "PPP2CA", "PPP2R1A", "DYRK1A", "WNT3A", "FZD1",
                      "DVL1", "AXIN1", "APC", "CTNNB1"),
  "MAPK cascade"   = c("MAPT", "MAPK1", "MAPK3", "MAPK8", "MAPK9", "MAPK14",
                      "MAP2K1", "MAP2K2", "MAP2K4", "BRAF", "RAF1", "RASGRF1",
                      "GRB2", "SOS1", "HRAS", "KRAS", "RPS6KA3", "RPS6KB1"),
  "CDK5 pathway"   = c("MAPT", "CDK5", "CDK5R1", "CDK5R2", "PPP1CA", "PPP1CB",
                      "PPP2CA", "PPP3CA", "PPP3R1", "PPP5C", "CAMK2A", "PRKCA",
                      "ABL1", "ABL2"),
  "SRC/FYN/SYK"    = c("MAPT", "SRC", "FYN", "SYK", "LYN", "HCK", "LCK",
                      "PTK2", "PTK2B", "PTPN11", "GRB2", "PIK3R1"),
  "PPP Phosphatases" = c("MAPT", "PPP1CA", "PPP1CB", "PPP1CC", "PPP2CA",
                         "PPP2CB", "PPP2R1A", "PPP2R5A", "PPP3CA", "PPP3CB",
                         "PPP3CC", "PPP3R1", "PPP5C", "PTEN")
)

# Find which axis genes are in our sig set
axis_summary <- data.frame()
for (ax_name in names(tau_axes)) {
  ax_genes <- tau_axes[[ax_name]]
  ax_sig <- intersect(ax_genes, sig_genes_all)
  axis_summary <- rbind(axis_summary, data.frame(
    Axis = ax_name, Total = length(ax_genes), In_Sig = length(ax_sig),
    Coverage = round(100 * length(ax_sig) / length(ax_genes), 1),
    Top_Hits = paste(head(ax_sig[order(abs(gene_rho$rho_median[match(ax_sig, gene_rho$EntrezGeneSymbol)]), decreasing = TRUE)], 5), collapse = ", "),
    stringsAsFactors = FALSE
  ))
}

cat("\nTau signaling axis coverage:\n")
print(axis_summary)

# Build simplified tau signaling map — MAPT + 1st-degree interactors from key axes
tau_map_genes <- unique(c("MAPT", unlist(lapply(tau_axes[1:3], function(x) intersect(x, sig_genes_all)))))
tau_map_genes <- tau_map_genes[tau_map_genes %in% names(V(g))]

if (length(tau_map_genes) > 2) {
  g_tau <- induced_subgraph(g, intersect(tau_map_genes, V(g)$name))

  # Color by axis membership
  V(g_tau)$axis <- "Other"
  for (ax_name in names(tau_axes)) {
    V(g_tau)$axis[V(g_tau)$name %in% tau_axes[[ax_name]]] <- ax_name
  }
  V(g_tau)$axis[V(g_tau)$name == "MAPT"] <- "MAPT"

  axis_colors <- c("MAPT" = "black", "GSK3B axis" = "#d7191c",
                   "MAPK cascade" = "#fdae61", "CDK5 pathway" = "#2c7bb6",
                   "SRC/FYN/SYK" = "#1b9e77", "PPP Phosphatases" = "#984ea3",
                   "Other" = "grey70")

  set.seed(123)
  p_tau_map <- ggraph(g_tau, layout = "stress") +
    geom_edge_link(color = "grey80", alpha = 0.5) +
    geom_node_point(aes(size = degree, fill = axis), shape = 21) +
    geom_node_text(aes(label = name, size = ifelse(name == "MAPT", 4.5, 2.8)),
                   repel = TRUE, max.overlaps = 50, box.padding = 0.3) +
    scale_fill_manual(values = axis_colors, name = "Signaling Axis") +
    scale_size_continuous(range = c(2, 8), guide = "none") +
    labs(
      title = "Tau Regulation Signaling Map",
      subtitle = sprintf("MAPT + 1st-degree interactors from GSK3B, MAPK, CDK5 axes | %d nodes",
                         vcount(g_tau)),
      caption = "OmniPath PPI | Colors = signaling module affiliation"
    ) +
    theme_void() +
    theme(plot.title = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(size = 9, color = "grey40"))

  ggsave(file.path(out_dir, "tau_signaling_map.png"),
         p_tau_map, width = 14, height = 12, dpi = 200)
  cat("Saved: tau_signaling_map.png\n")
}

write.csv(axis_summary, file.path(out_dir, "tau_signaling_axis_coverage.csv"), row.names = FALSE)

# ===========================================================================
# Step 5: Statistical enrichment of tau-related gene sets
# ===========================================================================
cat("\n========== Step 5: Tau Gene Set Enrichment ==========\n")

# Prepare ranked gene list for GSEA
ranked_genes <- gene_rho$rho_median
names(ranked_genes) <- gene_rho$EntrezGeneSymbol
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

# Convert to Entrez
gene_list_entrez <- tryCatch({
  bitr(names(ranked_genes), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
}, error = function(e) NULL)

if (!is.null(gene_list_entrez)) {
  ranked_entrez <- ranked_genes[gene_list_entrez$SYMBOL]
  names(ranked_entrez) <- gene_list_entrez$ENTREZID[match(names(ranked_entrez), gene_list_entrez$SYMBOL)]
  ranked_entrez <- sort(ranked_entrez, decreasing = TRUE)

  # Fisher test: are tau-related pathways over-represented?
  # Use Reactome pathways
  tau_reactome <- tryCatch(
    enrichPathway(gene = gene_list_entrez$ENTREZID[names(ranked_entrez) %in%
                   names(ranked_entrez[abs(ranked_entrez) > 0.3])],
                  organism = "human", pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                  readable = TRUE),
    error = function(e) NULL
  )

  if (!is.null(tau_reactome) && nrow(tau_reactome@result) > 0) {
    # Dotplot of tau-related Reactome terms
    tau_rx <- tau_reactome@result
    tau_rx$tau_related <- grepl(paste(tau_pathway_keywords, collapse = "|"),
                                tau_rx$Description, ignore.case = TRUE)
    tau_rx_sig <- tau_rx[tau_rx$tau_related & tau_rx$p.adjust < 0.05, ]

    if (nrow(tau_rx_sig) > 0) {
      tau_rx_sig <- tau_rx_sig[order(tau_rx_sig$p.adjust), ]
      tau_rx_plot <- head(tau_rx_sig, 20)
      tau_rx_plot <- tau_rx_plot[order(tau_rx_plot$p.adjust, decreasing = TRUE), ]
      tau_rx_plot$label <- factor(
        substr(tau_rx_plot$Description, 1, 70),
        levels = substr(tau_rx_plot$Description, 1, 70)
      )

      p_rx <- ggplot(tau_rx_plot, aes(x = -log10(p.adjust), y = label,
                                       size = Count, fill = -log10(p.adjust))) +
        geom_point(shape = 21, alpha = 0.85) +
        scale_fill_gradientn(colors = c("#2c7bb6", "#fdae61", "#d7191c"),
                              name = "-log10(P)") +
        scale_size_continuous(range = c(3, 8)) +
        labs(
          title = "Tau-Relevant Reactome Pathways Enriched in pTau-Associated Proteome",
          subtitle = sprintf("%d pathways | FDR < 0.05", nrow(tau_rx_sig)),
          x = "-log10(adjusted P)", y = ""
        ) +
        theme_bw(base_size = 11) +
        theme(plot.title = element_text(face = "bold"))

      ggsave(file.path(out_dir, "tau_reactome_enrichment.png"),
             p_rx, width = 12, height = 6, dpi = 200)
      cat("Saved: tau_reactome_enrichment.png\n")
    }
  }
}

# ===========================================================================
# Step 6: Candidate prioritization table
# ===========================================================================
cat("\n========== Step 6: Candidate Prioritization ==========\n")

# Score each gene by: |rho| × degree × (is_MAPT_interactor) × (is_in_tau_pathway)
sig_gene_info$network_degree <- V(g)$degree[match(sig_gene_info$EntrezGeneSymbol, V(g)$name)]
sig_gene_info$network_degree[is.na(sig_gene_info$network_degree)] <- 0

sig_gene_info$is_MAPT_interactor <- sig_gene_info$EntrezGeneSymbol %in% mapt_genes
sig_gene_info$in_tau_pathway <- sig_gene_info$Category1 != "Others"

# Composite priority score
sig_gene_info$priority_score <- with(sig_gene_info,
  abs(rho_median) * log1p(network_degree) *
  (1 + is_MAPT_interactor * 0.5) * (1 + in_tau_pathway * 0.3)
)
sig_gene_info <- sig_gene_info[order(sig_gene_info$priority_score, decreasing = TRUE), ]

candidates <- sig_gene_info[, c("EntrezGeneSymbol", "rho_median", "rho_mean",
                                 "network_degree", "Category1", "is_MAPT_interactor",
                                 "n_proteins", "priority_score")]
colnames(candidates) <- c("Gene", "PTAU_rho", "Mean_Rho", "Degree",
                           "Functional_Category", "MAPT_Interactor",
                           "N_Proteins", "Priority_Score")
candidates$Rank <- 1:nrow(candidates)
candidates <- candidates[, c("Rank", setdiff(colnames(candidates), "Rank"))]

cat("\nTop 20 priority candidates for mechanistic validation:\n")
print(head(candidates, 20))

write.csv(candidates, file.path(out_dir, "candidate_priority_ranking.csv"), row.names = FALSE)

# ===========================================================================
# 7. Final biological interpretation
# ===========================================================================
cat("\n")
cat("========== Biological Interpretation ==========\n\n")

cat("1. Dominant signaling pathways in the pTau-associated proteome:\n\n")

# Get top pathway categories
top_pathways <- module_enrich_summary %>%
  group_by(Source) %>%
  slice_min(p.adjust, n = 5) %>%
  ungroup()

cat("   Key enriched pathways:\n")
for (i in 1:min(nrow(top_pathways), 10)) {
  cat(sprintf("   - %s (%s, P=%.1e)\n", top_pathways$Description[i],
              top_pathways$Source[i], top_pathways$p.adjust[i]))
}

cat("\n2. Tau signaling axis coverage:\n\n")
for (i in 1:nrow(axis_summary)) {
  cat(sprintf("   %s: %d/%d genes (%.0f%% coverage)\n",
              axis_summary$Axis[i], axis_summary$In_Sig[i],
              axis_summary$Total[i], axis_summary$Coverage[i]))
}

# Identify multi-axis integrators
multi_axis_genes <- sig_gene_info$EntrezGeneSymbol[
  sig_gene_info$is_MAPT_interactor & sig_gene_info$in_tau_pathway
]
cat(sprintf("\n3. Multi-axis integrators (MAPT-interacting + in tau pathway): %d genes\n",
            length(multi_axis_genes)))
cat(sprintf("   Top: %s\n", paste(head(candidates$Gene[candidates$MAPT_Interactor & candidates$Functional_Category != "Others"], 10), collapse = ", ")))

cat(sprintf("\n4. Top candidates for mechanistic validation:\n   %s\n",
            paste(head(candidates$Gene, 10), collapse = ", ")))

cat(sprintf("\nOutput directory: %s\n", normalizePath(out_dir)))
cat("\n========== Done ==========\n")
