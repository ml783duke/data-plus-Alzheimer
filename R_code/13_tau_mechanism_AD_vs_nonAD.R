###############################################################################
# 12_tau_mechanism_CN_vs_nonCN.R
# Tau Mechanism Analysis: AD vs non-AD (EMCI+LMCI+AD)
# - Run full pathway + network analysis per group
# - Compare communities, enrichment, functional categories, candidates
###############################################################################

library(readxl); library(dplyr); library(igraph); library(ggraph)
library(ggplot2); library(clusterProfiler); library(org.Hs.eg.db)
library(enrichplot); library(patchwork)

set.seed(42)
out_dir <- "output/tau_mechanism_AD_vs_nonAD"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("========== Tau Mechanism: AD vs non-AD ==========\n")

# ---- 0. Load data ----
cat("\n--- Loading ---\n")
master <- read.csv("lifestyle_7var_ptau_proteomics.csv")
df_bl <- master[!is.na(master$ptau217_csf) & !is.na(df_bl$dx_entry), ]
protein_cols <- colnames(df_bl)[21:ncol(df_bl)]
convert_protein <- function(x) { x[x == "NA" | x == ""] <- NA; as.numeric(x) }
df_bl[protein_cols] <- lapply(df_bl[protein_cols], convert_protein)

# Download OmniPath once
omni_url <- "https://omnipathdb.org/interactions?format=tsv&fields=sources,references&genesymbols=1"
tmp <- tempfile(fileext = ".tsv"); download.file(omni_url, tmp, method = "auto")
omni <- read.delim(tmp, stringsAsFactors = FALSE)

# --- Tau categories (same as 10) ---
tau_categories <- list(
  "Tau Kinases" = c("GSK3B","GSK3A","CDK5","MAPT","MAPK1","MAPK3","MAPK8","MAPK9","MAPK10","MAPK11","MAPK12","MAPK13","MAPK14","DYRK1A","CSNK1D","CSNK1E","CSNK1G1","CSNK2A1","CSNK2A2","TTBK1","TTBK2","MARK1","MARK2","MARK3","MARK4","PRKAA1","PRKAA2","NUAK1","NUAK2","PHKG1","PHKG2","PRKCA","PRKCB","PRKCG","PRKCD","PRKCE","PRKCZ","PRKCI","PRKACA","PRKACB","PRKG1","CAMK2A","CAMK2B","CAMK2D","PINK1","LRRK2","SRPK1","SRPK2","EIF2AK2","PKN1","SGK1","TTK","CHEK2","RPS6KA3","RPS6KB1","PDPK1"),
  "Tau Phosphatases" = c("PPP1CA","PPP1CB","PPP1CC","PPP2CA","PPP2CB","PPP2R1A","PPP2R5A","PPP2R5B","PPP2R5C","PPP2R5D","PPP2R5E","PPP3CA","PPP3CB","PPP3CC","PPP3R1","PPP5C","PTEN","DUSP1","DUSP6","CDKN3"),
  "Tau Folding/Stability" = c("PIN1","HSP90AA1","HSP90AB1","HSPA1A","HSPA8","BAG2","BAG3","STUB1","FKBP5","FKBP4","PPID","PPIA","PPIB","PPIC","PPIF","DNAJA1","DNAJA2","DNAJB1","DNAJB6","HSPB1","HSPB8","CRYAB"),
  "Synaptic Proteins" = c("DLG4","DLG2","DLG1","DLG3","SYN1","SYN2","SYN3","SYP","SYNGR1","SYNGR3","SNAP25","SNAP91","STX1A","STXBP1","VAMP2","SYT1","SYT4","SYT5","GRIN1","GRIN2A","GRIN2B","GRIA1","GRIA2","GRIA3","GRIA4","GABRA1","GABRB2","GABRG2","SHANK2","SHANK3","HOMER1","HOMER2","HOMER3","ARC","BDNF","NGF","NTRK2","NGFR","GAP43","NRGN","STMN1","STMN2","CAMK2A","CAMK2B","PRKCG","YWHAH","YWHAE","YWHAG","YWHAB","SNCA","SNCB","DLGAP1"),
  "Cytoskeleton/Axonal" = c("MAP2","MAP1B","MAP1A","TUBB","TUBB2A","TUBB3","TUBB4A","TUBA1A","TUBA1B","TUBA4A","ACTB","ACTG1","ACTN1","ACTN2","ACTN4","SPTAN1","SPTBN1","SPTBN2","DCTN1","DCTN2","DYNC1H1","DYNC1I2","DYNC1LI1","DYNC1LI2","DYNLL1","DYNLL2","KIF1A","KIF1B","KIF2A","KIF3A","KIF5A","KIF5B","KIF5C","KLC1","KLC2","KLC3","KLC4","DNM1","DNM2","DNM3","CLTC","CLTA","AP2B1","AP2A1","AP2A2","AP2M1"),
  "Microglia/Immune" = c("CD33","TREM2","TYROBP","CSF1R","CX3CR1","ITGAM","ITGB2","CD68","CD14","TLR2","TLR4","MYD88","NFKB1","NFKB2","RELA","RELB","NFKBIA","TNF","IL1B","IL6","IL10","TGFB1","TGFB2","C1QA","C1QB","C1QC","C3","C4A","C4B","SRC","SYK","LYN","HCK","FGR","BTK"),
  "Growth Factor Signaling" = c("AKT1","AKT2","AKT3","MTOR","RPS6KB1","RPS6KA1","RPS6KA2","RPS6KA3","RPS6KA4","PIK3CA","PIK3CB","PIK3R1","PIK3R2","IGF1R","INSR","IRS1","IRS2","EGFR","ERBB2","ERBB3","ERBB4","FGFR1","FGFR2","FGFR3","NTRK1","NTRK2","NTRK3","PDGFRA","PDGFRB","VEGFA","FLT1","KDR","GRB2","SOS1","SOS2","HRAS","KRAS","NRAS","BRAF","RAF1","MAP2K1","MAP2K2")
)

tau_axes <- list(
  "GSK3B axis" = c("MAPT","GSK3B","GSK3A","AKT1","AKT2","PTEN","PDPK1","PIK3CA","PIK3R1","PRKCA","PRKCB","PRKCG","PRKCZ","PPP1CA","PPP2CA","PPP2R1A","DYRK1A","WNT3A","FZD1","DVL1","AXIN1","APC","CTNNB1"),
  "MAPK cascade" = c("MAPT","MAPK1","MAPK3","MAPK8","MAPK9","MAPK14","MAP2K1","MAP2K2","MAP2K4","BRAF","RAF1","RASGRF1","GRB2","SOS1","HRAS","KRAS","RPS6KA3","RPS6KB1"),
  "CDK5 pathway" = c("MAPT","CDK5","CDK5R1","CDK5R2","PPP1CA","PPP1CB","PPP2CA","PPP3CA","PPP3R1","PPP5C","CAMK2A","PRKCA","ABL1","ABL2"),
  "SRC/FYN/SYK" = c("MAPT","SRC","FYN","SYK","LYN","HCK","LCK","PTK2","PTK2B","PTPN11","GRB2","PIK3R1"),
  "PPP Phosphatases" = c("MAPT","PPP1CA","PPP1CB","PPP1CC","PPP2CA","PPP2CB","PPP2R1A","PPP2R5A","PPP3CA","PPP3CB","PPP3CC","PPP3R1","PPP5C","PTEN")
)

tau_pathway_keywords <- c("tau","MAPK","PI3K","AKT","GSK","CDK5","phosphorylation","dephosphorylation","synaptic","synapse","neuroinflammation","microglia","cytoskeleton","axon","Alzheimer","neurodegener","amyloid","APP","calcium","CAMK","PKC","SRC","FYN","phosphatase","kinase","ubiquitin","autophagy","apoptosis","Wnt","Notch","mTOR")

# ---- 1. Analysis function ----
run_tau_mechanism <- function(df_grp, grp_label) {
  cat(sprintf("\n========== %s (%d subjects) ==========\n", grp_label, nrow(df_grp)))

  # Correlation
  df_grp$PTAU_log2 <- log2(df_grp$ptau217_csf)
  missing_rate <- sapply(df_grp[protein_cols], function(x) mean(is.na(x))*100)
  pk <- names(missing_rate[missing_rate <= 20])
  pl2 <- log2(as.matrix(df_grp[pk]))
  pl2[is.infinite(pl2)|is.nan(pl2)] <- NA
  ptau_vals <- df_grp$PTAU_log2

  results <- data.frame(protein_id=colnames(pl2), n=NA_integer_, spearman_rho=NA_real_, p_value=NA_real_)
  for (i in seq_len(ncol(pl2))) { pv <- pl2[,i]; vi <- !is.na(pv); nv <- sum(vi); results$n[i] <- nv
    if (nv >= 10) { tst <- tryCatch(cor.test(pv[vi], ptau_vals[vi], method="spearman", exact=FALSE), error=function(e) NULL)
      if (!is.null(tst)) { results$spearman_rho[i] <- tst$estimate; results$p_value[i] <- tst$p.value } } }
  results <- results[!is.na(results$p_value), ]; results$fdr <- p.adjust(results$p_value, method="BH")
  results$EntrezGeneSymbol <- gsub("\\.[0-9]+\\.[0-9]+$", "", results$protein_id)

  sig <- results[results$fdr < 0.05, ]
  sig_genes <- unique(sig$EntrezGeneSymbol[!is.na(sig$EntrezGeneSymbol) & sig$EntrezGeneSymbol != ""])
  cat(sprintf("FDR<0.05: %d proteins -> %d genes\n", nrow(sig), length(sig_genes)))

  gene_rho <- sig %>% filter(!is.na(EntrezGeneSymbol)&EntrezGeneSymbol!="") %>%
    group_by(EntrezGeneSymbol) %>% summarise(rho_median=median(spearman_rho,na.rm=TRUE), .groups="drop")

  # Network
  inter_sig <- omni[omni$source_genesymbol %in% sig_genes & omni$target_genesymbol %in% sig_genes, ]
  edges <- inter_sig[,c("source_genesymbol","target_genesymbol")]; colnames(edges) <- c("from","to")
  g <- graph_from_data_frame(edges, directed=FALSE)
  V(g)$rho <- gene_rho$rho_median[match(V(g)$name, gene_rho$EntrezGeneSymbol)]
  V(g)$degree <- degree(g)
  cat(sprintf("Network: %d nodes, %d edges\n", vcount(g), ecount(g)))

  # Communities
  comm <- cluster_louvain(g); V(g)$module <- comm$membership; n_mod <- max(comm$membership)

  # Module stats
  mod_stats <- data.frame()
  for (m in 1:n_mod) { nodes_m <- V(g)$name[V(g)$module==m]; n_m <- length(nodes_m)
    avg_deg <- if(n_m>1) mean(degree(induced_subgraph(g,nodes_m))) else 0
    top_hubs <- paste(names(head(sort(degree(induced_subgraph(g,nodes_m)),decreasing=TRUE),3)),collapse=", ")
    mod_stats <- rbind(mod_stats, data.frame(Module=m, Size=n_m, Avg_Deg=round(avg_deg,1),
      Avg_Rho=round(mean(V(g)$rho[V(g)$module==m],na.rm=TRUE),4), Top_Hubs=top_hubs)) }
  mod_stats <- mod_stats[order(mod_stats$Size, decreasing=TRUE), ]

  # Functional classification
  sig_gene_info <- gene_rho %>% mutate(Category1=NA_character_)
  for (i in seq_len(nrow(sig_gene_info))) {
    gn <- sig_gene_info$EntrezGeneSymbol[i]
    cats <- names(tau_categories)[sapply(tau_categories, function(gs) gn %in% gs)]
    if (length(cats)>=1) sig_gene_info$Category1[i] <- cats[1] }
  sig_gene_info$Category1[is.na(sig_gene_info$Category1)] <- "Others"
  cat_summary <- sig_gene_info %>% group_by(Category1) %>%
    summarise(N=n(), Avg_Rho=round(mean(rho_median,na.rm=TRUE),4),
              Top_Genes=paste(head(EntrezGeneSymbol[order(abs(rho_median),decreasing=TRUE)],3),collapse=", "), .groups="drop") %>% arrange(desc(N))

  # MAPT interactors
  mapt_edges <- omni[(omni$source_genesymbol=="MAPT"|omni$target_genesymbol=="MAPT") &
                     (omni$source_genesymbol %in% sig_genes | omni$target_genesymbol %in% sig_genes), ]
  mapt_genes <- setdiff(unique(c(mapt_edges$source_genesymbol, mapt_edges$target_genesymbol)), "MAPT")

  # Tau signaling axis coverage
  axis_sum <- data.frame()
  for (ax in names(tau_axes)) {
    ax_sig <- intersect(tau_axes[[ax]], sig_genes)
    axis_sum <- rbind(axis_sum, data.frame(Axis=ax, Total=length(tau_axes[[ax]]),
      In_Sig=length(ax_sig), Coverage=round(100*length(ax_sig)/length(tau_axes[[ax]]),1),
      Top_Hits=paste(head(ax_sig[order(abs(gene_rho$rho_median[match(ax_sig,gene_rho$EntrezGeneSymbol)]),decreasing=TRUE)],5),collapse=", "))) }

  # Candidates
  sig_gene_info$network_degree <- V(g)$degree[match(sig_gene_info$EntrezGeneSymbol, V(g)$name)]
  sig_gene_info$network_degree[is.na(sig_gene_info$network_degree)] <- 0
  sig_gene_info$is_MAPT <- sig_gene_info$EntrezGeneSymbol %in% mapt_genes
  sig_gene_info$in_tau_pathway <- sig_gene_info$Category1 != "Others"
  sig_gene_info$priority_score <- with(sig_gene_info, abs(rho_median)*log1p(network_degree)*(1+is_MAPT*0.5)*(1+in_tau_pathway*0.3))
  sig_gene_info <- sig_gene_info[order(sig_gene_info$priority_score, decreasing=TRUE), ]

  list(grp=grp_label, n_subj=nrow(df_grp), n_sig=nrow(sig), n_genes=length(sig_genes),
       n_interactions=nrow(inter_sig), n_nodes=vcount(g), n_edges=ecount(g),
       n_modules=n_mod, mod_stats=mod_stats, cat_summary=cat_summary,
       mapt_genes=mapt_genes, axis_sum=axis_sum,
       candidates=head(sig_gene_info[,c("EntrezGeneSymbol","rho_median","Category1","is_MAPT","priority_score")],20),
       gene_rho=gene_rho, sig_genes=sig_genes, g=g)
}

# ---- 2. Run both groups ----
df_cn <- df_bl[df_bl$dx_entry == "AD", ]
df_noncn <- df_bl[df_bl$dx_entry %in% c("CN", "EMCI", "LMCI"), ]
res_cn <- run_tau_mechanism(df_cn, "CN")
res_nc <- run_tau_mechanism(df_noncn, "non-CN")

# ---- 3. Comparison ----
cat("\n\n========== COMPARISON ==========\n")

# --- Functional category comparison ---
cat_all <- rbind(
  res_cn$cat_summary %>% mutate(Group="CN"),
  res_nc$cat_summary %>% mutate(Group="non-CN"))
cat_all$Category1 <- factor(cat_all$Category1,
  levels = cat_all %>% group_by(Category1) %>% summarise(s=sum(N)) %>% arrange(desc(s)) %>% pull(Category1))

p_cat <- ggplot(cat_all[cat_all$Category1!="Others",], aes(x=Category1, y=N, fill=Group)) +
  geom_col(position="dodge", alpha=0.85) +
  scale_fill_manual(values=c("AD"="#2c7bb6","non-AD"="#d7191c")) +
  labs(title="Functional Categories: AD vs non-AD", subtitle=sprintf("AD: %d genes | non-AD: %d genes", res_cn$n_genes, res_nc$n_genes),
       x="", y="Number of Genes") +
  theme_bw(11) + theme(plot.title=element_text(face="bold"), axis.text.x=element_text(angle=30,hjust=1))

ggsave(file.path(out_dir, "functional_categories_comparison.png"), p_cat, width=10, height=5.5, dpi=200)

# --- Module comparison ---
mod_all <- rbind(
  res_cn$mod_stats %>% mutate(Group="CN"),
  res_nc$mod_stats %>% mutate(Group="non-CN"))

p_mod <- ggplot(mod_all, aes(x=Size, fill=Group)) +
  geom_density(alpha=0.5) +
  scale_fill_manual(values=c("AD"="#2c7bb6","non-AD"="#d7191c")) +
  labs(title="Module Size Distribution: AD vs non-AD", x="Module Size (# genes)", y="Density") +
  theme_bw(11) + theme(plot.title=element_text(face="bold"))
ggsave(file.path(out_dir, "module_size_distribution.png"), p_mod, width=7, height=4.5, dpi=200)

# --- Tau axis coverage comparison ---
axis_all <- rbind(res_cn$axis_sum %>% mutate(Group="CN"),
                  res_nc$axis_sum %>% mutate(Group="non-CN"))

p_axis <- ggplot(axis_all, aes(x=Axis, y=Coverage, fill=Group)) +
  geom_col(position="dodge", alpha=0.85) +
  geom_text(aes(label=sprintf("%.0f%%",Coverage)), position=position_dodge(0.9), vjust=-0.3, size=3) +
  scale_fill_manual(values=c("AD"="#2c7bb6","non-AD"="#d7191c")) +
  labs(title="Tau Signaling Axis Coverage: AD vs non-AD",
       subtitle=sprintf("%% of curated axis genes that are PTAU-significant"), x="", y="Coverage (%)") +
  ylim(0,105) + theme_bw(11) + theme(plot.title=element_text(face="bold"))
ggsave(file.path(out_dir, "axis_coverage_comparison.png"), p_axis, width=8, height=5, dpi=200)

# --- MAPT interactor overlap ---
cn_mapt <- res_cn$mapt_genes; nc_mapt <- res_nc$mapt_genes
shared_mapt <- intersect(cn_mapt, nc_mapt); cn_only <- setdiff(cn_mapt, nc_mapt); nc_only <- setdiff(nc_mapt, cn_mapt)

mapt_venn <- data.frame(
  Category=c("CN only","Shared","non-CN only"),
  Count=c(length(cn_only), length(shared_mapt), length(nc_only)),
  Genes=c(paste(cn_only,collapse=", "), paste(head(shared_mapt,10),collapse=", "), paste(nc_only,collapse=", ")))
cat("\nMAPT interactors:\n")
print(mapt_venn)

p_venn <- ggplot(mapt_venn, aes(x=Category, y=Count, fill=Category)) +
  geom_col(alpha=0.85) + geom_text(aes(label=Count), vjust=-0.3, size=5) +
  scale_fill_manual(values=c("CN only"="#2c7bb6","Shared"="#984ea3","non-CN only"="#d7191c")) +
  labs(title="MAPT Interactors: AD vs non-AD", subtitle=sprintf("Shared: %d genes", length(shared_mapt)),
       x="", y="Number of Genes") +
  theme_bw(11) + theme(plot.title=element_text(face="bold"), legend.position="none")
ggsave(file.path(out_dir, "mapt_venn.png"), p_venn, width=6, height=4.5, dpi=200)

# --- Top candidates comparison ---
cand_all <- rbind(res_cn$candidates %>% mutate(Group="CN"),
                  res_nc$candidates %>% mutate(Group="non-CN"))
write.csv(cand_all, file.path(out_dir, "candidates_comparison.csv"), row.names=FALSE)

cat("\nTop 10 CN candidates:\n")
print(head(res_cn$candidates, 10))
cat("\nTop 10 non-CN candidates:\n")
print(head(res_nc$candidates, 10))

# ---- 4. Network comparison figures ----
# CN + non-CN networks side by side (top 150 genes each for readability)
for (res in list(res_cn, res_nc)) {
  grp_name <- res$grp
  g <- res$g
  if (vcount(g) > 150) g <- induced_subgraph(g, V(g)$name[order(V(g)$degree, decreasing=TRUE)][1:150])
  comm <- cluster_louvain(g); V(g)$module <- comm$membership
  mod_colors <- hcl.colors(max(comm$membership), "Dynamic")
  names(mod_colors) <- sort(unique(comm$membership))

  p <- ggraph(g, layout="fr") +
    geom_edge_link(color="grey90", alpha=0.15) +
    geom_node_point(aes(size=degree, fill=factor(module)), shape=21, color="grey30", stroke=0.2, alpha=0.85) +
    geom_node_text(aes(label=ifelse(degree>quantile(degree,0.95)|abs(rho)>0.7,name,"")),
                   size=2.5, repel=TRUE, max.overlaps=30, box.padding=0.3) +
    scale_fill_manual(values=mod_colors, name="Module", guide="none") +
    scale_size_continuous(range=c(1,6), guide="none") +
    labs(title=sprintf("%s: %d nodes", grp_name, vcount(g)),
         subtitle=sprintf("%d modules | Top %d by degree", max(comm$membership), vcount(g))) +
    theme_void() + theme(plot.title=element_text(face="bold",size=13), plot.subtitle=element_text(size=8,color="grey40"))
  ggsave(file.path(out_dir, sprintf("community_network_%s.png", gsub(" ","_",grp_name))), p, width=12, height=10, dpi=200)
  cat(sprintf("Saved: community_network_%s\n", grp_name))
}

# ---- 5. Summary table ----
sum_tbl <- data.frame(
  Metric=c("Subjects","FDR-sig proteins","Sig genes","Network nodes","Network edges",
           "Louvain modules","MAPT interactors","GSK3B axis cov","MAPK cascade cov",
           "CDK5 pathway cov","SRC/FYN/SYK cov","PPP phosphatase cov"),
  CN=c(res_cn$n_subj,res_cn$n_sig,res_cn$n_genes,res_cn$n_nodes,res_cn$n_edges,
       res_cn$n_modules,length(res_cn$mapt_genes),
       res_cn$axis_sum$Coverage[1],res_cn$axis_sum$Coverage[2],res_cn$axis_sum$Coverage[3],
       res_cn$axis_sum$Coverage[4],res_cn$axis_sum$Coverage[5]),
  nonCN=c(res_nc$n_subj,res_nc$n_sig,res_nc$n_genes,res_nc$n_nodes,res_nc$n_edges,
          res_nc$n_modules,length(res_nc$mapt_genes),
          res_nc$axis_sum$Coverage[1],res_nc$axis_sum$Coverage[2],res_nc$axis_sum$Coverage[3],
          res_nc$axis_sum$Coverage[4],res_nc$axis_sum$Coverage[5])
)
cat("\nSummary:\n")
print(sum_tbl, row.names=FALSE)
write.csv(sum_tbl, file.path(out_dir, "comparison_summary.csv"), row.names=FALSE)

cat(sprintf("\nOutput: %s\n", normalizePath(out_dir)))
cat("========== Done ==========\n")
