###############################################################################
# 14_tau_mechanism_by_DX.R
# Tau Mechanism Analysis by Individual DX Groups (CN, EMCI, LMCI, AD)
###############################################################################

library(readxl); library(dplyr); library(igraph); library(ggraph)
library(ggplot2); library(clusterProfiler); library(org.Hs.eg.db)
library(enrichplot); library(patchwork)

set.seed(42)
out_dir <- "output/tau_mechanism_by_DX"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("========== Tau Mechanism by DX ==========\n")

# ---- Load ----
master <- read.csv("lifestyle_7var_ptau_proteomics.csv")
df_bl <- master[!is.na(master$ptau217_csf) & !is.na(master$dx_entry), ]
protein_cols <- colnames(df_bl)[21:ncol(df_bl)]
convert_protein <- function(x) { x[x == "NA" | x == ""] <- NA; as.numeric(x) }
df_bl[protein_cols] <- lapply(df_bl[protein_cols], convert_protein)

omni_url <- "https://omnipathdb.org/interactions?format=tsv&fields=sources,references&genesymbols=1"
tmp <- tempfile(fileext = ".tsv"); download.file(omni_url, tmp, method = "auto")
omni <- read.delim(tmp, stringsAsFactors = FALSE)

# Tau categories
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

DX_GROUPS <- c("CN", "EMCI", "LMCI", "AD")

# ---- Analysis function ----
run_tau_dx <- function(df_grp, grp_label) {
  cat(sprintf("\n--- %s (%d) ---\n", grp_label, nrow(df_grp)))
  df_grp$PTAU_log2 <- log2(df_grp$ptau217_csf)
  missing_rate <- sapply(df_grp[protein_cols], function(x) mean(is.na(x))*100)
  pk <- names(missing_rate[missing_rate <= 20])
  pl2 <- log2(as.matrix(df_grp[pk])); pl2[is.infinite(pl2)|is.nan(pl2)] <- NA
  ptau_vals <- df_grp$PTAU_log2

  results <- data.frame(protein_id=colnames(pl2), n=NA_integer_, spearman_rho=NA_real_, p_value=NA_real_)
  for (i in seq_len(ncol(pl2))) { pv <- pl2[,i]; vi <- !is.na(pv); nv <- sum(vi); results$n[i] <- nv
    if (nv >= 10) { tst <- tryCatch(cor.test(pv[vi], ptau_vals[vi], method="spearman", exact=FALSE), error=function(e) NULL)
      if (!is.null(tst)) { results$spearman_rho[i] <- tst$estimate; results$p_value[i] <- tst$p.value } } }
  results <- results[!is.na(results$p_value), ]; results$fdr <- p.adjust(results$p_value, method="BH")
  # Extract gene symbol from protein_id (strip SomaScan suffix like .10010.10)
  results$EntrezGeneSymbol <- gsub("\\.[0-9]+\\.[0-9]+$", "", results$protein_id)

  sig <- results[results$fdr < 0.05, ]
  sig_genes <- unique(sig$EntrezGeneSymbol[!is.na(sig$EntrezGeneSymbol) & sig$EntrezGeneSymbol != ""])
  cat(sprintf("Sig: %d proteins -> %d genes\n", nrow(sig), length(sig_genes)))

  gene_rho <- sig %>% filter(!is.na(EntrezGeneSymbol)&EntrezGeneSymbol!="") %>%
    group_by(EntrezGeneSymbol) %>% summarise(rho_median=median(spearman_rho,na.rm=TRUE), .groups="drop")

  inter_sig <- omni[omni$source_genesymbol %in% sig_genes & omni$target_genesymbol %in% sig_genes, ]
  edges <- inter_sig[,c("source_genesymbol","target_genesymbol")]; colnames(edges) <- c("from","to")
  g <- graph_from_data_frame(edges, directed=FALSE)
  V(g)$rho <- gene_rho$rho_median[match(V(g)$name, gene_rho$EntrezGeneSymbol)]
  V(g)$degree <- degree(g)

  mapt_edges <- omni[(omni$source_genesymbol=="MAPT"|omni$target_genesymbol=="MAPT") &
                     (omni$source_genesymbol %in% sig_genes | omni$target_genesymbol %in% sig_genes), ]
  mapt_genes <- setdiff(unique(c(mapt_edges$source_genesymbol, mapt_edges$target_genesymbol)), "MAPT")

  sig_gi <- gene_rho %>% mutate(Category1=NA_character_)
  for (i in seq_len(nrow(sig_gi))) {
    gn <- sig_gi$EntrezGeneSymbol[i]
    cats <- names(tau_categories)[sapply(tau_categories, function(gs) gn %in% gs)]
    if (length(cats)>=1) sig_gi$Category1[i] <- cats[1] }
  sig_gi$Category1[is.na(sig_gi$Category1)] <- "Others"
  cat_sum <- sig_gi %>% group_by(Category1) %>%
    summarise(N=n(), .groups="drop") %>% arrange(desc(N))

  axis_sum <- data.frame()
  for (ax in names(tau_axes)) {
    ax_sig <- intersect(tau_axes[[ax]], sig_genes)
    axis_sum <- rbind(axis_sum, data.frame(Axis=ax, Total=length(tau_axes[[ax]]),
      In_Sig=length(ax_sig), Coverage=round(100*length(ax_sig)/length(tau_axes[[ax]]),1))) }

  sig_gi$network_degree <- V(g)$degree[match(sig_gi$EntrezGeneSymbol, V(g)$name)]
  sig_gi$network_degree[is.na(sig_gi$network_degree)] <- 0
  sig_gi$is_MAPT <- sig_gi$EntrezGeneSymbol %in% mapt_genes
  sig_gi$in_tau_pathway <- sig_gi$Category1 != "Others"
  sig_gi$priority_score <- with(sig_gi, abs(rho_median)*log1p(network_degree)*(1+is_MAPT*0.5)*(1+in_tau_pathway*0.3))
  sig_gi <- sig_gi[order(sig_gi$priority_score, decreasing=TRUE), ]

  list(grp=grp_label, n_subj=nrow(df_grp), n_sig=nrow(sig), n_genes=length(sig_genes),
       n_nodes=vcount(g), n_edges=ecount(g), mapt_genes=mapt_genes,
       cat_sum=cat_sum, axis_sum=axis_sum,
       candidates=head(sig_gi[,c("EntrezGeneSymbol","rho_median","Category1","is_MAPT","priority_score")],10),
       gene_rho=gene_rho, sig_genes=sig_genes, g=g)
}

# ---- Run all 4 DX groups ----
res_all <- list()
for (g in DX_GROUPS) {
  df_g <- df_bl[df_bl$dx_entry == g, ]
  if (nrow(df_g) < 30) next
  res_all[[g]] <- run_tau_dx(df_g, g)
}

# ---- Comparison ----
cat("\n\n========== COMPARISON ==========\n")

# 1. Axis coverage
axis_comp <- do.call(rbind, lapply(names(res_all), function(g) {
  res_all[[g]]$axis_sum %>% mutate(DX = g) }))

p_axis <- ggplot(axis_comp, aes(x=Axis, y=Coverage, fill=DX)) +
  geom_col(position="dodge", alpha=0.85) +
  geom_text(aes(label=sprintf("%.0f%%",Coverage)), position=position_dodge(0.9), vjust=-0.3, size=2.8) +
  scale_fill_manual(values=c("CN"="#2c7bb6","EMCI"="#fdae61","LMCI"="#1b9e77","AD"="#d7191c")) +
  labs(title="Tau Signaling Axis Coverage by DX",
       subtitle="% of curated axis genes that are PTAU-significant", x="", y="Coverage (%)") +
  ylim(0,110) + theme_bw(11) + theme(plot.title=element_text(face="bold"))
ggsave(file.path(out_dir, "axis_coverage_by_DX.png"), p_axis, width=10, height=5.5, dpi=200)

# 2. Functional categories
cat_comp <- do.call(rbind, lapply(names(res_all), function(g) {
  res_all[[g]]$cat_sum %>% mutate(DX = g) }))
cat_comp <- cat_comp[cat_comp$Category1 != "Others", ]

p_cat <- ggplot(cat_comp, aes(x=Category1, y=N, fill=DX)) +
  geom_col(position="dodge", alpha=0.85) +
  scale_fill_manual(values=c("CN"="#2c7bb6","EMCI"="#fdae61","LMCI"="#1b9e77","AD"="#d7191c")) +
  labs(title="Functional Categories by DX", x="", y="Number of Genes") +
  theme_bw(11) + theme(plot.title=element_text(face="bold"), axis.text.x=element_text(angle=30,hjust=1))
ggsave(file.path(out_dir, "functional_categories_by_DX.png"), p_cat, width=10, height=5.5, dpi=200)

# 3. MAPT overlap
mapt_list <- lapply(res_all, function(r) r$mapt_genes)
mapt_all <- unique(unlist(mapt_list))
mapt_mat <- data.frame(Gene=mapt_all)
for (g in DX_GROUPS) mapt_mat[[g]] <- mapt_mat$Gene %in% mapt_list[[g]]
mapt_mat$n_DX <- rowSums(mapt_mat[,DX_GROUPS])
mapt_mat <- mapt_mat[order(mapt_mat$n_DX, decreasing=TRUE), ]
write.csv(mapt_mat, file.path(out_dir, "mapt_interactors_by_DX.csv"), row.names=FALSE)

p_mapt <- ggplot(mapt_mat, aes(x=factor(n_DX), fill=factor(n_DX))) +
  geom_bar(alpha=0.85) + geom_text(stat="count", aes(label=after_stat(count)), vjust=-0.3, size=4) +
  scale_fill_manual(values=c("1"="grey70","2"="#fdae61","3"="#1b9e77","4"="#d7191c")) +
  labs(title="MAPT Interactors: Conservation Across DX",
       subtitle=sprintf("Total unique: %d | In all 4 DX: %.0f genes", nrow(mapt_mat), sum(mapt_mat$n_DX==4)),
       x="# DX groups present", y="Count") +
  theme_bw(11) + theme(plot.title=element_text(face="bold"), legend.position="none")
ggsave(file.path(out_dir, "mapt_conservation_by_DX.png"), p_mapt, width=6, height=4.5, dpi=200)

# 4. Summary table
sum_tbl <- data.frame(Metric=c("Subjects","FDR-sig","Sig genes","Nodes","Edges","MAPT interactors",
  "GSK3B axis","MAPK cascade","CDK5 pathway","SRC/FYN/SYK","PPP phosphatase"))
for (g in DX_GROUPS) {
  r <- res_all[[g]]
  sum_tbl[[g]] <- c(r$n_subj, r$n_sig, r$n_genes, r$n_nodes, r$n_edges, length(r$mapt_genes),
    r$axis_sum$Coverage[1:5])
}
print(sum_tbl, row.names=FALSE)
write.csv(sum_tbl, file.path(out_dir, "comparison_summary.csv"), row.names=FALSE)

# 5. Candidates per DX
cand_all <- do.call(rbind, lapply(names(res_all), function(g) {
  res_all[[g]]$candidates %>% mutate(DX=g) }))
write.csv(cand_all, file.path(out_dir, "candidates_by_DX.csv"), row.names=FALSE)

cat(sprintf("\nOutput: %s\n", normalizePath(out_dir)))
cat("========== Done ==========\n")
