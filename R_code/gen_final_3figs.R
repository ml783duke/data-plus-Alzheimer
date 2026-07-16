###############################################################################
# Final 3 poster figures: Lifestyle + MAPT peak + Tau network
###############################################################################
library(ggplot2); library(dplyr); library(tidyr); library(igraph); library(ggraph)

out_dir <- "output_figure"
DUKE_NAVY  <- "#003087"; DUKE_BLUE <- "#00539B"; DUKE_LIGHT <- "#7288A0"; DUKE_ORANGE <- "#C84E00"
dx_colors  <- c("CN"="#8BA5BF","EMCI"=DUKE_NAVY,"LMCI"=DUKE_BLUE,"AD"=DUKE_ORANGE)

theme_post <- theme_minimal(base_size=16) +
  theme(plot.title=element_text(face="bold",size=22,color=DUKE_NAVY),
        plot.subtitle=element_text(size=14,color=DUKE_LIGHT),
        plot.margin=margin(12,12,12,12),
        legend.position="bottom",
        panel.grid.minor=element_blank(),
        panel.grid.major=element_line(color="grey92",linewidth=0.3),
        axis.title=element_text(size=16,color=DUKE_NAVY),
        axis.text=element_text(size=14))

cat("=== FIGURE 1: Lifestyle → pTau217 Forest Plot ===\n")

csv <- read.csv("lifestyle_7var_ptau.csv")
master <- readxl::read_excel("master_data.xlsx", sheet="Sheet1")
sl <- master[!duplicated(master$RID), c("RID","PTGENDER")]
csv <- merge(csv, sl, by="RID", all.x=TRUE)
csv$SEX <- ifelse(csv$PTGENDER=="Male", 1L, 0L)
df_a <- csv[!is.na(csv$ptau217_csf) & !is.na(csv$age_at_csf_ptau), ]
df_a$PTAU_log2 <- log2(df_a$ptau217_csf)

VARS <- c("DHA","EPA","HCys","NPIK","NPIKTOT","MH14ALCH","MH16SMOK")
VT   <- c("DHA"="cont","EPA"="cont","HCys"="cont","NPIK"="bin","NPIKTOT"="cont","MH14ALCH"="bin","MH16SMOK"="bin")
vl   <- c("DHA"="DHA (n-3 PUFA)","EPA"="EPA (n-3 PUFA)","HCys"="Homocysteine",
          "NPIK"="NPI Sleep Domain","NPIKTOT"="NPI Total Score",
          "MH14ALCH"="Alcohol Abuse History","MH16SMOK"="Smoking Status")

uni <- data.frame(variable=VARS, n=NA_integer_, beta=NA_real_, se=NA_real_,
                  ci_low=NA_real_, ci_high=NA_real_, p_value=NA_real_)
for (i in seq_along(VARS)) {
  v <- VARS[i]; ds <- df_a[!is.na(df_a[[v]]), ]; uni$n[i] <- nrow(ds)
  xv <- ds[[v]]
  xm <- if (VT[v]=="cont" && sd(xv,na.rm=TRUE) > 0) as.vector(scale(xv)) else xv
  fit <- lm(PTAU_log2 ~ xm + age_at_csf_ptau, data=ds); s <- summary(fit)
  uni$beta[i]  <- s$coefficients["xm","Estimate"]
  uni$se[i]    <- s$coefficients["xm","Std. Error"]
  uni$p_value[i] <- s$coefficients["xm","Pr(>|t|)"]
  uni$ci_low[i]  <- uni$beta[i] - 1.96*uni$se[i]
  uni$ci_high[i] <- uni$beta[i] + 1.96*uni$se[i]
}
uni$fdr <- p.adjust(uni$p_value, method="BH")
uni$significant <- uni$fdr < 0.05
uni$label <- vl[uni$variable]
uni <- uni[order(uni$beta), ]
uni$label <- factor(uni$label, levels=uni$label)

p1 <- ggplot(uni, aes(x=beta, y=label)) +
  geom_vline(xintercept=0, linetype="dashed", color="grey60", linewidth=0.7) +
  geom_point(aes(size=n), color=DUKE_BLUE, alpha=0.9) +
  geom_errorbarh(aes(xmin=ci_low, xmax=ci_high), height=0.3, color=DUKE_BLUE, linewidth=1) +
  scale_size_continuous(name="N", range=c(4, 8)) +
  labs(title="Lifestyle Factors Are Not Directly Associated with CSF pTau217",
       subtitle=sprintf("Univariate models: log2(pTau217) ~ variable + Age | N=%d | 0/7 FDR < 0.05 | All effects cross zero",
                        nrow(df_a)),
       x="Standardized Effect (per 1-SD or 1 vs 0)", y="") +
  theme_post + theme(panel.grid.major.y=element_blank(),
                      axis.text.y=element_text(size=16),
                      legend.text=element_text(size=12))

ggsave(file.path(out_dir,"Fig1_lifestyle_forest.png"), p1, width=12, height=8, dpi=300)
cat("Fig1 done\n")

cat("\n=== FIGURE 2: MAPT Interactors by DX (EMCI peak) ===\n")

mapt_mat <- read.csv("output/tau_mechanism_by_DX/mapt_interactors_by_DX.csv")
mapt_counts <- data.frame(
  DX=c("CN","EMCI","LMCI","AD"),
  Count=c(sum(mapt_mat$CN),sum(mapt_mat$EMCI),sum(mapt_mat$LMCI),sum(mapt_mat$AD)))
mapt_counts$DX <- factor(mapt_counts$DX, levels=c("CN","EMCI","LMCI","AD"))

p2 <- ggplot(mapt_counts, aes(x=DX, y=Count, fill=DX)) +
  geom_col(alpha=0.9, width=0.55) +
  geom_text(aes(label=Count), vjust=-0.5, size=10, color=DUKE_NAVY, fontface="bold") +
  scale_fill_manual(values=dx_colors, guide="none") +
  scale_y_continuous(limits=c(0, 28), expand=c(0,0)) +
  labs(title="Tau-Interacting Proteins Peak at EMCI",
       subtitle=sprintf("MAPT-interacting proteins significantly correlated with pTau217 | %d total unique",
                        nrow(mapt_mat)),
       x="", y="Number of MAPT-Interacting Proteins") +
  theme_post

ggsave(file.path(out_dir,"Fig2_mapt_EMCI_peak.png"), p2, width=10, height=8, dpi=300)
cat("Fig2 done\n")

cat("\n=== FIGURE 3: Tau Signaling Map ===\n")

omni_url <- "https://omnipathdb.org/interactions?format=tsv&fields=sources,references&genesymbols=1"
tmp <- tempfile(fileext=".tsv"); download.file(omni_url,tmp,method="auto")
omni <- read.delim(tmp,stringsAsFactors=FALSE)
sig <- read.csv("output/tau_mechanism_analysis/gene_functional_classification.csv")
gene_rho <- sig %>% filter(EntrezGeneSymbol!="",!is.na(rho_median)) %>%
  group_by(EntrezGeneSymbol) %>% summarise(rho=median(rho_median,na.rm=TRUE),.groups="drop")

tau_axes <- list(
  "GSK3B axis"=c("GSK3B","GSK3A","AKT1","AKT2","PTEN","PDPK1","PRKCA","PRKCB","PRKCG","PRKCZ","PPP2R1A","DYRK1A"),
  "MAPK cascade"=c("MAPK1","MAPK3","MAPK8","MAPK14","MAP2K1","MAP2K2","BRAF","RAF1","GRB2","HRAS","KRAS","RPS6KA3","RPS6KB1"),
  "CDK5 pathway"=c("CDK5","PPP3CA","PPP3R1","PPP5C","CAMK2A","ABL1","ABL2"),
  "SRC/FYN/SYK"=c("SRC","FYN","SYK","LYN","PTK2","PTK2B","PTPN11"),
  "PPP phosphatase"=c("PPP1CA","PPP1CB","PPP1CC","PPP2CA","PPP2CB","PPP2R5A","PPP3CB","PPP3CC"))

focus <- unique(c("MAPT", intersect(unlist(tau_axes), gene_rho$EntrezGeneSymbol)))
inter_focus <- omni[omni$source_genesymbol %in% focus & omni$target_genesymbol %in% focus, ]
edges <- data.frame(from=inter_focus$source_genesymbol, to=inter_focus$target_genesymbol)
nodes <- data.frame(gene=unique(c(edges$from, edges$to)))
nodes$rho <- gene_rho$rho[match(nodes$gene, gene_rho$EntrezGeneSymbol)]
nodes$axis <- "Other"
for (ax in names(tau_axes)) nodes$axis[nodes$gene %in% tau_axes[[ax]]] <- ax
nodes$axis[nodes$gene=="MAPT"] <- "MAPT"
g_tau <- graph_from_data_frame(edges, directed=FALSE, vertices=nodes)
V(g_tau)$degree <- degree(g_tau)

ax_col <- c("MAPT"="black","GSK3B axis"=DUKE_BLUE,"MAPK cascade"=DUKE_ORANGE,
            "CDK5 pathway"=DUKE_LIGHT,"SRC/FYN/SYK"="#1B6B3A",
            "PPP phosphatase"="#6B3A8A","Other"="grey85")

set.seed(42)
layout_df <- create_layout(g_tau, layout="stress")
mapt_xy <- layout_df[layout_df$name=="MAPT", c("x","y")]

p3 <- ggraph(g_tau, layout="stress") +
  geom_edge_link(color="grey55", alpha=0.7) +
  geom_node_point(aes(size=degree, fill=axis), shape=21, color="grey40", stroke=0.3) +
  geom_node_text(aes(label=ifelse(name=="MAPT","",name), size=3.5),
                 repel=TRUE, max.overlaps=60, box.padding=0.6, color=DUKE_NAVY, force=2) +
  annotate("text", x=mapt_xy$x, y=mapt_xy$y-0.3, label="MAPT", size=9, fontface="bold", color="black") +
  scale_fill_manual(values=ax_col, name="Signaling Axis") +
  scale_size_continuous(range=c(2,10), guide="none") +
  labs(title="Five Signaling Axes Converge on Tau",
       subtitle=sprintf("%d nodes, %d edges | Multi-pathway redundancy may explain single-target trial failures",
                        vcount(g_tau), ecount(g_tau))) +
  theme_void() +
  theme(plot.title=element_text(face="bold",size=22,color=DUKE_NAVY,hjust=0.5),
        plot.subtitle=element_text(size=14,color=DUKE_LIGHT,hjust=0.5),
        plot.margin=margin(3,3,3,3),
        legend.position="bottom", legend.text=element_text(size=14),
        legend.title=element_text(size=15))

ggsave(file.path(out_dir,"Fig3_tau_signaling_map.png"), p3, width=14, height=12, dpi=300)
cat("Fig3 done\n\nALL DONE\n")
