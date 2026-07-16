###############################################################################
# 3 standalone poster figures
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

cat("=== FIGURE 1: pTau217 LOESS by DX ===\n")

csv <- read.csv("lifestyle_7var_ptau_proteomics.csv")
df_l <- csv %>% mutate(p=as.numeric(ptau217_csf)) %>%
  filter(!is.na(p), !is.na(age_at_csf_ptau), !is.na(dx_entry),
         dx_entry %in% c("CN","EMCI","LMCI","AD"))
df_l$dx_entry <- factor(df_l$dx_entry, levels=c("CN","EMCI","LMCI","AD"))

fits <- df_l %>% group_by(dx_entry) %>% group_modify(~ {
  av <- .x$age_at_csf_ptau; pv <- .x$ptau217_csf; vi <- !is.na(av)&!is.na(pv)
  if (sum(vi)<10) return(tibble::tibble())
  fd <- data.frame(x=av[vi], y=pv[vi])
  lo <- tryCatch(loess(y~x,span=0.75,degree=1,data=fd), error=function(e)NULL)
  if (is.null(lo)) return(tibble::tibble())
  aseq <- seq(min(av,na.rm=TRUE), max(av,na.rm=TRUE), length.out=100)
  tibble::tibble(age=aseq, fitted=predict(lo, data.frame(x=aseq)))
}) %>% filter(n()>0)

label_pos <- fits %>% group_by(dx_entry) %>% slice_min(age,n=1) %>% ungroup()
label_pos$y_nudge <- 0
label_pos$y_nudge[label_pos$dx_entry=="EMCI"] <- -0.4
label_pos$y_nudge[label_pos$dx_entry=="LMCI"] <-  0.6
label_pos$y_nudge[label_pos$dx_entry=="CN"]   <- -0.5
label_pos$y_nudge[label_pos$dx_entry=="AD"]   <-  0.4

p1 <- ggplot() +
  geom_point(data=df_l, aes(x=age_at_csf_ptau,y=ptau217_csf,color=dx_entry), alpha=0.15, size=1) +
  geom_line(data=fits, aes(x=age,y=fitted,color=dx_entry), linewidth=1.6) +
  geom_text(data=label_pos, aes(x=age+0.8,y=fitted+y_nudge,label=dx_entry,color=dx_entry),
            hjust=0,vjust=0.5,size=7,fontface="bold") +
  scale_color_manual(values=dx_colors) +
  labs(title="CSF pTau217 Age Trajectories by Disease Stage",
       subtitle=sprintf("LOESS (degree=1, span=0.75) | %d subjects | AD shows paradoxical age-related decline",
                        sum(table(df_l$dx_entry))),
       x="Age (years)", y="CSF pTau217 (pg/mL)") +
  theme_post + theme(legend.position="none")

ggsave(file.path(out_dir,"Fig1_loess_pTau217.png"), p1, width=12, height=8, dpi=300)
cat("Fig1 done\n")

cat("\n=== FIGURE 2: MAPT Interactors by DX (EMCI peak) ===\n")

mapt_mat <- read.csv("output/tau_mechanism_by_DX/mapt_interactors_by_DX.csv")
mapt_counts <- data.frame(
  DX=c("CN","EMCI","LMCI","AD"),
  Count=c(sum(mapt_mat$CN),sum(mapt_mat$EMCI),sum(mapt_mat$LMCI),sum(mapt_mat$AD)))
mapt_counts$DX <- factor(mapt_counts$DX, levels=c("CN","EMCI","LMCI","AD"))

p2 <- ggplot(mapt_counts, aes(x=DX, y=Count, fill=DX)) +
  geom_col(alpha=0.9, width=0.55) +
  geom_text(aes(label=Count), vjust=-0.5, size=9, color=DUKE_NAVY, fontface="bold") +
  scale_fill_manual(values=dx_colors, guide="none") +
  scale_y_continuous(limits=c(0, 28), expand=c(0,0)) +
  labs(title="Tau-Interacting Proteins Peak at EMCI",
       subtitle=sprintf("MAPT interactors significantly correlated with pTau217 | %d total unique | CN=6, AD=6, EMCI=22",
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

mapt_sig <- omni[(omni$source_genesymbol=="MAPT"|omni$target_genesymbol=="MAPT") &
  (omni$source_genesymbol %in% gene_rho$EntrezGeneSymbol | omni$target_genesymbol %in% gene_rho$EntrezGeneSymbol), ]

tau_axes <- list(
  "GSK3B axis"=c("GSK3B","GSK3A","AKT1","AKT2","PTEN","PDPK1","PRKCA","PRKCB","PRKCG","PRKCZ","PPP2R1A","DYRK1A"),
  "MAPK cascade"=c("MAPK1","MAPK3","MAPK8","MAPK14","MAP2K1","MAP2K2","BRAF","RAF1","GRB2","HRAS","KRAS","RPS6KA3","RPS6KB1"),
  "CDK5 pathway"=c("CDK5","PPP3CA","PPP3R1","PPP5C","CAMK2A","ABL1","ABL2"),
  "SRC/FYN/SYK"=c("SRC","FYN","SYK","LYN","PTK2","PTK2B","PTPN11"),
  "PPP phosphatase"=c("PPP1CA","PPP1CB","PPP1CC","PPP2CA","PPP2CB","PPP2R5A","PPP3CB","PPP3CC"))

focus <- unique(c("MAPT", intersect(unlist(tau_axes), gene_rho$EntrezGeneSymbol)))

# Get ALL interactions among focus genes (not just MAPT-centric)
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
  labs(title="Multi-Axis Convergence on Tau",
       subtitle=sprintf("%d nodes, %d edges | 5 signaling axes converge on MAPT", vcount(g_tau), ecount(g_tau))) +
  theme_void() +
  theme(plot.title=element_text(face="bold",size=22,color=DUKE_NAVY,hjust=0.5),
        plot.subtitle=element_text(size=14,color=DUKE_LIGHT,hjust=0.5),
        plot.margin=margin(3,3,3,3),
        legend.position="bottom", legend.text=element_text(size=14),
        legend.title=element_text(size=15))

ggsave(file.path(out_dir,"Fig3_tau_signaling_map.png"), p3, width=14, height=12, dpi=300)
cat("Fig3 done\n\nALL DONE\n")
