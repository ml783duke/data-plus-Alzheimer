###############################################################################
# Regenerate all 3 poster figures with LARGER fonts
###############################################################################
library(ggplot2); library(dplyr); library(tidyr); library(igraph); library(ggraph); library(patchwork)

out_dir <- "output_figure"
DUKE_NAVY  <- "#003087"; DUKE_BLUE <- "#00539B"; DUKE_LIGHT <- "#7288A0"; DUKE_ORANGE <- "#C84E00"
dx_colors  <- c("CN"="#8BA5BF","EMCI"=DUKE_NAVY,"LMCI"=DUKE_BLUE,"AD"=DUKE_ORANGE)

# Base theme with larger fonts
theme_post <- theme_minimal(base_size = 14) +
  theme(plot.title = element_text(face="bold", size=18, color=DUKE_NAVY),
        plot.subtitle = element_text(size=12, color=DUKE_LIGHT),
        plot.margin = margin(8, 10, 8, 10),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color="grey92", linewidth=0.3),
        axis.title = element_text(size=14, color=DUKE_NAVY),
        axis.text = element_text(size=12),
        strip.text = element_text(face="bold", size=13, color=DUKE_NAVY))

cat("========== FIGURE 1: LOESS + Lifestyle ==========\n")

# ---- Left: LOESS ----
csv <- read.csv("lifestyle_7var_ptau_proteomics.csv")
df_l <- csv %>% mutate(p=as.numeric(ptau217_csf)) %>%
  filter(!is.na(p), !is.na(age_at_csf_ptau), !is.na(dx_entry),
         dx_entry %in% c("CN","EMCI","LMCI","AD"))
df_l$dx_entry <- factor(df_l$dx_entry, levels=c("CN","EMCI","LMCI","AD"))

fits <- df_l %>% group_by(dx_entry) %>% group_modify(~ {
  av <- .x$age_at_csf_ptau; pv <- .x$ptau217_csf; vi <- !is.na(av) & !is.na(pv)
  if (sum(vi) < 10) return(tibble::tibble())
  fd <- data.frame(x=av[vi], y=pv[vi])
  lo <- tryCatch(loess(y~x, span=0.75, degree=1, data=fd), error=function(e) NULL)
  if (is.null(lo)) return(tibble::tibble())
  aseq <- seq(min(av,na.rm=TRUE), max(av,na.rm=TRUE), length.out=100)
  tibble::tibble(age=aseq, fitted=predict(lo, data.frame(x=aseq)))
}) %>% filter(n() > 0)
label_pos <- fits %>% group_by(dx_entry) %>% slice_min(age, n=1) %>% ungroup()
# Place labels at curve start, nudged to avoid overlap
# Fine-tune y-placement so labels sit clearly ABOVE their lines
# Fine-tune y-placement: all nudged slightly down
label_pos$y_nudge <- 0
label_pos$y_nudge[label_pos$dx_entry == "EMCI"] <- -0.4
label_pos$y_nudge[label_pos$dx_entry == "LMCI"] <-  0.6
label_pos$y_nudge[label_pos$dx_entry == "CN"]   <- -0.5
label_pos$y_nudge[label_pos$dx_entry == "AD"]   <-  0.4

# Only use label_pos (not duplicate from df_l)
pL <- ggplot() +
  geom_point(data=df_l, aes(x=age_at_csf_ptau, y=ptau217_csf, color=dx_entry),
             alpha=0.15, size=0.9) +
  geom_line(data=fits, aes(x=age, y=fitted, color=dx_entry), linewidth=1.5) +
  geom_text(data=label_pos, aes(x=age+0.8, y=fitted + y_nudge, label=dx_entry, color=dx_entry),
            hjust=0, vjust=0.5, size=6, fontface="bold") +
  scale_color_manual(values=dx_colors) +
  labs(title="CSF pTau217 Age Trajectories by Disease Stage",
       x="Age (years)", y="CSF pTau217 (pg/mL)") +
  theme_post

# ---- Right: Forest ----
csv2 <- read.csv("lifestyle_7var_ptau.csv")
master <- readxl::read_excel("master_data.xlsx", sheet="Sheet1")
sl <- master[!duplicated(master$RID), c("RID","PTGENDER")]
csv2 <- merge(csv2, sl, by="RID", all.x=TRUE)
csv2$SEX <- ifelse(csv2$PTGENDER=="Male", 1L, 0L)
df_a <- csv2[!is.na(csv2$ptau217_csf) & !is.na(csv2$age_at_csf_ptau), ]
df_a$PTAU_log2 <- log2(df_a$ptau217_csf)

VARS <- c("DHA","EPA","HCys","NPIK","NPIKTOT","MH14ALCH","MH16SMOK")
VT   <- c("DHA"="cont","EPA"="cont","HCys"="cont","NPIK"="bin","NPIKTOT"="cont","MH14ALCH"="bin","MH16SMOK"="bin")
vl   <- c("DHA"="DHA (n-3 PUFA)","EPA"="EPA (n-3 PUFA)","HCys"="Homocysteine","NPIK"="NPI Sleep","NPIKTOT"="NPI Total","MH14ALCH"="Alcohol Abuse","MH16SMOK"="Smoking")

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
uni$label <- vl[uni$variable]
uni <- uni[order(uni$beta), ]
uni$label <- factor(uni$label, levels=uni$label)

pR <- ggplot(uni, aes(x=beta, y=label)) +
  geom_vline(xintercept=0, linetype="dashed", color="grey60", linewidth=0.5) +
  geom_point(aes(size=n), color=DUKE_BLUE, alpha=0.9) +
  geom_errorbarh(aes(xmin=ci_low, xmax=ci_high), height=0.25, color=DUKE_BLUE, linewidth=0.9) +
  scale_size_continuous(name="N", range=c(3,6)) +
  labs(title="Lifestyle Factors -> CSF pTau217",
       subtitle=sprintf("log2(pTau)~variable+Age | N=%d | 0/7 FDR<0.05", nrow(df_a)),
       x="Standardized Effect", y="") +
  theme_post + theme(legend.position="bottom", panel.grid.major.y=element_blank(),
                      axis.text.y=element_text(size=13), legend.text=element_text(size=11))

fig1 <- pL + pR + plot_layout(widths=c(1, 0.7)) +
  plot_annotation(title="pTau217 Age Trajectories & Lifestyle Associations",
    subtitle="Left: LOESS (degree=1, span=0.75) by DX | Right: Univariate models adjusted for age",
    theme=theme(plot.title=element_text(face="bold", size=20, color=DUKE_NAVY),
                plot.subtitle=element_text(size=13, color=DUKE_LIGHT)))
ggsave(file.path(out_dir, "Fig1_loess_lifestyle.png"), fig1, width=16, height=6.5, dpi=300)
cat("Fig1 done\n")

# ===========================================================================
cat("========== FIGURE 2: MAPT + Axis Coverage ==========\n")

mapt_mat <- read.csv("output/tau_mechanism_by_DX/mapt_interactors_by_DX.csv")
mapt_counts <- data.frame(
  DX = c("CN","EMCI","LMCI","AD"),
  Count = c(sum(mapt_mat$CN), sum(mapt_mat$EMCI), sum(mapt_mat$LMCI), sum(mapt_mat$AD)))
mapt_counts$DX <- factor(mapt_counts$DX, levels=c("CN","EMCI","LMCI","AD"))

p_mapt <- ggplot(mapt_counts, aes(x=DX, y=Count, fill=DX)) +
  geom_col(alpha=0.9, width=0.55) +
  geom_text(aes(label=Count), vjust=-0.5, size=7, color=DUKE_NAVY, fontface="bold") +
  scale_fill_manual(values=dx_colors, guide="none") +
  scale_y_continuous(limits=c(0, 28), expand=c(0,0)) +
  labs(title="MAPT-Interacting Proteins by Disease Stage",
       subtitle=sprintf("Proteins correlated with pTau217 AND interact with tau (OmniPath) | %d total unique", nrow(mapt_mat)),
       x="", y="Count") +
  theme_post

# Axis coverage
sum_tbl <- read.csv("output/tau_mechanism_by_DX/comparison_summary.csv")
axis_data <- data.frame()
for (i in 7:11) {
  axis_data <- rbind(axis_data, data.frame(
    Axis = c("GSK3B axis","MAPK cascade","CDK5 pathway","SRC/FYN/SYK","PPP phosphatase")[i-6],
    CN   = as.numeric(sum_tbl[i,"CN"]), EMCI = as.numeric(sum_tbl[i,"EMCI"]),
    LMCI = as.numeric(sum_tbl[i,"LMCI"]), AD = as.numeric(sum_tbl[i,"AD"])))
}
axis_long <- pivot_longer(axis_data, -Axis, names_to="DX", values_to="Coverage")
axis_long$DX <- factor(axis_long$DX, levels=c("CN","EMCI","LMCI","AD"))

p_axis <- ggplot(axis_long, aes(x=Axis, y=Coverage, fill=DX)) +
  geom_col(position="dodge", alpha=0.9, width=0.7) +
  geom_text(aes(label=sprintf("%.0f%%",Coverage)), position=position_dodge(0.7),
            vjust=-0.4, size=3.8, color=DUKE_NAVY) +
  scale_fill_manual(values=dx_colors) +
  scale_y_continuous(limits=c(0, 65), expand=c(0,0)) +
  labs(title="Tau Signaling Axis Coverage",
       subtitle="% curated axis genes significantly associated with pTau217 (FDR<0.05)",
       x="", y="Coverage (%)") +
  theme_post + theme(legend.position="bottom", axis.text.x=element_text(angle=25,hjust=1,size=13),
                      legend.text=element_text(size=11))

fig2 <- p_mapt / p_axis + plot_layout(heights=c(1, 1.2)) +
  plot_annotation(title="EMCI Is the Peak Activation Stage of the Tau Regulatory Network",
    subtitle="MAPT interactors and signaling axis coverage peak at EMCI-LMCI, decline in AD",
    theme=theme(plot.title=element_text(face="bold", size=20, color=DUKE_NAVY),
                plot.subtitle=element_text(size=13, color=DUKE_LIGHT)))
ggsave(file.path(out_dir, "Fig2_mapt_axis.png"), fig2, width=11, height=9, dpi=300)
cat("Fig2 done\n")

# ===========================================================================
cat("========== FIGURE 3: Enrichment + Network ==========\n")

go <- read.csv("output/pathway_enrichment_full/go_bp_enrichment.csv")
rx <- read.csv("output/pathway_enrichment_full/reactome_enrichment.csv")
kg <- read.csv("output/pathway_enrichment_full/kegg_enrichment.csv")
common <- c("Description","p.adjust","Count","label","Source")
mk <- function(x, src) {
  top <- head(x[order(x$p.adjust),], 10)
  top$label <- factor(stringr::str_wrap(top$Description, 45),
                      levels=rev(stringr::str_wrap(top$Description, 45)))
  top$Source <- src; top[, common]
}
ea <- mk(go,"GO BP")
ea$Source <- factor(ea$Source, levels="GO BP")

p_e <- ggplot(ea, aes(x=-log10(p.adjust), y=label, size=Count)) +
  geom_point(color=DUKE_BLUE, alpha=0.85) +
  scale_size_continuous(range=c(4, 10), name="Count") +
  labs(title="GO BP Enrichment of pTau-Associated Proteins",
       x=expression(-log[10](adjusted~P)), y="") +
  theme_minimal(base_size=16) +
  theme(plot.title=element_text(face="bold",size=20,color=DUKE_NAVY),
        legend.position="bottom", legend.text=element_text(size=14),
        legend.title=element_text(size=14),
        axis.text.y=element_text(size=18),
        axis.text.x=element_text(size=14),
        axis.title.x=element_text(size=17),
        panel.grid.minor=element_blank())

# Tau map
omni_url <- "https://omnipathdb.org/interactions?format=tsv&fields=sources,references&genesymbols=1"
tmp <- tempfile(fileext=".tsv"); download.file(omni_url, tmp, method="auto")
omni <- read.delim(tmp, stringsAsFactors=FALSE)
sig <- read.csv("output/tau_mechanism_analysis/gene_functional_classification.csv")
gene_rho <- sig %>% filter(EntrezGeneSymbol!="", !is.na(rho_median)) %>%
  group_by(EntrezGeneSymbol) %>% summarise(rho=median(rho_median,na.rm=TRUE), .groups="drop")

mapt_sig <- omni[(omni$source_genesymbol=="MAPT"|omni$target_genesymbol=="MAPT") &
                 (omni$source_genesymbol %in% gene_rho$EntrezGeneSymbol | omni$target_genesymbol %in% gene_rho$EntrezGeneSymbol), ]

tau_axes <- list(
  "GSK3B axis"=c("GSK3B","GSK3A","AKT1","AKT2","PTEN","PDPK1","PRKCA","PRKCB","PRKCG","PRKCZ","PPP2R1A","DYRK1A"),
  "MAPK cascade"=c("MAPK1","MAPK3","MAPK8","MAPK14","MAP2K1","MAP2K2","BRAF","RAF1","GRB2","HRAS","KRAS","RPS6KA3","RPS6KB1"),
  "CDK5 pathway"=c("CDK5","PPP3CA","PPP3R1","PPP5C","CAMK2A","ABL1","ABL2"),
  "SRC/FYN/SYK"=c("SRC","FYN","SYK","LYN","PTK2","PTK2B","PTPN11"),
  "PPP phosphatase"=c("PPP1CA","PPP1CB","PPP1CC","PPP2CA","PPP2CB","PPP2R5A","PPP3CB","PPP3CC"))

focus <- unique(c("MAPT", intersect(unlist(tau_axes), gene_rho$EntrezGeneSymbol)))
fe <- mapt_sig[mapt_sig$source_genesymbol %in% focus & mapt_sig$target_genesymbol %in% focus, ]
edges <- data.frame(from=fe$source_genesymbol, to=fe$target_genesymbol)
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
# Compute layout and get MAPT coordinates for manual nudge
layout_df <- create_layout(g_tau, layout="stress")
mapt_xy <- layout_df[layout_df$name=="MAPT", c("x","y")]
# Exclude MAPT from repel labels, add manually with offset
p_m <- ggraph(g_tau, layout="stress") +
  geom_edge_link(color="grey85", alpha=0.4) +
  geom_node_point(aes(size=degree, fill=axis), shape=21, color="grey40", stroke=0.3) +
  geom_node_text(aes(label=ifelse(name=="MAPT","",name), size=3),
                 repel=TRUE, max.overlaps=60, box.padding=0.3, color=DUKE_NAVY) +
  annotate("text", x=mapt_xy$x, y=mapt_xy$y - 0.3, label="MAPT",
           size=6, fontface="bold", color="black") +
  scale_fill_manual(values=ax_col, name="Signaling Axis") +
  scale_size_continuous(range=c(2,10), guide="none") +
  labs(title="Tau Signaling Map",
       subtitle=sprintf("%d nodes, %d edges | 5 axes converge on MAPT", vcount(g_tau), ecount(g_tau))) +
  theme_void() +
  theme(plot.title=element_text(face="bold",size=20,color=DUKE_NAVY,hjust=0.5),
        plot.subtitle=element_text(size=14,color=DUKE_LIGHT,hjust=0.5),
        plot.margin=margin(3,3,3,3),
        legend.position="bottom", legend.text=element_text(size=13),
        legend.title=element_text(size=14))

fig3 <- p_e + p_m + plot_layout(widths=c(1, 1.3)) +
  plot_annotation(title="Multi-Pathway Convergence on Tau Phosphorylation",
    subtitle="Left: GO Biological Process enrichment (ORA) | Right: MAPT-centered interaction network",
    theme=theme(plot.title=element_text(face="bold",size=22,color=DUKE_NAVY),
                plot.subtitle=element_text(size=15,color=DUKE_LIGHT)))
ggsave(file.path(out_dir, "Fig3_enrichment_network.png"), fig3, width=20, height=9, dpi=300)
cat("Fig3 done\n")

cat("\n========== ALL DONE ==========\n")
