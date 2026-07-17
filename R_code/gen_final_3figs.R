###############################################################################
# Final 3 poster figures: Lifestyle + MAPT peak + Tau network
###############################################################################
library(ggplot2); library(dplyr); library(tidyr); library(igraph); library(ggraph)

out_dir <- "output_figure"
DUKE_NAVY  <- "#003087"; DUKE_BLUE <- "#00539B"; DUKE_LIGHT <- "#7288A0"; DUKE_ORANGE <- "#C84E00"
# DX color palette: coordinated Duke-adjacent, clearly distinguishable across 4 groups
# CN = grey-blue | EMCI = Duke Navy | LMCI = dark green | AD = Duke Orange
dx_colors  <- c("CN"="#7288A0","EMCI"="#003087","LMCI"="#1B6B3A","AD"="#C84E00")

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

# Subtitle styling: find the subtitle's y-position relative to plot area
p2_base <- ggplot(mapt_counts, aes(x = DX, y = Count, fill = DX)) +
  geom_col(alpha = 0.9, width = 0.55) +
  geom_text(aes(label = Count), vjust = -0.5, size = 10, color = DUKE_NAVY, fontface = "bold") +
  scale_fill_manual(values = dx_colors, guide = "none") +
  scale_y_continuous(limits = c(0, 28), expand = c(0, 0)) +
  labs(title = "Tau-Interacting Proteins Peak at EMCI",
       subtitle = "", x = "", y = "Number of MAPT-Interacting Proteins") +
  theme_post

# Add boxed MAPT → pTau217 as annotation (replaces subtitle)
p2 <- p2_base +
  annotate("label", x = 2.05, y = 26.2, label = "MAPT",
           fill = "white", color = DUKE_NAVY, size = 4.5, fontface = "bold",
           label.padding = unit(0.2, "lines")) +
  annotate("text",  x = 2.20, y = 26.2, label = "→ pTau217",
           size = 5, color = DUKE_LIGHT, hjust = 0, fontface = "italic")

ggsave(file.path(out_dir,"Fig2_mapt_EMCI_peak.png"), p2, width=10, height=8, dpi=300)
cat("Fig2 done\n")

cat("\n=== FIGURE 3: GO BP Enrichment by DX (Top 5 per group) ===\n")

# Load enrichment data
enrich <- read.csv("output/pathway_enrichment_by_DX/enrichment_by_DX.csv")
go <- enrich[enrich$Source == "GO_BP", ]

# Top 5 per DX group -> union of all terms
go_top5 <- go %>%
  group_by(DX) %>%
  slice_min(p.adjust, n = 5) %>%
  ungroup()

# Union of top 5 terms across all DX groups
terms_all <- unique(go_top5$Description)
cat(sprintf("Union of top 5 terms: %d unique pathways\n", length(terms_all)))

# Keep only those terms from the full dataset
df <- go[go$Description %in% terms_all, ]
df$DX <- factor(df$DX, levels = c("CN", "EMCI", "LMCI", "AD"))

# Flag whether this term is in this DX's own top 5
df$in_top5 <- paste(df$Description, df$DX) %in% paste(go_top5$Description, go_top5$DX)

# Sort terms: 4-group terms first, then partial
term_groups <- table(go_top5$Description)
term_order <- names(sort(term_groups, decreasing = TRUE))
df$Description <- factor(df$Description, levels = rev(term_order))

# Short labels
make_label <- function(d) {
  ifelse(d == "neuron projection morphogenesis", "Neuron projection\nmorphogenesis",
  ifelse(d == "plasma membrane bounded cell projection morphogenesis", "Plasma membrane bounded\ncell projection morphogenesis",
  ifelse(d == "cell projection morphogenesis", "Cell projection\nmorphogenesis",
  ifelse(d == "axonogenesis", "Axonogenesis",
  ifelse(d == "cell morphogenesis involved in neuron differentiation", "Cell morphogenesis involved\nin neuron differentiation",
  ifelse(d == "axon development", "Axon development", d))))))
}
df$label <- factor(make_label(as.character(df$Description)),
                   levels = rev(make_label(term_order)))

# DX labels above first row
first_term <- names(sort(term_groups, decreasing = TRUE))[1]
df_first <- df[df$Description == first_term & df$in_top5, ]
df_first <- df_first[match(c("CN","EMCI","LMCI","AD"), df_first$DX), ]

# Only plot points that are in each group's top 5
df_plot <- df[df$in_top5, ]

p3 <- ggplot(df_plot, aes(x = RichFactor, y = label, color = DX)) +
  geom_point(size = 5, alpha = 0.9) +
  geom_text(data = df_first,
            aes(x = RichFactor, y = length(term_order) + 0.35, label = DX, color = DX),
            size = 5.5, fontface = "bold", vjust = 0) +
  scale_color_manual(values = dx_colors, guide = "none") +
  scale_x_continuous(limits = c(0.29, 0.74), expand = c(0, 0),
                     labels = scales::label_number(accuracy = 0.05)) +
  labs(title = "GO Biological Process Enrichment by Disease Stage",
       subtitle = "",
       x = "RichFactor", y = "") +
  theme_post

ggsave(file.path(out_dir, "Fig3_tau_signaling_map.png"), p3, width = 13, height = 8.5, dpi = 300)
cat("Fig3 done\n\nALL DONE\n")
