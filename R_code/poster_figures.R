###############################################################################
# Generate poster figures with Duke blue color scheme
# Duke palette: #003087 (navy), #00539B (blue), #7288A0 (light), #C84E00 (orange accent)
###############################################################################

library(ggplot2); library(dplyr); library(igraph); library(ggraph); library(patchwork)

out_dir <- "output_figure"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Duke color palette
DUKE_NAVY  <- "#003087"
DUKE_BLUE <- "#00539B"
DUKE_LIGHT <- "#7288A0"
DUKE_ORANGE <- "#C84E00"
DUKE_BG    <- "#F0F2F5"

# Poster theme: minimal, maximize data area
theme_poster <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, color = DUKE_NAVY),
    plot.subtitle = element_text(size = 9, color = DUKE_LIGHT),
    plot.margin = margin(10, 10, 10, 10),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    axis.title = element_text(size = 11, color = DUKE_NAVY),
    axis.text = element_text(size = 9)
  )

# Duke color gradient for 4 DX groups
dx_colors <- c("CN" = DUKE_LIGHT, "EMCI" = DUKE_BLUE, "LMCI" = DUKE_NAVY, "AD" = DUKE_ORANGE)
dx_colors_binary <- c("CN" = DUKE_BLUE, "non-CN" = DUKE_ORANGE,
                       "AD" = DUKE_ORANGE, "non-AD" = DUKE_BLUE)

cat("========== Generating Poster Figures ==========\n")

# ===========================================================================
# Fig 2A: GO BP + Reactome enrichment dotplots
# ===========================================================================
cat("\n--- Fig 2A: Pathway Enrichment ---\n")

# Load enrichment results
go <- read.csv("output/pathway_enrichment_full/go_bp_enrichment.csv")
rx <- read.csv("output/pathway_enrichment_full/reactome_enrichment.csv")

go_plot <- head(go[order(go$p.adjust), ], 10)
rx_plot <- head(rx[order(rx$p.adjust), ], 10)

go_plot$label <- substr(go_plot$Description, 1, 55)
go_plot$label <- factor(go_plot$label, levels = rev(go_plot$label))
go_plot$Source <- "GO BP"

rx_plot$label <- substr(rx_plot$Description, 1, 55)
rx_plot$label <- factor(rx_plot$label, levels = rev(rx_plot$label))
rx_plot$Source <- "Reactome"

enrich_all <- rbind(go_plot, rx_plot)
enrich_all$Source <- factor(enrich_all$Source, levels = c("GO BP", "Reactome"))

p_2a <- ggplot(enrich_all, aes(x = -log10(p.adjust), y = label, size = Count)) +
  geom_point(color = DUKE_BLUE, alpha = 0.85) +
  scale_size_continuous(range = c(3, 8)) +
  facet_wrap(~ Source, scales = "free_y", ncol = 1) +
  labs(title = "Pathway Enrichment of pTau-Associated Proteins",
       subtitle = "GO Biological Process + Reactome | ORA, FDR < 0.05 | Top 10 terms each",
       x = expression(-log[10](adjusted~P)), y = "") +
  theme_poster +
  theme(strip.text = element_text(face = "bold", size = 11, color = DUKE_NAVY))

ggsave(file.path(out_dir, "Fig2A_pathway_enrichment.png"), p_2a, width = 10, height = 7, dpi = 300)
cat("Saved: Fig2A\n")

# ===========================================================================
# Fig 2B: Tau signaling axis coverage by DX
# ===========================================================================
cat("\n--- Fig 2B: Tau Axis Coverage by DX ---\n")

sum_tbl <- read.csv("output/tau_mechanism_by_DX/comparison_summary.csv")
axis_data <- data.frame()
for (i in 7:11) {
  axis_data <- rbind(axis_data, data.frame(
    Axis = c("GSK3B axis","MAPK cascade","CDK5 pathway","SRC/FYN/SYK","PPP phosphatase")[i-6],
    CN = as.numeric(sum_tbl[i, "CN"]),
    EMCI = as.numeric(sum_tbl[i, "EMCI"]),
    LMCI = as.numeric(sum_tbl[i, "LMCI"]),
    AD = as.numeric(sum_tbl[i, "AD"])
  ))
}

axis_long <- tidyr::pivot_longer(axis_data, -Axis, names_to = "DX", values_to = "Coverage")
axis_long$DX <- factor(axis_long$DX, levels = c("CN","EMCI","LMCI","AD"))

p_2b <- ggplot(axis_long, aes(x = Axis, y = Coverage, fill = DX)) +
  geom_col(position = "dodge", alpha = 0.9, width = 0.7) +
  geom_text(aes(label = sprintf("%.0f%%", Coverage)),
            position = position_dodge(0.7), vjust = -0.4, size = 3, color = DUKE_NAVY) +
  scale_fill_manual(values = dx_colors, name = "") +
  scale_y_continuous(limits = c(0, 85), expand = c(0, 0)) +
  labs(title = "Tau Signaling Axis Coverage by Disease Stage",
       subtitle = "% curated axis genes significantly associated with pTau217 (FDR < 0.05)",
       x = "", y = "Coverage (%)") +
  theme_poster +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 10))

ggsave(file.path(out_dir, "Fig2B_tau_axis_by_DX.png"), p_2b, width = 10, height = 5.5, dpi = 300)
cat("Saved: Fig2B\n")

# ===========================================================================
# Fig 3A: MAPT interactors bar chart by DX
# ===========================================================================
cat("\n--- Fig 3A: MAPT Interactors by DX ---\n")

mapt_mat <- read.csv("output/tau_mechanism_by_DX/mapt_interactors_by_DX.csv")

# Count per DX
mapt_counts <- data.frame(
  DX = c("CN","EMCI","LMCI","AD"),
  Count = c(sum(mapt_mat$CN), sum(mapt_mat$EMCI), sum(mapt_mat$LMCI), sum(mapt_mat$AD))
)
mapt_counts$DX <- factor(mapt_counts$DX, levels = c("CN","EMCI","LMCI","AD"))

p_3a <- ggplot(mapt_counts, aes(x = DX, y = Count, fill = DX)) +
  geom_col(alpha = 0.9, width = 0.6) +
  geom_text(aes(label = Count), vjust = -0.5, size = 6, color = DUKE_NAVY, fontface = "bold") +
  scale_fill_manual(values = dx_colors, guide = "none") +
  scale_y_continuous(limits = c(0, max(mapt_counts$Count) * 1.15), expand = c(0, 0)) +
  labs(title = "MAPT-Interacting Proteins by Disease Stage",
       subtitle = sprintf("Proteins significantly correlated with pTau217 AND interact with tau (OmniPath) | %d total unique",
                          nrow(mapt_mat)),
       x = "", y = "Number of MAPT-Interacting Proteins") +
  theme_poster

ggsave(file.path(out_dir, "Fig3A_mapt_by_DX.png"), p_3a, width = 7, height = 5.5, dpi = 300)
cat("Saved: Fig3A\n")

# ===========================================================================
# Fig 3B: Tau signaling map (MAPT-centered network)
# ===========================================================================
cat("\n--- Fig 3B: Tau Signaling Map ---\n")

# Rebuild the tau signaling network with Duke colors
omni_url <- "https://omnipathdb.org/interactions?format=tsv&fields=sources,references&genesymbols=1"
tmp <- tempfile(fileext = ".tsv")
download.file(omni_url, tmp, method = "auto")
omni <- read.delim(tmp, stringsAsFactors = FALSE)

# Load gene rho
sig <- read.csv("output/tau_mechanism_analysis/gene_functional_classification.csv")
gene_rho <- sig %>% filter(EntrezGeneSymbol != "" & !is.na(rho_median)) %>%
  group_by(EntrezGeneSymbol) %>% summarise(rho = median(rho_median, na.rm=TRUE), .groups="drop")

# MAPT edges
mapt_all <- omni[omni$source_genesymbol == "MAPT" | omni$target_genesymbol == "MAPT", ]
mapt_sig <- mapt_all[mapt_all$source_genesymbol %in% gene_rho$EntrezGeneSymbol |
                     mapt_all$target_genesymbol %in% gene_rho$EntrezGeneSymbol, ]

tau_axes_map <- list(
  "GSK3B axis" = c("GSK3B","GSK3A","AKT1","AKT2","PTEN","PDPK1","PIK3CA","PIK3R1","PRKCA","PRKCB","PRKCG","PRKCZ","PPP2R1A","DYRK1A"),
  "MAPK cascade" = c("MAPK1","MAPK3","MAPK8","MAPK14","MAP2K1","MAP2K2","BRAF","RAF1","GRB2","HRAS","KRAS","RPS6KA3","RPS6KB1"),
  "CDK5 pathway" = c("CDK5","PPP3CA","PPP3R1","PPP5C","CAMK2A","ABL1","ABL2"),
  "SRC/FYN/SYK" = c("SRC","FYN","SYK","LYN","HCK","PTK2","PTK2B","PTPN11"),
  "PPP phosphatases" = c("PPP1CA","PPP1CB","PPP1CC","PPP2CA","PPP2CB","PPP2R5A","PPP3CB","PPP3CC")
)

all_axis_genes <- unique(unlist(tau_axes_map))
focus_genes <- unique(c("MAPT", intersect(all_axis_genes, gene_rho$EntrezGeneSymbol)))

focus_edges <- mapt_sig[
  (mapt_sig$source_genesymbol %in% focus_genes & mapt_sig$target_genesymbol %in% focus_genes) |
  mapt_sig$source_genesymbol == "MAPT" | mapt_sig$target_genesymbol == "MAPT",
]

edges <- data.frame(from = focus_edges$source_genesymbol, to = focus_edges$target_genesymbol)
nodes <- data.frame(gene = unique(c(edges$from, edges$to)))
nodes$rho <- gene_rho$rho[match(nodes$gene, gene_rho$EntrezGeneSymbol)]
nodes$axis <- "Other"
for (ax_name in names(tau_axes_map)) nodes$axis[nodes$gene %in% tau_axes_map[[ax_name]]] <- ax_name
nodes$axis[nodes$gene == "MAPT"] <- "MAPT"

g_tau <- graph_from_data_frame(edges, directed = FALSE, vertices = nodes)
V(g_tau)$degree <- degree(g_tau)

axis_colors <- c("MAPT" = "black", "GSK3B axis" = DUKE_BLUE,
                 "MAPK cascade" = DUKE_ORANGE, "CDK5 pathway" = DUKE_LIGHT,
                 "SRC/FYN/SYK" = "#1B6B3A", "PPP phosphatases" = "#6B3A8A",
                 "Other" = "grey85")

set.seed(123)
p_3b <- ggraph(g_tau, layout = "stress") +
  geom_edge_link(color = "grey85", alpha = 0.5) +
  geom_node_point(aes(size = degree, fill = axis), shape = 21, color = "grey40", stroke = 0.3) +
  geom_node_text(aes(label = name, size = ifelse(name == "MAPT", 4.5, 2.5)),
                 repel = TRUE, max.overlaps = 60, box.padding = 0.3, color = DUKE_NAVY) +
  scale_fill_manual(values = axis_colors, name = "Signaling Axis") +
  scale_size_continuous(range = c(2, 10), guide = "none") +
  labs(title = "Tau Signaling Map: Multi-Axis Convergence on MAPT",
       subtitle = sprintf("MAPT-interacting proteins in pTau-significant set | %d nodes, %d edges",
                          vcount(g_tau), ecount(g_tau))) +
  theme_void() +
  theme(plot.title = element_text(face = "bold", size = 14, color = DUKE_NAVY, hjust = 0.5),
        plot.subtitle = element_text(size = 9, color = DUKE_LIGHT, hjust = 0.5),
        plot.margin = margin(5, 5, 5, 5),
        legend.position = "bottom")

ggsave(file.path(out_dir, "Fig3B_tau_signaling_map.png"), p_3b, width = 14, height = 12, dpi = 300)
cat("Saved: Fig3B\n")

# ===========================================================================
# Fig 4: Functional categories + candidates
# ===========================================================================
cat("\n--- Fig 4: Functional Categories + Candidates ---\n")

cat_sum <- read.csv("output/tau_mechanism_analysis/functional_category_summary.csv")
cat_sum <- cat_sum[cat_sum$Category1 != "Others", ]
cat_sum$Category1 <- factor(cat_sum$Category1, levels = cat_sum$Category1[order(cat_sum$N)])

p_4a <- ggplot(cat_sum, aes(x = N, y = Category1, fill = Avg_Rho)) +
  geom_col(alpha = 0.9, width = 0.7) +
  geom_text(aes(label = N), hjust = -0.3, size = 3.5, color = DUKE_NAVY) +
  scale_fill_gradient2(low = DUKE_BLUE, mid = "white", high = DUKE_ORANGE, midpoint = 0) +
  scale_x_continuous(limits = c(0, max(cat_sum$N) * 1.2), expand = c(0, 0)) +
  labs(title = "Functional Classification of pTau-Associated Proteins",
       subtitle = sprintf("%d genes | Mean PTAU rho shown by fill color",
                          sum(cat_sum$N)),
       x = "Number of Genes", y = "", fill = "Mean rho") +
  theme_poster

ggsave(file.path(out_dir, "Fig4A_functional_categories.png"), p_4a, width = 9, height = 5, dpi = 300)
cat("Saved: Fig4A\n")

cat("\n========== Done ==========\n")
cat(sprintf("Output: %s\n", normalizePath(out_dir)))

# ===========================================================================
# Fig 1B: LOESS pTau217 trend by DX
# ===========================================================================
cat("\n--- Fig 1B: LOESS pTau217 by DX ---\n")

csv <- read.csv("lifestyle_7var_ptau_proteomics.csv")
df_loess <- csv %>%
  mutate(ptau217_csf = as.numeric(ptau217_csf)) %>%
  filter(!is.na(ptau217_csf) & !is.na(age_at_csf_ptau) & !is.na(dx_entry),
         dx_entry %in% c("CN", "EMCI", "LMCI", "AD"))

# Fit LOESS per DX
MIN_VALID_PAIRS <- 10
LOESS_SPAN <- 0.75

loess_fits <- df_loess %>%
  group_by(dx_entry) %>%
  group_modify(~ {
    age_vals <- .x$age_at_csf_ptau; ptau_vals <- .x$ptau217_csf
    valid_idx <- !is.na(age_vals) & !is.na(ptau_vals)
    if (sum(valid_idx) < MIN_VALID_PAIRS) return(tibble::tibble())
    fit_df <- data.frame(x = age_vals[valid_idx], y = ptau_vals[valid_idx])
    lo <- tryCatch(loess(y ~ x, span = LOESS_SPAN, degree = 1, data = fit_df),
                   error = function(e) NULL)
    if (is.null(lo)) return(tibble::tibble())
    age_seq <- seq(min(age_vals, na.rm = TRUE), max(age_vals, na.rm = TRUE), length.out = 100)
    tibble::tibble(age = age_seq, fitted = predict(lo, data.frame(x = age_seq)))
  }) %>% filter(n() > 0)

df_loess$dx_entry <- factor(df_loess$dx_entry, levels = c("CN", "EMCI", "LMCI", "AD"))

p_1b <- ggplot(df_loess, aes(x = age_at_csf_ptau, y = ptau217_csf, color = dx_entry)) +
  geom_point(alpha = 0.25, size = 1.2) +
  geom_line(data = loess_fits, aes(x = age, y = fitted, group = dx_entry, color = dx_entry),
            linewidth = 1.4) +
  scale_color_manual(values = dx_colors, name = "") +
  labs(
    title = expression(paste("CSF pTau217 Age Trajectories by Disease Stage")),
    subtitle = sprintf("LOESS (degree=1, span=%.2f) | %d subjects | CN+EMCI+LMCI+AD",
                       LOESS_SPAN, nrow(df_loess)),
    x = "Age at CSF Collection (years)",
    y = "CSF pTau217 (pg/mL)"
  ) +
  theme_poster +
  theme(legend.position = c(0.12, 0.75),
        legend.background = element_rect(fill = "white", color = "grey90"),
        legend.key.size = unit(0.4, "cm"))

ggsave(file.path(out_dir, "Fig1B_ptau217_loess_by_DX.png"), p_1b, width = 9, height = 6, dpi = 300)
cat("Saved: Fig1B\n")

# ===========================================================================
# Fig 5: Lifestyle forest plot (CSF pTau217, univariate)
# ===========================================================================
cat("\n--- Fig 5: Lifestyle Forest Plot ---\n")

# Re-create with Duke colors using the same analysis from script 11
csv <- read.csv("lifestyle_7var_ptau.csv")
master <- readxl::read_excel("master_data.xlsx", sheet = "Sheet1")
sex_lookup <- master[!duplicated(master$RID), c("RID", "PTGENDER")]
csv <- merge(csv, sex_lookup, by = "RID", all.x = TRUE)
csv$SEX <- ifelse(csv$PTGENDER == "Male", 1L, 0L)

df_a <- csv[!is.na(csv$ptau217_csf) & !is.na(csv$age_at_csf_ptau), ]
df_a$PTAU_log2 <- log2(df_a$ptau217_csf)

LIFESTYLE_VARS <- c("DHA", "EPA", "HCys", "NPIK", "NPIKTOT", "MH14ALCH", "MH16SMOK")
VAR_TYPE <- c("DHA"="cont","EPA"="cont","HCys"="cont","NPIK"="bin","NPIKTOT"="cont","MH14ALCH"="bin","MH16SMOK"="bin")
var_labels <- c("DHA"="DHA (n-3 PUFA)","EPA"="EPA (n-3 PUFA)","HCys"="Homocysteine",
  "NPIK"="NPI Sleep Domain","NPIKTOT"="NPI Total Score","MH14ALCH"="Alcohol Abuse History","MH16SMOK"="Smoking Status")

uni <- data.frame(variable=LIFESTYLE_VARS, n=NA_integer_, beta=NA_real_, se=NA_real_,
                  ci_low=NA_real_, ci_high=NA_real_, p_value=NA_real_)
for (i in seq_along(LIFESTYLE_VARS)) {
  v <- LIFESTYLE_VARS[i]
  df_sub <- df_a[!is.na(df_a[[v]]), ]
  uni$n[i] <- nrow(df_sub)
  x_val <- df_sub[[v]]
  x_model <- if (VAR_TYPE[v]=="cont" && sd(x_val,na.rm=TRUE)>0) as.vector(scale(x_val)) else x_val
  fit <- lm(PTAU_log2 ~ x_model + age_at_csf_ptau, data = df_sub)
  s <- summary(fit)
  uni$beta[i] <- s$coefficients["x_model","Estimate"]
  uni$se[i] <- s$coefficients["x_model","Std. Error"]
  uni$p_value[i] <- s$coefficients["x_model","Pr(>|t|)"]
  uni$ci_low[i] <- uni$beta[i] - 1.96 * uni$se[i]
  uni$ci_high[i] <- uni$beta[i] + 1.96 * uni$se[i]
}
uni$fdr <- p.adjust(uni$p_value, method="BH")
uni$significant <- uni$fdr < 0.05
uni$label <- var_labels[uni$variable]
uni <- uni[order(uni$beta), ]
uni$label <- factor(uni$label, levels = uni$label)
uni$sig_anno <- ifelse(uni$significant, "*", "")

p_5 <- ggplot(uni, aes(x = beta, y = label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.5) +
  geom_point(aes(size = n), color = DUKE_BLUE, alpha = 0.9) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.2, color = DUKE_BLUE, alpha = 0.8, linewidth = 0.8) +
  geom_text(aes(x = ifelse(beta > 0, ci_high + 0.01, ci_low - 0.01), label = sig_anno),
            size = 6, color = DUKE_ORANGE, vjust = 0.5) +
  scale_size_continuous(name = "N", range = c(2.5, 5)) +
  labs(
    title = "Lifestyle Factors: Association with CSF pTau217",
    subtitle = sprintf("Univariate: log2(pTau217) ~ variable + Age | %d subjects | * FDR < 0.05 | None significant",
                       nrow(df_a)),
    x = "Standardized Effect (per 1-SD or 1 vs 0)", y = ""
  ) +
  theme_poster +
  theme(panel.grid.major.y = element_blank())

ggsave(file.path(out_dir, "Fig5_lifestyle_forest.png"), p_5, width = 9, height = 4.5, dpi = 300)
cat("Saved: Fig5\n")

# ===========================================================================
# Fig 2A (FIXED): GO + Reactome + KEGG enrichment, highlight key terms
# ===========================================================================
cat("\n--- Fig 2A (fixed): 3-database enrichment ---\n")

go <- read.csv("output/pathway_enrichment_full/go_bp_enrichment.csv")
rx <- read.csv("output/pathway_enrichment_full/reactome_enrichment.csv")
kg <- read.csv("output/pathway_enrichment_full/kegg_enrichment.csv")

go_plot <- head(go[order(go$p.adjust), ], 8)
rx_plot <- head(rx[order(rx$p.adjust), ], 8)
kg_plot <- head(kg[order(kg$p.adjust), ], 8)

go_plot$label <- substr(go_plot$Description, 1, 55)
go_plot$label <- factor(go_plot$label, levels = rev(go_plot$label))
go_plot$Source <- "GO BP"

rx_plot$label <- substr(rx_plot$Description, 1, 55)
rx_plot$label <- factor(rx_plot$label, levels = rev(rx_plot$label))
rx_plot$Source <- "Reactome"

kg_plot$label <- substr(kg_plot$Description, 1, 55)
kg_plot$label <- factor(kg_plot$label, levels = rev(kg_plot$label))
kg_plot$Source <- "KEGG"

enrich_all <- rbind(go_plot, rx_plot, kg_plot)
enrich_all$Source <- factor(enrich_all$Source, levels = c("GO BP", "Reactome", "KEGG"))

p_2a_fix <- ggplot(enrich_all, aes(x = -log10(p.adjust), y = label, size = Count)) +
  geom_point(color = DUKE_BLUE, alpha = 0.85) +
  scale_size_continuous(range = c(2.5, 7), name = "Gene Count") +
  facet_wrap(~ Source, scales = "free_y", ncol = 1) +
  labs(title = "Pathway Enrichment of pTau-Associated Proteins",
       subtitle = expression(paste("GO BP + Reactome + KEGG | ORA, FDR < 0.05 | Top 8 terms each | ", italic("axonogenesis"), " (P=1e-8), neurotrophin signaling (P=3e-3), RTK signaling (P=3e-3)")),
       x = expression(-log[10](adjusted~P)), y = "") +
  theme_poster +
  theme(strip.text = element_text(face = "bold", size = 10, color = DUKE_NAVY))

ggsave(file.path(out_dir, "Fig2A_pathway_enrichment.png"), p_2a_fix, width = 11, height = 8, dpi = 300)
cat("Saved: Fig2A (fixed)\n")
