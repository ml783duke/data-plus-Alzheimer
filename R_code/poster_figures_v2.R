###############################################################################
# Poster Figures v2 — 3 merged panels
###############################################################################
library(ggplot2); library(dplyr); library(tidyr); library(igraph); library(ggraph); library(patchwork)

out_dir <- "output_figure"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

DUKE_NAVY  <- "#003087"; DUKE_BLUE <- "#00539B"; DUKE_LIGHT <- "#7288A0"
DUKE_ORANGE <- "#C84E00"; DUKE_BG <- "#F0F2F5"

dx_colors <- c("CN" = DUKE_LIGHT, "EMCI" = DUKE_BLUE, "LMCI" = DUKE_NAVY, "AD" = DUKE_ORANGE)

theme_poster <- theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 13, color = DUKE_NAVY),
        plot.subtitle = element_text(size = 8, color = DUKE_LIGHT),
        plot.margin = margin(5, 8, 5, 8),
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
        axis.title = element_text(size = 10, color = DUKE_NAVY),
        axis.text = element_text(size = 8),
        strip.text = element_text(face = "bold", size = 9, color = DUKE_NAVY))

cat("========== Generating 3 Poster Figures ==========\n")

# ===========================================================================
# FIGURE 1: LOESS pTau217 by DX (left) + Lifestyle forest plot (right)
# ===========================================================================
cat("\n--- Fig 1: pTau trajectories + Lifestyle ---\n")

# Left panel: LOESS
csv <- read.csv("lifestyle_7var_ptau_proteomics.csv")
df_loess <- csv %>%
  mutate(ptau217_csf = as.numeric(ptau217_csf)) %>%
  filter(!is.na(ptau217_csf) & !is.na(age_at_csf_ptau) & !is.na(dx_entry),
         dx_entry %in% c("CN", "EMCI", "LMCI", "AD"))

MIN_VALID_PAIRS <- 10; LOESS_SPAN <- 0.75
loess_fits <- df_loess %>%
  group_by(dx_entry) %>%
  group_modify(~ {
    av <- .x$age_at_csf_ptau; pv <- .x$ptau217_csf
    vi <- !is.na(av) & !is.na(pv)
    if (sum(vi) < MIN_VALID_PAIRS) return(tibble::tibble())
    fd <- data.frame(x = av[vi], y = pv[vi])
    lo <- tryCatch(loess(y ~ x, span = LOESS_SPAN, degree = 1, data = fd), error = function(e) NULL)
    if (is.null(lo)) return(tibble::tibble())
    aseq <- seq(min(av, na.rm = TRUE), max(av, na.rm = TRUE), length.out = 100)
    tibble::tibble(age = aseq, fitted = predict(lo, data.frame(x = aseq)))
  }) %>% filter(n() > 0)

df_loess$dx_entry <- factor(df_loess$dx_entry, levels = c("CN", "EMCI", "LMCI", "AD"))

p_loess <- ggplot(df_loess, aes(x = age_at_csf_ptau, y = ptau217_csf, color = dx_entry)) +
  geom_point(alpha = 0.2, size = 0.8) +
  geom_line(data = loess_fits, aes(x = age, y = fitted), linewidth = 1.2) +
  scale_color_manual(values = dx_colors, name = "") +
  labs(title = "CSF pTau217 Age Trajectories", x = "Age (years)", y = "CSF pTau217") +
  theme_poster + theme(legend.position = c(0.12, 0.78),
                       legend.background = element_rect(fill = "white", color = "grey90"),
                       legend.key.size = unit(0.3, "cm"))

# Right panel: Lifestyle forest
lcsv <- read.csv("lifestyle_7var_ptau.csv")
master <- readxl::read_excel("master_data.xlsx", sheet = "Sheet1")
sex_lookup <- master[!duplicated(master$RID), c("RID", "PTGENDER")]
lcsv <- merge(lcsv, sex_lookup, by = "RID", all.x = TRUE)
lcsv$SEX <- ifelse(lcsv$PTGENDER == "Male", 1L, 0L)
df_a <- lcsv[!is.na(lcsv$ptau217_csf) & !is.na(lcsv$age_at_csf_ptau), ]
df_a$PTAU_log2 <- log2(df_a$ptau217_csf)

LIFESTYLE_VARS <- c("DHA", "EPA", "HCys", "NPIK", "NPIKTOT", "MH14ALCH", "MH16SMOK")
VAR_TYPE <- c("DHA"="cont","EPA"="cont","HCys"="cont","NPIK"="bin","NPIKTOT"="cont","MH14ALCH"="bin","MH16SMOK"="bin")
var_labels <- c("DHA"="DHA (n-3 PUFA)","EPA"="EPA (n-3 PUFA)","HCys"="Homocysteine",
  "NPIK"="NPI Sleep","NPIKTOT"="NPI Total","MH14ALCH"="Alcohol Abuse","MH16SMOK"="Smoking")

uni <- data.frame(variable = LIFESTYLE_VARS, n = NA_integer_, beta = NA_real_, se = NA_real_,
                  ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_)
for (i in seq_along(LIFESTYLE_VARS)) {
  v <- LIFESTYLE_VARS[i]; df_sub <- df_a[!is.na(df_a[[v]]), ]; uni$n[i] <- nrow(df_sub)
  x_val <- df_sub[[v]]
  x_model <- if (VAR_TYPE[v] == "cont" && sd(x_val, na.rm = TRUE) > 0) as.vector(scale(x_val)) else x_val
  fit <- lm(PTAU_log2 ~ x_model + age_at_csf_ptau, data = df_sub); s <- summary(fit)
  uni$beta[i] <- s$coefficients["x_model", "Estimate"]
  uni$se[i] <- s$coefficients["x_model", "Std. Error"]
  uni$p_value[i] <- s$coefficients["x_model", "Pr(>|t|)"]
  uni$ci_low[i] <- uni$beta[i] - 1.96 * uni$se[i]
  uni$ci_high[i] <- uni$beta[i] + 1.96 * uni$se[i]
}
uni$fdr <- p.adjust(uni$p_value, method = "BH")
uni$significant <- uni$fdr < 0.05
uni$label <- var_labels[uni$variable]
uni <- uni[order(uni$beta), ]
uni$label <- factor(uni$label, levels = uni$label)

p_forest <- ggplot(uni, aes(x = beta, y = label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  geom_point(aes(size = n), color = DUKE_BLUE, alpha = 0.9) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.2, color = DUKE_BLUE, linewidth = 0.7) +
  scale_size_continuous(name = "N", range = c(2, 4)) +
  labs(title = "Lifestyle Factors → CSF pTau217",
       subtitle = sprintf("log2(pTau) ~ variable + Age | N=%d | 0/7 FDR < 0.05", nrow(df_a)),
       x = "Standardized Effect", y = "") +
  theme_poster + theme(panel.grid.major.y = element_blank())

# Merge Fig 1
fig1 <- p_loess + p_forest + plot_layout(widths = c(1, 0.7)) +
  plot_annotation(title = "pTau217 Age Trajectories & Lifestyle Associations",
                  subtitle = "Left: LOESS (degree=1, span=0.75) by DX | Right: Univariate models adjusted for age",
                  theme = theme(plot.title = element_text(face = "bold", size = 14, color = DUKE_NAVY),
                                plot.subtitle = element_text(size = 8, color = DUKE_LIGHT)))

ggsave(file.path(out_dir, "Fig1_loess_lifestyle.png"), fig1, width = 14, height = 5.5, dpi = 300)
cat("Saved: Fig1\n")

# ===========================================================================
# FIGURE 2: MAPT interactors (top) + Tau axis coverage (bottom)
# ===========================================================================
cat("\n--- Fig 2: EMCI peak activation ---\n")

# Top: MAPT interactors
mapt_mat <- read.csv("output/tau_mechanism_by_DX/mapt_interactors_by_DX.csv")
mapt_counts <- data.frame(
  DX = c("CN","EMCI","LMCI","AD"),
  Count = c(sum(mapt_mat$CN), sum(mapt_mat$EMCI), sum(mapt_mat$LMCI), sum(mapt_mat$AD))
)
mapt_counts$DX <- factor(mapt_counts$DX, levels = c("CN","EMCI","LMCI","AD"))

p_mapt <- ggplot(mapt_counts, aes(x = DX, y = Count, fill = DX)) +
  geom_col(alpha = 0.9, width = 0.55) +
  geom_text(aes(label = Count), vjust = -0.5, size = 5, color = DUKE_NAVY, fontface = "bold") +
  scale_fill_manual(values = dx_colors, guide = "none") +
  scale_y_continuous(limits = c(0, 28), expand = c(0, 0)) +
  labs(title = "MAPT-Interacting Proteins by Disease Stage",
       subtitle = sprintf("Proteins correlated with pTau217 AND interact with tau (OmniPath) | %d total unique", nrow(mapt_mat)),
       x = "", y = "Count") +
  theme_poster

# Bottom: Axis coverage
sum_tbl <- read.csv("output/tau_mechanism_by_DX/comparison_summary.csv")
axis_data <- data.frame()
for (i in 7:11) {
  axis_data <- rbind(axis_data, data.frame(
    Axis = c("GSK3B axis","MAPK cascade","CDK5 pathway","SRC/FYN/SYK","PPP phosphatase")[i-6],
    CN = as.numeric(sum_tbl[i, "CN"]), EMCI = as.numeric(sum_tbl[i, "EMCI"]),
    LMCI = as.numeric(sum_tbl[i, "LMCI"]), AD = as.numeric(sum_tbl[i, "AD"])))
}
axis_long <- pivot_longer(axis_data, -Axis, names_to = "DX", values_to = "Coverage")
axis_long$DX <- factor(axis_long$DX, levels = c("CN","EMCI","LMCI","AD"))

p_axis <- ggplot(axis_long, aes(x = Axis, y = Coverage, fill = DX)) +
  geom_col(position = "dodge", alpha = 0.9, width = 0.7) +
  geom_text(aes(label = sprintf("%.0f%%", Coverage)),
            position = position_dodge(0.7), vjust = -0.4, size = 2.8, color = DUKE_NAVY) +
  scale_fill_manual(values = dx_colors, name = "") +
  scale_y_continuous(limits = c(0, 65), expand = c(0, 0)) +
  labs(title = "Tau Signaling Axis Coverage",
       subtitle = "% curated axis genes significantly associated with pTau217 (FDR < 0.05)",
       x = "", y = "Coverage (%)") +
  theme_poster + theme(axis.text.x = element_text(angle = 25, hjust = 1))

# Merge Fig 2
fig2 <- p_mapt / p_axis +
  plot_layout(heights = c(1, 1.2)) +
  plot_annotation(title = "EMCI Is the Peak Activation Stage of the Tau Regulatory Network",
                  subtitle = "MAPT interactors and signaling axis coverage both peak at EMCI-LMCI and decline in AD",
                  theme = theme(plot.title = element_text(face = "bold", size = 14, color = DUKE_NAVY),
                                plot.subtitle = element_text(size = 8, color = DUKE_LIGHT)))

ggsave(file.path(out_dir, "Fig2_mapt_axis.png"), fig2, width = 11, height = 9, dpi = 300)
cat("Saved: Fig2\n")

# ===========================================================================
# FIGURE 3: Enrichment (left) + Tau signaling map (right)
# ===========================================================================
cat("\n--- Fig 3: Pathway convergence & signaling map ---\n")

# Left: Enrichment dotplot
go <- read.csv("output/pathway_enrichment_full/go_bp_enrichment.csv")
rx <- read.csv("output/pathway_enrichment_full/reactome_enrichment.csv")
kg <- read.csv("output/pathway_enrichment_full/kegg_enrichment.csv")

make_plot <- function(x, src) {
  top <- head(x[order(x$p.adjust), ], 7)
  top$label <- factor(substr(top$Description, 1, 50), levels = rev(substr(top$Description, 1, 50)))
  top$Source <- src; top
}
ea <- rbind(make_plot(go, "GO BP"), make_plot(rx, "Reactome"), make_plot(kg, "KEGG"))
ea$Source <- factor(ea$Source, levels = c("GO BP", "Reactome", "KEGG"))

p_enrich <- ggplot(ea, aes(x = -log10(p.adjust), y = label, size = Count)) +
  geom_point(color = DUKE_BLUE, alpha = 0.85) +
  scale_size_continuous(range = c(2, 5.5), name = "Count") +
  facet_wrap(~ Source, scales = "free_y", ncol = 1) +
  labs(title = "Pathway Enrichment (ORA)", x = expression(-log[10](adjusted~P)), y = "") +
  theme_poster + theme(strip.text = element_text(size = 8))

# Right: Tau signaling map (simplified for poster space)
omni_url <- "https://omnipathdb.org/interactions?format=tsv&fields=sources,references&genesymbols=1"
tmp <- tempfile(fileext = ".tsv"); download.file(omni_url, tmp, method = "auto")
omni <- read.delim(tmp, stringsAsFactors = FALSE)

sig <- read.csv("output/tau_mechanism_analysis/gene_functional_classification.csv")
gene_rho <- sig %>% filter(EntrezGeneSymbol != "" & !is.na(rho_median)) %>%
  group_by(EntrezGeneSymbol) %>% summarise(rho = median(rho_median, na.rm = TRUE), .groups = "drop")

mapt_all <- omni[omni$source_genesymbol == "MAPT" | omni$target_genesymbol == "MAPT", ]
mapt_sig <- mapt_all[mapt_all$source_genesymbol %in% gene_rho$EntrezGeneSymbol |
                     mapt_all$target_genesymbol %in% gene_rho$EntrezGeneSymbol, ]

tau_axes_map <- list(
  "GSK3B axis" = c("GSK3B","GSK3A","AKT1","AKT2","PTEN","PDPK1","PRKCA","PRKCB","PRKCG","PRKCZ","PPP2R1A","DYRK1A"),
  "MAPK cascade" = c("MAPK1","MAPK3","MAPK8","MAPK14","MAP2K1","MAP2K2","BRAF","RAF1","GRB2","HRAS","KRAS","RPS6KA3","RPS6KB1"),
  "CDK5 pathway" = c("CDK5","PPP3CA","PPP3R1","PPP5C","CAMK2A","ABL1","ABL2"),
  "SRC/FYN/SYK" = c("SRC","FYN","SYK","LYN","PTK2","PTK2B","PTPN11"),
  "PPP phosphatase" = c("PPP1CA","PPP1CB","PPP1CC","PPP2CA","PPP2CB","PPP2R5A","PPP3CB","PPP3CC"))

focus_genes <- unique(c("MAPT", intersect(unlist(tau_axes_map), gene_rho$EntrezGeneSymbol)))
focus_edges <- mapt_sig[mapt_sig$source_genesymbol %in% focus_genes & mapt_sig$target_genesymbol %in% focus_genes, ]

edges <- data.frame(from = focus_edges$source_genesymbol, to = focus_edges$target_genesymbol)
nodes <- data.frame(gene = unique(c(edges$from, edges$to)))
nodes$rho <- gene_rho$rho[match(nodes$gene, gene_rho$EntrezGeneSymbol)]
nodes$axis <- "Other"
for (ax in names(tau_axes_map)) nodes$axis[nodes$gene %in% tau_axes_map[[ax]]] <- ax
nodes$axis[nodes$gene == "MAPT"] <- "MAPT"

g_tau <- graph_from_data_frame(edges, directed = FALSE, vertices = nodes)
V(g_tau)$degree <- degree(g_tau)

axis_colors <- c("MAPT" = "black", "GSK3B axis" = DUKE_BLUE, "MAPK cascade" = DUKE_ORANGE,
                 "CDK5 pathway" = DUKE_LIGHT, "SRC/FYN/SYK" = "#1B6B3A",
                 "PPP phosphatase" = "#6B3A8A", "Other" = "grey85")

set.seed(42)
p_map <- ggraph(g_tau, layout = "stress") +
  geom_edge_link(color = "grey85", alpha = 0.4) +
  geom_node_point(aes(size = degree, fill = axis), shape = 21, color = "grey40", stroke = 0.2) +
  geom_node_text(aes(label = name, size = ifelse(name == "MAPT", 4, 2.2)),
                 repel = TRUE, max.overlaps = 60, box.padding = 0.2, color = DUKE_NAVY) +
  scale_fill_manual(values = axis_colors, name = "Signaling Axis") +
  scale_size_continuous(range = c(1.5, 8), guide = "none") +
  labs(title = "Tau Signaling Map",
       subtitle = sprintf("%d nodes, %d edges | 5 axes converge on MAPT", vcount(g_tau), ecount(g_tau))) +
  theme_void() +
  theme(plot.title = element_text(face = "bold", size = 12, color = DUKE_NAVY, hjust = 0.5),
        plot.subtitle = element_text(size = 7, color = DUKE_LIGHT, hjust = 0.5),
        plot.margin = margin(2, 2, 2, 2),
        legend.position = "bottom", legend.text = element_text(size = 7))

# Merge Fig 3
fig3 <- p_enrich + p_map + plot_layout(widths = c(1, 1.3)) +
  plot_annotation(title = "Multi-Pathway Convergence on Tau Phosphorylation",
                  subtitle = "Left: Database-driven enrichment (GO+Reactome+KEGG) | Right: MAPT-centered interaction network, colored by signaling axis",
                  theme = theme(plot.title = element_text(face = "bold", size = 14, color = DUKE_NAVY),
                                plot.subtitle = element_text(size = 8, color = DUKE_LIGHT)))

ggsave(file.path(out_dir, "Fig3_enrichment_network.png"), fig3, width = 16, height = 7.5, dpi = 300)
cat("Saved: Fig3\n")

cat("\n========== Done ==========\n")
