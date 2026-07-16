library(ggplot2); library(dplyr); library(igraph); library(ggraph); library(patchwork)
DUKE_NAVY <- "#003087"; DUKE_BLUE <- "#00539B"; DUKE_LIGHT <- "#7288A0"; DUKE_ORANGE <- "#C84E00"
out_dir <- "output_figure"

# Enrichment
go <- read.csv("output/pathway_enrichment_full/go_bp_enrichment.csv")
rx <- read.csv("output/pathway_enrichment_full/reactome_enrichment.csv")
kg <- read.csv("output/pathway_enrichment_full/kegg_enrichment.csv")
make_plot <- function(x, src) {
  top <- head(x[order(x$p.adjust), ], 7)
  top$label <- factor(substr(top$Description, 1, 50), levels = rev(substr(top$Description, 1, 50)))
  top$Source <- src; top
}
common <- c("Description", "p.adjust", "Count", "label", "Source")
ea <- rbind(make_plot(go, "GO BP")[, common],
            make_plot(rx, "Reactome")[, common],
            make_plot(kg, "KEGG")[, common])
ea$Source <- factor(ea$Source, levels = c("GO BP", "Reactome", "KEGG"))

p_e <- ggplot(ea, aes(x = -log10(p.adjust), y = label, size = Count)) +
  geom_point(color = DUKE_BLUE, alpha = 0.85) +
  scale_size_continuous(range = c(2, 5.5), name = "Count") +
  facet_wrap(~ Source, scales = "free_y", ncol = 1) +
  labs(title = "Pathway Enrichment (ORA)", x = expression(-log[10](adjusted~P)), y = "") +
  theme_minimal(11) +
  theme(strip.text = element_text(face = "bold", size = 8),
        plot.title = element_text(face = "bold", color = DUKE_NAVY),
        legend.position = "bottom", panel.grid.minor = element_blank())

# Tau map
cat("Downloading OmniPath...\n")
omni_url <- "https://omnipathdb.org/interactions?format=tsv&fields=sources,references&genesymbols=1"
tmp <- tempfile(fileext = ".tsv"); download.file(omni_url, tmp, method = "auto")
omni <- read.delim(tmp, stringsAsFactors = FALSE)
cat(sprintf("OmniPath: %d rows\n", nrow(omni)))

sig <- read.csv("output/tau_mechanism_analysis/gene_functional_classification.csv")
gene_rho <- sig %>%
  filter(EntrezGeneSymbol != "" & !is.na(rho_median)) %>%
  group_by(EntrezGeneSymbol) %>%
  summarise(rho = median(rho_median, na.rm = TRUE), .groups = "drop")
cat(sprintf("Gene rho: %d genes\n", nrow(gene_rho)))

mapt_sig <- omni[
  (omni$source_genesymbol == "MAPT" | omni$target_genesymbol == "MAPT") &
  (omni$source_genesymbol %in% gene_rho$EntrezGeneSymbol |
   omni$target_genesymbol %in% gene_rho$EntrezGeneSymbol), ]
cat(sprintf("MAPT-sig edges: %d\n", nrow(mapt_sig)))

tau_axes <- list(
  "GSK3B axis" = c("GSK3B","GSK3A","AKT1","AKT2","PTEN","PDPK1","PRKCA","PRKCB","PRKCG","PRKCZ","PPP2R1A","DYRK1A"),
  "MAPK cascade" = c("MAPK1","MAPK3","MAPK8","MAPK14","MAP2K1","MAP2K2","BRAF","RAF1","GRB2","HRAS","KRAS","RPS6KA3","RPS6KB1"),
  "CDK5 pathway" = c("CDK5","PPP3CA","PPP3R1","PPP5C","CAMK2A","ABL1","ABL2"),
  "SRC/FYN/SYK" = c("SRC","FYN","SYK","LYN","PTK2","PTK2B","PTPN11"),
  "PPP phos" = c("PPP1CA","PPP1CB","PPP1CC","PPP2CA","PPP2CB","PPP2R5A","PPP3CB","PPP3CC"))

focus <- unique(c("MAPT", intersect(unlist(tau_axes), gene_rho$EntrezGeneSymbol)))
fe <- mapt_sig[mapt_sig$source_genesymbol %in% focus & mapt_sig$target_genesymbol %in% focus, ]
cat(sprintf("Focus edges: %d\n", nrow(fe)))

edges <- data.frame(from = fe$source_genesymbol, to = fe$target_genesymbol)
nodes <- data.frame(gene = unique(c(edges$from, edges$to)))
nodes$rho <- gene_rho$rho[match(nodes$gene, gene_rho$EntrezGeneSymbol)]
nodes$axis <- "Other"
for (ax in names(tau_axes)) nodes$axis[nodes$gene %in% tau_axes[[ax]]] <- ax
nodes$axis[nodes$gene == "MAPT"] <- "MAPT"
cat(sprintf("Nodes: %d, Axis distribution: %s\n", nrow(nodes),
            paste(names(table(nodes$axis)), table(nodes$axis), sep = "=", collapse = ", ")))

g_tau <- graph_from_data_frame(edges, directed = FALSE, vertices = nodes)
V(g_tau)$degree <- degree(g_tau)

ax_col <- c("MAPT" = "black", "GSK3B axis" = DUKE_BLUE, "MAPK cascade" = DUKE_ORANGE,
            "CDK5 pathway" = DUKE_LIGHT, "SRC/FYN/SYK" = "#1B6B3A",
            "PPP phos" = "#6B3A8A", "Other" = "grey85")

set.seed(42)
p_m <- ggraph(g_tau, layout = "stress") +
  geom_edge_link(color = "grey85", alpha = 0.4) +
  geom_node_point(aes(size = degree, fill = axis), shape = 21, color = "grey40", stroke = 0.2) +
  geom_node_text(aes(label = name, size = ifelse(name == "MAPT", 4, 2.2)),
                 repel = TRUE, max.overlaps = 60, box.padding = 0.2, color = DUKE_NAVY) +
  scale_fill_manual(values = ax_col, name = "Signaling Axis") +
  scale_size_continuous(range = c(1.5, 8), guide = "none") +
  labs(title = "Tau Signaling Map",
       subtitle = sprintf("%d nodes, %d edges | 5 axes converge on MAPT", vcount(g_tau), ecount(g_tau))) +
  theme_void() +
  theme(plot.title = element_text(face = "bold", size = 12, color = DUKE_NAVY, hjust = 0.5),
        plot.subtitle = element_text(size = 7, color = DUKE_LIGHT, hjust = 0.5),
        plot.margin = margin(2, 2, 2, 2),
        legend.position = "bottom", legend.text = element_text(size = 7))

cat("Building combined figure...\n")
fig3 <- p_e + p_m + plot_layout(widths = c(1, 1.3)) +
  plot_annotation(
    title = "Multi-Pathway Convergence on Tau Phosphorylation",
    subtitle = "Left: Database-driven enrichment (GO+Reactome+KEGG) | Right: MAPT-centered network, colored by axis",
    theme = theme(plot.title = element_text(face = "bold", size = 14, color = DUKE_NAVY),
                  plot.subtitle = element_text(size = 8, color = DUKE_LIGHT)))

ggsave(file.path(out_dir, "Fig3_enrichment_network.png"), fig3, width = 16, height = 7.5, dpi = 300)
cat("Saved: Fig3_enrichment_network.png\n")
