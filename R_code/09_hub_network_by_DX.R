###############################################################################
# Generate per-DX hub protein interaction networks
###############################################################################

library(readxl)
library(dplyr)
library(igraph)
library(ggraph)
library(ggplot2)
library(patchwork)

set.seed(42)

out_dir <- "output/ptau_protein_pathway_by_DX"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load data ----
cat("Loading data...\n")
omni_url <- "https://omnipathdb.org/interactions?format=tsv&fields=sources,references&genesymbols=1"
tmp <- tempfile(fileext = ".tsv")
download.file(omni_url, tmp, method = "auto")
omni <- read.delim(tmp, stringsAsFactors = FALSE)

master <- read_excel("master_data.xlsx", sheet = "Sheet1")
dict <- read.csv("protein_raw_data/protein dict.csv", stringsAsFactors = FALSE)

df_bl <- master[master$VISCODE2 == "bl" & !is.na(master$PTAU) & !is.na(master$DX), ]
protein_cols <- grep("^X[0-9]+\\.[0-9]+$", colnames(df_bl), value = TRUE)

convert_protein <- function(x) { x[x == "NA" | x == ""] <- NA; as.numeric(x) }
df_bl[protein_cols] <- lapply(df_bl[protein_cols], convert_protein)

DX_GROUPS <- c("CN", "EMCI", "LMCI", "AD")
plot_list <- list()

for (g in DX_GROUPS) {
  cat(sprintf("\n========== DX=%s ==========\n", g))
  df_dx <- df_bl[df_bl$DX == g, ]
  if (nrow(df_dx) < 30) next

  # --- Correlation ---
  missing_rate <- sapply(df_dx[protein_cols], function(x) mean(is.na(x)) * 100)
  proteins_kept <- names(missing_rate[missing_rate <= 20])

  df_dx$PTAU_log2 <- log2(df_dx$PTAU)
  pm <- as.matrix(df_dx[proteins_kept])
  pl2 <- log2(pm)
  pl2[is.infinite(pl2) | is.nan(pl2)] <- NA

  ptau_vals <- df_dx$PTAU_log2
  results <- data.frame(protein_id = colnames(pl2), n = NA_integer_,
                        spearman_rho = NA_real_, p_value = NA_real_)
  for (i in seq_len(ncol(pl2))) {
    pv <- pl2[, i]; vi <- !is.na(pv); nv <- sum(vi)
    results$n[i] <- nv
    if (nv >= 10) {
      tst <- tryCatch(cor.test(pv[vi], ptau_vals[vi], method = "spearman", exact = FALSE),
                      error = function(e) NULL)
      if (!is.null(tst)) { results$spearman_rho[i] <- tst$estimate; results$p_value[i] <- tst$p.value }
    }
  }
  results <- results[!is.na(results$p_value), ]
  results$fdr <- p.adjust(results$p_value, method = "BH")
  results <- merge(results, dict[, c("Analytes", "EntrezGeneSymbol")],
                   by.x = "protein_id", by.y = "Analytes", all.x = TRUE)

  sig <- results[results$fdr < 0.05, ]
  sig_genes <- unique(sig$EntrezGeneSymbol[!is.na(sig$EntrezGeneSymbol) & sig$EntrezGeneSymbol != ""])
  cat(sprintf("Sig genes: %d\n", length(sig_genes)))

  gene_rho <- sig %>%
    filter(!is.na(EntrezGeneSymbol) & EntrezGeneSymbol != "") %>%
    group_by(EntrezGeneSymbol) %>%
    summarise(rho_median = median(spearman_rho, na.rm = TRUE), .groups = "drop")

  # --- Sig-sig interactions ---
  inter_sig <- omni[omni$source_genesymbol %in% sig_genes & omni$target_genesymbol %in% sig_genes, ]
  cat(sprintf("Sig-sig interactions: %d\n", nrow(inter_sig)))

  if (nrow(inter_sig) < 5) { cat("  Too few edges, skipping.\n"); next }

  # Top 25 hubs
  hc <- table(c(inter_sig$source_genesymbol, inter_sig$target_genesymbol))
  hub_df <- data.frame(gene = names(hc), n_int = as.integer(hc))
  hub_df <- hub_df[order(hub_df$n_int, decreasing = TRUE), ]
  top_hubs <- head(hub_df, 25)

  # Sub-network
  hub_edges <- inter_sig[
    inter_sig$source_genesymbol %in% top_hubs$gene &
    inter_sig$target_genesymbol %in% top_hubs$gene,
  ]
  hub_genes <- unique(c(hub_edges$source_genesymbol, hub_edges$target_genesymbol))

  hub_nodes <- data.frame(gene = hub_genes)
  hub_nodes$rho <- gene_rho$rho_median[match(hub_nodes$gene, gene_rho$EntrezGeneSymbol)]
  hub_nodes$n_int <- top_hubs$n_int[match(hub_nodes$gene, top_hubs$gene)]

  edges <- data.frame(from = hub_edges$source_genesymbol, to = hub_edges$target_genesymbol)
  gg <- graph_from_data_frame(edges, directed = FALSE, vertices = hub_nodes)
  V(gg)$degree <- degree(gg)

  cat(sprintf("Network: %d nodes, %d edges\n", length(hub_genes), nrow(edges)))

  # --- Individual plot ---
  p <- ggraph(gg, layout = "fr") +
    geom_edge_link(color = "grey85", alpha = 0.4) +
    geom_node_point(aes(size = degree, fill = rho), shape = 21, color = "grey40", stroke = 0.4) +
    geom_node_text(aes(label = name), size = 3, repel = TRUE, max.overlaps = 40, box.padding = 0.3) +
    scale_fill_gradient2(low = "#2c7bb6", mid = "white", high = "#d7191c",
                          midpoint = 0, na.value = "grey85", name = "PTAU rho") +
    scale_size_continuous(range = c(3, 10), name = "Connections") +
    labs(
      title = sprintf("Hub Protein Network — %s", g),
      subtitle = sprintf("Top 25 connected proteins | %d nodes, %d edges | %d subjects",
                         length(hub_genes), nrow(edges), nrow(df_dx)),
      caption = "OmniPath PPI | Color = PTAU rho | Size = degree"
    ) +
    theme_void() +
    theme(plot.title = element_text(face = "bold", size = 15),
          plot.subtitle = element_text(size = 9, color = "grey40"))

  ggsave(file.path(out_dir, sprintf("hub_network_%s.png", g)),
         p, width = 12, height = 10, dpi = 200)
  cat(sprintf("Saved: hub_network_%s.png\n", g))

  # --- For combined plot ---
  plot_list[[g]] <- ggraph(gg, layout = "fr") +
    geom_edge_link(color = "grey90", alpha = 0.3) +
    geom_node_point(aes(size = degree, fill = rho), shape = 21, color = "grey50", stroke = 0.3) +
    geom_node_text(aes(label = name), size = 2.5, repel = TRUE, max.overlaps = 40, box.padding = 0.2) +
    scale_fill_gradient2(low = "#2c7bb6", mid = "white", high = "#d7191c",
                          midpoint = 0, na.value = "grey85", guide = "none") +
    scale_size_continuous(range = c(2, 8), guide = "none") +
    labs(title = sprintf("%s (%d nodes)", g, length(hub_genes))) +
    theme_void() +
    theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
          plot.margin = margin(3, 3, 3, 3))
}

# ---- Combined 2x2 ----
cat("\nBuilding combined figure...\n")
combined <- wrap_plots(plot_list, ncol = 2) +
  plot_annotation(
    title = "Hub Protein Interaction Networks Across Disease Stages",
    subtitle = "Top 25 most-connected PTAU-significant proteins | OmniPath curated PPI",
    theme = theme(plot.title = element_text(face = "bold", size = 16),
                  plot.subtitle = element_text(size = 10, color = "grey40"))
  )
ggsave(file.path(out_dir, "hub_network_combined_2x2.png"),
       combined, width = 18, height = 16, dpi = 200)
cat("Saved: hub_network_combined_2x2.png\nDone!\n")
