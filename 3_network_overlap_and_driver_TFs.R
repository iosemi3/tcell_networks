#!/usr/bin/env Rscript
## ============================================================================
## 03_network_overlap_and_driver_TFs.R
##
## Standalone script:
##   - consumes Pando/xgboost graph objects (tbl_graphs.RDS)
##   - recomputes consensus vertex centrality (vertexatrs)
##   - computes global network overlap (Jaccard over gene sets)
##   - identifies:
##       * TOP TFs by consensus centrality
##       * subtype-specific non-backbone driver TFs (high degree, not in consensus)
##   - generates:
##       * global overlap heatmap
##       * driver TF heatmap (top 3 per subtype)
##
## This script complements:
##   01: graph-based consensus centralities
##   02: frozen-matrix backbone/activation/effector analysis
## ============================================================================

library(igraph)
library(reshape2)
library(ggplot2)
library(dplyr)
library(purrr)

## -------------------- Load graph objects ------------------------------------

objects <- readRDS(
  "/Users/sebas/Desktop/tcell network paper submission/rds_files/tbl_graphs.RDS"
)
# 'objects' is expected to be a named list of graph-like objects:
#   names(objects) ~ "<subtype>_graph_obj"

## -------------------- Recompute consensus centralities (vertexatrs) ---------

# intersection of all graphs
graph_itersection <- do.call(igraph::graph.intersection, objects)

# vertex attributes (should include centrality measures from Pando)
vertexatrs <- data.frame(igraph::vertex_attr(graph_itersection))
rownames(vertexatrs) <- vertexatrs$name

# keep only centrality columns and compute consensus metrics
vertexatrs <- vertexatrs[grep("centrality", colnames(vertexatrs))]
vertexatrs$median_centrality <- apply(vertexatrs, 1, median, na.rm = TRUE)
vertexatrs$mean_centrality   <- apply(vertexatrs, 1, mean,   na.rm = TRUE)
vertexatrs <- na.omit(vertexatrs)

# optional export of consensus table
# write.csv(vertexatrs,
#           "/corgi/sebas/tcell_multi/backbone_degcent_all/vertexatrs_consensus_centralities.csv")

## -------------------- Global network overlap (Jaccard) ----------------------

# Szymkiewicz–Simpson coefficient (if you want it later)
overlap_coefficient <- function(set1, set2) {
  intersection_size <- length(intersect(set1, set2))
  set1_size <- length(set1)
  set2_size <- length(set2)
  intersection_size / sqrt(set1_size * set2_size)
}

jaccard_index <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))
  union_size   <- length(union(set1, set2))
  intersection / union_size
}

n_networks <- length(objects)
global_overlap_networks <- matrix(NA_real_, nrow = n_networks, ncol = n_networks)

for (i in seq_len(n_networks)) {
  for (j in seq_len(n_networks)) {
    set1 <- names(V(objects[[i]]))
    set2 <- names(V(objects[[j]]))
    global_overlap_networks[i, j] <- jaccard_index(set1, set2)
  }
}

dimnames(global_overlap_networks) <- list(names(objects), names(objects))
print(global_overlap_networks)

# heatmap of global overlap
df <- melt(global_overlap_networks, varnames = c("Var1", "Var2"), value.name = "value")
df$Var1 <- gsub("_graph_obj", "", df$Var1)
df$Var2 <- gsub("_graph_obj", "", df$Var2)

ggplot(df, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "cornsilk", high = "brown3") +
  geom_text(aes(label = sprintf("%.2f", value)), size = 6, color = "black", fontface = "bold") +
  labs(
    x     = "Group",
    y     = "Group",
    title = "Global Overlap (Jaccard index)"
  ) +
  theme_minimal() +
  theme(
    legend.text   = element_text(face = "bold"),
    legend.title  = element_text(face = "bold"),
    legend.position = "none",
    panel.grid    = element_blank(),
    axis.text.x   = element_text(face = "bold", size = 15, angle = 45, hjust = 1),
    axis.text.y   = element_text(face = "bold", size = 15)
  )

## Optional: export overlap matrix
# write.csv(global_overlap_networks,
#           "/corgi/sebas/tcell_multi/backbone_degcent_all/global_overlap_networks_jaccard.csv")

## -------------------- TOP vs LOW consensus TFs ------------------------------

col_mean <- mean(vertexatrs$mean_centrality, na.rm = TRUE)
morethanmean_centrality <- subset(vertexatrs, mean_centrality >  col_mean)
lessthanmean_centrality <- subset(vertexatrs, mean_centrality <= col_mean)

# Optional exports:
# write.csv(morethanmean_centrality,
#           "/corgi/sebas/tcell_multi/backbone_degcent_all/top_consensus_TFs.csv")
# write.csv(lessthanmean_centrality,
#           "/corgi/sebas/tcell_multi/backbone_degcent_all/low_consensus_TFs.csv")

## -------------------- Subtype-specific non-backbone drivers -----------------
## For each network:
##   - remove consensus TFs (rows in vertexatrs)
##   - compute degree centrality
##   - rank TFs and keep all sorted → result_list
##   - then take top 3 per subtype for plotting

result_list <- list()

for (i in seq_along(objects)) {
  g <- objects[[i]]
  
  # vertex names
  vertex_names <- names(V(g))
  if (is.null(vertex_names)) {
    vertex_names <- as.character(seq_len(vcount(g)))
  }
  
  # remove consensus backbone/common TFs (present in vertexatrs)
  idx_remove <- which(vertex_names %in% rownames(vertexatrs))
  idx_keep   <- setdiff(seq_len(vcount(g)), idx_remove)
  
  if (length(idx_keep) == 0) {
    next
  }
  
  g_filtered <- induced_subgraph(g, vids = idx_keep)
  degree_centrality_filtered <- degree(g_filtered)
  
  result <- data.frame(
    Vertex          = names(degree_centrality_filtered),
    CentralityScore = as.numeric(degree_centrality_filtered),
    stringsAsFactors = FALSE
  )
  
  result <- result[order(result$CentralityScore, decreasing = TRUE), ]
  result_list[[i]] <- result
  names(result_list)[i] <- names(objects)[i]
}

# Optional: export each subtype's ordered driver list
# for (nm in names(result_list)) {
#   write.csv(
#     result_list[[nm]],
#     file = paste0("/corgi/sebas/tcell_multi/backbone_degcent_all/",
#                   nm, "_non_backbone_drivers.csv"),
#     row.names = FALSE
#   )
# }

## -------------------- Combine into single data.frame ------------------------

result_df <- bind_rows(
  lapply(names(result_list), function(nm) {
    df <- result_list[[nm]]
    df$Subtype <- nm
    df
  })
)

# Ensure consistent column naming
colnames(result_df) <- c("Vertex", "CentralityScore", "Subtype")

# Keep only top 3 drivers per subtype (as in your original comment)
df_list_top5 <- result_df %>%
  group_by(Subtype) %>%
  arrange(desc(CentralityScore)) %>%
  slice_head(n = 3) %>%
  ungroup()

# simplify subtype labels (drop "_graph_obj")
df_list_top5$Subtype <- gsub("_graph_obj", "", df_list_top5$Subtype)

## -------------------- Driver TF heatmap -------------------------------------

ggplot(df_list_top5, aes(x = Subtype, y = Vertex)) +
  geom_tile(aes(fill = CentralityScore), color = "black", size = 0.5) +
  scale_fill_gradient(low = "cornsilk", high = "brown3") +
  labs(
    y    = "Driver gene",
    x    = "Subtype",
    fill = "Degree centrality"
  ) +
  guides(fill = guide_legend(override.aes = list(size = 3))) +
  theme(
    panel.grid  = element_line(size = 1, linewidth = 0.5, colour = "black"),
    axis.text   = element_text(face = "bold", size = 14),
    axis.title  = element_text(face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(face = "bold"),
    legend.text  = element_text(face = "bold"),
    legend.key.size = unit(1, "lines")
  )

## Optional: export combined driver table
# write.csv(df_list_top5,
#           "/corgi/sebas/tcell_multi/backbone_degcent_all/non_backbone_top3_drivers_per_subtype.csv",
#           row.names = FALSE)
