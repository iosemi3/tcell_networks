## ============================================================================
## BACKBONE TF consensus centrality (main analysis script)
## Requires:
##  - objects: list of igraph/tbl_graph graphs, one per subtype (e.g. *_graph_obj)
##  - Output:
##     * vertexatrs (consensus centrality on intersection graph)
##     * degree_centrality_networks (TF x subtype degree matrix)
##     * degree heatmap (top 20 TFs by naive degree)
##     * eigenvector centrality stacked bar (top 25 TFs by max eigen)
## ============================================================================

library(igraph)
library(reshape2)
library(ggplot2)
library(ggsci)
library(viridis)
library(ggridges)
library(dplyr)

## ----------------------- Load graphs ----------------------------------------

objects <- readRDS(
  "/Users/sebas/Desktop/tcell network paper submission/rds_files/tbl_graphs.RDS"
)

## ----------------------- Intersection graph & consensus centrality ----------

# degree-based intersection across all subtype-specific graphs
graph_itersection <- do.call(igraph::graph.intersection, objects)

# vertex attributes (should include Pando centrality measures)
vertexatrs <- data.frame(igraph::vertex_attr(graph_itersection))
rownames(vertexatrs) <- vertexatrs$name

# keep only centrality columns and compute consensus metrics
vertexatrs <- vertexatrs[grep("centrality", colnames(vertexatrs))]
vertexatrs$median_centrality <- apply(vertexatrs, 1, median, na.rm = TRUE)
vertexatrs$mean_centrality   <- apply(vertexatrs, 1, mean,   na.rm = TRUE)
vertexatrs <- na.omit(vertexatrs)

common.tfs <- rownames(vertexatrs)

## ----------------------- Degree centrality matrix per network ---------------

## Build degree_centrality_networks from 'objects'
## rows: genes (TFs)
## cols: one column per network, named <name>_deg.cent
## where <name> is names(objects), e.g. naive_graph_obj → naive_graph_obj_deg.cent

# collect all TFs across all graphs
all_genes <- sort(unique(unlist(
  lapply(objects, function(obj) {
    g <- as.igraph(obj)
    names(V(g))
  })
)))

# initialize matrix
degree_centrality_networks <- matrix(
  NA_real_,
  nrow = length(all_genes),
  ncol = length(objects)
)
rownames(degree_centrality_networks) <- all_genes
colnames(degree_centrality_networks) <- paste0(names(objects), "_deg.cent")

# fill with degree centrality per network (unnormalized degree)
for (nm in names(objects)) {
  g   <- as.igraph(objects[[nm]])
  deg <- igraph::degree(g, v = V(g))
  
  col_name <- paste0(nm, "_deg.cent")  # e.g. "naive_graph_obj_deg.cent"
  degree_centrality_networks[names(deg), col_name] <- as.numeric(deg)
}

degree_centrality_networks <- as.data.frame(degree_centrality_networks)

## Optional: export degree matrix used downstream
# write.csv(
#   degree_centrality_networks,
#   "/corgi/martin/multiome_T/csv_files/Sebastian_network/merged_graph_obj_degree_centrality_from_tbl_graphs.csv"
# )

## ----------------------- Degree heatmap (top 20 TFs) ------------------------

## Choose the column that defines the ordering for TFs on the y-axis
## (here: naive network; adjust if naming differs)
order_col <- "naive_graph_obj_deg.cent"

deg_ord <- degree_centrality_networks %>%
  filter(!is.na(.data[[order_col]])) %>%
  arrange(desc(.data[[order_col]]))

deg_ord$TF <- rownames(deg_ord)

# melt to long format for plotting
deg_melt <- reshape2::melt(deg_ord, id.vars = "TF")

# order TF factor levels by decreasing degree in 'order_col'
deg_melt$TF <- factor(deg_melt$TF, levels = rev(unique(deg_ord$TF)))

# keep only top 20 TFs by naive degree
deg_melt_top <- deg_melt[deg_melt$TF %in% head(deg_ord$TF, 20), ]

ggplot(deg_melt_top, aes(x = variable, y = TF)) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradient(low = "blue", high = "red") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10, face = "bold")
  )

## ----------------------- Per-network degree export (optional) ---------------

for (i in names(objects)) {
  colname <- i
  g <- objects[[i]]
  degree_centrality <- igraph::degree(g, v = igraph::V(g))
  df <- data.frame(degree_centrality)
  print(head(df))
  # write.csv(
  #   df,
  #   file = paste0("/corgi/sebas/tcell_multi/backbone_degcent_all/",
  #                 colname, "_degree_centrality.csv")
  # )
}

## ----------------------- Eigenvector centrality across networks -------------


eigen_centrality <- matrix(
  NA_real_,
  nrow = length(common.tfs),
  ncol = length(objects)
)
rownames(eigen_centrality) <- common.tfs
colnames(eigen_centrality) <- names(objects)

for (i in names(objects)) {
  current_object <- objects[[i]]
  eig_vals <- igraph::eigen_centrality(current_object)$vector[common.tfs]
  eigen_centrality[, i] <- eig_vals
}

eigen_centrality <- as.data.frame(eigen_centrality)

eigen_centrality_melted <- eigen_centrality
eigen_centrality_melted$gene <- rownames(eigen_centrality)
eigen_centrality_melted <- melt(eigen_centrality_melted, id.vars = "gene")

## Select top 25 genes by maximum eigen-centrality across any subtype
genes25 <- eigen_centrality_melted %>%
  group_by(gene) %>%
  summarize(maxval = max(value, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(maxval)) %>%
  slice(1:25) %>%
  pull(gene)

eig_top <- eigen_centrality_melted %>%
  filter(gene %in% genes25)

ggplot(eig_top, aes(x = reorder(gene, value), y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  scale_fill_npg() +
  labs(x = "gene", y = "eigenvector centrality") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 20, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 20, face = "bold"),
    legend.title = element_text(size = 10, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
