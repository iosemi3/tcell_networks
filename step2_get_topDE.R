library(Seurat)
library(dplyr)
library(ggplot2)

# Load your Seurat object
load(snakemake@input[["seurat_object"]])

# Ensure that the loaded object has the expected name
if (!exists("pbmc_1103")) {
  stop("The Seurat object 'pbmc_1103' was not found in the loaded RData file.")
}

# Perform differential expression analysis
all_markers <- FindAllMarkers(pbmc_1103, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Get the top 10 markers for each cluster
top10_markers <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

# Save to CSV
write.csv(top10_markers, file = snakemake@output[["cluster_markers"]], row.names = FALSE)

# Plot the top 10 markers for each cluster
pdf("results/top10_markers_per_cluster.pdf")
for (cluster in unique(top10_markers$cluster)) {
  cluster_markers <- top10_markers %>% filter(cluster == !!cluster)
  plot <- VlnPlot(pbmc_1103, features = cluster_markers$gene, pt.size = 0.1) + ggtitle(paste("Cluster", cluster))
  print(plot)
}
dev.off()

# Perform the analysis based on specific identities instead of clusters
Idents(pbmc_1103) <- pbmc_1103$labels.DICE
all_markers_identity <- FindAllMarkers(pbmc_1103, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Get the top 10 markers for each identity
top10_markers_identity <- all_markers_identity %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

# Save to CSV
write.csv(top10_markers_identity, file = snakemake@output[["identity_markers"]], row.names = FALSE)

# Plot the top 10 markers for each identity
pdf("results/top10_markers_per_identity.pdf")
for (identity in unique(top10_markers_identity$cluster)) {
  identity_markers <- top10_markers_identity %>% filter(cluster == !!identity)
  plot <- VlnPlot(pbmc_1103, features = identity_markers$gene, pt.size = 0.1) + ggtitle(paste("Identity", identity))
  print(plot)
}
dev.off()
