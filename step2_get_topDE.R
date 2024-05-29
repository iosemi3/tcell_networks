library(Seurat)
library(dplyr)

# Load your Seurat object
pbmc_1103 <- readRDS(snakemake@input[["seurat_object"]])

# Perform differential expression analysis
all_markers <- FindAllMarkers(pbmc_1103, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Get the top 10 markers for each cluster
top10_markers <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

# Save to CSV
write.csv(top10_markers, file = snakemake@output[["cluster_markers"]], row.names = FALSE)

# Perform the analysis based on specific identities instead of clusters
Idents(pbmc_1103) <- pbmc_1103$labels.DICE
all_markers_identity <- FindAllMarkers(pbmc_1103, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Get the top 10 markers for each identity
top10_markers_identity <- all_markers_identity %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

# Save to CSV
write.csv(top10_markers_identity, file = snakemake@output[["identity_markers"]], row.names = FALSE)
