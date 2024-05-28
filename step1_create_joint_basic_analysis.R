# Load necessary libraries
library(Seurat)
library(SingleR)
library(celldex)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(dplyr)

# Set seed for reproducibility
set.seed(1234)

# Define input and output paths
input_feature_matrix <- snakemake@input$input_feature_matrix
input_fragments <- snakemake@input$input_fragments
umap_plot_file <- "results/umap_plot.pdf"
piechart_plot_file <- 'results/piechart.pdf'
pmbc_rdata_file <- 'results/pmbc_1103.rds'

# Load the RNA and ATAC data
counts <- Read10X_h5(input_feature_matrix)
fragpath <- input_fragments

# Get gene annotations for hg38
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"

# Create a Seurat object containing the RNA data
pbmc <- CreateSeuratObject(counts = counts$`Gene Expression`, assay = "RNA")

# Create ATAC assay and add it to the object
pbmc[["ATAC"]] <- CreateChromatinAssay(counts = counts$Peaks, sep = c(":", "-"), fragments = fragpath, annotation = annotation)

# Subset to 1000 cells for faster processing
subset_cells <- sample(colnames(pbmc), 1000)
pbmc <- subset(pbmc, cells = subset_cells)

# RNA data processing
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT")
pbmc[['percent.ribo']] <- PercentageFeatureSet(pbmc, pattern = "^RP[SL]")

# Filtering cells
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 7.5 & nCount_RNA > 500 & nCount_RNA < 11000)

# Normalize data and perform cell cycle scoring
pbmc <- NormalizeData(pbmc, normalization.method = 'LogNormalize')
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
pbmc <- CellCycleScoring(pbmc, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
pbmc <- ScaleData(pbmc, vars.to.regress = c('S.Score', 'G2M.Score'))

# Find highly variable features and perform PCA
pbmc <- FindVariableFeatures(pbmc, selection.method = 'mean.var.plot', verbose = TRUE)
pbmc <- RunPCA(pbmc, npcs = 30, verbose = TRUE)

# Find neighbors and clusters
pbmc <- FindNeighbors(pbmc, dims = 1:30)
pbmc <- FindClusters(pbmc, resolution = 0.2)

# Automatic cell annotation by SingleR
dice <- celldex::DatabaseImmuneCellExpressionData()
pred.dice <- SingleR(test = as.SingleCellExperiment(pbmc), ref = dice, assay.type.test = 1, labels = dice$label.fine)
pbmc@meta.data <- cbind(pbmc@meta.data, pred.dice)
colnames(pbmc@meta.data)[which(names(pbmc@meta.data) == "labels")] <- "labels.DICE"
Idents(pbmc) <- pbmc@meta.data$labels.DICE

# Dropping unwanted cells
unwanted_cells <- c('Monocytes, CD16+', 'B cells, naive', 'T cells, CD4+, TFH', 'T cells, CD8+, naive', 'T cells, CD8+, naive, stimulated')
pbmc <- subset(pbmc, idents = unwanted_cells, invert = TRUE)

# Run UMAP
pbmc <- RunUMAP(pbmc, dims = 1:30, reduction = 'pca', reduction.name = 'umap')

# Export UMAP plot
umap_plot <- DimPlot(pbmc, reduction = "umap", pt.size = 1.5)
pdf(umap_plot_file)
print(umap_plot)
dev.off()

# Create and export piechart
piechart <- pbmc@meta.data %>%
  group_by(seurat_clusters, labels.DICE) %>%
  summarise(count = n()) %>%
  mutate(percent = count / sum(count) * 100) %>%
  ggplot(aes(x = "", y = percent, fill = labels.DICE)) +
  geom_col(width = 1) +
  coord_polar("y", start = 0) +
  facet_wrap(~ seurat_clusters) +
  theme(axis.text.x = element_blank()) +
  labs(title = "Celltype abundance per cluster RNA")

pdf(piechart_plot_file)
print(piechart)
dev.off()

# Save the Seurat object
saveRDS(pbmc, pmbc_rdata_file)
