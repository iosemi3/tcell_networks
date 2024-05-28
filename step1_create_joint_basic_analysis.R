
library(ggplot2)
library(dplyr)
library(celldex)
library(Signac)
library(Seurat)
library(chromVAR)
library(motifmatchr)
library(Matrix)
library(SummarizedExperiment)
library(BiocParallel)
library(hdf5r)
#library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(EnsDb.Hsapiens.v86)
library(JASPAR2020)
library(TFBSTools)
#library(rhdf5)
library(plotly)
set.seed(1234)
################-----SNAKEMAKE
input_feature_matrix <- snakemake@input$input_feature_matrix
input_fragments <- snakemake@input$input_fragments
#input_feature_matrix <- Read10X_h5("/Users/sebas/documents/atac_tcell/from_server/filtered_feature_bc_matrix.h5")
#input_fragments <- "/Users/sebas/documents/atac_tcell/from_server/atac_fragments.tsv.gz"

umap_plot_file <- "results/umap_plot.pdf"
piechart_plot_file<- 'results/piechart.pdf'
pmbc_rdata_file<- 'results/pmbc_1103.rds'


################
################
################-----BLOCK 1: DATA PROCESSING 
################
################

#-=-=-=-=-=-=-=-=-=
# load the RNA and ATAC data
counts <- Read10X_h5(input_feature_matrix)
fragpath <- paste(input_fragments)
# get gene annotations for hg38
annotation <-GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"

# create a Seurat object containing the RNA adata
pbmc <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)
# create ATAC assay and add it to the object
pbmc[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)
pmbc_1103 <-pbmc
pmbc_1103
remove(counts)
remove(annotation)

############---RNA-------#############

pmbc_1103[["percent.mt"]] <- PercentageFeatureSet(pmbc_1103, pattern = "^MT")
pmbc_1103[['percent.ribo']]<- PercentageFeatureSet(pmbc_1103, "^RP[SL]")

###-filtering by nr of features and prct mt-###
pmbc_1103<- subset(pmbc_1103, 
                   subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 7.5 &nCount_RNA > 500 &nCount_RNA<11000)

pmbc_1103<- NormalizeData(pmbc_1103,normalization.method = 'LogNormalize')
#cell cycle scoring
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
pmbc_1103<- CellCycleScoring(pmbc_1103, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
pmbc_1103<- ScaleData(pmbc_1103, vars.to.regress = c('S.Score','G2M.Score'))

###-HIGHLY VARIABLE FEATURE ANALYSIS-###

pmbc_1103<-FindVariableFeatures(pmbc_1103, selection.method = 'mean.var.plot',binning.method = 'equal_frequency', verbose = T)
#pmbc_1103<- ScaleData(pmbc_1103, features = rownames(pmbc_1103))

#####
pmbc_1103<- RunPCA(pmbc_1103,ndims.print = c(1,3,5), verbose = T)
DimPlot(pmbc_1103)
DimHeatmap(pmbc_1103,dims = c(1,2))
pmbc_1103<- FindNeighbors(pmbc_1103, dims = 1:30,reduction = 'PCA_RNA')
pmbc_1103<- FindClusters(pmbc_1103, resolution = 0.2)
# marker genes
cluster_mar_gene=FindAllMarkers(pmbc_1103, only.pos = F, test.use = 'negbinom')
write.table(cluster_markers,'marker_genes_cluster_tcellmultiome.tsv')

#########
ElbowPlot(pmbc_1103,reduction = 'PCA_RNA')
pmbc_1103 <- RunPCA(pmbc_1103)

###cell anotation
### automatic cell annotation by SingleR
##-PREP
#####
dice<-celldex::DatabaseImmuneCellExpressionData()
#####
##Dice
#####
pred.dice <- SingleR(test = as.SingleCellExperiment(pmbc_1103), ref = dice, assay.type.test=1,
                      labels = dice$label.fine)
##### integrating
pmbc_1103@meta.data<-cbind(pmbc_1103@meta.data,pred.dice)
colnames(pmbc_1103@meta.data)[which(names(pmbc_1103@meta.data) == "labels")] <- "labels.DICE"
Idents(pmbc_1103)=pmbc_1103@meta.data$labels.DICE
#dropping cells that belong to: 
#annotated categories with less than 150 cells
#categories that are unexpected from this experiment and might indicate contamination
lessthan150=c('Monocytes, CD16+',
              'B cells, naive',
              'T cells, CD4+, TFH',
              'T cells, CD8+, naive',
              'T cells, CD8+, naive, stimulated ')
pmbc_1103=subset(x = pmbc_1103, idents = c('Monocytes, CD16+',
                                 'B cells, naive',
                                 'T cells, CD4+, TFH',
                                 'T cells, CD8+, naive',
                                 'T cells, CD8+, naive, stimulated'),
       invert = TRUE)

pmbc_1103<- RunUMAP(pmbc_1103, dims = 1:30,reduction = 'PCA_RNA',
                    reduction.name ='umap' )
#export as pdf
umap_plot<-DimPlot(pmbc_1103, reduction = "umap",pt.size = 1.5)
pdf(umap_plot_file)
print(umap_plot)
dev.off()

#piechart
piechart <- tibble(
  cluster = pmbc_1103$seurat_clusters,
  cell_type = pmbc_1103$labels.DICE
) %>%
  group_by(cluster,cell_type) %>%
  count() %>%
  group_by(cluster) %>%
  mutate(
    percent=(100*n)/sum(n)
  ) %>%
  ungroup() %>%
  mutate(
    cluster=paste("Cluster",cluster)
  ) %>%
  ggplot(aes(x="",y=percent, fill=cell_type)) +
  geom_col(width=1) +scale_color_manual(values = c( "T cells, CD4+, naive, stimulated"="darkgreen",
                                                    "T cells, CD4+, naive"="red",
                                                    "T cells, CD4+, naive TREG"="black",
                                                    "T cells, CD4+, Th1"="yellow4",
                                                    "T cells, CD4+, Th1_17"="royalblue1",
                                                    "T cells, CD4+, Th2"="darkmagenta")) +
  coord_polar("y", start=0) +
  facet_wrap(vars(cluster)) +  
  theme(axis.text.x=element_blank()) +
  xlab(NULL) +
  ylab(NULL)+ggtitle('Celltype abundance per cluster RNA')

pdf(piechart)
print(piechart)
dev.off()

#export pmbc_1103
saveRDS(pmbc_1103, pmbc_rdata_file)
