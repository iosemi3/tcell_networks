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
library(SingleR)
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
subset_cells <- sample(colnames(pmbc_1103),3000)
pmbc_1103 <- subset(pmbc_1103, cells = subset_cells)

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
pmbc_1103<- RunPCA(pmbc_1103,ndims.print = c(1,3,5), verbose = T, reduction.name = "PCA_RNA")
pmbc_1103<- FindNeighbors(pmbc_1103, dims = 1:30, reduction = "PCA_RNA" )
pmbc_1103<- FindClusters(pmbc_1103, resolution = 0.2, reduction = "PCA_RNA")
# marker genes
cluster_mar_gene=FindAllMarkers(pmbc_1103, only.pos = F, test.use = 'negbinom')

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
pmbc_1103<- RunUMAP(pmbc_1103, dims = 1:30,reduction = 'PCA_RNA',
                    reduction.name ='umap' )

saveRDS(pmbc_1103, pmbc_rdata_file)
