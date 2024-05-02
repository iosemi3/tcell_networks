if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
devtools::install_github("stuart-lab/signac", ref = "develop")



BiocManager::install("hdf5r", force = T)
BiocManager::install("biovizBase", force = T)
BiocManager::install("EnsDb.Hsapiens.v86", force = T)
BiocManager::install("JASPAR2020", force = T)
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38", force = T)
BiocManager::install("celldex", force = T)
BiocManager::install('SingleR', force = T)
BiocManager::install('SingleCellExperiment', force=T)

install.packages('grr')
install.packages('JASPAR2020')
install.packages('rhdf5')
install.packages('plotly')
install.packages('enrichR')
install.packages("scran")
install.packages('rsvd')
install.packages('SeuratWrappers')
install.packages('Seurat')
install.packages('celldex')
install.packages('SingleR')

library(ggplot2)
library(celldex)
library(scran)
library (SingleR)
library(celldex)
library(dplyr)
library(enrichR)
setRepositories(ind=1:3) # needed to automatically install Bioconductor dependencies
install.packages("Signac")
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

setwd('/corgi/sebas/tcell_multi/multiome_local_seb_20221207')


#saving session info for martin
install.packages('pander')
library(pander)
multiome_tcell_sebastian_nov2022= data.frame(installed.packages()[,1])
rownames(multiome_tcell_sebastian_nov2022)=NULL
write.csv(x = multiome_tcell_sebastian_nov2022, 'libraries_multi_seb_nov2022.csv', sep = '')


################
################
################-----BLOCK 1: DATA PROCESSING 
################
################

#-=-=-=-=-=-=-=-=-=
# load the RNA and ATAC data
counts <- Read10X_h5("/corgi/sebas/tcell_multi/aggregated_tlibs/AGG6789_good/outs/filtered_feature_bc_matrix.h5")

fragpath <- "/corgi/sebas/tcell_multi/aggregated_tlibs/AGG6789_good/outs/atac_fragments.tsv.gz"



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

DefaultAssay(pmbc_1103) <- "RNA"

mat <- GetAssayData(object = pmbc_1103, assay = "RNA", slot = "data")
mat2<-GetAssayData(object = pmbc_1103, assay = "ATAC", slot = "data")

write.csv(mat, "rna_features.csv")




############---RNA-------#############


pmbc_1103[["percent.mt"]] <- PercentageFeatureSet(pmbc_1103, pattern = "^MT")
#pmbc_1103[['percent.ribo']]<- PercentageFeatureSet(pmbc_1103, "^RP[SL]")

# Visualize QC metrics as a violin plot
VlnPlot(pmbc_1103, features = c("nFeature_RNA", "nCount_RNA", "percent.mt",'percent.ribo'), ncol = 3,pt.size = 0)
VlnPlot(pmbc_1103,features = 'nCount_RNA',pt.size = 0)
plot1 <- FeatureScatter(pmbc_1103, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pmbc_1103, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(pmbc_1103,feature1 = 'percent.ribo',feature2 = 'nCount_RNA')
CombinePlots(plots = list(plot1, plot2))
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
pmbc_1103<- ScaleData(pmbc_1103, features = rownames(pmbc_1103))


#####

pmbc_1103<- RunPCA(pmbc_1103,reduction.name = 'PCA_RNA',ndims.print = c(1,3,5), verbose = T)
DimPlot(pmbc_1103)
DimHeatmap(pmbc_1103,dims = c(1,2))


pbmc <- JackStraw(pmbc_1103, red.num.replicate = 100)
pbmc <- ScoreJackStraw(pmbc_1103, dims = 1:20)

pmbc_1103<- FindNeighbors(pmbc_1103, dims = 1:30,reduction = 'PCA_RNA')
pmbc_1103<- FindClusters(pmbc_1103, resolution = 0.2)






# marker genes
cluster_mar_gene=FindAllMarkers(pmbc_1103, only.pos = F, test.use = 'negbinom')
write.table(cluster_markers,'marker_genes_cluster_tcellmultiome.tsv')

#### top 10 genes per cluster
cluster_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

DoHeatmap(pmbc_1103, features = top10$gene) + NoLegend()
remove(cluster_markers)



#########
ElbowPlot(pmbc_1103,reduction = 'PCA_RNA')
pmbc_1103 <- RunPCA(pmbc_1103)
VizDimLoadings(pmbc_1103, dims = 1:4, reduction = "PCA_RNA")
DimHeatmap(pmbc_1103, dims = 1:10, cells = 500, balanced = TRUE, reduction ='PCA_RNA')
########Whats behind my first 5 components? CellMarker and Reactome on genes for PCs ENRICHMENT
############# tO be done



#Dimensionality of DS
#pmbc_1103<- JackStraw(pmbc_1103, num.replicate = 100)
#pmbc_1103<- ScoreJackStraw(pmbc_1103, dims = 1:5)
#JackStrawPlot(pmbc_1103, dims = 1:5)
#orElbowplot
#ElbowPlot(pmbc_1103)


###cell anotation

VlnPlot(pmbc_1103, features = c('NKG7'))
FeaturePlot(pmbc_1103,features = 'RORC',pt.size = 1.5)
FeaturePlot(alra.out,features = 'RORC',pt.size = 1.5, cols = c("lightgrey", "red"))

DoHeatmap(pmbc_1103,
          features = c('BCL6','CD3D',
                       'ICOS','CXCR5', #TFH
                       'STAT3','IRF4','IL17A','RUNX1','BATF',
                       'RORA','APOE','PIK3C2B','RORC',#TH17
                       'TBX21','STAT1',
                       'STAT4','ANXA3',#Th1
                       'GATA3','SMAD2',
                       'STAT6','RUNX2',#Th2
                       'FOXP3',
                       'IKZF2',
                       'CCL5','CTLA4',#Treg
                       'TRGV9',
                       'KLRD1',
                       'GNLY',
                       'NKG7'), #NK
          size = 5) 


#PATHWAY ENRICHKMENT ANALYSIS. Multiple dbs available
DEenrichRPlot(pmbc_1103,ident.1 = '4',max.genes = 20,enrich.database = 'CORUM')


FeatureScatter(pmbc_1103, feature1 = 'CD4', feature2 = 'CD3D') #useful for visualizing facs like plots



###enrichr
setEnrichrSite("Enrichr")
websiteLive <- TRUE
dbs <- listEnrichrDbs()
if (is.null(dbs)) websiteLive <- FALSE
if (websiteLive) head(dbs)



dbs <- c("CellMarker_Augmented_2021", "Azimuth_Cell_Types_2021", 
  'KEGG_2021_Human','Reactome_2016')
if (websiteLive) {enriched <- enrichr(c(paste(pmbc.markers$gene[pmbc.markers$cluster=='2'])), dbs)
  }

#if (websiteLive) enriched[["CellMarker_Augmented_2021"]]
if (websiteLive) plotEnrich(enriched[[1]], showTerms = 10, numChar = 40, y = "Count", orderBy = "P.value")
###enrichr


####### RUNNING TOP20 GENES ON ENRICHR AND SEEING WASSUP

enriched <- enrichr((top20), dbs)
if (websiteLive) plotEnrich(enriched[[2]], showTerms = 10, numChar = 40, y = "Count", orderBy = "P.value")



### automatic cell annotation by SingleR
##-PREP
#####

dice<-celldex::DatabaseImmuneCellExpressionData()
monaco<-celldex::MonacoImmuneData()
hpca<- celldex::HumanPrimaryCellAtlasData()
#####
##Dice
#####
pred.dice <- SingleR(test = as.SingleCellExperiment(pmbc_1103), ref = dice, assay.type.test=1,
                      labels = dice$label.fine)
table(pred.dice$labels)


#####
##HPCA
#####
pred.HPCA <- SingleR(test = as.SingleCellExperiment(pmbc_1103), ref = hpca, assay.type.test=1,
                     labels = hpca$label.fine)
table(pred.HPCA$labels)
#####

##Monaco
#####
pred.monaco <- SingleR(test = as.SingleCellExperiment(pmbc_1103), ref = monaco, assay.type.test=1,
                     labels = monaco$label.fine)
table(pred.monaco$labels)
##Plotting
#####

par(mai=c(1,2,1,1))
barplot(table(pred.monaco$labels),horiz = T, las=2, main = 'Monaco REF')
barplot(table(pred.dice$labels),horiz = T, las=2, main = 'DICE REF')
barplot(table(pred.HPCA$labels),horiz = T, las=2, main = 'HPCA')
#####

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

DimPlot(pmbc_1103, group.by = 'labels.DICE', pt.size =1.5)+ggtitle('scRNA seq')
#2D UMAP
pmbc_1103<- RunUMAP(pmbc_1103, dims = 1:30,reduction = 'PCA_RNA',
                    reduction.name ='X2D_umap_RNA' )


#3D UMAP
pmbc_1103<- RunUMAP(pmbc_1103, dims = 1:30,reduction = 'PCA_RNA',n.components = 3L,
                    reduction.name ='3D_umap' )

DimPlot(pmbc_1103, reduction = "X2D_umap_RNA",pt.size = 1.5)

# looks like our data has a 3d structure for the UMAP
#using plotly for 3d umap (thanks to: https://github.com/Dragonmasterx87/Interactive-3D-Plotting-in-Seurat-3.0.0/blob/master/3D%20UMAP%20Plotting%20v1.3.R)
# Visualize what headings are called so that you can extract them to form a dataframe
Embeddings(object = pmbc_1103, reduction = 'X3D_umap')

# Prepare a dataframe for cell plotting
plot.data <- FetchData(object = pmbc_1103, vars = c("x3d_umap_1", "x3d_umap_2", "x3d_umap_3",
                                                    'labels.DICE'))

# Make a column of row name identities (these will be your cell/barcode names)
plot.data$label <- paste(rownames(plot.data))

# Plot your data, in this example my Seurat object had 21 clusters (0-20)
fig <- plot_ly(data = plot.data, 
               x = ~x3d_umap_1, y = ~x3d_umap_2, z = ~x3d_umap_3, 
               color = ~labels.DICE, 
               colors = c(
                 
                 "darkgreen",
                 
                 "red",
                 
                 "black",
                 "yellow4",
                 "royalblue1",
                 
                 "darkmagenta"),
               type = "scatter3d", 
               mode = "markers", 
               marker = list(size = 3, width=2), # controls size of points
               text=~label, #This is that extra column we made earlier for which we will use for cell ID
               hoverinfo="text") #When you visualize your plotly object, hovering your mouse pointer over a point shows cell names


#piechart
tibble(
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


remove(dice)
remove(pred.dice)
remove(g2m.genes)
remove(s.genes)
remove(lessthan150)


########################################################################

############------ATAC----------###############

###QC 
DefaultAssay(pmbc_1103) <- "ATAC"

pmbc_1103 <- NucleosomeSignal(pmbc_1103)

pmbc_1103 <- TSSEnrichment(pmbc_1103)

VlnPlot(
  object = pmbc_1103,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 4,
  pt.size = 0
)
#qc plots
pmbc_1103$high.tss <- ifelse(pmbc_1103$TSS.enrichment > 2, 'High', 'Low')
TSSPlot(pmbc_1103, group.by = 'high.tss') + NoLegend()

pmbc_1103$nucleosome_group <- ifelse(pmbc_1103$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = pmbc_1103, group.by = 'nucleosome_group')
# filter out low quality cells
VlnPlot(
  object = pmbc_1103,
  features = c("nCount_RNA","nCount_ATAC"),
  ncol = 4,
  pt.size = 0
)
pmbc_1103 <- subset(
  x = pmbc_1103,
  subset = nCount_ATAC < 39000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)

##peak calling per RNA cluster using DICE ds
peaks <- CallPeaks(pmbc_1103, macs2.path = '/home/mahogny/miniconda3/bin/macs2',
                   group.by = 'labels.DICE')
# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(pmbc_1103),
  features = peaks,
  cells = colnames(pmbc_1103), 
  verbose = T
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
pmbc_1103[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)

#####ATAC processing
DefaultAssay(pmbc_1103) <- "peaks"
pmbc_1103 <- RunTFIDF(pmbc_1103)
pmbc_1103 <- FindTopFeatures(pmbc_1103, min.cutoff = 3, verbose = T )
pmbc_1103 <- RunSVD(pmbc_1103, verbose = T, reduction.name = 'SVD_ATAC')
#calc dimensionality reduction
DepthCor(pmbc_1103, reduction = 'SVD_ATAC',n = 20)


####--Clustering --- Removing first component due to strong correlation of depth
pmbc_1103 <- RunUMAP(object = pmbc_1103, reduction = 'SVD_ATAC', dims = 2:30, reduction.name = 'UMAP_ATAC')
pmbc_1103<- FindNeighbors(object = pmbc_1103, reduction = 'SVD_ATAC', dims = 2:30)
pmbc_1103 <- FindClusters(object = pmbc_1103, verbose = FALSE, algorithm = 3)
DimPlot(object = pmbc_1103, pt.size = 1,group.by = 'labels.DICE')+ ggtitle('scATAC seq')
####3D UMAP more informative
pmbc_1103 <- RunUMAP(object = pmbc_1103, reduction = 'SVD_ATAC', dims = 2:30,n.components = 3L, reduction.name = 'X3D_UMAP_ATAC')
# Visualize what headings are called so that you can extract them to form a dataframe
Embeddings(object = pmbc_1103, reduction = 'X3D_UMAP_ATAC')

# Prepare a dataframe for cell plotting
plot.data <- FetchData(object = pmbc_1103, vars = c("x3d_umap_atac_1", "x3d_umap_atac_2", "x3d_umap_atac_3",
                                                    'labels.DICE'))

# Make a column of row name identities (these will be your cell/barcode names)
plot.data$label <- paste(rownames(plot.data))

# Plot 
fig_atac <- plot_ly(data = plot.data, 
               x = ~x3d_umap_atac_1, y = ~x3d_umap_atac_2, z = ~x3d_umap_atac_3, 
               color = ~labels.DICE, 
               colors = c(
                 
                 "darkgreen",
                 
                 "red",
                 
                 "black",
                 "yellow4",
                 "royalblue1",
                 
                 "darkmagenta"),
               type = "scatter3d", 
               mode = "markers", 
               marker = list(size = 3, width=2), # controls size of points
               text=~label, #This is that extra column we made earlier for which we will use for cell ID
               hoverinfo="text") #When you visualize your plotly object, hovering your mouse pointer over a point shows cell names
fig_atac
#piechart---not very informative
tibble(
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
  ylab(NULL)+ggtitle('Celltype abundance per cluster ATAC')





#####---- Joint RNA+ATAC umap
# build a joint neighbor graph using both assays
pmbc_1103 <- FindMultiModalNeighbors(
  object = pmbc_1103,
  reduction.list = list('SVD_ATAC', 'PCA_RNA'), 
  dims.list = list(2:30, 1:30),
  verbose = TRUE
)

# build a joint UMAP visualization
pmbc_1103<- RunUMAP(
  object = pmbc_1103,
  nn.name = "weighted.nn",
  assay = 'ATAC',
  verbose = TRUE
)

DimPlot(pmbc_1103, reduction = "umap", pt.size = 1.5, group.by = 'labels.DICE')+ ggtitle('Joint sc UMAP')
#3D dimension UMAP required.


pmbc_1103<- RunUMAP(
  object = pmbc_1103,
  nn.name = "weighted.nn",
  n.components =  3L,
  assay = 'ATAC',
  verbose = TRUE,
  reduction.name = 'X3D_UMAP_JOINT'
)

# Visualize what headings are called so that you can extract them to form a dataframe
Embeddings(object = pmbc_1103, reduction = 'X3D_UMAP_JOINT')

# Prepare a dataframe for cell plotting
plot.data <- FetchData(object = pmbc_1103, vars = c("x3d_umap_joint_1", "x3d_umap_joint_2", "x3d_umap_joint_3",
                                                    'labels.DICE'))

# Make a column of row name identities (these will be your cell/barcode names)
plot.data$label <- paste(rownames(plot.data))

# Plot 
fig_joint <- plot_ly(data = plot.data, 
                    x = ~x3d_umap_joint_1, y = ~x3d_umap_joint_2, z = ~x3d_umap_joint_3, 
                    color = ~labels.DICE, 
                    colors = c(
                      
                      "darkgreen",
                      
                      "red",
                      
                      "black",
                      "yellow4",
                      "royalblue1",
                      
                      "darkmagenta"),
                    type = "scatter3d", 
                    mode = "markers", 
                    marker = list(size = 3, width=2), # controls size of points
                    text=~label, #This is that extra column we made earlier for which we will use for cell ID
                    hoverinfo="text") #When you visualize your plotly object, hovering your mouse pointer over a point shows cell names

fig_joint


#######----
#######----
#######----
#######----
#######----END BLOCK 1
#######----
#######----
#######----

#Activation BLOCK
####How many genes overlap between naive, naive stim, naive Treg
naive_stim_markers=FindMarkers(pmbc_1103, 
            ident.1 = 'T cells, CD4+, naive, stimulated',
            assay = 'RNA',
            test.use = 'negbinom',
            logfc.threshold = 0.3,
            min.cells.feature = 3,
            verbose = T)

naive_markers=FindMarkers(pmbc_1103, 
                          ident.1 = 'T cells, CD4+, naive',
                          assay = 'RNA',
                          test.use = 'negbinom',
                          logfc.threshold = 0.3,
                          min.cells.feature = 3,
                          verbose = T)
naive_treg_markers=FindMarkers(pmbc_1103, 
                               ident.1 = 'T cells, CD4+, naive TREG',
                               assay = 'RNA',
                               test.use = 'negbinom',
                               logfc.threshold = 0.3,
                               min.cells.feature = 3,
                               verbose = T)

venn.diagram(
       x = list(rownames(naive_markers), row.names(naive_stim_markers), rownames(naive_treg_markers)),
       category.names = c("Naive" , "Naive stimulated" , "TREG naive"),
       filename = 'NaiveVennGenemarkers.jpeg')



naive_markers$genes=rownames(naive_markers)
naive_stim_markers$genes=rownames(naive_stim_markers)
shared1=merge.data.frame(naive_markers,naive_stim_markers,by = 'genes')
#enrichment of cluster genes

###enrichr
setEnrichrSite("Enrichr")
websiteLive <- TRUE
dbs <- listEnrichrDbs()
if (is.null(dbs)) websiteLive <- FALSE
if (websiteLive) head(dbs)



dbs <- c("KEGG_2019_Human", "GO_Cellular_Component_2021", 
         'KEGG_2021_Human','Reactome_2016','GO_Molecular_Function_2021',
         'GO_Biological_Process_2021')
if (websiteLive) {enriched <- enrichr(intersect(rownames(naive_markers),rownames(naive_treg_markers)), dbs)
}

#taking the go biological processes 2021 ds
if (websiteLive) plotEnrich(enriched[[6]], 
                            showTerms = 5,numChar = 60, 
                            y = "Count", orderBy = "P.value",
                            title = 'GO enrichment Naive & Naive Treg intersect')
###enrichr





####---ChromVAR
DefaultAssay(pmbc_1103) <- "peaks"
#Get motifs from JASPAR
# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

pfm_test=getMatrixSet(
  x= JASPAR2020,
  opts= list(collect = 'FAM', tax_group= 'vertebrates', all_versions = F)
)


#add motif info
pmbc_1103<- AddMotifs(pmbc_1103, 
                 genome = BSgenome.Hsapiens.UCSC.hg38,pfm = pfm,verbose = T)

#COMPUTE MOTIF ACTIVITY
pbmc_chromvared<- RunChromVAR(
  object = pmbc_1103,
  genome = BSgenome.Hsapiens.UCSC.hg38, assay = 'peaks'
)

DefaultAssay(pbmc_chromvared) <- 'chromvar'

#------------------correlation BW chromvar output TFs and expression
#gene jaspar and uniprot dictionary
jaspmotifdict= read.csv('/Users/sebas/Documents/multimodal Tcell/jaspar_uniprot.csv', header = T, sep = ',')
jaspmotifdict=jaspmotifdict[,c(2,3)]
colnames(jaspmotifdict)= c('uni_id','gene')
# gene ensembl and uniprot dictionary
unidict=read.csv('/Users/sebas/Documents/multimodal Tcell/genes_uniprot_mart_export.txt', header = T, sep = '\t')
colnames(unidict)= c('uni_id','gene')
unidict=unidict[!unidict$uni_id=='',]
#motif jaspar and gene jaspar dictionary
motifjaspar=read.csv('/Users/sebas/Documents/multimodal Tcell/JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.txt', header = F, sep= '\t')
motifjaspar=motifjaspar[grep('MA',motifjaspar$V1),]
colnames(motifjaspar)=c('motif_jaspar','gene_jaspar')
#final gene motif dictionary TO USE
jaspdict= merge.data.frame(jaspmotifdict,unidict,by = 'uni_id')
colnames(jaspdict)= c('uni_id','gene_jaspar','gene_ensembl')
jaspdict= merge.data.frame(jaspdict,motifjaspar, by= 'gene_jaspar')
remove(motifjaspar)
remove(unidict)
remove(jaspmotifdict)

motif.counts=as.data.frame(GetAssayData(pbmc_chromvared, slot = 'data'))
motif.counts$motif_jaspar=rownames(motif.counts)
motif.counts=merge.data.frame(motif.counts,jaspdict,by = 'motif_jaspar')

#second try with discrete gene expression and discrete zscore (per cell)
test.melt=melt(motif.counts)
test.melt=test.melt[,c('gene_ensembl','variable','value')]
DefaultAssay(pmbc_1103) <- 'RNA'
gene.counts=as.data.frame(GetAssayData(pmbc_1103,slot = 'data'))
gene.counts$gene_ensembl=rownames(gene.counts)
gene.counts_jaspar.motif= merge.data.frame(gene.counts, jaspdict, by = 'gene_ensembl')

test.melt2=melt(gene.counts_jaspar.motif)
test.melt2$imp.gene.exp=test.melt2$value
test.melt2= within.data.frame(test.melt2,rm(value))
test.melt2=test.melt2[,c('gene_ensembl','variable','imp.gene.exp')]



common.barcode=intersect(test.melt$variable,test.melt2$variable)
remove(test.melt)
remove(test.melt2)
remove(gene.counts)
remove(gene.counts_jaspar.motif)




gene.rna=as.data.frame(GetAssayData(pmbc_1103,slot = 'data'))
DefaultAssay(pbmc_chromvared)='chromvar'
motif.counts=as.data.frame(GetAssayData(pbmc_chromvared, slot = 'data'))

gene.rna = gene.rna[,common.barcode]
motif.counts = motif.counts[,common.barcode]

jaspdict = merge(
  data.frame(genei = 1:nrow(gene.rna), gene_ensembl=rownames(gene.rna)),
  jaspdict)
jaspdict = merge(
  data.frame(motifi = 1:nrow(motif.counts), motif_jaspar=rownames(motif.counts)),
  unique(jaspdict))

jaspdict$corr = NA
for(i in 1:nrow(jaspdict)){
  jaspdict$corr[i] = cor(
    as.double(gene.rna[jaspdict$genei[i],]),
    as.double(motif.counts[jaspdict$motifi[i],]))
}
head(jaspdict)

jaspdict = jaspdict[!is.na(jaspdict$corr),]


jaspdict = jaspdict[order(jaspdict$corr),]
plot(sort(jaspdict$corr))

jaspdict
tail(jaspdict, n=30)

jaspdict[jaspdict$gene_ensembl=="GATA3",]
jaspdict[jaspdict$gene_ensembl=="FOXP3",]
jaspdict[jaspdict$gene_ensembl=="TBX21",]
jaspdict[jaspdict$gene_ensembl=="BCL6",]

jaspdict = jaspdict[,c("gene_ensembl","gene_jaspar","corr")]


#optional BELOW

jaspdict1 = data.frame(motif=jaspdict$gene_jaspar, symbol1=jaspdict$gene_ensembl, corr1=jaspdict$corr)
jaspdict2 = data.frame(motif=jaspdict$gene_jaspar, symbol2=jaspdict$gene_ensembl, corr2=jaspdict$corr)

jasphetro = merge(jaspdict1, jaspdict2)
jasphetro[jasphetro$symbol1!=jasphetro$symbol2,]

##############------- CHECKING DIFFERENTIALLY OPEN PEAKS
#Get motifs from JASPAR
# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

pfm_test=getMatrixSet(
  x= JASPAR2020,
  opts= list(collect = 'FAM', tax_group= 'vertebrates', all_versions = F)
)


#add motif info
pbmc<- AddMotifs(pbmc, 
                 genome = BSgenome.Hsapiens.UCSC.hg38,pfm = pfm)



#Naive vs naivestimulated
da_peaks_naiveall <- FindMarkers(
  object = pmbc_1103,
  ident.1 = 'T cells, CD4+, naive, stimulated',
  ident.2 = 'T cells, CD4+, naive',
  logfc.threshold = 0.3,
  test.use = 'negbinom',
  min.cells.feature = 3,
  latent.vars = 'nCount_peaks'
)
# get top differentially accessible peaks
top.da.peak <- rownames(da_peaks_naivestim[da_peaks_naivestim$p_val_adj< 0.005,])


# find peaks open in  cells.
#Chosing set of background peaks.
open.peaks <- AccessiblePeaks(pmbc_1103,idents = c('T cells, CD4+, naive, stimulated','T cells, CD4+, naive'))
meta.feature <- GetAssayData(pmbc_1103, assay = "peaks", slot = "meta.features")
peaks.matched <- MatchRegionStats(
  meta.feature = meta.feature[open.peaks, ],
  query.feature = meta.feature[top.da.peak, ],
  verbose = T
)
# test enriched motifs
enriched.motifs.naiveall <- FindMotifs(
  object = pmbc_1103,
  background = peaks.matched,
  features = top.da.peak
)
MotifPlot(
  object = pmbc_1103,
  motifs = head(rownames(enriched.motifs.naiveall))
)


#Naive vs naive T reg
da_peaks_naivetregcomp <- FindMarkers(
  object = pmbc_1103,
  ident.1 = 'T cells, CD4+, naive TREG',
  logfc.threshold = 0.2,
  test.use = 'negbinom',
  min.cells.feature = 3,
  latent.vars = 'nCount_peaks'
)
# get top differentially accessible peaks
top.da.peak <- rownames(da_peaks_naivetregcomp[da_peaks_naivetregcomp$p_val_adj< 0.005,])
# find peaks open in  cells.
#Chosing set of background peaks.
open.peaks <- AccessiblePeaks(pmbc_1103,idents = c('T cells, CD4+, naive TREG'))
meta.feature <- GetAssayData(pmbc_1103, assay = "peaks", slot = "meta.features")
peaks.matched <- MatchRegionStats(
  meta.feature = meta.feature[open.peaks, ],
  query.feature = meta.feature[top.da.peak, ],
  verbose = T
)
# test enriched motifs
enriched.motifs.naivetregcomp <- FindMotifs(
  object = pmbc_1103,
  background = peaks.matched,
  features = top.da.peak
)
MotifPlot(
  object = pmbc_1103,
  motifs = head(rownames(enriched.motifs.naivetregcomp))
)

####creating heatmap with correlation values between genexp and motif abundance in open region
##for the enriched motifs for each subclass
####### naive all comparison
commons.filter= intersect(jaspdict$gene_jaspar,enriched.motifs.naiveall$motif.name[1:10])
naive_jaspdict= jaspdict[jaspdict$gene_jaspar %in% commons.filter,] 

ggheatmap <- ggplot(naive_jaspdict, aes(gene_ensembl, gene_jaspar, fill = signif(corr,digits = 2)))+
  geom_tile(color = "white",)+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()+
  geom_text(aes(gene_ensembl, gene_jaspar, label = signif(corr,digits = 2)), color = "black", size = 3) +
  theme(
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.8),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5)) +
  labs(y= "Motif", x = "Gene") 
# Print the heatmap
print(ggheatmap)

###### naive treg
commons.filter= intersect(jaspdict$gene_jaspar,enriched.motifs.naivetregcomp$motif.name[enriched.motifs.naivetregcomp$pvalue<0.05])
naivetreg_jaspdict= jaspdict[jaspdict$gene_jaspar %in% commons.filter,] 

ggheatmap <- ggplot(naivetreg_jaspdict, aes(gene_ensembl, gene_jaspar, fill = signif(corr,digits = 2)))+
  geom_tile(color = "white",)+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()+
  geom_text(aes(gene_ensembl, gene_jaspar, label = signif(corr,digits = 2)), color = "black", size = 3) +
  theme(
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.8),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5)) +
  labs(y= "Motif", x = "Gene") 
# Print the heatmap
print(ggheatmap)

#######----
#######----
#######----
#######----
#######----END BLOCK 2: activation
#######----
#######----
#######----


######### BlOCK 3: CYTOKINE STIMULATION
#finding markers via intersecting with naive
#finding overrepresented/ enriched motifs for each category
#checking correlation between motif and gene expression
######------TH1
th1_markers=FindMarkers(pmbc_1103, 
                               ident.1 = 'T cells, CD4+, Th1',
                              #ident.2 = 'T cells, CD4+, naive',
                               assay = 'RNA',
                               test.use = 'negbinom',
                               logfc.threshold = 0.3,
                               min.cells.feature = 3,
                               verbose = T)
a=c(intersect(rownames((naive_markers)),rownames((th1_markers)))[1:10])

#motifs
da_peaks_th1 <- FindMarkers(
  object = pmbc_1103,
  ident.1 = 'T cells, CD4+, naive',
  ident.2 = 'T cells, CD4+, Th1',
  logfc.threshold = 0.3,
  test.use = 'negbinom',
  min.cells.feature = 3,
  latent.vars = 'nCount_peaks'
)
# get top differentially accessible peaks
top.da.peak <- rownames(da_peaks_th1[da_peaks_th1$p_val_adj< 0.005,])
# find peaks open in  cells.
#Chosing set of background peaks.
open.peaks <- AccessiblePeaks(pmbc_1103,idents = c('T cells, CD4+, naive','T cells, CD4+, Th1'))
meta.feature <- GetAssayData(pmbc_1103, assay = "peaks", slot = "meta.features")
peaks.matched <- MatchRegionStats(
  meta.feature = meta.feature[open.peaks, ],
  query.feature = meta.feature[top.da.peak, ],
  verbose = T
)
# test enriched motifs
enriched.motifs.th1 <- FindMotifs(
  object = pmbc_1103,
  background = peaks.matched,
  features = top.da.peak
)
enriched.motifs.th1$motif.name[enriched.motifs.th1$pvalue<0.05]
MotifPlot(
  object = pmbc_1103,
  motifs = head(rownames(enriched.motifs.th1))
)
####creating heatmap with correlation values between genexp and motif abundance in open region
##for the enriched motifs for each subclass
corrheatmap6=function(x){
commons.filter= intersect(jaspdict$gene_jaspar,x$motif.name[1:6])
x_jaspdict= jaspdict[jaspdict$gene_jaspar %in% commons.filter,] 

ggheatmap <- ggplot(x_jaspdict, aes(gene_ensembl, gene_jaspar, fill = signif(corr,digits = 2)))+
  geom_tile(color = "white",)+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()+
  geom_text(aes(gene_ensembl, gene_jaspar, label = signif(corr,digits = 2)), color = "black", size = 3) +
  theme(
    legend.justification = c(1, 0),
    legend.position = c(1.35, 0.8),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5)) +
  labs(y= "Motif", x = "Gene") 
# Print the heatmap
return (print(ggheatmap))
}
corrheatmap6(enriched.motifs.th1)

######------TH1_17
th1_17_markers=FindMarkers(pmbc_1103, 
                        ident.1 = 'T cells, CD4+, Th1_17',
                        assay = 'RNA',
                        test.use = 'negbinom',
                        logfc.threshold = 0.3,
                        min.cells.feature = 3,
                        verbose = T)
b=intersect(rownames((naive_markers)),rownames((th1_17_markers)))[1:10]
#motifs
da_peaks_th1_17 <- FindMarkers(
  object = pmbc_1103,
  ident.1 = 'T cells, CD4+, naive',
  ident.2 = 'T cells, CD4+, Th1_17',
  logfc.threshold = 0.3,
  test.use = 'negbinom',
  min.cells.feature = 3,
  latent.vars = 'nCount_peaks'
)
# get top differentially accessible peaks
top.da.peak <- rownames(da_peaks_th1_17[da_peaks_th1_17$p_val_adj< 0.005,])
# find peaks open in  cells.
#Chosing set of background peaks.
open.peaks <- AccessiblePeaks(pmbc_1103,idents = c('T cells, CD4+, naive','T cells, CD4+, Th1_17'))
meta.feature <- GetAssayData(pmbc_1103, assay = "peaks", slot = "meta.features")
peaks.matched <- MatchRegionStats(
  meta.feature = meta.feature[open.peaks, ],
  query.feature = meta.feature[top.da.peak, ],
  verbose = T
)
# test enriched motifs
enriched.motifs.th1_17 <- FindMotifs(
  object = pmbc_1103,
  background = peaks.matched,
  features = top.da.peak
)
enriched.motifs.th1_17$motif.name[enriched.motifs.th1_17$pvalue<0.05]
MotifPlot(
  object = pmbc_1103,
  motifs = head(rownames(enriched.motifs.th1_17))
)
#correlation heatmap
corrheatmap6(enriched.motifs.th1_17)

######------TH2
th2_markers=FindMarkers(pmbc_1103, 
                           ident.1 = 'T cells, CD4+, Th2',
                           assay = 'RNA',
                           test.use = 'negbinom',
                           logfc.threshold = 0.3,
                           min.cells.feature = 3,
                           verbose = T)
c=intersect(rownames((naive_markers)),rownames((th2_markers)))[1:10]
da_peaks_th2 <- FindMarkers(
  object = pmbc_1103,
  ident.1 = 'T cells, CD4+, naive',
  ident.2 = 'T cells, CD4+, Th2',
  logfc.threshold = 0.3,
  test.use = 'negbinom',
  min.cells.feature = 3,
  latent.vars = 'nCount_peaks'
)
# get top differentially accessible peaks
top.da.peak <- rownames(da_peaks_th2[da_peaks_th2$p_val_adj< 0.005,])
# find peaks open in  cells.
#Chosing set of background peaks.
open.peaks <- AccessiblePeaks(pmbc_1103,idents = c('T cells, CD4+, naive','T cells, CD4+, Th2'))
meta.feature <- GetAssayData(pmbc_1103, assay = "peaks", slot = "meta.features")
peaks.matched <- MatchRegionStats(
  meta.feature = meta.feature[open.peaks, ],
  query.feature = meta.feature[top.da.peak, ],
  verbose = T
)
# test enriched motifs
enriched.motifs.th2 <- FindMotifs(
  object = pmbc_1103,
  background = peaks.matched,
  features = top.da.peak
)
enriched.motifs.th2$motif.name[enriched.motifs.th2$pvalue<0.05]
MotifPlot(
  object = pmbc_1103,
  motifs = head(rownames(enriched.motifs.th2))
)
#correlation heatmap
corrheatmap6(enriched.motifs.th2)


######------TH17
th17_markers=FindMarkers(pmbc_1103, 
                        ident.1 = 'T cells, CD4+, Th17',
                        assay = 'RNA',
                        test.use = 'negbinom',
                        logfc.threshold = 0.3,
                        min.cells.feature = 3,
                        verbose = T)
d=intersect(rownames((naive_markers)),rownames((th17_markers)))[1:10]

#motif
da_peaks_th17 <- FindMarkers(
  object = pmbc_1103,
  ident.1 = 'T cells, CD4+, naive',
  ident.2 = 'T cells, CD4+, Th17',
  logfc.threshold = 0.3,
  test.use = 'negbinom',
  min.cells.feature = 3,
  latent.vars = 'nCount_peaks'
)
# get top differentially accessible peaks
top.da.peak <- rownames(da_peaks_th17[da_peaks_th17$p_val_adj< 0.005,])
# find peaks open in  cells.
#Chosing set of background peaks.
open.peaks <- AccessiblePeaks(pmbc_1103,idents = c('T cells, CD4+, naive','T cells, CD4+, Th17'))
meta.feature <- GetAssayData(pmbc_1103, assay = "peaks", slot = "meta.features")
peaks.matched <- MatchRegionStats(
  meta.feature = meta.feature[open.peaks, ],
  query.feature = meta.feature[top.da.peak, ],
  verbose = T
)
# test enriched motifs
enriched.motifs.th17 <- FindMotifs(
  object = pmbc_1103,
  background = peaks.matched,
  features = top.da.peak
)
enriched.motifs.th17$motif.name[enriched.motifs.th17$pvalue<0.05]
MotifPlot(
  object = pmbc_1103,
  motifs = head(rownames(enriched.motifs.th17))
)
#corrheatmap
corrheatmap6(enriched.motifs.th17)
######------memory T reg
memory_treg_markers=FindMarkers(pmbc_1103, 
                         ident.1 = 'T cells, CD4+, memory TREG',
                         assay = 'RNA',
                         test.use = 'negbinom',
                         logfc.threshold = 0.3,
                         min.cells.feature = 3,
                         verbose = T)
e=intersect(rownames((naive_markers)),rownames((memory_treg_markers)))[1:10]

#motif
da_peaks_mtreg <- FindMarkers(
  object = pmbc_1103,
  ident.1 = 'T cells, CD4+, naive',
  ident.2 = 'T cells, CD4+, memory TREG',
  logfc.threshold = 0.3,
  test.use = 'negbinom',
  min.cells.feature = 3,
  latent.vars = 'nCount_peaks'
)
# get top differentially accessible peaks
top.da.peak <- rownames(da_peaks_mtreg[da_peaks_mtreg$p_val_adj< 0.005,])
# find peaks open in  cells.
#Chosing set of background peaks.
open.peaks <- AccessiblePeaks(pmbc_1103,idents = c('T cells, CD4+, naive','T cells, CD4+, memory TREG'))
meta.feature <- GetAssayData(pmbc_1103, assay = "peaks", slot = "meta.features")
peaks.matched <- MatchRegionStats(
  meta.feature = meta.feature[open.peaks, ],
  query.feature = meta.feature[top.da.peak, ],
  verbose = T
)
# test enriched motifs
enriched.motifs.mtreg <- FindMotifs(
  object = pmbc_1103,
  background = peaks.matched,
  features = top.da.peak
)
enriched.motifs.mtreg$motif.name[enriched.motifs.mtreg$pvalue<0.05]
MotifPlot(
  object = pmbc_1103,
  motifs = head(rownames(enriched.motifs.mtreg))
)
##corr heatmap
corrheatmap6(enriched.motifs.mtreg)
######------NK
NK_markers=FindMarkers(pmbc_1103, 
                                ident.1 = 'NK cells',
                                assay = 'RNA',
                                test.use = 'negbinom',
                                logfc.threshold = 0.3,
                                min.cells.feature = 3,
                                verbose = T)
f=intersect(rownames((naive_markers)),rownames((NK_markers)))[1:10]

#motif
da_peaks_nk <- FindMarkers(
  object = pmbc_1103,
  ident.1 = 'T cells, CD4+, naive',
  ident.2 = 'NK cells',
  logfc.threshold = 0.3,
  test.use = 'negbinom',
  min.cells.feature = 3,
  latent.vars = 'nCount_peaks'
)
# get top differentially accessible peaks
top.da.peak <- rownames(da_peaks_nk[da_peaks_nk$p_val_adj< 0.005,])
# find peaks open in  cells.
#Chosing set of background peaks.
open.peaks <- AccessiblePeaks(pmbc_1103,idents = c('T cells, CD4+, naive','NK cells'))
meta.feature <- GetAssayData(pmbc_1103, assay = "peaks", slot = "meta.features")
peaks.matched <- MatchRegionStats(
  meta.feature = meta.feature[open.peaks, ],
  query.feature = meta.feature[top.da.peak, ],
  verbose = T
)
# test enriched motifs
enriched.motifs.nk <- FindMotifs(
  object = pmbc_1103,
  background = peaks.matched,
  features = top.da.peak
)
enriched.motifs.nk$motif.name[enriched.motifs.nk$pvalue<0.05]
MotifPlot(
  object = pmbc_1103,
  motifs = head(rownames(enriched.motifs.nk))
)
##corr heatmap
corrheatmap6(enriched.motifs.nk)



#common genes expressed bw comparisons
g=list(a,b,c,d,e,f)

intersect(unique(c(a,b,c,d,e)),f)

#### extracting top enriched motifs and making a heatmap with the p values as fill, per cell type

loop.elements= noquote(ls(pattern = 'enriched.motifs.'))
lapply(loop.elements,head(loop.elements))

heatmap.data=data.frame(stringsAsFactors = F)
for(i in loop.elements){
  celltype=gsub("^.*\\.","",i)
  motif.name=head(get(i)$motif.name)
  pval=head(get(i)$pvalue)
  pval=formatC(pval, format = "e", digits = 2)
  onelevel=cbind(celltype,motif.name,pval)
  heatmap.data=rbind(heatmap.data,onelevel)
}
heatmap.data=as.matrix(acast(heatmap.data, heatmap.data$motif.name~heatmap.data$celltype,fill = 0))
class(heatmap.data)='numeric'

heatmap.data=-log10(heatmap.data)
heatmap.data[sapply(heatmap.data, is.infinite)] <- NA
heatmap.data[sapply(heatmap.data, is.na)] <- 0
Heatmap(heatmap.data, rect_gp = gpar(col='white', lwd=2))


ggplot(as.data.frame(heatmap.data),aes(rownames(heatmap.data),colnames(heatmap.data)))+geom_tile()+scale_fill_gradient()

g

#######----
#######----
#######----
#######----
#######----END BLOCK 3
#######----
#######----
#######----

######### BlOCK 4:Regulon status (Network analysis )
#Gene regulatory network analysis using SCENIC and RCIS target db.
#checking correlation between motif and gene expression


###--------Gene regulatory network analysis using SCENIC

BiocManager::install(c("AUCell", "RcisTarget",'Seurat','Signac'))

if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("aertslab/SCENIC") 
devtools::install_github("aertslab/SCopeLoomR")
library(AUCell)
library(RcisTarget)
library(GENIE3)
library(SCENIC)
library(Seurat)
library(Signac)
#files for scenic motif ranking #download db files manually, goes faster
setwd('cisTarget_databases/')
dbFiles <- c("https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-7species.mc9nr.feather",
             "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-7species.mc9nr.feather")
dir.create("cisTarget_databases");setwd("cisTarget_databases")

## Get data from sce object:
#----------------RNA mine
#----------------
cellInfo <- data.frame(seuratCluster=levels(pmbc_1103))
scenicOptions <- initializeScenic(org="hgnc", dbDir="/Users/sebas/Documents/multimodal Tcell/cisTarget_databases/", nCores=10)

### Co-expression network

gene.dgc= as.matrix(GetAssayData(pmbc_1103,slot = 'counts',assay = 'RNA'))
genesKept <- geneFiltering(gene.dgc, scenicOptions, minCountsPerGene = 500, minSamples = 3) #same filtering criteria as for seurat
#keep only genes in db
exprMat_filtered <- gene.dgc[genesKept, ] #### used also for csdr
dim(exprMat_filtered)
rm(gene.dgc)
####exporting for running on server
write.table(exprMat_filtered,'expression_matrix_filtered.txt')

#computing correlations (Spearman) (ON SERVER from here on)
runCorrelation(exprMat_filtered, scenicOptions)
#run genie3(ON SERVER)
runGenie3(exprMat_filtered, scenicOptions)


#import

### Build and score the GRN

scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] 
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) 
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions, exprMat_filtered)





#### exploring output of regulon targets
regulontargets= read.table('output/Step2_regulonTargetsInfo.tsv', header = T)
regulontargets=regulontargets[regulontargets$highConfAnnot=='TRUE',]

master_regulons=unique(regulontargets$TF)

#checking regulon activity per ct

cellInfo <- data.frame(seuratCluster=Idents(pmbc_1103))

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$seuratCluster),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- data.frame(t(scale(t(regulonActivity_byCellType), center = T, scale=T)))
regulonActivity_byCellType_Scaled$regulon=rownames(regulonActivity_byCellType_Scaled)
#selecting only non extended
regulonActivity_byCellType_Scaled=regulonActivity_byCellType_Scaled[!grepl('_extended',regulonActivity_byCellType_Scaled$regulon),]
regulonActivity_byCellType_Scaled=subset(regulonActivity_byCellType_Scaled, select= -regulon)


ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity")
#updating list of master regulons for only regulons that passed AUC 
master_regulons=master_regulons[-c(7,9)]
#checking correlation with motifs
corrheatregulon=function(x){
  commons.filter= intersect(jaspdict$gene_jaspar,x)
  x_jaspdict= jaspdict[jaspdict$gene_jaspar %in% commons.filter,] 
  
  ggheatmap <- ggplot(x_jaspdict, aes(gene_ensembl, gene_jaspar, fill = signif(corr,digits = 2)))+
    geom_tile(color = "white",)+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Pearson\nCorrelation") +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 12, hjust = 1))+
    coord_fixed()+
    geom_text(aes(gene_ensembl, gene_jaspar, label = signif(corr,digits = 2)), color = "black", size = 3) +
    theme(
      legend.justification = c(1, 0),
      legend.position = c(0.4, 0.8),
      legend.direction = "horizontal")+
    guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                                 title.position = "top", title.hjust = 0.5)) +
    labs(y= "Motif", x = "Gene") 
  # Print the heatmap
  return (print(ggheatmap))
}
corrheatregulon(master_regulons)


#exporting to loom for visualizing in scope and cytoscape
export2loom(scenicOptions, exprMat_filtered)
colnames(regulontargets)
cytoscape_export= regulontargets[,c(1:2)]
write.csv(cytoscape_export,'output/cytoscape_export.csv')


###--------Gene regulatory network analysis using PANDO
#note: the implementation is on python, but the first steps are made in R
#creating cistopic object
devtools::install_github('quadbiolab/Pando')
library(Pando)
library(grr)
library(JASPAR2020)
pmbc_1103 <- initiate_grn(pmbc_1103)
# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

pmbc_1103 <- find_motifs(
  pmbc_1103,
  pfm = pfm,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

pmbc_1103 <- infer_grn(
  pmbc_1103,
  peak_to_gene_method = 'Signac',
  method = 'glm'
)

pmbc_1103 <- find_modules(pmbc_1103)
modules <- NetworkModules(pmbc_1103)
modules@meta
table(modules@meta$tf)

#1.......generalized linear model
glm_modules=names(table(modules@meta$tf))
plot_gof(pmbc_1103, point_size=3)

plot_module_metrics(pmbc_1103) #size of modules

pmbc_1103 <- get_network_graph(pmbc_1103)
plot_network_graph(pmbc_1103) 
write.csv(glm_modules,"non_regressedcc_glm_pando_sebas.csv")
#checking canonical TFS/GENES present in the module
canon_bois=c('BCL6','CD3D',
             'ICOS','CXCR5', #TFH
             'STAT3','IRF4','IL17A','RUNX1','BATF',
             'RORA','APOE','PIK3C2B','RORC',#TH17
             'TBX21','STAT1',
             'STAT4','ANXA3',#Th1
             'GATA3','SMAD2',
             'STAT6','RUNX2',#Th2
             'FOXP3',
             'IKZF2',
             'CCL5','CTLA4',#Treg
             'TRGV9',
             'KLRD1',
             'GNLY',
             'NKG7') #NK

master_regs= c('ATF4',  
'CEBPG',
'DDIT3',
'ELF1',
'GATA3',
'JUNB',
'NFKB1',
'REL',
'RUNX3',
'SREBF2',
'YY1')
canon_bois[which(canon_bois %in% glm_modules)]
master_regs[which(master_regs %in% glm_modules)]

#
#repeating the same with a bayesian regression model
pmbc_1103 <- infer_grn(
  pmbc_1103,
  peak_to_gene_method = 'Signac',
  method = 'brms'
)

pmbc_1103 <- find_modules(pmbc_1103)
modules_brm <- NetworkModules(pmbc_1103)
modules_brm@meta
#2.......bayesian regression model
brm_modules=names(table(modules_brm@meta$tf))

brm_modules=names(table(modules@meta$tf))
plot_gof(pmbc_1103, point_size=3)

plot_module_metrics(pmbc_1103) #size of modules

pmbc_1103 <- get_network_graph(pmbc_1103)
plot_network_graph(pmbc_1103) 
write.csv(brm_modules,"non_regressedcc_brm_pando_sebas.csv")

canon_bois[which(canon_bois %in% brm_modules)]
master_regs[which(master_regs %in% brm_modules)]


#repeating the same with a gradient boost regression model (SAME AS IN SCENIC)

pmbc_1103 <- infer_grn(
  pmbc_1103,
  peak_to_gene_method = 'Signac',
  method = 'xgb'
)

pmbc_1103 <- find_modules(pmbc_1103)
modules_xgb <- NetworkModules(pmbc_1103)
#3.......gradient boost regression model
modules_xgb@meta$tf

grm_modules=names(table(modules_xgb@meta$tf))

table(modules_xgb@meta$tf)

write.csv(grm_modules,"non_regressedcc_xgboost_pando_sebas.csv")
#which canonicals are within the modules
canon_bois[which(canon_bois %in% grm_modules)]
master_regs[which(master_regs %in% grm_modules)]
###plotting
plot_gof(pmbc_1103, point_size=3)

plot_module_metrics(pmbc_1103) #size of modules

pmbc_1103 <- get_network_graph(pmbc_1103)
plot_network_graph(pmbc_1103) #grn


#checking centrality
pagerank_centrality_info=data.frame(NetworkGraph(pmbc_1103))

pagerank_centrality_info$name[order(-pagerank_centrality_info$centrality)][1:20]


GetGRN(pmbc_1103)
GetNetwork(pmbc_1103)
NetworkFeatures(pmbc_1103)
NetworkGraph(pmbc_1103)
NetworkModules(pmbc_1103)
NetworkParams(pmbc_1103)
NetworkTFs(pmbc_1103)
Params(pmbc_1103)
gof(pmbc_1103)
GetAssaySummary(pmbc_1103)








##########GRN ANALYSIS
######test to validate.
library(tidyr)
library(Seurat)
library(Signac)
library(SeuratWrappers)
library(EnsDb.Hsapiens.v86)
library(GenomicRanges)
fragpath <- "/corgi/sebas/tcell_multi/aggregated_tlibs/AGG6789_good/outs/atac_fragments.tsv.gz"
# get gene annotations for hg38
annotation <-GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- "UCSC"

genes_NKT <- readRDS("/corgi/martin/multiome_T/RDS_final/T_cells_noCC_DICE_ATAC_final_peaks to genes_NKT.RDS")
DefaultAssay(genes_NKT)='RNA'

# Create a copy of the original Seurat object for modification
modified_seurat_obj <- genes_NKT

# Set the seed for reproducibility
set.seed(123)

# Introduce random changes to gene expression and ATAC-seq peak values
gene_expression <- GetAssayData(modified_seurat_obj, assay = "RNA")
atac_seq_peaks <- GetAssayData(modified_seurat_obj, assay = "ATAC")

gene_expression_modified <- gene_expression + rnorm(n = length(gene_expression), mean = 0, sd = 0.1)
atac_seq_peaks_modified <- atac_seq_peaks + rnorm(n = length(atac_seq_peaks), mean = 0, sd = 0.1)


# create a Seurat object containing the RNA adata
modified_seurat_obj <- CreateSeuratObject(
  counts = gene_expression_modified,
  assay = "RNA"
)

# create ATAC assay and add it to the object

modified_seurat_obj[["ATAC"]] <- CreateChromatinAssay(
  counts = atac_seq_peaks_modified,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)


# Get the cell identities from genes_NKT
cell_identities <- Idents(genes_NKT)

# Assign the cell identities to the modified_seurat_obj
Idents(modified_seurat_obj) <- cell_identities


#testing for pando run

library(grr)
library(TFBSTools)
library(JASPAR2020)
library(Signac)
library(Pando)
#note: the VARIABLE FEATURES FOR RNA and VARIABLE PEAKS FOR ATAC need to be run with the same parameters as for the regular network
modified_seurat_obj<-FindVariableFeatures(modified_seurat_obj, 
                                          selection.method = 'mean.var.plot',
                                          binning.method = 'equal_frequency', 
                                          verbose = T,
                                          assay = 'RNA')
modified_seurat_obj<-FindVariableFeatures(modified_seurat_obj, 
                                          selection.method = 'mean.var.plot',
                                          binning.method = 'equal_frequency', 
                                          verbose = T,
                                          assay = 'ATAC')

genes_NKT<-FindVariableFeatures(genes_NKT, 
                                selection.method = 'mean.var.plot',
                                binning.method = 'equal_frequency', 
                                verbose = T,
                                assay = 'RNA')


# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

peaks <- CallPeaks(modified_seurat_obj, macs2.path = '/home/mahogny/miniconda3/bin/macs2'
                   )
# remove peaks on nonstandard chromosomes and in genomic blacklist regions
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38_unified, invert = TRUE)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(pmbc_1103),
  features = peaks,
  cells = colnames(pmbc_1103), 
  verbose = T
)

# create a new assay using the MACS2 peak set and add it to the Seurat object
pmbc_1103[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = fragpath,
  annotation = annotation
)

#####ATAC processing
DefaultAssay(pmbc_1103) <- "peaks"
pmbc_1103 <- RunTFIDF(pmbc_1103)
pmbc_1103 <- FindTopFeatures(pmbc_1103, min.cutoff = 3, verbose = T )
pmbc_1103 <- RunSVD(pmbc_1103, verbose = T, reduction.name = 'SVD_ATAC')
#calc dimensionality reduction
DepthCor(pmbc_1103, reduction = 'SVD_ATAC',n = 20)


####--Clustering --- Removing first component due to strong correlation of depth
pmbc_1103 <- RunUMAP(object = pmbc_1103, reduction = 'SVD_ATAC', dims = 2:30, reduction.name = 'UMAP_ATAC')
pmbc_1103<- FindNeighbors(object = pmbc_1103, reduction = 'SVD_ATAC', dims = 2:30)
pmbc_1103 <- FindClusters(object = pmbc_1103, verbose = FALSE, algorithm = 3)
DimPlot(object = pmbc_1103, pt.size = 1,group.by = 'labels.DICE')+ ggtitle('scATAC seq')
####3D UMAP more informative
pmbc_1103 <- RunUMAP(object = pmbc_1103, reduction = 'SVD_ATAC', dims = 2:30,n.components = 3L, reduction.name = 'X3D_UMAP_ATAC')

Th1_random<- subset(modified_seurat_obj,idents="T cells, CD4+, Th1")

modified_seurat_obj <- find_motifs(
  modified_seurat_obj,
  pfm = pfm,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

modified_seurat_obj<-initiate_grn(modified_seurat_obj)



Th1_random <- infer_grn(
  Th1_random,
  peak_to_gene_method = 'Signac',
  method = 'glm'
)


# Print inferred coefficients
coef(seurat_object)

# Find gene and regulatory modules 
test_srt <- find_modules(test_srt)

# Print modules
NetworkModules(test_srt)













install.packages('reticulate', dependencies = TRUE, INSTALL_opts = '--no-lock')
BiocManager::install('reticulate', force = T)
BiocManager::install("GenomeInfoDb", force = T)
install.packages("remotes")
library(remotes)
remotes::install_github("quadbio/Pando", ref = "v1.0.2")

#####----------------BACKBONE OF GENERAL 
library(Pando)
library(tidyr)
#### import networks
path_to_files="/corgi/martin/multiome_T/RDS_final/Pando_noCC/"
list.files('/corgi/martin/multiome_T/RDS_final/Pando_noCC/')
set.seed(666)
#----------import data
for (i in list.files(path_to_files)) {
  var_name <- paste0(sub("^[^_]+_([^.]+).*", "\\1", i)) # creates a variable name based on the current iteration
  file_path <- paste0("/corgi/martin/multiome_T/RDS_final/Pando_noCC/", i) # creates the path to the file based on the current iteration
  value <- readRDS(file_path) # reads the RDS file and assigns it to the value variable
  var_name=paste0(var_name,"_graph")
  assign(var_name, value) # assigns the value to the variable with the generated name
  remove(value)
}
#----------do things with data
# list all files in the current directory that contain "_graph"
file_list <- grep("_graph", ls(), value = TRUE)
library(dplyr)
library(Seurat)
# iterate through each file in the list
for (file in file_list) {
  #create it as a graph object
  var_name=paste0(file,"_obj")
  print(var_name)
  value=find_modules(get(file)) #default parameters: pval=0.05, rsq-threhs=0.1,.. 50 targets per
  value=get_network_graph(value)
  value=GetNetwork(value)
  value=NetworkGraph(value)
  print(value)
  assign(var_name,value)
  remove(var_name)
  remove(value)
  #optionally here we can input a command to extract the centrality for each type
}

value=find_modules(naive_graph) #default parameters: pval=0.05, rsq-threhs=0.1,.. 50 targets per
value=get_network_graph(value)
value=GetNetwork(value)
value=NetworkGraph(value)


#apply intersection
# create a list of objects using mget
library(igraph)
BiocManager::install("RCy3", force = T)
library(RCy3)

objects <- mget(grep("_obj", ls(), value = TRUE, fixed = T)) #SKIP TO COMMUNITIES HERE
objects<-readRDS('/corgi/martin/multiome_T/Pando_for_Sebastian_XGB/tbl_graphs.RDS') #as of 13 nov 2023 pando stopped working


###### random networks: check whether the centrality is due to random chance (p=0.5) or not, using watts-strogatz model (better modelling small networks like naive!). 


######NAIVE!
for (i in names(objects)) {
  stp = i #subtype name
  current_object = get(i)
  g=as.igraph(current_object)
  degree_distribution_values_g <- degree_distribution(g)
  degree_g<-degree(g)
  namecnt=sub("(.*)_graph_obj", "\\1", stp)
  print(namecnt)
  name1=paste0(namecnt,'Wdegdist_real')
  name2=paste0(namecnt,'Wdegcent_real')
  assign(name1,degree_distribution_values_g)
  assign(name2,degree_g)
  remove(g_degree_distribution_values_real)
  remove(name2)
  #done
  #simulate 1000 random graphs and get the matrix with that simulated deg distribution
  
  library(igraph)
  
  # Get the number of nodes in graph "g"
  n <- vcount(g)
  
  # Calculate the average number of first-degree neighbors in graph "g"
  avg_k <- mean(neighborhood.size(g, order=1,mode = "all"))
  
  # Calculate the average edge density
  avg_edge_density <- edge_density(g)
  
  
  # Create a list to store the random networks
  gl <- vector('list', 10000)
  
  for (i in 1:10000) {
    # Generate a random network based on the Watts-Strogatz model
    random_network <- watts.strogatz.game(1, n, avg_k, avg_edge_density)
    
    # Rename vertices to match the original network
    V(random_network)$name <- V(g)$name
    
    # Store the random network in the list
    gl[[i]] <- random_network
  }
  
  
  # Loop through the list of graphs and calculate degree dist for each
  gl_degdists<- unlist(
    lapply(gl, degree_distribution, mode=c("total"))
  ) 
  name2=paste0(namecnt,'Wdegdist_sim')
  assign(name2,gl_degdists)
  print(name2)
  remove(gl_degdists)
  remove(name2)
  
  gl_degcents<- unlist(
    lapply(gl, degree, mode=c("total"))
  ) 
  name2=paste0(namecnt,'Wdegcent_sim')
  assign(name2,gl_degcents)
  print(name2)
  remove(gl_degcents)
}

plot_density_overlay <- function(data1, data2) {
  # Create the density plot for data1 in green
  plot(density(data1), main = "Density Plot", col = "green", lwd=2)
  
  # Overlay the density plot for data2 in red
  lines(density(data2), col = "red", lwd=2)
}

plot_density_overlay(NKdegdist_real, NKdegdist_sim)

#perform permutation
permutation_test_and_plot <- function(real_deg_dist, sim_deg_dist, n_permutations = 100) {
  set.seed(666)
  real_var_name <- deparse(substitute(real_deg_dist))
  # Calculate observed test statistic (difference of means)
  observed_diff <- mean(real_deg_dist) - mean(sim_deg_dist)
  
  # Pool the data
  pooled_data <- c(real_deg_dist, sim_deg_dist)
  
  # Initialize variable to store permuted differences
  permuted_diffs <- numeric(n_permutations)
  
  # Perform permutations
  for (i in 1:n_permutations) {
    permuted_data <- sample(pooled_data)
    permuted_real <- permuted_data[1:length(real_deg_dist)]
    permuted_sim <- permuted_data[(length(real_deg_dist) + 1):length(pooled_data)]
    permuted_diffs[i] <- mean(permuted_real) - mean(permuted_sim)
  }
  
  # Calculate p-value
  p_value <- sum(abs(permuted_diffs) >= abs(observed_diff)) / n_permutations
  
  # Find common range for both datasets
  common_xlim = range(c(real_deg_dist, sim_deg_dist))
  
  # Create the histogram for simulated data first
  hist(sim_deg_dist, breaks = 20, col = rgb(1,0,0,0.5), xlim = common_xlim, freq = FALSE,
       xlab = "Degree", ylab = "Density", main = real_var_name)
  
  # Add histogram for real data with some transparency
  hist(real_deg_dist, breaks = 20, col = rgb(0,1,0,0.7), add = TRUE, freq = FALSE)
  
  # Add lines for the average measures with thicker lines (lwd = 2)
  abline(v = mean(real_deg_dist), col = "green", lty = 2, lwd = 2)
  abline(v = mean(sim_deg_dist), col = "red", lty = 2, lwd = 2)
  
  # Add legend with p-value and number of permutations
  legend_text <- c("Real", "Simulated", "Real Avg", "Sim Avg", paste0("p-value: ", round(p_value, 4)), paste0("Permutations: ", n_permutations))
  legend("topright", legend = legend_text, 
         col = c(rgb(0,1,0,0.5), rgb(1,0,0,0.5), "green", "red", "black", "black"), 
         pch = c(15, 15, NA, NA, NA, NA), lty = c(NA, NA, 2, 2, NA, NA),
         text.font = 2)
  
  return(p_value)
}
#plot
p_value <- permutation_test_and_plot(
  naivedegdist_real, naivedegdist_sim, 1000)
print(paste0("P-value: ", p_value))

######ALL OTHERs, using BARABASI!
library(igraph)

for (i in names(objects)) {
  stp = i #subtype name
  current_object = get(i)
  g = as.igraph(current_object)
  degree_distribution_values_g <- degree_distribution(g)
  degree_g <- degree(g)
  
  namecnt = sub("(.*)_graph_obj", "\\1", stp)
  print(namecnt)
  
  name1 = paste0(namecnt, 'degdist_real')
  name2 = paste0(namecnt, 'degcent_real')
  
  assign(name1, degree_distribution_values_g)
  assign(name2, degree_g)
  
  # Simulate 10000 random graphs using the Barabási-Albert model
  
  # Create a list to store the random networks
  gl <- vector('list', 10000)
  
  for (i in 1:10000) {
    gl[[i]] <- barabasi.game(n = gorder(g), m = floor(ecount(g) / gorder(g)), directed = TRUE)
    V(gl[[i]])$name <- V(g)$name
  }
  
  # Loop through the list of graphs and calculate degree dist for each
  gl_degdists <- unlist(lapply(gl, degree_distribution, mode=c("total")))
  name2 = paste0(namecnt, 'degdist_sim')
  assign(name2, gl_degdists)
  print(name2)
  
  gl_degcents <- unlist(lapply(gl, degree, mode=c("total")))
  name2 = paste0(namecnt, 'degcent_sim')
  assign(name2, gl_degcents)
  print(name2)
}


plot_density_overlay(NKdegdist_real, NKdegdist_sim)

permutation_test_and_plot(Th17degdist_real,Th17degdist_sim)

wilcox.test(naiveWdegdist_real, naiveWdegdist_sim, alternative = "two.sided")
wilcox.test(naive_stimdegdist_real, naive_stimdegcent_sim, alternative = "two.sided")
wilcox.test(Th1degdist_real, Th1degdist_sim, alternative = "two.sided")
wilcox.test(Th2degdist_real, Th2degdist_sim, alternative = "two.sided")
wilcox.test(Th17degdist_real, Th17degdist_sim, alternative = "two.sided")
wilcox.test(Th1_17degdist_real, Th1_17degdist_sim, alternative = "two.sided")
wilcox.test(Treg_naivedegdist_real, Treg_naivedegdist_sim, alternative = "two.sided")
wilcox.test(Treg_memdegdist_real, Treg_memdegdist_sim, alternative = "two.sided")


wilcoxon_test_and_plot <- function(real_deg_dist, sim_deg_dist) {
  # Perform the Wilcoxon rank sum test
  test_result <- wilcox.test(real_deg_dist, sim_deg_dist, alternative = "two.sided")
  p_value <- test_result$p.value
  
  real_var_name <- deparse(substitute(real_deg_dist))
  
  # Find common range for both datasets
  common_xlim = range(c(real_deg_dist, sim_deg_dist))
  
  # Create the histogram for simulated data first
  hist(sim_deg_dist, breaks = 10, col = rgb(1,0,0,0.5), xlim = common_xlim, freq = FALSE,
       xlab = "Degree", ylab = "Density", main = real_var_name)
  
  # Add histogram for real data with some transparency
  hist(real_deg_dist, breaks = 10, col = rgb(0,1,0,0.7), add = TRUE, freq = FALSE)
  
  # Add lines for the median values with thicker lines (lwd = 2)
  abline(v = median(real_deg_dist), col = "green", lty = 2, lwd = 2)
  abline(v = median(sim_deg_dist), col = "red", lty = 2, lwd = 2)
  
  # Add legend with p-value
  legend_text <- c("Real", "Simulated", "Real Median", "Sim Median", paste0("p-value: ", signif(p_value, digits = 4)))
  legend("topright", legend = legend_text, 
         col = c(rgb(0,1,0,0.5), rgb(1,0,0,0.5), "green", "red", "black"), 
         pch = c(15, 15, NA, NA, NA), lty = c(NA, NA, 2, 2, NA),
         text.font = 2)
  
  return(p_value)
}




#### DEGREE DISTRIBUTIONS

wilcoxon_test_and_plot(naive_stimdegdist_real,naive_stimdegdist_sim)
wilcoxon_test_and_plot(naiveWdegdist_real, naiveWdegdist_sim)
wilcoxon_test_and_plot(Th1degdist_real,Th1degdist_sim)
wilcoxon_test_and_plot(Th2degdist_real,Th2degdist_sim)
wilcoxon_test_and_plot(Th17degdist_real,Th17degdist_sim)
wilcoxon_test_and_plot(Th1_17degdist_real,Th1_17degdist_sim)
wilcoxon_test_and_plot(NKdegdist_real,NKdegdist_sim)
wilcoxon_test_and_plot(Treg_memdegdist_real,Treg_memdegdist_sim)
wilcoxon_test_and_plot(Treg_naivedegdist_real,Treg_naivedegdist_sim)

t.test(naiveWdegdist_real, naiveWdegdist_sim, alternative = "two.sided")
t.test(naive_stimdegdist_real,naive_stimdegdist_sim, alternative = "two.sided")
t.test(Th1degdist_real,Th1degdist_sim, alternative = "two.sided")
t.test(Th2degdist_real,Th2degdist_sim, alternative = "two.sided")
t.test(Th17degdist_real,Th17degdist_sim, alternative = "two.sided")
t.test(Th1_17degdist_real,Th1_17degdist_sim, alternative = "two.sided")
t.test(NKdegdist_real,NKdegdist_sim, alternative = "two.sided")
t.test(Treg_memdegdist_real,Treg_memdegdist_sim, alternative = "two.sided")
t.test(Treg_naivedegdist_real,Treg_naivedegdist_sim, alternative = "two.sided")

#### CENTRALITIES

wilcoxon_test_and_plot(naive_stimdegcent_real,naive_stimdegcent_sim)
wilcoxon_test_and_plot(naiveWdegcent_real, naiveWdegcent_sim)
wilcoxon_test_and_plot(Th1degcent_real,Th1degcent_sim)
wilcoxon_test_and_plot(Th2degcent_real,Th2degcent_sim)
wilcoxon_test_and_plot(Th17degcent_real,Th17degcent_sim)
wilcoxon_test_and_plot(Th1_17degcent_real,Th1_17degcent_sim)
wilcoxon_test_and_plot(NKdegcent_real,NKdegcent_sim)
wilcoxon_test_and_plot(Treg_memdegcent_real,Treg_memdegcent_sim)
wilcoxon_test_and_plot(Treg_naivedegcent_real,Treg_naivedegcent_sim)

wilcox.test(naiveWdegcent_real, naiveWdegcent_sim, alternative = "two.sided")
wilcox.test(naive_stimdegcent_real,naive_stimdegcent_sim, alternative = "two.sided")
wilcox.test(Th1degcent_real,Th1degcent_sim, alternative = "two.sided")
wilcox.test(Th2degcent_real,Th2degcent_sim, alternative = "two.sided")
wilcox.test(Th17degcent_real,Th17degcent_sim, alternative = "two.sided")
wilcox.test(Th1_17degcent_real,Th1_17degcent_sim, alternative = "two.sided")
wilcox.test(NKdegcent_real,NKdegcent_sim, alternative = "two.sided")
wilcox.test(Treg_memdegcent_real,Treg_memdegcent_sim, alternative = "two.sided")
wilcox.test(Treg_naivedegcent_real,Treg_naivedegcent_sim, alternative = "two.sided")




t.test(naiveWdegcent_real, naiveWdegcent_sim, alternative = "two.sided")
t.test(naive_stimdegcent_real,naive_stimdegcent_sim, alternative = "two.sided")
t.test(Th1degcent_real,Th1degcent_sim, alternative = "two.sided")
t.test(Th2degcent_real,Th2degcent_sim, alternative = "two.sided")
t.test(Th17degcent_real,Th17degcent_sim, alternative = "two.sided")
t.test(Th1_17degcent_real,Th1_17degcent_sim, alternative = "two.sided")
t.test(NKdegcent_real,NKdegcent_sim, alternative = "two.sided")
t.test(Treg_memdegcent_real,Treg_memdegcent_sim, alternative = "two.sided")
t.test(Treg_naivedegcent_real,Treg_naivedegcent_sim, alternative = "two.sided")


t.test(naiveWdegcent_real, naiveWdegcent_sim, alternative = "two.sided")
t.test(naive_stimWdegcent_real,naive_stimWdegcent_sim, alternative = "two.sided")
t.test(Th1Wdegcent_real,Th1Wdegcent_sim, alternative = "two.sided")
t.test(Th2Wdegcent_real,Th2Wdegcent_sim, alternative = "two.sided")
t.test(Th17Wdegcent_real,Th17Wdegcent_sim, alternative = "two.sided")
t.test(Th1_17Wdegcent_real,Th1_17Wdegcent_sim, alternative = "two.sided")
t.test(NKWdegcent_real,NKWdegcent_sim, alternative = "two.sided")
t.test(Treg_memWdegcent_real,Treg_memWdegcent_sim, alternative = "two.sided")
t.test(Treg_naiveWdegcent_real,Treg_naiveWdegcent_sim, alternative = "two.sided")

#t test works

t_test_and_plot <- function(real_deg_dist, sim_deg_dist) {
  # Perform the two-sided t-test
  test_result <- t.test(real_deg_dist, sim_deg_dist, alternative = "two.sided")
  p_value <- test_result$p.value
  
  # Calculate means and format them with 3 decimal places
  real_mean <- formatC(mean(real_deg_dist), format = "f", digits = 3)
  sim_mean <- formatC(mean(sim_deg_dist), format = "f", digits = 3)
  
  real_var_name <- deparse(substitute(real_deg_dist))
  
  # Find common range for both datasets
  common_xlim = range(c(real_deg_dist, sim_deg_dist))
  
  # Create the histogram for simulated data first
  hist(sim_deg_dist, breaks = 5, col = rgb(1,0,0,0.5), xlim = common_xlim, freq = FALSE,
       xlab = "Degree", ylab = "Density", main = real_var_name)
  
  # Add histogram for real data with some transparency
  hist(real_deg_dist, breaks = 20, col = rgb(0,1,0,0.7), add = TRUE, freq = FALSE)
  
  # Add lines for the mean values with thicker lines (lwd = 2)
  abline(v = mean(real_deg_dist), col = "green", lty = 2, lwd = 2)
  abline(v = mean(sim_deg_dist), col = "red", lty = 2, lwd = 2)
  
  # Add legend with p-value and means
  legend_text <- c("Real", "Simulated", 
                   paste("Real Mean:", real_mean), 
                   paste("Sim Mean:", sim_mean), 
                   paste0("p-value: ", signif(p_value, digits = 4)))
  legend("topright", legend = legend_text, 
         col = c(rgb(0,1,0,0.5), rgb(1,0,0,0.5), "green", "red", "black"), 
         pch = c(15, 15, NA, NA, NA), lty = c(NA, NA, 2, 2, NA),
         text.font = 2)
  
  return(p_value)
}



t_test_and_plot(naiveWdegcent_real, naiveWdegcent_sim)
t_test_and_plot(naive_stimWdegcent_real,naive_stimWdegcent_sim)
t_test_and_plot(Th1Wdegcent_real,Th1Wdegcent_sim)
t_test_and_plot(Th2Wdegcent_real,Th2Wdegcent_sim)
t_test_and_plot(Th17Wdegcent_real,Th17Wdegcent_sim)
t_test_and_plot(Th1_17Wdegcent_real,Th1_17Wdegcent_sim)
t_test_and_plot(NKWdegcent_real,NKWdegcent_sim)
t_test_and_plot(Treg_memWdegcent_real,Treg_memWdegcent_sim)
t_test_and_plot(Treg_naiveWdegcent_real,Treg_naiveWdegcent_sim)

t_test_and_plot(naiveWdegcent_real, naiveWdegcent_sim)
t_test_and_plot(naive_stimdegcent_real,naive_stimdegcent_sim)
t_test_and_plot(Th1degcent_real,Th1degcent_sim)
t_test_and_plot(Th2Wdegcent_real,Th2Wdegcent_sim)
t_test_and_plot(Th17degcent_real,Th17degcent_sim)
t_test_and_plot(Th1_17degcent_real,Th1_17degcent_sim)
t_test_and_plot(NKdegcent_real,NKdegcent_sim)
t_test_and_plot(Treg_memdegcent_real,Treg_memdegcent_sim)
t_test_and_plot(Treg_naivedegcent_real,Treg_naivedegcent_sim)









plot_density_overlay(naiveWdegdist_real, naiveWdegdist_sim)
plot_density_overlay(naive_stimWdegdist_real,naive_stimWdegdist_sim)
plot_density_overlay(Th1Wdegdist_real,Th1Wdegdist_sim)
plot_density_overlay(Th2Wdegdist_real,Th2Wdegdist_sim)
plot_density_overlay(Th17Wdegdist_real,Th17Wdegdist_sim)
plot_density_overlay(Th1_17Wdegdist_real,Th1_17Wdegdist_sim)
plot_density_overlay(NKWdegdist_real,NKWdegdist_sim)
plot_density_overlay(Treg_memWdegdist_real,Treg_memWdegdist_sim)
plot_density_overlay(Treg_naiveWdegdist_real,Treg_naiveWdegdist_sim)













### hypothesis testing formula: vertex_name, gl_centralities, g_centralities_real, p_value_threshold

library(ggplot2)

perform_hypothesis_test <- function(simulated, real) {
  # Calculate the global mean of gl_centralities
  global_mean_gl <- mean(gl_centralities)
  
  # Perform a one-sample t-test comparing the global mean to g_centralities_real
  t_test_result <- t.test(g_centralities_real, mu = global_mean_gl, alternative="two.sided", conf.level=0.95)
  
  # Create a histogram of g_centralities_real
  hist_plot <- ggplot(data.frame(Value = g_centralities_real), aes(x = g_centralities_real)) +
    geom_histogram(binwidth = 0.02, fill = "blue", color = "black") +
    
    # Add a vertical line for the global mean
    geom_vline(xintercept = global_mean_gl, color = "red", linetype = "dashed", size = 1) +
    
    # Add a vertical line for the p-value threshold
    geom_vline(xintercept = quantile(g_centralities_real, 1 - p_value_threshold), color = "green", linetype = "dashed", size = 1) +
    
    labs(x = "Eigenvector Centrality") +
    theme_minimal()
  
  # Print the plot
  print(hist_plot)
  
  # Print the hypothesis test result with p-value
  if (t_test_result$p.value <= p_value_threshold) {
    cat("Reject the null hypothesis: The means are significantly different (p =", t_test_result$p.value, ").\n")
  } else {
    cat("Fail to reject the null hypothesis: The means are not significantly different (p =",t_test_result$p.value, ").\n")
  }
}


perform_hypothesis_test(gl_centralities = Th17centralities_sim, g_centralities_real = Th17centralities_real, p_value_threshold = 0.05)

# Load the igraph library if not already loaded
library(igraph)

# Define a layout for the graph (e.g., Fruchterman-Reingold layout)
layout <- layout_with_fr(g)

# Calculate degree centrality
degree_centrality <- degree(g)

# Identify the top 8 most highly central vertices based on degree centrality
top_vertices <- names(sort(degree_centrality, decreasing = TRUE)[1:8])

# Identify the first-degree neighbors of the top vertices
neighbor_vertices <- unique(unlist(neighbors(g, top_vertices)))

# Create a plot of the graph with customizations
plot(g, 
     layout = layout,             
     main = "Customized Graph", 
     vertex.label = NA,  # Hide all vertex labels
     vertex.size = 10,            
     edge.arrow.size = 0.5,       
     edge.curved = 0.2,
     vertex.color = ifelse(V(g)$name %in% top_vertices | V(g)$name %in% neighbor_vertices, "red", "blue")) 

# Create a legend with corresponding vertex names
legend("topright", legend = top_vertices, col = "red", pch = 1, cex = 0.8)

# Add text to display the edge density value at the bottom of the figure
text(x = layout[, 1], 
     y = min(layout[, 2]) - 1, 
     labels = paste("Edge Density:", round(edge_density(g), 4)), 
     pos = 1, 
     col = "black")

# Adjust the plot layout
par(mar = c(1, 1, 3, 1))  # Increase bottom margin for edge density text

# Save or display the plot










####### backbone bois, degree and eigenvectors

###degree
graph_itersection= do.call(graph.intersection,objects)
#get vertex attributes
vertexatrs= data.frame(vertex_attr(graph_itersection))
rownames(vertexatrs)=vertexatrs$name
#select only centrality
vertexatrs=vertexatrs[grep('centrality',colnames(vertexatrs))]
vertexatrs$median_centrality=apply(vertexatrs, 1, median, na.rm=T)
vertexatrs$mean_centrality=apply(vertexatrs, 1, mean, na.rm=T)
#remove TFS with any NA occurence for centrality
vertexatrs<- na.omit(vertexatrs)

# creating a df with the centralities for each TF and each network, based on the consensus when aggregating all the networks together
common.tfs= rownames(vertexatrs)
#calculating the degree_centrality: why? it reflects the number of connections(edges) it has. Other options available.Try pagerank





colnames(degree_centrality_networks)[9]='naive_graph_obj_deg.cent'
write.csv(degree_centrality_networks,'/corgi/sebas/tcell_multi/backbone_degcent_all/all_backbonebois.csv')
degree_centrality_networks=degree_centrality_networks[ , !colnames(degree_centrality_networks) %in% c("degree_centrality")]
library(viridis)
library(reshape2)
library(ggridges)

degcent_2=degree_centrality_networks
degcent_2$gene=rownames(degcent_2)
degcent_2=degcent_2[1:8,] #plot for supplementary
degcent_2=melt(degcent_2)
degcent2_ordered <- degcent_2[order(-degcent_2$value),]


# Take the top 100 entries
degcent2_top <- head(degcent2_ordered, 100)

# Create a stacked bar plot for figure
ggplot(degcent2_top, aes(x = reorder(gene, value), y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  scale_fill_npg() +
  labs(x = "gene", y = "value") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 20, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 10, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )



###### plot for backbone
molten_degcent=melt(degcent_2, id.vars = 'gene')
molten_degcent=molten_degcent[order(molten_degcent$value, decreasing = T),]
molten_degcent$gene=factor(molten_degcent$gene,  levels=unique(molten_degcent$gene))
molten_degcent$variable <- gsub("_graph_obj_deg\\.cent", "", molten_degcent$variable)
####heatmap from figure
ggplot(molten_degcent, aes(x = gene, y = variable)) +
  geom_tile(aes(fill = value), color = "black", size = 0.5) +
  scale_fill_gradient(low = "cornsilk", high = "brown3") +
  labs(x = "Backbone gene", y = "Subtype") +
  guides(fill = guide_legend(override.aes = list(size = 3))) +
  theme(
    panel.grid = element_line(size = 1, linewidth = 0.5, colour = "black"),
    axis.text = element_text(face = "bold", size = 14),
    axis.title = element_text(face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_text(face = "bold"),
    legend.key.size = unit(1, "lines")
  )





degree_centrality_networks$TF<-rownames(degree_centrality_networks)
degree_centrality_networks=reshape2::melt(degree_centrality_networks)
degree_centrality_networks['variable']=gsub("_graph_obj", "", degree_centrality_networks$variable)
degree_centrality_networks= degree_centrality_networks %>% arrange(desc(degree_centrality_networks$value))
#write.csv(degree_centrality_networks,'degcent_networks.csv')
library(ggplot2)
ggplot(degree_centrality_networks[1:20,], aes(x = variable, y = TF), filled.contour()) +
  geom_tile(aes(fill=value)) + scale_fill_gradient(low = "blue", high = "red")


#degree centrality for each one of them. 

for (i in names(objects)) {
  colname <- i
  i <- get(i)
  degree_centrality <- degree(i, v=V(i))
  df <- data.frame(degree_centrality)
  print(df)
  write.csv(df, file = paste0('/corgi/sebas/tcell_multi/backbone_degcent_all/',colname,"_degree_centrality.csv"))
}


########eigenvector centrality
# Create an empty matrix to store eigenvector centrality values
eigen_centrality <- matrix(NA, nrow = length(common.tfs), ncol = length(objects))
rownames(eigen_centrality) <- common.tfs
colnames(eigen_centrality) <- names(objects)

# Loop through each graph and compute eigenvector centrality for common.tfs
for (i in names(objects)) {
  current_object <- objects[[i]]
  centrality_values <- eigen_centrality(current_object)$vector[common.tfs]
  eigen_centrality[, i] <- centrality_values
}
eigen_centrality=as.data.frame(eigen_centrality)
# Take the top 100 entries
eigen_centrality_melted=eigen_centrality
eigen_centrality_melted$gene=rownames(eigen_centrality)
eigen_centrality_melted <- melt(eigen_centrality_melted,id.vars = 'gene')
eigen_centrality_top <- eigen_centrality_melted[order(-eigen_centrality_melted$value), ]
eigen_centrality_top=head(eigen_centrality_top,75)
eigen_centrality_top$variable <- sub("_graph_obj$", "", eigen_centrality_top$variable)

# Create a stacked bar plot
library(ggsci)
ggplot(eigen_centrality_melted, aes(x = reorder(gene, value), y = value, fill = variable)) +
  geom_bar(stat = "identity") +
  scale_fill_npg() +
  labs(x = "gene", y = "value") +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 20, face = "bold", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 20, face = "bold"),  # Increase the size of legend labels
    legend.title = element_text(size = 10, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )


#backbones
#not backbones



merged_cents= read.csv('/corgi/martin/multiome_T/csv_files/Sebastian_network/merged_graph_obj_degree_centrality.csv', header = T)
merged_cents=merged_cents[,-1]
rownames(merged_cents)=merged_cents[,1]
merged_cents=merged_cents[,-1]

backbones=merged_cents[rownames(merged_cents) %in% backbone.bois, ]
activation= merged_cents[merged_cents$naive==0,]
stimulation= merged_cents[merged_cents$naive==0 & merged_cents$naive.stim==0,]
effector=stimulation[rowSums(stimulation != 0) == 1, ]

####---- for cytoscape plotting
#Activation
sorted_stim <- activation[order(-activation$naive.stim),]
sorted_stim<-sorted_stim[1:6,]
sorted_activation<-rownames(sorted_activation)
write.csv(sorted_activation,'/corgi/sebas/tcell_multi/cytoscape/activation_genes.csv', row.names = F, quote = F)

#stimulation
sorted_stim <- stimulation[order(-activation$Th2),]
sorted_stim<-sorted_stim[1:6,]
sorted_stim= rownames(sorted_stim)
write.csv(sorted_stim,'/corgi/sebas/tcell_multi/cytoscape/stimulation_genes.csv', row.names = F, quote = F)

#effector

sorted_effectors_all_st <- character()

for (i in 1:9) {###overlapping of effectors, take top 3 for each subtype
  effectors <- head(rownames(effector)[which(effector[, i] > 0)][order(effector[effector[, i] > 0, i], decreasing = TRUE)], 3)
  sorted_effectors_all_st <<- c(sorted_effectors_all_st, effectors)
}

sorted_effectors_all_st ### input for overlapping of effector

          
write.csv(head(rownames(effector)[which(effector[, 1] > 0)][order(effector[effector[, 1] > 0, 1], decreasing = T)], 6),'/corgi/sebas/tcell_multi/cytoscape/effector_th1_genes.csv', row.names = F, quote = F)
write.csv(head(rownames(effector)[which(effector[, 2] > 0)][order(effector[effector[, 2] > 0, 2], decreasing = T)], 6),'/corgi/sebas/tcell_multi/cytoscape/effector_th2_genes.csv', row.names = F, quote = F)
write.csv(head(rownames(effector)[which(effector[, 3] > 0)][order(effector[effector[, 3] > 0, 3], decreasing = T)], 6),'/corgi/sebas/tcell_multi/cytoscape/effector_th1_17_genes.csv', row.names = F, quote = F)
write.csv(head(rownames(effector)[which(effector[, 4] > 0)][order(effector[effector[, 4] > 0, 4], decreasing = T)], 6),'/corgi/sebas/tcell_multi/cytoscape/effector_th17_genes.csv', row.names = F, quote = F)
write.csv(head(rownames(effector)[which(effector[, 5] > 0)][order(effector[effector[, 5] > 0, 5], decreasing = T)], 6),'/corgi/sebas/tcell_multi/cytoscape/effector_NKT_genes.csv', row.names = F, quote = F)
write.csv(head(rownames(effector)[which(effector[, 6] > 0)][order(effector[effector[, 6] > 0, 6], decreasing = T)], 6),'/corgi/sebas/tcell_multi/cytoscape/effector_naive.stim_genes.csv', row.names = F, quote = F)
write.csv(head(rownames(effector)[which(effector[, 7] > 0)][order(effector[effector[, 7] > 0, 7], decreasing = T)], 6),'/corgi/sebas/tcell_multi/cytoscape/effector_Treg.mem_genes.csv', row.names = F, quote = F)
write.csv(head(rownames(effector)[which(effector[, 8] > 0)][order(effector[effector[, 8] > 0, 8], decreasing = T)], 6),'/corgi/sebas/tcell_multi/cytoscape/effector_Treg.naive_genes.csv', row.names = F, quote = F)
write.csv(head(rownames(effector)[which(effector[, 9] > 0)][order(effector[effector[, 9] > 0, 9], decreasing = T)], 6),'/corgi/sebas/tcell_multi/cytoscape/effector_naive_genes.csv', row.names = F, quote = F)



df_top <- stimulation %>% filter(rowSums(.) > mean(rowSums(.))) %>%
  arrange(desc(mean(rowSums(.)))) %>% arrange(desc(rowSums(.))) %>%
  slice(1:20)
df_top$gene=rownames(df_top)

df_top=melt(backbones)
backbones$gene=rownames(backbones)
df_top=melt(backbones)

ggplot(df_top, aes(x = variable, y = gene)) +
  geom_tile(aes(fill = value), color = "black", size = 0.5) +
  scale_fill_gradient(low = "cornsilk", high = "brown3") +
  labs(y = "Backbone TF", x = 'Subtype', fill = "Degree Centrality") +
  guides(fill = guide_legend(override.aes = list(size = 3))) +
  theme(panel.grid = element_line(size = 1, linewidth = 0.5, colour = 'black'),
        axis.text = element_text(face = "bold", size = 14),
        axis.title = element_text(face = "bold", size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold")  # Set legend text to be bold
  )







#########compute global overlap
library(igraph)

# Function to calculate overlap coefficient using Szymkiewicz-Simpson similarity
overlap_coefficient <- function(set1, set2) {
  intersection_size <- length(intersect(set1, set2))
  set1_size <- length(set1)
  set2_size <- length(set2)
  return(intersection_size / sqrt(set1_size * set2_size))
}

jaccard_index <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))
  union <- length(union(set1, set2))
  
  jaccard <- intersection / union
  return(jaccard)
}




# Create an empty matrix to store overlap coefficients
n_networks <- length(objects)
global_overlap_networks <- matrix(NA, nrow = n_networks, ncol = n_networks)

# Loop through each pair of networks
for (i in 1:n_networks) {
  for (j in 1:n_networks) {
    # Get the gene sets for the current pair of networks
    set1 <- V(objects[[i]])
    set2 <- V(objects[[j]])
    
    # Calculate overlap coefficient and store in the matrix
    global_overlap_networks[i, j] <- jaccard_index(set1, set2)
  }
}

# Assign row and column names to overlap matrix
dimnames(global_overlap_networks) <- list(names(objects), names(objects))

# Print the overlap matrix
print(global_overlap_networks)

library(ggplot2)
library(reshape2)

# Convert the matrix into a data frame
df <- melt(global_overlap_networks)
df <- melt(global_overlap_networks)
df$Var1 <- gsub("_graph_obj", "", rownames(global_overlap_networks))

ggplot(df, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "cornsilk", high = "brown3") +
  geom_text(aes(label = sprintf("%.2f", value)), size = 6, color = "black", fontface = "bold") +
  labs(x = "Group", y = "Group", title = "Global Overlap Net") +
  theme_minimal() +
  theme(
    legend.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text.x = element_text(face = "bold", size = 15),
    axis.text.y = element_text(face = "bold", size = 15))








#----------export .csv for:
#A) TOP is defined as more than the mean centrality
#extract common TFs.coms should be filtered by mean centrality
col_mean <- mean(vertexatrs$mean_centrality)
morethanmean_centrality= subset(vertexatrs, vertexatrs$mean_centrality>col_mean)
lessthanmean_centrality= subset(vertexatrs, vertexatrs$mean_centrality<col_mean)
#export as csv to be done

#B) TOP SPECIFIC FOR EACH NETWORK 

#take the ones absent in vertexatrs object, for each subtype
#create result_list ordered by degree centrality
result_list <- list()
for (i in 1:length(objects)) {
  g <- objects[[i]]
  degree_centrality <- degree(g)
  vertex_names <- V(g)$name
  if (is.na(vertex_names[1])) vertex_names <- 1:vcount(g)
  vertex_names_to_remove_index <- which(vertex_names %in% rownames(vertexatrs))
  vertex_names_to_keep_index <- setdiff(1:vcount(g), vertex_names_to_remove_index)
  g_filtered <- induced.subgraph(g, vertex_names_to_keep_index)
  degree_centrality_filtered <- degree(g_filtered)
  result <- data.frame(vertex_name=vertex_names[vertex_names_to_keep_index], centrality_score=degree_centrality_filtered)
  result_list[[i]] <- result[order(result$centrality_score, decreasing=TRUE), ]
  names(result_list)[i] <- names(objects)[i]
}
#export result_list as csv if export needed

####$$$$
# Initialize an empty data frame
result_df <- data.frame(ListName = character(),
                        VertexName = character(),
                        CentralityScore = numeric(),
                        stringsAsFactors = FALSE)

# Loop through the list and extract the data
for (i in 1:length(result_list)) {
  list_name <- names(result_list)[i]
  sub_df <- result_list[[i]]
  sub_df$ListName <- list_name
  
  # Combine the data into the result data frame
  result_df <- rbind(result_df, sub_df)
}

# Rename the columns
colnames(result_df) <- c("Vertex", "CentralityScore", "Subtype")

# Display the result data frame
head(result_df)


####$$$$


library(purrr)
library(dplyr)
library(ggplot2)
# extract the top 5 rows of each dataframe
df_list_top5 <-  result_df %>%
  group_by(Subtype) %>%
  arrange(desc(CentralityScore)) %>%
  slice_head(n = 3)

df_list_top5 <- result_df[result_df$Vertex %in% df_list_top5$Vertex,]



# Display the result
df_list_top5

df_list_top5m <- dcast((df_list_top5), Subtype  ~ Vertex, value.var = 'CentralityScore', fill = 0)
df_list_top5 <- melt(df_list_top5m, value.name = "CentralityScore")
df_list_top5$logcentrality=log(1+df_list_top5$CentralityScore)
# Find maximum value within each category
max_vals <- tapply(df_list_top5$CentralityScore, max)

#plot
ggplot(df_list_top5, aes(x = Subtype, y = variable)) +
  geom_tile(aes(fill=CentralityScore), color = "black", size = 0.5) + scale_fill_gradient(low = "cornsilk", high = "brown3")+
  labs(y = "Driver gene")+
  guides(fill = guide_legend(override.aes = list(size = 3)))+
  theme(panel.grid = element_line(size = 1,linewidth = 0.5, colour = 'black'),
        axis.text = element_text(face = "bold", size = 14),
        axis.title = element_text(face = "bold", size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(face = "bold"),
        legend.key.size = unit(1, "lines")
  )


#####--------------------------------------MODULARITY

#### Function to test out all possible community detection algorithms of the igraph package
library(igraph)

detect_communities <- function(graph, algorithm) {
  communities <- NULL
  
  if (algorithm == "walktrap") {
    # Walktrap algorithm
    wt <- walktrap.community(graph)
    communities <- as.list(wt)
  } else if (algorithm == "fastgreedy") {
    # Fastgreedy algorithm
    fg <- fastgreedy.community(graph)
    communities <- as.list(fg)
  } else if (algorithm == "louvain") {
    # Louvain algorithm
    lv <- cluster_louvain(graph)
    communities <- as.list(lv)
  } else if (algorithm == "labelpropagation") {
    # Label Propagation algorithm
    lp <- label.propagation.community(graph)
    communities <- as.list(lp)
  } else if (algorithm == "infomap") {
    # Infomap algorithm
    im <- cluster_infomap(graph)
    communities <- as.list(im)
  } else if (algorithm == "leadingeigenvector") {
    # Leading Eigenvector algorithm
    le <- leading.eigenvector.community(graph)
    communities <- as.list(le)
  } else {
    # Invalid algorithm
    stop("Invalid community detection algorithm specified.")
  }
  
  return(communities)
}



#Find all objects with "CLUSTER_WALKTRAP" in their name
graphs#graph object
info_cluster_list <- grep("CLUSTER_INFOMAP", ls(), value = TRUE) #communities object

#INFOMAP
# Assuming 'graphs' is a list of your graph objects
par(mfrow = c(3, 3))  # Set up a 3x3 plot grid

for (i in names(graphs)) {
  # Get the graph object by its name
  g <- graphs[[i]]
  
  # Compute the infomap community of the graph
  com <- infomap.community(g)
  
  # Get the sizes of the communities
  community_sizes <- sizes(com)
  
  # Plot the size distribution of the communities
  hist(community_sizes, main = paste("Community Sizes of", i), 
       xlab = "Size", ylab = "Frequency", col = "lightblue", breaks = 7)
  
  # Optionally assign the community list to a variable in the global environment
  assign(paste0('INFOMP_', i), as.list(com))
  print(table(community_sizes))
}
#LABELPROP
# Assuming 'graphs' is a list of your graph objects
par(mfrow = c(3, 3))  # Set up a 3x3 plot grid
library(igraph)
for (i in names(graphs)) {
  # Get the graph object by its name
  print(i)
  g <- graphs[[i]]
  
 
  
  # Compute communities using label propagation with weights
  com <- label.propagation.community(g, weights = E(g)$weight)
  
  # Get the sizes of the communities
  community_sizes <- sizes(com)
  
  # Plot the size distribution of the communities
  hist(community_sizes, main = paste("Community Sizes of", i), 
       xlab = "Size", ylab = "Frequency", col = "lightblue", breaks = 7)
  
  # Optionally assign the community list to a variable in the global environment
  assign(paste0('LABELPROP_', i), as.list(com))
  print(table(community_sizes))
}

#calculate random walk trap community
#the assumption is that better connected node clusters will retain the random generated walk within
#that cluster (ITS A TRAP!), just like a good ol' Ikea building is designed.
library(igraph)
graphs=objects
remove(objects)
pula_list <- names(graphs)
com=NULL
pula=NULL
# Iterate over the list of objects
for (i in names(graphs)) {
  # Get the object by its name
  pula <- graphs[[i]]
  # Compute the walktrap community of the pula object
  com <-walktrap.community(pula)
  com=as.list(com)
  #print(com)
  assign(paste0('CLUSTER_WALKTRAP_',i),com)
  
}


#Find all objects with "CLUSTER_WALKTRAP" in their name
graphs#graph object
#cluster_list <- grep("LABELPROP_", ls(), value = TRUE) #communities object
cluster_list <- grep("CLUSTER_WALKTRAP", ls(), value = TRUE)
library(gprofiler2)

###GLOBAL ENRICHMENT
perform_GO_enrichment <- function(community_obj_name, organism = "hsapiens") {
  community_obj <- get(community_obj_name)
  membership <- membership(community_obj)
  
  results_list <- list()
  
  for (group_id in unique(membership)) {
    group_members <- names(which(membership == group_id))
    
    if (length(group_members) > 10) {
      gostres <- try(gost(query = group_members,
                          exclude_iea = TRUE,
                          user_threshold = 0.05,
                          correction_method = "fdr",
                          domain_scope = 'annotated',
                          sources = c('GO:BP'),
                          organism = organism), silent = TRUE)
      
      if (inherits(gostres, "try-error")) {
        next
      } else {
        gprofout <- gostres$result[c("p_value", 'term_id', "term_name")]
        
        if (!is.null(gprofout) && nrow(gprofout) > 0) {
          gprofout <- gprofout[gprofout$p_value < 0.05,]
          if (nrow(gprofout) > 0) {
            gprofout <- gprofout[1, ]
            results_list[[paste0(community_obj_name, "_Group", group_id)]] <- gprofout
          }
        }
      }
    }
  }
  
  return(results_list)
}

all_results_GOBP <- lapply(cluster_list, perform_GO_enrichment)
all_results_GOMF <- lapply(cluster_list, perform_GO_enrichment)



############### Backbones!!!!
# Backbone genes
backbone_bois <- c('BCL6','CREM', 'HMGB2', 'MAF', 'PBX4','RBPJ', 'STAT4', 'ZNF292')
perform_GO_enrichment_with_backbone_presence <- function(community_obj_name, backbone_genes, organism = "hsapiens") {
  # If backbone_genes is a list of object names, fetch and combine genes from all driver objects
  if (is.list(backbone_genes)) {
    backbone_genes <- unique(unlist(lapply(backbone_genes, function(obj_name) get(obj_name))))
  }
  
  # Fetch the community object
  community_obj <- get(community_obj_name)
  membership <- membership(community_obj)
  
  results_list <- list()
  
  for (group_id in unique(membership)) {
    group_members <- names(which(membership == group_id))
    
    # Identify backbone/driver genes present in the group
    present_backbone_genes <- intersect(group_members, backbone_genes)
    
    # Check if the community contains any backbone genes
    if (length(present_backbone_genes) > 0) {
      if (length(group_members) > 10) {
        gostres <- try(gost(query = group_members,
                            exclude_iea = TRUE,
                            user_threshold = 0.05,
                            correction_method = "fdr",
                            domain_scope = 'annotated',
                            sources = c('GO:MF'),
                            organism = organism), silent = TRUE)
        
        if (inherits(gostres, "try-error")) {
          next
        } else {
          gprofout <- gostres$result[c("p_value", 'term_id', "term_name")]
          
          if (!is.null(gprofout) && nrow(gprofout) > 0) {
            gprofout <- gprofout[gprofout$p_value < 0.05,]
            if (nrow(gprofout) > 0) {
              gprofout <- gprofout[1, ]
              results_list[[paste0(community_obj_name, "_Group", group_id)]] <- list(
                GO_Results = gprofout,
                Genes = group_members,
                Present_Backbone_Drivers = present_backbone_genes
              )
            }
          }
        }
      }
    }
  }
  
  return(results_list)
}
all_results_with_backbone_presence<-lapply(cluster_list, perform_GO_enrichment_with_backbone_presence, backbone_bois)


############### Drivers!!!!
#driver genes
NKT_drivers <- c("ZNF732", "RPS8", "KLRB1", "GNPTAB")
Th1_drivers <- c("EOMES", "CCL5", "CENPE", "SERPINB9", "PAM", "SSBP2", "AHI1", "DTNB")
Th2_drivers <- c("KLF2", "ABI1", "ARHGEF3")
Th17_drivers <- c("HLF", "PVT1", "HMGCR", "MAST4", "LRRFIP2", "EDEM3", "RANBP9", "MAP3K4", "GABPB1", "SLC16A7", "SLC9A9", "PHLPP1", "PBX1")
Th1_17_drivers <- c("HLF", "EOMES", "GMNN", "ZC3HAV1")
mTreg_drivers <- c("SOX13", "PLIN2", "TXNIP", "BCAS3", "TMEM260")
nTreg_drivers <- c("TFEC", "DACH1", "ESR1", "DUSP4", "ST6GALNAC3", "CARMIL1")
naive.stim_drivers <- c("MYB", "CFLAR", "MBD5", "STAMBPL1", "HEATR5A", "KCNQ5")

list_drivers<-grep("_drivers", ls(), value = TRUE)
perform_GO_enrichment_with_driver_presence <- function(community_obj_name, driver_obj_names, organism = "hsapiens") {
  # Fetch and combine genes from all driver objects
  backbone_genes <- unique(unlist(lapply(driver_obj_names, function(obj_name) get(obj_name))))
  
  # Rest of the function remains the same
  community_obj <- get(community_obj_name)
  membership <- membership(community_obj)
  
  results_list <- list()
  
  for (group_id in unique(membership)) {
    group_members <- names(which(membership == group_id))
    
    # Identify backbone/driver genes present in the group
    present_backbone_genes <- intersect(group_members, backbone_genes)
    
    # Check if the community contains any backbone genes
    if (length(present_backbone_genes) > 0) {
      if (length(group_members) > 10) {
        gostres <- try(gost(query = group_members,
                            exclude_iea = TRUE,
                            user_threshold = 0.05,
                            correction_method = "fdr",
                            domain_scope = 'annotated',
                            sources = c('GO:MF'),
                            organism = organism), silent = TRUE)
        
        if (inherits(gostres, "try-error")) {
          next
        } else {
          gprofout <- gostres$result[c("p_value", 'term_id', "term_name")]
          
          if (!is.null(gprofout) && nrow(gprofout) > 0) {
            gprofout <- gprofout[gprofout$p_value < 0.05,]
            if (nrow(gprofout) > 0) {
              gprofout <- gprofout[1, ]
              results_list[[paste0(community_obj_name, "_Group", group_id)]] <- list(
                GO_Results = gprofout,
                Genes = group_members,
                Present_Backbone_Drivers = present_backbone_genes
              )
            }
          }
        }
      }
    }
  }
  
  return(results_list)
}
# Applying the function to each community object
all_results_with_driver_presence <- lapply(cluster_list, perform_GO_enrichment_with_driver_presence, list_drivers)



###### function for martin's sanity
extract_info_for_gene <- function(nested_list, target_gene) {
  extracted_info <- data.frame(Group = character(),
                               GO_Results = I(list()),
                               Genes = I(list()),
                               Present_Backbone_Drivers = I(list()),
                               stringsAsFactors = FALSE)
  
  for (outer_list in nested_list) {
    for (inner_list_name in names(outer_list)) {
      inner_list <- outer_list[[inner_list_name]]
      if (target_gene %in% inner_list$Genes) {
        # Append to the dataframe
        extracted_info <- rbind(extracted_info, data.frame(Group = inner_list_name,
                                                           GO_Results = I(list(inner_list$GO_Results)),
                                                           Genes = I(list(inner_list$Genes)),
                                                           Present_Backbone_Drivers = I(list(inner_list$Present_Backbone_Drivers))))
      }
    }
  }
  
  return(extracted_info)
}




a=extract_info_for_gene(all_results_with_backbone_presence,'HMGB2')
a











#### SS VALUE
# Function to calculate Szymkiewicz-Simpson coefficient
szymkiewicz_simpson <- function(set1, set2) {
  if (length(set1) == 0 || length(set2) == 0) return(0)
  return(length(intersect(set1, set2)) / min(length(set1), length(set2)))
}

# Extract and compare genes from all groups
extract_and_compare_genes_v1 <- function(all_results) {
  group_genes <- list()
  group_drivers <- list()
  
  # Extract genes and backbone drivers from each group
  for (result_list in all_results) {
    for (group_name in names(result_list)) {
      group_genes[[group_name]] <- result_list[[group_name]]$Genes
      group_drivers[[group_name]] <- result_list[[group_name]]$Present_Backbone_Drivers
    }
  }
  
  group_names <- names(group_genes)
  # Modify group names and add backbone drivers as suffix
  modified_group_names <- sapply(group_names, function(name) {
    # Remove specific prefixes and suffixes
    short_name <- gsub("^CLUSTER_WALKTRAP_", "", name)
    short_name <- sub("_Group[0-9]+.*$", "", short_name)
    # Concatenate multiple drivers with underscores and add as suffix
    driver_suffix <- ifelse(is.null(group_drivers[[name]]) || group_drivers[[name]] == "", "", paste0("_", paste(group_drivers[[name]], collapse = "_")))
    return(paste0(short_name, driver_suffix))
  })
  
  # Compute overlap matrix
  overlap_matrix <- matrix(nrow = length(group_genes), ncol = length(group_genes), dimnames = list(modified_group_names, modified_group_names))
  for (i in seq_along(group_genes)) {
    for (j in seq_along(group_genes)) {
      overlap_matrix[i, j] <- szymkiewicz_simpson(group_genes[[group_names[i]]], group_genes[[group_names[j]]])
    }
  }
  
  return(overlap_matrix)
}

overlap_matrix <- extract_and_compare_genes_v1(all_results_with_backbone_presence)


##### function for martins sanity v2- compute overlap on all the communities, irrespective of molecular function annotated or not
extract_groups <- function(community_obj_name, backbone_genes) {
  community_obj <- get(community_obj_name)
  membership <- membership(community_obj)
  group_genes <- list()
  
  for (group_id in unique(membership)) {
    group_members <- names(which(membership == group_id))
    if (length(group_members) > 10 && length(intersect(group_members, backbone_genes)) > 0) {
      group_genes[[paste0(community_obj_name, "_Group", group_id)]] <- group_members
    }
  }
  
  return(group_genes)
}

# Apply this to each community object and combine the results
all_group_genes <- lapply(cluster_list, extract_groups, backbone_bois)
all_group_genes <- do.call(c, all_group_genes)

group_names <- names(all_group_genes)
overlap_matrix <- matrix(0, nrow = length(group_names), ncol = length(group_names), dimnames = list(group_names, group_names))

for (i in 1:length(group_names)) {
  for (j in 1:length(group_names)) {
    if (i != j) {
      overlap_matrix[i, j] <- szymkiewicz_simpson(all_group_genes[[i]], all_group_genes[[j]])
    }
  }
}

create_new_name_with_drivers <- function(group_name, all_group_genes, backbone_genes, list_drivers) {
  group_backbones <- intersect(all_group_genes[[group_name]], backbone_genes)
  group_size <- length(all_group_genes[[group_name]])
  
  # Extracting the subtype from the group name
  subtype <- sub("CLUSTER_WALKTRAP_", "", sub("_graph_obj_Group[0-9]+", "", group_name))
  
  # Fetching the drivers for the subtype
  driver_obj_name <- paste0(subtype, "_drivers")
  if (driver_obj_name %in% list_drivers) {
    group_drivers <- get(driver_obj_name)
  } else {
    group_drivers <- NULL
  }
  
  # Finding the intersection of group genes with drivers
  present_drivers <- intersect(all_group_genes[[group_name]], group_drivers)
  
  # Creating the name
  combined_genes <- unique(c(group_backbones, present_drivers))
  if (length(combined_genes) > 0) {
    return(paste0(subtype, "_", paste(combined_genes, collapse = "_"), "_", group_size))
  } else {
    return(paste0(subtype, "_", group_size))
  }
}

# Assuming list_drivers is a character vector with driver names
new_names_with_drivers <- sapply(names(all_group_genes), create_new_name_with_drivers, all_group_genes, backbone_bois, list_drivers)

# Assuming list_drivers is a character vector with driver names
new_names_with_drivers <- sapply(names(all_group_genes), create_new_name_with_drivers, all_group_genes, backbone_bois, list_drivers)
rownames(overlap_matrix) <- new_names_with_drivers
colnames(overlap_matrix) <- new_names_with_drivers








library(ggplot2)
library(reshape2)

# Transform the matrix into a long format
melted_overlap_matrix <- melt(overlap_matrix)

ggplot(melted_overlap_matrix, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0.5, limit = c(0, 1), space = "Lab", 
                       name="Overlap\nCoefficient") +
  theme_minimal() +
  theme(axis.text.x = element_text( angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text( vjust = 1, face = "bold"),
        axis.title = element_blank(),
        legend.title.align = 0.5) +
  coord_fixed()










#check if drivers and backbones are in the same community, if yes, get GO results and components. List1 is drivers, list 2 is backbones

####Overlapping values, community info: which drivers, which backbones, what molecular function, what pvalue
extract_and_compare_genes <- function(all_results_with_driver_presence, all_results_with_backbone_presence) {
  group_genes <- list()
  group_drivers <- list()
  group_go_terms <- list()
  
  # Extract genes, backbone drivers, and GO results from each group in backbone presence
  for (result_list in all_results_with_backbone_presence) {
    for (group_name in names(result_list)) {
      group_genes[[group_name]] <- result_list[[group_name]]$Genes
      group_drivers[[group_name]] <- result_list[[group_name]]$Present_Backbone_Drivers
      go_results <- result_list[[group_name]]$GO_Results
      if(nrow(go_results) > 0){
        group_go_terms[[group_name]] <- paste(go_results$term_id[1])
      } else {
        group_go_terms[[group_name]] <- NA
      }
    }
  }
  
  # Extract genes of interest from driver presence
  genes_of_interest <- unique(unlist(lapply(all_results_with_driver_presence, function(x) unlist(lapply(x, function(y) y$Present_Backbone_Drivers)))))
  
  group_names <- names(group_genes)
  # Modify group names to include backbone drivers, genes of interest, and GO terms
  modified_group_names <- sapply(group_names, function(name) {
    # Remove specific prefixes and suffixes
    short_name <- gsub("^CLUSTER_WALKTRAP_", "", name)
    short_name <- sub("_Group[0-9]+.*$", "", short_name)
    
    # Concatenate multiple drivers with underscores and add as suffix
    driver_suffix <- paste0("_", paste(unique(c(group_drivers[[name]], genes_of_interest[genes_of_interest %in% group_genes[[name]]])), collapse = "_"))
    
    # Add GO term and p-value to the suffix
    go_suffix <- ifelse(!is.na(group_go_terms[[name]]), paste0("_", group_go_terms[[name]]), "")
    
    return(paste0(short_name, driver_suffix, go_suffix))
  })
  
  # Compute overlap matrix
  overlap_matrix <- matrix(nrow = length(group_genes), ncol = length(group_genes), dimnames = list(modified_group_names, modified_group_names))
  for (i in seq_along(group_genes)) {
    for (j in seq_along(group_genes)) {
      overlap_matrix[i, j] <- szymkiewicz_simpson(group_genes[[group_names[i]]], group_genes[[group_names[j]]])
    }
  }
  
  return(overlap_matrix)
}

# Use the function with your data
overlap_matrix2 <- extract_and_compare_genes(all_results_with_driver_presence, all_results_with_backbone_presence)


overlap_matrix2 <- extract_and_compare_genes(all_results_with_driver_presence, all_results_with_backbone_presence)

library(ggplot2)
library(reshape2)

# Transform the matrix into a long format
melted_overlap_matrix <- melt(overlap_matrix2)

write.csv(melted_overlap_matrix,'/corgi/sebas/tcell_multi/overlap_matrix_backbones_drivers/overlap_mtx_driv_bb_MF.csv')

ggplot(melted_overlap_matrix, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0.5, limit = c(0, 1), space = "Lab", 
                       name="Overlap\nCoefficient") +
  theme_minimal() +
  theme(axis.text.x = element_text( angle = 45, hjust = 1, face = "bold"),
        axis.text.y = element_text( vjust = 1, face = "bold"),
        axis.title = element_blank(),
        legend.title.align = 0.5) +
  coord_fixed()


##alluvial
library(ggalluvial)

create_alluvial_plot <- function(all_results) {
  plot_data <- data.frame(CellSubtype = character(),
                          BackboneGene = character(),
                          GOFunction = character(),
                          freq = numeric())
  
  for (result_list in all_results) {
    for (group_name in names(result_list)) {
      if (length(result_list[[group_name]]$Present_Backbone_Drivers) > 0) {
        cell_subtype <- gsub("^CLUSTER_WALKTRAP_", "", group_name)
        cell_subtype <- sub("_Group[0-9]+.*$", "", cell_subtype)
        
        for (backbone_gene in result_list[[group_name]]$Present_Backbone_Drivers) {
          for (go_function in result_list[[group_name]]$GO_Results$term_id) {
            plot_data <- rbind(plot_data, data.frame(CellSubtype = cell_subtype,
                                                     BackboneGene = backbone_gene,
                                                     GOFunction = go_function,
                                                     freq = 1))
          }
        }
      }
    }
  }
  plot_data$freq <- ave(plot_data$freq, plot_data$CellSubtype, plot_data$BackboneGene, plot_data$GOFunction, FUN = length)
  return(plot_data)
}

# Example usage
alluvial_plot_data <- create_alluvial_plot(all_results_with_backbone_presence)
alluvial_plot_data$CellSubtype <- gsub("_graph_obj", "", alluvial_plot_data$CellSubtype)

# Custom color palette (modify colors as needed)
my_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")


svg("/corgi/sebas/tcell_multi/alluvial/alluvial.svg", width = 11, height = 8.5)
ggplot(alluvial_plot_data, aes(axis1 = CellSubtype, axis2 = BackboneGene, axis3 = GOFunction, y = freq)) +
  geom_alluvium(aes(fill = BackboneGene)) +
  geom_stratum(color = "black", size= 1) + # Add gray borders to the stratum
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  ggtitle("Alluvial Plot: Cell Subtype, Backbone Gene, and GO Function") +
  scale_fill_manual(values = my_colors) # Apply custom color palette

dev.off()
my_colors <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", 
               "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")






#####_--------------------------$$$$$$$$$$$$$$$$    REWIRING METRIC

### local edge probability/ density for FDNs  * SS

### SS
extract_and_compare_genes_v2 <- function(all_results) {
  group_genes <- list()
  group_drivers <- list()
  
  # Extract genes and backbone drivers from each group
  for (result_list in all_results) {
    for (group_name in names(result_list)) {
      simple_name <- gsub("^CLUSTER_WALKTRAP_", "", group_name)
      simple_name <- sub("_Group[0-9]+.*$", "", simple_name)
      simple_name <- gsub("_graph_obj", "", simple_name)
      
      if (!is.null(group_genes[[simple_name]])) {
        group_genes[[simple_name]] <- unique(c(group_genes[[simple_name]], result_list[[group_name]]$Genes))
      } else {
        group_genes[[simple_name]] <- result_list[[group_name]]$Genes
      }
      
      if (!is.null(group_drivers[[simple_name]])) {
        group_drivers[[simple_name]] <- unique(c(group_drivers[[simple_name]], result_list[[group_name]]$Present_Backbone_Drivers))
      } else {
        group_drivers[[simple_name]] <- result_list[[group_name]]$Present_Backbone_Drivers
      }
    }
  }
  
  group_names <- names(group_genes)
  
  # Initialize dataframe
  final_df <- data.frame(Subtype = character(), 
                         Comparison = character(), 
                         SS_Value = numeric(), 
                         Backbone = character(), 
                         row.names = NULL)
  
  # Compute overlap and create dataframe
  for (i in group_names) {
    for (j in group_names) {
      if (i != j) {  # Exclude self vs. self comparisons
        ss_value <- szymkiewicz_simpson(group_genes[[i]], group_genes[[j]])
        
        if (!is.null(group_drivers[[i]]) && length(group_drivers[[i]]) > 0) {
          for (backbone in group_drivers[[i]]) {
            final_df <- rbind(final_df, data.frame(Subtype = i, Comparison = j, SS_Value = ss_value, Backbone = backbone))
          }
        } else {
          final_df <- rbind(final_df, data.frame(Subtype = i, Comparison = j, SS_Value = ss_value, Backbone = NA))
        }
      }
    }
  }
  
  return(final_df)
} #excluding self vs self interactions

overlaps_melteddf<- extract_and_compare_genes_v2(all_results_with_backbone_presence)


### local edge density
library(igraph)

library(igraph)

# Function to calculate local edge density in a directed graph
local_edge_density_directed <- function(graph, node) {
  # Get the neighbors of the node (both in and out)
  neighbors <- unique(c(neighbors(graph, node, mode = "out"),
                        neighbors(graph, node, mode = "in")))
  
  # Create a subgraph with these neighbors
  subgraph <- induced_subgraph(graph, neighbors)
  
  # Calculate the number of actual and possible edges
  actual_edges <- gsize(subgraph)
  total_nodes <- vcount(subgraph)
  possible_edges <- total_nodes * (total_nodes - 1)
  
  # Calculate and return the edge density
  edge_density <- actual_edges / possible_edges
  return(edge_density)
}


# Function to calculate local edge density for multiple graphs and nodes, and store in a matrix
compute_local_edge_density_matrix <- function(graphs, nodes) {
  # Prepare matrix with appropriate dimensions and names
  graph_names <- sapply(names(graphs), function(name) sub("_graph_obj$", "", name))
  density_matrix <- matrix(NA, nrow = length(graphs), ncol = length(nodes), 
                           dimnames = list(graph_names, nodes))
  
  # Iterate over each graph
  for (i in seq_along(graphs)) {
    graph <- graphs[[i]]
    
    # Apply the function to each node in the node list
    for (node in nodes) {
      if(node %in% V(graph)$name) { # Check if the node exists in the graph
        density_matrix[i, node] <- local_edge_density_directed(graph, node)
      }
    }
  }
  
  return(density_matrix)
}


local_edge_dens_result_matrix <- compute_local_edge_density_matrix(graphs, backbone_bois)

overlaps_melteddf

# Assuming you have overlaps_melteddf and local_edge_density_matrix ready
# Add a new column for Rewiring_score
overlaps_melteddf$Rewiring_score <- NA

for (i in 1:nrow(overlaps_melteddf)) {
  subtype <- overlaps_melteddf$Subtype[i]
  backbone_gene <- overlaps_melteddf$Backbone[i]
  ss_value <- overlaps_melteddf$SS_Value[i]
  
  # Find the corresponding local edge density
  if (subtype %in% rownames(local_edge_dens_result_matrix) && backbone_gene %in% colnames(local_edge_dens_result_matrix)) {
    local_density <- local_edge_dens_result_matrix[subtype, backbone_gene]
    # Calculate Rewiring_score
    overlaps_melteddf$Rewiring_score[i] <- ss_value * local_density
  }
}

# View the updated dataframe
print(overlaps_melteddf)

#minmax normalization
overlaps_melteddf$Norm_Rewiring_score <- (overlaps_melteddf$Rewiring_score - min(overlaps_melteddf$Rewiring_score)) / (max(overlaps_melteddf$Rewiring_score) - min(overlaps_melteddf$Rewiring_score))




## getting pseudotime per celltype

library(SeuratWrappers)
library(nlme)
library(sf)
library(monocle3)
library(Seurat)
genes_NKT <- readRDS("/corgi/martin/multiome_T/RDS_final/T_cells_noCC_DICE_ATAC_final_peaks to genes_NKT.RDS")
DefaultAssay(genes_NKT)<-'RNA'
# Extract data from Seurat object
integrated <- genes_NKT

integrated <- ScaleData(integrated)
integrated <- FindVariableFeatures(integrated)
integrated <- RunPCA(integrated)
integrated <- RunUMAP(integrated, dims = 1:30, reduction.name = "UMAP")
integrated <- FindNeighbors(integrated, dims = 1:30)
integrated <- FindClusters(integrated)
Idents(integrated)=integrated@meta.data$SingleR.labels
cds <- as.cell_data_set(integrated)
cds <- cluster_cells(cds)
cds<- learn_graph(cds)
Tcell.cds=cds
#cluster info
list.cluster=genes_NKT@active.ident
Tcell.cds@clusters$UMAP$clusters=list.cluster
#umap coords from seurat
Tcell.cds@int_colData@listData$reducedDims$UMAP= genes_NKT@reductions$umap_RNA@cell.embeddings


#plot
plot_cells(Tcell.cds,reduction_method = 'UMAP',color_cells_by = 'ident',label_groups_by_cluster = F,
           group_label_size = 5)

Tcell.cds=learn_graph(Tcell.cds,use_partition = F)


plot_cells(Tcell.cds,
           color_cells_by = 'ident',
           label_groups_by_cluster = F,
           label_branch_points = F,
           label_roots = F,
           label_leaves = T,
           group_label_size = 5,
          )

cellInfo <- data.frame(seuratCluster=Idents(genes_NKT))



#order
Tcell.cds=order_cells(Tcell.cds,reduction_method = 'UMAP', root_cells = c(rownames(cellInfo)[cellInfo$seuratCluster=='T cells, CD4+, naive']))
plot_cells(Tcell.cds,
           color_cells_by = 'pseudotime',
           label_groups_by_cluster = F,
           label_branch_points = F,
           label_roots = F,
           label_leaves = F)



#pseudotime values
Tcell.cds$monocle3_pseudotime=pseudotime(Tcell.cds)
pseudotime.RNA=data.frame(colData(Tcell.cds))
library(ggplot2)
library(ggridges)
ggplot(pseudotime.RNA,aes(monocle3_pseudotime,reorder(ident,monocle3_pseudotime,median),fill= ident)) + geom_boxplot()
ggplot(pseudotime.RNA, aes(x = reorder(ident, monocle3_pseudotime, median), y = monocle3_pseudotime, fill = ident)) +
  geom_violin()


pseudotime.cellsRNA=pseudotime.RNA
pseudotime.cellsRNA$barcode=rownames(pseudotime.cellsRNA)
pseudotime.cellsRNA=pseudotime.cellsRNA[,c('monocle3_pseudotime','ident')]
library(dplyr)

avg_pseudo_df <- pseudotime.cellsRNA %>%
  group_by(ident) %>%
  summarize(avg_pseudo = mean(monocle3_pseudotime))

# Display the resulting dataframe
print(avg_pseudo_df)

remove(pseudotime.cellsRNA)
remove(integrated)
remove(cds)

# Remove the prefix from the ident column with a more robust regular expression
avg_pseudo_df$ident <- gsub("^T cells, CD4\\+\\,\\s*", "", avg_pseudo_df$ident)
avg_pseudo_df <- avg_pseudo_df %>%
  mutate(ident = case_when(
    ident == "naive TREG" ~ "Treg_naive",
    ident == "NKT cells" ~ "NK",
    ident == "naive, stimulated" ~ "naive_stim",
    ident == "memory TREG" ~ "Treg_mem",
    TRUE ~ ident # Default case to keep other values as they are
  ))

avg_pseudo_df <- avg_pseudo_df %>% 
  rename(Subtype = ident)

# Display the updated dataframe
print(avg_pseudo_df)


# Display the updated dataframe
print(avg_pseudo_df)


# Merge the dataframes
pseudo_rewscore_df <- merge(avg_pseudo_df, overlaps_melteddf, by = "Subtype")

pseudo_rewscore_df <- pseudo_rewscore_df[order(pseudo_rewscore_df$avg_pseudo), ]
pseudo_rewscore_df$Subtype <- factor(pseudo_rewscore_df$Subtype, levels = unique(pseudo_rewscore_df$Subtype[order(pseudo_rewscore_df$avg_pseudo)]))

ggplot(merged_df, aes(x = Subtype, y = Norm_Rewiring_score, group = Backbone, color = Backbone)) +
  stat_smooth(method = "loess", se = FALSE, aes(group = Backbone), size = 1.2) +  # Add smooth lines
  theme_minimal() +
  labs(title = "Line Plot of Normalized Rewiring Score by Backbone",
       x = "Subtype (Ordered by Avg Pseudotime)",
       y = "Normalized Rewiring Score",
       color = "Backbone") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) # Rotate x-axis labels for better readability

ggplot(merged_df, aes(x = Subtype, y = Norm_Rewiring_score, group = Backbone, color = Backbone)) +
  stat_smooth(method = "loess", se = FALSE, aes(group = Backbone), size = 1.2) +  # Only smooth lines
  facet_wrap(~ Backbone, nrow = 4, ncol = 2) +  # Facet by Backbone
  theme_minimal() +
  labs(title = "Line Plot of Normalized Rewiring Score by Backbone",
       x = "Subtype (Ordered by Avg Pseudotime)",
       y = "Normalized Rewiring Score",
       color = "Backbone") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))





####custom ordering
# Define the custom order for subtypes
custom_order <- c("naive", "Treg_naive", "naive_stim", "Th2", "Th17", "Treg_mem", "Th1_17", "Th1", "NK")

# Set the levels of Subtype in pseudo_rewscore_df according to the custom order
pseudo_rewscore_df$Subtype <- factor(pseudo_rewscore_df$Subtype, levels = custom_order)

# Define the backbones and the subtype you want to check
backbones_missing_naive <- c("STAT4", "MAF", "CREM", "RBPJ", "ZNF292")
subtype_to_check <- "naive"

# Function to add a row with zeros if a specific backbone-subtype combination is missing
add_missing_zeros <- function(df, backbone, subtype) {
  if (!any(df$Backbone == backbone & df$Subtype == subtype)) {
    # Create a row with zeros
    zero_row <- data.frame(
      Subtype = factor(subtype, levels = levels(df$Subtype)),
      avg_pseudo = NA,  # Assuming avg_pseudo is not relevant
      Comparison = NA,  # Assuming Comparison is not relevant
      SS_Value = 0,
      Backbone = backbone,
      Rewiring_score = 0,
      Norm_Rewiring_score = 0,
      Inverse = NA  # Assuming Inverse is not relevant
    )
    # Add the row to the dataframe
    df <- rbind(df, zero_row)
  }
  return(df)
}

# Apply the function to each backbone for the "naive" subtype
for (backbone in backbones_missing_naive) {
  pseudo_rewscore_df <- add_missing_zeros(pseudo_rewscore_df, backbone, subtype_to_check)
}

ggplot(pseudo_rewscore_df, aes(x = Subtype, y = Rewiring_score, group = Backbone, color = Backbone)) +
  stat_smooth(method = "loess", se = FALSE, aes(group = Backbone), size = 1.2) +  # Only smooth lines
  facet_wrap(~ Backbone, nrow = 4, ncol = 2) +  # Facet by Backbone
  theme_minimal() +
  labs(title = "Line Plot of Normalized Rewiring Score by Backbone",
       x = "Subtype (Ordered by Avg Pseudotime)",
       y = "Normalized Rewiring Score",
       color = "Backbone") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))












### plot using the average
avg_rewiring_score <- pseudo_rewscore_df %>%
  group_by(Subtype, Backbone) %>%
  summarize(Avg_Norm_Rewiring_score = mean(Norm_Rewiring_score, na.rm = TRUE)) %>%
  ungroup()

# Reorder the Subtype factor based on avg_pseudo
avg_rewiring_score$Subtype <- factor(avg_rewiring_score$Subtype, levels = unique(avg_rewiring_score$Subtype[order(pseudo_rewscore_df$avg_pseudo[pseudo_rewscore_df$Subtype %in% avg_rewiring_score$Subtype])]))

ggplot(avg_rewiring_score, aes(x = Subtype, y = Avg_Norm_Rewiring_score, group = Backbone, color = Backbone)) +
  geom_line(size = 1.2) +  # Increase line thickness
  geom_point() + 
  geom_smooth()+
  stat_smooth(method = "loess", se = FALSE, aes(group = Backbone), size = 1.2) +  # Add smooth lines
  theme_minimal() +
  labs(title = "Line Plot of Average Normalized Rewiring Score by Backbone",
       x = "Subtype (Ordered by Avg Pseudotime)",
       y = "Average Normalized Rewiring Score",
       color = "Backbone") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) # Rotate x-axis labels for better readability






##### function to extract dataframe with both GO:BP and GO:MF # same function as before, just making it look for go_bp
perform_BP_enrichment_with_backbone_presence <- function(community_obj_name, backbone_genes, organism = "hsapiens") {
  # If backbone_genes is a list of object names, fetch and combine genes from all driver objects
  if (is.list(backbone_genes)) {
    backbone_genes <- unique(unlist(lapply(backbone_genes, function(obj_name) get(obj_name))))
  }
  
  # Fetch the community object
  community_obj <- get(community_obj_name)
  membership <- membership(community_obj)
  
  results_list <- list()
  
  for (group_id in unique(membership)) {
    group_members <- names(which(membership == group_id))
    
    # Identify backbone/driver genes present in the group
    present_backbone_genes <- intersect(group_members, backbone_genes)
    
    # Check if the community contains any backbone genes
    if (length(present_backbone_genes) > 0) {
      if (length(group_members) > 10) {
        gostres <- try(gost(query = group_members,
                            exclude_iea = TRUE,
                            user_threshold = 0.05,
                            correction_method = "fdr",
                            domain_scope = 'annotated',
                            sources = c('GO:BP'),
                            organism = organism), silent = TRUE)
        
        if (inherits(gostres, "try-error")) {
          next
        } else {
          gprofout <- gostres$result[c("p_value", 'term_id', "term_name")]
          
          if (!is.null(gprofout) && nrow(gprofout) > 0) {
            gprofout <- gprofout[gprofout$p_value < 0.05,]
            if (nrow(gprofout) > 0) {
              gprofout <- gprofout[1, ]
              results_list[[paste0(community_obj_name, "_Group", group_id)]] <- list(
                GO_Results = gprofout,
                Genes = group_members,
                Present_Backbone_Drivers = present_backbone_genes
              )
            }
          }
        }
      }
    }
  }
  
  return(results_list)
}
all_results_with_backbone_presence_MF<-lapply(cluster_list, perform_GO_enrichment_with_backbone_presence, backbone_bois)
all_results_with_backbone_presence_BP<-lapply(cluster_list, perform_BP_enrichment_with_backbone_presence, backbone_bois)


extract_backbone_GO_info <- function(results_list) {
  all_data <- do.call(rbind, lapply(results_list, function(group_list) {
    do.call(rbind, lapply(names(group_list), function(group_name) {
      if (!is.null(group_list[[group_name]]$Present_Backbone_Drivers)) {
        backbone_genes <- group_list[[group_name]]$Present_Backbone_Drivers
        go_results <- group_list[[group_name]]$GO_Results
        gene_count <- length(group_list[[group_name]]$Genes)
        subtype <- gsub("CLUSTER_WALKTRAP_(.*)_graph_obj.*", "\\1", group_name)
        
        data.frame(
          GO_pval = go_results$p_value,
          GO_term_id = go_results$term_id,
          term_name = go_results$term_name,
          size_group = gene_count,
          Subtype = subtype,
          Present_Backbone_Driver = backbone_genes,
          stringsAsFactors = FALSE
        )
      } else {
        NULL
      }
    }))
  }))
  all_data
}

# Apply the function to your nested list
backbone_GOMF_info_df <- extract_backbone_GO_info(all_results_with_backbone_presence_MF)

# Assuming your dataframe is named backbone_GO_info_df
colnames(backbone_GOMF_info_df)[colnames(backbone_GOMF_info_df) %in% c("GO_pval", "GO_term_id", "term_name", "size_group")] <- 
  paste0(colnames(backbone_GOMF_info_df)[colnames(backbone_GOMF_info_df) %in% c("GO_pval", "GO_term_id", "term_name", "size_group")], "_MF")

# View the updated dataframe
str(backbone_GO_info_df)
str(backbone_GOMF_info_df)

write.csv(backbone_GO_info_df, '//corgi/sebas/tcell_multi/go_results_community/backbone_BP_info_df.csv')
write.csv(backbone_GOMF_info_df, '//corgi/sebas/tcell_multi/go_results_community/backbone_MF_info_df.csv')


##### Remaking the rewiring metric: substitute local edge density with the degree centrality

actual_edges_to_degree_centrality <- function(graph, node) {
  # Get the in-degree and out-degree of the node
  in_degree <- degree(graph, node, mode = "in")
  out_degree <- degree(graph, node, mode = "out")
  total_degree <- in_degree + out_degree
  
  # Get the neighbors of the node (both in and out)
  neighbors <- unique(c(neighbors(graph, node, mode = "out"),
                        neighbors(graph, node, mode = "in")))
  
  # Create a subgraph with these neighbors
  subgraph <- induced_subgraph(graph, neighbors)
  
  # Calculate the number of actual edges in the subgraph
  actual_edges <- gsize(subgraph)
  
  # Calculate and return the ratio of actual edges to total degree
  if (total_degree > 0) {
    ratio <-  total_degree
  } else {
    ratio <- 0  # Avoid division by zero
  }
  
  return(ratio)
}

# Function to calculate the ratio of actual edges to degree centrality for multiple graphs and nodes, and store in a matrix
compute_edges_to_degree_centrality_matrix <- function(graphs, nodes) {
  # Prepare matrix with appropriate dimensions and names
  graph_names <- sapply(names(graphs), function(name) sub("_graph_obj$", "", name))
  ratio_matrix <- matrix(NA, nrow = length(graphs), ncol = length(nodes),
                         dimnames = list(graph_names, nodes))
  
  # Iterate over each graph
  for (i in seq_along(graphs)) {
    graph <- graphs[[i]]
    
    # Apply the function to each node in the node list
    for (node in nodes) {
      if(node %in% V(graph)$name) { # Check if the node exists in the graph
        ratio_matrix[i, node] <- actual_edges_to_degree_centrality(graph, node)
      }
    }
  }
  
  return(ratio_matrix)
}

edges_to_degree_centrality_matrix <- compute_edges_to_degree_centrality_matrix(graphs, backbone_bois)

overlaps_melteddf

# Add a new column for Rewiring_score
overlaps_melteddf$Rewiring_score <- NA

for (i in 1:nrow(overlaps_melteddf)) {
  subtype <- overlaps_melteddf$Subtype[i]
  backbone_gene <- overlaps_melteddf$Backbone[i]
  ss_value <- overlaps_melteddf$SS_Value[i]
  
  # Find the corresponding actual edges to degree centrality ratio
  if (subtype %in% rownames(edges_to_degree_centrality_matrix) && backbone_gene %in% colnames(edges_to_degree_centrality_matrix)) {
    edge_to_centrality_ratio <- edges_to_degree_centrality_matrix[subtype, backbone_gene]
    # Calculate Rewiring_score using the ratio
    overlaps_melteddf$Rewiring_score[i] <- ss_value * edge_to_centrality_ratio
  }
}

# View the updated dataframe
print(overlaps_melteddf)

#minmax normalization
overlaps_melteddf$Norm_Rewiring_score <- (overlaps_melteddf$Rewiring_score - min(overlaps_melteddf$Rewiring_score)) / (max(overlaps_melteddf$Rewiring_score) - min(overlaps_melteddf$Rewiring_score))




#####plotting
# Merge the dataframes
pseudo_rewscore_df <- merge(avg_pseudo_df, overlaps_melteddf, by = "Subtype")

pseudo_rewscore_df <- pseudo_rewscore_df[order(pseudo_rewscore_df$avg_pseudo), ]
pseudo_rewscore_df$Subtype <- factor(pseudo_rewscore_df$Subtype, levels = unique(pseudo_rewscore_df$Subtype[order(pseudo_rewscore_df$avg_pseudo)]))

ggplot(pseudo_rewscore_df, aes(x = Subtype, y = Norm_Rewiring_score, group = Backbone, color = Backbone)) +
  stat_smooth(method = "loess", se = FALSE, aes(group = Backbone), size = 1.2) +  # Add smooth lines
  theme_minimal() +
  labs(title = "Line Plot of Normalized Rewiring Score by Backbone",
       x = "Subtype (Ordered by Avg Pseudotime)",
       y = "Normalized Rewiring Score",
       color = "Backbone") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) # Rotate x-axis labels for better readability

ggplot(pseudo_rewscore_df, aes(x = Subtype, y = Norm_Rewiring_score, group = Backbone, color = Backbone)) +
  stat_smooth(method = "loess", se = FALSE, aes(group = Backbone), size = 1.2) +  # Only smooth lines
  facet_wrap(~ Backbone, nrow = 4, ncol = 2) +  # Facet by Backbone
  theme_minimal() +
  labs(title = "Line Plot of Normalized Rewiring Score by Backbone",
       x = "Subtype (Ordered by Avg Pseudotime)",
       y = "Normalized Rewiring Score",
       color = "Backbone") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))





















#### part2: overlapping for specific genes with specific communities
paired_elements_modularity1=c()
paired_elements_modularity2=c()

for(i in modularity_groups_ls){
  stp=gsub("_graph_obj_wlktrp_out.*", "",i)
  subtype=stp #subtype
  
  j=get(i)
  for(grp in j){
  if(length(grp)>10 & any(backbone.bois %in% grp)){
    backbone.presence=backbone.bois[which(backbone.bois %in% grp)]
    if(length(backbone.presence)>1){
      for(boi in backbone.presence){
        name=paste0(stp,'_',boi)
        paired_elements_modularity1=c(paired_elements_modularity1,rep(subtype))
        paired_elements_modularity2=c(paired_elements_modularity2,boi)
        assign(name,grp)
      }
    }else{
    print(backbone.presence)
    #print(grp)
    name=paste0(stp,'_',backbone.presence)
    paired_elements_modularity1=c(paired_elements_modularity1,rep(subtype))
    paired_elements_modularity2=c(paired_elements_modularity2,backbone.presence)
    assign(name,grp)
    }
    }
  }
}
paired_elements_modularity=data.frame(cbind(paired_elements_modularity1,paired_elements_modularity2))
paired_elements_modularity

install.packages("bios2mds")
install.packages('Bios2cor')
library(bios2mds)
library(Bios2cor)

# Custom overlap coefficient function for each specific set
overlap_coefficient <- function(set1, set2) {
  if (length(set1) == 0 || length(set2) == 0) {
    return(0)
  }
  intersection_size <- length(intersect(set1, set2))
  min_set_size <- min(length(set1), length(set2))
  return(intersection_size / min_set_size)
}




overlap_backbone=matrix(0, nrow = 0, ncol = 0)




#####$$$$$
for (i in 1:nrow(paired_elements_modularity)) {
  ecs <- paired_elements_modularity$paired_elements_modularity1[i]
  ecsprime <- paired_elements_modularity$paired_elements_modularity2[i]
  ecsterm <- paste0(ecs, '_', ecsprime)
  
  for (j in 1:nrow(paired_elements_modularity)) {
    if (i != j) {
      wai <- paired_elements_modularity$paired_elements_modularity1[j]
      waiprime <- paired_elements_modularity$paired_elements_modularity2[j]
      waiterm <- paste0(wai, '_', waiprime)
      
      # Assuming you have data for each element in your workspace
      set1 <- get(ecsterm)
      set2 <- get(waiterm)
      
      # Expand the result matrix if necessary
      if (!ecsterm %in% rownames(overlap_backbone)) {
        overlap_backbone <- rbind(overlap_backbone, rep(0, ncol(overlap_backbone)))
        rownames(overlap_backbone)[nrow(overlap_backbone)] <- ecsterm
      }
      if (!waiterm %in% colnames(overlap_backbone)) {
        overlap_backbone <- cbind(overlap_backbone, rep(0, nrow(overlap_backbone)))
        colnames(overlap_backbone)[ncol(overlap_backbone)] <- waiterm
      }
      
      overlap_backbone[rownames(overlap_backbone) == ecsterm, colnames(overlap_backbone) == waiterm] <- overlap_coefficient(set1, set2)
    }
  }
}

library(ribiosUtils)
overlap_coefficient(naive_BCL6,naive_stim_BCL6)
library(pheatmap)

overlap_backbone #export as supplemental csv
overlap_backbone_df=data.frame(overlap_backbone)
##Figure for overlap analysis
# Calculate the range of values in the matrix
value_range <- range(overlap_backbone)

# Create a data frame for the histogram values
hist_values <- hist(overlap_backbone, breaks = "FD", plot = FALSE)
hist_df <- data.frame(values = hist_values$mids, counts = hist_values$counts)

# Set the color palette from ggsci
npg_palette <- pal_npg()(3)

# Create the barplot using ggplot
ggplot(hist_df, aes(x = values, y = counts, fill = values)) +
  geom_bar(stat = "identity") +
  scale_fill_gradientn(colors = npg_palette[c(1:3)]) +
  labs(x = "Overlap coeff", y = "Frequency") +
  theme_minimal() +
  theme(
    axis.text = element_text(face = "bold", size = 20),
    axis.title = element_text(face = "bold", size = 20),
    panel.grid = element_blank(),
    legend.position = "none"
  )

molten_overlap_backbone= melt(overlap_backbone)
### subset matrix
####values above 0.9 but not 1 that have overlap
high_overlap_backbone=data.frame(
       Var2 = unique(molten_overlap_backbone$Var2[molten_overlap_backbone$value > 0.9 & molten_overlap_backbone$value != 1]),
       Var1 = unique(molten_overlap_backbone$Var1[molten_overlap_backbone$value > 0.9 & molten_overlap_backbone$value != 1])
   )

high_overlap_backbone <- data.frame(
  Var2 = unique(molten_overlap_backbone$Var2[molten_overlap_backbone$value > 0.9 & molten_overlap_backbone$value != 1]),
  Var1 = unique(molten_overlap_backbone$Var1[molten_overlap_backbone$value > 0.9 & molten_overlap_backbone$value != 1])
)

# Extract the last term after the last "_"
high_overlap_backbone$LastTerm_Var2 <- sapply(strsplit(as.character(high_overlap_backbone$Var2), "_"), function(x) tail(x, n = 1))
high_overlap_backbone$LastTerm_Var1 <- sapply(strsplit(as.character(high_overlap_backbone$Var1), "_"), function(x) tail(x, n = 1))

# Remove the last term from the original strings
high_overlap_backbone$Var2 <- sapply(strsplit(as.character(high_overlap_backbone$Var2), "_"), function(x) paste(head(x, -1), collapse = "_"))
high_overlap_backbone$Var1 <- sapply(strsplit(as.character(high_overlap_backbone$Var1), "_"), function(x) paste(head(x, -1), collapse = "_"))

#### taking from backbone_functions only the high overlapping values, for alluvial plotting, FIGURE 3
high_overlap_backbone_functions= backbone_functions[backbone_functions$subtype == 'naive' | backbone_functions$subtype == 'Treg_naive', ]
high_overlap_backbone_functions=high_overlap_backbone_functions[!high_overlap_backbone_functions$gene %in%  c('HMGB2','ZNF292'),]
top=ggplot(data = high_overlap_backbone_functions,
           aes(axis1 = subtype, axis2 = termGO, axis3 =gene ,
               y = loggrsize)) +
  scale_x_discrete(limits = c("Subtype", "Function", "Gene"), expand = c(0.2, 0.02)) +
  geom_alluvium(aes(fill = gene),width=1/3, linewidth=1.5) +
  scale_fill_npg()+
  geom_stratum(aes(fill = gene),width=1/3, linewidth = 1.5)+
  coord_flip()+
  theme_void() +
  theme(legend.position = "none")
top+ggrepel::geom_text_repel(
  aes(label = ifelse(after_stat(x) == 1, as.character(after_stat(stratum)), NA), fontface='bold'),
  stat = "stratum", size = 6, direction = "y", nudge_x = +.3
) +  ggrepel::geom_text_repel(
  aes(label = ifelse(after_stat(x) == 2, as.character(after_stat(stratum)), NA), fontface='bold'),
  stat = "stratum", size = 6, direction = "y", nudge_x = +.4
) +  ggrepel::geom_text_repel(
  aes(label = ifelse(after_stat(x) == 3, as.character(after_stat(stratum)), NA), fontface='bold'),
  stat = "stratum", size = 6, direction = "y", nudge_x = +.3
)










low_overlap_backbone <- data.frame(
  Var2 = unique(molten_overlap_backbone$Var2[molten_overlap_backbone$value < 0.1 & molten_overlap_backbone$value != 0]),
  Var1 = unique(molten_overlap_backbone$Var1[molten_overlap_backbone$value < 0.1 & molten_overlap_backbone$value != 0])
)

# Extract the last term after the last "_"
low_overlap_backbone$LastTerm_Var2 <- sapply(strsplit(as.character(low_overlap_backbone$Var2), "_"), function(x) tail(x, n = 1))
low_overlap_backbone$LastTerm_Var1 <- sapply(strsplit(as.character(low_overlap_backbone$Var1), "_"), function(x) tail(x, n = 1))

# Remove the last term from the original strings
low_overlap_backbone$Var2 <- sapply(strsplit(as.character(low_overlap_backbone$Var2), "_"), function(x) paste(head(x, -1), collapse = "_"))
low_overlap_backbone$Var1 <- sapply(strsplit(as.character(low_overlap_backbone$Var1), "_"), function(x) paste(head(x, -1), collapse = "_"))

low_overlap_backbone_functions= backbone_functions[!rownames(backbone_functions) %in% rownames(high_overlap_backbone_functions),]
low=ggplot(data = low_overlap_backbone_functions,
           aes(axis1 = subtype, axis2 = termGO, axis3 =gene ,
               y = loggrsize)) +
  scale_x_discrete(limits = c("Subtype", "Function", "Gene"), expand = c(0.2, 0.02)) +
  geom_alluvium(aes(fill = gene),width=1/3, linewidth=1.5) +
  scale_fill_npg()+
  geom_stratum(aes(fill = gene),width=1/3, linewidth = 1.5)+
  coord_flip()+
  theme_void() +
  theme(legend.position = "none")
low+ggrepel::geom_text_repel(
  aes(label = ifelse(after_stat(x) == 1, as.character(after_stat(stratum)), NA), fontface='bold'),
  stat = "stratum", size = 6, direction = "y", nudge_x = +.3
) +  ggrepel::geom_text_repel(
  aes(label = ifelse(after_stat(x) == 2, as.character(after_stat(stratum)), NA), fontface='bold'),
  stat = "stratum", size = 6, direction = "y", nudge_x = +.4
) +  ggrepel::geom_text_repel(
  aes(label = ifelse(after_stat(x) == 3, as.character(after_stat(stratum)), NA), fontface='bold'),
  stat = "stratum", size = 6, direction = "y", nudge_x = +.3
)













### exporting all communities for all subtypes and all backbones
counter <- 1

for (i in modularity_groups_ls) {
  stp <- gsub("_graph_obj_wlktrp_out.*", "", i)
  subtype <- stp
  
  j <- get(i)
  for (grp in j) {
    if (length(grp) > 10 & any(backbone.bois %in% grp)) {
      to.write <- grp
      name1 <- subtype
      nameout <- paste0(name1, '_', counter, '.csv')
      write.csv(to.write, file = paste0('/corgi/sebas/tcell_multi/multiome_local_seb_20221207/backbone_communities/', nameout))
      
      counter <- counter + 1
    }
  }
}



#heatmap for overlap, all vs all
pheatmap::pheatmap(overlap_backbone,main = 'All vs ALL community overlapping') #looks crappy
#function to plot whichever cell type against all
create_heatmap <- function(term) {
  library(pheatmap)
  
  # Filter rows containing the specified term and exclude rows not containing the term
  subset_df <- overlap_backbone_df[grepl(term, rownames(overlap_backbone_df)), ]
  
  # Create heatmap using pheatmap
  pheatmap(subset_df, clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean")
}
create_heatmap('Th1_17')









#### Overlapping coefficient between communities with driver bearers
paired_elements_modularity1=c()
paired_elements_modularity2=c()

for(i in modularity_groups_ls){
  stp=gsub("_graph_obj_wlktrp_out.*", "",i)
  subtype=stp #subtype
  
  j=get(i)
  for(grp in j){
    if(length(grp)>10 & any(driver.bois %in% grp)){
      driver.presence=driver.bois[which(driver.bois %in% grp)]
      if(length(driver.presence)>1){
        for(boi in driver.presence){
          name=paste0(stp,'_',boi)
          paired_elements_modularity1=c(paired_elements_modularity1,rep(subtype))
          paired_elements_modularity2=c(paired_elements_modularity2,boi)
          assign(name,grp)
        }
      }else{
        print(driver.presence)
        #print(grp)
        name=paste0(stp,'_',driver.presence)
        paired_elements_modularity1=c(paired_elements_modularity1,rep(subtype))
        paired_elements_modularity2=c(paired_elements_modularity2,driver.presence)
        assign(name,grp)
      }
    }
  }
}

paired_elements_modularity_driver=data.frame(cbind(paired_elements_modularity1,paired_elements_modularity2))
paired_elements_modularity_driver



# Custom overlap coefficient function
overlap_coefficient <- function(set1, set2) {
  if (length(set1) == 0 || length(set2) == 0) {
    return(0)
  }
  intersection_size <- length(intersect(set1, set2))
  min_set_size <- min(length(set1), length(set2))
  return(intersection_size / min_set_size)
}


overlap_driver=matrix(0, nrow = 0, ncol = 0)



#####$$$$$
for (i in 1:nrow(paired_elements_modularity_driver)) {
  ecs <- paired_elements_modularity_driver$paired_elements_modularity1[i]
  ecsprime <- paired_elements_modularity_driver$paired_elements_modularity2[i]
  ecsterm <- paste0(ecs, '_', ecsprime)
  
  for (j in 1:nrow(paired_elements_modularity_driver)) {
    if (i != j) {
      wai <- paired_elements_modularity_driver$paired_elements_modularity1[j]
      waiprime <- paired_elements_modularity_driver$paired_elements_modularity2[j]
      waiterm <- paste0(wai, '_', waiprime)
      
      # Assuming you have data for each element in your workspace
      set1 <- get(ecsterm)
      set2 <- get(waiterm)
      
      # Expand the result matrix if necessary
      if (!ecsterm %in% rownames(overlap_driver)) {
        overlap_driver <- rbind(overlap_driver, rep(0, ncol(overlap_driver)))
        rownames(overlap_driver)[nrow(overlap_driver)] <- ecsterm
      }
      if (!waiterm %in% colnames(overlap_driver)) {
        overlap_driver <- cbind(overlap_driver, rep(0, nrow(overlap_driver)))
        colnames(overlap_driver)[ncol(overlap_driver)] <- waiterm
      }
      
      overlap_driver[rownames(overlap_driver) == ecsterm, colnames(overlap_driver) == waiterm] <- overlap_coefficient(set1, set2)
    }
  }
}

overlap_driver_df=data.frame(overlap_driver)
pheatmap::pheatmap(overlap_driver, cluster_rows = F, cluster_cols = F)
overlap_driver_df_m=melt(overlap_driver_df)

#naive sitm
# Filter rows and columns containing 'naive'
rows <- grep('naive', rownames(overlap_driver), invert = T)

# Get columns NOT containing 'naive'
cols <- grep('naive', colnames(overlap_driver), invert = TRUE)

# Create a sub-matrix with filtered rows and columns
filtered_overlap_driver <- overlap_driver[rows, cols]

# Print the filtered matrix

pheatmap::pheatmap(filtered_overlap_driver)

hist(filtered_overlap_driver)




############- ORIGIN OF THE REWIRING

library(tidygraph)
library(tibble)

edge_list <-
  Th1_graph_obj %>%
  activate(edges) %>%
  data.frame()

edge_list$corr

headnode_list<-
  Th17_graph_obj %>%
  activate(nodes) %>%
  data.frame()
node_list

library(ggplot2)

ggplot(edge_list[edge_list$from_node=='AHR',], aes(x = regions, fill = from_node)) +
  geom_bar() +
  geom_text(stat = 'count', aes(label = after_stat(count)), vjust = -0.5, color = 'black') +
  labs(x = "Regions", y = "Frequency", fill = "From Node")





#CART cell relevant genes, their possition in the network. Targets or drivers?>
cart_genes<-read.table('/corgi/sebas/NCBI_data/genes_ens_ids.tsv', header = F) #68 long
backbone.bois %in% cart_genes$V1 #no backbone bois in the list
driver.bois %in% cart_genes$V1 #no driver bois in the list

edge_list$from_node[which(cart_genes$V1 %in% edge_list$to_node)]

#list of _graphs
#objects
barp_size_cart_targets=NULL
for(i in names(objects)){
  stp= i #subtype name
  i=get(i)
  edge_list <-
    i %>%
    activate(edges) %>%
    data.frame()
  
  cartgenes_tf=edge_list$from_node[which(cart_genes$V1 %in% edge_list$to_node)] #which TFS target cart relevant genes present in the targets
  cartgenes_targets=edge_list$to_node[which(edge_list$to_node %in% cart_genes$V1)]
  percent_of_target=(length(unique(cartgenes_targets))/length(unique(cart_genes$V1)))*100
  
  tablegenes=data.frame(sort(table(edge_list$from_node), decreasing = T))
  barp_size_cart_targets=rbind(barp_size_cart_targets, data.frame(subtype=stp, percent_of_target=sprintf("%.2f", percent_of_target), TF=unique(cartgenes_tf)))
}

# Use wesanderson color palette
my_palette <- wes_palette("Darjeeling2", 4)

# Plot bar chart
ggplot(barp_size_cart_targets, aes(x=subtype, y=as.numeric(percent_of_target), fill=TF)) +
  geom_bar(stat="identity", position="dodge") +
  scale_fill_manual(values=my_palette) +
  labs(x="Subtype", y="Percent of Targets", fill="TF") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1, face = 'bold', size=10), legend.position="bottom",)

library(ggplot2)
library(dplyr)
library(tidygraph)
library(gridExtra)

#making plots of TF and how many targets it has, per cell type
# define the folder path where you want to save the plot
folder_path <- "/corgi/sebas/tcell_multi/figures"

# create a list to store the plots
plots <- list()

for(i in names(objects)){
  stp= i #subtype name
  i=get(i)
  edge_list <-
    i %>%
    activate(edges) %>%
    data.frame()
  cartgenes_targets=edge_list$from_node[which(cart_genes$V1 %in% edge_list$to_node)] #which TFS target cart genes present in the targets
  tablegenes=data.frame(sort(table(edge_list$from_node), decreasing = T))
  tablegenes$Var1 <- factor(tablegenes$Var1, levels = tablegenes$Var1[order(tablegenes$Freq, decreasing = TRUE)])
  print(stp)
  print(tablegenes)
  p <- ggplot(data = tablegenes, aes(x = Var1, y = Freq)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    labs(y = "# of Targets", x = "Gene") + 
    theme(axis.text.x = element_text(face = "bold"))+
    coord_flip() +
    ggtitle(stp) # add the title to the plot
  plots[[stp]] <- p # add the plot to the list
}

# combine the plots into a multi-page pdf
pdf(file = paste0(folder_path, "/plots.pdf"), width = 10, height = 14)
for (i in names(objects)) {
  print(plots[[i]])
}
dev.off()


##########---------------ENHANCERS
library(tidygraph)
edge_list <-
  Th1_graph_obj %>%
  activate(edges) %>%
  data.frame()

extract_regions <- function(edge_list) {
  
  result_df <- data.frame(chromosome = character(),
                          start = integer(),
                          end = integer(),
                          bb_from = character(),
                          bb_to = character(),
                          stringsAsFactors = FALSE)
  
  for (i in 1:nrow(edge_list)) {
    
    regions <- strsplit(edge_list[i, "regions"], ";")[[1]]
    bb_from <- edge_list[i, "from_node"]  # Extract from_node value
    bb_to <- edge_list[i, "to_node"]  # Extract to_node value
    
    for (j in 1:length(regions)) {
      
      region_info <- strsplit(regions[j], "-")
      chromosome <- region_info[[1]][1]
      start <- as.integer(region_info[[1]][2])
      end <- as.integer(region_info[[1]][3])
      
      result_df <- rbind(result_df, data.frame(chromosome = chromosome,
                                               start = start,
                                               end = end,
                                               bb_from = bb_from,
                                               bb_to = bb_to,
                                               stringsAsFactors = FALSE))
      
    }
    
  }
  
  result_df <- unique(result_df)
  return(result_df)
  
}


extract_regions(edge_list) #works

for(i in names(objects)){
  stp=gsub("_graph_obj_wlktrp_out.*", "",i)
  subtype=stp #subtype
  j=get(i) #access i as object
  
  edge_list <-
    j %>%
    activate(edges) %>%
    data.frame()
  out=extract_regions(edge_list)
  nameoftheobject=paste0(stp,'_regions')
  assign(nameoftheobject,out)
}

regions <- mget(grep("_regions", ls(), value = TRUE, fixed = TRUE), envir = .GlobalEnv)
regions <- Filter(function(x) !is.function(x) && !is.numeric(x), regions)




for (i in names(regions)) {
  print(i)
  i = regions[[i]]
  duplicated_start = duplicated(i$start)
  duplicated_end = duplicated(i$end)
  
  if (any(duplicated_start) || any(duplicated_end)) {
    # Find rows that appear more than once
    duplicated_rows = c(which(duplicated_start), which(duplicated_end))
    duplicated_rows = sort(unique(duplicated_rows))
    # Print the duplicated rows
    print(i[duplicated_rows, ])
  } else {
    # Print "Unique all over"
    print("Unique all over")
  }
}

### export
for (i in names(regions)) {
  thename=i
  i = regions[[i]]
  print(i)
  write.csv(i,paste0('/corgi/sebas/tcell_multi/hotspot/regions_with_bb_genes/',thename,'.csv'),col.names = T)
}







cart_genes$V1
----------------------
#plot generic network for workflow description
# Create a graph with 5 nodes and edges between them
graph_edges <- data.frame(from = c('a', 'a','a', 'b','b', 'c','c'),
                          to   = c('b', 'c', 'g','c','f', 'd','e'))

g <- graph_from_data_frame(graph_edges, directed = FALSE)

# Define the layout
layout <- layout_with_kk(g)

# Define the node colors
V(g)$color <- ifelse(V(g)$name %in% c('a', 'b', 'c'), "lightgreen", "cyan")

# Plot the graph with bold vertex labels
plot(g, layout = layout, vertex.label.color = "black", vertex.label.font = 2, vertex.size = 20)


--------------

######- USEFUL FUNCTIONS
#rename crap
objects_to_rename <- grep("CLUSTER_WALKTRAP_", ls(), value = TRUE) # finds objects that contain "_graph"
for (i in objects_to_rename){
  new_name <- gsub('_obj$','',i) # remove those objects from the environment
  assign(new_name, get(i), envir = .GlobalEnv)
  rm(i)
}

#removing crap
objects_to_remove <- grep("^CLUSTER_WALKTRAP.*_obj$", ls(), value = TRUE) # finds objects that contain "_graph"
if(length(objects_to_remove)>0){
  rm(list = objects_to_remove) # remove those objects from the environment
}
remove(objects_to_remove)
















#### repeating XGBOOST and GLM but with a GREAT approach.
#GLM
pmbc_1103 <- infer_grn(
  pmbc_1103,
  peak_to_gene_method = 'GREAT',
  method = 'glm'
)

pmbc_1103 <- find_modules(pmbc_1103)
great_modules <- NetworkModules(pmbc_1103)
great_modules@meta
table(great_modules@meta$tf)
great_glm_modules=names(table(great_modules@meta$tf))
plot_gof(pmbc_1103, point_size=3)

plot_module_metrics(pmbc_1103) #size of modules

pmbc_1103 <- get_network_graph(pmbc_1103)
plot_network_graph(pmbc_1103) 
write.csv(great_glm_modules,"great_glm_pando_sebas.csv")

canon_bois[which(canon_bois %in% great_glm_modules)]
master_regs[which(master_regs %in% great_glm_modules)]

#XGBOOST
pmbc_1103 <- infer_grn(
  pmbc_1103,
  peak_to_gene_method = 'GREAT',
  method = 'xgb'
)

pmbc_1103 <- find_modules(pmbc_1103)
great_modules_xgb <- NetworkModules(pmbc_1103)
great_modules_xgb@meta$tf
great_grm_modules=names(table(great_modules_xgb@meta$tf))
table(great_modules_xgb@meta$tf)

write.csv(great_grm_modules,"great_xgboost_pando_sebas.csv")

canon_bois[which(canon_bois %in% great_grm_modules)]
master_regs[which(master_regs %in% great_grm_modules)]
###plotting
plot_gof(pmbc_1103, point_size=3)

plot_module_metrics(pmbc_1103) #size of modules

pmbc_1103 <- get_network_graph(pmbc_1103)
plot_network_graph(pmbc_1103) #grn














####clustering based on regulon activity
nPcs <- c(10) # For toy dataset
# nPcs <- c(5,15,50)
scenicOptions@settings$seed <- 123 # same seed for all of them
# Run t-SNE with different settings:
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50))
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50), onlyHighConf=TRUE, filePrefix="int/tSNE_oHC")
# Plot as pdf (individual files in int/):
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_", list.files("int"), value=T), value=T))
par(mfrow=c(length(nPcs), 3))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, varName="seuratCluster", cex=.5)



#checking regulon activity per ct

cellInfo <- data.frame(seuratCluster=Idents(pmbc_1103))

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$seuratCluster),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))


ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity")

topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
viewTable(topRegulators)



###cell type specific regulators based on Regulon specificty score
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "seuratCluster"])



rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)
plotRSS_oneSet(rss, setName = "T cells, CD4+, Th2")







#same but binarized with the min prc of cells with that regulon
minPerc <- .1
binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")
cellInfo_binarizedCells <- cellInfo[which(rownames(cellInfo)%in% colnames(binaryRegulonActivity)),, drop=FALSE]
regulonActivity_byCellType_Binarized <- sapply(split(rownames(cellInfo_binarizedCells), cellInfo_binarizedCells$seuratCluster), 
                                               function(cells) rowMeans(binaryRegulonActivity[,cells, drop=FALSE]))
binaryActPerc_subset <- regulonActivity_byCellType_Binarized[which(rowSums(regulonActivity_byCellType_Binarized>minPerc)>0),]
ComplexHeatmap::Heatmap(binaryActPerc_subset, name="Regulon activity (%)", col = c("white","pink","red"))



#regulon activities on UMAP from seurat

dr_coords <- Embeddings(pmbc_1103, reduction="umap")

tfs <- c("FOXP1",'GATA3','NFKB1')
par(mfrow=c(2,2))
AUCell::AUCell_plotTSNE(dr_coords, cellsAUC=selectRegulons(regulonAUC, tfs), plots = "AUC")

######SCENIC+ ON PYTHON. EXPORTING PMBC FILE TO ANNDATA TO IMPORT INTO PYTHON...
install.packages('renv')
renv::init()
library(reticulate)
renv::use_python()
py_pkgs <- c(
  "scanpy",
  "python-igraph",
  "louvain"
)
reticulate::py_install(py_pkgs)
pkgs <- c(
  "renv",
  "reticulate",
  "png",
  "ggplot2",
  "BiocManager",
  "Seurat"
)
reticulate::py_install(pkgs)



DefaultAssay(pmbc_1103)='RNA'
exprs <- GetAssayData(seurat)
meta <- seurat[[]]
feature_meta <- GetAssay(seurat)[[]]
embedding <- Embeddings(seurat, "umap")

#######----
#######----
#######----
#######----
#######----END BLOCK 4
#######----
#######----
#######----


######### BlOCK 5:TF family abundance analysis
#data from http://humantfs.ccbr.utoronto.ca/download.php
#

TFpaper= read.csv('family info analysis/TF_names_v_1.01.txt',header = F)
TFfampaper= read.csv('family info analysis/TableS2.csv', header = T)
TFfampaper=TFfampaper[,c(1,3,4,5,6)]
TFfampaper$TFname=sub('\\_.*','',TFfampaper$Motif.ID)

#function to create list of tf based on the marker expression genes for each cat,
#then create df of TF family enrichment per ct
#then create barplot
#input used: curated gene marker expression list per ct from seurat

familyabundance_ct= function(x){
  x$genes=rownames(x)
  TFmarkersofX=intersect(x$genes,TFpaper$V1)
  TFfamofX=intersect(TFmarkersofX,TFfampaper$TFname)
  df_fam=TFfampaper[TFfampaper$TFname %in% TFfamofX,]
  df_fam=data.frame(table(df_fam$TF.Family))
  rownames(df_fam)=df_fam$Var1
  df_fam=df_fam[2]
  #colnames(df_fam)=c('Family','Freq')
  return(df_fam)
}

markerdict= noquote(c(ls(pattern = '_markers')))

tf_fam_ct=list()
for( i in markerdict){
  tf_fam_ct[[i]]=familyabundance_ct(get(noquote(i)))
}
tf_fam_ct= do.call(rbind,tf_fam_ct)
library(stringr)
tf_fam_ct$fam=str_split_fixed(rownames(tf_fam_ct), "\\.", 2)[,2]
tf_fam_ct$subtype=str_split_fixed(rownames(tf_fam_ct), "\\.", 2)[,1]
tf_fam_ct$subtype=gsub('_markers','',tf_fam_ct$subtype)


#checking frequency sums
for(i in unique(tf_fam_ct$subtype)){
  print(i)
  print(sum(tf_fam_ct$Freq[tf_fam_ct$subtype==i]))
}
for( i in markerdict){
  print(familyabundance_ct(get(noquote(i))))
}

###plotting
library(hrbrthemes)
library(viridis)
# Small multiple
ggplot(tf_fam_ct, aes(fill=fam, y=Freq, x=subtype, )) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_viridis(discrete = T) +
  ggtitle("Transcription factor family abundance per subtype") +
  theme_ipsum() +
  xlab("")

######Maybe do the same with oepn regions?



###########checking for biases of Regulon targets to TF
regulons_bias= regulontargets[regulontargets$TF %in% master_regulons,]

tf_fam_regulon=list()
for(i in unique(regulons_bias$TF)){
  present=intersect(regulons_bias$gene[regulons_bias$TF== i],TFfampaper$TFname)
  #presentpct= (length(intersect(regulons_bias$gene[regulons_bias$TF== i],TFfampaper$TFname)) 
  #/ length(regulons_bias$gene[regulons_bias$TF==i])) * 100
  #missingpct= 100-presentpct
  df_fam=TFfampaper[TFfampaper$TFname %in% present,]
  #table(df_fam$TF.Family)
  #print(presentpct)
  #print(missingpct)
  print(data.frame(table(df_fam$TF.Family)))
  tf_fam_regulon[[i]]=data.frame(table(df_fam$TF.Family))
}

tf_fam_regulon= do.call(rbind,tf_fam_regulon)
library(stringr)
tf_fam_regulon$TF=str_split_fixed(rownames(tf_fam_regulon), "\\.", 2)[,1]
colnames(tf_fam_regulon)[1]='fam'

ggplot(tf_fam_regulon, aes(fill=fam, y=Freq, x=TF, )) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_viridis(discrete = T) +
  ggtitle("Transcription factor family abundance per regulon") +
  theme_ipsum() +
  xlab("")

#enrichment of regulons
setEnrichrSite("Enrichr")
websiteLive <- TRUE
dbs <- listEnrichrDbs()
if (is.null(dbs)) websiteLive <- FALSE
if (websiteLive) head(dbs)



dbs <- c("KEGG_2019_Human", "GO_Cellular_Component_2021", 
         'KEGG_2021_Human','Reactome_2016','GO_Molecular_Function_2021',
         'GO_Biological_Process_2021')
if (websiteLive) {enriched <- enrichr(regulontargets$gene[regulontargets$TF=='NFKB1'], dbs)
}

#taking the go biological processes 2021 ds
if (websiteLive) plotEnrich(enriched[[6]], 
                            showTerms = 5,numChar = 60, 
                            y = "Count", orderBy = "P.value",
                            title = 'x')



####-------------------------------------------------

#getting all expressed genes #Johans idea
allexpressed=gene.rna
allexpressed=data.frame(rowMeans(gene.rna))
allexpressed$gene=rownames(allexpressed)
allexpressed=allexpressed[allexpressed$rowMeans.gene.rna. > 0,]
colnames(allexpressed)[1]='Row mean expression'

hist(log(allexpressed$`Row mean expression`),breaks = 250)
abline(v = median(log(allexpressed)),                       # Add line for mean
       col = "red",
       lwd = 2)
text(x = median(log(allexpressed)) * 1.7,                   # Add text for mean
     y = median(log(allexpressed)) * 1.7,
     paste("Median =", median(log(allexpressed))),
     col = "red",
     cex = 1,
     )

#merging expressed genes and tf family info
cluster_mar_gene
genedictionary= merge(cluster_mar_gene,allexpressed,by='gene',all=TRUE)
#Merging TF names
TF_frompaper=TFpaper
colnames(TF_frompaper)='TF'
TF_frompaper$gene=TF_frompaper$TF
genedictionary= merge(genedictionary,TF_frompaper,by='gene',all=TRUE)
remove(TF_frompaper)
#merging with correlation matrix bw motif and gene (jaspdict)
jaspdict2=jaspdict
colnames(jaspdict2)=c('gene','motif_jaspar','pearson_corr_diffgenexp_vs_predictedmotif')
genedictionary= merge(genedictionary,jaspdict2,by='gene',all=TRUE)
remove(jaspdict2)
#no clear cutoff. taking everyhint thats more than 0

hist(allexpressed$`Row mean expression`,breaks = 40)

####-------------------------------------------------

#######----
#######----
#######----
#######----
#######----END BLOCK 5
#######----
#######----
#######----


######### BlOCK 6: 
#RNA VELOCITY USING VELOCYTO
#cite!
#Doing it in python. Package does not work in R
cellinfotoexport=cellInfo
cellinfotoexport$barcode=rownames(cellinfotoexport)
write.csv(cellinfotoexport,'Documents/multimodal Tcell/cell_annotation_seurat_tcell.csv')
#Potentially useless: The assumption of calculating the ratio of spliced/unspliced as being a good way of estimating terminal states is incomplete. 
#It does not take into account any biological context. 

#markov montecarlo chain using CellRank on python for pseudotime analysis.
#cellrank is shit

#monocle for pseudotime
devtools::install_github('satijalab/seurat-data')
install.packages('terra', repos='https://rspatial.r-universe.dev')
devtools::install_github('cole-trapnell-lab/monocle3', force = T)
remotes::install_github('satijalab/seurat-wrappers')
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils',
                       'HDF5Array', 'terra', 'ggrastr'))
library(SeuratWrappers)
library(nlme)
library(sf)
library(monocle3)
library(Seurat)
genes_NKT <- readRDS("/corgi/martin/multiome_T/RDS_final/T_cells_noCC_DICE_ATAC_final_peaks to genes_NKT.RDS")
DefaultAssay(genes_NKT)<-'RNA'
# Extract data from Seurat object
integrated <- genes_NKT

integrated <- ScaleData(integrated)
integrated <- FindVariableFeatures(integrated)
integrated <- RunPCA(integrated)
integrated <- RunUMAP(integrated, dims = 1:30, reduction.name = "UMAP")
integrated <- FindNeighbors(integrated, dims = 1:30)
integrated <- FindClusters(integrated)
Idents(integrated)=integrated@meta.data$SingleR.labels
cds <- as.cell_data_set(integrated)
cds <- cluster_cells(cds)
cds<- learn_graph(cds)
Tcell.cds=cds
#cluster info
list.cluster=genes_NKT@active.ident
Tcell.cds@clusters$UMAP$clusters=list.cluster
#umap coords from seurat
Tcell.cds@int_colData@listData$reducedDims$UMAP= genes_NKT@reductions$umap_RNA@cell.embeddings


#plot
plot_cells(Tcell.cds,reduction_method = 'UMAP',color_cells_by = 'ident',label_groups_by_cluster = F,
           group_label_size = 5)

Tcell.cds=learn_graph(Tcell.cds,use_partition = F)


plot_cells(Tcell.cds,
           color_cells_by = 'ident',
           label_groups_by_cluster = F,
           label_branch_points = F,
           label_roots = F,
           label_leaves = F,
           group_label_size = 5)

cellInfo <- data.frame(seuratCluster=Idents(genes_NKT))



#order
Tcell.cds=order_cells(Tcell.cds,reduction_method = 'UMAP', root_cells = c(rownames(cellInfo)[cellInfo$seuratCluster=='T cells, CD4+, naive']))
plot_cells(Tcell.cds,
           color_cells_by = 'pseudotime',
           label_groups_by_cluster = F,
           label_branch_points = F,
           label_roots = F,
           label_leaves = F)
#pseudotime values
Tcell.cds$monocle3_pseudotime=pseudotime(Tcell.cds)
pseudotime.RNA=data.frame(colData(Tcell.cds))
library(ggplot2)
library(ggridges)
ggplot(pseudotime.RNA,aes(monocle3_pseudotime,reorder(ident,monocle3_pseudotime,median),fill= ident)) + geom_boxplot()
ggplot(pseudotime.RNA, aes(x = reorder(ident, monocle3_pseudotime, median), y = monocle3_pseudotime, fill = ident)) +
  geom_violin()


pseudotime.cellsRNA=pseudotime.RNA
pseudotime.cellsRNA$barcode=rownames(pseudotime.cellsRNA)
pseudotime.cellsRNA=pseudotime.cellsRNA[,c('monocle3_pseudotime','ident')]
library(dplyr)

avg_pseudo_df <- pseudotime.cellsRNA %>%
  group_by(ident) %>%
  summarize(avg_pseudo = mean(monocle3_pseudotime))

# Display the resulting dataframe
print(avg_pseudo_df)










###genes on pseudotime
pseudotime.genes.RNA=graph_test(Tcell.cds,neighbor_graph = 'principal_graph',cores = 6)
#filtering
pseudotime.genes.RNA<-pseudotime.genes.RNA %>%
  arrange(q_value) %>%
  filter(status=='OK') %>%
  filter(q_value<0.05)
pseudo.rna.gene.plot=pseudotime.genes.RNA
pseudo.rna.gene.plot$gene=rownames(pseudo.rna.gene.plot)
pseudo.rna.gene.plot=pseudo.rna.gene.plot[,c('gene','q_value','p_value')]

FeaturePlot(pmbc_1103,features = master_regulons)

#CHECK pseudotime in seurat and do heatmap of genes~cells

pmbc_1103$pseudotime= pseudotime(Tcell.cds)

gene.zsc=as.data.frame(GetAssayData(pmbc_1103,slot = 'scale.data'))
gene.zsc= gene.zsc[,intersect(colnames(gene.zsc),pseudotime.cellsRNA$barcode)]
gene.zsc= gene.zsc[intersect(rownames(gene.zsc), pseudo.rna.gene.plot$gene),]
gene.zsc$gene=rownames(gene.zsc)
gene.zsc=melt(gene.zsc)
colnames(gene.zsc)[2:3]=c('barcode','zscore')
gene.zsc=merge(gene.zsc,pseudotime.cellsRNA)

ggplot(gene.zsc[gene.zsc$gene==intersect(gene.zsc$gene, master_regulons),],aes(reorder(barcode,monocle3_pseudotime), gene, fill=zscore))+geom_tile()

FeaturePlot(pmbc_1103,features = 'pseudotime',label = T,repel = T)+ ggtitle('RNA pseudotime')
#better on 3d

plot.data <- FetchData(object = pmbc_1103, vars = c("x3dumap_1", "x3dumap_2", "x3dumap_3",'pseudotime'))
plot.data$label <- paste(rownames(plot.data))

plot_ly(data = plot.data, 
               x = ~x3dumap_1, y = ~x3dumap_2, z = ~x3dumap_3, 
               color = ~pseudotime, 
               mode = "markers", 
               marker = list(size = 3, width=2), # controls size of points
               text=~label, #This is that extra column we made earlier for which we will use for cell ID
               hoverinfo="text") #When you visualize your plotly object, hovering your mouse pointer over a point shows cell names


#############REPEATING THE SAME BUT WITH ATAC- note does not make much sense
DefaultAssay(pmbc_1103)='ATAC'
Tcell.cds <- as.cell_data_set(pmbc_1103)
#### getting clustering info from seurat
#partitions
rec.part= c(rep(1,length(Tcell.cds@colData@rownames)))
names(rec.part)=Tcell.cds@colData@rownames
rec.part=as.factor(rec.part)
Tcell.cds@clusters$UMAP$partitions=rec.part
#cluster info
list.cluster=pmbc_1103@active.ident
Tcell.cds@clusters$UMAP$clusters=list.cluster
#umap coords from seurat
Tcell.cds@int_colData@listData$reducedDims$UMAP= pmbc_1103@reductions$UMAP_ATAC@cell.embeddings
#plot
plot_cells(Tcell.cds,reduction_method = 'UMAP',color_cells_by = 'ident',label_groups_by_cluster = F,
           group_label_size = 5) +theme(legend.position = 'right')

Tcell.cds=learn_graph(Tcell.cds,use_partition = F,)

plot_cells(Tcell.cds,
           color_cells_by = 'ident',
           label_groups_by_cluster = F,
           label_branch_points = F,
           label_roots = F,
           label_leaves = F,
           group_label_size = 5)
#order
Tcell.cds=order_cells(Tcell.cds,reduction_method = 'UMAP', root_cells = rownames(cellInfo)[cellInfo$seuratCluster=='T cells, CD4+, naive'])
plot_cells(Tcell.cds,
           color_cells_by = 'pseudotime',
           label_groups_by_cluster = F,
           label_branch_points = F,
           label_roots = F,
           label_leaves = F)
#pseudotime values
Tcell.cds$monocle3_pseudotime=pseudotime(Tcell.cds)
pseudotime.ATAC=data.frame(colData(Tcell.cds))
ggplot(pseudotime.ATAC,aes(monocle3_pseudotime,reorder(ident,monocle3_pseudotime,median),fill= ident)) + geom_boxplot()
###genes on pseudotime
pseudotime.genes.ATAC=graph_test(Tcell.cds,neighbor_graph = 'principal_graph',cores = 6)

pseudotime.genes.ATAC %>%
  arrange(q_value) %>%
  filter(status=='OK')%>%
  head()

FeaturePlot(pmbc_1103,features = master_regulons)

#CHECK pseudotime in seurat

pmbc_1103$pseudotime_atac= pseudotime(Tcell.cds)

FeaturePlot(pmbc_1103,features = 'pseudotime_atac',label = T,repel = T)+ ggtitle('ATAC pseudotime')
#better on 3d

plot.data <- FetchData(object = pmbc_1103, vars = c("x3d_umap_atac_1", "x3d_umap_atac_2", "x3d_umap_atac_3",'pseudotime_atac'))
plot.data$label <- paste(rownames(plot.data))

plot_ly(data = plot.data, 
        x = ~x3d_umap_atac_1, y = ~x3d_umap_atac_2, z = ~x3d_umap_atac_3, 
        color = ~pseudotime_atac, 
        mode = "markers", 
        marker = list(size = 3, width=2), # controls size of points
        text=~label, #This is that extra column we made earlier for which we will use for cell ID
        hoverinfo="text") #When you visualize your plotly object, hovering your mouse pointer over a point shows cell names














# stopping with monocle. Dynverse seems more useful






#### Dynverse : NOTE. Running it on beagle. does not run well here
devtools::install_github("dynverse/dyno")
library(dyno)
###preparing the data
cellInfo2 <- data.frame(seuratCluster=Idents(pmbc_1103))
cellInfo2$cell_id=rownames(cellInfo2)
colnames(cellInfo2)[1]='group_id'

cellInfo3 <- as.character(cellInfo2$group_id)
names(cellInfo3) <- cellInfo2$cell_id



gene.rna.dynverse=gene.rna[rownames(gene.rna) %in% cluster_mar_gene$gene,]
#gene.rna.dynverse=melt(gene.rna.dynverse)
#colnames(gene.rna.dynverse)[1]='barcode'
#gene.rna.dynverse_dict=merge(gene.rna.dynverse, cellInfo2, by='barcode')
#gene.rna.dynverse= data.frame(gene.rna.dynverse_dict[,2:3])
#gene.rna.dynverse=pivot_wider(gene.rna.dynverse_dict,names_from = 'seuratCluster',
#                              values_from = 'value')


gene.rna.dynverse=wrap_expression(counts = pmbc_1103@assays$RNA@counts[rownames(pmbc_1103@assays$RNA@counts) %in% rownames(gene.rna.dynverse),],
                                  expression = pmbc_1103@assays$RNA@scale.data[rownames(pmbc_1103@assays$RNA@scale.data) %in% rownames(gene.rna.dynverse),])

mycounts=pmbc_1103@assays$RNA@counts[rownames(pmbc_1103@assays$RNA@counts) %in% rownames(gene.rna.dynverse),]
myexpression=pmbc_1103@assays$RNA@scale.data[rownames(pmbc_1103@assays$RNA@scale.data) %in% rownames(gene.rna.dynverse),]
gene.rna.dynverse <- 
  wrap_data(
    id = "pmbc_1103",
    cell_ids = colnames(pmbc_1103@assays$RNA@counts[rownames(pmbc_1103@assays$RNA@counts) %in% rownames(gene.rna.dynverse),])
  ) %>% 
  add_expression(
    counts = mycounts,
    expression = myexpression
  )




dynguidelines::guidelines_shiny(gene.rna.dynverse)

gene.rna.dynverse=add_prior_information(dataset = gene.rna.dynverse,
                                        groups_id = cellInfo3)

# Reproduces the guidelines as created in the shiny app
# Reproduces the guidelines as created in the shiny app
answers <- dynguidelines::answer_questions(
  multiple_disconnected = FALSE, 
  expect_topology = TRUE, 
  expected_topology = "multifurcation", 
  n_cells = 2537, 
  n_features = 6036, 
  memory = "3GB", 
  prior_information = c("start_id", "groups_id", "dimred"), 
  docker = TRUE
)
guidelines <- dynguidelines::guidelines(answers = answers)

#running the trajectory inferrence
model=infer_trajectory(gene.rna.dynverse,"raceid_stemid", verbose = T)

#observing
plot_dimred(model) 
plot_dendro(model,'pseudotime')
plot_graph(model)
plot_onedim(add_root(model))
plot_heatmap(model, expression_source = pmbc_1103@assays$RNA@scale.data)



#######----
#######----
#######----
#######----
#######---- 6
#######----
#######----
#######----

######## BlOCK 7: 
#TF footprinting
#https://stuartlab.org/signac/articles/footprint.html

DefaultAssay(pmbc_1103) <- "peaks"

pwm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information
pcmbc_1103 <- AddMotifs(pmbc_1103, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pwm)
#run signac's footprinting
pmbc_1103_footprint <-Footprint(
  object = pmbc_1103,
  motif.name = c("ATF4","CEBPG",#"DDIT3",
                 "ELF1",
                 "GATA3","JUNB","NFKB1","REL",
                 "RUNX3","SREBF2","YY1"),
  genome = BSgenome.Hsapiens.UCSC.hg38,
  in.peaks = T,
  upstream = 500,
  downstream=500
)

PlotFootprint(pmbc_1103_footprint, features = c("ATF4","CEBPG",#"DDIT3",
                                       "ELF1",
                                       "GATA3","JUNB","NFKB1","REL",
                                       "RUNX3","SREBF2","YY1"))


PlotFootprint(pmbc_1103_footprint, features = c("JUNB"), label = T, 
              label.idents = c('T cells, CD4+, naive','T cells, CD4+, naive TREG',
                               'T cells, CD4+, naive, stimulated'))
PlotFootprint(pmbc_1103_footprint, features = c("ATF4"), label = T, 
              label.idents = c('T cells, CD4+, naive','T cells, CD4+, naive TREG',
                               'T cells, CD4+, naive, stimulated'))



#TF footprinting for canonical TFs


pmbc_1103_footprint_CANON<-Footprint(
  object = pmbc_1103,
  motif.name = c('BCL6',#'CD3D',
                 #'ICOS'
                 #'CXCR5', #TFH
                 'STAT3','IRF4',#'IL17A'
                 'RUNX1','BATF',
                 'RORA',#'APOE'
                 #'PIK3C2B'
                 'RORC',#TH17
                 'TBX21','STAT1',
                 #'STAT4',
                 #'ANXA3',#Th1
                 'GATA3',#'SMAD2',
                 #'STAT6',
                 'RUNX2',#Th2
                 'FOXP3',
                 #'IKZF2',
                 #'CCL5','CTLA4',#Treg
                 #'TRGV9',
                 #'KLRD1',
                 #'GNLY',
                 #'NKG7',
                 'EOMES'), #NKT 
  genome = BSgenome.Hsapiens.UCSC.hg38,
  in.peaks = T,
  upstream = 500,
  downstream=500
)

PlotFootprint(pmbc_1103_footprint_CANON, features = c("BATF"), label = T, 
              label.idents = c('T cells, CD4+, naive','T cells, CD4+, naive TREG',
                               'T cells, CD4+, naive, stimulated','T cells, CD4+, Th17'))
PlotFootprint(pmbc_1103_footprint_CANON, features = c("RORC"), label = T, 
              label.idents = c('T cells, CD4+, naive','T cells, CD4+, naive TREG',
                               'T cells, CD4+, naive, stimulated','T cells, CD4+, Th17')) ###TH17

PlotFootprint(pmbc_1103_footprint_CANON, features = c("EOMES"), label = T, 
              label.idents = c('T cells, CD4+, naive','T cells, CD4+, naive TREG',
                               'T cells, CD4+, naive, stimulated','NK cells'))#NKT

topie=data.frame(table(Idents(pmbc_1103_footprint_CANON)))
topie$pct=round(topie$Freq/sum(topie$Freq)*100)
topie$varpct=paste(topie$Var1,topie$pct)
topie$varpct=paste(topie$varpct, '%', sep = '')
pie(topie$Freq,labels = topie$varpct)

#######----
#######----
#######----
#######----
#######---- END BLOCK 7
#######----
#######----
#######----







####### RUNNING TOP20 GENES ON ENRICHR AND SEEING WASSUP

enriched <- enrichr((top20), dbs)
if (websiteLive) plotEnrich(enriched[[2]], showTerms = 10, numChar = 50, y = "Count", orderBy = "P.value")









#######----Linking peaks to genes
DefaultAssay(pmbc_1103) <- "peaks"

# first compute the GC content for each peak
pmbc_1103<- RegionStats(pmbc_1103, genome = BSgenome.Hsapiens.UCSC.hg38)

# link peaks to genes #TAKES A WHILE TO COMPUTE
pmbc_1103 <- LinkPeaks(
  object = pmbc_1103,
  method = 'spearman',
  peak.assay = "peaks",
  expression.assay = "RNA",
  min.cells = 3, pvalue_cutoff = 0.05,
  verbose = T
)
#saving

####coverage plots
DefaultAssay(pmbc_1103) <- "peaks"
CoveragePlot(
  object = pmbc_1103,
  region = 'TCF7',
  features = 'PDE3B',
  expression.assay = 'RNA',
  extend.upstream = 500,
  extend.downstream = 500
)

DefaultAssay(pbmc_chromvared) <- "RNA"
VlnPlot(pbmc_chromvared, features = 'CTCF',pt.size = 0)
FeaturePlot(
  object = pbmc_chromvared,
  features = 'TOP2A'
  ,pt.size = 1
)



#####---- DA PEAKS between cell types--- THIS CAN BE USED FOR ANNOTATION
DefaultAssay(pmbc_1103) <- 'peaks'

da_peaks_naive <- FindMarkers(
  object = pmbc_1103,assay = 'peaks', 
)

head(da_peaks_test)

plot1 <- VlnPlot(
  object = pmbc_1103,
  features = rownames(da_peaks_test)[1],
  pt.size = 0.1
)
plot2 <- FeaturePlot(
  object = pmbc_1103,
  features = rownames(da_peaks_test)[1],
  pt.size = 3
)

plot1 | plot2





### Gene activity: checking chromatin accessibility linked to each gene.

#note, on geneactivity, there is possibility to extend upstream or downstream of TSS
gene.activities <- GeneActivity(pmbc_1103,verbose = T)
# add the gene activity matrix to the Seurat object as a new assay and normalize it
pmbc_1103[['Gene.activity']] <- CreateAssayObject(counts = gene.activities)
pmbc_1103 <- NormalizeData(
  object = pmbc_1103,
  assay = 'Gene.activity',
  normalization.method = 'LogNormalize',
  scale.factor = median(pmbc_1103$nCount_RNA)
)

geneacts=data.frame(as.matrix(GetAssayData(pmbc_1103,assay = 'Gene.activity')))



#plotting SOME marker genes on atac: NOTE, SHIT NOT WORKN
DefaultAssay(pmbc_1103) <- 'Gene.activity'
FeaturePlot(
  object = pmbc_1103, pt.size = .4,
  features = 'FOXP3'
)
#####Plotting genes of g2m on atac umap


FeaturePlot(
  object = pbmc,
  features = 'HSF1',
  pt.size = 2,
  min.cutoff = 'q10',
  max.cutoff = 'q90')


#vackup
pbmc2=pbmc



##Another way of finding DA regions is checking FOLD CHANGE
###fc <- FoldChange(pbmc, ident.1 = "CD4 Naive", ident.2 = "CD14 Mono")
##head(fc)


####- check closest gene to each peak
open_Treg <- rownames(da_peaks_test[da_peaks_test$avg_log2FC > 0.3, ])
open_Th1 <- rownames(da_peaks_test[da_peaks_test$avg_log2FC < -0.3, ])

closest_genes_treg <- ClosestFeature(pbmc, regions = open_Treg)
closest_genes_th1 <- ClosestFeature(pbmc, regions = open_Th1)
head(closest_genes_treg)
head(closest_genes_th1)

CoveragePlot(
  object = pbmc,
  region = rownames(da_peaks_test)[1],
  extend.upstream = 40000,
  extend.downstream = 20000
)






####---ChromVAR

#Get motifs from JASPAR
# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

pfm_test=getMatrixSet(
  x= JASPAR2020,
  opts= list(collect = 'FAM', tax_group= 'vertebrates', all_versions = F)
)


#add motif info
pbmc<- AddMotifs(pbmc, 
                 genome = BSgenome.Hsapiens.UCSC.hg38,pfm = pfm)


#####-make function out of this goddamit

da_peaks_naive <- FindMarkers(
  object = pbmc,
  ident.1 = 'T cells, CD4+, naive',
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)
# get top differentially accessible peaks
top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005, ])
# find peaks open in  cells.
#Chosing set of background peaks.
open.peaks <- AccessiblePeaks(pbmc,idents = 'gDT')
meta.feature <- GetAssayData(pbmc, assay = "peaks", slot = "meta.features")
peaks.matched <- MatchRegionStats(
  meta.feature = meta.feature[open.peaks, ],
  query.feature = meta.feature[top.da.peak, ],
  n = 50000
)


# test enriched motifs
enriched.motifs <- FindMotifs(
  object = pbmc,
  background = peaks.matched,
  features = top.da.peak
)
MotifPlot(
  object = pbmc,
  motifs = head(rownames(enriched.motifs))
)

#COMPUTE MOTIF ACTIVITY
pbmc_chromvared<- RunChromVAR(
  object = pmbc_1103,
  genome = BSgenome.Hsapiens.UCSC.hg38
)

DefaultAssay(pbmc_chromvared) <- 'chromvar'


# look at the activity of MOTIF X
FeaturePlot(
  object = pbmc_chromvared,
  features = 'MA0486.1',
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.6)


###Differential activity scores bw CTs. USES LIMMA
##supposedly, similar results as enrichment test 

differential.activity <- FindMarkers(
  object = pbmc_chromvared,
  ident.1 = 'TfH',
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff")

MotifPlot(
  object = pbmc_chromvared,
  motifs = head(rownames(differential.activity)),
  assay = 'peaks'
)
############------ATAC----------###############
#####



### Correlation TF expression and motifs

#gene jaspar and uniprot dictionary
jaspmotifdict= read.csv('/Users/sebas/Documents/multimodal Tcell/jaspar_uniprot.csv', header = T, sep = ',')
jaspmotifdict=jaspmotifdict[,c(2,3)]
colnames(jaspmotifdict)= c('uni_id','gene')
# gene ensembl and uniprot dictionary
unidict=read.csv('/Users/sebas/Documents/multimodal Tcell/genes_uniprot_mart_export.txt', header = T, sep = '\t')
colnames(unidict)= c('uni_id','gene')
unidict=unidict[!unidict$uni_id=='',]
#motif jaspar and gene jaspar dictionary
motifjaspar=read.csv('/Users/sebas/Documents/multimodal Tcell/JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.txt', header = F, sep= '\t')
motifjaspar=motifjaspar[grep('MA',motifjaspar$V1),]
colnames(motifjaspar)=c('motif_jaspar','gene_jaspar')
#final gene motif dictionary TO USE
jaspdict= merge.data.frame(jaspmotifdict,unidict,by = 'uni_id')
colnames(jaspdict)= c('uni_id','gene_jaspar','gene_ensembl')
jaspdict= merge.data.frame(jaspdict,motifjaspar, by= 'gene_jaspar')
remove(motifjaspar)
remove(unidict)
remove(jaspmotifdict)



DefaultAssay(pmbc_1103)='RNA'
#VARIABLE GENES DF


#-------------------------------------------#

#scRNAseq imputation using ALRA 

remotes::install_github('satijalab/seurat-wrappers')
library(rsvd)
library(SeuratWrappers)
#get assay data for RNA
alra.out=SeuratWrappers::RunALRA(pmbc_1103)

gene.counts=as.data.frame(GetAssayData(alra.out,slot = 'data'))
gene.counts$gene_ensembl=rownames(gene.counts)
gene.counts_jaspar.motif= merge.data.frame(gene.counts, jaspdict, by = 'gene_ensembl')



#highly variable genes
varfeatures=as.data.frame(VariableFeatures(alra.out))
colnames(varfeatures)='gene_ensembl'
variable_tf_expressed = merge.data.frame(varfeatures,jaspdict,by = 'gene_ensembl')
gene.counts.variable_tf=merge.data.frame(gene.counts,variable_tf_expressed,by = 'gene_ensembl')




#get assay data for CHROMVAR
DefaultAssay(pbmc_chromvared)='chromvar'
motif.counts=as.data.frame(GetAssayData(pbmc_chromvared, slot = 'data'))
motif.counts$motif_jaspar=rownames(motif.counts)
motif.counts=merge.data.frame(motif.counts,jaspdict,by = 'motif_jaspar')

#checking common genes bw motifs and hvg
length(intersect(motif.counts$gene_ensembl,gene.counts.variable_tf$gene_ensembl))

#checking common bw motifs and all genes
length(intersect(motif.counts$gene_ensembl,gene.counts$gene_ensembl))


#------------------correlation

##### only gives global estimate

#second try with discrete gene expression and discrete zscore (per cell)
test.melt=melt(motif.counts)
test.melt=test.melt[,c('gene_ensembl','variable','value')]


test.melt2=melt(gene.counts_jaspar.motif)
test.melt2$imp.gene.exp=test.melt2$value
test.melt2= within.data.frame(test.melt2,rm(value))
test.melt2=test.melt2[,c('gene_ensembl','variable','imp.gene.exp')]


common.barcode=intersect(test.melt$variable,test.melt2$variable)

############ NEW ATTEMPT

gene.rna=as.data.frame(GetAssayData(pmbc_1103,slot = 'scale.data'))
DefaultAssay(pbmc_chromvared)='chromvar'
motif.counts=as.data.frame(GetAssayData(pbmc_chromvared, slot = 'data'))

gene.rna = gene.rna[,common.barcode]
motif.counts = motif.counts[,common.barcode]

jaspdict = merge(
  data.frame(genei = 1:nrow(gene.rna), gene_ensembl=rownames(gene.rna)),
  jaspdict)
jaspdict = merge(
  data.frame(motifi = 1:nrow(motif.counts), motif_jaspar=rownames(motif.counts)),
  unique(jaspdict))

jaspdict$corr = NA
for(i in 1:nrow(jaspdict)){
  jaspdict$corr[i] = cor(
    as.double(gene.rna[jaspdict$genei[i],]),
    as.double(motif.counts[jaspdict$motifi[i],]))
}
head(jaspdict)

jaspdict = jaspdict[!is.na(jaspdict$corr),]

jaspdict = jaspdict[order(jaspdict$corr),]
plot(sort(jaspdict$corr))

jaspdict
tail(jaspdict, n=30)

jaspdict[jaspdict$gene_ensembl=="GATA3",]
jaspdict[jaspdict$gene_ensembl=="FOXP3",]
jaspdict[jaspdict$gene_ensembl=="TBX21",]
jaspdict[jaspdict$gene_ensembl=="BCL6",]

jaspdict = jaspdict[,c("gene_ensembl","gene_jaspar","corr")]

jaspdict1 = data.frame(motif=jaspdict$gene_jaspar, symbol1=jaspdict$gene_ensembl, corr1=jaspdict$corr)
jaspdict2 = data.frame(motif=jaspdict$gene_jaspar, symbol2=jaspdict$gene_ensembl, corr2=jaspdict$corr)

jasphetro = merge(jaspdict1, jaspdict2)
jasphetro[jasphetro$symbol1!=jasphetro$symbol2,]


###--------Gene regulatory network analysis using SCENIC

BiocManager::install(c("AUCell", "RcisTarget",'Seurat','Signac'))

if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("aertslab/SCENIC") 
devtools::install_github("aertslab/SCopeLoomR")
library(AUCell)
library(RcisTarget)
library(GENIE3)
library(SCENIC)
library(Seurat)
library(Signac)
#files for scenic motif ranking #download db files manually, goes faster
setwd('Documents/multimodal Tcell/cisTarget_databases/')
setwd("/Users/sebas")
dbFiles <- c("https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-7species.mc9nr.feather",
             "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-7species.mc9nr.feather")
dir.create("cisTarget_databases");setwd("cisTarget_databases")


vignetteFile <- "https://raw.githubusercontent.com/aertslab/SCENIC/master/vignettes/SCENIC_Running.Rmd"
download.file(vignetteFile, "SCENIC_myRun.Rmd")



## Get data from sce object:

#----------------RNA mine
#----------------

cellInfo <- data.frame(seuratCluster=Idents(pmbc_1103))
scenicOptions <- initializeScenic(org="hgnc", dbDir="/Users/sebas/Documents/multimodal Tcell/cisTarget_databases/", nCores=10)

### Co-expression network
gene.dgc= as.matrix(GetAssayData(pmbc_1103,slot = 'data'))
genesKept <- geneFiltering(gene.dgc, scenicOptions, minCountsPerGene = 500, minSamples = 3) #same filtering criteria as for seurat
#keep only genes in db
exprMat_filtered <- gene.dgc[genesKept, ]
dim(exprMat_filtered)
rm(gene.dgc)

#computing correlations (Spearman)
runCorrelation(exprMat_filtered, scenicOptions)
#run genie3
runGenie3(exprMat_filtered, scenicOptions)

### Build and score the GRN

scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) # Toy run settings
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)


regulonsmine<- loadInt(scenicOptions, "regulons")
regulonsmine_AUC <- loadInt(scenicOptions, "aucell_regulons")
head(cbind(onlyNonDuplicatedExtended(names(regulons))))



###regulon on off

### select AUC threshold
aucellApp <- plotTsne_AUCellApp(scenicOptions, exprMat_filtered)
savedSelections <- shiny::runApp(aucellApp)
# Save the modified thresholds:
newThresholds <- savedSelections$thresholds
scenicOptions@fileNames$int["aucell_thresholds",1] <- "int/newThresholds.Rds"
saveRDS(newThresholds, file=getIntName(scenicOptions, "aucell_thresholds"))
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
# scenicOptions@settings$devType="png"
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions, exprMat_filtered)





#####Visualization and blah
# Export:
scenicOptions@fileNames$output["loomFile",] <- "output/multiometcell_SCENIC.loom"
export2loom(scenicOptions, exprMat_filtered)



####clustering based on regulon activity
nPcs <- c(10) # For toy dataset
# nPcs <- c(5,15,50)
scenicOptions@settings$seed <- 123 # same seed for all of them
# Run t-SNE with different settings:
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50))
fileNames <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50), onlyHighConf=TRUE, filePrefix="int/tSNE_oHC")
# Plot as pdf (individual files in int/):
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_", list.files("int"), value=T), value=T))
par(mfrow=c(length(nPcs), 3))
fileNames <- paste0("int/",grep(".Rds", grep("tSNE_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(fileNames, scenicOptions, showLegend=FALSE, varName="seuratCluster", cex=.5)



#checking regulon activity per ct

cellInfo <- data.frame(seuratCluster=Idents(pmbc_1103))

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$seuratCluster),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity")

topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
viewTable(topRegulators)



###cell type specific regulators based on Regulon specificty score
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=cellInfo[colnames(regulonAUC), "seuratCluster"])
rssPlot <- plotRSS(rss)
plotly::ggplotly(rssPlot$plot)
plotRSS_oneSet(rss, setName = "T cells, CD4+, Th2")
#same but binarized with the min prc of cells with that regulon
minPerc <- .1
binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")
cellInfo_binarizedCells <- cellInfo[which(rownames(cellInfo)%in% colnames(binaryRegulonActivity)),, drop=FALSE]
regulonActivity_byCellType_Binarized <- sapply(split(rownames(cellInfo_binarizedCells), cellInfo_binarizedCells$seuratCluster), 
                                               function(cells) rowMeans(binaryRegulonActivity[,cells, drop=FALSE]))
binaryActPerc_subset <- regulonActivity_byCellType_Binarized[which(rowSums(regulonActivity_byCellType_Binarized>minPerc)>0),]
ComplexHeatmap::Heatmap(binaryActPerc_subset, name="Regulon activity (%)", col = c("white","pink","red"))



#regulon activities on UMAP from seurat

dr_coords <- Embeddings(pmbc_1103, reduction="umap")

tfs <- c("FOXP1",'BCL6','NFKB1')
par(mfrow=c(2,2))
AUCell::AUCell_plotTSNE(dr_coords, cellsAUC=selectRegulons(regulonAUC, tfs), plots = "AUC")
#----------------
#----------------RNA Gosia 2022 5 day stim
#----------------
remotes::install_github("mojaveazure/seurat-disk")
SeuratDisk::Convert("/Users/sebas/Documents/otherfiles/trynka dataset 2022 5 days/stimulatedCells_highlyActiveCD4_5d_HVGs_processed.h5ad", dest = "h5seurat", overwrite = TRUE)
gosia <- SeuratDisk::LoadH5Seurat("/Users/sebas/Documents/otherfiles/trynka dataset 2022 5 days/stimulatedCells_highlyActiveCD4_5d_HVGs_processed.h5seurat")
#reduced sample size
gosia=CreateSeuratObject(counts = as.matrix(GetAssayData(gosia,slot = 'counts')), project = 'gosia',
                         assay = 'RNA',min.cells = 3,min.features = 400)

scenicOptions <- initializeScenic(org="hgnc", dbDir="/Users/sebas/Documents/multimodal Tcell/cisTarget_databases/", nCores=10)

### Co-expression network
gosia.dgc= as.matrix(GetAssayData(gosia,slot = 'counts'))
gosia.dgc= exp(gosia.dgc)

genesKept <- geneFiltering(gosia.dgc, scenicOptions,minSamples = 3, minCountsPerGene = 500)
#keep only genes in db
exprMat_filtered <- gosia.dgc[genesKept, ]
#exprMat_filtered= data.frame(exprMat_filtered)
dim(exprMat_filtered)
rm(gosia.dgc)

#computing correlations (Spearman)
runCorrelation(exprMat_filtered, scenicOptions)
#run genie3
runGenie3(exprMat_filtered, scenicOptions)

### Build and score the GRN

scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) # Toy run settings
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)


regulons <- loadInt(scenicOptions, "regulons")
regulons=data.frame(unlist(regulons))

regulons_auc <- loadInt(scenicOptions, "aucell_regulons")
head(cbind(onlyNonDuplicatedExtended(names(regulons))))
table(intersect(names(regulonsmine),names(regulonsmine_AUC)))

#----------------
#----------------RNA Gosia 2021

mtx=ReadMtx(features = '/Users/sebas/Documents/otherfiles/trynka dataset 2021 5 days effectorness cano gamez/NCOMMS-19-7936188_scRNAseq_genes.tsv',
            feature.sep = '\t', feature.column = 1,
            cells = '/Users/sebas/Documents/otherfiles/trynka dataset 2021 5 days effectorness cano gamez/NCOMMS-19-7936188_scRNAseq_barcodes.tsv',
            cell.sep = '\t',
            mtx = '/Users/sebas/Documents/otherfiles/trynka dataset 2021 5 days effectorness cano gamez/NCOMMS-19-7936188_scRNAseq_raw_UMIs.mtx')
gosia_2021=CreateSeuratObject(mtx)

meta=read.csv('/Users/sebas/Documents/otherfiles/trynka dataset 2021 5 days effectorness cano gamez//NCOMMS-19-7936188_metadata.txt', header = T,sep = '\t')


gosia_2021@meta.data=cbind(gosia_2021@meta.data,meta)
remove(meta)
remove(mtx)

#filtering using my criteria
gosia_2021[["percent.mt"]] <- PercentageFeatureSet(gosia_2021, pattern = "^MT")
VlnPlot(gosia_2021, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
plot1 <- FeatureScatter(gosia_2021, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(gosia_2021, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

gosia_2021<- subset(gosia_2021, 
                   subset = nFeature_RNA > 500 & nFeature_RNA < 5500 & percent.mt < 3 &nCount_RNA > 500 &nCount_RNA<30000)
gosia_2021<- NormalizeData(gosia_2021,normalization.method = 'LogNormalize')

gosia_2021





scenicOptions <- initializeScenic(org="hgnc", dbDir="/Users/sebas/Documents/multimodal Tcell/cisTarget_databases/", nCores=10)


### Co-expression network
gosia.dgc= as.matrix(GetAssayData(gosia_2021,slot = 'counts'))

genesKept <- geneFiltering(gosia.dgc, scenicOptions,minSamples = 3, minCountsPerGene = 500)
#keep only genes in db
exprMat_filtered <- gosia.dgc[genesKept, ]
#exprMat_filtered= data.frame(exprMat_filtered)
dim(exprMat_filtered)
rm(gosia.dgc)

#computing correlations (Spearman)
runCorrelation(exprMat_filtered, scenicOptions)
#run genie3
runGenie3(exprMat_filtered, scenicOptions)

### Build and score the GRN

scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) # Toy run settings
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)


regulons <- loadInt(scenicOptions, "regulons")
regulons=data.frame(unlist(regulons))

regulons_auc <- loadInt(scenicOptions, "aucell_regulons")
head(cbind(onlyNonDuplicatedExtended(names(regulons))))
table(intersect(names(regulonsmine),names(regulonsmine_AUC)))
