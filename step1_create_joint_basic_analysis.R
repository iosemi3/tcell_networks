if (!require("BiocManager", quietly = TRUE))

library(ggplot2)
library(scran)
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
counts <- Read10X_h5("/Users/sebas/documents/atac_tcell/from_server/filtered_feature_bc_matrix.h5")
fragpath <- "/Users/sebas/documents/atac_tcell/from_server/atac_fragments.tsv.gz"
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
pbmc <- JackStraw(pmbc_1103, red,num.replicate = 100)
pbmc <- ScoreJackStraw(pmbc_1103, dims = 1:20)
pmbc_1103<- FindNeighbors(pmbc_1103, dims = 1:30,reduction = 'PCA_RNA')
pmbc_1103<- FindClusters(pmbc_1103, resolution = 0.2)
# marker genes
cluster_mar_gene=FindAllMarkers(pmbc_1103, only.pos = F, test.use = 'negbinom')
write.table(cluster_markers,'marker_genes_cluster_tcellmultiome.tsv')

#########
ElbowPlot(pmbc_1103,reduction = 'PCA_RNA')
pmbc_1103 <- RunPCA(pmbc_1103)
VizDimLoadings(pmbc_1103, dims = 1:10, reduction = "PCA_RNA")

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










########################################################################
############------ATAC----------###############

###QC 
DefaultAssay(pmbc_1103) <- "ATAC"

pmbc_1103 <- NucleosomeSignal(pmbc_1103)

pmbc_1103 <- TSSEnrichment(pmbc_1103)

#qc plots
pmbc_1103$high.tss <- ifelse(pmbc_1103$TSS.enrichment > 2, 'High', 'Low')
pmbc_1103$nucleosome_group <- ifelse(pmbc_1103$nucleosome_signal > 4, 'NS > 4', 'NS < 4')

pmbc_1103 <- subset(
  x = pmbc_1103,
  subset = nCount_ATAC < 39000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)

##peak calling per RNA cluster using DICE ds
peaks <- CallPeaks(pmbc_1103, macs2.path = '/Users/sebas/opt/anaconda3/envs/macs2/bin/macs2',
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
DepthCor(pmbc_1103, reduction = 'SVD_ATAC')


####--Clustering --- Removing first component due to strong correlation of depth
pmbc_1103 <- RunUMAP(object = pmbc_1103, reduction = 'SVD_ATAC', dims = 2:30, reduction.name = 'UMAP_ATAC')
pmbc_1103<- FindNeighbors(object = pmbc_1103, reduction = 'SVD_ATAC', dims = 2:30)
pmbc_1103 <- FindClusters(object = pmbc_1103, verbose = FALSE, algorithm = 3, graph.name = 'ATAC_nn')
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
#1.......generalized linear model
glm_modules=names(table(modules@meta$tf))
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

#repeating the same with a gradient boost regression model (SAME AS IN SCENIC)

pmbc_1103 <- infer_grn(
  pmbc_1103,
  peak_to_gene_method = 'Signac',
  method = 'xgb'
)

pmbc_1103 <- find_modules(pmbc_1103)
modules_xgb <- NetworkModules(pmbc_1103)
#3.......gradient boost regression model
modules_xgb@meta






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
Tcell.cds@int_colData@listData$reducedDims$UMAP= pmbc_1103@reductions$umap@cell.embeddings
#plot
plot_cells(Tcell.cds,reduction_method = 'UMAP',color_cells_by = 'ident',label_groups_by_cluster = F,
           group_label_size = 5) +theme(legend.position = 'right')

Tcell.cds=learn_graph(Tcell.cds,use_partition = F)

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
pseudotime.RNA=data.frame(colData(Tcell.cds))
ggplot(pseudotime.RNA,aes(monocle3_pseudotime,reorder(ident,monocle3_pseudotime,median),fill= ident)) + geom_boxplot()
pseudotime.cellsRNA=pseudotime.RNA
pseudotime.cellsRNA$barcode=rownames(pseudotime.cellsRNA)
pseudotime.cellsRNA=pseudotime.cellsRNA[,c('barcode','monocle3_pseudotime')]
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
  object = pmbc_1103,reduction = "", pt.size = .4,
  features = c("ID2","RBPJ",
               "RNASEH2B","GYPC","CD38","ELL2")
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
