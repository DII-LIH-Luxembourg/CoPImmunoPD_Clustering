###############################################################################
##
## Single Cell Temra analysis
##
## This Workflow is based on:
## https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_01_qc.html
##
## Oliver Hunewald
## 10.05.2023
##
###############################################################################

rm(list=ls())
dev.off()

library(Seurat)
library(tidyverse)
library(Matrix)
library(RCurl)
library(scales)
library(cowplot)
library(SingleCellExperiment)
library(ensembldb)
library(fs)
library(here)
library(stringr)
library(clustree)
library(dplyr)
library(readxl)
library(pheatmap)
library(rafalib)
library("data.table")

setwd("~/CoPImmunoPD")

md <- read_xlsx("metadata/metadata.xlsx") 
# load the external metadata
head(md)

# Create a Seurat object for each sample
for (i in 1:dim(md)[1]) {

  datapath <- paste0("reads/",md$file_id[i])
    
  seurat_data <- Read10X(data.dir = datapath)
  seurat_obj <- CreateSeuratObject(counts = seurat_data, 
                                   min.features = 200)
  
  metadata <- seurat_obj@meta.data
  metadata$condition <- md$condition[i]
  metadata$sample <- md$sample_id[i]
  seurat_obj@meta.data <- metadata
  
  assign(paste0("id_", md$sample_id[i]), seurat_obj)
}

table(id_HC_CCR7pCD45ROm@meta.data$nCount_RNA)
dim(id_HC_CCR7mCD45ROm)
dim(id_HC_CCR7mCD45ROp)
dim(id_PD_CCR7mCD45ROm)
dim(id_PD_CCR7mCD45ROp)

## remove data with positive ccr7 expression
# Subset on the expression level of a gene/feature
VlnPlot(id_HC_CCR7pCD45ROm, group.by = "sample", features = c("CCR7"))

id_HC_CCR7pCD45ROm <- subset(x = id_HC_CCR7pCD45ROm, subset = CCR7 >= 1)
id_HC_CCR7pCD45ROp <- subset(x = id_HC_CCR7pCD45ROp, subset = CCR7 >= 1)
id_PD_CCR7pCD45ROm <- subset(x = id_PD_CCR7pCD45ROm, subset = CCR7 >= 1)
id_PD_CCR7pCD45ROp <- subset(x = id_PD_CCR7pCD45ROp, subset = CCR7 >= 1)

id_HC_CCR7mCD45ROm <- subset(x = id_HC_CCR7mCD45ROm, subset = CCR7 < 1)
id_HC_CCR7mCD45ROp <- subset(x = id_HC_CCR7mCD45ROp, subset = CCR7 < 1)
id_PD_CCR7mCD45ROm <- subset(x = id_PD_CCR7mCD45ROm, subset = CCR7 < 1)
id_PD_CCR7mCD45ROp <- subset(x = id_PD_CCR7mCD45ROp, subset = CCR7 < 1)

VlnPlot(id_HC_CCR7pCD45ROm, group.by = "sample", features = c("CCR7"))

## ----------------------------------------------------------------------------
alldata <- merge(id_HC_CCR7mCD45ROm, c(id_HC_CCR7mCD45ROp, id_HC_CCR7pCD45ROm, 
                                          id_HC_CCR7pCD45ROp, id_PD_CCR7mCD45ROm,
                                          id_PD_CCR7mCD45ROp, id_PD_CCR7pCD45ROm, 
                                          id_PD_CCR7pCD45ROp),
                    add.cell.id = c("HC_m_m","HC_m_p","HC_p_m","HC_p_p","PD_m_m",
                                    "PD_m_p","PD_p_m","PD_p_p"))


rm(id_HC_CCR7mCD45ROm,id_HC_CCR7mCD45ROp,id_HC_CCR7pCD45ROm,id_HC_CCR7pCD45ROp,
   id_PD_CCR7mCD45ROm,id_PD_CCR7mCD45ROp,id_PD_CCR7pCD45ROm,id_PD_CCR7pCD45ROp)

gc()

as.data.frame(alldata@assays$RNA@counts[1:10, 1:2])
head(alldata@meta.data, 10)
dim(alldata)

###############################################################################
## Qualtiy Control
## Calculate QC
###############################################################################
# Manually
total_counts_per_cell <- colSums(alldata@assays$RNA@counts)
mito_genes <- rownames(alldata)[grep("^MT-", rownames(alldata))]
alldata$percent_mito <- colSums(alldata@assays$RNA@counts[mito_genes, ])/total_counts_per_cell

head(mito_genes, 10)

# Manually
ribo_genes <- rownames(alldata)[grep("^RP[SL]", rownames(alldata))]
head(ribo_genes, 10)
alldata$percent_ribo <- colSums(alldata@assays$RNA@counts[ribo_genes, ])/total_counts_per_cell

# Percentage hemoglobin genes - includes all genes starting with HB except HBP.
alldata <- PercentageFeatureSet(alldata, "^HB[^(P)]", col.name = "percent_hb")

alldata <- PercentageFeatureSet(alldata, "PECAM1|PF4", col.name = "percent_plat")

feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb")
VlnPlot(alldata, group.by = "sample", features = feats, pt.size = 0.1, ncol = 3) +
  NoLegend()

dim(alldata)

table(alldata$sample)

###############################################################################
## Filtering
## Detection-based filtering
###############################################################################
selected_c <- WhichCells(alldata, expression = nFeature_RNA > 200)
selected_f <- rownames(alldata)[Matrix::rowSums(alldata) > 3]

data.filt <- subset(alldata, features = selected_f, cells = selected_c)
dim(alldata)
dim(data.filt)
###############################################################################
## Mito/Ribo filtering
###############################################################################
selected_mito <- WhichCells(data.filt, expression = percent_mito < 0.2)
selected_ribo <- WhichCells(data.filt, expression = percent_ribo > 0.05)

# and subset the object to only keep those cells
data.filt <- subset(data.filt, cells = selected_mito)
data.filt <- subset(data.filt, cells = selected_ribo)
dim(alldata)
dim(data.filt)

table(data.filt$sample)

feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hb")

VlnPlot(data.filt, group.by = "sample", features = feats, pt.size = 0.1, ncol = 3) +
  NoLegend()

###############################################################################
# Filter Genes
###############################################################################
genes.use <- grep(pattern = "^RP[SL][[:digit:]]|^RP[[:digit:]]|^RPSA",
                  rownames(data.filt),
                  value=TRUE, invert=TRUE) 

data.filt <- data.filt[genes.use, ]
# Filter MALAT1
data.filt <- data.filt[!grepl("MALAT1", rownames(data.filt)), ]
# Filter Mitocondrial
data.filt <- data.filt[!grepl("^MT-", rownames(data.filt)), ]
# Filter Hemoglobin gene (optional if that is a problem on your data)
data.filt <- data.filt[!grepl("^HB[^(P)]", rownames(data.filt)), ]

dim(alldata)
dim(data.filt)

table(data.filt$sample)

###############################################################################
# Calculate cell-cycle scores
###############################################################################
# Before running CellCycleScoring the data need to be normalized and
# logtransformed.
data.filt = NormalizeData(data.filt)

cc.genes$s.genes <- cc.genes$s.genes[-which(cc.genes$s.genes == "MLF1IP")]
cc.genes$g2m.genes <- cc.genes$g2m.genes[-which(cc.genes$g2m.genes == "FAM64A")]
cc.genes$g2m.genes <- cc.genes$g2m.genes[-which(cc.genes$g2m.genes == "HN1")]

data.filt <- CellCycleScoring(object = data.filt, g2m.features = cc.genes$g2m.genes,
                              s.features = cc.genes$s.genes)

VlnPlot(data.filt, features = c("S.Score", "G2M.Score"), group.by = "sample",
        ncol = 4, pt.size = 0.1)

###############################################################################
# Predict doublets
###############################################################################
# library(remotes)
suppressMessages(require(DoubletFinder))

data.filt = FindVariableFeatures(data.filt, verbose = F)
data.filt = ScaleData(data.filt, vars.to.regress = c("nFeature_RNA", "percent_mito"),
                      verbose = F)
data.filt = RunPCA(data.filt, verbose = F, npcs = 20)
data.filt = RunUMAP(data.filt, dims = 1:10, verbose = F)

# Can run parameter optimization with paramSweep
sweep.res <- paramSweep_v3(data.filt) 
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE) 
bcmvn <- find.pK(sweep.stats)
barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)
# Define the expected number of doublet cellscells.
nExp <- round(ncol(data.filt) * 0.04)  # expect 4% doublets
data.filt <- doubletFinder_v3(data.filt, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)

# Name of the DF prediction can change, so extract the correct column name.
DF.name = colnames(data.filt@meta.data)[grepl("DF.classification", colnames(data.filt@meta.data))]

cowplot::plot_grid(ncol = 2, DimPlot(data.filt, group.by = "sample") + NoAxes(),
                   DimPlot(data.filt, group.by = DF.name) + NoAxes())


VlnPlot(data.filt, features = "nFeature_RNA", group.by = DF.name, pt.size = 0.1)

dim(data.filt)
data.filt = data.filt[, data.filt@meta.data[, DF.name] == "Singlet"]
dim(data.filt)

VlnPlot(data.filt, group.by = "sample", features = c("CCR7"))

saveRDS(data.filt, "./data/seurat_all_qc.rds")

ElbowPlot(data.filt, reduction = "pca", ndims = 20)

VizDimLoadings(data.filt, dims = 1:5, reduction = "pca", ncol = 5, balanced = T)

###############################################################################
## Dimension reduction
###############################################################################
suppressWarnings(suppressMessages(data.filt <- FindVariableFeatures(data.filt, 
                                                                    selection.method = "vst",
                                                                    nfeatures = 2000, 
                                                                    verbose = FALSE, 
                                                                    assay = "RNA")))

top20 <- head(VariableFeatures(data.filt), 20)
LabelPoints(plot = VariableFeaturePlot(data.filt), points = top20, repel = TRUE)

plot_grid(ncol = 3, 
          DimPlot(data.filt, reduction = "pca", group.by = "sample",dims = 1:2), 
          DimPlot(data.filt, reduction = "pca", group.by = "sample", dims = 3:4),
          DimPlot(data.filt, reduction = "pca", group.by = "sample", dims = 5:6))

VizDimLoadings(data.filt, dims = 1:5, reduction = "pca", ncol = 5, balanced = T)

ElbowPlot(data.filt, reduction = "pca", ndims = 20)

VlnPlot(data.filt, group.by = "sample", features = c("CCR7"))

###############################################################################
## UMAP
###############################################################################
data.filt <- RunUMAP(data.filt, reduction = "pca", dims = 1:20, n.components = 2, 
                     n.neighbors = 30, n.epochs = 200, min.dist = 0.3, 
                     learning.rate = 1, spread = 1, seed.use = 42)

DimPlot(data.filt, reduction = "umap", group.by = "sample")

###############################################################################
## Function for pairwise Integration
###############################################################################
do_integration <- function(sublist) {
  
  for (i in 1:length(sublist)) {
    sublist[[i]] <- NormalizeData(sublist[[i]], verbose = FALSE)
    sublist[[i]] <- FindVariableFeatures(sublist[[i]], selection.method = "vst",
                                         nfeatures = 2000, verbose = FALSE)
  }
  
  hvgs_per_dataset <- lapply(sublist, function(x) {
    x@assays$RNA@var.features
  })
  
  temp <- unique(unlist(hvgs_per_dataset))
  overlap <- sapply(hvgs_per_dataset, function(x) {
    temp %in% x
  })
  
  alldata.anchors <- FindIntegrationAnchors(object.list = sublist, dims = 1:20,
                                            reduction = "cca")
  
  alldata.int <- IntegrateData(anchorset = alldata.anchors, dims = 1:20, 
                               new.assay.name = "CCA")
  
  return(alldata.int)
  
}

###############################################################################
## Function integartion based on sub-type
###############################################################################
integrate_pairwise <- function(myDataList, subtype) {

  pht <- subtype
  
  if(pht == "m_m") {
    sublist <- c(myDataList$HC_CCR7mCD45ROm, myDataList$PD_CCR7mCD45ROm)
    alldata.int <- do_integration(sublist)
    
    return(alldata.int)
  
  } else if (pht == "p_p") {
    sublist <- c(myDataList$HC_CCR7pCD45ROp, myDataList$PD_CCR7pCD45ROp)
    alldata.int <- do_integration(sublist)
    
    return(alldata.int)
    
  } else if (pht == "p_m") {
    sublist <- c(myDataList$HC_CCR7pCD45ROm, myDataList$PD_CCR7pCD45ROm)
    alldata.int <- do_integration(sublist)
    
    return(alldata.int)
    
  } else if (pht == "m_p") {
    sublist <- c(myDataList$HC_CCR7mCD45ROp, myDataList$PD_CCR7mCD45ROp)
    alldata.int <- do_integration(sublist)
    
    return(alldata.int)
    
  } else {
    
    print("Wrong Argument!")
    
  }
}

###############################################################################
## Integration
###############################################################################
print(names(data.filt@reductions))

data.filt.list <- SplitObject(data.filt, split.by = "sample")

pht <- "m_m" # phenotype
alldata.int <- integrate_pairwise(data.filt.list, pht)
##-----------------------------------------------------------------------------
names(alldata.int@assays)
# by default, Seurat now sets the integrated assay as the default assay, so any
# operation you now perform will be on the ingegrated data.
alldata.int@active.assay
##-----------------------------------------------------------------------------
# Run Dimensionality reduction on integrated space
alldata.int <- ScaleData(alldata.int, verbose = FALSE)
alldata.int <- RunPCA(alldata.int, npcs = 20, verbose = FALSE, seed.use = 42)
alldata.int <- RunUMAP(alldata.int, dims = 1:20, seed.use = 42)

plot_grid(ncol = 2,
          DimPlot(data.filt, reduction = "umap", group.by = "sample"),
          DimPlot(alldata.int, reduction = "umap", group.by = "sample")
          )

##-----------------------------------------------------------------------------
hvgs <- unique(unlist(hvgs_per_dataset))

assaylist <- list()
genelist <- list()
for (i in 1:length(data.filt.list)) {
  assaylist[[i]] <- t(as.matrix(GetAssayData(data.filt.list[[i]], "data")[hvgs, ]))
  genelist[[i]] <- hvgs
}

lapply(assaylist, dim)

###############################################################################
## Building kNN / SNN graph
###############################################################################
# check that CCA is still the active assay
alldata.int@active.assay
alldata.int <- FindNeighbors(alldata.int, dims = 1:20, k.param = 60, prune.SNN = 1/15)

# check the names for graphs in the object.
names(alldata.int@graphs)

pheatmap(alldata.int@graphs$CCA_nn[1:200, 1:200], col = c("white", "black"), 
         border_color = "grey90", main = "KNN graph", legend = F, cluster_rows = F, 
         cluster_cols = F, fontsize = 2)


alldata.int <- FindClusters(object = alldata.int,
                            resolution = c(0.1, 0.2, 0.3, 0.4), random.seed = 42)

saveRDS(alldata.int, paste0("./data/sce_",pht,".rds"))









