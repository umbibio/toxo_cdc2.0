library(Seurat)
library(openxlsx)
library(ggplot2)
library(ggfortify)
library(jcolors)
require(gridExtra)
library(grid)
library(matrixStats)
library(tidyverse)
library(RColorBrewer)
library(parallel)
library(openxlsx)
library(plotly)


source('./util_funcs.R')


num.cores <- detectCores(all.tests = FALSE, logical = TRUE)

## Tw of the datasets are mutant lines.
## ID of the KO genes
Crk2 <- 'TGGT1-218220'
Ark3 <- 'TGGT1-203010'

## Count files
intra.file.csv <- "../Input/toxo_scRNA_MJ/RH.intra.expr.csv"
extra.file.csv <- "../Input/toxo_scRNA_MJ/RH.extra.expr.csv"
crk2.file.csv  <- "../Input/toxo_scRNA_MJ/RH.crk2.expr.csv"
ark3.file.csv  <- "../Input/toxo_scRNA_MJ/RH.ark3.expr.csv"

## IDs
prod.desc  <- read.xlsx('../Input/toxo_genomics/genes/ProductDescription_GT1.xlsx')
TGGT1_ME49 <- read.xlsx('../Input/toxo_genomics/Orthologs/TGGT1_ME49 Orthologs.xlsx')


getExpr <- function(in.file, TGGT1_ME49){
  file.counts <- read.csv(in.file)
  genes <- file.counts$X
  ind <- which(genes %in% TGGT1_ME49$TGME49)
  file.counts <- file.counts[ind, ]
  genes <- genes[ind]
  genes <- TGGT1_ME49$TGGT1[match(genes, TGGT1_ME49$TGME49)]
  expr <- file.counts[,-1]
  rownames(expr) <- genes
  
  return(expr)
}

intra.counts <- getExpr(intra.file.csv, TGGT1_ME49)
extra.counts <- getExpr(extra.file.csv, TGGT1_ME49)
crk2.counts  <- getExpr(crk2.file.csv, TGGT1_ME49)
ark3.counts  <- getExpr(ark3.file.csv, TGGT1_ME49)


## individual Seurat objects
feats <- c("nFeature_RNA","nCount_RNA")

# Intra
S.O.intra <- CreateSeuratObject(counts = intra.counts)
S.O.intra$orig.ident <- 'intra'
VlnPlot(S.O.intra, features = feats, pt.size = 0.1,ncol = 2) + NoLegend()
FeatureScatter(S.O.intra, "nCount_RNA", "nFeature_RNA",  pt.size = 0.5)
selected_c <- WhichCells(S.O.intra, expression = nFeature_RNA > 100 & nFeature_RNA < 950)
selected_f <- rownames(S.O.intra)[ Matrix::rowSums(S.O.intra) > 5]
selected_f <- unique(c(selected_f, Crk2, Ark3)) ## Make sure Crk2 and Ark 3 are selected
S.O.intra  <- subset(S.O.intra, features=selected_f, cells=selected_c)

saveRDS(S.O.intra, '../Input/toxo_cdc/rds/S.O_intra_not_down_sample.rds')

# Extra
S.O.extra <- CreateSeuratObject(counts = extra.counts)
S.O.extra$orig.ident <- 'extra'
VlnPlot(S.O.extra, features = feats, pt.size = 0.1,ncol = 2) + NoLegend()
FeatureScatter(S.O.extra, "nCount_RNA", "nFeature_RNA",  pt.size = 0.5)
selected_c <- WhichCells(S.O.extra, expression = nFeature_RNA > 100 & nFeature_RNA < 950)
selected_f <- rownames(S.O.extra)[ Matrix::rowSums(S.O.extra) > 5]
selected_f <- unique(c(selected_f, Crk2, Ark3)) ## Make sure Crk2 and Ark 3 are selected
S.O.extra  <- subset(S.O.extra, features=selected_f, cells=selected_c)

# CRK2
S.O.crk2  <- CreateSeuratObject(counts = crk2.counts)
S.O.crk2$orig.ident <- 'crk2'
VlnPlot(S.O.crk2, features = feats, pt.size = 0.1,ncol = 2) + NoLegend()
FeatureScatter(S.O.crk2, "nCount_RNA", "nFeature_RNA",  pt.size = 0.5)
selected_c <- WhichCells(S.O.crk2, expression = nFeature_RNA > 100 & nFeature_RNA < 900)
selected_f <- rownames(S.O.crk2)[ Matrix::rowSums(S.O.crk2) > 5]
selected_f <- unique(c(selected_f, Crk2, Ark3)) ## Make sure Crk2 and Ark 3 are selected
S.O.crk2   <- subset(S.O.crk2, features=selected_f, cells=selected_c)

# ARK3
S.O.ark3  <- CreateSeuratObject(counts = ark3.counts)
S.O.ark3$orig.ident <- 'ark3'
VlnPlot(S.O.ark3, features = feats, pt.size = 0.1,ncol = 2) + NoLegend()
FeatureScatter(S.O.ark3, "nCount_RNA", "nFeature_RNA",  pt.size = 0.5)
selected_c <- WhichCells(S.O.ark3, expression = nFeature_RNA > 100 & nFeature_RNA < 1000)
selected_f <- rownames(S.O.ark3)[ Matrix::rowSums(S.O.ark3) > 5]
selected_f <- unique(c(selected_f, Crk2, Ark3)) ## Make sure Crk2 and Ark 3 are selected
S.O.ark3   <- subset(S.O.ark3, features=selected_f, cells=selected_c)


S.O.list <- list(intra = S.O.intra, extra = S.O.extra, crk2 = S.O.crk2, ark3 = S.O.ark3)

## Downsample to 6000 cells
set.seed(100)
S.O.list <- mclapply(S.O.list, function(S.O){
  S.O <- subset(x = S.O, downsample = 8000)
  
}, mc.cores = num.cores)



## Individually process the data and transfer labels
## Boothroyed data
S.O.tg.boothroyd <- readRDS('../Input/boothroyd_sc_all_data/rds/S.O.tg_RH_boothroyd.rds')

## split the data, process each, transfer the lables
S.Os <- mclapply(S.O.list, function(S.O){
  S.O <- prep_S.O(S.O, res = 0.4)
  anchors <- FindTransferAnchors(reference = S.O.tg.boothroyd, query = S.O, dims = 1:30)
  predictions <- TransferData(anchorset = anchors, refdata = S.O.tg.boothroyd@meta.data$phase,dims = 1:30)
  predictions$phase <- predictions$predicted.id
  #predictions$phase[which(predictions$prediction.score.max < 0.7)] <- 'NA'
  S.O <- AddMetaData(object = S.O, metadata = predictions)
  return(S.O)
}, mc.cores = num.cores)

spps <- names(S.Os)

S.Os <- lapply(1:length(S.Os), function(i){
  S.Os[[i]]@meta.data$spp <- spps[i]
  S.Os[[i]]
})

saveRDS(S.Os, '../Input/toxo_cdc/rds/S.O.intra_extra_crk2_ark3_lables_not_anchored_list.rds')
saveRDS(S.Os[[1]], '../Input/toxo_cdc/rds/S.O.intra_lables.rds')


## Test plots
Idents(S.Os[[1]]) <- 'phase'
p <- DimPlot(S.Os[[1]], reduction = "pca", 
             #group.by = "cell", 
             #split.by = 'spp',
             pt.size = 1,
             #shape.by='spp',
             label = F, label.size = 4) + #NoLegend() + 
  theme(panel.spacing = unit(0.5, "lines")) + 
  theme(axis.text.x = element_text(face="bold", size=12, angle=0)) +
  theme(axis.text.y = element_text(face="bold", size=12, angle=0)) +
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )


plot(p)


### Merging data
alldata.lab <- merge(S.Os[[1]], S.Os[2:4], add.cell.ids=c("intra","extra","crk2","ark3"))
alldata.lab$spp <- alldata.lab$orig.ident

alldata.lab <- prep_S.O(alldata.lab, res = 0.4)


## Test plots
Idents(alldata.lab) <- 'spp'
p <- DimPlot(alldata.lab, reduction = "umap", 
             #group.by = "cell", 
             #split.by = 'spp',
             pt.size = 1,
             #shape.by='spp',
             label = F, label.size = 4) + #NoLegend() + 
  theme(panel.spacing = unit(0.5, "lines")) + 
  theme(axis.text.x = element_text(face="bold", size=12, angle=0)) +
  theme(axis.text.y = element_text(face="bold", size=12, angle=0)) +
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )


plot(p)


saveRDS(alldata.lab, '../Input/toxo_cdc/rds/S.O.intra_extra_crk2_ark3_lables_not_anchored_merged.rds')

### Anchoring data to intra
## shared PCA/UMAP

S.O.list <- lapply(X = S.Os, FUN = function(x) {
  ## Extract the count data
  
  ## extract the count data from each as.matrix(S.O.list[[1]][["RNA"]]@data)
  ## Replace genes with Bdiv orthologous when needed
  ## recreate the new Seurat object.
  DefaultAssay(x) <- 'RNA'
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 6000)
})

features <- SelectIntegrationFeatures(object.list = S.O.list, nfeatures = 6000)
spps <- names(S.Os)
reference_dataset <- 1
anchors <- FindIntegrationAnchors(object.list = S.O.list, 
                                  anchor.features = features, reference = reference_dataset)
S.O.integrated <- IntegrateData(anchorset = anchors)
# switch to integrated assay. Make sure to set to RNA for Differential Expression
DefaultAssay(S.O.integrated) <- "integrated"
S.O.integrated <- ScaleData(object = S.O.integrated, verbose = FALSE)
S.O.integrated <- RunPCA(S.O.integrated, features = VariableFeatures(object = S.O.integrated))
S.O.integrated <- FindNeighbors(S.O.integrated, dims = 1:10, reduction = 'pca')
S.O.integrated <- FindClusters(S.O.integrated, resolution = 0.4)
S.O.integrated <- RunUMAP(S.O.integrated, dims = 1:13)

## Test plots
Idents(S.O.integrated) <- 'phase'
p <- DimPlot(S.O.integrated, reduction = "umap", 
             #group.by = "cell", 
             split.by = 'spp',
             pt.size = 1,
             #shape.by='spp',
             label = F, label.size = 4) + #NoLegend() + 
  theme(panel.spacing = unit(0.5, "lines")) + 
  theme(axis.text.x = element_text(face="bold", size=12, angle=0)) +
  theme(axis.text.y = element_text(face="bold", size=12, angle=0)) +
  theme(
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )


plot(p)

saveRDS(S.O.integrated, '../Input/toxo_cdc/rds/S.O.intra_extra_crk2_ark3_lables_anchored_integrated.rds')
