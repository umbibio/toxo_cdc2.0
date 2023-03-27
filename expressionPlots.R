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


## Tw of the datasets are mutant lines.
## ID of the KO genes
Crk2 <- 'TGGT1-218220'
Ark3 <- 'TGGT1-203010'

prod.desc  <- read.xlsx('../Input/toxo_genomics/genes/ProductDescription_GT1.xlsx')
TGGT1_ME49 <- read.xlsx('../Input/toxo_genomics/Orthologs/TGGT1_ME49 Orthologs.xlsx')


S.O.list <- readRDS('../Input/toxo_cdc/rds/S.O.intra_extra_crk2_ark3_lables_not_anchored_list.rds')
S.O.merged <- readRDS('../Input/toxo_cdc/rds/S.O.intra_extra_crk2_ark3_lables_not_anchored_merged.rds')
S.O.integrated <- readRDS('../Input/toxo_cdc/rds/S.O.intra_extra_crk2_ark3_lables_anchored_integrated.rds')



## Test plots
S.Os <- S.O.list
Idents(S.Os[[4]]) <- 'phase'
p <- DimPlot(S.Os[[4]], reduction = "pca", 
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

## Experssion plots
prod.desc[grep('AP2', prod.desc$ProductDescription),]
DefaultAssay(S.O.merged) <- 'RNA'
Idents(S.O.merged) <- 'phase'
gene.id <- 'TGGT1-233680'
gene.id <- 'TGGT1-216880'
FeaturePlot(S.O.merged, features = gene.id, split.by = 'spp', label = T)

Idents(S.O.merged) <- 'spp'
S.O.crk2 <- subset(S.O.merged, ident = 'crk2')
S.O.crk2$phase <- factor(S.O.crk2$phase, levels = c('G1.a', 'G1.b', 'S', 'M', 'C'))
Idents(S.O.crk2) <- 'phase'
S.O.ark3 <- subset(S.O.merged, ident = 'ark3')
S.O.ark3$phase <- factor(S.O.ark3$phase, levels = c('G1.a', 'G1.b', 'S', 'M', 'C'))
Idents(S.O.ark3) <- 'phase'
S.O.extra <- subset(S.O.merged, ident = 'extra')
S.O.extra$phase <- factor(S.O.extra$phase, levels = c('G1.a', 'G1.b', 'S', 'M', 'C'))
Idents(S.O.extra) <- 'phase'
S.O.intra <- subset(S.O.merged, ident = 'intra')
S.O.intra$phase <- factor(S.O.intra$phase, levels = c('G1.a', 'G1.b', 'S', 'M', 'C'))
Idents(S.O.intra) <- 'phase'

p1 <- VlnPlot(S.O.ark3, features = gene.id)
p2 <- VlnPlot(S.O.crk2, features = gene.id)
p3 <- VlnPlot(S.O.extra, features = gene.id)
p4 <- VlnPlot(S.O.intra, features = gene.id)
p1|p2|p3|p4

Idents(S.O.merged) <- 'phase'
VlnPlot(S.O.merged, features = gene.id, split.by = 'spp')

