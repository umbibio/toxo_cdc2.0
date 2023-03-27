library(Seurat)
library(openxlsx)
library(ggplot2)
library(tidyverse)


source('util_funcs.R')

## Estimate empirical covariance matrix from expression data.
## Network smoothing is applied to each expression matrix with parameters alpha = 0.3 and max.iter = 3
## All data are down sampled to include 800 cells/cluster max and 2000 most variable features.
## The network in k-nn extracted from the Seurat object.

sc.rna.genes.expr.pt <- readRDS('../Input/toxo_cdc/rds/ME49_59/sc_rna_genes_expr_pt.rds')
#sc.atac.genes.expr.pt <- readRDS('../Input/toxo_cdc/rds/ME49_59/sc_atac_genes_expr_pt.rds')


## Calculate empirical covariance for scRNA
Ft.tg <- sc.rna.genes.expr.pt %>% dplyr::select(c('Sample', "GeneID", "log2.expr")) %>% 
  pivot_wider(names_from = 'Sample', values_from = 'log2.expr')

## Remove genes with 0 variance
X <- as.matrix(Ft.tg[,-1])
rownames(X) <- Ft.tg$GeneID
sds <- apply(X, 1, sd)
rm.ind <- which(sds == 0)

if(length(rm.ind) > 0){
  X <- X[-rm.ind,]
}


X.scale <- scale(t(X))
S.tg <- cov(X.scale)

saveRDS('../Input/toxo_cdc/rds/ME49_59/S_tg.rds')
