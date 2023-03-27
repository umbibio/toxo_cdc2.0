library(Seurat)
library(ggplot2)
library(ggfortify)
library(jcolors)
require(gridExtra)
library(grid)
library(gam)
library(princurve)
library(parallel)
library(tidyverse)
library(MyEllipsefit)
library(sctransform)
library(openxlsx)
library(doParallel)
library(tidytext)
library(ggrepel)
library(dtwclust)
library(geomtextpath)


source('./util_funcs.R')
num.cores <- detectCores(all.tests = FALSE, logical = TRUE)

## Splines
sc.rna.spline.fits <- readRDS('../Input/toxo_cdc/rds/sc_rna_spline_fits_all_genes.rds')
sc.atac.spline.fits <- readRDS('../Input/toxo_cdc/rds/sc_atac_spline_fits_all_genes.rds')

marker.genes <- readRDS('../Input/toxo_cdc/rds/Intra_markers_sig.rds')
marker.genes.phase <- marker.genes %>% transmute(GeneID = gene, phase = cluster) %>% distinct()


## Turn the data into wide format (time by gene) and center & scale each gene
sc.rna.dtw.wide <- sc.rna.spline.fits %>% 
  pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
  mutate_at(vars(matches("TGME")), scale) %>%
  as.data.frame()

sc.atac.dtw.wide <- sc.atac.spline.fits %>% 
  pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
  mutate_at(vars(matches("TGME")), scale) %>%
  as.data.frame()


sc.rna.mu.scale <- sc.rna.dtw.wide %>% 
  pivot_longer(-x, names_to = 'GeneID', values_to = 'expr')

sc.atac.mu.scale <- sc.atac.dtw.wide %>% 
  pivot_longer(-x, names_to = 'GeneID', values_to = 'expr')


## Filter to include markers only
sc.rna.mu.scale <- sc.rna.mu.scale %>% dplyr::filter(GeneID %in% marker.genes$gene)
sc.atac.mu.scale <- sc.atac.mu.scale %>% dplyr::filter(GeneID %in% marker.genes$gene)



genes <- unique(sc.atac.mu.scale$GeneID)

cc.dat <- mclapply(1:length(genes), function(i){
  atac.g <- sc.atac.mu.scale %>% dplyr::filter(GeneID == genes[i]) %>% 
    transmute(t = x, y = expr) %>% arrange(t)
  rna.g <- sc.rna.mu.scale %>% dplyr::filter(GeneID == genes[i]) %>% 
    transmute(t = x, y = expr) %>% arrange(t)
  tmp <- ccf(c(atac.g$y), c(rna.g$y), plot = F)
  L <- list(Lag = tmp$lag[which.max(tmp$acf)], cc = max(tmp$acf))
  return(L)
}, mc.cores = num.cores)

cc.dat <- data.frame(GeneID = genes, Lags = unlist(lapply(cc.dat, `[[`, 1)), ccs = unlist(lapply(cc.dat, `[[`, 2)))

saveRDS(cc.dat, '../Input/toxo_cdc/rds/sc_rna_sc_atac_cross_cor_lag.rds')
ggplot(cc.dat, aes(x = ccs)) + 
  geom_histogram(aes(y = ..density..),
                 colour = 1, fill = "white") +
  geom_density(lwd = 1.2,
               linetype = 2,
               colour = 2)

