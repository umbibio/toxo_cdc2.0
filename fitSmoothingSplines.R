library(tidyverse)
library(openxlsx)
library(doParallel)
library(npreg)

## For parallel calculations
num.cores <- detectCores(all.tests = FALSE, logical = TRUE)


## Read in the data.

sc.rna.genes.expr.pt <- readRDS('../Input/Ap2Review/sc_rna_genes_expr_pt_boothroyd.rds')

## Fit smoothing splines and sample at regular time-points

genes <- unique(sc.rna.genes.expr.pt$GeneID)
## Expression
sc.rna.spline.fits <- mclapply(1:length(genes), function(i){
  tmp <- sc.rna.genes.expr.pt %>% dplyr::filter(GeneID == genes[i]) %>%
    transmute(GeneID = GeneID, x = pt.shifted.scaled, y = log2.expr)

  y <- tmp$y
  t <- tmp$x
  
   
  ### Enforcing preiodicity
  ntimes <- length(t)
  y.ext <- c(y, y, y)
  t.ext <- c(t, t+6, t+12) 
  w.ext <- rep(1, length(y.ext))
  w.ext[which(y.ext == 0)] <- 1/2
  sc.rna.sp.ext <- smooth.spline(t.ext, y.ext, lambda = 0.1, w = w.ext)
  sc.rna.sp.ext <- predict(sc.rna.sp.ext, seq(6, 12, by = 1/3)) 
  sc.rna.sp <- sc.rna.sp.ext
  sc.rna.sp$x <- sc.rna.sp$x - 6
  ######
  
  #w <- rep(1, length(y))
  #w[which(y == 0)] <- 1/2
  #sc.rna.sp <- smooth.spline(t, y, lambda = 0.1, w = w)
  #sc.rna.sp <- predict(sc.rna.sp, seq(0, 6, by = 1/3)) 
  
  #plot(t, y, col = 'red')
  #points(sc.rna.sp.ext$x, sc.rna.sp.ext$y, type = 'l', col = 'black', lwd = 2)
  #points(sc.rna.sp$x, sc.rna.sp$y, type = 'l', col = 'blue', lwd = 2)
  mu <- data.frame(x = sc.rna.sp$x, y = sc.rna.sp$y) 
  mu <- data.frame(GeneID = rep(tmp$GeneID[1], length(mu[,1])), x = mu[,1], y = mu[,2])
  return(mu)
}, mc.cores = num.cores)

sc.rna.spline.fits <- bind_rows(sc.rna.spline.fits)

saveRDS(sc.rna.spline.fits, '../Input/Ap2Review/sc_rna_spline_fits_boothroyd.rds')

i = which(genes == 'TGGT1-250800')
