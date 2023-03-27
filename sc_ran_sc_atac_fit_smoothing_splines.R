library(tidyverse)
library(openxlsx)
library(doParallel)
library(npreg)

## For parallel calculations
num.cores <- detectCores(all.tests = FALSE, logical = TRUE)


## Read in the data.

sc.rna.genes.expr.pt <- readRDS('../Input/toxo_cdc/rds/sc_rna_genes_expr_pt.rds')
sc.atac.genes.expr.pt <- readRDS('../Input/toxo_cdc/rds/sc_atac_genes_expr_pt.rds')

#Intra.markers.sig <- readRDS('../Input/toxo_cdc/rds/Intra_markers_sig.rds')
#DEG.sig <- readRDS('../Input/toxo_cdc/rds/Intra_DEGs_up_and_down_sig.rds') ## both up & down


## Common genes between the data sets
comm.genes <- intersect(unique(sc.rna.genes.expr.pt$GeneID), unique(sc.atac.genes.expr.pt$GeneID))

#comm.genes <- comm.genes[comm.genes %in% Intra.markers.sig$gene]
#comm.genes <- comm.genes[comm.genes %in% DEG.sig$gene]


# smoothFilter <- function(y){
#   window.length <- floor(length(y) / 10)
#   filt.y <- rep(T, length(y))
#   segs <- list()
#   prop.non.zeros <- {}
#   point.density <- {}
#   for(i in 1:10){
#     #vd <- density(y[((i - 1) * window.length + 1):(i * window.length)])
#     #vplot(d)
#     # 
#     segs <- c(segs, list(((i - 1) * window.length + 1):(i * window.length)))
#     prop.non.zero <- sum(y[segs[[i]]] > 0) / 
#       length(segs[[i]])
#     
#     prop.non.zeros <- c(prop.non.zeros, prop.non.zero)
#     point.density <- c(point.density, length(y[segs[[i]]])/length(y))
#     #print(paste(t[seg][1], t[seg][length(seg)]))
#     #print(prop.non.zero)
#     #Sys.sleep(0.8)
#   }
#   
#   keep <- quantile(prop.non.zeros, prob = 0.10)
#   for(i in 1:10){
#     if(prop.non.zeros[i] <= keep){
#       filt.y[segs[[i]]] <- F
#     }
#   }
#   
#   
#   return(which(y == 0 & filt.y))
# }
# 
# 
# density.fit <- function(t, y){
#   total.points <- length(t)
#   
#   ## 20 min long windows
#   window.length <- floor(20 * total.points / (6 * 60))
#   window.time <- 20 / 60
#   ## 5 min stride
#   stride <- floor(5 * total.points / (6 * 60))
#   stride.time <- 5 / 60
#   total.slides <- floor(((total.points - window.length) / stride) + 1)
#   
#   segs <- list()
#   prop.non.zeros <- {}
#   point.density <- {}
#   w <- matrix(0, nrow =total.points, ncol = total.points)
#   for(i in 1:total.slides){
#     strt <- (i - 1) * stride + 1
#     stp <- (i - 1) * stride + 1 + window.length
#     
#     strt.t <- (i - 1) * stride.time 
#     stp.t <- (i - 1) * stride.time + window.time
#     
#     segs <- c(segs, list(strt:stp))
#     
#     ind <- t >=  strt.t & t <= stp.t
#     prop.non.zero <- sum(y[ind] > 0) / sum(ind)
#     prop.non.zeros <- c(prop.non.zeros, prop.non.zero)
#     point.density <- c(point.density, sum(ind)/total.points)
#     w[, ind] <- point.density[i] 
#   }
#   
#   #weights[, stp:ncol(weights)] <- weights[stp, stp]
#   
#   w <- colSums(w) / nrow(w)
#   w <- min(y) + ((max(y) - min(y))/(max(w) - min(w))) * w
#   return(w)
# }

do.filt <- F
## Fit smoothing splines and sample at regular time-points

## Expression
sc.rna.spline.fits <- mclapply(1:length(comm.genes), function(i){
  tmp <- sc.rna.genes.expr.pt %>% dplyr::filter(GeneID == comm.genes[i]) %>%
    transmute(GeneID = GeneID, x = pt.shifted.scaled, y = log2.expr)
 
  if(do.filt){
    ind <- smoothFilter(tmp$y)
    tmp2 <- tmp[-ind,]
    y <- tmp2$y
    t <- tmp2$x
  }else{
    y <- tmp$y
    t <- tmp$x
    
  }
  #sc.rna.sp <- ss(t, y, periodic = T, lambda = 1e-4)
  w <- rep(1, length(y))
  w[which(y == 0)] <- 1/2
  sc.rna.sp <- smooth.spline(t, y, lambda = 0.1, w = w)
  
  sc.rna.sp <- predict(sc.rna.sp, seq(0, 6, by = 1/3)) 
  #sc.rna.pp <- fitPsplines(t, y)
  #plot(tmp$x, tmp$y)
  #points(t, y, col = 'red')
  #points(sc.rna.sp$x, sc.rna.sp$y, type = 'l', col = 'red')
  #points(sc.rna.sp$x, sc.rna.sp$y, type = 'l', col = 'blue')
  #points(sc.rna.pp$x, sc.rna.sp$y, type = 'l', col = 'green')
  mu <- data.frame(x = sc.rna.sp$x, y = sc.rna.sp$y) 
  mu <- data.frame(GeneID = rep(tmp$GeneID[1], length(mu[,1])), x = mu[,1], y = mu[,2])
  return(mu)
}, mc.cores = num.cores)

sc.rna.spline.fits <- bind_rows(sc.rna.spline.fits)

#saveRDS(sc.rna.spline.fits, '../Input/toxo_cdc/rds/sc_rna_spline_fits_cell_cycle_genes.rds')
saveRDS(sc.rna.spline.fits, '../Input/toxo_cdc/rds/sc_rna_spline_fits_all_genes.rds')

## Access
sc.atac.spline.fits <- mclapply(1:length(comm.genes), function(i){
  tmp <- sc.atac.genes.expr.pt %>% dplyr::filter(GeneID == comm.genes[i]) %>%
    transmute(GeneID = GeneID, x = pt.shifted.scaled, y = log2.expr)
  
  if(do.filt){
    ind <- smoothFilter(tmp$y)
    tmp2 <- tmp[-ind,]
    y <- tmp2$y
    t <- tmp2$x
  }else{
    y <- tmp$y
    t <- tmp$x
    
  }
  #sc.atac.sp <- ss(t, y, periodic = T, lambda = 1e-4)
  w <- rep(1, length(y))
  w[which(y == 0)] <- 1/2
  sc.atac.sp <- smooth.spline(t, y, lambda = 0.1, w = w)
  sc.atac.sp <- predict(sc.atac.sp, seq(0, 6, by = 1/3)) 
  
  mu <- data.frame(x = sc.atac.sp$x, y = sc.atac.sp$y)
  mu <- data.frame(GeneID = rep(tmp$GeneID[1], length(mu[,1])), x = mu[,1], y = mu[,2])
  return(mu)
}, mc.cores = num.cores)

sc.atac.spline.fits <- bind_rows(sc.atac.spline.fits)

#saveRDS(sc.atac.spline.fits, '../Input/toxo_cdc/rds/sc_atac_spline_fits_cell_cycle_genes.rds')
saveRDS(sc.atac.spline.fits, '../Input/toxo_cdc/rds/sc_atac_spline_fits_all_genes.rds')

