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
#library(MyEllipsefit)
library(sctransform)
library(openxlsx)
library(doParallel)
library(tidytext)
library(ggrepel)
library(dtwclust)
library(geomtextpath)
library(openxlsx)


source('./util_funcs.R')

## IDs
prod.desc  <- read.xlsx('../Input/toxo_genomics/genes/ProductDescription_GT1.xlsx')
TGGT1_ME49 <- read.xlsx('../Input/toxo_genomics/Orthologs/TGGT1_ME49 Orthologs.xlsx')

prod.desc <- left_join(prod.desc, TGGT1_ME49, by = c('GeneID' = 'TGGT1'))


IMC_markers <- read.xlsx('../Input/compScBdTgPb/gene_function/IMC genes Toxoplasma.xlsx', sheet = 1)
ME49_TGGT1  <- read.xlsx('../Input/compScBdTgPb/Orthologs/TGGT1_ME49 Orthologs.xlsx')
IMC_markers <- left_join(IMC_markers, ME49_TGGT1, by = c('Gene.ID' = 'TGGT1'))

## scDATA

rna_sub <- readRDS('../Input/toxo_cdc/rds/S.O_intra_lables_pt.rds')
atac_sub <- readRDS('../Input/toxo_cdc/rds/S.O_intra_atac_lables_pt.rds')


## Splines
sc.rna.spline.fits <- readRDS('../Input/toxo_cdc/rds/sc_rna_spline_fits_all_genes.rds')
sc.atac.spline.fits <- readRDS('../Input/toxo_cdc/rds/sc_atac_spline_fits_all_genes.rds')

## Turn the data into wide format (time by gene) and center & scale each gene
sc.rna.dtw.wide <- sc.rna.spline.fits %>% 
  pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
  mutate_at(vars(matches("TGME")), ~scale(., center = T, scale = T)) %>%
  as.data.frame()

## Remove NAs (genes with 0 expression)
na.ind <- which(apply(sc.rna.dtw.wide, 2, function(x) any(is.na(x))))
if(length(na.ind)){
  sc.rna.dtw.wide <- sc.rna.dtw.wide[,-na.ind]
}

sc.atac.dtw.wide <- sc.atac.spline.fits %>% 
  pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
  mutate_at(vars(matches("TGME")), ~scale(., center = T, scale = T)) %>%
  as.data.frame()

na.ind <- which(apply(sc.atac.dtw.wide, 2, function(x) any(is.na(x))))
if(length(na.ind)){
  sc.atac.dtw.wide <- sc.atac.dtw.wide[,-na.ind]
}


sc.rna.mu.scale <- sc.rna.dtw.wide %>% 
  pivot_longer(-x, names_to = 'GeneID', values_to = 'expr')

sc.atac.mu.scale <- sc.atac.dtw.wide %>% 
  pivot_longer(-x, names_to = 'GeneID', values_to = 'expr')

sc.rna.atac.mu.scale <- inner_join(sc.rna.mu.scale, sc.atac.mu.scale , by = c('x', 'GeneID'))
colnames(sc.rna.atac.mu.scale) <- c('time', 'GeneID', 'scRNA', 'scATAC')
prod.desc$TGME49 <- gsub('_', '-', prod.desc$TGME49)
sc.rna.atac.mu.scale <- left_join(sc.rna.atac.mu.scale, prod.desc, by = c('GeneID' = 'TGME49'))

#sc.rna.hc_dtw <- dtwClustCurves(sc.rna.dtw.wide[,-1], nclust = 8L)
#sc.atac.hc_dtw <- dtwClustCurves(sc.atac.dtw.wide[,-1], nclust = 8L)

#saveRDS(sc.rna.hc_dtw, "../Input/toxo_cdc/rds/sc_rna_hc_dtw.rds")
#saveRDS(sc.atac.hc_dtw, "../Input/toxo_cdc/rds/sc_atac_hc_dtw.rds")

sc.rna.d <- as.matrix(dist(t(as.matrix(sc.rna.dtw.wide[,-1]))))
sc.atac.d <- as.matrix(dist(t(as.matrix(sc.atac.dtw.wide[,-1]))))

## Identify matching genes.
Crk2 <- 'TGME49-218220'
Ark3 <- 'TGME49-203010'
CycP1 <- 'TGME49-267580'
CycP2 <- 'TGME49-305330'
Cyc5 <- 'TGME49-293280'
my.features <- c(Crk2, Ark3, CycP1, CycP2, Cyc5)


conserved_TFs <- c('TGME49-286710',
                   'TGME49-206650',
                   'TGME49-236840',
                   'TGME49-252310',
                   'TGME49-203380',
                   'TGME49-201220',
                   'TGME49-242320')


my.features <- gsub('_', '-', IMC_markers$TGME49)

my.features <- 'TGME49-203380'
my.features <- 'TGME49-250800'


## New list from MJ, MNK1
MNK1  <- 'TGME49-275610'
IMC29 <- 'TGME49-243200'
IMC32 <- 'TGME49-232150'
BCC0  <- 'TGME49-294860'
BCC3  <- 'TGME49-311770'
FBXO1 <- 'TGME49-310930'
AC9   <- 'TGME49-246950'

my.features <- c(MNK1, IMC29, IMC32, BCC0, BCC3, FBXO1, AC9)
my.names <- c('MNK1', 'IMC29', 'IMC32', 'BCC0', 'BCC3', 'FBXO1', 'AC9')
sim.genes.rna <- lapply(1:length(my.features), function(i){
  ind <- which(colnames(sc.rna.d) == my.features[i])
  sim.genes <- quantile(sc.rna.d[ind, ], probs = 0.01)
  sim.genes.ind <- which(sc.rna.d[ind, ] < sim.genes)
  colnames(sc.rna.d)[sim.genes.ind]
})


sim.genes.atac <- lapply(1:length(my.features), function(i){
  ind <- which(colnames(sc.atac.d) == my.features[i])
  sim.genes <- quantile(sc.atac.d[ind, ], probs = 0.01)
  sim.genes.ind <- which(sc.atac.d[ind, ] < sim.genes)
  colnames(sc.atac.d)[sim.genes.ind]
})


sim.genes.rna.atac <- lapply(1:length(my.features), function(i){
  intersect(sim.genes.rna[[i]], sim.genes.atac[[i]])
})

plot_rna_atac_trends_V2 <- function(genes, title = ''){
  sc.rna.atac.mu.scale.sub <- sc.rna.atac.mu.scale %>% dplyr::filter(GeneID %in% genes) %>% 
    pivot_longer(-c('time', "GeneID", "GeneID.y", "ProductDescription"), 
                 names_to = 'data', values_to = 'normExpr') 
  
  sc.rna.atac.mu.scale.sub$label <- NA
  sc.rna.atac.mu.scale.sub$label[which(sc.rna.atac.mu.scale.sub$time == 3)] <-
    sc.rna.atac.mu.scale.sub$GeneID[which(sc.rna.atac.mu.scale.sub$time == 3)]
  
  p  <- ggplot(sc.rna.atac.mu.scale.sub, aes(x= time,y=normExpr)) +
    geom_path(aes(color = GeneID),alpha = 0.8, size = 0.8)+ 
    theme_bw(base_size = 14) +
    theme(legend.position = "right") +
    ylab('normExpr') + xlab('Time') +
    theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
    theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
    theme(strip.background = element_rect(colour="black", fill="white",
                                          size=0.5, linetype="solid")) +
    theme(strip.text = element_text(size = 14, face="bold", angle = 0)) + 
    
    coord_cartesian(xlim = c(0,6.5)) + 
    geom_text_repel(aes(label = label), size = 2.5, fontface = "bold",
                    box.padding = unit(0.6, "lines"),
                    max.overlaps = 300,
                    #segment.angle = 180,
                    nudge_x = 0.25, 
                    nudge_y = 0.25,
                    hjust=0.25,
                    #nudge_x=0.25, 
                    segment.size = 0.1,
                    na.rm = TRUE)+ 
    facet_grid(.~data, scales = 'free', space = 'free') +
    
    
    ggtitle(title) +
    theme(
      plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
      axis.title.x = element_text(size=14, face="bold", hjust = 1),
      axis.title.y = element_text(size=14, face="bold")
    ) + 
    theme(#legend.position = c(0.15, 0.85),
      legend.position = 'none',
      legend.title = element_text(colour="black", size=12, 
                                  face="bold"),
      legend.text = element_text(colour="black", size=12, 
                                 face="bold"))
  
  
  return(p)
  
}


for(i in 1:length(sim.genes.rna.atac)){
  if(length(sim.genes.rna.atac[[i]]) == 0){
    next
  }
  
  p <- plot_rna_atac_trends_V2(sim.genes.rna.atac[[i]], title = my.names[i])
  f.n <- paste("../Output/toxo_cdc/figures/MNK1/", my.features[i], '.pdf', sep = '')
  ggsave(filename=f.n,
         plot=p,
         width = 8, height = 6,
         units = "in"
  )
  
  tab <- prod.desc %>% dplyr::filter(TGME49 %in% sim.genes.rna.atac[[i]])
  f.n <- paste("../Output/toxo_cdc/figures/MNK1/",my.features[i], '.xlsx', sep = '')
  write.xlsx(tab, f.n)
}

p1 <- plot_rna_atac_trends_V2(sim.genes.rna.atac[[1]])
plot(p1)
prod.desc %>% dplyr::filter(TGME49 %in% sim.genes.rna.atac[[1]])


ggsave(filename="../Output/toxo_cdc/figures/AP2_rna_atac_clusters.pdf",
       plot=p1,
       width = 10, height = 8,
       units = "in"
)


## Filter out hypotheticals
non.hype <- prod.desc %>% dplyr::filter(TGME49 %in% sim.genes.rna.atac[[1]]) %>% 
  dplyr::filter(ProductDescription != "hypothetical protein")

sim.sim.genes.rna.atac.no.hyp <- non.hype$TGME49
p1 <- plot_rna_atac_trends_V2(sim.sim.genes.rna.atac.no.hyp)
plot(p1)
prod.desc %>% dplyr::filter(TGME49 %in% sim.genes.rna.atac[[1]])


ggsave(filename="../Output/toxo_cdc/figures/AP2_rna_atac_clusters.pdf",
       plot=p1,
       width = 10, height = 8,
       units = "in"
)



p <- plot_trends("TGME49-216660")
plot(p)
plot_trends("TGME49-267580")



#sig.AP2s <- inner_join(DEG.sig, AP2s, by = c('gene' = 'TGME49'))
## Manually add AP2XI-5
sig.AP2s <- AP2s[AP2s$TGME49 %in% c(DEG.sig$gene, "TGME49-216220"),]

saveRDS(sig.AP2s, '../Input/toxo_cdc/rds/sig_AP2s.rds')

## Clustering AP2s

sc.rna.AP2 <- sc.rna.dtw.wide[,colnames(sc.rna.dtw.wide) %in% sig.AP2s$TGME49 ]
sc.atac.AP2 <- sc.atac.dtw.wide[,colnames(sc.atac.dtw.wide) %in% sig.AP2s$TGME49 ]

sc.rna.AP2.markers.hc_dtw <- dtwClustCurves(sc.rna.AP2, nclust = 3L)
sc.atac.AP2.markers.hc_dtw <- dtwClustCurves(sc.atac.AP2, nclust = 3L)

plot(sc.rna.AP2.markers.hc_dtw, type = 'sc')
plot(sc.atac.AP2.markers.hc_dtw, type = 'sc')


## GGplot cluster graphs

sc.rna.long.AP2 <- inner_join(sc.rna.mu.scale, sig.AP2s, by = c('GeneID' = 'TGME49')) 
sc.rna.long.AP2 <- sc.rna.long.AP2 %>% 
  transmute(time = x, GeneID = GeneID, normExpr = expr, Name = Ap2Name)  %>% distinct()

sc.atac.long.AP2 <- inner_join(sc.atac.mu.scale, sig.AP2s, by = c('GeneID' = 'TGME49')) 
sc.atac.long.AP2 <- sc.atac.long.AP2 %>% 
  transmute(time = x, GeneID = GeneID, normExpr = expr, Name = Ap2Name) %>% distinct()


sc.rna.clust.info <- data.frame(GeneID = colnames(sc.rna.AP2), cluster = cutree(sc.rna.AP2.markers.hc_dtw, k = 3))
sc.atac.clust.info <- data.frame(GeneID = colnames(sc.atac.AP2), cluster = cutree(sc.atac.AP2.markers.hc_dtw, k = 3))

sc.rna.long.clust <- inner_join(sc.rna.long.AP2, sc.rna.clust.info, by = 'GeneID')
sc.atac.long.clust <- inner_join(sc.atac.long.AP2, sc.atac.clust.info, by = 'GeneID')

sc.rna.sc.atac.joint <- inner_join(sc.rna.long.clust, sc.atac.long.clust, 
                                   by = c("time", "GeneID", "Name"))
colnames(sc.rna.sc.atac.joint) <- c("time", "GeneID", "scRNA", "Name", "cluster.RNA", "scATAC", "cluster.ATAC")
sc.rna.sc.atac.joint.long <- sc.rna.sc.atac.joint %>% 
  pivot_longer(-c('time', "GeneID", "Name", "cluster.RNA", "cluster.ATAC"), 
               names_to = 'data', values_to = 'normExpr') 


sc.rna.sc.atac.joint.long$label <- NA
sc.rna.sc.atac.joint.long$label[which(sc.rna.sc.atac.joint.long$time == 3)] <-
  sc.rna.sc.atac.joint.long$Name[which(sc.rna.sc.atac.joint.long$time == 3)]

# sc.rna.sc.atac.joint.long$label <- NA
# AP2s.clusts <- sc.rna.sc.atac.joint.long %>% group_by(Name) %>% summarise(cluster.RNA = cluster.RNA[1]) %>% 
#   ungroup() %>% group_by(cluster.RNA) %>% mutate(lab.time = 1:n() %% 6)
# 
# sc.rna.sc.atac.joint.long$lab.time <- floor(sc.rna.sc.atac.joint.long$time) 
# sc.rna.sc.atac.joint.long <- left_join(sc.rna.sc.atac.joint.long, AP2s.clusts, by = c('Name', 'lab.time', 'cluster.RNA'))
# sc.rna.sc.atac.joint.long$label[sc.rna.sc.atac.joint.long$time == sc.rna.sc.atac.joint.long$lab.time] <-
#   sc.rna.sc.atac.joint.long$Name[sc.rna.sc.atac.joint.long$time == sc.rna.sc.atac.joint.long$lab.time]

sc.rna.sc.atac.joint.long$cluster.RNA <- paste('C', sc.rna.sc.atac.joint.long$cluster.RNA)

saveRDS(sc.rna.sc.atac.joint.long, '../Input/toxo_cdc/rds/AP2_sc_rna_sc_atac_joint_dtw_clust.rds')


plot_rna_atac_trends <- function(sc.rna.sc.atac.joint.long.sub){
  p  <- ggplot(sc.rna.sc.atac.joint.long.sub, aes(x= time,y=normExpr)) +
    geom_path(aes(color = Name,),alpha = 0.8, size = 0.8)+ 
    theme_bw(base_size = 14) +
    theme(legend.position = "right") +
    ylab('normExpr') + xlab('Time') +
    theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
    theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
    theme(strip.background = element_rect(colour="black", fill="white",
                                          size=0.5, linetype="solid")) +
    theme(strip.text = element_text(size = 14, face="bold", angle = 0)) + 
    
    coord_cartesian(xlim = c(0,6.5)) + 
    geom_text_repel(aes(label = label), size = 2.5, fontface = "bold",
                    box.padding = unit(0.6, "lines"),
                    max.overlaps = 300,
                    #segment.angle = 180,
                    nudge_x = 0.25, 
                    nudge_y = 0.25,
                    hjust=0.25,
                    #nudge_x=0.25, 
                    segment.size = 0.1,
                    na.rm = TRUE)+ 
    facet_grid(cluster.RNA~data, scales = 'free', space = 'free') +
    
    
    #ggtitle(titles[i]) +
    theme(
      plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
      axis.title.x = element_text(size=14, face="bold", hjust = 1),
      axis.title.y = element_text(size=14, face="bold")
    ) + 
    theme(#legend.position = c(0.15, 0.85),
      legend.position = 'none',
      legend.title = element_text(colour="black", size=12, 
                                  face="bold"),
      legend.text = element_text(colour="black", size=12, 
                                 face="bold"))
  
  
  return(p)
  
}

p1 <- plot_rna_atac_trends(sc.rna.sc.atac.joint.long)

plot(p1)

ggsave(filename="../Output/toxo_cdc/figures/AP2_rna_atac_clusters.pdf",
       plot=p1,
       width = 10, height = 8,
       units = "in"
)


### Subsets of AP2s
sig.AP2s$Ap2Name
my.AP2s <- c("AP2VIIb-3", "AP2IX-4" ) ## Know to interact with XII-2
sc.rna.sc.atac.joint.long.sub <- sc.rna.sc.atac.joint.long %>% dplyr::filter(Name %in% my.AP2s)
p1 <- plot_rna_atac_trends(sc.rna.sc.atac.joint.long.sub)
plot(p1)


my.AP2s <- c("AP2XI-3" ) ## 
sc.rna.sc.atac.joint.long.sub <- sc.rna.sc.atac.joint.long %>% dplyr::filter(Name %in% my.AP2s)
p1 <- plot_rna_atac_trends(sc.rna.sc.atac.joint.long.sub)
plot(p1)

## TGME49_303050 SWI5
