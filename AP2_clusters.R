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

## IDs
prod.desc  <- read.xlsx('../Input/toxo_genomics/genes/ProductDescription_GT1.xlsx')
TGGT1_ME49 <- read.xlsx('../Input/toxo_genomics/Orthologs/TGGT1_ME49 Orthologs.xlsx')

## Up & Down DEGs
DEG.sig <- readRDS('../Input/toxo_cdc/rds/Intra_DEGs_up_and_down_sig.rds') ## both up & down


## AP2s
AP2s <- read.xlsx('../Input/compScBdTgPb/genes/TF_Info_Updated_kz.xlsx')
AP2s <- left_join(AP2s, TGGT1_ME49, by = c('GeneID' = 'TGGT1'))
AP2s$TGME49 <- gsub('_', '-', AP2s$TGME49)

AP2s <- AP2s[grep("AP2", AP2s$Ap2Name),]

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

sc.atac.dtw.wide <- sc.atac.spline.fits %>% 
  pivot_wider(names_from = 'GeneID', values_from = 'y') %>% 
  mutate_at(vars(matches("TGME")), ~scale(., center = T, scale = T)) %>%
  as.data.frame()


sc.rna.mu.scale <- sc.rna.dtw.wide %>% 
  pivot_longer(-x, names_to = 'GeneID', values_to = 'expr')

sc.atac.mu.scale <- sc.atac.dtw.wide %>% 
  pivot_longer(-x, names_to = 'GeneID', values_to = 'expr')



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
