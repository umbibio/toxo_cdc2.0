library(RColorBrewer)
library(Seurat)
library(tidyverse)
library(openxlsx)



#library(sctransform)

source('./util_funcs.R')

num.cores <- detectCores(all.tests = FALSE, logical = TRUE)

## IDs
prod.desc  <- read.xlsx('../Input/toxo_genomics/genes/ProductDescription_GT1.xlsx')
TGGT1_ME49 <- read.xlsx('../Input/toxo_genomics/Orthologs/TGGT1_ME49 Orthologs.xlsx')
prod.desc <- left_join(prod.desc, TGGT1_ME49, by = c('GeneID' = 'TGGT1'))

## AP2s
AP2s <- read.xlsx('../Input/compScBdTgPb/genes/TF_Info_Updated_kz.xlsx')
AP2s <- left_join(AP2s, TGGT1_ME49, by = c('GeneID' = 'TGGT1'))
AP2s$TGME49 <- gsub('_', '-', AP2s$TGME49)

AP2s <- AP2s[grep("AP2", AP2s$Ap2Name),]


atac_sub <- readRDS('../Input/toxo_cdc/rds/S.O_intra_atac_lables_pt.rds')
rna_sub <- readRDS('../Input/toxo_cdc/rds/S.O_intra_lables_pt.rds')

# ## Differential gene expression
Idents(rna_sub) <- 'phase'
Intra.markers <- FindAllMarkers(object = rna_sub, only.pos = T, min.pct = 0)

Intra.markers$GeneID <- gsub('-', '_', Intra.markers$gene)
Intra.markers.top <- Intra.markers %>% group_by(cluster) %>% top_n(2, avg_log2FC)
FeaturePlot(object = rna_sub, 
            features = Intra.markers.top$gene, 
            cols = c("grey", "blue"), reduction = "pca")

Intra.markers.sig <- Intra.markers %>% dplyr::filter(avg_log2FC > log2(1.5) & p_val_adj < 0.05)

saveRDS(Intra.markers.sig, '../Input/toxo_cdc/rds/Intra_markers_sig.rds')

write.xlsx(Intra.markers.sig, '../Output/toxo_cdc/tabels/Intra_markers_sig.xlsx')
ss <- Intra.markers.sig %>% group_by(cluster) %>% summarise(num.DEG = n())
ss$cluster <- factor(ss$cluster, levels = c('G1.a', 'G1.b', 'S', 'M', 'C'))
p <- ggplot(data=ss, aes(x=cluster, y=num.DEG)) +
  geom_bar(stat="identity", fill="steelblue")+
  geom_text(aes(label=num.DEG), vjust=1.6, color="black", size=4, fontface="bold")+
  theme_minimal() + 
  ylab('DEGs') + xlab('') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  theme(
    plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
    axis.title.x = element_text(size=14, face="bold", hjust = 1),
    axis.title.y = element_text(size=14, face="bold")
  ) 

plot(p)

ggsave(filename="../Output/toxo_cdc/figures/tg_Intra_deg_numbers.pdf",
       plot=p,
       width = 8, height = 6,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)

print(ss)



## Differential accessibility analysis: ATAC
Idents(atac_sub) <- 'phase'
ATAC.markers <- FindAllMarkers(object = atac_sub, only.pos = TRUE, min.pct = 0)

ATAC.markers$GeneID <- gsub('-', '_', ATAC.markers$gene)
ATAC.markers.top <- ATAC.markers %>% group_by(cluster) %>% top_n(2, avg_log2FC)
FeaturePlot(object = atac_sub, 
            features = ATAC.markers.top$gene, 
            cols = c("grey", "blue"), reduction = "pca")

ATAC.markers.sig <- ATAC.markers %>% dplyr::filter(avg_log2FC > log2(1.5) & p_val_adj < 0.05)

saveRDS(ATAC.markers.sig, '../Input/toxo_cdc/rds/ATAC_markers_sig.rds')

ss <- ATAC.markers.sig %>% group_by(cluster) %>% summarise(num.DEG = n())
ss$cluster <- factor(ss$cluster, levels = c('G1.a', 'G1.b', 'S', 'M', 'C'))
p <- ggplot(data=ss, aes(x=cluster, y=num.DEG)) +
  geom_bar(stat="identity", fill="steelblue")+
  geom_text(aes(label=num.DEG), vjust=1.6, color="black", size=4, fontface="bold")+
  theme_minimal() + 
  ylab('DEGs') + xlab('') +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(axis.text.y = element_text(angle = 0, hjust = 1, size = 12, face="bold")) +
  theme(strip.background = element_rect(colour="black", fill="white",
                                        size=0.5, linetype="solid")) +
  theme(
    plot.title = element_text(size=14, face = "bold.italic", color = 'red'),
    axis.title.x = element_text(size=14, face="bold", hjust = 1),
    axis.title.y = element_text(size=14, face="bold")
  ) 

plot(p)

ggsave(filename="../Output/toxo_cdc/figures/tg_atac_deg_numbers.pdf",
       plot=p,
       width = 8, height = 6,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)

print(ss)


## DEGs (both Up & Down to identify depleated genes)
Idents(rna_sub) <- 'phase'
DEGs <- FindAllMarkers(object = rna_sub, only.pos = F, min.pct = 0)

DEGs$GeneID <- gsub('-', '_', DEGs$gene)
DEGs.sig <- DEGs %>% dplyr::filter(abs(avg_log2FC) > 0.5 & p_val_adj < 0.05)

DEGs.sig <- left_join(DEGs.sig, prod.desc, by = c('GeneID' = 'TGME49'))
saveRDS(DEGs.sig, '../Input/toxo_cdc/rds/Intra_DEGs_up_and_down_sig.rds')
write.xlsx(DEGs.sig, '../Output/toxo_cdc/tabels/Intra_DEGs_up_and_down_sig.xlsx')


### ATAC region-gene assignment phases

## Identify phase of the atach regions assigned to genes
Intra.markers.sig <- readRDS('../Input/toxo_cdc/rds/Intra_markers_sig.rds')
region_gene_assignment <- readRDS('../Input/toxo_cdc/rds/sc_atac_regions_gene_assignment.rds')
cell.cycle.markers <- Intra.markers.sig %>% transmute(GeneID = GeneID, phase = cluster)
region_cell_clycle_gene_assignment <- inner_join(region_gene_assignment, cell.cycle.markers, by = 'GeneID')

## IDs
prod.desc  <- read.xlsx('../Input/toxo_genomics/genes/ProductDescription_GT1.xlsx')
TGGT1_ME49 <- read.xlsx('../Input/toxo_genomics/Orthologs/TGGT1_ME49 Orthologs.xlsx')
TGGT1_ME49 <- inner_join(TGGT1_ME49, prod.desc, by = c('TGGT1' = 'GeneID'))
region_cell_clycle_gene_assignment <- left_join(region_cell_clycle_gene_assignment, TGGT1_ME49, by = c('GeneID' = 'TGME49'))
saveRDS(region_cell_clycle_gene_assignment, '../Input/toxo_cdc/rds/sc_atac_regions_cell_cycle_gene_assignment.rds')

write.table(region_cell_clycle_gene_assignment, 
            '../Input/toxo_cdc/tables/sc_atac_regions_cell_cycle_gene_assignment.txt', quote = F, 
            sep = '\t', row.names = F, col.names = F)

region_cell_clycle_gene_assignment_G1.a <- region_cell_clycle_gene_assignment %>% 
  dplyr::filter(phase == 'G1.a')
write.table(region_cell_clycle_gene_assignment_G1.a, 
            '../Input/toxo_cdc/tables/sc_atac_regions_cell_cycle_gene_assignment_G1a.txt', quote = F, 
            sep = '\t', row.names = F, col.names = F)


region_cell_clycle_gene_assignment_G1.b <- region_cell_clycle_gene_assignment %>% 
  dplyr::filter(phase == 'G1.b')
write.table(region_cell_clycle_gene_assignment_G1.b, 
            '../Input/toxo_cdc/tables/sc_atac_regions_cell_cycle_gene_assignment_G1b.txt', quote = F, 
            sep = '\t', row.names = F, col.names = F)

region_cell_clycle_gene_assignment_S <- region_cell_clycle_gene_assignment %>% 
  dplyr::filter(phase == 'S')
write.table(region_cell_clycle_gene_assignment_S, 
            '../Input/toxo_cdc/tables/sc_atac_regions_cell_cycle_gene_assignment_S.txt', quote = F, 
            sep = '\t', row.names = F, col.names = F)

region_cell_clycle_gene_assignment_M <- region_cell_clycle_gene_assignment %>% 
  dplyr::filter(phase == 'M')
write.table(region_cell_clycle_gene_assignment_M, 
            '../Input/toxo_cdc/tables/sc_atac_regions_cell_cycle_gene_assignment_M.txt', quote = F, 
            sep = '\t', row.names = F, col.names = F)

region_cell_clycle_gene_assignment_C <- region_cell_clycle_gene_assignment %>% 
  dplyr::filter(phase == 'C')
write.table(region_cell_clycle_gene_assignment_C, 
            '../Input/toxo_cdc/tables/sc_atac_regions_cell_cycle_gene_assignment_C.txt', quote = F, 
            sep = '\t', row.names = F, col.names = F)

