library(tidyverse)
library(bedtoolsr)
library(openxlsx)
library(grid)
library(matrixStats)
library(tidyverse)
library(tidytext)
library(RColorBrewer)
library(parallel)
library(ComplexHeatmap)
library(circlize)
library(doParallel)
library(edgeR)



#######################################################
############# Peak Gene Assignment ####################
############ Raw counts per region ####################
#######################################################

prod.desc <- read.xlsx('../Input/toxo_genomics/genes/ProductDescription_GT1.xlsx')
gtf.file <- "../Input/toxo_genomics/genome/ToxoDB-59_TgondiiGT1.gtf"
gtf <- read.table(gtf.file, header = F, sep = '\t', quote = NULL)
gtf <- gtf %>% dplyr::filter(grepl('TGGT1', V1))
## using qval = 0.01 and gene size 69350000 bp in macs2 calling 
## narrow peaks called by macs2
narrow.peaks.file <-  "../Input/cut_and_run_toxo/bam/macs2/AP2XII-8_120_150_peaks.narrowPeak"
narrow.peaks <- read.table(narrow.peaks.file, header = F, sep = '\t', quote = NULL)
narrow.peaks <- narrow.peaks %>% dplyr::filter(grepl('TGGT1', V1))

#Unique ID for peaks
narrow.peaks$V4 <- paste(narrow.peaks$V1 , paste(narrow.peaks$V2, narrow.peaks$V3, sep= "-"), sep = ":")

peaks.all.sort <- narrow.peaks %>% arrange(V1, V2, V3)

## Get the transcripts
gtf.trans <- gtf %>% dplyr::filter(V3 == 'mRNA')
parse.str <- strsplit(gtf.trans$V9, split = ';')
inds <- unlist(lapply(parse.str , function(x) which(grepl("geneID", x))))
gtf.trans$V9 <- gsub("geneID=", "", unlist(lapply(1:length(inds), function(i) parse.str[[i]][inds[[i]]])))


## overlap the genes with peaks
genic.peaks <- bedtoolsr::bt.intersect(a = peaks.all.sort, b = gtf.trans, wo = T)


## remaining peaks are intergenic
intergenic.peaks <- peaks.all.sort %>% dplyr::filter(!(V4 %in% genic.peaks$V4))
intergenic.peaks <- intergenic.peaks %>% arrange(V1, V2, V3) 

## Exons
gtf.exon <- gtf %>% dplyr::filter(V3 == 'exon')
gtf.exon.sort <- gtf.exon %>% arrange(V1, V4, V5)
gtf.exon$V9 <- gsub("Parent=", "", gtf.exon$V9)
gtf.exon$V9 <- gsub("-t.*", "", gtf.exon$V9)
gtf.exon <- gtf.exon %>% group_by(V9) %>% mutate(exon.ord = ifelse(V7 == '+', 1:n(), seq(n(), 1, by = -1)),
                                                 multiple.exon = ifelse(n() > 1, T, F))

## Filter for first exon coordinates (exon1 coordinates)
tmp.neg <- gtf.exon %>% filter(V7 == "-") %>% group_by(V9) %>%  dplyr::slice(which.max(V5))
tmp.pos <- gtf.exon %>% filter(V7 == "+") %>% group_by(V9) %>%  dplyr::slice(which.min(V5))
gtf.exon1 <- bind_rows(tmp.pos, tmp.neg)
gtf.exon1.sort <- gtf.exon1 %>% arrange(V1, V4, V5)

## Assign the peaks to nearest upstream gene (look at 5 closest in case of bi-directional)

intergenic.peaks.genes.dist <- bedtoolsr::bt.closest(a = intergenic.peaks, b = gtf.exon1.sort, D = "b", k = 2)

## bedtools pots a 1 at the and of gene name!
intergenic.peaks.genes.dist$V19 <- unlist(lapply(strsplit(intergenic.peaks.genes.dist$V19, split = ''), function(x){
  n = length(x)
  paste(x[1:(n-1)], collapse = '')
}))

## V16 <= 0 means the peak is at upstream 
## Find closest gene among top 5 that is upstreaam (min V16)
intergenic.peaks.genes.dist.upstream <- intergenic.peaks.genes.dist %>% dplyr::filter(V21 <= 0)
intergenic.peaks.genes.dist.upstream <- intergenic.peaks.genes.dist.upstream %>% group_by(V4) %>% 
  mutate(V22 = V21[which.min(abs(V21))])


## Filter the rest
intergenic.peaks.genes.dist.upstream <- intergenic.peaks.genes.dist.upstream %>% dplyr::filter(V21 == V22)

## filter the ones that are too far (1000 bp)
intergenic.peaks.genes.dist.upstream <- intergenic.peaks.genes.dist.upstream %>% dplyr::filter(abs(V21) < 1000)

intergenic.narrow.peaks <- narrow.peaks %>% dplyr::filter(V4 %in% intergenic.peaks.genes.dist.upstream$V4)
write.table(intergenic.narrow.peaks,'../Input/cut_and_run_toxo/bam/macs2/AP2XII-8_120_150_intergenic_peaks.narrowPeak',
            quote = F, row.names = F, col.names = F, sep = '\t')

intergenic.peaks <- intergenic.peaks.genes.dist.upstream %>% ungroup() %>%
  transmute(Chr = V1, strt = V2, stp = V3, peak_id = V4, GeneID = V19)

intergenic.peaks.bed <- intergenic.peaks
intergenic.peaks.bed$strand <- gtf.trans$V7[match(intergenic.peaks.bed$GeneID, gtf.trans$V9)]

write.table(intergenic.peaks.bed,'../Input/cut_and_run_toxo/bam/macs2/AP2XII-8_120_150_intergenic_peaks.bed',
            quote = F, row.names = F, col.names = F, sep = '\t')


intergenic.peaks <- left_join(intergenic.peaks.bed, prod.desc, by = 'GeneID')

write.xlsx(intergenic.peaks, '../Input/cut_and_run_toxo/bam/macs2/AP2XII-8_120_150_intergenic_peaks_gene.xlsx')


