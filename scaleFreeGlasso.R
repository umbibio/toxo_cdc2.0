library(openxlsx)
library(ggplot2)
library(ggfortify)
library(jcolors)
require(gridExtra)
library(grid)
library(matrixStats)
library(tidyverse)
library(RColorBrewer)
library(sctransform)
library(glassoFast)
library(igraph)
library(ggraph)
library(graphlayouts)
library(fdrtool)
library(parallel)

num.cores <- detectCores()

getNetwork <- function(gl, genes, prod.desc){
  
  par.cor <- -cov2cor(gl$wi)
  tmp <- par.cor
  diag(tmp) <- 0
  tmp <- c(tmp)
  ind <- which(tmp != 0)
  if(length(ind) > 0){
    cutoff.cor <- quantile(abs(tmp[ind]), probs = 0.9)
  }else{
    cutoff.cor <- 0
  }
  #cutoff.cor <- par.cor %>% as.data.frame() %>% summarise(across(everything(), ~ quantile(abs(.x), prob = 0.25))) %>% 
  #  quantile(prob = 0.75)
  
  
  recovered.network <- par.cor
  recovered.network[which(abs(par.cor) <=  as.numeric(cutoff.cor))] <- 0
  #recovered.network[which(abs(par.cor) >  as.numeric(cutoff.cor))] <- 1
  diag(recovered.network) <- 0
  recovered.network[upper.tri(recovered.network)] <- 0
  colnames(recovered.network) <- genes
  
  recovered.network.adj <- recovered.network
  recovered.network.adj[which(abs(par.cor) >  as.numeric(cutoff.cor))] <- 1
  
  node.con <- colSums(recovered.network.adj)
  node.con <- data.frame(id = names(node.con), size = as.numeric(node.con))
  
  edge.list <- recovered.network %>% as.data.frame() %>% mutate(from = genes) %>% 
    pivot_longer(-from, names_to = 'to', values_to = 'weight') %>% dplyr::filter(weight != 0) %>%
    transmute(from = from, to = to, weight = abs(weight), sgn = sign(weight))
  
  
  nodes <- data.frame(GeneID = gsub("-", "_", genes), id = genes) %>% left_join(prod.desc, by = 'GeneID')
  
  nodes.list <- nodes %>% dplyr::filter(id %in% unique(c(edge.list$from, edge.list$to))) %>%
    transmute(id = id, type = ProductDescription) %>% left_join(node.con, by = 'id')
  
  nodes.list$class <- 'Protein'
  nodes.list$class[grepl('hypothetical', nodes.list$type)] <- 'hypo'
  nodes.list$class[grepl('conserved', nodes.list$type)] <- 'cons'
  nodes.list$class[grepl('histone', nodes.list$type)] <- 'hist'
  nodes.list$class[grepl('domain', nodes.list$type)] <- 'domain'
  nodes.list$class[grepl('AP2', nodes.list$type)] <- 'AP2'
  
  network <- graph_from_data_frame(d = edge.list, vertices = nodes.list, directed = F)
  
  return(network)
  
}


networksScaleFreeStats <- function(genes, input.dir, prod.desc){
  in.files <- list.files(input.dir)
  
  networks <- mclapply(1:length(in.files), function(i){
    tryCatch(
      expr = {
        gl <- readRDS(paste(input.dir, in.files[i], sep = ''))
        network <- getNetwork(gl, genes, prod.desc)
        d <- degree(network, mode="all")
        fit <- fit_power_law(d+1)
        
        tmp <- rle(sort(d))
        r2 <- summary(lm(log(tmp$lengths)~log(tmp$values)))$adj.r.squared
        
        net.stats <- c(as.numeric(gsub('\\.RData', '', strsplit(in.files[i], split = '_')[[1]][2])),
                       gsize(network),
                       fit$xmin, fit$KS.p, r2)
        return(list(network = network, net.stats = net.stats))
      },
      error = function(v){
        message(paste('error:', v))
      },
      warning = function(w){
        message(paste('error:', w))
      }
    )  }, mc.cores = num.cores)
  
  
  err.ind <- which(unlist(mclapply(networks, function(x) is.null(x$net.stats),  mc.cores = num.cores)))
  
  if(length(err.ind) > 0){
    file.info <- as.data.frame(do.call('rbind', mclapply(networks[-err.ind], `[[`,  2, mc.cores = num.cores)))
  }else{
    file.info <- as.data.frame(do.call('rbind', mclapply(networks, `[[`,  2, mc.cores = num.cores)))
  }
  
  colnames(file.info) <- c('rho', 'size', 'xmin', 'pval', 'r2')
  
  return(list(file.info = file.info, networks = networks))
}

genOptNetwork <- function(genes, file.info, networks, prod.desc){
  file.info.filt <- file.info %>% dplyr::filter(!is.na(r2))
  ind.high.r2 <- which(file.info.filt$r2 > 0.7)
  edge.high.r2 <- mclapply(ind.high.r2, function(i){
    tmp <- as_data_frame(simplify(networks[[i]]$network))[,1:2]
    colnames(tmp) = c('src', 'trg')
    tmp.reverse <- data.frame(src = tmp[,2], trg = tmp[,1])
    tmp <- rbind(tmp, tmp.reverse)
    edge.hash <- paste(tmp$src, tmp$trg, sep = '_')
    return(edge.hash)
  }, mc.cores = num.cores)
  
  edge.rle <- rle(sort(unlist(edge.high.r2)))
  const.edge <- edge.rle$values[edge.rle$lengths >= quantile(edge.rle$lengths, probs = 0.10)]
  
  split.edges <- strsplit(const.edge, split = '_')
  cons.edge.list <- data.frame(from = unlist(lapply(split.edges, `[[`, 1)), to = unlist(lapply(split.edges, `[[`, 2)))
  
  nodes <- data.frame(GeneID = gsub("-", "_", genes), id = genes) %>% left_join(prod.desc, by = 'GeneID')
  nodes.list <- nodes %>% dplyr::filter(id %in% unique(c(cons.edge.list$from, cons.edge.list$to))) %>%
    transmute(id = id, type = ProductDescription)
  
  nodes.list$class <- 'Protein'
  nodes.list$class[grepl('hypothetical', nodes.list$type)] <- 'hypo'
  nodes.list$class[grepl('conserved', nodes.list$type)] <- 'cons'
  nodes.list$class[grepl('histone', nodes.list$type)] <- 'hist'
  nodes.list$class[grepl('domain', nodes.list$type)] <- 'domain'
  nodes.list$class[grepl('AP2', nodes.list$type)] <- 'AP2'
  
  network.opt <- graph_from_data_frame(d = cons.edge.list, vertices = nodes.list, directed = F)
  network.opt <- set_vertex_attr(network.opt, 'size', index = V(network.opt), degree(network.opt))
  network.opt <- simplify(network.opt)
  
  return(network.opt)
}

getClustSubNet <- function(network.opt, phase, genes.clust){
  all.edges <- as_data_frame(simplify(network.opt))[,1:2]
  clust.ind <- which(all.edges[,1] %in% genes.clust$GeneID[genes.clust$phase %in% phase] & 
                       all.edges[,2] %in% genes.clust$GeneID[genes.clust$phase %in% phase])
  
  nodes <- data.frame(GeneID = unique(c(all.edges[clust.ind, 1], all.edges[clust.ind, 2])))
  nodes <- left_join(nodes, genes.clust, by = 'GeneID')
  
  clust.edges <- data.frame(all.edges[clust.ind, ])
  
  g.nds <- graph_from_data_frame(d = clust.edges, vertices = nodes, directed = F)
  node.size <- degree(g.nds, v = V(g.nds))
  
  g.nds <- set_vertex_attr(g.nds, 'size', index = V(g.nds), node.size)
  g.nds <- set_vertex_attr(g.nds, 'clust', index = V(g.nds), nodes$phase)
  
  return(g.nds)
  
}

rankedNetDeg <- function(net, prod.desc, genes.clust){
  
  opt.d <- degree(net)
  opt.d <- data.frame(GeneID = names(opt.d), degree = opt.d)
  opt.d$GeneID <- gsub('-', '_', opt.d$GeneID)
  opt.d <- left_join(opt.d, prod.desc, by = 'GeneID') %>% arrange(desc(degree))
  tmp <- genes.clust
  tmp$GeneID <- gsub('-', '_', tmp$GeneID)
  opt.d <- left_join(opt.d, tmp, by = 'GeneID')
  
  return(opt.d)
  
}

getEdgeList <- function(clust.net, prod.desc){
  e.l <- as_data_frame(simplify(clust.net))[,1:2]
  
  colnames(e.l) <- c('src', 'trg')
  e.l.reverse <- data.frame(src = e.l[,2], trg = e.l[,1])
  e.l <- rbind(e.l, e.l.reverse)
  e.l$src <- gsub('-', '_', e.l$src)
  e.l$trg <- gsub('-', '_', e.l$trg)
  e.l <- left_join(e.l, prod.desc, by = c('src' = 'GeneID'))
  e.l <- left_join(e.l, prod.desc, by = c('trg' = 'GeneID'))
  colnames(e.l) <- c('src', 'trg', 'srcDesc', 'trgDesc')
  
  return(e.l)
}

getNetFromVert <- function(g.obj, nds){
  nds <- gsub('_', '-', nds)
  all.edges <- as_data_frame(simplify(g.obj))[,1:2]
  ind <- which(all.edges[,1] %in% nds | 
                       all.edges[,2] %in% nds)
  
  nodes <- data.frame(GeneID = unique(c(all.edges[ind, 1], all.edges[ind, 2])))

  node.phase <- data.frame(GeneID = V(g.obj)$name, phase = vertex_attr(g.obj, 'phase', index = V(g.obj)))
  nodes <- left_join(nodes, node.phase, by = 'GeneID')
  edges <- data.frame(all.edges[ind, ])
  
  g.nds <- graph_from_data_frame(d = edges, vertices = nodes, directed = F)
  node.size <- degree(g.nds, v = V(g.nds))
  
  g.nds <- set_vertex_attr(g.nds, 'size', index = V(g.nds), node.size)

  #node.phase <- vertex_attr(g.obj, 'phase', index = V(g.obj))
  #g.nds <- set_vertex_attr(g.nds, 'phase', index = V(g.nds), node.phase[match(V(g.nds), V(g.obj))])
  
  return(g.nds)
}

highlightNetwork <- function(g.obj, nds){
  nds <- gsub('_', '-', nds)

  g.obj <- set_edge_attr(g.obj, 'subnet', E(g.obj)[from(nds)], nds)
  nodes <- data.frame(GeneID = unique(c(all.edges[ind, 1], all.edges[ind, 2])))
  
  edges <- data.frame(all.edges[ind, ])
  
  g.nds <- graph_from_data_frame(d = edges, vertices = nodes, directed = F)
  node.size <- degree(g.nds, v = V(g.nds))
  
  g.nds <- set_vertex_attr(g.nds, 'size', index = V(g.nds), node.size)
  
  node.phase <- vertex_attr(g.obj, 'phase', index = V(g.obj))
  g.nds <- set_vertex_attr(g.nds, 'phase', index = V(g.nds), node.phase[match(V(g.nds), V(g.obj))])
  
  return(g.nds)
}


tg.prod.desc <- read.xlsx('../Input/compScBdTgPb/genes/ProductDescription_GT1.xlsx')
#tg.in.Cov <- readRDS('../Input/compScBdTgPb/glasso/S_tg_smooth_gam.RData')
tg.in.Cov <- readRDS('../Input/compScBdTgPb/glasso/S_tg_smooth_gam_de_trended.RData')
#tg.input.dir <- "../Input/compScBdTgPb/glasso/Tg_out/"
tg.input.dir <- "../Input/compScBdTgPb/glasso/Tg_de_trended_out/"
tg.sc.tc.mus.scale <- readRDS('../Input/compScBdTgPb/RData/tg_sc.tc.mus.scale.RData')
tg.genes <- rownames(tg.in.Cov)

tg.L <- networksScaleFreeStats(tg.genes, tg.input.dir, tg.prod.desc)

tg.opt.network <- genOptNetwork(tg.genes, tg.L$file.info, tg.L$networks, tg.prod.desc)
saveRDS(tg.opt.network, '../Input/compScBdTgPb/RData/tg_opt_net.RData')

tg.genes.clust <- tg.sc.tc.mus.scale %>% dplyr::select(GeneID, cluster, phase) %>% distinct()
tg.genes.clust$GeneID <- gsub('_', '-', tg.genes.clust$GeneID)

tg.clust.net <- getClustSubNet(tg.opt.network, c('G1', 'S', 'M', 'C'), tg.genes.clust)
saveRDS(tg.clust.net, '../Input/compScBdTgPb/RData/tg_clust_net.RData')

# define a custom color palette
got_palette <- c("#1A5878", "#C44237", "#AD8941", "#E99093", "#50594B")

node.size <- degree(tg.clust.net, v = V(tg.clust.net))
g.nds <- set_vertex_attr(tg.clust.net, 'size', index = V(tg.clust.net), node.size)



g.nds.filt <- induced_subgraph(g.nds, V(g.nds)[node.size > 1])

p <- ggraph(g.nds.filt,layout = "stress")+
  #geom_edge_link0(aes(edge_width = weight),edge_colour = "grey66")+
  geom_edge_link0(edge_colour = "grey66")+
  geom_node_point(aes(fill = phase, size = size),shape=21)+
  geom_node_text(aes(filter = size > 20, label = name),family="serif")+
  scale_fill_manual(values = got_palette)+
  scale_edge_width(range = c(0.2,3))+
  scale_size(range = c(1,6))+
  theme_graph(base_family = 'Helvetica') +
  theme(legend.position = "right")

plot(p)

ggsave(filename="../Output/compScBdTgPb/figs/tg_cyclic_genes_network.pdf",
       plot=p,
       width = 10, height = 8,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)


tg.clust.net.d <- rankedNetDeg(tg.clust.net, tg.prod.desc, tg.genes.clust)
saveRDS(tg.clust.net.d, '../Input/compScBdTgPb/RData/tg_clust_net_d.RData')
write.xlsx(tg.clust.net.d, '../Output/compScBdTgPb/tables/tg_genes_cyclic_network_nodes.xlsx')


tg.opt.net.d <- rankedNetDeg(tg.opt.network, tg.prod.desc, tg.genes.clust)
saveRDS(tg.opt.net.d, '../Input/compScBdTgPb/RData/tg_opt_net_d.RData')
write.xlsx(tg.opt.net.d, '../Output/compScBdTgPb/tables/tg_genes_opt_network_nodes.xlsx')


tg.e.l <- getEdgeList(tg.clust.net, tg.prod.desc)
saveRDS(tg.e.l, '../Input/compScBdTgPb/RData/tg_e_l.RData')
write.xlsx(tg.e.l, '../Output/compScBdTgPb/tables/tg_cyclic_network_edge_list.xlsx')


tg.opt.e.l <- getEdgeList(tg.opt.network, tg.prod.desc)
saveRDS(tg.opt.e.l, '../Input/compScBdTgPb/RData/tg_opt_e_l.RData')
write.xlsx(tg.opt.e.l, '../Output/compScBdTgPb/tables/tg_opt_network_edge_list.xlsx')



### Interactome of selected gens

## B.d egress Markers
bd.sc.egress.markers.tc.mus.long <- readRDS('../Input/compScBdTgPb/RData/bd_sc_egress._markers_tc_mus_long.RData')
egress.clust1 <- bd.sc.egress.markers.tc.mus.long %>% dplyr::select(GeneID, cluster) %>% dplyr::filter(cluster == 1) %>%
  distinct()

egress.clust2 <- bd.sc.egress.markers.tc.mus.long %>% dplyr::select(GeneID, cluster) %>% dplyr::filter(cluster == 2) %>%
  distinct()

egress.clust3 <- bd.sc.egress.markers.tc.mus.long %>% dplyr::select(GeneID, cluster) %>% dplyr::filter(cluster == 3) %>%
  distinct()

egress.clust4 <- bd.sc.egress.markers.tc.mus.long %>% dplyr::select(GeneID, cluster) %>% dplyr::filter(cluster == 4) %>%
  distinct()

egress.clust1.net <- getNetFromVert(bd.clust.net, egress.clust1$GeneID)
egress.clust2.net <- getNetFromVert(bd.clust.net, egress.clust2$GeneID)
egress.clust3.net <- getNetFromVert(bd.clust.net, egress.clust3$GeneID)
egress.clust4.net <- getNetFromVert(bd.clust.net, egress.clust4$GeneID)

egress.clust1.edge.list <- getEdgeList(egress.clust1.net, bd.prod.desc)
egress.clust2.edge.list <- getEdgeList(egress.clust2.net, bd.prod.desc)
egress.clust3.edge.list <- getEdgeList(egress.clust3.net, bd.prod.desc)
egress.clust4.edge.list <- getEdgeList(egress.clust4.net, bd.prod.desc)

write.xlsx(egress.clust1.edge.list, '../Output/compScBdTgPb/tables/bd_egress_clust1_edges.xlsx')
write.xlsx(egress.clust2.edge.list, '../Output/compScBdTgPb/tables/bd_egress_clust2_edges.xlsx')
write.xlsx(egress.clust3.edge.list, '../Output/compScBdTgPb/tables/bd_egress_clust3_edges.xlsx')
write.xlsx(egress.clust4.edge.list, '../Output/compScBdTgPb/tables/bd_egress_clust4_edges.xlsx')

p4 <- ggraph(egress.clust4.net,layout = "stress") +
  #geom_edge_link0(aes(edge_width = weight),edge_colour = "grey66")+
  geom_edge_link0(edge_colour = "grey66")+
  geom_node_point(aes(fill = phase, size = 3),shape=21)+
  geom_node_text(aes(filter = size > 1, label = name),family="serif")+
  scale_fill_manual(values = got_palette)+
  scale_edge_width(range = c(0.2,3))+
  scale_size(range = c(1,6))+
  theme_graph()+
  theme_graph(base_family = 'Helvetica') +
  theme(legend.position = "right")

plot(p4)

ggsave(filename="../Output/compScBdTgPb/figs/bd_egress_clust4_network.pdf",
       plot=p4,
       width = 6, height = 6,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)




### T. gondii AP2 network

AP2s <- tg.prod.desc$GeneID[grep('AP2', tg.prod.desc$ProductDescription)]

AP2.net <- getNetFromVert(tg.clust.net, AP2s)
AP2.edge.list <- getEdgeList(AP2.net, tg.prod.desc)

write.xlsx(AP2.edge.list, '../Output/compScBdTgPb/tables/tg_AP2_edges.xlsx')

p5 <- ggraph(AP2.net,layout = "stress") +
  #geom_edge_link0(aes(edge_width = weight),edge_colour = "grey66")+
  geom_edge_link0(edge_colour = "grey66")+
  geom_node_point(aes(fill = phase, size = 3),shape=21)+
  geom_node_text(aes(filter = size > 1, label = name),family="serif")+
  scale_fill_manual(values = got_palette)+
  scale_edge_width(range = c(0.2,3))+
  scale_size(range = c(1,6))+
  theme_graph()+
  theme_graph(base_family = 'Helvetica') +
  theme(legend.position = "right")

plot(p5)

ggsave(filename="../Output/compScBdTgPb/figs/bd_egress_clust4_network.pdf",
       plot=p4,
       width = 6, height = 6,
       units = "in", # other options are "in", "cm", "mm"
       dpi = 300
)



## Node networks
nds <- 'TGGT1-214840'
tg.clust.net.nds <- set_edge_attr(tg.clust.net, 'subnet', E(tg.clust.net)[from(nds)], nds)
tg.clust.net.nds <- set_edge_attr(tg.clust.net.nds, 'subnet.weight', E(tg.clust.net.nds)[from(nds)], 2)
tg.clust.net.nds <- set_edge_attr(tg.clust.net.nds, 'subnet.weight', E(tg.clust.net.nds)[!from(nds)], 1)

ggraph(tg.clust.net.nds,layout = "stress") +
  #geom_edge_link(aes(colour = subnet, width = subnet.weight, alpha = subnet.weight)) + 
  geom_edge_link(aes(colour = subnet)) + 
  #geom_edge_link0(aes(colour = subnet))+
  #geom_edge_link0(edge_colour = "grey66")+
  geom_node_point(aes(fill = phase, size = size), shape=21)+
  geom_node_text(aes(filter = name == nds, label = name),family="serif")+
  scale_fill_manual(values = got_palette)+
  scale_edge_width(range = c(0.8,2))+
  scale_size(range = c(1,6))+
  theme_graph()+
  theme(legend.position = "top")


node.net <- getNetFromVert(tg.clust.net, 'TGGT1_271930')

ggraph(node.net,layout = "stress") +
  #geom_edge_link0(aes(edge_width = weight),edge_colour = "grey66")+
  geom_edge_link0(edge_colour = "grey66")+
  geom_node_point(aes(fill = phase, size = 3),shape=21)+
  geom_node_text(aes(filter = size > 1, label = name),family="serif")+
  scale_fill_manual(values = got_palette)+
  scale_edge_width(range = c(0.2,3))+
  scale_size(range = c(1,6))+
  theme_graph()+
  theme(legend.position = "right")

#saveRDS(file.info, '../Input/BdivCellCycle/RDS/file_info.RData')

## Fitting the rho/r2 curve with smooth splines
file.info.filt <- file.info %>% dplyr::filter(!is.na(r2))
spline.fit <- smooth.spline(x = file.info.filt$rho, 
                            y = file.info.filt$r2,lambda = 0.0005)

## Compute the derivatives of the fitted splines
s.0 <- predict(spline.fit, spline.fit$x, deriv=0)
s.1 <- predict(spline.fit, spline.fit$x, deriv=1)
s.derv <- data.frame(s0=s.0$y, s1=s.1$y)

## Get the location of the extrema
locs <- rle(den.sign <- sign(s.derv$s1))
plot(file.info.filt$rho, file.info.filt$r2, pch = 20, xlab = 'rho', ylab = 'r^2')
points(s.0$x, s.0$y, type = 'l', lwd=2,  col = 'blue')
points(file.info.filt$rho[ind.high.r2], file.info.filt$r2[ind.high.r2], pch = 21, col = 'red', )

# inc.ind <- which(locs$values == 1)
# maxima.ind = {}
# for(i in inc.ind){
#   maxima.ind = c(maxima.ind, sum(locs$lengths[1:i]))
# }
# 
# ## Interpolate a point between the location where derivative changes sign
# maxima = (spline.fit$x[maxima.ind] + spline.fit$x[(maxima.ind + 1)]) / 2
# maxima = maxima[!is.na(maxima)]
# ## Get the maximum values
# maxval = predict(spline.fit, maxima)
# 
# ## Get the outliers
# maxima.outliers = which(maxval$y >= quantile(maxval$y, prob = 0.8))
# 
# ## Peaks for entities of interest
# entity.x = maxval$x[maxima.outliers]
# entity.y = maxval$y[maxima.outliers]
# 
# plot(file.info.filt$rho, file.info.filt$r2, pch = 20, col = 'blue')
# points(s.0$x, s.0$y, type = 'l')
# points(entity.x, entity.y, col = 'red', pch = 8)
# 
# network.opt <- networks[[which(file.info$rho >= entity.x[length(entity.x)])[1]]]$network

#### Network similarity/consistency

saveRDS(network.opt, '../Input/BdivCellCycle/RDS/bd_network_opt.RData')


clust.net <- getClustSubNet(network.opt, c('G1b', 'G1c'), tc_clust)

# define a custom color palette
got_palette <- c("#1A5878", "#C44237", "#AD8941", "#E99093", "#50594B")

ggraph(clust.net,layout = "stress")+
  #geom_edge_link0(aes(edge_width = weight),edge_colour = "grey66")+
  geom_edge_link0(edge_colour = "grey66")+
  geom_node_point(aes(fill = clust, size = size),shape=21)+
  geom_node_text(aes(filter = size>40, label = name),family="serif")+
  scale_fill_manual(values = got_palette)+
  scale_edge_width(range = c(0.2,3))+
  scale_size(range = c(1,6))+
  theme_graph()+
  theme(legend.position = "right")


clust.net.d <- rankedNetDeg(clust.net, prod.desc, tc_clusts)



########### OLDE STUFF

## calculate the jaccard distance
# getJD <- function(edge.high.r2, i, j){
#   return(length(intersect(edge.high.r2[[i]], edge.high.r2[[j]]))/length(union(edge.high.r2[[i]], edge.high.r2[[j]])))
# }
# 
# distMat.list <- mclapply(1:(length(edge.high.r2) - 1), function(i){
#   mclapply((i+1):length(edge.high.r2), function(j){
#     return(getJD(edge.high.r2, i, j))
#     }, mc.cores = num.cores - 1)
#   }, mc.cores = 1)
# 
# dist.mat <- do.call('rbind', distMat.list)
# dist.mat <- as.dist(dist.mat[upper.tri(dist.mat)] <- dist.mat[lower.tri(dist.mat)])
# hc_eucledian <- hclust(dist.mat, method = "ward.D" )


####

d <- degree(bd_network.opt, mode="all")
fit <- fit_power_law(d+1)
fit
hist(d, nclass = 200)

tmp <- rle(sort(d))
ff <- lm(log(tmp$lengths)~log(tmp$values))
ff.sum <- summary(ff)

plot(log(tmp$lengths)~log(tmp$values))

abline(b = ff.sum$coefficients[2], a =  ff.sum$coefficients[1], col='red')


nodes.degree <- dplyr::bind_rows(vertex_attr(bd_network.opt, index = V(bd_network.opt)))
#putative.regs <- nodes.degree$name[which(nodes.degree$size > quantile(nodes.degree$size, probs = 0.3))]
putative.regs <- nodes.degree$name
e.l <- as.data.frame(get.edgelist(bd_network.opt))
colnames(e.l) <- c('src', 'trg')
e.l.reverse <- data.frame(src = e.l[,2], trg = e.l[,1])
e.l <- rbind(e.l, e.l.reverse)

e.l.regs <- e.l %>% dplyr::filter(src %in% putative.regs) %>% group_by(src) %>% 
  summarise(targs = list(trg), num.trg = length(trg))

tc_clusts <- read.xlsx('../Output/scClockOut/sc_time_course_gene_clusters.xlsx')
tc_clusts <- tc_clusts %>% dplyr::select(GeneID, cluster, phase) %>% distinct()
tc_clusts$GeneID <- gsub('_', '-', gsub('___.*', '', tc_clusts$GeneID))
clust.nums <- tc_clusts %>% group_by(phase) %>% summarise(clust.genes = list(GeneID), clust.gene.num = n())

clust.nums$dummy = 1
e.l.regs$dummy <- 1
e.l.regs.clust <- left_join(e.l.regs, clust.nums, by = 'dummy')

e.l.regs.clust$all <- length(genes)

e.l.regs.clust <- e.l.regs.clust %>% rowwise() %>% 
  mutate(overlap = length(intersect(c(targs), c(clust.genes))),
         overlap.genes = list(intersect(c(targs), c(clust.genes)))) %>% 
  mutate(pvalue = fisher.test(matrix(c(overlap,  num.trg- overlap, clust.gene.num - overlap, 
                                       all - (num.trg + clust.gene.num - overlap) ), byrow = T, ncol = 2, nrow = 2),
                              alternative = "greater")$p.value)
e.l.regs.clust$qval <- fdrtool(e.l.regs.clust$pvalue, statistic = 'pvalue', plot = F)$qval

e.l.regs.clust$proportion <- e.l.regs.clust$overlap / e.l.regs.clust$clust.gene.num
con.stat <- e.l.regs.clust %>% dplyr::filter(qval < 0.01) %>% group_by(phase) %>% summarise(quant = quantile(proportion, prob = 0.5))
e.l.regs.clust <- left_join(e.l.regs.clust, con.stat, by = 'phase')

## matching with Caroline's list
TFs <- read.xlsx('../Input/BdivCellCycle/regulators/gene_regulators_kz.xlsx', sheet = 1)
colnames(TFs) <- c('Pf', 'Bdiv','annotation','notes')
TFs <- TFs %>% mutate(Bdiv = strsplit(Bdiv, ',')) %>% unnest(Bdiv)
TFs$Bdiv <- as.character(gsub('\\t', '', gsub('.*best', '', gsub(' ', '', TFs$Bdiv ))))
TFs <- TFs %>% dplyr::filter(!grepl('noo', Bdiv))
TFs$classification <- 'TF'
TFs$domain <- 'NA'
TFs$GeneName <- 'NA'


epigenetic.regs <- read.xlsx('../Input/BdivCellCycle/regulators/gene_regulators_kz.xlsx', sheet = 2)
colnames(epigenetic.regs) <- c('classification', 'Pf', 'GeneName', 'annotation', 'domain', 'Bdiv', 'notes')
epigenetic.regs <- epigenetic.regs %>% mutate(Bdiv = strsplit(Bdiv, ',')) %>% unnest(Bdiv)
epigenetic.regs$Bdiv <- as.character(gsub('\\t', '', gsub('.*best', '', gsub(' ', '', epigenetic.regs$Bdiv ))))

tmp1 <- TFs %>% transmute(Bdiv = Bdiv, Pf = Pf, GeneName = GeneName, domain = domain, 
                          classification = classification, annotation = annotation, notes = notes)

tmp2 <- epigenetic.regs %>% transmute(Bdiv = Bdiv, Pf = Pf, GeneName = GeneName, domain = domain,  
                                      classification = classification, annotation = annotation, notes = notes)

known.regs <- rbind(tmp1, tmp2)
known.regs$Bdiv <- gsub('_', '-', known.regs$Bdiv)

e.l.regs.clust <- left_join(e.l.regs.clust, known.regs, by = c('src' = 'Bdiv'))

saveRDS(e.l.regs.clust, '../Input/BdivCellCycle/RDS/e_l_regs_clust.RData')
#putative.clust.regs <- e.l.regs.clust %>% 
#  dplyr::filter(qval < 0.01, proportion >= quant) %>% arrange(desc(overlap))
putative.clust.regs <- e.l.regs.clust %>% 
  dplyr::filter(qval < 0.01) %>% arrange(desc(overlap))

putative.clust.regs$src <- gsub('-', '_', putative.clust.regs$src)
putative.clust.regs <- left_join(putative.clust.regs, prod.desc, by = c('src' = 'GeneID'))

saveRDS(putative.clust.regs, '../Input/BdivCellCycle/RDS/putative_clust_regs.RData')
write.xlsx(putative.clust.regs, '../Output/BdivCellCycle/cluster_putative_regulatrs.xlsx')



cat('top 2 overall regulators\n')
putative.clust.regs %>% ungroup() %>% 
  dplyr::select(src, num.trg, phase, clust.gene.num, overlap, qval, proportion, ProductDescription) %>% 
  group_by(phase) %>% slice_max(order_by = proportion, n = 2)

cat('\ntop 2 know regulators\n')
putative.clust.regs %>% ungroup() %>% dplyr::filter(!is.na(Pf)) %>%
  dplyr::select(src, num.trg, phase, clust.gene.num, overlap, qval, 
                proportion, Pf, GeneName, ProductDescription) %>% 
  group_by(phase) %>% slice_max(order_by = proportion, n = 2)



### Graph regulators



# define a custom color palette
got_palette <- c("#1A5878", "#C44237", "#AD8941", "#E99093", "#50594B")

getNetFromVert <- function(g.obj,nds, reg.info){
  all.edges <- get.edgelist(g.obj)
  nds.edges <- unique(c(which(all.edges[,1] %in% nds), which(all.edges[,2] %in% nds)))
  nodes <- data.frame(GeneID = unique(c(all.edges[nds.edges, ])))
  reg.info$src <- gsub('_', '-', reg.info$src)
  
  nodes <- cbind(nodes, reg.info[match(nodes$GeneID, reg.info$src),])
  
  nodes$reg <- ifelse(nodes$GeneID %in% nds, 'Y', 'N')
  g.nds <- graph_from_data_frame(d = data.frame(all.edges[nds.edges, ]), vertices = nodes, directed = F)
  node.size <- degree(g.nds, v = V(g.nds))
  g.nds <- set_vertex_attr(g.nds, 'size', index = V(g.nds), node.size)
  return(g.nds)
}

putative.clust.regs.known <- putative.clust.regs %>% dplyr::filter(!is.na(Pf))
nds <- gsub('_', '-',unique(putative.clust.regs.known$src))
g.nds <- getNetFromVert(bd_network.opt,nds,putative.clust.regs.known)

ggraph(g.nds,layout = "stress")+
  #geom_edge_link0(aes(edge_width = weight),edge_colour = "grey66")+
  geom_edge_link0(edge_colour = "grey66")+
  geom_node_point(aes(fill = phase, size = size),shape=21)+
  geom_node_text(aes(filter = reg == 'Y', label = name),family="serif")+
  scale_fill_manual(values = got_palette)+
  scale_edge_width(range = c(0.2,3))+
  scale_size(range = c(1,6))+
  theme_graph()+
  theme(legend.position = "right")


## Plot expression curves of regulators
## Plot a few curves to check the alignments

library(sme)
source('util_funcs.R')
sc.tc.df.adj <- readRDS('../Input/scClock/sc.tc.df.adj.RData')
sc.tc.fits <- readRDS("../Input/scClock/sme_fits_sc_tc_20min.RData")
vs = gsub('-', '_', nds)


pdf(file = "../Output/scClockFigs/putative_regulators_sme_fits_sc.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches

par(mfrow = c(4,4))
for(v in vs){
  ind <- which(unique(sc.tc.df.adj$variable) == v)
  plot.sme(sc.tc.fits[[ind]], v)
}

dev.off()



## Regulators
tb <- putative.clust.regs %>% ungroup() %>% dplyr::filter(!is.na(Pf)) %>%
  dplyr::select(src, phase,Pf, GeneName, annotation, ProductDescription) %>% 
  group_by(src) %>% summarise(phase = phase[1], Pf = Pf[1], 
                              GeneName = gsub('NA', '', paste(unique(GeneName), collapse  = ',')), 
                              annotation = gsub('NA', '', paste(unique(annotation), collapse  = ',')), 
                              Product.Description = Product.Description[1])

print(tb, width = Inf)


######################################## OLDER stuff







#######################################
## Examine sub network connectivity
sub.net <- list()
ph <- unique(tc_clusts$phase)
sub.net.info <- data.frame(phase = rep(0, length(ph)), 
                           num.nodes = rep(0, length(ph)), 
                           num.edges = rep(0, length(ph)))
for(i in 1:length(ph)){
  ind <- tc_clusts$GeneID[tc_clusts$phase == ph[i]] %in% names(V(bd_network.opt))
  putative.clust.regs.ph <- putative.clust.regs %>% dplyr::filter(phase == ph[i])
  regs <- putative.clust.regs.ph %>% dplyr::filter(!is.na(Pf))
  g.sub <- induced_subgraph(bd_network.opt, unique(c(tc_clusts$GeneID[tc_clusts$phase == ph[i]][ind], 
                                              gsub('_', '-',unique(regs$src)))))
  reg <- ifelse(names(V(g.sub)) %in% gsub('_', '-',unique(regs$src)), 'Y', 'N')
  g.sub <- set_vertex_attr(g.sub, 'reg', index = V(g.sub), reg)
  sub.net <- c(sub.net, list(g.sub))
  d <- degree(g.sub, mode="all")
  sub.net.info$phase[i] <- ph[i]
  sub.net.info$num.nodes[i] <- length(V(g.sub))
  sub.net.info$num.edges[i] <- gsize(g.sub)
}



# define a custom color palette
got_palette <- c("#1A5878", "#C44237", "#AD8941", "#E99093", "#50594B")

my.graph <- sub.net[[3]]
Isolated = which(degree(my.graph) < 3)
my.graph = delete.vertices(my.graph, Isolated)

ggraph(my.graph,layout = "stress")+
  #geom_edge_link0(aes(edge_width = weight),edge_colour = "grey66")+
  geom_edge_link0(edge_colour = "grey66")+
  geom_node_point(aes(fill = class, size = size),shape=21)+
  geom_node_text(aes(filter = reg == 'Y', label = name),family="serif")+
  scale_fill_manual(values = got_palette)+
  scale_edge_width(range = c(0.2,3))+
  scale_size(range = c(1,6))+
  theme_graph()+
  theme(legend.position = "right")




###########
file.info <- data.frame(rho = rep(0, length(in.files)), size = rep(0, length(in.files)), 
                        xmin = rep(0, length(in.files)), pval = rep(0, length(in.files)))


plot(file.info$rho, file.info$r2)
which(file.info$rho==0.1376)
bd_network <- bd_networks[[which(file.info$rho==0.1376)]]$bd_network

d <- degree(bd_network, mode="all")
fit <- fit_power_law(d+1)
fit
hist(d, nclass = 200)

tmp <- rle(sort(d))
ff <- lm(log(tmp$lengths)~log(tmp$values))

plot(log(tmp$lengths)~log(tmp$values))

abline(b = -1.122650, a =  6.335302)

plot(file.info$rho,file.info$r2)

for(i in 1:length(in.files)){
  cat(paste('processing file', i))
  cat('\n')
  gl <- readRDS(paste(input.dir, in.files[i], sep = ''))
  bd_network <- getNetwork(gl, genes, prod.desc)
  bd_networks <- c(bd_networks, list(bd_network))
  d <- degree(bd_network, mode="all")
  ## Find best fitting xmin
  #qs <- quantile(d+1, probs = c(0.75, 0.99))
  #fits <- lapply(unique(sort(d[d+1 >qs[1] & d+1 <qs[2]] + 1)), function(xmin) fit_power_law(d+1, xmin))
  #xmins <- unlist(lapply(fits, `[[`, 3))
  #pvals <- unlist(lapply(fits, `[[`, 6))
  #best.ind <- which.min(pvals)
  fit <- fit_power_law(d+1)

  file.info$rho[i] <- as.numeric(gsub('\\.RData', '', strsplit(in.files[i], split = '_')[[1]][2]))
  file.info$size[i] <- gsize(bd_network)
  #file.info$xmin[i] <- xmins[best.ind]
  #file.info$pval[i] <- pvals[best.ind]
  file.info$xmin[i] <- fit$xmin
  file.info$pval[i] <- fit$KS.p
  
}

#saveRDS(bd_networks, '../Input/scClock/glasso/bd_networks.RData')
##ff <- fdrtool(file.info$pval,statistic = 'pvalue', color.figure = F)
##file.info$qval <- ff$qval

file.info.filt <- file.info %>% dplyr::filter(size > 0 & size < quantile(file.info$size, probs = 0.95) & pval < 0.05)
plot(file.info.filt$size, file.info.filt$pval)

rho.opt <- file.info.filt$rho[which.min(file.info.filt$pval)]
opt.network <- bd_networks[[which(file.info$rho == rho.opt)]] 

## identify hubs
d <- degree(opt.network, mode="all")
pl <- power.law.fit(d+1)
hist(d)
xx <- as.data.frame(do.call('cbind',vertex_attr(opt.network, index = V(opt.network))))
xx[xx$name %in% names(d[order(d,decreasing=TRUE)][1:20]),] %>% arrange(desc(size))




## Get the gene clusters
tc_clusts <- read.xlsx('../Output/scClockOut/sc_time_course_gene_clusters.xlsx')
tc_clusts <- tc_clusts %>% dplyr::select(GeneID, cluster, phase) %>% distinct()
tc_clusts$GeneID <- gsub('_', '-', gsub('___.*', '', tc_clusts$GeneID))
clust.nums <- tc_clusts %>% group_by(phase) %>% summarise(total.genes = n())

e.l <- as.data.frame(get.edgelist(opt.network))
colnames(e.l) <- c('src', 'trg')

## sub network containing gene clusters
cluster.src.all <- left_join(tc_clusts, e.l, by = c('GeneID' = 'trg')) %>% na.omit() %>%
  group_by(src) %>% summarise(targs = list(GeneID), num.targs = length(GeneID)) %>% ungroup() 
cluster.src.all <- cluster.src.all%>% mutate(proportion = num.targs/nrow(tc_clusts)) 
cluster.src.all <- cluster.src.all %>% arrange(desc(proportion))

hubs <- cluster.src.all$src[cluster.src.all$num.targs > quantile(cluster.src.all$num.targs, probs = 0.9)]

orthologs <- read.xlsx('../Input/orthologs/all_reciprocal_bet_hits.xlsx')
orthologs.Bd.Tg <- read.xlsx('../Input/orthologs/GT1_BDiv.xlsx')

tmp <- cluster.src.all %>% dplyr::filter(num.targs > quantile(cluster.src.all$num.targs, probs = 0.9)) %>%
  dplyr::select(src, num.targs)
tmp$src <- gsub('-', '_', tmp$src)
putative.regs <- inner_join(orthologs.Bd.Tg, tmp, by = c('BDiv' = 'src')) %>% arrange(desc(num.targs))

ind <- tc_clusts$GeneID %in% names(V(opt.network))
g.sub <- induced_subgraph(opt.network, c(hubs, tc_clusts$GeneID[ind]))


# define a custom color palette
got_palette <- c("#1A5878", "#C44237", "#AD8941", "#E99093", "#50594B")

my.graph <- g.sub
Isolated = which(degree(g.sub) <= 3)
my.graph = delete.vertices(my.graph, Isolated)

ggraph(my.graph,layout = "stress")+
  #geom_edge_link0(aes(edge_width = weight),edge_colour = "grey66")+
  geom_edge_link0(edge_colour = "grey66")+
  geom_node_point(aes(fill = class, size = size),shape=21)+
  geom_node_text(aes(filter = size >= 150, label = name),family="serif")+
  scale_fill_manual(values = got_palette)+
  scale_edge_width(range = c(0.2,3))+
  scale_size(range = c(1,6))+
  theme_graph()+
  theme(legend.position = "right")



## per phase sub network containing gene clusters
cluster.src <- left_join(tc_clusts, e.l, by = c('GeneID' = 'trg')) %>% na.omit() %>%
  group_by(src, phase) %>% summarise(targs = list(GeneID), num.targs = length(GeneID)) %>% ungroup() 
cluster.src <- left_join(cluster.src , clust.nums, by = 'phase') %>% mutate(proportion = num.targs/total.genes) 
cluster.src <- cluster.src %>% arrange(desc(proportion))



## Examine sub network connectivity
sub.net <- list()
ph <- unique(tc_clusts$phase)
sub.net.info <- data.frame(phase = rep(0, length(ph)), 
                           num.nodes = rep(0, length(ph)), 
                           num.edges = rep(0, length(ph)))
for(i in 1:length(ph)){
  ind <- tc_clusts$GeneID[tc_clusts$phase == ph[i]] %in% names(V(opt.network))
  g.sub <- induced_subgraph(opt.network, tc_clusts$GeneID[tc_clusts$phase == ph[i]][ind])
  sub.net <- c(sub.net, list(g.sub))
  d <- degree(g.sub, mode="all")
  sub.net.info$phase[i] <- ph[i]
  sub.net.info$num.nodes[i] <- length(V(g.sub))
  sub.net.info$num.edges[i] <- gsize(g.sub)
}

## Get hubs of each modules

regs <- list()
for(i in 1:length(ph)){
  ph.reg <- cluster.src %>% dplyr::filter(phase == ph[i] &  proportion > 0.2) %>% dplyr::select(src, num.targs)
  regs <- c(regs, ph.reg)
}

S.M.reg <- cluster.src %>% dplyr::filter(phase == 'S/M' &  proportion > 0.3) %>% dplyr::select(src, num.targs)
G1a.reg <- cluster.src %>% dplyr::filter(phase == 'G1c' &  proportion > 0.2) %>% dplyr::select(src, num.targs)


# define a custom color palette
got_palette <- c("#1A5878", "#C44237", "#AD8941", "#E99093", "#50594B")

my.graph <- sub.net[[1]]
Isolated = which(degree(my.graph) < 3)
my.graph = delete.vertices(my.graph, Isolated)

ggraph(my.graph,layout = "stress")+
  #geom_edge_link0(aes(edge_width = weight),edge_colour = "grey66")+
  geom_edge_link0(edge_colour = "grey66")+
  geom_node_point(aes(fill = class, size = size),shape=21)+
  geom_node_text(aes(filter = size >= 150, label = name),family="serif")+
  scale_fill_manual(values = got_palette)+
  scale_edge_width(range = c(0.2,3))+
  scale_size(range = c(1,6))+
  theme_graph()+
  theme(legend.position = "right")


xx <- as.data.frame(do.call('cbind',vertex_attr(bd_networks[[1440]], index = V(bd_networks[[1440]]))))
d <- degree(bd_networks[[1440]], mode="all")

xx[xx$name %in% names(d[order(d,decreasing=TRUE)][1:20]),] %>% arrange(desc(size))
pl <- power.law.fit(d+1)

d = degree(bd_network, mode = "all")
hist(d)


## Get the gene clusters
tc_clusts <- read.xlsx('../Output/scClockOut/sc_time_course_gene_clusters.xlsx')

tc_clusts <- tc_clusts %>% dplyr::select(GeneID, cluster, phase) %>% distinct()
tc_clusts$GeneID <- gsub('_', '-', gsub('___.*', '', tc_clusts$GeneID))
clust.nums <- tc_clusts %>% group_by(phase) %>% summarise(total.genes = n())

e.l <- as.data.frame(get.edgelist(bd_networks[[1440]]))
colnames(e.l) <- c('src', 'trg')
cluster.src <- left_join(tc_clusts, e.l, by = c('GeneID' = 'trg')) %>% na.omit() %>%
  group_by(src, phase) %>% summarise(targs = list(GeneID), num.targs = length(GeneID)) %>% ungroup() 
cluster.src <- left_join(cluster.src , clust.nums, by = 'phase') %>% mutate(proportion = num.targs/total.genes) 
cluster.src <- cluster.src %>% arrange(desc(proportion))

S.M.reg <- cluster.src %>% dplyr::filter(phase == 'S/M' &  proportion > 0.3) %>% dplyr::select(src, num.targs)
G1a.reg <- cluster.src %>% dplyr::filter(phase == 'G1c' &  proportion > 0.2) %>% dplyr::select(src, num.targs)

g.sub <- induced_subgraph(bd_networks[[1440]], c(S.M.reg$src, tc_clusts$GeneID[tc_clusts$phase == 'S/M']))
ind <- c(G1a.reg$src, tc_clusts$GeneID[tc_clusts$phase == 'G1c']) %in% names(V(bd_networks[[1440]]))
g.sub <- induced_subgraph(bd_networks[[1440]], c(G1a.reg$src, tc_clusts$GeneID[tc_clusts$phase == 'G1c'])[ind])

got_palette <- c("#1A5878", "#C44237", "#AD8941", "#E99093", "#50594B")
vertex_attr(g.sub, "reg", S.M.reg$src) <- 'reg'
vertex_attr(g.sub, "reg", G1a.reg$src) <- 'reg'

ggraph(g.sub,layout = "stress")+
  #geom_edge_link0(aes(edge_width = weight),edge_colour = "grey66")+
  geom_edge_link0(edge_colour = "grey66")+
  geom_node_point(aes(fill = class,size = size),shape=21)+
  geom_node_text(aes(filter = !is.na(reg), label = name),family="serif")+
  scale_fill_manual(values = got_palette)+
  scale_edge_width(range = c(0.2,3))+
  scale_size(range = c(1,6))+
  theme_graph()+
  theme(legend.position = "right")

orthologs <- read.xlsx('../Input/orthologs/all_reciprocal_bet_hits.xlsx')
orthologs.Bd.Tg <- read.xlsx('../Input/orthologs/GT1_BDiv.xlsx')

orthologs[which(orthologs$Bdiv %in% gsub('-', '_', S.M.reg$src)),]
orthologs[which(orthologs$Bdiv %in% gsub('-', '_', G1a.reg$src)),]

xx <- as.data.frame(do.call('cbind',vertex_attr(g.sub, index = V(g.sub))))
d <- degree(g.sub, mode="all")

xx[xx$name %in% names(d[order(d,decreasing=TRUE)][1:20]),]

