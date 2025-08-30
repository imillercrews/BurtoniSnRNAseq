#### Burtoni snseq seurat analysis
### neuropeptides comparison
###Note: Seurat requires R version > 4
## use lambcomp1 to run R with command 
# > R-4.0.3
###DEsingle
##https://miaozhun.github.io/DEsingle/


### set working directory
setwd("/stor/work/Hofmann/All_projects/A_burtoni_snseq/seurat/")

#### load libraries ####
##install libraries
# install.packages("tidyverse",
#                  repos = "https://cloud.r-project.org")
# install.packages("Seurat",
#                  repos = "https://cloud.r-project.org")
# install.packages("patchwork",
#                  repos = "https://cloud.r-project.org")
#install.packages('BiocManager')
#BiocManager::install('multtest')
# install.packages('metap',
#                  repos = "https://cloud.r-project.org")
# install.packages('clustree')
# if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")
# BiocManager::install("DEsingle")
#load libraries
library(Seurat)
library(patchwork)
library(clustree)
library(pheatmap)
library(DEsingle)
library(igraph)
library(ggridges)
library(forcats)
library(tidyverse)
library(ggrepel)


#### load data ####
### single cell data
# load('burtoni.snseq.combined.sct.all.RData')
load('burtoni.snseq.combined.sct.RData')

### scsorter data
load("burtoni.scsorter.data.scsort.output.RData")

### load souporcell data
burtoni.souporcell.filtered = read.csv('../souporcell/burtoni.souporcell.filtered.csv')

### load sctype.hypo data
burtoni.sctypemarkers.hypo = read.csv('./sctype.hypo/sctypemarkers.hypo.unknown.csv')

### Add metadata
## cell type
burtoni.snseq.combined.sct = AddMetaData(
  object = burtoni.snseq.combined.sct,
  metadata = burtoni.scsorter.data.scsort.output %>% 
    select(Cell.type,
           Cell.id) %>% 
    column_to_rownames(var = "Cell.id"),
  col.name = 'Cell.type'
)

## genotype
burtoni.snseq.combined.sct = AddMetaData(
  object = burtoni.snseq.combined.sct,
  metadata = burtoni.souporcell.filtered %>% 
    select(Genotype.id,
           Cell.id) %>% 
    column_to_rownames(var = "Cell.id"),
  col.name = 'Genotype.id'
)

## sctype.hypo
burtoni.snseq.combined.sct = AddMetaData(
  object = burtoni.snseq.combined.sct,
  metadata = burtoni.sctypemarkers.hypo %>% 
    select(sctypemarkers.hypo,
           Cell.id) %>% 
    column_to_rownames(var = "Cell.id"),
  col.name = 'sctypemarkers.hypo'
)

# resolution range
resolution.range <- seq(from = 0, to = 2, by = 0.2)

#reduced range
resolution.range.reduced <- seq(from = 0, to = 1.2, by = 0.2)

#reduced range 2
resolution.range.reduced.2 <- seq(from = 0.2, to = 1, by = 0.1)

### neuropeptide list
neuropeptides.list = read_csv("../Gene.lists/neuropeptides.list_orthologs.csv")



# #### recluster without vascular ####
# 
# burtoni.snseq.combined.sct.all.subset = burtoni.snseq.combined.sct.all
# 
# #set idents to broad cell class
# Idents(object = burtoni.snseq.combined.sct.all.subset) <- "Cell.class.broad"
# 
# #subset to neurons
# burtoni.snseq.combined.sct.all.subset = subset(burtoni.snseq.combined.sct.all.subset,
#                                                idents = c("neurons",
#                                                           "glia"))
# 
# #set idents to cell type
# Idents(object = burtoni.snseq.combined.sct.all.subset) <- "Cell.type"
# 
# ## run PCA, UMAP, and cluster 
# #use 0.4 resolution
# burtoni.snseq.combined.sct.all.subset.recluster = burtoni.snseq.combined.sct.all.subset %>% 
#   RunPCA() %>%
#   FindNeighbors(dims = 1:15) %>%
#   RunUMAP(dims = 1:15) %>%
#   FindClusters(resolution = 0.8)
# 
# 
# ### clustree
# # cluster across resolutions
# burtoni.snseq.combined.sct.all.subset.clustree <- Seurat::FindClusters(object = burtoni.snseq.combined.sct.all.subset.recluster, 
#                                                                        resolution = resolution.range.reduced)
# 
# #clustree
# clustree(burtoni.snseq.combined.sct.all.subset.clustree, 
#          prefix = "integrated_snn_res.",
#          scale_node_text = TRUE) +
#   scale_edge_color_continuous(low = "black", 
#                               high = "black") +
#   theme(legend.position = "bottom")
# ggsave('cell.subset/cell.subset.clustree.png',
#        width = 10,
#        height = 10)
# 
# 
# ### graph 
# # # idents to new clusters
# # Idents(object = burtoni.snseq.combined.sct.all.subset.recluster) <- "integrated_snn_res.0.8"
# # 
# # ## hover locator
# # DimPlot(burtoni.snseq.combined.sct.all.subset.recluster, 
# #         reduction = "umap", 
# #         label = TRUE,
# #         repel = TRUE, 
# #         group.by = "SCT_snn_res.0.8") %>% 
# #   HoverLocator(information = FetchData(object = burtoni.snseq.combined.sct.all.subset.recluster, 
# #                                        vars = 'SCT_snn_res.0.8'))
# 
# ##dom vs sub
# DimPlot(burtoni.snseq.combined.sct.all.subset.recluster, 
#         reduction = "umap", 
#         group.by = "orig.ident")
# ggsave('cell.subset/DomVsSub.dimplot.cell.subset.png')
# 
# ## clusters
# DimPlot(burtoni.snseq.combined.sct.all.subset.recluster, 
#         reduction = "umap", 
#         label = TRUE,
#         repel = TRUE, 
#         group.by = "integrated_snn_res.0.8")
# ggsave('cell.subset/Clusters.dimplot.cell.subset.png',
#        width = 10,
#        height = 10)
# 
# ## cell type
# DimPlot(burtoni.snseq.combined.sct.all.subset.recluster, 
#         reduction = "umap", 
#         label = TRUE,
#         repel = TRUE, 
#         group.by = "Cell.type")
# ggsave('cell.subset/Cell.type.dimplot.cell.subset.png',
#        width = 10,
#        height = 10)
# 
# ## broad cell class
# DimPlot(burtoni.snseq.combined.sct.all.subset.recluster, 
#         reduction = "umap", 
#         label = TRUE,
#         repel = TRUE, 
#         group.by = "Cell.class.broad")
# ggsave('cell.subset/Cell.class.broad.dimplot.cell.subset.png',
#        width = 10,
#        height = 10)
# 
# 
# ##expression level
# burtoni.snseq.combined.sct.all.subset.recluster.expression = full_join(full_join(burtoni.snseq.combined.sct.all.subset.recluster@reductions$umap@cell.embeddings %>% 
#                                                                                    as.data.frame() %>% 
#                                                                                    rownames_to_column("Cell.id"),
#                                                                                  burtoni.snseq.combined.sct.all.subset.recluster@meta.data %>% 
#                                                                                    rownames_to_column("Cell.id")),
#                                                                        burtoni.snseq.combined.sct.all.subset.recluster@assays$integrated@scale.data %>% 
#                                                                          as.data.frame()  %>% 
#                                                                          filter(rownames(burtoni.snseq.combined.sct.all.nucb2.neurons.recluster@assays$integrated@scale.data) %in% c('avp')) %>% 
#                                                                          t() %>% 
#                                                                          as.data.frame() %>% 
#                                                                          rownames_to_column('Cell.id')) %>% 
#   full_join(burtoni.snseq.combined.sct.all.subset.recluster@reductions$pca@cell.embeddings %>% 
#               as.data.frame() %>% 
#               select(PC_1,
#                      PC_2,
#                      PC_3,
#                      PC_4,
#                      PC_5) %>% 
#               rownames_to_column("Cell.id"))
# 
# # seperation by cell type
# table(burtoni.snseq.combined.sct.all.subset.recluster.expression %>% 
#         select(Cell.type,
#                integrated_snn_res.0.8))
# 
# # seperation by cell class broad
# table(burtoni.snseq.combined.sct.all.subset.recluster.expression %>% 
#         select(Cell.class.broad,
#                integrated_snn_res.0.8))
# 
# # difference social status
# #bias
# # 12500 to 10004
# table(burtoni.snseq.combined.sct.all.subset.recluster.expression %>% 
#         select(orig.ident))
# 
# #compare across clusters
# table(burtoni.snseq.combined.sct.all.subset.recluster.expression %>% 
#         select(integrated_snn_res.0.8,
#                orig.ident)) %>% 
#   as.data.frame() %>% 
#   pivot_wider(names_from = orig.ident,
#               values_from = Freq) %>% 
#   mutate(Dom.sub.ratio = dom_burtoni_snseq/sub_burtoni_snseq) %>% 
#   filter(Dom.sub.ratio > 1.25) 
# 
# #graph
# table(burtoni.snseq.combined.sct.all.subset.recluster.expression %>% 
#         select(integrated_snn_res.0.8,
#                orig.ident)) %>% 
#   as.data.frame() %>% 
#   pivot_wider(names_from = orig.ident,
#               values_from = Freq) %>% 
#   mutate(Dom.sub.ratio = dom_burtoni_snseq/sub_burtoni_snseq,
#          Total = dom_burtoni_snseq + sub_burtoni_snseq) %>% 
#   ggplot(aes(x = Dom.sub.ratio,
#              y=Total,
#              label = integrated_snn_res.0.8)) +
#   geom_label() +
#   theme_classic() +
#   geom_vline(xintercept = 1.25)
# ggsave('cell.subset/Cluster.bias.DomvsSub.cell.subset.png',
#        width = 10,
#        height = 10)
# 
# ## cell class across clusters
# table(burtoni.snseq.combined.sct.all.subset.recluster.expression %>% 
#         select(integrated_snn_res.0.8,
#                Cell.class.broad)) %>% 
#   as.data.frame() %>% 
#   ggplot(aes(x = Freq,
#              y = integrated_snn_res.0.8,
#              fill = Cell.class.broad)) +
#   geom_bar(stat = 'identity') +
#   theme_classic() 
# ggsave('cell.subset/Cluster.cell.class.broad.cell.subset.png',
#        width = 10,
#        height = 10)
# 
# # cell class across clusters
# table(burtoni.snseq.combined.sct.all.subset.recluster.expression %>% 
#         select(integrated_snn_res.0.8,
#                Cell.class.broad)) %>% 
#   as.data.frame() %>% 
#   group_by(integrated_snn_res.0.8) %>% 
#   mutate(Total.cluster = sum(Freq)) %>% 
#   ungroup() %>% 
#   mutate(Percent = 100*Freq/Total.cluster) %>% 
#   ggplot(aes(x = Percent,
#              y = integrated_snn_res.0.8,
#              fill = Cell.class.broad)) +
#   geom_bar(stat = 'identity') +
#   theme_classic() 
# ggsave('cell.subset/Cluster.cell.class.broad.cell.subset.percent.png',
#        width = 10,
#        height = 10)
# 
# 
# 
# #### neuropeptides ####
# ## neuropeptides
# #avp vs oxt
# ## umap plots
# FeaturePlot(object = burtoni.snseq.combined.sct.all, 
#             features = c("oxt",
#                          "avp"),
#             pt.size = 0.5,
#             min.cutoff = 4,
#             max.cutoff = 5,
#             order = TRUE,
#             blend = T)
# ggsave('avp.oxt/Avp.Oxt.umap.png',
#        width = 20,
#        height = 20)
# #across samples
# FeaturePlot(object = burtoni.snseq.combined.sct.all, 
#             features = c("oxt",
#                          "avp"),
#             pt.size = 0.5,
#             min.cutoff = 4,
#             max.cutoff = 6,
#             order = TRUE,
#             blend = T,
#             split.by = "orig.ident")
# ggsave('avp.oxt/Avp.Oxt.subvsdom.umap.png',
#        width = 10,
#        height = 10)
# 
# 
# 
# ##vinplot
# #avp
# VlnPlot(burtoni.snseq.combined.sct.all, 
#         features = c("avp"),
#         split.by = "orig.ident") +
#   ylim(4,25)
# ggsave('avp.oxt/Avp.subvsdom.vlnplot.png',
#        width = 10,
#        height = 10)
# #oxt
# VlnPlot(burtoni.snseq.combined.sct.all, 
#         features = c("oxt"),
#         split.by = "orig.ident") +
#   ylim(4,25)
# ggsave('avp.oxt/Oxt.subvsdom.vlnplot.png',
#        width = 10,
#        height = 10)
# 
# 
# ##check number of cells expressing gene above threshold
# #oxt
# #141
# sum(GetAssayData(object = burtoni.snseq.combined.sct.all, 
#                  slot = "data")["oxt",]>4)
# #avp
# #157
# sum(GetAssayData(object = burtoni.snseq.combined.sct.all, 
#                  slot = "data")["avp",]>4)
# 
# 
# #create list of oxt cells?
# # oxt.cell.ids = WhichCells(burtoni.snseq.combined.sct.all,
# #                           subset.name = )
# 
# ### subset data to relevant cells
# burtoni.snseq.combined.sct.all.avp.oxt = subset(burtoni.snseq.combined.sct.all,
#                                                 subset = avp > 4 | oxt > 4)
# 
# ## run PCA, UMAP, and cluster 
# #use 0.6 resolution
# burtoni.snseq.combined.sct.all.avp.oxt.recluster = burtoni.snseq.combined.sct.all.avp.oxt %>% 
#   RunPCA() %>%
#   FindNeighbors(dims = 1:15) %>%
#   RunUMAP(dims = 1:15) %>%
#   FindClusters(resolution = 1)
# 
# ## check cells per group per cluster
# #get table
# Avp.cells.per.cluster.table = table(burtoni.snseq.combined.sct.all.avp.oxt.recluster@active.ident, 
#                                     burtoni.snseq.combined.sct.all.avp.oxt.recluster@meta.data$orig.ident) %>% 
#   as.data.frame.matrix() %>% 
#   rownames_to_column("cluster.id") %>% 
#   pivot_longer(cols = c("dom_burtoni_snseq",
#                         "sub_burtoni_snseq"),
#                names_to = "orig.ident",
#                values_to = "count") %>% 
#   group_by(orig.ident) %>% 
#   mutate(total = sum(count)) %>% 
#   ungroup() %>% 
#   mutate(percentage = 100*count/total)
# 
# #graph
# Avp.cells.per.cluster.table %>% 
#   as.data.frame() %>% 
#   mutate(cluster.id = as.numeric(cluster.id)) %>% 
#   ggplot(aes(x = cluster.id,
#              y = percentage,
#              group = orig.ident,
#              color = orig.ident)) +
#   geom_point() +
#   theme_bw()
# ggsave('avp.oxt/Cells.per.cluster.DomvsSub.avp.oxt.png',
#        width = 10,
#        height = 10)
# 
# ### graph 
# ##dom vs sub
# DimPlot(burtoni.snseq.combined.sct.all.avp.oxt.recluster, 
#         reduction = "umap", 
#         group.by = "orig.ident")
# # ggsave('avp.oxt/DomVsSub.dimplot.avp.oxt.all.png')
# 
# ## clusters
# DimPlot(burtoni.snseq.combined.sct.all.avp.oxt.recluster, 
#         reduction = "umap", 
#         label = TRUE,
#         repel = TRUE)
# ggsave('avp.oxt/Clusters.dimplot.avp.oxt.all.png',
#        width = 10,
#        height = 10)
# 
# #avp vs oxt
# ## umap plots
# FeaturePlot(object = burtoni.snseq.combined.sct.all.avp.oxt.recluster, 
#             features = c("oxt",
#                          "avp"),
#             pt.size = 0.5,
#             max.cutoff = 6,
#             order = TRUE,
#             label = T,
#             blend = T)
# ggsave('avp.oxt/Avp.Oxt.clustering.umap.png',
#        width = 15,
#        height = 7)
# 
# #across samples
# DimPlot(burtoni.snseq.combined.sct.all.avp.oxt.recluster, 
#         reduction = "umap",
#         group.by = 'orig.ident',
#         pt.size = 1) +
#   theme(legend.position = 'none') 
# ggsave('avp.oxt/Avp.Oxt.clustering.umap.DomvsSub.png',
#        width = 10,
#        height = 10)
# 
# 
# ##pca
# # clusters
# DimPlot(burtoni.snseq.combined.sct.all.avp.oxt.recluster, 
#         reduction = "pca", 
#         label = TRUE,
#         repel = TRUE)
# ggsave('avp.oxt/Avp.Oxt.pca.PC1.PC2.png',
#        width = 10,
#        height = 10)
# #PC 2 and 3
# DimPlot(burtoni.snseq.combined.sct.all.avp.oxt.recluster, 
#         reduction = "pca", 
#         label = TRUE,
#         repel = TRUE,
#         dims = c(3,2))
# ggsave('avp.oxt/Avp.Oxt.pca.PC2.PC3.png',
#        width = 10,
#        height = 10)
# #avp and oxt
# FeaturePlot(object = burtoni.snseq.combined.sct.all.avp.oxt.recluster, 
#             features = c("oxt"),
#             pt.size = 3,
#             max.cutoff = 6,
#             order = TRUE,
#             label = T,
#             reduction = 'pca',
#             cols = c('blue',
#                      'red'),
#             dims = c(3,2))+
#   theme(text = element_text(size = 30))
# ggsave('avp.oxt/Avp.Oxt.pca.PC2.PC3.expression.png',
#        width = 10,
#        height = 10)
# 
# #avp and oxt
# #no label
# FeaturePlot(object = burtoni.snseq.combined.sct.all.avp.oxt.recluster, 
#             features = c("avp"),
#             pt.size = 3,
#             max.cutoff = 6,
#             order = TRUE,
#             reduction = 'pca',
#             cols = c('blue',
#                      'red'),
#             dims = c(3,2))+
#   theme(text = element_text(size = 30))
# ggsave('avp.oxt/Avp.Oxt.pca.PC2.PC3.expression.nolabel.png',
#        width = 10,
#        height = 10)
# 
# #heatmap PCA
# DimHeatmap(burtoni.snseq.combined.sct.all.avp.oxt.recluster, 
#            dims = 1:4, 
#            cells = 500, 
#            balanced = TRUE)
# ggsave('avp.oxt/Avp.Oxt.pca.heatmap.png',
#        width = 10,
#        height = 10)
# #loadings
# VizDimLoadings(burtoni.snseq.combined.sct.all.avp.oxt.recluster, 
#                dims = 1:4, 
#                reduction = "pca")
# ggsave('avp.oxt/Avp.Oxt.pca.loadings.png',
#        width = 10,
#        height = 10)
# 
# #calculate variance
# pca <- burtoni.snseq.combined.sct.all.avp.oxt.recluster[["pca"]]
# total_variance <- slot(burtoni.snseq.combined.sct.all.avp.oxt.recluster[["pca"]], "misc")$total.variance
# eigValues = (pca@stdev)^2  ## EigenValues
# varExplained = eigValues / total_variance
# 
# 
# 
# ## marker genes clusters
# burtoni.snseq.combined.sct.all.avp.oxt.recluster.markers <- FindAllMarkers(burtoni.snseq.combined.sct.all.avp.oxt.recluster, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# 
# burtoni.snseq.combined.sct.all.avp.oxt.recluster.markers.top10 = burtoni.snseq.combined.sct.all.avp.oxt.recluster.markers %>%
#   group_by(cluster) %>%
#   top_n(n = 10, wt = avg_log2FC) %>% 
#   arrange(avg_log2FC,
#           .by_group = TRUE)
# 
# DoHeatmap(burtoni.snseq.combined.sct.all.avp.oxt.recluster, features = burtoni.snseq.combined.sct.all.avp.oxt.recluster.markers.top10$gene) + NoLegend()
# ggsave('avp.oxt/Avp.Oxt.marker.heatmap.top10.png',
#        width = 10,
#        height = 10)
# 
# 
# DoHeatmap(burtoni.snseq.combined.sct.all.avp.oxt.recluster, features = burtoni.snseq.combined.sct.all.avp.oxt.recluster.markers %>%
#             group_by(cluster) %>%
#             arrange(avg_log2FC,
#                     .by_group = TRUE) %>% 
#             pull(gene)) + NoLegend()
# ggsave('avp.oxt/Avp.Oxt.marker.heatmap.png',
#        width = 10,
#        height = 10)
# 
# 
# # #create dendrogram
# # BuildClusterTree(burtoni.snseq.combined.sct.all.avp.oxt.recluster) %>% PlotClusterTree()
# # ggsave('avp.oxt/Avp.Oxt.dendrogram.png',
# #        width = 10,
# #        height = 10)
# 
# 
# ##create cell type dendrogram
# #remove rows with all zeros 
# burtoni.snseq.combined.sct.all.avp.oxt.recluster.data <- as.data.frame(log1p(AverageExpression(burtoni.snseq.combined.sct.all.avp.oxt.recluster, verbose = FALSE)$RNA))
# burtoni.snseq.combined.sct.all.avp.oxt.recluster.data$gene <- rownames(burtoni.snseq.combined.sct.all.avp.oxt.recluster.data)
# #create dendrogram
# burtoni.snseq.combined.sct.all.avp.oxt.recluster.tree = pvclust(burtoni.snseq.combined.sct.all.avp.oxt.recluster.data %>% 
#                                                                   select(-c(gene)) %>%                                                            mutate(sum = rowSums(.)) %>% 
#                                                                   filter(sum > 0) %>% 
#                                                                   select(-c(sum)))
# 
# #graph
# png(file = 'avp.oxt/Avp.Oxt.dendrogram.png',
#     width = 10,
#     height = 10,
#     units = 'in',
#     res = 480)
# plot(burtoni.snseq.combined.sct.all.avp.oxt.recluster.tree)
# dev.off()
# 
# ## gene expression between nodes
# burtoni.snseq.combined.sct.all.avp.oxt.recluster.tree.domvssub = burtoni.snseq.combined.sct.all.avp.oxt.recluster
# 
# Idents(burtoni.snseq.combined.sct.all.avp.oxt.recluster.tree.domvssub) <- "orig.ident"
# 
# burtoni.snseq.combined.sct.all.avp.oxt.recluster.tree.domvssub <- as.data.frame(log1p(AverageExpression(burtoni.snseq.combined.sct.all.avp.oxt.recluster.tree.domvssub, verbose = FALSE)$RNA))
# burtoni.snseq.combined.sct.all.avp.oxt.recluster.tree.domvssub$gene <- rownames(burtoni.snseq.combined.sct.all.avp.oxt.recluster.tree.domvssub)
# 
# ggplot(burtoni.snseq.combined.sct.all.avp.oxt.recluster.tree.domvssub, 
#        aes(sub_burtoni_snseq, 
#            dom_burtoni_snseq)) + 
#   geom_point() +
#   geom_abline(intercept = 0,
#               slope = 1) +
#   theme_bw() +
#   xlim(0,200) +
#   ylim(0,200)
# # ggsave('avp.oxt/Avp.domvssub.expression.png',
# #        width = 10,
# #        height = 10)
# 
# ## compare gene expression
# #avp 1
# burtoni.snseq.combined.sct.all.avp.subset = subset(burtoni.snseq.combined.sct.all.avp.oxt.recluster,
#                                                    idents = c("1"))
# 
# Idents(burtoni.snseq.combined.sct.all.avp.subset) <- "orig.ident"
# 
# burtoni.snseq.combined.sct.all.avp.subset <- as.data.frame(log1p(AverageExpression(burtoni.snseq.combined.sct.all.avp.subset, verbose = FALSE)$RNA))
# burtoni.snseq.combined.sct.all.avp.subset$gene <- rownames(burtoni.snseq.combined.sct.all.avp.subset)
# 
# ggplot(burtoni.snseq.combined.sct.all.avp.subset, 
#        aes(sub_burtoni_snseq, 
#            dom_burtoni_snseq)) + 
#   geom_point() +
#   geom_abline(intercept = 0,
#               slope = 1) +
#   theme_bw()
# ggsave('avp.oxt/Avp.1.domvssub.expression.png',
#        width = 10,
#        height = 10)
# 
# 
# #avp 3, 4, 5
# burtoni.snseq.combined.sct.all.avp.subset.multi = subset(burtoni.snseq.combined.sct.all.avp.oxt.recluster,
#                                                          idents = c("3",
#                                                                     "4",
#                                                                     "5"))
# 
# Idents(burtoni.snseq.combined.sct.all.avp.subset.multi) <- "orig.ident"
# 
# burtoni.snseq.combined.sct.all.avp.subset.multi <- as.data.frame(log1p(AverageExpression(burtoni.snseq.combined.sct.all.avp.subset.multi, verbose = FALSE)$RNA))
# burtoni.snseq.combined.sct.all.avp.subset.multi$gene <- rownames(burtoni.snseq.combined.sct.all.avp.subset.multi)
# 
# ggplot(burtoni.snseq.combined.sct.all.avp.subset.multi, 
#        aes(sub_burtoni_snseq, 
#            dom_burtoni_snseq)) + 
#   geom_point() +
#   geom_abline(intercept = 0,
#               slope = 1) +
#   theme_bw() 
# ggsave('avp.oxt/Avp.3.4.5.domvssub.expression.png',
#        width = 10,
#        height = 10)
# 
# # LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
# 
# 
# ## create volcano plots
# #across social status
# #AVP 3 
# burtoni.snseq.combined.sct.all.avp.oxt.recluster.3 = subset(burtoni.snseq.combined.sct.all.avp.oxt.recluster,
#                                                             idents = 3)
# 
# burtoni.snseq.combined.sct.all.avp.oxt.recluster.3 = subset(burtoni.snseq.combined.sct.all.avp.oxt.recluster.3,
#                                                             subset = avp > 4 )
# #change ident
# Idents(burtoni.snseq.combined.sct.all.avp.oxt.recluster.3) <- "orig.ident"
# #find markers across samples
# burtoni.snseq.combined.sct.all.avp.oxt.recluster.3.markers = FindMarkers(object = burtoni.snseq.combined.sct.all.avp.oxt.recluster.3,
#                                                                          ident.1 = "dom_burtoni_snseq", 
#                                                                          ident.2 = "sub_burtoni_snseq",
#                                                                          test.use = "DESeq2",
#                                                                          assay = 'RNA') 
# 
# #graph volcano plot
# burtoni.snseq.combined.sct.all.avp.oxt.recluster.3.markers %>% 
#   mutate(Sig = ifelse(p_val < 0.01 & avg_log2FC > 0.5,
#                       "significant dom",
#                       ifelse(p_val < 0.01 & avg_log2FC < -0.5,
#                              "significant sub",
#                              "not significant")))%>% 
#   rownames_to_column(var = 'gene') %>% 
#   mutate(gene.id = ifelse(Sig != "not significant",
#                           gene,
#                           NA)) %>%
#   ggplot(aes(x = avg_log2FC,
#              y = -log10(p_val),
#              color = Sig)) +
#   geom_point(size = 3) +
#   theme_bw() +
#   xlim(-2.5, 2.5) +
#   ggtitle("Volcano plot AVP.3") +
#   scale_color_manual(values = c('black',
#                                 "red",
#                                 "blue")) +
#   geom_vline(xintercept = -0.5,
#              linetype="dotted")+
#   geom_vline(xintercept = 0.5,
#              linetype="dotted") +
#   geom_hline(yintercept = -log10(0.01),
#              linetype="dotted")+ 
#   geom_text_repel(aes(label = gene.id),
#                   size = 5) +
#   theme(text = element_text(size = 30),)   
# ggsave('avp.oxt/Avp.3.domvssub.volcano.png',
#        width = 10,
#        height = 10)
# 
# 
# 
# #oxt 0 
# burtoni.snseq.combined.sct.all.avp.oxt.recluster.0 = subset(burtoni.snseq.combined.sct.all.avp.oxt.recluster,
#                                                             idents = 0)
# 
# burtoni.snseq.combined.sct.all.avp.oxt.recluster.0 = subset(burtoni.snseq.combined.sct.all.avp.oxt.recluster.0,
#                                                             subset = oxt > 4 )
# #change ident
# Idents(burtoni.snseq.combined.sct.all.avp.oxt.recluster.0) <- "orig.ident"
# #find markers across samples
# burtoni.snseq.combined.sct.all.avp.oxt.recluster.0.markers = FindMarkers(object = burtoni.snseq.combined.sct.all.avp.oxt.recluster.0,
#                                                                          ident.1 = "dom_burtoni_snseq", 
#                                                                          ident.2 = "sub_burtoni_snseq") 
# 
# burtoni.snseq.combined.sct.all.avp.oxt.recluster.0.markers = FindMarkers(object = burtoni.snseq.combined.sct.all.avp.oxt.recluster.0,
#                                                                          ident.1 = "dom_burtoni_snseq", 
#                                                                          ident.2 = "sub_burtoni_snseq",
#                                                                          test.use = "DESeq2",
#                                                                          assay = 'RNA') 
# 
# #graph volcano plot
# burtoni.snseq.combined.sct.all.avp.oxt.recluster.0.markers %>% 
#   mutate(Sig = ifelse(p_val < 0.01 & avg_log2FC > 0.5,
#                       "significant dom",
#                       ifelse(p_val < 0.01 & avg_log2FC < -0.5,
#                              "significant sub",
#                              "not significant"))) %>% 
#   rownames_to_column(var = 'gene') %>% 
#   mutate(gene.id = ifelse(Sig != "not significant",
#                           gene,
#                           NA)) %>% 
#   ggplot(aes(x = avg_log2FC,
#              y = -log10(p_val),
#              color = Sig)) +
#   geom_point(size=3) +
#   theme_bw() +
#   xlim(-2.5, 2.5) +
#   ggtitle("Volcano plot oxt.0") +
#   scale_color_manual(values = c('black',
#                                 "red",
#                                 "blue")) +
#   geom_vline(xintercept = -0.5,
#              linetype="dotted")+
#   geom_vline(xintercept = 0.5,
#              linetype="dotted") +
#   geom_hline(yintercept = -log10(0.01),
#              linetype="dotted") + 
#   geom_text_repel(aes(label = gene.id),
#                   size = 5) +
#   theme(text = element_text(size = 30))    
# ggsave('avp.oxt/oxt.0.domvssub.volcano.png',
#        width = 10,
#        height = 10)
# 
# 
# 
# #AVP 1
# burtoni.snseq.combined.sct.all.avp.oxt.recluster.1 = subset(burtoni.snseq.combined.sct.all.avp.oxt.recluster,
#                                                             idents = 1)
# 
# burtoni.snseq.combined.sct.all.avp.oxt.recluster.1 = subset(burtoni.snseq.combined.sct.all.avp.oxt.recluster.1,
#                                                             subset = avp > 4 )
# #change ident
# Idents(burtoni.snseq.combined.sct.all.avp.oxt.recluster.1) <- "orig.ident"
# #find markers across samples
# burtoni.snseq.combined.sct.all.avp.oxt.recluster.1.markers = FindMarkers(object = burtoni.snseq.combined.sct.all.avp.oxt.recluster.1,
#                                                                          ident.1 = "dom_burtoni_snseq", 
#                                                                          ident.2 = "sub_burtoni_snseq",
#                                                                          test.use = "DESeq2",
#                                                                          assay = 'RNA') 
# 
# #graph volcano plot
# burtoni.snseq.combined.sct.all.avp.oxt.recluster.1.markers %>% 
#   mutate(Sig = ifelse(p_val < 0.01 & avg_log2FC > 0.5,
#                       "significant dom",
#                       ifelse(p_val < 0.01 & avg_log2FC < -0.5,
#                              "significant sub",
#                              "not significant")))%>% 
#   rownames_to_column(var = 'gene') %>% 
#   mutate(gene.id = ifelse(Sig != "not significant",
#                           gene,
#                           NA)) %>%
#   ggplot(aes(x = avg_log2FC,
#              y = -log10(p_val),
#              color = Sig)) +
#   geom_point() +
#   theme_bw() +
#   xlim(-2.5, 2.5) +
#   ggtitle("Volcano plot AVP.1") +
#   scale_color_manual(values = c('black',
#                                 "red",
#                                 "blue")) +
#   geom_vline(xintercept = -0.5,
#              linetype="dotted")+
#   geom_vline(xintercept = 0.5,
#              linetype="dotted") +
#   geom_hline(yintercept = -log10(0.01),
#              linetype="dotted")+ 
#   geom_text_repel(aes(label = gene.id)) 
# ggsave('avp.oxt/Avp.1.domvssub.volcano.png',
#        width = 10,
#        height = 10)
# 
# 
# 
# #oxt 2 
# burtoni.snseq.combined.sct.all.avp.oxt.recluster.2 = subset(burtoni.snseq.combined.sct.all.avp.oxt.recluster,
#                                                             idents = 2)
# 
# burtoni.snseq.combined.sct.all.avp.oxt.recluster.2 = subset(burtoni.snseq.combined.sct.all.avp.oxt.recluster.2,
#                                                             subset = oxt > 4 )
# #change ident
# Idents(burtoni.snseq.combined.sct.all.avp.oxt.recluster.2) <- "orig.ident"
# #find markers across samples
# burtoni.snseq.combined.sct.all.avp.oxt.recluster.2.markers = FindMarkers(object = burtoni.snseq.combined.sct.all.avp.oxt.recluster.2,
#                                                                          ident.1 = "dom_burtoni_snseq", 
#                                                                          ident.2 = "sub_burtoni_snseq") 
# 
# burtoni.snseq.combined.sct.all.avp.oxt.recluster.2.markers = FindMarkers(object = burtoni.snseq.combined.sct.all.avp.oxt.recluster.2,
#                                                                          ident.1 = "dom_burtoni_snseq", 
#                                                                          ident.2 = "sub_burtoni_snseq",
#                                                                          test.use = "DESeq2",
#                                                                          assay = 'RNA') 
# 
# #graph volcano plot
# burtoni.snseq.combined.sct.all.avp.oxt.recluster.2.markers %>% 
#   mutate(Sig = ifelse(p_val < 0.01 & avg_log2FC > 0.5,
#                       "significant dom",
#                       ifelse(p_val < 0.01 & avg_log2FC < -0.5,
#                              "significant sub",
#                              "not significant"))) %>% 
#   rownames_to_column(var = 'gene') %>% 
#   mutate(gene.id = ifelse(Sig != "not significant",
#                           gene,
#                           NA)) %>% 
#   ggplot(aes(x = avg_log2FC,
#              y = -log10(p_val),
#              color = Sig)) +
#   geom_point() +
#   theme_bw() +
#   xlim(-0.5, 0.5) +
#   ggtitle("Volcano plot oxt.2") +
#   scale_color_manual(values = c('black',
#                                 "red",
#                                 "blue")) +
#   geom_vline(xintercept = -0.5,
#              linetype="dotted")+
#   geom_vline(xintercept = 0.5,
#              linetype="dotted") +
#   geom_hline(yintercept = -log10(0.01),
#              linetype="dotted") + 
#   geom_text_repel(aes(label = gene.id)) 
# ggsave('avp.oxt/oxt.2.domvssub.volcano.png',
#        width = 10,
#        height = 10)
# 
# 
# ###dotplot
# ##oxt and avp
# DotPlot(burtoni.snseq.combined.sct.all.avp.oxt.recluster, 
#         features = c("avp",
#                      "oxt"), 
#         cols = c("grey", 
#                  "red"), 
#         dot.scale = 8,
#         col.max = 1,
#         dot.min = .5) + 
#   RotatedAxis() +
#   ggtitle("Avp.Oxt")
# ggsave('avp.oxt/Avp.Oxt.clustering.dotplot.png',
#        width = 10,
#        height = 10)
# 
# ## cell type
# DotPlot(burtoni.snseq.combined.sct.all.avp.oxt.recluster, 
#         features = burtoni.snseq.combined.sct.all.avp.oxt.recluster@meta.data %>% 
#           select(ends_with(".score1")) %>% 
#           colnames(), 
#         cols = c("grey", 
#                  "red"), 
#         dot.scale = 8,
#         col.min = 1,
#         dot.min = .4) + 
#   RotatedAxis()
# ggsave('avp.oxt/Avp.Oxt.celltype.dotplot.png',
#        width = 10,
#        height = 10)
# 
# ##markers parvo vs magno
# DotPlot(burtoni.snseq.combined.sct.all.avp.oxt.recluster, 
#         features = c('KCNMB4',
#                      "RELN",
#                      "cnr1",
#                      "CNR1"), 
#         cols = c("grey", 
#                  "red"),
#         col.min =  1,
#         dot.min = .1) + 
#   RotatedAxis() +
#   ggtitle("Avp.Oxt")
# ggsave('avp.oxt/Avp.Oxt.magno.parvo.markers.dotplot.png',
#        width = 10,
#        height = 10)
# 
# ### clustree
# # cluster across resolutions
# burtoni.snseq.combined.sct.all.avp.oxt.recluster.clustree <- Seurat::FindClusters(object = burtoni.snseq.combined.sct.all.avp.oxt.recluster, 
#                                                                                   resolution = resolution.range.reduced)
# #check data
# head(burtoni.snseq.combined.sct.all.avp.oxt.recluster.clustree[[]])
# #set presentation colors
# presentation.color <- c('#66c2a5',
#                         '#fc8d62',
#                         '#8da0cb',
#                         '#e78ac3',
#                         '#a6d854',
#                         '#ffd92f',
#                         '#e5c494',
#                         '#b3b3b3')
# #clustree
# clustree(burtoni.snseq.combined.sct.all.avp.oxt.recluster.clustree, 
#          prefix = "integrated_snn_res.",
#          node_colour = 'cluster',
#          node_size_range = c(10,20),
#          scale_node_text = TRUE) +
#   scale_edge_color_continuous(low = "black", 
#                               high = "black") +
#   scale_color_manual(values = presentation.color)+
#   theme(legend.position = "bottom")
# ggsave('avp.oxt/avp.oxt.clustree.png',
#        width = 10,
#        height = 10)
# 
# ## loop through resolutions for presentation
# for (x in resolution.range.reduced) {
#   DimPlot(burtoni.snseq.combined.sct.all.avp.oxt.recluster.clustree, 
#           reduction = "umap",
#           group.by = paste('integrated_snn_res',
#                            x,
#                            sep = '.'),
#           label = T,
#           pt.size = 2,
#           label.size = 10,
#           repel = T) +
#     theme(legend.position = 'none') +
#     scale_color_manual(values = presentation.color)
#   ggsave(paste('avp.oxt/resolution/Avp.oxt.resolution',
#                x,
#                'dimplot.avp.oxt.all.png',
#                sep = '.'),
#          width = 10,
#          height = 10)
# }
# 
# ### identify marker genes
# ## all
# burtoni.snseq.combined.sct.all.avp.oxt.recluster.clustree.markers <- FindAllMarkers(burtoni.snseq.combined.sct.all.avp.oxt.recluster.clustree, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# 
# ##compare clusters
# #AVP
# #1 vs 3
# burtoni.snseq.combined.sct.all.avp.oxt.recluster.clustree.markers.1v3 <- FindMarkers(burtoni.snseq.combined.sct.all.avp.oxt.recluster.clustree, 
#                                                                                      ident.1 = 1, 
#                                                                                      ident.2 = 3, 
#                                                                                      min.pct = 0.5, 
#                                                                                      logfc.threshold = 1)
# 
# #OXT
# #2 vs 0
# burtoni.snseq.combined.sct.all.avp.oxt.recluster.clustree.markers.2v0 <- FindMarkers(burtoni.snseq.combined.sct.all.avp.oxt.recluster.clustree, 
#                                                                                      ident.1 = 2, 
#                                                                                      ident.2 = 0, 
#                                                                                      min.pct = 0.55, 
#                                                                                      logfc.threshold = 1)
# 
# ## compare concordance across groups
# burtoni.snseq.combined.sct.all.avp.oxt.recluster.clustree.markers.compare = full_join(burtoni.snseq.combined.sct.all.avp.oxt.recluster.clustree.markers.2v0 %>% 
#                                                                                         rownames_to_column("gene"),
#                                                                                       burtoni.snseq.combined.sct.all.avp.oxt.recluster.clustree.markers.1v3 %>% 
#                                                                                         rownames_to_column("gene"),
#                                                                                       by = "gene",
#                                                                                       suffix = c(".2v0",
#                                                                                                  ".1v3"))
# 
# ## graph concordance
# burtoni.snseq.combined.sct.all.avp.oxt.recluster.clustree.markers.compare %>% 
#   ggplot(aes(x = avg_log2FC.2v0,
#              y = avg_log2FC.1v3)) +
#   geom_point() +
#   theme_classic() +
#   geom_hline(yintercept = 0) +
#   geom_vline(xintercept = 0)
# 
# 
# ## 2 v 0 
# burtoni.snseq.combined.sct.all.avp.oxt.recluster.clustree.markers.compare %>% 
#   mutate(empty.1v3 = ifelse(is.na(avg_log2FC.1v3),
#                             "NA",
#                             "Overlap")) %>% 
#   ggplot(aes(x = avg_log2FC.2v0,
#              y = -log(p_val_adj.2v0),
#              color = empty.1v3)) +
#   geom_point() +
#   theme_classic()
# 
# ## 1 v 3 
# burtoni.snseq.combined.sct.all.avp.oxt.recluster.clustree.markers.compare %>% 
#   mutate(empty.2v0 = ifelse(is.na(avg_log2FC.2v0),
#                             "NA",
#                             "Overlap")) %>% 
#   ggplot(aes(x = avg_log2FC.1v3,
#              y = -log(p_val_adj.1v3),
#              color = empty.2v0)) +
#   geom_point() +
#   theme_classic()
# 
# # remove na
# burtoni.snseq.combined.sct.all.avp.oxt.recluster.clustree.markers.compare %>% 
#   na.omit() %>% 
#   View()
# 
# # remove pval adjust
# # remove na
# burtoni.snseq.combined.sct.all.avp.oxt.recluster.clustree.markers.compare %>% 
#   na.omit() %>% 
#   filter(p_val_adj.1v3<0.05) %>% 
#   filter(p_val_adj.2v0<0.05) %>% 
#   View()
# 
# ### combine with scsorter cell type
# ## subset highly variable genes
# ## identify top variable genes
# burtoni.snseq.combined.sct.all.avp.oxt.recluster.filter <- FindVariableFeatures(burtoni.snseq.combined.sct.all.avp.oxt.recluster, 
#                                                                                 selection.method = "vst", 
#                                                                                 nfeatures = 2000, 
#                                                                                 verbose = F)
# 
# # identify top 2000 variable genes
# topgenes <- head(VariableFeatures(burtoni.snseq.combined.sct.all.avp.oxt.recluster.filter), 
#                  2000)
# 
# #get assay data
# #as matrix
# burtoni.snseq.combined.sct.all.avp.oxt.recluster.filter.extract = GetAssayData(burtoni.snseq.combined.sct.all.avp.oxt.recluster.filter) %>% 
#   as.matrix()
# 
# # filter genes
# topgene_filter = rowSums(as.matrix(burtoni.snseq.combined.sct.all.avp.oxt.recluster.filter.extract)[topgenes, ]!=0) > ncol(burtoni.snseq.combined.sct.all.avp.oxt.recluster.filter.extract)*.1
# #create filter gene list
# topgenes = topgenes[topgene_filter]
# 
# # subset data by gene list
# burtoni.snseq.combined.sct.all.avp.oxt.recluster.filter.extract = burtoni.snseq.combined.sct.all.avp.oxt.recluster.filter.extract[rownames(burtoni.snseq.combined.sct.all.avp.oxt.recluster.filter.extract) %in% topgenes, ]
# 
# #subset annotation file
# burtoni.scsorter.data.scsort.output.avp.oxt = inner_join(burtoni.scsorter.data.scsort.output,
#                                                          burtoni.snseq.combined.sct.all.avp.oxt.recluster.filter.extract %>% 
#                                                            colnames() %>% 
#                                                            as.data.frame() %>% 
#                                                            rename('Cell.id' = '.'))
# 
# ##pheatmap
# # create heatmap using pheatmap
# pheatmap(burtoni.snseq.combined.sct.all.avp.oxt.recluster.filter.extract %>% 
#            t(),
#          annotation_row = burtoni.scsorter.data.scsort.output.avp.oxt %>% 
#            select(Cell.type,
#                   Cell.id) %>% 
#            column_to_rownames("Cell.id"),
#          scale = "column") +
#   ggtitle('heatmap top 2000 variable genes')
# ggsave('avp.oxt/heatmap cell types avp oxt top 2000 variable genes.png',
#        width = 10,
#        height = 10)
# 
# ## create UMAP plot 
# #extract umap data
# ## combine umap and metadata
# burtoni.snseq.combined.sct.all.avp.oxt.recluster.umap = full_join(burtoni.snseq.combined.sct.all.avp.oxt.recluster@reductions$umap@cell.embeddings %>% 
#                                                                     as.data.frame() %>% 
#                                                                     rownames_to_column("Cell.id"),
#                                                                   burtoni.snseq.combined.sct.all.avp.oxt.recluster@meta.data %>% 
#                                                                     rownames_to_column("Cell.id"))
# 
# #graph umap
# burtoni.snseq.combined.sct.all.avp.oxt.recluster.umap %>% 
#   ggplot(aes(x=UMAP_1,
#              y= UMAP_2,
#              color = Cell.type)) +
#   geom_point() +
#   theme_bw()
# ggsave('avp.oxt/scsorter/UMAP scsorter cell types avp oxt.png',
#        width = 10,
#        height = 10)
# 
# #check cell count
# table(burtoni.snseq.combined.sct.all.avp.oxt.recluster.umap$Cell.type)
# 
# #check by orig.ident
# table(burtoni.snseq.combined.sct.all.avp.oxt.recluster.umap %>% 
#         select(Cell.type,
#                orig.ident))
# 
# 
# #neurons vs non-neurons
# burtoni.snseq.combined.sct.all.avp.oxt.recluster.umap %>% 
#   mutate(Cell.type.neuron = ifelse(Cell.type %in% c('neurons',
#                                                     'excitatory',
#                                                     'inhibitory'),
#                                    "neurons",
#                                    "other")) %>% 
#   ggplot(aes(x=UMAP_1,
#              y= UMAP_2,
#              color = Cell.type.neuron,
#              shape = integrated_snn_res.1)) +
#   geom_point() +
#   theme_bw()
# ggsave('avp.oxt/scsorter/UMAP scsorter cell types neurons vs other avp oxt.png',
#        width = 10,
#        height = 10)
# 
# 
# 
# 
# 
# ##expression level
# burtoni.snseq.combined.sct.all.avp.oxt.recluster.expression = full_join(burtoni.snseq.combined.sct.all.avp.oxt.recluster.umap,
#                                                                         burtoni.snseq.combined.sct.all.avp.oxt.recluster@assays$integrated@scale.data %>% 
#                                                                           as.data.frame() %>% 
#                                                                           filter(rownames(burtoni.snseq.combined.sct.all.avp.oxt.recluster@assays$integrated@scale.data) %in% c('avp',
#                                                                                                                                                                                 'oxt')) %>% 
#                                                                           t() %>% as.data.frame() %>% 
#                                                                           rownames_to_column('Cell.id'))
# 
# #graph
# #oxt
# burtoni.snseq.combined.sct.all.avp.oxt.recluster.expression %>%
#   filter(oxt > 4) %>% 
#   ggplot(aes(x=Cell.type ,
#              y=avp)) + 
#   geom_violin() +
#   theme_bw()
# ggsave('avp.oxt/scsorter/oxt scsorter cell types.png',
#        width = 10,
#        height = 10)
# 
# #avp
# burtoni.snseq.combined.sct.all.avp.oxt.recluster.expression %>%
#   filter(avp > 4) %>% 
#   ggplot(aes(x=Cell.type ,
#              y=avp)) + 
#   geom_violin() +
#   theme_bw()
# ggsave('avp.oxt/scsorter/avp scsorter cell types.png',
#        width = 10,
#        height = 10)
# 
# #cell counts of oxt avp 
# table(burtoni.snseq.combined.sct.all.avp.oxt.recluster.expression %>% 
#         mutate(Oxt.cell = ifelse(oxt > 4 & avp > 4,
#                                  'both',
#                                  ifelse(oxt > 4,
#                                         'oxt',
#                                         'avp'))) %>% 
#         select(Cell.type,
#                Oxt.cell))
# 
# #cell counts of oxt avp cell neurons
# table(burtoni.snseq.combined.sct.all.avp.oxt.recluster.expression %>% 
#         filter(Cell.type %in% c("neurons", 
#                                 "excitatory",
#                                 "inhibitory")) %>% 
#         mutate(Oxt.cell = ifelse(oxt > 4 & avp > 4,
#                                  'both',
#                                  ifelse(oxt > 4,
#                                         'oxt',
#                                         'avp'))) %>% 
#         select(Cell.type,
#                Oxt.cell))
# 
# 
# 
# #### neuropeptides neurons ####
# ### avp oxt neurons
# ### subset data to relevant cells
# burtoni.snseq.combined.sct.all.avp.oxt.neurons = burtoni.snseq.combined.sct.all
# 
# #set idents
# Idents(object = burtoni.snseq.combined.sct.all.avp.oxt.neurons) <- "Cell.type"
# 
# #subset neurons
# burtoni.snseq.combined.sct.all.avp.oxt.neurons = subset(burtoni.snseq.combined.sct.all.avp.oxt.neurons,
#                                                         idents = c("neurons",
#                                                                    "excitatory",
#                                                                    "inhibitory"))
# 
# #subset expression
# burtoni.snseq.combined.sct.all.avp.oxt.neurons = subset(burtoni.snseq.combined.sct.all.avp.oxt.neurons,
#                                                         subset = avp > 4 | oxt > 4)
# 
# 
# ## run PCA, UMAP, and cluster 
# #use 0.4 resolution
# burtoni.snseq.combined.sct.all.avp.oxt.neurons.recluster = burtoni.snseq.combined.sct.all.avp.oxt.neurons %>% 
#   RunPCA() %>%
#   FindNeighbors(dims = 1:15) %>%
#   RunUMAP(dims = 1:15) %>%
#   FindClusters(resolution = 0.8)
# 
# ### clustree
# # cluster across resolutions
# burtoni.snseq.combined.sct.all.avp.oxt.neurons.clustree <- Seurat::FindClusters(object = burtoni.snseq.combined.sct.all.avp.oxt.neurons.recluster, 
#                                                                                 resolution = resolution.range)
# #check data
# head(burtoni.snseq.combined.sct.all.avp.oxt.neurons.clustree[[]])
# #set presentation colors
# presentation.color <- c('#66c2a5',
#                         '#fc8d62',
#                         '#8da0cb',
#                         '#e78ac3',
#                         '#a6d854',
#                         '#ffd92f',
#                         '#e5c494',
#                         '#b3b3b3')
# #clustree
# clustree(burtoni.snseq.combined.sct.all.avp.oxt.neurons.clustree, 
#          prefix = "integrated_snn_res.",
#          node_colour = 'cluster',
#          node_size_range = c(10,20),
#          scale_node_text = TRUE) +
#   scale_edge_color_continuous(low = "black", 
#                               high = "black") +
#   scale_color_manual(values = presentation.color)+
#   theme(legend.position = "bottom")
# ggsave('avp.oxt/avp/avp.oxt.neurons.clustree.png',
#        width = 10,
#        height = 10)
# 
# 
# ### graph 
# ##dom vs sub
# DimPlot(burtoni.snseq.combined.sct.all.avp.oxt.neurons.recluster, 
#         reduction = "umap", 
#         group.by = "orig.ident")
# ggsave('avp.oxt/avp/DomVsSub.dimplot.avp.oxt.neurons.all.png')
# 
# ## clusters
# DimPlot(burtoni.snseq.combined.sct.all.avp.oxt.neurons.recluster, 
#         reduction = "umap", 
#         label = TRUE,
#         repel = TRUE)
# ggsave('avp.oxt/avp/Clusters.dimplot.avp.oxt.neurons.all.png',
#        width = 10,
#        height = 10)
# 
# ##expression level
# burtoni.snseq.combined.sct.all.avp.oxt.neurons.recluster.expression = full_join(full_join(burtoni.snseq.combined.sct.all.avp.oxt.neurons.recluster@reductions$umap@cell.embeddings %>% 
#                                                                                             as.data.frame() %>% 
#                                                                                             rownames_to_column("Cell.id"),
#                                                                                           burtoni.snseq.combined.sct.all.avp.oxt.neurons.recluster@meta.data %>% 
#                                                                                             rownames_to_column("Cell.id")),
#                                                                                 burtoni.snseq.combined.sct.all.avp.oxt.neurons.recluster@assays$integrated@scale.data %>% 
#                                                                                   as.data.frame() %>% 
#                                                                                   filter(rownames(burtoni.snseq.combined.sct.all.avp.oxt.neurons.recluster@assays$integrated@scale.data) %in% c('avp',
#                                                                                                                                                                                                 'oxt')) %>% 
#                                                                                   t() %>% as.data.frame() %>% 
#                                                                                   rownames_to_column('Cell.id')) %>% 
#   full_join(burtoni.snseq.combined.sct.all.avp.oxt.neurons.recluster@reductions$pca@cell.embeddings %>% 
#               as.data.frame() %>% 
#               select(PC_1,
#                      PC_2,
#                      PC_3,
#                      PC_4,
#                      PC_5) %>% 
#               rownames_to_column("Cell.id")) %>% 
#   mutate(Oxt.cell = ifelse(oxt > 4 & avp > 4,
#                            'both',
#                            ifelse(oxt > 4,
#                                   'oxt',
#                                   'avp')))
# 
# 
# ##graph
# burtoni.snseq.combined.sct.all.avp.oxt.neurons.recluster.expression %>% 
#   filter(Oxt.cell != 'both') %>% 
#   ggplot(aes(x = UMAP_1,
#              y = UMAP_2,
#              color = integrated_snn_res.0.8,
#              shape = Oxt.cell)) +
#   geom_point() +
#   theme_bw()
# ggsave('avp.oxt/avp/Umap.avp.oxt.neurons.all.png',
#        width = 10,
#        height = 10)
# 
# burtoni.snseq.combined.sct.all.avp.oxt.neurons.recluster.expression %>% 
#   filter(Oxt.cell != 'both') %>% 
#   ggplot(aes(x = UMAP_1,
#              y = UMAP_2,
#              shape = integrated_snn_res.0.8,
#              color = Oxt.cell)) +
#   geom_point() +
#   theme_bw()
# ggsave('avp.oxt/avp/Umap.expression.avp.oxt.neurons.all.png',
#        width = 10,
#        height = 10)
# 
# ## cell count per cluster
# table(burtoni.snseq.combined.sct.all.avp.oxt.neurons.recluster.expression %>% 
#         filter(Oxt.cell != 'both') %>% 
#         select(integrated_snn_res.0.8,
#                Oxt.cell))
# 
# ### PCA
# #1 vs 2
# burtoni.snseq.combined.sct.all.avp.oxt.neurons.recluster.expression %>% 
#   filter(Oxt.cell != 'both') %>% 
#   ggplot(aes(x = PC_1,
#              y = PC_2,
#              color = integrated_snn_res.0.8,
#              shape = Oxt.cell)) +
#   geom_point() +
#   theme_bw()
# ggsave('avp.oxt/avp/PCA.1.2.avp.oxt.neurons.all.png',
#        width = 10,
#        height = 10)
# #2 vs 3
# burtoni.snseq.combined.sct.all.avp.oxt.neurons.recluster.expression %>% 
#   filter(Oxt.cell != 'both') %>% 
#   ggplot(aes(x = PC_3,
#              y = PC_2,
#              shape = integrated_snn_res.0.8,
#              color = Oxt.cell)) +
#   geom_point() +
#   theme_bw()
# ggsave('avp.oxt/avp/PCA.2.3.expression.avp.oxt.neurons.all.png',
#        width = 10,
#        height = 10)
# 
# burtoni.snseq.combined.sct.all.avp.oxt.neurons.recluster.expression %>% 
#   filter(Oxt.cell != 'both') %>% 
#   ggplot(aes(x = PC_3,
#              y = PC_2,
#              color = integrated_snn_res.0.8,
#              shape = Oxt.cell)) +
#   geom_point() +
#   theme_bw()
# ggsave('avp.oxt/avp/PCA.2.3.avp.oxt.neurons.all.png',
#        width = 10,
#        height = 10)
# 
# 
# #### all neuropeptides neurons ####
# ### avp, oxt, gnrh, sst, neurons
# ### subset data to relevant cells
# burtoni.snseq.combined.sct.all.neuropeptides.neurons = burtoni.snseq.combined.sct.all
# 
# #set idents
# Idents(object = burtoni.snseq.combined.sct.all.neuropeptides.neurons) <- "Cell.type"
# 
# #subset neurons
# burtoni.snseq.combined.sct.all.neuropeptides.neurons = subset(burtoni.snseq.combined.sct.all.neuropeptides.neurons,
#                                                               idents = c("neurons",
#                                                                          "excitatory",
#                                                                          "inhibitory"))
# 
# #subset expression
# burtoni.snseq.combined.sct.all.neuropeptides.neurons = subset(burtoni.snseq.combined.sct.all.neuropeptides.neurons,
#                                                               subset = avp > 4 | oxt > 4 | ENSONIG00000011023 > 4 | ENSONIG00000033642 > 4)
# 
# 
# ## run PCA, UMAP, and cluster 
# #use 0.4 resolution
# burtoni.snseq.combined.sct.all.neuropeptides.neurons.recluster = burtoni.snseq.combined.sct.all.neuropeptides.neurons %>% 
#   RunPCA() %>%
#   FindNeighbors(dims = 1:15) %>%
#   RunUMAP(dims = 1:15) %>%
#   FindClusters(resolution = 0.4)
# 
# ### clustree
# # cluster across resolutions
# burtoni.snseq.combined.sct.all.neuropeptides.neurons.clustree <- Seurat::FindClusters(object = burtoni.snseq.combined.sct.all.neuropeptides.neurons.recluster, 
#                                                                                       resolution = resolution.range)
# #check data
# head(burtoni.snseq.combined.sct.all.neuropeptides.neurons.clustree[[]])
# 
# #clustree
# clustree(burtoni.snseq.combined.sct.all.neuropeptides.neurons.clustree, 
#          prefix = "integrated_snn_res.",
#          node_colour = 'cluster',
#          node_size_range = c(10,20),
#          scale_node_text = TRUE) +
#   scale_edge_color_continuous(low = "black", 
#                               high = "black") +
#   theme(legend.position = "bottom")
# ggsave('neuropeptides/neuropeptides.neurons.clustree.png',
#        width = 10,
#        height = 10)
# 
# 
# ### graph 
# ##dom vs sub
# DimPlot(burtoni.snseq.combined.sct.all.neuropeptides.neurons.recluster, 
#         reduction = "umap", 
#         group.by = "orig.ident")
# ggsave('neuropeptides/DomVsSub.dimplot.neuropeptides.neurons.all.png')
# 
# ## clusters
# DimPlot(burtoni.snseq.combined.sct.all.neuropeptides.neurons.recluster, 
#         reduction = "umap", 
#         label = TRUE,
#         repel = TRUE)
# ggsave('neuropeptides/Clusters.dimplot.neuropeptides.neurons.all.png',
#        width = 10,
#        height = 10)
# 
# ##expression level
# burtoni.snseq.combined.sct.all.neuropeptides.neurons.recluster.expression = full_join(full_join(burtoni.snseq.combined.sct.all.neuropeptides.neurons.recluster@reductions$umap@cell.embeddings %>% 
#                                                                                                   as.data.frame() %>% 
#                                                                                                   rownames_to_column("Cell.id"),
#                                                                                                 burtoni.snseq.combined.sct.all.neuropeptides.neurons.recluster@meta.data %>% 
#                                                                                                   rownames_to_column("Cell.id")),
#                                                                                       burtoni.snseq.combined.sct.all.neuropeptides.neurons.recluster@assays$integrated@scale.data %>% 
#                                                                                         as.data.frame() %>% 
#                                                                                         filter(rownames(burtoni.snseq.combined.sct.all.neuropeptides.neurons.recluster@assays$integrated@scale.data) %in% c('avp',
#                                                                                                                                                                                                             'oxt',
#                                                                                                                                                                                                             'ENSONIG00000011023',
#                                                                                                                                                                                                             'ENSONIG00000033642')) %>% 
#                                                                                         t() %>% as.data.frame() %>% 
#                                                                                         rownames_to_column('Cell.id')) %>% 
#   full_join(burtoni.snseq.combined.sct.all.neuropeptides.neurons.recluster@reductions$pca@cell.embeddings %>% 
#               as.data.frame() %>% 
#               select(PC_1,
#                      PC_2,
#                      PC_3,
#                      PC_4,
#                      PC_5) %>% 
#               rownames_to_column("Cell.id")) %>% 
#   mutate(Cell.type.neuropeptides = case_when(avp > 4 ~ 'avp',
#                                              oxt > 4 ~ 'oxt',
#                                              ENSONIG00000011023 > 3.9 ~ 'gnrh',
#                                              ENSONIG00000033642 > 4 ~ 'sst',
#                                              TRUE ~ "both"))
# 
# 
# ##graph
# burtoni.snseq.combined.sct.all.neuropeptides.neurons.recluster.expression %>% 
#   ggplot(aes(x = UMAP_1,
#              y = UMAP_2,
#              color = integrated_snn_res.0.4,
#              shape = Cell.type.neuropeptides)) +
#   geom_point() +
#   theme_bw()
# ggsave('neuropeptides/Umap.neuropeptides.neurons.all.png',
#        width = 10,
#        height = 10)
# 
# #poster
# burtoni.snseq.combined.sct.all.neuropeptides.neurons.recluster.expression %>% 
#   ggplot(aes(x = UMAP_1,
#              y = UMAP_2,
#              color = Cell.type.neuropeptides)) +
#   geom_point(size = 3) +
#   theme_classic() +
#   scale_color_manual(values = c('orange3',
#                                 "purple",
#                                 "21B9CA",
#                                 "yellow3"))+ 
#   labs(color='Neuropeptides') +
#   theme(legend.text=element_text(size=20),
#         legend.title=element_text(size=20),
#         legend.position = c(0.8, 0.2)) +
#   guides(colour = guide_legend(override.aes = list(size=10)))
# ggsave('neuropeptides/Umap.neuropeptides.neurons.all.celltype.pdf',
#        width = 12,
#        height = 10)
# 
# ## cell count per cluster
# table(burtoni.snseq.combined.sct.all.neuropeptides.neurons.recluster.expression %>% 
#         select(integrated_snn_res.0.4,
#                Cell.type.neuropeptides))
# 
# table(burtoni.snseq.combined.sct.all.neuropeptides.neurons.recluster.expression %>% 
#         select(integrated_snn_res.0.4,
#                Cell.type.neuropeptides,
#                orig.ident))
# 
# table(burtoni.snseq.combined.sct.all.neuropeptides.neurons.recluster.expression %>% 
#         select(integrated_snn_res.0.4,
#                orig.ident))
# 
# table(burtoni.snseq.combined.sct.all.neuropeptides.neurons.recluster.expression %>% 
#         select(Cell.type.neuropeptides,
#                orig.ident))
# 
# ##chi squared test
# chisq.test.nueropeptides.neurons = chisq.test(table(burtoni.snseq.combined.sct.all.neuropeptides.neurons.recluster.expression %>% 
#                                                       select(Cell.type.neuropeptides,
#                                                              orig.ident)) %>% 
#                                                 t())
# #X-squared = 26.909, df = 3, p-value = 6.152e-06
# 
# ## test if 50% for each
# chisq.test(table(burtoni.snseq.combined.sct.all.neuropeptides.neurons.recluster.expression %>% 
#                    select(Cell.type.neuropeptides,
#                           orig.ident) %>% 
#                    filter(Cell.type.neuropeptides == 'oxt')),
#            p =c(0.5,0.5))
# 
# 
# 
# ## graph
# 
# table(burtoni.snseq.combined.sct.all.neuropeptides.neurons.recluster.expression %>% 
#         select(Cell.type.neuropeptides,
#                orig.ident)) %>% 
#   as_tibble() %>% 
#   ggplot(aes(x = Cell.type.neuropeptides,
#              y = n,
#              group = orig.ident,
#              fill = orig.ident)) +
#   geom_bar(stat = 'identity',
#            position=position_dodge(width=0.8),
#            width = 0.7) +
#   theme_classic() +
#   scale_fill_manual(values = c('orange3',
#                                "21B9CA"))
# ggsave('neuropeptides/Bargraph.neuropeptides.neurons.all.celltype.count.pdf',
#        width = 5,
#        height = 5)
# 
# ### PCA
# #1 vs 2
# burtoni.snseq.combined.sct.all.neuropeptides.neurons.recluster.expression %>% 
#   ggplot(aes(x = PC_1,
#              y = PC_2,
#              color = integrated_snn_res.0.4,
#              shape = Cell.type.neuropeptides)) +
#   geom_point() +
#   theme_bw()
# ggsave('neuropeptides/PCA.1.2.neuropeptides.neurons.all.png',
#        width = 10,
#        height = 10)
# #2 vs 3
# burtoni.snseq.combined.sct.all.neuropeptides.neurons.recluster.expression %>% 
#   ggplot(aes(x = PC_3,
#              y = PC_2,
#              shape = integrated_snn_res.0.4,
#              color = Cell.type.neuropeptides)) +
#   geom_point() +
#   theme_bw()
# ggsave('neuropeptides/PCA.2.3.expression.neuropeptides.neurons.all.png',
#        width = 10,
#        height = 10)
# 
# 
# 
# #### compare neuropeptide membership ####
# ### avp.oxt vs avp.oxt.neurons
# png('avp.oxt/avp/heatmap avp.oxt vs avp.oxt.neurons.png',
#     width = 10,
#     height = 10,
#     units = 'in',
#     res = 300)
# table(left_join(burtoni.snseq.combined.sct.all.avp.oxt.neurons.recluster.expression %>% 
#                   select(Cell.id,
#                          integrated_snn_res.0.8,
#                          Oxt.cell) %>% 
#                   rename(Avp.oxt.neuron.cluster = integrated_snn_res.0.8),
#                 burtoni.snseq.combined.sct.all.avp.oxt.recluster.expression %>% 
#                   select(Cell.id,
#                          integrated_snn_res.1) %>% 
#                   rename(Avp.oxt.cluster = integrated_snn_res.1)) %>% 
#         filter(Oxt.cell != 'both') %>% 
#         select(Avp.oxt.neuron.cluster,
#                Avp.oxt.cluster,
#                Oxt.cell))  %>% 
#   unclass() %>% 
#   as.data.frame() %>% 
#   as.matrix() %>% 
#   pheatmap(main = 'avp.oxt vs avp.oxt.neurons',
#            color = colorRampPalette(c("white", "red"))(50))
# dev.off()
# 
# ### avp.neurons vs avp.oxt.neurons
# png('avp.oxt/avp/heatmap avp.neurons vs avp.oxt.neurons.png',
#     width = 10,
#     height = 10,
#     units = 'in',
#     res = 300)
# table(right_join(burtoni.snseq.combined.sct.all.avp.oxt.neurons.recluster.expression %>% 
#                    select(Cell.id,
#                           integrated_snn_res.0.8,
#                           Oxt.cell) %>% 
#                    rename(Avp.oxt.neuron.cluster = integrated_snn_res.0.8),
#                  burtoni.snseq.combined.sct.all.avp.neurons.recluster.expression %>% 
#                    select(Cell.id,
#                           integrated_snn_res.1) %>% 
#                    rename(Avp.neuron.cluster = integrated_snn_res.1)) %>% 
#         filter(Oxt.cell != 'both') %>% 
#         select(Avp.oxt.neuron.cluster,
#                Avp.neuron.cluster,
#                Oxt.cell))  %>% 
#   unclass() %>% 
#   as.data.frame() %>% 
#   as.matrix() %>% 
#   pheatmap(main = 'avp.neurons vs avp.oxt.neurons',
#            color = colorRampPalette(c("white", "red"))(50))
# dev.off()
# 
# #better color pallette
# png('avp.oxt/avp/heatmap avp.neurons vs avp.oxt.neurons orange.png',
#     width = 10,
#     height = 10,
#     units = 'in',
#     res = 300)
# table(right_join(burtoni.snseq.combined.sct.all.avp.oxt.neurons.recluster.expression %>% 
#                    select(Cell.id,
#                           integrated_snn_res.0.8,
#                           Oxt.cell) %>% 
#                    rename(Avp.oxt.neuron.cluster = integrated_snn_res.0.8),
#                  burtoni.snseq.combined.sct.all.avp.neurons.recluster.expression %>% 
#                    select(Cell.id,
#                           integrated_snn_res.1) %>% 
#                    rename(Avp.neuron.cluster = integrated_snn_res.1)) %>% 
#         filter(Oxt.cell != 'both') %>% 
#         select(Avp.oxt.neuron.cluster,
#                Avp.neuron.cluster,
#                Oxt.cell))  %>% 
#   unclass() %>% 
#   as.data.frame() %>% 
#   as.matrix() %>% 
#   pheatmap(main = 'avp.neurons vs avp.oxt.neurons',
#            color = colorRampPalette(c("white", "orange3"))(30))
# dev.off()
# 
# 
# ### avp.neurons vs avp.oxt
# png('avp.oxt/avp/heatmap avp.oxt vs avp.neurons.png',
#     width = 10,
#     height = 10,
#     units = 'in',
#     res = 300)
# table(right_join(burtoni.snseq.combined.sct.all.avp.oxt.recluster.expression %>% 
#                    select(Cell.id,
#                           integrated_snn_res.1) %>% 
#                    rename(Avp.oxt.cluster = integrated_snn_res.1),
#                  burtoni.snseq.combined.sct.all.avp.neurons.recluster.expression %>% 
#                    select(Cell.id,
#                           integrated_snn_res.1,
#                           Oxt.cell) %>% 
#                    rename(Avp.cluster.neuron = integrated_snn_res.1)) %>% 
#         filter(Oxt.cell != 'both') %>% 
#         select(Avp.oxt.cluster,
#                Avp.cluster.neuron,
#                Oxt.cell)) %>% 
#   unclass() %>% 
#   as.data.frame() %>% 
#   as.matrix() %>% 
#   pheatmap(main = 'avp.oxt vs avp.neurons',
#            color = colorRampPalette(c("white", "red"))(50))
# dev.off()
# 
# ### avp.neurons vs neuropeptides all
# png('avp.oxt/avp/heatmap all neuropeptides vs avp.neurons.png',
#     width = 5,
#     height = 5,
#     units = 'in',
#     res = 720)
# table(right_join(burtoni.snseq.combined.sct.all.neuropeptides.neurons.recluster.expression %>% 
#                    select(Cell.id,
#                           integrated_snn_res.0.4) %>% 
#                    rename(neuropeptides.cluster = integrated_snn_res.0.4),
#                  burtoni.snseq.combined.sct.all.avp.neurons.recluster.expression %>% 
#                    select(Cell.id,
#                           integrated_snn_res.1,
#                           Oxt.cell) %>% 
#                    rename(Avp.cluster.neuron = integrated_snn_res.1)) %>% 
#         select(neuropeptides.cluster,
#                Avp.cluster.neuron)) %>% 
#   unclass() %>% 
#   as.data.frame() %>% 
#   as.matrix() %>% 
#   pheatmap(color = colorRampPalette(c("white", "orange3"))(50),
#            fontsize = 20,
#            cex = 1,
#            legend = FALSE)
# dev.off()
# 
# ### avp.neurons across threshold
# png('avp.oxt/avp/heatmap avp.neurons across threshold.png',
#     width = 10,
#     height = 10,
#     units = 'in',
#     res = 300)
# table(left_join(burtoni.snseq.combined.sct.all.avp.neurons.4.recluster.expression %>% 
#                   select(Cell.id,
#                          integrated_snn_res.1) %>% 
#                   rename(Avp.neuron.cluster.4 = integrated_snn_res.1),
#                 burtoni.snseq.combined.sct.all.avp.neurons.recluster.expression %>% 
#                   select(Cell.id,
#                          integrated_snn_res.0.8) %>% 
#                   rename(Avp.neuron.cluster.2 = integrated_snn_res.0.8)) %>%  
#         select(Avp.neuron.cluster.4,
#                Avp.neuron.cluster.2))  %>% 
#   unclass() %>% 
#   as.data.frame() %>% 
#   as.matrix() %>% 
#   pheatmap(main = 'avp.neurons across threshold',
#            color = colorRampPalette(c("white", "red"))(50))
# dev.off()
# 
# ### all neurons vs neuropeptides all
# png('avp.oxt/avp/heatmap all neurons vs avp.neurons.png',
#     width = 5,
#     height = 5,
#     units = 'in',
#     res = 720)
# table(right_join(burtoni.snseq.combined.sct.all.neurons.recluster.expression %>% 
#                    select(Cell.id,
#                           integrated_snn_res.0.8) %>% 
#                    rename(neuron.cluster = integrated_snn_res.0.8),
#                  burtoni.snseq.combined.sct.all.avp.neurons.recluster.expression %>% 
#                    select(Cell.id,
#                           integrated_snn_res.0.8) %>% 
#                    rename(Avp.cluster.neuron = integrated_snn_res.0.8)) %>% 
#         select(neuron.cluster,
#                Avp.cluster.neuron)) %>% 
#   unclass() %>% 
#   as.data.frame() %>% 
#   as.matrix() %>% 
#   pheatmap(color = colorRampPalette(c("white", "orange3"))(50),
#            fontsize = 20,
#            cex = 1,
#            display_numbers = T,
#            number_format = "%.0f") 
# dev.off()
# 

#### neurons ####
## use new dataframe
# burtoni.snseq.combined.sct.all.neurons = burtoni.snseq.combined.sct.all
burtoni.snseq.combined.sct.all.neurons = burtoni.snseq.combined.sct

# remove
# rm(burtoni.snseq.combined.sct.all)

#set idents
Idents(object = burtoni.snseq.combined.sct.all.neurons) <- "sctypemarkers.hypo"

#subset to neurons
burtoni.snseq.combined.sct.all.neurons = subset(burtoni.snseq.combined.sct.all.neurons,
                                            idents = c("C7-1: GLU",
                                                       "C7-2: GABA"))

# #set idents
# Idents(object = burtoni.snseq.combined.sct.all.neurons) <- "Cell.type"
# 
# #subset to neurons
# burtoni.snseq.combined.sct.all.neurons = subset(burtoni.snseq.combined.sct.all.neurons,
#                                                 idents = c("neurons",
#                                                            "excitatory",
#                                                            "inhibitory"))

# set assay integrated for clustering
DefaultAssay(burtoni.snseq.combined.sct.all.neurons) = 'integrated'

## run PCA, UMAP, and cluster 
#use 0.4 resolution
burtoni.snseq.combined.sct.all.neurons.recluster = burtoni.snseq.combined.sct.all.neurons %>% 
  RunPCA() %>%
  FindNeighbors(dims = 1:15) %>%
  RunUMAP(dims = 1:15) %>%
  FindClusters(resolution = 0.8)

# remove
# rm(burtoni.snseq.combined.sct.all.neurons)

# 
# ### clustree
# # cluster across resolutions
# burtoni.snseq.combined.sct.all.neurons.clustree <- Seurat::FindClusters(object = burtoni.snseq.combined.sct.all.neurons.recluster, 
#                                                                         resolution = resolution.range.reduced.2)
# 
# #clustree
# clustree(burtoni.snseq.combined.sct.all.neurons.clustree, 
#          prefix = "integrated_snn_res.",
#          scale_node_text = TRUE) +
#   scale_edge_color_continuous(low = "black", 
#                               high = "black") +
#   theme(legend.position = "bottom")
# ggsave('neurons/neurons.clustree.png',
#        width = 10,
#        height = 10)
# 
# 
# ### graph 
# # idents to new clusters
# Idents(object = burtoni.snseq.combined.sct.all.neurons.recluster) <- "integrated_snn_res.0.8"
# 
# ## hover locator
# DimPlot(burtoni.snseq.combined.sct.all.neurons.recluster, 
#         reduction = "umap", 
#         label = TRUE,
#         repel = TRUE, 
#         group.by = "SCT_snn_res.0.8") %>% 
#   HoverLocator(information = FetchData(object = burtoni.snseq.combined.sct.all.neurons.recluster, 
#                                        vars = 'SCT_snn_res.0.8'))
# 
# ##dom vs sub
# DimPlot(burtoni.snseq.combined.sct.all.neurons.recluster, 
#         reduction = "umap", 
#         group.by = "orig.ident")
# ggsave('neurons/DomVsSub.dimplot.neurons.all.png')
# 
# ## clusters
# DimPlot(burtoni.snseq.combined.sct.all.neurons.recluster, 
#         reduction = "umap", 
#         label = TRUE,
#         repel = TRUE, 
#         group.by = "integrated_snn_res.0.8")
# ggsave('neurons/Clusters.dimplot.avp.neurons.all.png',
#        width = 10,
#        height = 10)
# 
# 
# DimPlot(burtoni.snseq.combined.sct.all.neurons.recluster, 
#         reduction = "umap", 
#         label = TRUE,
#         repel = TRUE, 
#         group.by = "Cell.type")



##expression level
burtoni.snseq.combined.sct.all.neurons.recluster.expression = full_join(burtoni.snseq.combined.sct.all.neurons.recluster@reductions$umap@cell.embeddings %>% 
                                                                                    as.data.frame() %>% 
                                                                                    rownames_to_column("Cell.id"),
                                                                                  burtoni.snseq.combined.sct.all.neurons.recluster@meta.data %>% 
                                                                                    rownames_to_column("Cell.id")) %>% full_join(burtoni.snseq.combined.sct.all.neurons.recluster@reductions$pca@cell.embeddings %>% 
              as.data.frame() %>% 
              select(PC_1,
                     PC_2,
                     PC_3,
                     PC_4,
                     PC_5) %>% 
              rownames_to_column("Cell.id")) 






#### neuropeptides comparison ####
### create matrix of neuropeptide expression level 
## use neuropeptides.list
neuropeptides.list.gene.names =  neuropeptides.list %>% 
  filter(!is.na(Gene.name.nile.tilapia)) %>%
  pull(Gene.name.nile.tilapia) %>% 
  unique()

## use all neurons expression recluster
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides = burtoni.snseq.combined.sct.all.neurons.recluster@assays$SCT@data %>% 
  as.data.frame() %>% 
  filter(rownames(burtoni.snseq.combined.sct.all.neurons.recluster@assays$SCT@data) %in% neuropeptides.list.gene.names) %>% 
  t() %>% 
  as.data.frame()

## reduce columns to neuropeptides
# create presence absence matrix
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix = burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides %>% 
  select(any_of(neuropeptides.list.gene.names)) %>% 
  mutate(across(where(is.numeric), 
                function(x) ifelse(x < 1, 0, x))) %>% 
  mutate(across(where(is.numeric), 
                function(x) ifelse(x >= 1, 1, x))) 

###rename to actual gene name
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix = burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix %>% 
  dplyr::rename(#ADM = ENSONIG00000036130,
    CHGA = ENSONIG00000000725,
    GHRH = ENSONIG00000020252,
    GNRH1 = ENSONIG00000011023,
    #KNG1 = ENSONIG00000008378,
    NMB = ENSONIG00000002537,
    NPFF = ENSONIG00000005229,
    SST = ENSONIG00000033642)

## remove genes with no expression
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix = burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix %>% 
  select_if(colSums(.) != 0)

#how many neurons have neuropeptides
#total 14592
#~77% have at least one
#~47% have at least two
# 
# #graph number of neuropeptides per cell
# burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix %>%
#   mutate(Total.neuropeptides = rowSums(.)) %>% 
#   ggplot(aes(Total.neuropeptides)) +
#   geom_histogram(binwidth = 0.3) +
#   theme_classic() +
#   scale_x_continuous(breaks=seq(0,11,1)) +
#   xlab('Neuropeptides per neuron') +
#   theme(axis.text = element_text(size = 20),
#         axis.title = element_text(size = 20),
#         title =element_text(size=20, face='bold'))
# ggsave('neuropeptides/comparison/Neuropeptides per neuron.png',
#        width = 10,
#        height = 10)
# 
# ## graph gene per number of neuropeptide cell
# burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix %>%
#   mutate(Total.neuropeptides = rowSums(.)) %>% 
#   pivot_longer(cols = -c(Total.neuropeptides),
#                names_to = "gene",
#               values_to = "present") %>% 
#   filter(Total.neuropeptides > 0) %>% 
#   group_by(Total.neuropeptides,
#            gene) %>% 
#   summarise(Present = sum(present)) %>% 
#   group_by(gene) %>% 
#   mutate(Total.count = sum(Present)) %>% 
#   ungroup() %>% 
#   mutate(Percentage = 100*Present/Total.count,
#          Percentage = ifelse(is.na(Percentage),
#                              0,
#                              Percentage)) %>% 
#   ggplot(aes(x = Total.neuropeptides,
#              y= reorder(gene, 
#                         Total.count))) +
#   geom_tile(aes(fill = Percentage)) +
#   scale_fill_gradientn(colours=c("white",
#                                  "black"), 
#                        limits = c(0,100)) +
#   theme_classic() +
#   ylab('Celltype') +
#   xlab('Neuropeptides per cell') +
#   ggtitle('Percentage of neuropeptide cells expressing multiple neuropeptides')
# ggsave('neuropeptides/comparison/Neuropeptides per neuron celltype.png',
#        width = 10,
#        height = 10)
# 
# #create ridge lines
# burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix %>%
#   mutate(Total.neuropeptides = rowSums(.)) %>% 
#   pivot_longer(cols = -c(Total.neuropeptides),
#                names_to = "gene",
#                values_to = "present") %>% 
#   filter(Total.neuropeptides > 0) %>% 
#   group_by(Total.neuropeptides,
#            gene) %>% 
#   summarise(Present = sum(present)) %>% 
#   group_by(gene) %>% 
#   mutate(Total.count = sum(Present)) %>% 
#   ungroup() %>% 
#   mutate(Percentage = 100*Present/Total.count) %>% 
#   ggplot(aes(x = Total.neuropeptides,
#              y= reorder(gene, 
#                         Total.count),
#              height = Percentage)) +
#   geom_density_ridges(stat = "identity",
#                       scale = 10,
#                       alpha=0.75) +
#   theme_ridges() +
#   ylab('Celltype') +
#   xlab('Neuropeptides per cell') +
#   ggtitle('Percentage of neuropeptide cells expressing multiple neuropeptides')
# ggsave('neuropeptides/comparison/Neuropeptides per neuron celltype ridgeline.png',
#        width = 10,
#        height = 10)

## calculate the overlap 
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap <- crossprod(as.matrix(burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix))

##calculate total count
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.count = diag(burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap) %>% 
  as.data.frame()
#add colname
colnames(burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.count) = 'Total.count'
#add gene names
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.count = burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.count %>% 
  rownames_to_column(var = 'Gene.name')

# ##graph total count
# burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.count %>%
#   ggplot(aes(Total.count)) +
#   geom_histogram(binwidth = 100) +
#   theme_classic()
# ggsave('neuropeptides/comparison/Neuropeptides neuron counts.png',
#        width = 10,
#        height = 10)

# calculate the percentage
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent = t(burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap * 100 / diag(burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap))     
#convert data to long format
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent = burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent %>% 
  as.data.frame.table() 

# ##graph
# burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent %>% 
#   filter(Freq < 100) %>% 
#   ggplot(aes(x=Var1, 
#              y=Var2)) +
#   geom_tile(aes(fill=Freq)) +
#   theme(axis.text.x = element_text(angle = 90)) +
#   scale_fill_gradientn(colours=c("blue","red")) 
# ggsave('neuropeptides/comparison/Neuropeptides percentage overlap.png',
#        width = 10,
#        height = 10)

#add gene count data
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder = left_join(burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent,
                                                                                                                     burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.count,
                                                                                                                     by = c('Var1' = 'Gene.name')) %>% 
  dplyr::rename(Var1.total.count = Total.count) %>% 
  left_join(burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.count,
            by = c('Var2' = 'Gene.name')) %>% 
  dplyr::rename(Var2.total.count = Total.count)

#reorder 
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder <- burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder %>%
  arrange(Var1.total.count,
          Var2.total.count) %>%
  mutate(Var1 = fct_inorder(factor(Var1, ordered=TRUE)),
         Var2 = fct_inorder(factor(Var2, ordered=TRUE)))

# #graph reorder
# burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder %>% 
#   filter(Freq < 100) %>% 
#   ggplot(aes(x=Var1, 
#              y=Var2)) +
#   geom_tile(aes(fill=Freq)) +
#   theme(axis.text.x = element_text(angle = 90)) +
#   scale_fill_gradientn(colours=c("blue","red")) 
# ggsave('neuropeptides/comparison/Neuropeptides percentage overlap reorder.png',
#        width = 10,
#        height = 10)
# 
# #graph reorder
# #filter cell count
# burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder %>% 
#   filter(Freq < 100) %>% 
#   filter(Var2.total.count > 250) %>% 
#   filter(Var1.total.count > 250) %>% 
#   ggplot(aes(x=Var1, 
#              y=Var2)) +
#   geom_tile(aes(fill=Freq)) +
#   theme(axis.text.x = element_text(angle = 90)) +
#   scale_fill_gradientn(colours=c("blue","red")) 
# ggsave('neuropeptides/comparison/Neuropeptides percentage overlap reorder filter.png',
#        width = 10,
#        height = 10)



# calculate the percentage
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.df = burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap %>% 
  as.data.frame.table() 

# ##graph
# burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.df %>% 
#   ggplot(aes(x=Var1, 
#              y=Var2)) +
#   geom_tile(aes(fill=Freq)) +
#   theme(axis.text.x = element_text(angle = 90)) +
#   scale_fill_gradientn(colours=c("blue","red")) 
# ggsave('neuropeptides/comparison/Neuropeptides overlap.png',
#        width = 10,
#        height = 10)

#add gene count data
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.df.reorder = left_join(burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.df,
                                                                                                                burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.count,
                                                                                                                by = c('Var1' = 'Gene.name')) %>% 
  dplyr::rename(Var1.total.count = Total.count) %>% 
  left_join(burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.count,
            by = c('Var2' = 'Gene.name')) %>% 
  dplyr::rename(Var2.total.count = Total.count)

#reorder 
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.df.reorder <- burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.df.reorder %>%
  arrange(Var1.total.count,
          Var2.total.count) %>%
  mutate(Var1 = fct_inorder(factor(Var1, ordered=TRUE)),
         Var2 = fct_inorder(factor(Var2, ordered=TRUE)))

# remove individual overlap
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.df.reorder = burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.df.reorder %>% 
  mutate(Freq.na = ifelse(Var1.total.count == Var2.total.count,
                          NA,
                          Freq))

# #graph reorder
# burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.df.reorder %>% 
#   ggplot(aes(x=Var1, 
#              y=Var2)) +
#   geom_tile(aes(fill=Freq.na)) +
#   theme(axis.text.x = element_text(angle = 90)) +
#   scale_fill_gradientn(colours=c("blue","red")) 
# ggsave('neuropeptides/comparison/Neuropeptides overlap rescaled.png',
#        width = 10,
#        height = 10)






#### neuropeptides comparison across social status ####
## set idents
Idents(object = burtoni.snseq.combined.sct.all.neurons.recluster) = "orig.ident"
### create matrix of neuropeptide expression level 
## use neuropeptides.list
## reduce columns to neuropeptides
## seperate into dom and sub
# create presence absence matrix
#dom
dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix = GetAssayData(subset(burtoni.snseq.combined.sct.all.neurons.recluster,
                    idents = 'dom_burtoni_snseq'),
             assay = 'SCT',
             slot = 'data') %>% 
  as.data.frame() %>% 
  filter(rownames(burtoni.snseq.combined.sct.all.neurons.recluster@assays$SCT@data) %in% neuropeptides.list.gene.names) %>% 
  t() %>% 
  as.data.frame() %>% 
  select(any_of(neuropeptides.list.gene.names)) %>% 
  mutate(across(where(is.numeric), 
                function(x) ifelse(x < 1, 0, x))) %>% 
  mutate(across(where(is.numeric), 
                function(x) ifelse(x >= 1, 1, x))) 

#sub
sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix = GetAssayData(subset(burtoni.snseq.combined.sct.all.neurons.recluster,
                                                                                                           idents = 'sub_burtoni_snseq'),
                                                                                                    assay = 'SCT',
                                                                                                    slot = 'data') %>% 
  as.data.frame() %>% 
  filter(rownames(burtoni.snseq.combined.sct.all.neurons.recluster@assays$SCT@data) %in% neuropeptides.list.gene.names) %>% 
  t() %>% 
  as.data.frame() %>% 
  select(any_of(neuropeptides.list.gene.names)) %>% 
  mutate(across(where(is.numeric), 
                function(x) ifelse(x < 1, 0, x))) %>% 
  mutate(across(where(is.numeric), 
                function(x) ifelse(x >= 1, 1, x))) 

##rename to actual gene name
#dom
dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix = dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix %>% 
  dplyr::rename(#ADM = ENSONIG00000036130,
    CHGA = ENSONIG00000000725,
    GHRH = ENSONIG00000020252,
    GNRH1 = ENSONIG00000011023,
    #KNG1 = ENSONIG00000008378,
    NMB = ENSONIG00000002537,
    NPFF = ENSONIG00000005229,
    SST = ENSONIG00000033642)

#sub
sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix = sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix %>% 
  dplyr::rename(#ADM = ENSONIG00000036130,
    CHGA = ENSONIG00000000725,
    GHRH = ENSONIG00000020252,
    GNRH1 = ENSONIG00000011023,
    #KNG1 = ENSONIG00000008378,
    NMB = ENSONIG00000002537,
    NPFF = ENSONIG00000005229,
    SST = ENSONIG00000033642)

#how many neurons have neuropeptides
#dom has 6454 neurons and 6061 expressing neuropeptide
#dom has 8893 neurons and 8328 expressing neuropeptide
dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix %>%
  mutate(Total.neuropeptides = rowSums(.)) %>% 
  filter(Total.neuropeptides != 0) %>% 
  nrow()
# total neurons
dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix %>% 
  nrow()


#sub has 6187 neurons and 3608 expressing neuropeptide
#sub has 7820 neurons and 4498 expressing neuropeptide
sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix %>%
  mutate(Total.neuropeptides = rowSums(.)) %>% 
  filter(Total.neuropeptides != 0) %>% 
  nrow()

# total neurons
sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix %>% 
  nrow()

#dom has 8893 neurons and 8328 expressing neuropeptide
#sub has 7820 neurons and 4498 expressing neuropeptide
# 
# #graph number of neuropeptides per cell
# #dom
# dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix %>%
#   mutate(Total.neuropeptides = rowSums(.)) %>% 
#   ggplot(aes(Total.neuropeptides)) +
#   geom_histogram(binwidth = 0.3) +
#   theme_classic() +
#   scale_x_continuous(breaks=seq(0,11,1)) +
#   ggtitle('Neuropeptides per neuron dom')
# ggsave('neuropeptides/comparison/social.status/Neuropeptides per neuron dom.png',
#        width = 10,
#        height = 10)
# 
# #sub
# sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix %>%
#   mutate(Total.neuropeptides = rowSums(.)) %>% 
#   ggplot(aes(Total.neuropeptides)) +
#   geom_histogram(binwidth = 0.3) +
#   theme_classic() +
#   scale_x_continuous(breaks=seq(0,11,1)) +
#   ggtitle('Neuropeptides per neuron sub')
# ggsave('neuropeptides/comparison/social.status/Neuropeptides per neuron sub.png',
#        width = 10,
#        height = 10)
# 
# #graph together
# dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix %>%
#   mutate(Total.neuropeptides = rowSums(.)) %>%
#   select(Total.neuropeptides) %>%
#   mutate(Status = 'dom') %>%
#   rbind(sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix %>%
#           mutate(Total.neuropeptides = rowSums(.)) %>%
#           select(Total.neuropeptides) %>%
#           mutate(Status = 'sub')) %>%
#   mutate(count = 1) %>%
#   group_by(Status,
#            Total.neuropeptides) %>%
#   summarise(count.total = sum(count)) %>%
#   group_by(Status) %>%
#   mutate(Total = sum(count.total)) %>%
#   ungroup() %>%
#   mutate(Percent = 100*count.total/Total) %>% 
#   ggplot(aes(y = Percent,
#              x = Total.neuropeptides,
#              group = Status,
#              color = Status)) +
#   geom_line(size = 4) +
#   geom_point(size = 10) +
#   theme_classic() +
#   scale_x_continuous(breaks=seq(0,11,1)) +
#   theme(axis.text = element_text(size = 20),
#          axis.title = element_text(size = 20),
#         title =element_text(size=20, face='bold')) +
#   scale_color_manual(values = c("#4e499e",
#                                 "#60bb46"))+
#   xlab('Neuropeptides per neuron')
# ggsave('neuropeptides/comparison_new/social.status/Neuropeptides per neuron across status.png',
#        width = 10,
#        height = 10)

# 
# 
# ## graph gene per number of neuropeptide cell
# #dom
# dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix %>%
#   select_if(colSums(.) != 0) %>% 
#   mutate(Total.neuropeptides = rowSums(.)) %>% 
#   pivot_longer(cols = -c(Total.neuropeptides),
#                names_to = "gene",
#                values_to = "present") %>% 
#   filter(Total.neuropeptides > 0) %>% 
#   group_by(Total.neuropeptides,
#            gene) %>% 
#   summarise(Present = sum(present)) %>% 
#   group_by(gene) %>% 
#   mutate(Total.count = sum(Present)) %>% 
#   ungroup() %>% 
#   mutate(Percentage = 100*Present/Total.count) %>%
#   ggplot(aes(x = Total.neuropeptides,
#              y= reorder(gene, 
#                         Total.count))) +
#   geom_tile(aes(fill = Percentage)) +
#   scale_fill_gradientn(colours=c("white",
#                                  "black"), 
#                        limits = c(0, 100)) +
#   theme_classic() +
#   ylab('Celltype') +
#   xlab('Neuropeptides per cell') +
#   ggtitle('Percentage of neuropeptide cells expressing multiple neuropeptides dom')
# ggsave('neuropeptides/comparison/social.status/Neuropeptides per neuron celltype dom.png',
#        width = 10,
#        height = 10)
# 
# #create ridge lines
# dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix %>%
#   select_if(colSums(.) != 0) %>% 
#   mutate(Total.neuropeptides = rowSums(.)) %>% 
#   pivot_longer(cols = -c(Total.neuropeptides),
#                names_to = "gene",
#                values_to = "present") %>% 
#   filter(Total.neuropeptides > 0) %>% 
#   group_by(Total.neuropeptides,
#            gene) %>% 
#   summarise(Present = sum(present)) %>% 
#   group_by(gene) %>% 
#   mutate(Total.count = sum(Present)) %>% 
#   ungroup() %>% 
#   mutate(Percentage = 100*Present/Total.count) %>% 
#   ggplot(aes(x = Total.neuropeptides,
#              y= reorder(gene, 
#                         Total.count),
#              height = Percentage)) +
#   geom_density_ridges(stat = "identity",
#                  scale = 10,
#                  alpha=0.75) +
#   theme_ridges() +
#   ylab('Celltype') +
#   xlab('Neuropeptides per cell') +
#   ggtitle('Percentage of neuropeptide cells expressing multiple neuropeptides dom')
# ggsave('neuropeptides/comparison/social.status/Neuropeptides per neuron celltype dom ridgeline.png',
#        width = 10,
#        height = 10)
# 
# 
# 
# #sub
# sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix %>%
#   select_if(colSums(.) != 0) %>% 
#   mutate(Total.neuropeptides = rowSums(.)) %>% 
#   pivot_longer(cols = -c(Total.neuropeptides),
#                names_to = "gene",
#                values_to = "present") %>% 
#   filter(Total.neuropeptides > 0) %>% 
#   group_by(Total.neuropeptides,
#            gene) %>% 
#   summarise(Present = sum(present)) %>% 
#   group_by(gene) %>% 
#   mutate(Total.count = sum(Present)) %>% 
#   ungroup() %>% 
#   mutate(Percentage = 100*Present/Total.count) %>% 
#   ggplot(aes(x = Total.neuropeptides,
#              y= reorder(gene, 
#                         Total.count))) +
#   geom_tile(aes(fill = Percentage)) +
#   scale_fill_gradientn(colours=c("white",
#                                  "black"), 
#                        limits = c(0, 100)) +
#   theme_classic() +
#   ylab('Celltype') +
#   xlab('Neuropeptides per cell') +
#   ggtitle('Percentage of neuropeptide cells expressing multiple neuropeptides sub')
# ggsave('neuropeptides/comparison/social.status/Neuropeptides per neuron celltype sub.png',
#        width = 10,
#        height = 10)
# 
# #create ridge lines
# sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix %>%
#   select_if(colSums(.) != 0) %>% 
#   mutate(Total.neuropeptides = rowSums(.)) %>% 
#   pivot_longer(cols = -c(Total.neuropeptides),
#                names_to = "gene",
#                values_to = "present") %>% 
#   filter(Total.neuropeptides > 0) %>% 
#   group_by(Total.neuropeptides,
#            gene) %>% 
#   summarise(Present = sum(present)) %>% 
#   group_by(gene) %>% 
#   mutate(Total.count = sum(Present)) %>% 
#   ungroup() %>% 
#   mutate(Percentage = 100*Present/Total.count) %>% 
#   ggplot(aes(x = Total.neuropeptides,
#              y= reorder(gene, 
#                         Total.count),
#              height = Percentage)) +
#   geom_density_ridges(stat = "identity",
#                       scale = 10,
#                       alpha=0.75) +
#   theme_ridges() +
#   ylab('Celltype') +
#   xlab('Neuropeptides per cell') +
#   ggtitle('Percentage of neuropeptide cells expressing multiple neuropeptides sub')
# ggsave('neuropeptides/comparison/social.status/Neuropeptides per neuron celltype sub ridgeline.png',
#        width = 10,
#        height = 10)
# 
# #difference 
# sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix %>%
#   select_if(colSums(.) != 0) %>% 
#   mutate(Total.neuropeptides = rowSums(.)) %>% 
#   pivot_longer(cols = -c(Total.neuropeptides),
#                names_to = "gene",
#                values_to = "present") %>% 
#   filter(Total.neuropeptides > 0) %>% 
#   group_by(Total.neuropeptides,
#            gene) %>% 
#   summarise(Present = sum(present)) %>% 
#   group_by(gene) %>% 
#   mutate(Total.count = sum(Present)) %>% 
#   ungroup() %>% 
#   mutate(Percentage = 100*Present/Total.count) %>% 
#   full_join(dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix %>%
#               select_if(colSums(.) != 0) %>% 
#               mutate(Total.neuropeptides = rowSums(.)) %>% 
#               pivot_longer(cols = -c(Total.neuropeptides),
#                            names_to = "gene",
#                            values_to = "present") %>% 
#               filter(Total.neuropeptides > 0) %>% 
#               group_by(Total.neuropeptides,
#                        gene) %>% 
#               summarise(Present = sum(present)) %>% 
#               group_by(gene) %>% 
#               mutate(Total.count = sum(Present)) %>% 
#               ungroup() %>% 
#               mutate(Percentage = 100*Present/Total.count),
#             by = c('Total.neuropeptides',
#                    'gene'),
#             suffix = c(".sub", 
#                       ".dom")) %>% 
#   mutate(Percentage.diff = Percentage.dom - Percentage.sub,
#          Ratio = Total.count.dom/Total.count.sub) %>% 
#   drop_na() %>%
#   ggplot(aes(x = Total.neuropeptides,
#              y= reorder(gene, 
#                         Ratio,
#                         FUN = mean))) +
#   geom_tile(aes(fill = Percentage.diff)) +
#   scale_fill_gradientn(colours=c("dark blue",
#                                  "blue",
#                                  "white",
#                                  "red",
#                                  "dark red"), 
#                        limits = c(-100, 100)) +
#   theme_classic() +
#   ylab('Celltype (order by bias)') +
#   xlab('Neuropeptides per cell') +
#   ggtitle('Percentage of neuropeptide cells expressing multiple neuropeptides across social status')
# ggsave('neuropeptides/comparison/social.status/Neuropeptides per neuron celltype difference.png',
#        width = 10,
#        height = 10)
#   
# #difference 
# #filter
# sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix %>%
#   mutate(Total.neuropeptides = rowSums(.)) %>% 
#   pivot_longer(cols = -c(Total.neuropeptides),
#                names_to = "gene",
#                values_to = "present") %>% 
#   filter(Total.neuropeptides > 0) %>% 
#   group_by(Total.neuropeptides,
#            gene) %>% 
#   summarise(Present = sum(present)) %>% 
#   group_by(gene) %>% 
#   mutate(Total.count = sum(Present)) %>% 
#   ungroup() %>% 
#   mutate(Percentage = 100*Present/Total.count) %>% 
#   full_join(dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix %>%
#               mutate(Total.neuropeptides = rowSums(.)) %>% 
#               pivot_longer(cols = -c(Total.neuropeptides),
#                            names_to = "gene",
#                            values_to = "present") %>% 
#               filter(Total.neuropeptides > 0) %>% 
#               group_by(Total.neuropeptides,
#                        gene) %>% 
#               summarise(Present = sum(present)) %>% 
#               group_by(gene) %>% 
#               mutate(Total.count = sum(Present)) %>% 
#               ungroup() %>% 
#               mutate(Percentage = 100*Present/Total.count),
#             by = c('Total.neuropeptides',
#                    'gene'),
#             suffix = c(".sub", 
#                        ".dom")) %>% 
#   mutate(Percentage.diff = Percentage.dom - Percentage.sub,
#          Ratio = Total.count.dom/Total.count.sub,
#          Total.count = Total.count.dom + Total.count.sub) %>% 
#   drop_na() %>% 
#   filter(Total.count >= 200) %>% 
#   filter(Ratio >= 0.5,
#          Ratio <= 2) %>% 
#   ggplot(aes(x = Total.neuropeptides,
#              y= reorder(gene, 
#                         Total.count,
#                         FUN = mean))) +
#   geom_tile(aes(fill = Percentage.diff)) +
#   scale_fill_gradientn(colours=c("dark blue",
#                                  "blue",
#                                  "white",
#                                  "red",
#                                  "dark red"), 
#                        limits = c(-100, 100)) +
#   theme_classic() +
#   ylab('Celltype (order by total)') +
#   xlab('Neuropeptides per cell') +
#   ggtitle('Percentage of neuropeptide cells expressing multiple neuropeptides across social status filter')
# ggsave('neuropeptides/comparison/social.status/Neuropeptides per neuron celltype difference filter.png',
#        width = 10,
#        height = 10)
#   
  
## calculate the overlap 
#dom
dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap <- crossprod(as.matrix(dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix))
#sub
sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap <- crossprod(as.matrix(sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix))

##calculate total count
#dom
dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.count = diag(dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap) %>% 
  as.data.frame()
#sub
sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.count = diag(sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap) %>% 
  as.data.frame()

#add colname
#dom
colnames(dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.count) = 'Total.count'
#sub
colnames(sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.count) = 'Total.count'
#add gene names
#dom
dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.count = dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.count %>% 
  rownames_to_column(var = 'Gene.name')
#sub
sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.count = sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.count %>% 
  rownames_to_column(var = 'Gene.name')
# 
# ##graph total count
# #dom
# dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.count %>%
#   ggplot(aes(Total.count)) +
#   geom_histogram(binwidth = 100) +
#   theme_classic() +
#   ggtitle('Neuropeptides neuron counts dom')
# ggsave('neuropeptides/comparison/social.status/Neuropeptides neuron counts dom.png',
#        width = 10,
#        height = 10)
# #sub
# sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.count %>%
#   ggplot(aes(Total.count)) +
#   geom_histogram(binwidth = 100) +
#   theme_classic() +
#   ggtitle('Neuropeptides neuron counts sub')
# ggsave('neuropeptides/comparison/social.status/Neuropeptides neuron counts sub.png',
#        width = 10,
#        height = 10)

# calculate the percentage
#dom
dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent = t(dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap * 100 / diag(dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap))   
#sub
sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent = t(sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap * 100 / diag(sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap))   
#convert data to long format
#dom
dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent = dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent %>% 
  as.data.frame.table()
#sub
sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent = sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent %>% 
  as.data.frame.table()
# 
# ##graph
# #dom
# dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent %>% 
#   filter(Freq < 100) %>% 
#   ggplot(aes(x=Var1, 
#              y=Var2)) +
#   geom_tile(aes(fill=Freq)) +
#   theme(axis.text.x = element_text(angle = 90)) +
#   scale_fill_gradientn(colours=c("blue","red")) +
#   ggtitle('Neuropeptides percentage overlap dom')
# ggsave('neuropeptides/comparison/social.status/Neuropeptides percentage overlap dom.png',
#        width = 10,
#        height = 10)
# #sub
# sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent %>% 
#   filter(Freq < 100) %>% 
#   ggplot(aes(x=Var1, 
#              y=Var2)) +
#   geom_tile(aes(fill=Freq)) +
#   theme(axis.text.x = element_text(angle = 90)) +
#   scale_fill_gradientn(colours=c("blue","red")) +
#   ggtitle('Neuropeptides percentage overlap sub')
# ggsave('neuropeptides/comparison/social.status/Neuropeptides percentage overlap sub.png',
#        width = 10,
#        height = 10)

#add gene count data
#dom
dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder = left_join(dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent,
                                                                                                                         dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.count,
                                                                                                                         by = c('Var1' = 'Gene.name')) %>% 
  dplyr::rename(Var1.total.count = Total.count) %>% 
  left_join(dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.count,
            by = c('Var2' = 'Gene.name')) %>% 
  dplyr::rename(Var2.total.count = Total.count)

#sub
sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder = left_join(sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent,
                                                                                                                         sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.count,
                                                                                                                         by = c('Var1' = 'Gene.name')) %>% 
  dplyr::rename(Var1.total.count = Total.count) %>% 
  left_join(sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.count,
            by = c('Var2' = 'Gene.name')) %>% 
  dplyr::rename(Var2.total.count = Total.count)

#reorder 
#dom
dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder <- dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder %>%
  arrange(Var1.total.count,
          Var2.total.count) %>%
  mutate(Var1 = fct_inorder(factor(Var1, ordered=TRUE)),
         Var2 = fct_inorder(factor(Var2, ordered=TRUE)))
#sub
sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder <- sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder %>%
  arrange(Var1.total.count,
          Var2.total.count) %>%
  mutate(Var1 = fct_inorder(factor(Var1, ordered=TRUE)),
         Var2 = fct_inorder(factor(Var2, ordered=TRUE)))
# 
# #graph reorder
# #dom
# dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder %>% 
#   filter(Freq < 100) %>% 
#   ggplot(aes(x=Var1, 
#              y=Var2)) +
#   geom_tile(aes(fill=Freq)) +
#   theme(axis.text.x = element_text(angle = 90)) +
#   scale_fill_gradientn(colours=c("blue","red")) +
#   ggtitle('Neuropeptides percentage overlap reorder dom')
# ggsave('neuropeptides/comparison/social.status/Neuropeptides percentage overlap reorder dom.png',
#        width = 10,
#        height = 10)
# #sub
# sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder %>% 
#   filter(Freq < 100) %>% 
#   ggplot(aes(x=Var1, 
#              y=Var2)) +
#   geom_tile(aes(fill=Freq)) +
#   theme(axis.text.x = element_text(angle = 90)) +
#   scale_fill_gradientn(colours=c("blue","red")) +
#   ggtitle('Neuropeptides percentage overlap reorder sub')
# ggsave('neuropeptides/comparison/social.status/Neuropeptides percentage overlap reorder sub.png',
#        width = 10,
#        height = 10)
# 
# #graph reorder
# #filter cell count
# #dom
# dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder %>% 
#   filter(Freq < 100) %>% 
#   filter(Var2.total.count > 250) %>% 
#   filter(Var1.total.count > 250) %>% 
#   ggplot(aes(x=Var1, 
#              y=Var2)) +
#   geom_tile(aes(fill=Freq)) +
#   theme(axis.text.x = element_text(angle = 90)) +
#   scale_fill_gradientn(colours=c("blue","red")) +
#   ggtitle('Neuropeptides percentage overlap reorder filter dom')
# ggsave('neuropeptides/comparison/social.status/Neuropeptides percentage overlap reorder filter dom.png',
#        width = 10,
#        height = 10)
# #sub
# sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder %>% 
#   filter(Freq < 100) %>% 
#   filter(Var2.total.count > 250) %>% 
#   filter(Var1.total.count > 250) %>% 
#   ggplot(aes(x=Var1, 
#              y=Var2)) +
#   geom_tile(aes(fill=Freq)) +
#   theme(axis.text.x = element_text(angle = 90)) +
#   scale_fill_gradientn(colours=c("blue","red")) +
#   ggtitle('Neuropeptides percentage overlap reorder filter sub')
# ggsave('neuropeptides/comparison/social.status/Neuropeptides percentage overlap reorder filter sub.png',
#        width = 10,
#        height = 10)
# 


### calculate  count
#dom
dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.df = dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap %>% 
  as.data.frame.table() 
#sub
sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.df = sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap %>% 
  as.data.frame.table() 
# 
# ##graph
# #dom
# dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.df %>% 
#   ggplot(aes(x=Var1, 
#              y=Var2)) +
#   geom_tile(aes(fill=Freq)) +
#   theme(axis.text.x = element_text(angle = 90)) +
#   scale_fill_gradientn(colours=c("blue","red")) +
#   ggtitle('Neuropeptides overlap dom')
# ggsave('neuropeptides/comparison/social.status/Neuropeptides overlap.png',
#        width = 10,
#        height = 10)
# #sub
# sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.df %>% 
#   ggplot(aes(x=Var1, 
#              y=Var2)) +
#   geom_tile(aes(fill=Freq)) +
#   theme(axis.text.x = element_text(angle = 90)) +
#   scale_fill_gradientn(colours=c("blue","red")) +
#   ggtitle('Neuropeptides overlap sub')
# ggsave('neuropeptides/comparison/social.status/Neuropeptides overlap.png',
#        width = 10,
#        height = 10)

#add gene count data
#dom
dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.df.reorder = left_join(dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.df,
                                                                                                                    dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.count,
                                                                                                                    by = c('Var1' = 'Gene.name')) %>% 
  dplyr::rename(Var1.total.count = Total.count) %>% 
  left_join(dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.count,
            by = c('Var2' = 'Gene.name')) %>% 
  dplyr::rename(Var2.total.count = Total.count)
#sub
sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.df.reorder = left_join(sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.df,
                                                                                                                    sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.count,
                                                                                                                    by = c('Var1' = 'Gene.name')) %>% 
  dplyr::rename(Var1.total.count = Total.count) %>% 
  left_join(sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.count,
            by = c('Var2' = 'Gene.name')) %>% 
  dplyr::rename(Var2.total.count = Total.count)

#reorder
#dom
dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.df.reorder <- dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.df.reorder %>%
  arrange(Var1.total.count,
          Var2.total.count) %>%
  mutate(Var1 = fct_inorder(factor(Var1, ordered=TRUE)),
         Var2 = fct_inorder(factor(Var2, ordered=TRUE)))
#sub
sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.df.reorder <- sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.df.reorder %>%
  arrange(Var1.total.count,
          Var2.total.count) %>%
  mutate(Var1 = fct_inorder(factor(Var1, ordered=TRUE)),
         Var2 = fct_inorder(factor(Var2, ordered=TRUE)))

# remove individual overlap
#dom
dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.df.reorder = dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.df.reorder %>% 
  mutate(Freq.na = ifelse(Var1.total.count == Var2.total.count,
                          NA,
                          Freq))
#sub
sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.df.reorder = sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.df.reorder %>% 
  mutate(Freq.na = ifelse(Var1.total.count == Var2.total.count,
                          NA,
                          Freq))
# 
# #graph reorder
# #dom
# dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.df.reorder %>% 
#   ggplot(aes(x=Var1, 
#              y=Var2)) +
#   geom_tile(aes(fill=Freq.na)) +
#   theme(axis.text.x = element_text(angle = 90)) +
#   scale_fill_gradientn(colours=c("blue","red")) +
#   ggtitle('Neuropeptides overlap rescaled dom')
# ggsave('neuropeptides/comparison/social.status/Neuropeptides overlap rescaled dom.png',
#        width = 10,
#        height = 10)
# #sub
# sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.df.reorder %>% 
#   ggplot(aes(x=Var1, 
#              y=Var2)) +
#   geom_tile(aes(fill=Freq.na)) +
#   theme(axis.text.x = element_text(angle = 90)) +
#   scale_fill_gradientn(colours=c("blue","red")) +
#   ggtitle('Neuropeptides overlap rescaled sub')
# ggsave('neuropeptides/comparison/social.status/Neuropeptides overlap rescaled sub.png',
#        width = 10,
#        height = 10)

### compare percentage overlap
## between doms and subs
#status
#set Var1 and Var2 to characters
#rename variables
status.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder = full_join(dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder %>% 
                                                                                                                              mutate(Var1 = as.character(Var1),
                                                                                                                                     Var2 = as.character(Var2)) %>% 
                                                                                                                              dplyr::rename(Freq.dom = Freq,
                                                                                                                                     Var1.total.count.dom = Var1.total.count,
                                                                                                                                     Var2.total.count.dom = Var2.total.count),
                                                                                                                            sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder %>% 
                                                                                                                              mutate(Var1 = as.character(Var1),
                                                                                                                                     Var2 = as.character(Var2))%>% 
                                                                                                                              dplyr::rename(Freq.sub = Freq,
                                                                                                                                     Var1.total.count.sub = Var1.total.count,
                                                                                                                                     Var2.total.count.sub = Var2.total.count))

#find difference dom-sub percent
status.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder = status.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder %>% 
  mutate(Freq.diff = Freq.dom - Freq.sub)

#reorder by 
#status
status.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder <- status.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder %>%
  arrange(Var1.total.count.dom,
          Var2.total.count.dom) %>%
  mutate(Var1 = fct_inorder(factor(Var1, ordered=TRUE)),
         Var2 = fct_inorder(factor(Var2, ordered=TRUE)))
# 
# #graph reorder
# #status
# status.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder %>% 
#   ggplot(aes(x=Var1, 
#              y=Var2)) +
#   geom_tile(aes(fill=Freq.diff)) +
#   theme(axis.text.x = element_text(angle = 90)) +
#   scale_fill_gradientn(colours=c('dark blue', "blue",'white',"red", 'dark red')) +
#   ggtitle('Neuropeptides percentage overlap reorder difference')
# ggsave('neuropeptides/comparison/social.status/Neuropeptides percentage difference overlap reorder status.png',
#        width = 10,
#        height = 10)
# 
# #graph reorder
# #status
# #filter
# status.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder %>% 
#   filter(Var2.total.count.dom > 200 & Var1.total.count.dom > 200) %>% 
#   filter(Var2.total.count.sub > 200 & Var1.total.count.sub >200) %>%
#   ggplot(aes(x=Var1, 
#              y=Var2)) +
#   geom_tile(aes(fill=Freq.diff)) +
#   theme(axis.text.x = element_text(angle = 90)) +
#   scale_fill_gradientn(colours=c('dark blue', "blue",'white',"red", 'dark red'), limits = c(-10, 10)) +
#   ggtitle('Neuropeptides percentage overlap reorder difference filter')
# ggsave('neuropeptides/comparison/social.status/Neuropeptides percentage difference overlap reorder status filter.png',
#        width = 10,
#        height = 10)

##graph nueropeptide counts
#status

#count neuropeptides across social status
status.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder %>% 
  select(Var2, 
         Var2.total.count.dom,
         Var2.total.count.sub) %>% 
  distinct() %>% 
  select(Var2.total.count.dom, 
         Var2.total.count.sub) %>% 
  droplevels %>% 
  colSums()


#dom has 8893 neurons and 8328 expressing neuropeptide
#sub has 7820 neurons and 4498 expressing neuropeptide
# neuron variable (x) and neuropeptide variable (x)
#Var2.total.count.dom=19273, Var2.total.count.sub =6534 
#ratio = 1.469328
#ratio = 0.6805833
# scale dom by: 4498/8328 = 0.3731322

#dom has 6454 neurons and 6061 expressing neuropeptide (0.9391075)
#sub has 6187 neurons and 3608 expressing neuropeptide (0.5715177)
# neuron variable (0.7836) and neuropeptide variable (0.4806)
#Var2.total.count.dom=14188, Var2.total.count.sub =5294
#ratio = 1.469328
#ratio = 0.6805833
# scale dom by: 3608/6061 = 0.3731322
## make list of biased genes
neuropeptides.biased.df = status.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder %>%
  mutate(Dom.count.scaled_0.5401057 = Var2.total.count.dom*0.5401057,
         Sub.count = Var2.total.count.sub, 
         Ratio.scaled = Dom.count.scaled_0.5401057/Sub.count,
         Keep = case_when(Ratio.scaled >= 1.5 & Var2.total.count.dom >= 25 ~ as.character(Var2),
                          Ratio.scaled <= 0.6667 & Var2.total.count.sub >= 25 ~ as.character(Var2),
                          TRUE ~ as.character(NA)),
         Keep.name =  case_when(Var2 == "avp"~ "avp",
                                Var2 == "oxt"~ "oxt",
                                Var2 == "SST"~ "SST",
                                Var2 == "GNRH1"~ "GNRH1",
                                TRUE ~ Keep),
         Bias = case_when(Ratio.scaled >= 1.5 & Var2.total.count.dom >=25 ~ "Dominant bias",
                                Ratio.scaled <= 0.667 & Var2.total.count.sub >=25 ~ "Subordinate bias",
                                TRUE ~ "No bias")) %>% 
  select(Keep, 
         Ratio.scaled,
         Dom.count.scaled_0.5401057,
         Sub.count,
         Bias) %>% 
  filter(!is.na(Keep)) %>% 
  distinct() 


# #ratio = 0.5952813
# status.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder %>%
#   mutate(Dom.count.scaled_0.5952813 = Var2.total.count.dom*0.5952813,
#          Sub.count = Var2.total.count.sub,
#          Ratio.scaled = Dom.count.scaled_0.5952813/Sub.count,
#          Keep = case_when(Ratio.scaled >= 1.5 & Var2.total.count.dom >= 100 ~ as.character(Var2),
#                           Ratio.scaled <= 0.667 & Var2.total.count.sub >= 100 ~ as.character(Var2),
#                           TRUE ~ as.character(NA)),
#          Keep.name =  case_when(Var2 == "avp"~ "avp",
#                                 Var2 == "oxt"~ "oxt",
#                                 Var2 == "GNRH1"~ "GNRH1",
#                                 TRUE ~ Keep),
#          Keep.color = case_when(Ratio.scaled >= 1.5 & Var2.total.count.dom >=100 ~ "Dominant bias",
#                                 Ratio.scaled <= 0.667 & Var2.total.count.sub >=100 ~ "Subordinate bias",
#                                 TRUE ~ "No bias")) %>%
#   select(c(Dom.count.scaled_0.5952813,
#            Sub.count,
#            Keep.name,
#            Keep.color)) %>%
#   distinct() %>%
#   ggplot(aes(x=Dom.count.scaled_0.5952813,
#              y=Sub.count,
#              label = Keep.name)) +
#   geom_abline(slope = 1,
#               intercept = 0) +
#   geom_abline(slope = 0.667,
#               intercept = 0,
#               linetype = 'dashed') +
#   geom_abline(slope = 1.5,
#               intercept = 0,
#               linetype = 'dashed') +
#   geom_point(aes(fill = Keep.color),
#              shape = 21,
#              size = 6) +
#   geom_label_repel()+
#   theme_classic() +
#   ggtitle('Neuropeptide neuron count across social status') +
#   xlab("Dominant neuron count scaled (~60%)") +
#   ylab("Subordinate neuron count") +
#   scale_fill_manual(breaks = c("Dominant bias",
#                                "Subordinate bias",
#                                "No bias"),
#                     values = c("#4e499e",
#                                "#60bb46",
#                                "black")) +
#   theme(legend.position = c(0.8,
#                             0.25),
#         legend.box.background = element_rect(colour = "black"),
#         legend.background = element_blank(),
#         legend.text=element_text(size=20),
#         legend.title = element_blank())+
#   theme(axis.text = element_text(size = 20),
#         axis.title = element_text(size = 20),
#         title =element_text(size=20, face='bold')) + 
#   ylim(0,3000) + 
#   xlim(0,3000)
# ggsave('neuropeptides/comparison_new/social.status/Neuropeptides count comparison across social status scaled.png',
#        width = 10,
#        height = 10)
# 
# #log scaled
# status.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder %>%
#   mutate(Dom.count.scaled_0.5952813 = Var2.total.count.dom*0.5952813,
#          Sub.count = Var2.total.count.sub,
#          Ratio.scaled = Dom.count.scaled_0.5952813/Sub.count,
#          Keep = case_when(Ratio.scaled >= 2 & Var2.total.count.dom >= 25 ~ as.character(Var2),
#                           Ratio.scaled <= 0.5 & Var2.total.count.sub >= 25 ~ as.character(Var2),
#                           TRUE ~ as.character(NA)),
#          Keep.name =  case_when(Var2 == "avp"~ "avp",
#                                 Var2 == "oxt"~ "oxt",
#                                 Var2 == "GNRH1"~ "GNRH1",
#                                 TRUE ~ Keep),
#          Keep.color = case_when(Ratio.scaled >= 2 & Var2.total.count.dom >=25 ~ "Dominant bias",
#                                 Ratio.scaled <= 0.5 & Var2.total.count.sub >=25 ~ "Subordinate bias",
#                                 TRUE ~ "No bias")) %>% 
#   select(c(Dom.count.scaled_0.5952813,
#            Sub.count,
#            Keep.name,
#            Keep.color)) %>%
#   distinct() %>%
#   ggplot(aes(x=Dom.count.scaled_0.5952813,
#              y=Sub.count,
#              label = Keep.name)) +
#   geom_abline(slope = 1,
#               intercept = 0) +
#   geom_point(aes(fill = Keep.color),
#              shape = 21,
#              size = 10) +
#   geom_label_repel(size = 10)+
#   theme_classic() +
#   ggtitle('Neuropeptide neuron count across social status') +
#   xlab("Dominant neuron count scaled (~60%)") +
#   ylab("Subordinate neuron count") +
#   scale_fill_manual(breaks = c("Dominant bias",
#                                "Subordinate bias",
#                                "No bias"),
#                     values = c("#4e499e",
#                                "#60bb46",
#                                "grey")) +
#   theme(legend.position = c(0.8,
#                             0.50),
#         legend.box.background = element_rect(colour = "black"),
#         legend.background = element_blank(),
#         legend.text=element_text(size=20),
#         legend.title = element_blank())+
#   theme(axis.text = element_text(size = 20),
#         axis.title = element_text(size = 20),
#         title =element_text(size=20, face='bold'))  +
#   scale_x_continuous(trans='log10') +
#   scale_y_continuous(trans='log10')
# ggsave('neuropeptides/comparison_new/social.status/Neuropeptides count comparison across social status scaled log.png',
#        width = 10,
#        height = 10)

# #log scaled
# #reduce
status.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder %>%
  mutate(Dom.count.scaled_0.5401057 = Var2.total.count.dom*0.5401057,
         Sub.count = Var2.total.count.sub,
         Ratio.scaled = Dom.count.scaled_0.5401057/Sub.count,
         Keep = case_when(Ratio.scaled >= 2 & Var2.total.count.dom >= 25 ~ as.character(Var2),
                          Ratio.scaled <= 0.5 & Var2.total.count.sub >= 25 ~ as.character(Var2),
                          TRUE ~ as.character(NA)),
         Keep.name =  case_when(Var2 == "avp"~ "avp",
                                Var2 == "oxt"~ "oxt",
                                Var2 == "GNRH1"~ "GNRH1",
                                TRUE ~ Keep),
         Keep.color = case_when(Ratio.scaled >= 2 & Var2.total.count.dom >=25 ~ "Dominant bias",
                                Ratio.scaled <= 0.5 & Var2.total.count.sub >=25 ~ "Subordinate bias",
                                TRUE ~ "No bias")) %>% 
  select(c(Dom.count.scaled_0.5401057,
           Sub.count,
           Keep.name,
           Keep.color)) %>%
  distinct() %>%
  filter(Dom.count.scaled_0.5401057 > 10 | Sub.count > 10) %>% 
  ggplot(aes(x=Dom.count.scaled_0.5401057,
             y=Sub.count,
             label = Keep.name)) +
  geom_abline(slope = 1,
              intercept = 0) +
  geom_point(aes(fill = Keep.color),
             shape = 21,
             size = 10) +
  geom_label_repel(size = 10)+
  theme_classic() +
  ggtitle('Neuropeptide neuron count across social status') +
  xlab("Dominant neuron count scaled (~54%)") +
  ylab("Subordinate neuron count") +
  scale_fill_manual(breaks = c("Dominant bias",
                               "Subordinate bias",
                               "No bias"),
                    values = c("#4e499e",
                               "#60bb46",
                               "grey")) +
  theme(legend.position = c(0.8,
                            0.50),
        legend.box.background = element_rect(colour = "black"),
        legend.background = element_blank(),
        legend.text=element_text(size=20),
        legend.title = element_blank())+
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        title =element_text(size=20, face='bold'))  +
  scale_x_continuous(trans='log10') +
  scale_y_continuous(trans='log10') 
ggsave('neuropeptides/comparison_new/social.status/Neuropeptides count comparison across social status scaled log reduce.png',
       width = 10,
       height = 10)

#ratio = 0.483645
# status.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder %>%
#   mutate(Dom.count.scaled_0.5952813 = Var2.total.count.dom*0.5952813,
#          Sub.count = Var2.total.count.sub,
#          Ratio.scaled = Dom.count.scaled_0.5952813/Sub.count,
#          Keep = case_when(Ratio.scaled >= 1.5 & Var2.total.count.dom >= 100 ~ as.character(Var2),
#                           Ratio.scaled <= 0.667 & Var2.total.count.sub >= 100 ~ as.character(Var2),
#                           TRUE ~ as.character(NA)),
#          Keep.name =  case_when(Var2 == "avp"~ "avp",
#                                 Var2 == "oxt"~ "oxt",
#                                 Var2 == "GNRH1"~ "GNRH1",
#                                 TRUE ~ Keep),
#          Keep.color = case_when(Ratio.scaled >= 1.5 & Var2.total.count.dom >=100 ~ "Dominant bias",
#                                 Ratio.scaled <= 0.667 & Var2.total.count.sub >=100 ~ "Subordinate bias",
#                                 TRUE ~ "No bias")) %>%
#   filter(Dom.count.scaled_0.5952813 < 2000) %>% 
#   select(c(Dom.count.scaled_0.5952813,
#            Sub.count,
#            Keep.name,
#            Keep.color)) %>%
#   distinct() %>%
#   ggplot(aes(x=Dom.count.scaled_0.5952813,
#              y=Sub.count,
#              label = Keep.name)) +
#   geom_abline(slope = 1,
#               intercept = 0) +
#   geom_abline(slope = 0.667,
#               intercept = 0,
#               linetype = 'dashed') +
#   geom_abline(slope = 1.5,
#               intercept = 0,
#               linetype = 'dashed') +
#   geom_point(aes(fill = Keep.color),
#              shape = 21,
#              size = 6) +
#   geom_label_repel()+
#   theme_classic() +
#   ggtitle('Neuropeptide neuron count across social status') +
#   xlab("Dominant neuron count scaled (~60%)") +
#   ylab("Subordinate neuron count") +
#   scale_fill_manual(breaks = c("Dominant bias",
#                                "Subordinate bias",
#                                "No bias"),
#                     values = c("#4e499e",
#                                "#60bb46",
#                                "black")) +
#   theme(legend.position = c(0.8,
#                             0.25),
#         legend.box.background = element_rect(colour = "black"),
#         legend.background = element_blank(),
#         legend.text=element_text(size=20),
#         legend.title = element_blank())+
#   theme(axis.text = element_text(size = 20),
#         axis.title = element_text(size = 20),
#         title =element_text(size=20, face='bold')) + 
#   ylim(0,1000) + 
#   xlim(0,1000)
# ggsave('neuropeptides/comparison_new/social.status/Neuropeptides count comparison across social status scaled outliers.png',
#        width = 10,
#        height = 10)


# #true count
# status.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder %>%
#   mutate(Dom.count.scaled_0.5952813 = Var2.total.count.dom*0.5952813,
#          Sub.count = Var2.total.count.sub,
#          Dom.count = Var2.total.count.dom,
#          Ratio.scaled = Dom.count.scaled_0.5952813/Sub.count,
#          Keep = case_when(Ratio.scaled >= 1.5 & Var2.total.count.dom >= 100 ~ as.character(Var2),
#                           Ratio.scaled <= 0.667 & Var2.total.count.sub >= 100 ~ as.character(Var2),
#                           TRUE ~ as.character(NA)),
#          Keep.name =  case_when(Var2 == "avp"~ "avp",
#                                 Var2 == "oxt"~ "oxt",
#                                 Var2 == "SST"~ "SST",
#                                 Var2 == "GNRH1"~ "GNRH1",
#                                 TRUE ~ Keep),
#          Keep.color = case_when(Ratio.scaled >= 1.5 & Var2.total.count.dom >=100 ~ "Dominant bias",
#                                 Ratio.scaled <= 0.667 & Var2.total.count.sub >=100 ~ "Subordinate bias",
#                                 TRUE ~ "No bias")) %>%
#   select(c(Dom.count,
#            Sub.count,
#            Keep.name,
#            Keep.color)) %>%
#   distinct() %>%
#   ggplot(aes(x=Dom.count,
#              y=Sub.count,
#              label = Keep.name)) +
#   geom_abline(slope = 1,
#               intercept = 0) +
#   geom_abline(slope = 0.667,
#               intercept = 0,
#               linetype = 'dashed') +
#   geom_abline(slope = 1.5,
#               intercept = 0,
#               linetype = 'dashed') +
#   geom_point(aes(fill = Keep.color),
#              shape = 21,
#              size = 4) +
#   geom_label_repel()+
#   theme_classic() +
#   ggtitle('Neuropeptide neuron count across social status') +
#   xlab("Dominant neuron count") +
#   ylab("Subordinate neuron count") +
#   scale_fill_manual(breaks = c("Dominant bias",
#                                "Subordinate bias",
#                                "No bias"),
#                     values = c("#4e499e",
#                                "#60bb46",
#                                "black")) +
#   theme(legend.position = c(0.8,
#                             0.25),
#         legend.box.background = element_rect(colour = "black"),
#         legend.background = element_blank(),
#         legend.text=element_text(size=10),
#         legend.title = element_blank())+
#   theme(axis.text = element_text(size = 10),
#         axis.title = element_text(size = 10),
#         title =element_text(size=10, face='bold'))
# ggsave('neuropeptides/comparison_new/social.status/Neuropeptides count comparison across social status.png',
#        width = 10,
#        height = 10)

#### network ####
### create network of all neuropeptides
### use percentages from both dom and sub
### all neuropeptides
## create igraph dataframes
# create edge list
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist = burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder %>%
  filter(Freq != 100) %>% 
  filter(Var2.total.count >= 50) %>%
  dplyr::rename(From = Var1,
                To = Var2,
                From.total.count = Var1.total.count,
                To.total.count = Var2.total.count)

# create vertices
# need to get all nodes listed
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.vertices = full_join(burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist %>% 
                                                                                                                                 dplyr::select(c(To,
                                                                                                                                                 To.total.count)) %>% 
                                                                                                                                 distinct() %>% 
                                                                                                                                 mutate(To = as.character(To)),
                                                                                                                               burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist %>% 
                                                                                                                                 dplyr::select(c(From,
                                                                                                                                                 From.total.count)) %>% 
                                                                                                                                 distinct() %>% 
                                                                                                                                 mutate(From = as.character(From)) %>% 
                                                                                                                                 dplyr::rename(To = From,
                                                                                                                                               To.total.count = From.total.count))

# create igraph
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder <- graph_from_data_frame(d = burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist,
                                                                                                                                   vertices = burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.vertices,
                                                                                                                                   directed = TRUE)

#ggraph
ggraph(burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder,
       layout = 'fr') +
  geom_edge_parallel(aes(width = Freq,
                         color = Freq,
                         start_cap = label_rect(node1.name),
                         end_cap = label_rect(node2.name)), 
                     arrow = arrow(length = unit(4, 
                                                 'mm'))) +
  geom_node_point(aes(size = To.total.count)) +
  geom_node_text(aes(label = name),
                 nudge_x = 0.15) +
  ggtitle('All neuropeptides percentage, filter > 50') +
  theme_classic()
ggsave('neuropeptides/comparison/network/Neuropeptides network all.png',
       width = 10,
       height = 10)

### all neuropeptides, filter 10%
## create igraph dataframes
# create edge list
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.filter = burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder %>%
  filter(Freq != 100) %>% 
  filter(Freq >= 10) %>% 
  filter(Var2.total.count >= 50) %>% 
  dplyr::rename(From = Var1,
                To = Var2,
                From.total.count = Var1.total.count,
                To.total.count = Var2.total.count)

# create vertices
# need to get all nodes listed
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.vertices.filter = full_join(burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.filter %>% 
  dplyr::select(c(To,
                  To.total.count)) %>% 
  distinct() %>% 
  mutate(To = as.character(To)),
  burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.filter %>% 
    dplyr::select(c(From,
                    From.total.count)) %>% 
    distinct() %>% 
    mutate(From = as.character(From)) %>% 
             dplyr::rename(To = From,
                           To.total.count = From.total.count)) %>% 
  left_join(neuropeptides.biased.df %>% 
              dplyr::select(c(Keep,
                              Bias)) %>% 
              dplyr::rename(To = Keep)) %>% 
  mutate(Bias = ifelse(is.na(Bias),
                       'None',
                       Bias))



# create igraph
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.filter <- graph_from_data_frame(d = burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.filter,
                                                      vertices = burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.vertices.filter,
                                                      directed = TRUE)

#ggraph
ggraph(burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.filter,
       layout = 'fr') +
  geom_edge_parallel(aes(width = Freq,
                    color = Freq,
                    start_cap = label_rect(node1.name),
                    end_cap = label_rect(node2.name)), 
                arrow = arrow(length = unit(4, 
                                            'mm'))) +
  geom_node_point(aes(size = To.total.count)) +
  geom_node_text(aes(label = name),
                 repel = T) +
  ggtitle('All neuropeptides percentage, filter > 10% & > 50') +
  theme_classic()
ggsave('neuropeptides/comparison/network/Neuropeptides network all filter.png',
       width = 10,
       height = 10)

#ggraph
#add bias
ggraph(burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.filter,
       layout = 'fr') +
  geom_edge_parallel(aes(width = Freq,
                         color = Freq,
                         start_cap = label_rect(node1.name),
                         end_cap = label_rect(node2.name)), 
                     arrow = arrow(length = unit(4, 
                                                 'mm'))) +
  geom_node_point(aes(size = To.total.count,
                      color = Bias)) +
  geom_node_text(aes(label = name),
                 repel = T) +
  ggtitle('All neuropeptides percentage, filter > 10% & > 50 bias') +
  theme_classic()
ggsave('neuropeptides/comparison/network/Neuropeptides network all filter bias.png',
       width = 10,
       height = 10)


### dom neuropeptides, filter 10%
## create igraph dataframes
# create edge list
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.dom = dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder %>%
  filter(Freq != 100) %>% 
  filter(Freq >= 10) %>% 
  filter(Var2.total.count >= 50) %>% 
  dplyr::rename(From = Var1,
                To = Var2,
                From.total.count = Var1.total.count,
                To.total.count = Var2.total.count)

# create vertices
# need to get all nodes listed
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.vertices.dom = full_join(burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.dom %>% 
                                                                                                                                        dplyr::select(c(To,
                                                                                                                                                        To.total.count)) %>% 
                                                                                                                                        distinct() %>% 
                                                                                                                                        mutate(To = as.character(To)),
                                                                                                                                      burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.dom %>% 
                                                                                                                                        dplyr::select(c(From,
                                                                                                                                                        From.total.count)) %>% 
                                                                                                                                        distinct() %>% 
                                                                                                                                        mutate(From = as.character(From)) %>% 
                                                                                                                                        dplyr::rename(To = From,
                                                                                                                                                      To.total.count = From.total.count)) %>% 
  left_join(neuropeptides.biased.df %>% 
              dplyr::select(c(Keep,
                              Bias)) %>% 
              dplyr::rename(To = Keep)) %>% 
  mutate(Bias = ifelse(is.na(Bias),
                       'None',
                       Bias))

# create igraph
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.dom <- graph_from_data_frame(d = burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.dom,
                                                                                                                                          vertices = burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.vertices.dom,
                                                                                                                                          directed = TRUE)

#ggraph
ggraph(burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.dom,
       layout = 'fr') +
  geom_edge_parallel(aes(width = Freq,
                         color = Freq,
                         start_cap = label_rect(node1.name),
                         end_cap = label_rect(node2.name)), 
                     arrow = arrow(length = unit(4, 
                                                 'mm'))) +
  geom_node_point(aes(size = To.total.count)) +
  geom_node_text(aes(label = name),
                 repel = T) +
  ggtitle('Dom neuropeptides percentage, filter > 10% & > 50') +
  theme_classic()
ggsave('neuropeptides/comparison/network/Neuropeptides network all filter dom.png',
       width = 10,
       height = 10)

#ggraph
#bias
ggraph(burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.dom,
       layout = 'fr') +
  geom_edge_parallel(aes(width = Freq,
                         color = Freq,
                         start_cap = label_rect(node1.name),
                         end_cap = label_rect(node2.name)), 
                     arrow = arrow(length = unit(4, 
                                                 'mm'))) +
  geom_node_point(aes(size = To.total.count,
                      color = Bias)) +
  geom_node_text(aes(label = name),
                 repel = T) +
  ggtitle('Dom neuropeptides percentage, filter > 10% & > 50 bias') +
  theme_classic()
ggsave('neuropeptides/comparison/network/Neuropeptides network all filter dom bias.png',
       width = 10,
       height = 10)

### dom neuropeptides, filter 10%, scaled
## create igraph dataframes
# create edge list
# 50*2.080786 = 63.8044 104.0393
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.dom.scaled = dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder %>%
  filter(Freq != 100) %>% 
  filter(Freq >= 10) %>% 
  filter(Var2.total.count >= 104) %>% 
  dplyr::rename(From = Var1,
                To = Var2,
                From.total.count = Var1.total.count,
                To.total.count = Var2.total.count)

# create vertices
# need to get all nodes listed
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.vertices.dom.scaled = full_join(burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.dom.scaled %>% 
                                                                                                                                     dplyr::select(c(To,
                                                                                                                                                     To.total.count)) %>% 
                                                                                                                                     distinct() %>% 
                                                                                                                                     mutate(To = as.character(To)),
                                                                                                                                   burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.dom.scaled %>% 
                                                                                                                                     dplyr::select(c(From,
                                                                                                                                                     From.total.count)) %>% 
                                                                                                                                     distinct() %>% 
                                                                                                                                     mutate(From = as.character(From)) %>% 
                                                                                                                                     dplyr::rename(To = From,
                                                                                                                                                   To.total.count = From.total.count)) %>% 
  left_join(neuropeptides.biased.df %>% 
              dplyr::select(c(Keep,
                              Bias)) %>% 
              dplyr::rename(To = Keep)) %>% 
  mutate(Bias = ifelse(is.na(Bias),
                       'None',
                       Bias))

# create igraph
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.dom.scaled <- graph_from_data_frame(d = burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.dom.scaled,
                                                                                                                                       vertices = burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.vertices.dom.scaled,
                                                                                                                                       directed = TRUE)

#ggraph
ggraph(burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.dom.scaled,
       layout = 'fr') +
  geom_edge_parallel(aes(width = Freq,
                         color = Freq,
                         start_cap = label_rect(node1.name),
                         end_cap = label_rect(node2.name)), 
                     arrow = arrow(length = unit(4, 
                                                 'mm'))) +
  geom_node_point(aes(size = To.total.count)) +
  geom_node_text(aes(label = name),
                 repel = T) +
  ggtitle('Dom neuropeptides percentage, filter > 10% & > 104, scaled') +
  theme_classic()
ggsave('neuropeptides/comparison/network/Neuropeptides network all filter dom scaled.png',
       width = 10,
       height = 10)

#ggraph
#bias
ggraph(burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.dom.scaled,
       layout = 'fr') +
  geom_edge_parallel(aes(width = Freq,
                         color = Freq,
                         start_cap = label_rect(node1.name),
                         end_cap = label_rect(node2.name)), 
                     arrow = arrow(length = unit(4, 
                                                 'mm'))) +
  geom_node_point(aes(size = To.total.count,
                      color = Bias)) +
  geom_node_text(aes(label = name),
                 repel = T) +
  ggtitle('Dom neuropeptides percentage, filter > 10% & > 74, scaled bias') +
  theme_classic()
ggsave('neuropeptides/comparison/network/Neuropeptides network all filter dom scaled bias.png',
       width = 10,
       height = 10)

### sub neuropeptides, filter 10%
## create igraph dataframes
# create edge list
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.sub = sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder %>%
  filter(Freq != 100) %>% 
  filter(Freq >= 10) %>% 
  filter(Var2.total.count >= 50) %>% 
  dplyr::rename(From = Var1,
                To = Var2,
                From.total.count = Var1.total.count,
                To.total.count = Var2.total.count)

# create vertices
# need to get all nodes listed
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.vertices.sub = full_join(burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.sub %>% 
                                                                                                                                     dplyr::select(c(To,
                                                                                                                                                     To.total.count)) %>% 
                                                                                                                                     distinct() %>% 
                                                                                                                                     mutate(To = as.character(To)),
                                                                                                                                   burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.sub %>% 
                                                                                                                                     dplyr::select(c(From,
                                                                                                                                                     From.total.count)) %>% 
                                                                                                                                     distinct() %>% 
                                                                                                                                     mutate(From = as.character(From)) %>% 
                                                                                                                                     dplyr::rename(To = From,
                                                                                                                                                   To.total.count = From.total.count)) %>% 
  left_join(neuropeptides.biased.df %>% 
              dplyr::select(c(Keep,
                              Bias)) %>% 
              dplyr::rename(To = Keep)) %>% 
  mutate(Bias = ifelse(is.na(Bias),
                       'None',
                       Bias))

# create igraph
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.sub <- graph_from_data_frame(d = burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.sub,
                                                                                                                                       vertices = burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.vertices.sub,
                                                                                                                                       directed = TRUE)

#ggraph
ggraph(burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.sub,
       layout = 'fr') +
  geom_edge_parallel(aes(width = Freq,
                         color = Freq,
                         start_cap = label_rect(node1.name),
                         end_cap = label_rect(node2.name)), 
                     arrow = arrow(length = unit(4, 
                                                 'mm'))) +
  geom_node_point(aes(size = To.total.count)) +
  geom_node_text(aes(label = name),
                 repel = T) +
  ggtitle('sub neuropeptides percentage, filter > 10% & > 50') +
  theme_classic()
ggsave('neuropeptides/comparison/network/Neuropeptides network all filter sub.png',
       width = 10,
       height = 10)

#ggraph
#bias
ggraph(burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.sub,
       layout = 'fr') +
  geom_edge_parallel(aes(width = Freq,
                         color = Freq,
                         start_cap = label_rect(node1.name),
                         end_cap = label_rect(node2.name)), 
                     arrow = arrow(length = unit(4, 
                                                 'mm'))) +
  geom_node_point(aes(size = To.total.count,
                      color = Bias)) +
  geom_node_text(aes(label = name),
                 repel = T) +
  ggtitle('sub neuropeptides percentage, filter > 10% & > 50 bias') +
  theme_classic()
ggsave('neuropeptides/comparison/network/Neuropeptides network all filter sub bias.png',
       width = 10,
       height = 10)

### dom neuropeptides, filter 10%
## create igraph dataframes
# create edge list
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.dom.bias = dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder %>%
  filter(Freq != 100) %>% 
  filter(Freq >= 10) %>% 
  filter(Var2.total.count >= 50) %>% 
  dplyr::rename(From = Var1,
                To = Var2,
                From.total.count = Var1.total.count,
                To.total.count = Var2.total.count) %>% 
  filter(!From %in% c(neuropeptides.biased.df$Keep)) %>% 
  filter(!To %in% c(neuropeptides.biased.df$Keep))


# create vertices
# need to get all nodes listed
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.vertices.dom.bias = full_join(burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.dom.bias %>% 
                                                                                                                                          dplyr::select(c(To,
                                                                                                                                                          To.total.count)) %>% 
                                                                                                                                          distinct() %>% 
                                                                                                                                          mutate(To = as.character(To)),
                                                                                                                                        burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.dom.bias %>% 
                                                                                                                                          dplyr::select(c(From,
                                                                                                                                                          From.total.count)) %>% 
                                                                                                                                          distinct() %>% 
                                                                                                                                          mutate(From = as.character(From)) %>% 
                                                                                                                                          dplyr::rename(To = From,
                                                                                                                                                        To.total.count = From.total.count)) %>% 
  left_join(neuropeptides.biased.df %>% 
              dplyr::select(c(Keep,
                              Bias)) %>% 
              dplyr::rename(To = Keep)) %>% 
  mutate(Bias = ifelse(is.na(Bias),
                       'None',
                       Bias))

# create igraph
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.dom.bias <- graph_from_data_frame(d = burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.dom.bias,
                                                                                                                                            vertices = burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.vertices.dom.bias,
                                                                                                                                            directed = TRUE)

#ggraph
ggraph(burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.dom.bias,
       layout = 'fr') +
  geom_edge_parallel(aes(width = Freq,
                         color = Freq,
                         start_cap = label_rect(node1.name),
                         end_cap = label_rect(node2.name)), 
                     arrow = arrow(length = unit(4, 
                                                 'mm'))) +
  geom_node_point(aes(size = To.total.count)) +
  geom_node_text(aes(label = name),
                 repel = T) +
  ggtitle('dom neuropeptides percentage, filter > 10% & > 50, bias removed') +
  theme_classic()
ggsave('neuropeptides/comparison/network/Neuropeptides network all filter dom bias removed.png',
       width = 10,
       height = 10)

### sub neuropeptides, filter 10%
## create igraph dataframes
# create edge list
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.sub.bias = sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder %>%
  filter(Freq != 100) %>% 
  filter(Freq >= 10) %>% 
  filter(Var2.total.count >= 50) %>% 
  dplyr::rename(From = Var1,
                To = Var2,
                From.total.count = Var1.total.count,
                To.total.count = Var2.total.count) %>% 
  filter(!From %in% c(neuropeptides.biased.df$Keep)) %>% 
  filter(!To %in% c(neuropeptides.biased.df$Keep))


# create vertices
# need to get all nodes listed
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.vertices.sub.bias = full_join(burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.sub.bias %>% 
                                                                                                                                     dplyr::select(c(To,
                                                                                                                                                     To.total.count)) %>% 
                                                                                                                                     distinct() %>% 
                                                                                                                                     mutate(To = as.character(To)),
                                                                                                                                   burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.sub.bias %>% 
                                                                                                                                     dplyr::select(c(From,
                                                                                                                                                     From.total.count)) %>% 
                                                                                                                                     distinct() %>% 
                                                                                                                                     mutate(From = as.character(From)) %>% 
                                                                                                                                     dplyr::rename(To = From,
                                                                                                                                                   To.total.count = From.total.count)) %>% 
  left_join(neuropeptides.biased.df %>% 
              dplyr::select(c(Keep,
                              Bias)) %>% 
              dplyr::rename(To = Keep)) %>% 
  mutate(Bias = ifelse(is.na(Bias),
                       'None',
                       Bias))

# create igraph
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.sub.bias <- graph_from_data_frame(d = burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.sub.bias,
                                                                                                                                       vertices = burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.vertices.sub.bias,
                                                                                                                                       directed = TRUE)

#ggraph
ggraph(burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.sub.bias,
       layout = 'fr') +
  geom_edge_parallel(aes(width = Freq,
                         color = Freq,
                         start_cap = label_rect(node1.name),
                         end_cap = label_rect(node2.name)), 
                     arrow = arrow(length = unit(4, 
                                                 'mm'))) +
  geom_node_point(aes(size = To.total.count)) +
  geom_node_text(aes(label = name),
                 repel = T) +
  ggtitle('sub neuropeptides percentage, filter > 10% & > 50, bias removed') +
  theme_classic()
ggsave('neuropeptides/comparison/network/Neuropeptides network all filter sub bias removed.png',
       width = 10,
       height = 10)


#### clustering test ####

###pvclust
library(pvclust)

test <- pvclust(dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap %>%  
                  scale(),
                parallel = T)

# Create dendro
plot(test)
pvrect(test, alpha=0.89, pv="au")

# Create dendro
#equal height
plot(test, 
     hang = -1)
pvrect(test, alpha=0.89, pv="au")

#order
test.order = test[["hclust"]][["order"]] %>% 
  as.data.frame() %>% 
  dplyr::rename(position = ".") %>% 
  rownames_to_column('order') %>%
  mutate(order = as.integer(order)) 


#reorder by columnames
test.gene.order = dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap %>% 
  colnames() %>% 
  as.data.frame() %>% 
  rownames_to_column('position') %>% 
  dplyr::rename(gene = ".") %>%  
  mutate(position = as.integer(position)) %>%
  full_join(test.order) %>% 
  arrange(order) %>% 
  select(-c(position))

#


dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.df.reorder %>%
  full_join(test.gene.order %>% 
              dplyr::rename(Var1 = gene,
                            order1 = order)) %>% 
  full_join(test.gene.order %>% 
              dplyr::rename(Var2 = gene,
                            order2 = order)) %>% 
  reorder(order1,
          order2) %>% 
  ggplot(aes(x=Var1, 
             y=Var2)) +
  geom_tile(aes(fill=Freq.na)) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradientn(colours=c("blue","red")) +
  ggtitle('Neuropeptides overlap rescaled dom')

dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.df.reorder %>% 
  ggplot(aes(x=Var1, 
             y=Var2)) +
  geom_tile(aes(fill=Freq.na)) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradientn(colours=c("blue","red")) +
  ggtitle('Neuropeptides overlap rescaled dom')

## sub
test2 <- pvclust(sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap %>%  
                   scale(),
                 parallel = T)

# Create dendro
plot(test2)
pvrect(test2, alpha=0.89, pv="au")

test2.order = test2[["hclust"]][["order"]] %>% 
  as.data.frame() %>% 
  dplyr::rename(position = ".")

#### for poster ####
### dom neuropeptides, filter 10%, scaled
## create igraph dataframes
# create edge list
# 100*2.080786 = 208.0786
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.dom.scaled = dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder %>%
  filter(Freq != 100) %>%
  filter(Freq >= 10) %>% 
  filter(Var2.total.count >= 208) %>% 
  dplyr::rename(From = Var1,
                To = Var2,
                From.total.count = Var1.total.count,
                To.total.count = Var2.total.count)

# create vertices
# need to get all nodes listed
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.vertices.dom.scaled = full_join(burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.dom.scaled %>% 
                                                                                                                                            dplyr::select(c(To,
                                                                                                                                                            To.total.count)) %>% 
                                                                                                                                            distinct() %>% 
                                                                                                                                            mutate(To = as.character(To)),
                                                                                                                                          burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.dom.scaled %>% 
                                                                                                                                            dplyr::select(c(From,
                                                                                                                                                            From.total.count)) %>% 
                                                                                                                                            distinct() %>% 
                                                                                                                                            mutate(From = as.character(From)) %>% 
                                                                                                                                            dplyr::rename(To = From,
                                                                                                                                                          To.total.count = From.total.count)) %>% 
  left_join(neuropeptides.biased.df %>% 
              dplyr::select(c(Keep,
                              Bias)) %>% 
              dplyr::rename(To = Keep)) %>% 
  mutate(Bias = ifelse(is.na(Bias),
                       'No bias',
                       Bias))

# create igraph
#circle
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.dom.scaled <- graph_from_data_frame(d = burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.dom.scaled,
                                                                                                                                              vertices = burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.vertices.dom.scaled,
                                                                                                                                              directed = TRUE)

#ggraph
ggraph(burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.dom.scaled,
       layout = 'linear',
       circular = TRUE,
       sort.by = To.total.count) +
  geom_edge_parallel(aes(width = Freq,
                         alpha = Freq,
                         start_cap = label_rect(node1.name),
                         end_cap = label_rect(node2.name)), 
                     arrow = arrow(length = unit(4, 
                                                 'mm')),
                     color = 'dimgrey') +
  geom_node_point(aes(size = To.total.count,
                      fill = Bias),
                  shape = 21) +
  geom_node_label(aes(label = name),
                 repel = T) +
  ggtitle('Dominant') +
  theme_classic()  +
  scale_fill_manual(breaks = c("Dominant bias",
                               "Subordinate bias",
                               "No bias"),
                    values = c("#4e499e",
                               "#60bb46",
                               "black")) +
  theme(plot.title = element_text(hjust = 0.5),
        title =element_text(size=10, face='bold'),
        axis.line =element_blank(),
        axis.title =element_blank(),
        axis.ticks  =element_blank(),
        axis.text =element_blank())
ggsave('neuropeptides/comparison/network/Neuropeptides network all filter dom scaled bias circle poster.pdf',
       width = 5,
       height = 5)

# #ggraph
# # lgl
# ggraph(burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.dom.scaled,
#        layout = 'lgl') +
#   geom_edge_parallel(aes(width = Freq,
#                          alpha = Freq,
#                          start_cap = label_rect(node1.name),
#                          end_cap = label_rect(node2.name)), 
#                      arrow = arrow(length = unit(4, 
#                                                  'mm')),
#                      color = 'dimgrey') +
#   geom_node_point(aes(size = To.total.count,
#                       fill = Bias),
#                   shape = 21) +
#   geom_node_label(aes(label = name),
#                   repel = T) +
#   scale_size(range = c(2, 8)) +
#   ggtitle('Dominant') +
#   theme_classic()  +
#   scale_fill_manual(breaks = c("Dominant bias",
#                                "Subordinate bias",
#                                "No bias"),
#                     values = c("#4e499e",
#                                "#60bb46",
#                                "black")) +
#   theme(legend.position = 'none',
#         plot.title = element_text(hjust = 0.5),
#         title =element_text(size=10, face='bold'),
#         axis.line =element_blank(),
#         axis.title =element_blank(),
#         axis.ticks  =element_blank(),
#         axis.text =element_blank())
# ggsave('neuropeptides/comparison/network/Neuropeptides network all filter dom scaled bias poster.pdf',
#        width = 5,
#        height = 5)

### sub neuropeptides, filter 10%
## create igraph dataframes
# create edge list
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.sub.scaled = sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder %>%
  filter(Freq != 100) %>% 
  filter(Freq >= 10) %>% 
  filter(Var2.total.count >= 100) %>% 
  dplyr::rename(From = Var1,
                To = Var2,
                From.total.count = Var1.total.count,
                To.total.count = Var2.total.count)

# create vertices
# need to get all nodes listed
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.vertices.sub.scaled = full_join(burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.sub.scaled %>% 
                                                                                                                                            dplyr::select(c(To,
                                                                                                                                                            To.total.count)) %>% 
                                                                                                                                            distinct() %>% 
                                                                                                                                            mutate(To = as.character(To)),
                                                                                                                                          burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.sub.scaled %>% 
                                                                                                                                            dplyr::select(c(From,
                                                                                                                                                            From.total.count)) %>% 
                                                                                                                                            distinct() %>% 
                                                                                                                                            mutate(From = as.character(From)) %>% 
                                                                                                                                            dplyr::rename(To = From,
                                                                                                                                                          To.total.count = From.total.count)) %>% 
  left_join(neuropeptides.biased.df %>% 
              dplyr::select(c(Keep,
                              Bias)) %>% 
              dplyr::rename(To = Keep)) %>% 
  mutate(Bias = ifelse(is.na(Bias),
                       'No bias',
                       Bias))

# create igraph
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.sub.scaled <- graph_from_data_frame(d = burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.sub.scaled,
                                                                                                                                              vertices = burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.vertices.sub.scaled,
                                                                                                                                              directed = TRUE)

#ggraph
# circle
ggraph(burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.sub.scaled,
       layout = 'linear',
       circular = TRUE,
       sort.by = To.total.count) +
  geom_edge_parallel(aes(width = Freq,
                         alpha = Freq,
                         start_cap = label_rect(node1.name),
                         end_cap = label_rect(node2.name)), 
                     arrow = arrow(length = unit(4, 
                                                 'mm')),
                     color = 'dimgrey') +
  geom_node_point(aes(size = To.total.count,
                      fill = Bias),
                  shape = 21) +
  geom_node_label(aes(label = name),
                  repel = T) +
  ggtitle('Subordinate') +
  theme_classic()  +
  scale_fill_manual(breaks = c("Dominant bias",
                               "Subordinate bias",
                               "No bias"),
                    values = c("#4e499e",
                               "#60bb46",
                               "black")) +
  theme(plot.title = element_text(hjust = 0.5),
        title =element_text(size=10, face='bold'),
        axis.line =element_blank(),
        axis.title =element_blank(),
        axis.ticks  =element_blank(),
        axis.text =element_blank())
ggsave('neuropeptides/comparison/network/Neuropeptides network all filter sub scaled bias circle poster.pdf',
       width = 5,
       height = 5)

#ggraph
# lgl
ggraph(burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.sub.scaled,
       layout = 'fr') +
  geom_edge_parallel(aes(width = Freq,
                         alpha = Freq,
                         start_cap = label_rect(node1.name),
                         end_cap = label_rect(node2.name)), 
                     arrow = arrow(length = unit(4, 
                                                 'mm')),
                     color = 'dimgrey') +
  geom_node_point(aes(size = To.total.count,
                      fill = Bias),
                  shape = 21) +
  geom_node_label(aes(label = name),
                  repel = T) +
  scale_size(range = c(2, 8)) +
  ggtitle('Subordinate') +
  theme_classic()  +
  scale_fill_manual(breaks = c("Dominant bias",
                               "Subordinate bias",
                               "No bias"),
                    values = c("#4e499e",
                               "#60bb46",
                               "black")) +
  theme(legend.position = 'none',
        plot.title = element_text(hjust = 0.5),
        title =element_text(size=10, face='bold'),
        axis.line =element_blank(),
        axis.title =element_blank(),
        axis.ticks  =element_blank(),
        axis.text =element_blank())
ggsave('neuropeptides/comparison/network/Neuropeptides network all filter sub scaled bias poster.pdf',
       width = 5,
       height = 5)


#### for presentation ####
### dom neuropeptides, filter 10%, scaled
## create igraph dataframes
# create edge list
# filter out Freq less than 10
# filter out if cell count is below 10
# filter out bias nodes Oxt, prl, gnrh1
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.dom.scaled = dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder %>% 
  filter(Freq != 100) %>%
  dplyr::rename(From = Var1,
                To = Var2,
                From.total.count = Var1.total.count,
                To.total.count = Var2.total.count) %>% 
  mutate(Freq = ifelse(Freq < 10,
                       NA,
                       Freq)) %>% 
  mutate(Freq = ifelse(From.total.count < 20 | To.total.count < 20,
                       NA,
                       Freq)) %>% 
  filter(To != 'oxt',
         To != 'prl',
         To != 'GNRH1',
         From != 'oxt',
         From != 'prl',
         From != 'GNRH1')

# create vertices
# need to get all nodes listed
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.vertices.dom.scaled = full_join(burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.dom.scaled %>% 
                                                                                                                                            dplyr::select(c(To,
                                                                                                                                                            To.total.count)) %>% 
                                                                                                                                            distinct() %>% 
                                                                                                                                            mutate(To = as.character(To)),
                                                                                                                                          burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.dom.scaled %>% 
                                                                                                                                            dplyr::select(c(From,
                                                                                                                                                            From.total.count)) %>% 
                                                                                                                                            distinct() %>% 
                                                                                                                                            mutate(From = as.character(From)) %>% 
                                                                                                                                            dplyr::rename(To = From,
                                                                                                                                                          To.total.count = From.total.count)) %>% 
  left_join(neuropeptides.biased.df %>% 
              dplyr::select(c(Keep,
                              Bias)) %>% 
              dplyr::rename(To = Keep)) %>% 
  mutate(Bias = ifelse(is.na(Bias),
                       'No bias',
                       Bias)) %>% 
  mutate(Label = ifelse(To.total.count < 20,
                        NA,
                        To))

# create igraph
#circle
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.dom.scaled <- graph_from_data_frame(d = burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.dom.scaled,
                                                                                                                                              vertices = burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.vertices.dom.scaled,
                                                                                                                                              directed = TRUE)

#ggraph
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.dom.scaled %>% 
  ggraph(layout = 'lgl') +
  geom_edge_parallel(aes(width = Freq,
                         alpha = Freq,
                         start_cap = label_rect(node1.name),
                         end_cap = label_rect(node2.name)), 
                     arrow = arrow(length = unit(4, 
                                                 'mm')),
                     color = 'dimgrey') +
  geom_node_point(aes(size = To.total.count),
                  shape = 21,
                  fill = "#4e499e") +
  geom_node_label(aes(label = Label),
                  repel = T) +
  ggtitle('Dominant') +
  theme_classic()  +
  theme(plot.title = element_text(hjust = 0.5),
        title =element_text(size=10, face='bold'),
        axis.line =element_blank(),
        axis.title =element_blank(),
        axis.ticks  =element_blank(),
        axis.text =element_blank()) +
  scale_edge_width(range = c(1, 6))+ 
  scale_size(range = c(1,7.5)) +
  theme(legend.position="none")
ggsave('neuropeptides/comparison/network/Neuropeptides network all filter dom scaled bias circle presentation outlier.pdf',
       width = 5,
       height = 5)


### sub neuropeptides, filter 10%, scaled
## create igraph dataframes
# create edge list
# filter out Freq less than 10
# filter out if cell count is below 10
# filter out bias nodes Oxt, prl, gnrh1
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.sub.scaled = sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder %>% 
  filter(Freq != 100) %>%
  dplyr::rename(From = Var1,
                To = Var2,
                From.total.count = Var1.total.count,
                To.total.count = Var2.total.count) %>% 
  mutate(Freq = ifelse(Freq < 10,
                       NA,
                       Freq)) %>% 
  mutate(Freq = ifelse(From.total.count < 20 | To.total.count < 20,
                       NA,
                       Freq)) %>% 
  filter(To != 'oxt',
         To != 'prl',
         To != 'GNRH1',
         From != 'oxt',
         From != 'prl',
         From != 'GNRH1')

# create vertices
# need to get all nodes listed
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.vertices.sub.scaled = full_join(burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.sub.scaled %>% 
                                                                                                                                            dplyr::select(c(To,
                                                                                                                                                            To.total.count)) %>% 
                                                                                                                                            distinct() %>% 
                                                                                                                                            mutate(To = as.character(To)),
                                                                                                                                          burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.sub.scaled %>% 
                                                                                                                                            dplyr::select(c(From,
                                                                                                                                                            From.total.count)) %>% 
                                                                                                                                            distinct() %>% 
                                                                                                                                            mutate(From = as.character(From)) %>% 
                                                                                                                                            dplyr::rename(To = From,
                                                                                                                                                          To.total.count = From.total.count)) %>% 
  left_join(neuropeptides.biased.df %>% 
              dplyr::select(c(Keep,
                              Bias)) %>% 
              dplyr::rename(To = Keep)) %>% 
  mutate(Bias = ifelse(is.na(Bias),
                       'No bias',
                       Bias)) %>% 
  mutate(Label = ifelse(To.total.count < 20,
                        NA,
                        To))

# create igraph
#circle
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.sub.scaled <- graph_from_data_frame(d = burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.sub.scaled,
                                                                                                                                              vertices = burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.vertices.sub.scaled,
                                                                                                                                              directed = TRUE)

#ggraph
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.sub.scaled %>% 
  ggraph(layout = 'lgl') +
  geom_edge_parallel(aes(width = Freq,
                         alpha = Freq,
                         start_cap = label_rect(node1.name),
                         end_cap = label_rect(node2.name)), 
                     arrow = arrow(length = unit(4, 
                                                 'mm')),
                     color = 'dimgrey') +
  geom_node_point(aes(size = To.total.count),
                  shape = 21,
                  fill = "#60bb46") +
  geom_node_label(aes(label = Label),
                  repel = T) +
  ggtitle('Subordinate') +
  theme_classic()  +
  theme(plot.title = element_text(hjust = 0.5),
        title =element_text(size=10, face='bold'),
        axis.line =element_blank(),
        axis.title =element_blank(),
        axis.ticks  =element_blank(),
        axis.text =element_blank()) +
  scale_edge_width(range = c(1, 4))+ 
  scale_size(range = c(1,5)) +
  theme(legend.position="none")
ggsave('neuropeptides/comparison/network/Neuropeptides network all filter sub scaled bias circle presentation outlier.pdf',
       width = 5,
       height = 5)

### dom neuropeptides, filter 10%, scaled
## create igraph dataframes
# create edge list
# filter out Freq less than 10
# filter out if cell count is below 10
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.dom.scaled = dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder %>% 
  filter(Freq != 100) %>%
  dplyr::rename(From = Var1,
                To = Var2,
                From.total.count = Var1.total.count,
                To.total.count = Var2.total.count) %>% 
  mutate(Freq = ifelse(Freq < 10,
                       NA,
                       Freq)) %>% 
  mutate(Freq = ifelse(From.total.count < 20 | To.total.count < 20,
                       NA,
                       Freq)) 

# create vertices
# need to get all nodes listed
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.vertices.dom.scaled = full_join(burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.dom.scaled %>% 
                                                                                                                                            dplyr::select(c(To,
                                                                                                                                                            To.total.count)) %>% 
                                                                                                                                            distinct() %>% 
                                                                                                                                            mutate(To = as.character(To)),
                                                                                                                                          burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.dom.scaled %>% 
                                                                                                                                            dplyr::select(c(From,
                                                                                                                                                            From.total.count)) %>% 
                                                                                                                                            distinct() %>% 
                                                                                                                                            mutate(From = as.character(From)) %>% 
                                                                                                                                            dplyr::rename(To = From,
                                                                                                                                                          To.total.count = From.total.count)) %>% 
  left_join(neuropeptides.biased.df %>% 
              dplyr::select(c(Keep,
                              Bias)) %>% 
              dplyr::rename(To = Keep)) %>% 
  mutate(Bias = ifelse(is.na(Bias),
                       'No bias',
                       Bias)) %>% 
  mutate(Label = ifelse(To.total.count < 20,
                        NA,
                        To))

# create igraph
#circle
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.dom.scaled <- graph_from_data_frame(d = burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.dom.scaled,
                                                                                                                                              vertices = burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.vertices.dom.scaled,
                                                                                                                                              directed = TRUE)

#ggraph
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.dom.scaled %>% 
  ggraph(layout = 'lgl') +
  geom_edge_parallel(aes(width = Freq,
                         alpha = Freq,
                         start_cap = label_rect(node1.name),
                         end_cap = label_rect(node2.name)), 
                     arrow = arrow(length = unit(4, 
                                                 'mm')),
                     color = 'dimgrey') +
  geom_node_point(aes(size = To.total.count),
                  shape = 21,
                  fill = "#4e499e") +
  geom_node_label(aes(label = Label),
                  repel = T) +
  ggtitle('Dominant') +
  theme_classic()  +
  theme(plot.title = element_text(hjust = 0.5),
        title =element_text(size=10, face='bold'),
        axis.line =element_blank(),
        axis.title =element_blank(),
        axis.ticks  =element_blank(),
        axis.text =element_blank()) +
  scale_edge_width(range = c(1, 8))+ 
  scale_size(range = c(1,9)) +
  theme(legend.position="none")
ggsave('neuropeptides/comparison/network/Neuropeptides network all filter dom scaled bias circle presentation.pdf',
       width = 5,
       height = 5)

### sub neuropeptides, filter 10%, scaled
## create igraph dataframes
# create edge list
# filter out Freq less than 10
# filter out if cell count is below 10
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.sub.scaled = sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder %>% 
  filter(Freq != 100) %>%
  dplyr::rename(From = Var1,
                To = Var2,
                From.total.count = Var1.total.count,
                To.total.count = Var2.total.count) %>% 
  mutate(Freq = ifelse(Freq < 10,
                       NA,
                       Freq)) %>% 
  mutate(Freq = ifelse(From.total.count < 20 | To.total.count < 20,
                       NA,
                       Freq)) 

# create vertices
# need to get all nodes listed
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.vertices.sub.scaled = full_join(burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.sub.scaled %>% 
                                                                                                                                            dplyr::select(c(To,
                                                                                                                                                            To.total.count)) %>% 
                                                                                                                                            distinct() %>% 
                                                                                                                                            mutate(To = as.character(To)),
                                                                                                                                          burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.sub.scaled %>% 
                                                                                                                                            dplyr::select(c(From,
                                                                                                                                                            From.total.count)) %>% 
                                                                                                                                            distinct() %>% 
                                                                                                                                            mutate(From = as.character(From)) %>% 
                                                                                                                                            dplyr::rename(To = From,
                                                                                                                                                          To.total.count = From.total.count)) %>% 
  left_join(neuropeptides.biased.df %>% 
              dplyr::select(c(Keep,
                              Bias)) %>% 
              dplyr::rename(To = Keep)) %>% 
  mutate(Bias = ifelse(is.na(Bias),
                       'No bias',
                       Bias)) %>% 
  mutate(Label = ifelse(To.total.count < 20,
                        NA,
                        To))

# create igraph
#circle
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.sub.scaled <- graph_from_data_frame(d = burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.sub.scaled,
                                                                                                                                              vertices = burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.vertices.sub.scaled,
                                                                                                                                              directed = TRUE)

#ggraph
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.sub.scaled %>% 
  ggraph(layout = 'lgl') +
  geom_edge_parallel(aes(width = Freq,
                         alpha = Freq,
                         start_cap = label_rect(node1.name),
                         end_cap = label_rect(node2.name)), 
                     arrow = arrow(length = unit(4, 
                                                 'mm')),
                     color = 'dimgrey') +
  geom_node_point(aes(size = To.total.count),
                  shape = 21,
                  fill = "#60bb46") +
  geom_node_label(aes(label = Label),
                  repel = T) +
  ggtitle('Subordinate') +
  theme_classic()  +
  theme(plot.title = element_text(hjust = 0.5),
        title =element_text(size=10, face='bold'),
        axis.line =element_blank(),
        axis.title =element_blank(),
        axis.ticks  =element_blank(),
        axis.text =element_blank()) +
  scale_edge_width(range = c(1, 4))+ 
  scale_size(range = c(1,5)) +
  theme(legend.position="none")
ggsave('neuropeptides/comparison/network/Neuropeptides network all filter sub scaled bias circle presentation.pdf',
       width = 5,
       height = 5)
  
  
  
  
  
  
#### networksfor presentation  ####  
### dom neuropeptides, filter 10%, scaled
## create igraph dataframes
# create edge list
# filter out Freq less than 10
# filter out if cell count is below 10
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.dom.scaled = dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder %>% 
  filter(Freq != 100) %>%
  dplyr::rename(From = Var1,
                To = Var2,
                From.total.count = Var1.total.count,
                To.total.count = Var2.total.count) %>% 
  mutate(Freq = ifelse(Freq < 10,
                       NA,
                       Freq)) %>% 
  mutate(Freq = ifelse(From.total.count < 20 | To.total.count < 20,
                       NA,
                       Freq)) 

# create vertices
# need to get all nodes listed
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.vertices.dom.scaled = full_join(burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.dom.scaled %>% 
                                                                                                                                            dplyr::select(c(To,
                                                                                                                                                            To.total.count)) %>% 
                                                                                                                                            distinct() %>% 
                                                                                                                                            mutate(To = as.character(To)),
                                                                                                                                          burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.dom.scaled %>% 
                                                                                                                                            dplyr::select(c(From,
                                                                                                                                                            From.total.count)) %>% 
                                                                                                                                            distinct() %>% 
                                                                                                                                            mutate(From = as.character(From)) %>% 
                                                                                                                                            dplyr::rename(To = From,
                                                                                                                                                          To.total.count = From.total.count)) %>% 
  left_join(neuropeptides.biased.df %>% 
              dplyr::select(c(Keep,
                              Bias)) %>% 
              dplyr::rename(To = Keep)) %>% 
  mutate(Bias = ifelse(is.na(Bias),
                       'No bias',
                       Bias)) %>% 
  mutate(Label = ifelse(To.total.count < 20,
                        NA,
                        To))

# create igraph
#circle
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.dom.scaled <- graph_from_data_frame(d = burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.dom.scaled %>%                                                                                                             filter(!is.na(Freq)),
                                                                                                                                              vertices = burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.vertices.dom.scaled %>% 
                                                                                                                                                filter(!is.na(Label)),
                                                                                                                                              directed = TRUE)

#ggraph
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.dom.scaled %>% 
  ggraph(layout = 'linear',
         circular = TRUE,
         sort.by = To.total.count) +
  geom_edge_parallel(aes(width = Freq,
                         alpha = Freq,
                         start_cap = label_rect(node1.name),
                         end_cap = label_rect(node2.name)), 
                     arrow = arrow(length = unit(4, 
                                                 'mm')),
                     color = 'dimgrey') +
  # geom_node_point(aes(size = To.total.count),
  #                 shape = 21,
  #                 fill = "#60bb46") +
  geom_node_label(aes(label = Label,
                      size = To.total.count),
                  repel = F) +
  ggtitle('Dominant') +
  theme_classic()  +
  theme(plot.title = element_text(hjust = 0.5),
        title =element_text(size=10, face='bold'),
        axis.line =element_blank(),
        axis.title =element_blank(),
        axis.ticks  =element_blank(),
        axis.text =element_blank()) +
  scale_edge_width(range = c(1, 4))+ 
  scale_size(range = c(2,10)) +
  theme(legend.position="none")+
  xlim(-1.15, 1.15) +
  ylim(-1.15, 1.15)
ggsave('neuropeptides/comparison_new/network/Neuropeptides network all filter dom scaled bias circle presentation.png',
       width = 6.65,
       height = 6.65)

### sub neuropeptides, filter 10%, scaled
## create igraph dataframes
# create edge list
# filter out Freq less than 10
# filter out if cell count is below 10
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.sub.scaled = sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder %>% 
  filter(Freq != 100) %>%
  dplyr::rename(From = Var1,
                To = Var2,
                From.total.count = Var1.total.count,
                To.total.count = Var2.total.count) %>% 
  mutate(Freq = ifelse(Freq < 10,
                       NA,
                       Freq)) %>% 
  mutate(Freq = ifelse(From.total.count < 20 | To.total.count < 20,
                       NA,
                       Freq)) 

# create vertices
# need to get all nodes listed
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.vertices.sub.scaled = full_join(burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.sub.scaled %>% 
                                                                                                                                            dplyr::select(c(To,
                                                                                                                                                            To.total.count)) %>% 
                                                                                                                                            distinct() %>% 
                                                                                                                                            mutate(To = as.character(To)),
                                                                                                                                          burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.sub.scaled %>% 
                                                                                                                                            dplyr::select(c(From,
                                                                                                                                                            From.total.count)) %>% 
                                                                                                                                            distinct() %>% 
                                                                                                                                            mutate(From = as.character(From)) %>% 
                                                                                                                                            dplyr::rename(To = From,
                                                                                                                                                          To.total.count = From.total.count)) %>% 
  left_join(neuropeptides.biased.df %>% 
              dplyr::select(c(Keep,
                              Bias)) %>% 
              dplyr::rename(To = Keep)) %>% 
  mutate(Bias = ifelse(is.na(Bias),
                       'No bias',
                       Bias)) %>% 
  mutate(Label = ifelse(To.total.count < 20,
                        NA,
                        To))

# create igraph
#circle
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.sub.scaled  <- graph_from_data_frame(d = burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.edgelist.sub.scaled %>%                                                                                                             filter(!is.na(Freq)) %>% 
                                mutate(Status = 'sub'),
                              vertices = burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.vertices.sub.scaled %>% 
                                filter(!is.na(Label)) %>%
                                dplyr::rename(To.total.count.sub = To.total.count) %>% 
                                full_join(burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.vertices.dom.scaled %>% 
                                            filter(!is.na(Label))),
                              directed = TRUE)


#ggraph
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.network.overlap.percent.reorder.sub.scaled %>% 
  ggraph(layout = 'linear',
         circular = TRUE,
         sort.by = To.total.count) +
  geom_edge_parallel(aes(width = Freq,
                         alpha = Freq,
                         start_cap = label_rect(node1.name),
                         end_cap = label_rect(node2.name)), 
                     arrow = arrow(length = unit(4, 
                                                 'mm')),
                     color = 'dimgrey') +
  # geom_node_point(aes(size = To.total.count),
  #                 shape = 21,
  #                 fill = "#60bb46") +
  geom_node_label(aes(label = Label,
                      size = To.total.count.sub,
                      alpha = is.na(To.total.count.sub)),
                  repel = F) +
  ggtitle('Subordinate') +
  theme_classic()  +
  theme(plot.title = element_text(hjust = 0.5),
        title =element_text(size=10, face='bold'),
        axis.line =element_blank(),
        axis.title =element_blank(),
        axis.ticks  =element_blank(),
        axis.text =element_blank()) +
  scale_edge_width(range = c(1, 4))+ 
  scale_size(range = c(2,10)) +
  scale_alpha_manual(values = c(1,0)) +
  theme(legend.position="none") +
  xlim(-1.25, 1.25) +
  ylim(-1.25, 1.25)
ggsave('neuropeptides/comparison_new/network/Neuropeptides network all filter sub scaled bias circle presentation.png',
       width = 6.65,
       height = 6.65)
















#### covariance matrix ####
### compare neuropeptide cov matrix for doms and subs
#https://cran.r-project.org/web/packages/equalCovs/index.html
#https://cran.r-project.org/web/packages/HDtest/index.html

## could try?: https://doi.org/10.1080/01621459.2012.758041

## use presence absence matrix
## doms
dom.neuropeptides.cov = cov(dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix)

## sub
sub.neuropeptides.cov = cov(sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix)

## filter out genes with 0 across both matrices
## get list of genes to remove
# dom
dom.neuropeptides.cov.remove.genes = dom.neuropeptides.cov %>%
  as.data.frame() %>% 
  filter(rowSums(across(where(is.numeric)))==0) %>%
  droplevels() %>% 
  rownames()

# sub
sub.neuropeptides.cov.remove.genes = sub.neuropeptides.cov %>%
  as.data.frame() %>% 
  filter(rowSums(across(where(is.numeric)))==0) %>%
  droplevels() %>% 
  rownames()

# find genes to remove from both datasets
all.neuropeptides.cov.remove.genes = intersect(dom.neuropeptides.cov.remove.genes,
      sub.neuropeptides.cov.remove.genes)

## subset matrices to genes with data
# dom
dom.neuropeptides.cov = dom.neuropeptides.cov %>% 
  as.data.frame() %>% 
  rownames_to_column('genes') %>% 
  filter(!(genes %in% all.neuropeptides.cov.remove.genes)) %>% 
  select(-all_of(c(all.neuropeptides.cov.remove.genes,
            'genes'))) %>% 
  as.matrix()

# sub
sub.neuropeptides.cov = sub.neuropeptides.cov %>% 
  as.data.frame() %>% 
  rownames_to_column('genes') %>% 
  filter(!(genes %in% all.neuropeptides.cov.remove.genes)) %>% 
  select(-all_of(c(all.neuropeptides.cov.remove.genes,
                   'genes'))) %>% 
  as.matrix()


## stats
library(equalCovs)
equalCovs(dom.neuropeptides.cov, sub.neuropeptides.cov, 8181, 6411)

# try function from HDtest?
library(HDtest)
equalCovs(dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix %>% 
            as.matrix(), 
          sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix %>% 
            as.matrix(), 
          0.05, 
          'DomvsSub.presence')

## use count values not presence/absence? 
dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.value = GetAssayData(subset(burtoni.snseq.combined.sct.all.neurons.recluster,
                    idents = 'dom_burtoni_snseq'),
             assay = 'SCT',
             slot = 'counts') %>% 
  as.data.frame() %>% 
  filter(rownames(burtoni.snseq.combined.sct.all.neurons.recluster@assays$SCT@counts) %in% neuropeptides.list.gene.names) %>% 
  t()
