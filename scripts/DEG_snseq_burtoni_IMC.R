#### Burtoni snseq seurat analysis
### DEG analysis
###Note: Seurat requires R version > 4
## use lambcomp1 to run R with command 
# > R-4.0.3

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
# install.packages("ggalluvial")

#load libraries
library(Seurat)
library(patchwork)
library(clustree)
library(pheatmap)
library(DEsingle)
library(dendextend)
library(tidyverse)
library(ggalluvial)


#### functions ####
### create function to get dotplot data
DotPlot.data = function (object, assay = NULL, features, cols = c("lightgrey", 
                                                                  "blue"), col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, 
                         idents = NULL, group.by = NULL, split.by = NULL, cluster.idents = FALSE, 
                         scale = TRUE, scale.by = "radius", scale.min = NA, scale.max = NA) 
{
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  scale.func <- switch(EXPR = scale.by, size = scale_size, 
                       radius = scale_radius, stop("'scale.by' must be either 'size' or 'radius'"))
  feature.groups <- NULL
  if (is.list(features) | any(!is.na(names(features)))) {
    feature.groups <- unlist(x = sapply(X = 1:length(features), 
                                        FUN = function(x) {
                                          return(rep(x = names(x = features)[x], each = length(features[[x]])))
                                        }))
    if (any(is.na(x = feature.groups))) {
      warning("Some feature groups are unnamed.", call. = FALSE, 
              immediate. = TRUE)
    }
    features <- unlist(x = features)
    names(x = feature.groups) <- features
  }
  cells <- unlist(x = CellsByIdentities(object = object, idents = idents))
  data.features <- FetchData(object = object, vars = features, 
                             cells = cells)
  data.features$id <- if (is.null(x = group.by)) {
    Idents(object = object)[cells, drop = TRUE]
  }
  else {
    object[[group.by, drop = TRUE]][cells, drop = TRUE]
  }
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  if (!is.null(x = split.by)) {
    splits <- object[[split.by, drop = TRUE]][cells, drop = TRUE]
    data.features$id <- paste(data.features$id, splits, sep = "_")
    unique.splits <- unique(x = splits)
    id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), 
                        "_", rep(x = unique(x = splits), times = length(x = id.levels)))
  }
  data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
    data.use <- data.features[data.features$id == ident, 
                              1:(ncol(x = data.features) - 1), drop = FALSE]
    avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
      return(mean(x = expm1(x = x)))
    })
    pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, 
                     threshold = 0)
    return(list(avg.exp = avg.exp, pct.exp = pct.exp))
  })
  names(x = data.plot) <- unique(x = data.features$id)
  if (cluster.idents) {
    mat <- do.call(what = rbind, args = lapply(X = data.plot, 
                                               FUN = unlist))
    mat <- scale(x = mat)
    id.levels <- id.levels[hclust(d = dist(x = mat))$order]
  }
  data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
    data.use <- as.data.frame(x = data.plot[[x]])
    data.use$features.plot <- rownames(x = data.use)
    data.use$id <- x
    return(data.use)
  })
  data.plot <- do.call(what = "rbind", args = data.plot)
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }
  ngroup <- length(x = levels(x = data.plot$id))
  if (ngroup == 1) {
    scale <- FALSE
    warning("Only one identity present, the expression values will be not scaled", 
            call. = FALSE, immediate. = TRUE)
  }
  else if (ngroup < 5 & scale) {
    warning("Scaling data with a low number of groups may produce misleading results", 
            call. = FALSE, immediate. = TRUE)
  }
  avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot), 
                           FUN = function(x) {
                             data.use <- data.plot[data.plot$features.plot == 
                                                     x, "avg.exp"]
                             if (scale) {
                               data.use <- scale(x = data.use)
                               data.use <- MinMax(data = data.use, min = col.min, 
                                                  max = col.max)
                             }
                             else {
                               data.use <- log1p(x = data.use)
                             }
                             return(data.use)
                           })
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))

  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(x = data.plot$features.plot, 
                                    levels = features)
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  
  return(data.plot)
}
#### load data ####
### single cell data
#load('burtoni.snseq.combined.sct.all.RData')
load('burtoni.snseq.combined.sct.RData')

### scsorter data
load("burtoni.scsorter.data.scsort.output.RData")

### load souporcell data
burtoni.souporcell.filtered = read.csv('../souporcell/burtoni.souporcell.filtered.csv')

### load sctype.hypo data
burtoni.sctypemarkers.hypo = read.csv('./sctype.hypo/sctypemarkers.hypo.csv')

### Add metadata
## cell type
burtoni.snseq.combined.sct.all = AddMetaData(
  object = burtoni.snseq.combined.sct.all,
  metadata = burtoni.scsorter.data.scsort.output %>% 
    select(Cell.type,
           Cell.id) %>% 
    column_to_rownames(var = "Cell.id"),
  col.name = 'Cell.type'
)

## genotype
burtoni.snseq.combined.sct.all = AddMetaData(
  object = burtoni.snseq.combined.sct.all,
  metadata = burtoni.souporcell.filtered %>% 
    select(Genotype.id,
           Cell.id) %>% 
    column_to_rownames(var = "Cell.id"),
  col.name = 'Genotype.id'
)

## sctype.hypo
burtoni.snseq.combined.sct.all = AddMetaData(
  object = burtoni.snseq.combined.sct.all,
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

# #### AR ####
# #ENSONIG00000012854
# #ara
# VlnPlot(burtoni.snseq.combined.sct.all,
#         features = c("ENSONIG00000012854"),
#         split.by = "orig.ident", 
#         group.by = 'Cell.type') +
#   ylim(2,25)
# 
# VlnPlot(burtoni.snseq.combined.sct.all,
#         features = c("ara"),
#         split.by = "orig.ident", 
#         group.by = 'Cell.type') +
#   ylim(2,25)
# 
# 
# FeaturePlot(object = burtoni.snseq.combined.sct.all, 
#             features = c("ara", 
#                          "ENSONIG00000012854"),
#             pt.size = 0.5,
#             min.cutoff = 2,
#             max.cutoff = 3,
#             order = TRUE,
#             blend = T)
# 
# ##expression level
# burtoni.snseq.combined.sct.all.ara = full_join(full_join(burtoni.snseq.combined.sct.all@reductions$umap@cell.embeddings %>% 
#                                                                                         as.data.frame() %>% 
#                                                                                         rownames_to_column("Cell.id"),
#                                                      burtoni.snseq.combined.sct.all@meta.data %>% 
#                                                                                         rownames_to_column("Cell.id")),
#                                            burtoni.snseq.combined.sct.all@assays$integrated@scale.data %>% 
#                                                                               as.data.frame() %>% 
#                                                                               filter(rownames(burtoni.snseq.combined.sct.all@assays$integrated@scale.data) %in% c('ENSONIG00000012854',
#                                                                                                                                                                                         'ara')) %>% 
#                                                                               t() %>% as.data.frame() %>% 
#                                                                               rownames_to_column('Cell.id')) %>% 
#   full_join(burtoni.snseq.combined.sct.all@reductions$pca@cell.embeddings %>% 
#               as.data.frame() %>% 
#               select(PC_1,
#                      PC_2,
#                      PC_3,
#                      PC_4,
#                      PC_5) %>% 
#               rownames_to_column("Cell.id")) %>% 
#   mutate(AR.cell = ifelse(ENSONIG00000012854 > 2 & ara > 2,
#                            'both',
#                            ifelse(ara > 2,
#                                   'ara',
#                                   ifelse(ENSONIG00000012854 > 2,
#                                          'ENSONIG00000012854',
#                                          NA))))
# 
# # seperation by cell type
# table(burtoni.snseq.combined.sct.all.ara %>% 
#         filter(!is.na(AR.cell)) %>% 
#         select(Cell.type,
#                orig.ident,
#                AR.cell))
# 
# table(burtoni.snseq.combined.sct.all.ara %>% 
#         filter(!is.na(AR.cell)) %>% 
#         select(              orig.ident,
#                AR.cell))
# 
# 
# #### gnrh ####
# ## gnrh
# #avp vs oxt
# ## umap plots
# FeaturePlot(object = burtoni.snseq.combined.sct.all, 
#             features = c("ENSONIG00000011023"),
#             pt.size = 0.5,
#             min.cutoff = 4,
#             max.cutoff = 5,
#             order = TRUE)
# ggsave('gnrh/gnrh.umap.png',
#        width = 20,
#        height = 20)
# #across samples
# FeaturePlot(object = burtoni.snseq.combined.sct.all, 
#             features = c("ENSONIG00000011023"),
#             pt.size = 0.5,
#             min.cutoff = 4,
#             max.cutoff = 6,
#             order = TRUE,
#             split.by = "orig.ident")
# ggsave('gnrh/gnrh.subvsdom.umap.png',
#        width = 10,
#        height = 10)
# 
# ##vinplot
# #gnrh
# VlnPlot(burtoni.snseq.combined.sct.all, 
#         features = c("ENSONIG00000011023"),
#         split.by = "orig.ident") +
#   ylim(4,25)
# ggsave('gnrh/gnrh.subvsdom.vlnplot.png',
#        width = 10,
#        height = 10)
# 
# ##check number of cells expressing gene above threshold
# #gnrh
# #351
# sum(GetAssayData(object = burtoni.snseq.combined.sct.all, 
#                  slot = "data")["ENSONIG00000011023",]>4)
# 
# ### subset data to relevant cells
# burtoni.snseq.combined.sct.all.gnrh = subset(burtoni.snseq.combined.sct.all,
#                                              subset = ENSONIG00000011023 > 4 )
# 
# ## run PCA, UMAP, and cluster 
# #use 0.4 resolution
# burtoni.snseq.combined.sct.all.gnrh.recluster = burtoni.snseq.combined.sct.all.gnrh %>% 
#   RunPCA() %>%
#   FindNeighbors(dims = 1:15) %>%
#   RunUMAP(dims = 1:15) %>%
#   FindClusters(resolution = 0.6)
# 
# ## check cells per group per cluster
# #get table
# Gnrh.cells.per.cluster.table = table(burtoni.snseq.combined.sct.all.gnrh.recluster@active.ident, 
#                                      burtoni.snseq.combined.sct.all.gnrh.recluster@meta.data$orig.ident) %>% 
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
# Gnrh.cells.per.cluster.table %>% 
#   as.data.frame() %>% 
#   mutate(cluster.id = as.numeric(cluster.id)) %>% 
#   ggplot(aes(x = cluster.id,
#              y = count,
#              group = orig.ident,
#              color = orig.ident)) +
#   geom_point() +
#   theme_bw()
# ggsave('gnrh/Cells.per.cluster.DomvsSub.gnrh.png',
#        width = 10,
#        height = 10)
# 
# #percent
# Gnrh.cells.per.cluster.table %>% 
#   as.data.frame() %>% 
#   mutate(cluster.id = as.numeric(cluster.id)) %>% 
#   ggplot(aes(x = cluster.id,
#              y = percentage,
#              group = orig.ident,
#              color = orig.ident)) +
#   geom_point() +
#   theme_bw()
# ggsave('gnrh/Cells.per.cluster.DomvsSub.gnrh.percentage.png',
#        width = 10,
#        height = 10)
# 
# ### graph 
# ##dom vs sub
# DimPlot(burtoni.snseq.combined.sct.all.gnrh.recluster, 
#         reduction = "umap", 
#         group.by = "orig.ident")
# ggsave('gnrh/DomVsSub.dimplot.gnrh.all.png')
# 
# ## clusters
# DimPlot(burtoni.snseq.combined.sct.all.gnrh.recluster, 
#         reduction = "umap", 
#         label = TRUE,
#         repel = TRUE)
# ggsave('gnrh/Clusters.dimplot.gnrh.all.png',
#        width = 10,
#        height = 10)
# 
# ###dotplot
# ## cell type
# DotPlot(burtoni.snseq.combined.sct.all.gnrh.recluster, 
#         features = burtoni.snseq.combined.sct.all.gnrh.recluster@meta.data %>% 
#           select(ends_with(".score1")) %>% 
#           colnames(), 
#         cols = c("grey", 
#                  "red"), 
#         dot.scale = 8,
#         col.min = 1,
#         dot.min = .4) + 
#   RotatedAxis()
# ggsave('gnrh/gnrh.celltype.dotplot.png',
#        width = 10,
#        height = 10)
# 
# 
# 
# ### clustree
# # cluster across resolutions
# burtoni.snseq.combined.sct.all.gnrh.recluster.clustree <- Seurat::FindClusters(object = burtoni.snseq.combined.sct.all.gnrh.recluster, 
#                                                                                resolution = resolution.range.reduced)
# #check data
# head(burtoni.snseq.combined.sct.all.gnrh.recluster.clustree[[]])
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
# clustree(burtoni.snseq.combined.sct.all.gnrh.recluster.clustree, 
#          prefix = "integrated_snn_res.",
#          node_colour = 'cluster',
#          node_size_range = c(10,20),
#          scale_node_text = TRUE) +
#   scale_edge_color_continuous(low = "black", 
#                               high = "black") +
#   scale_color_manual(values = presentation.color)+
#   theme(legend.position = "bottom")
# ggsave('gnrh/gnrh.clustree.png',
#        width = 10,
#        height = 10)
# 
# ###pca
# # clusters
# DimPlot(burtoni.snseq.combined.sct.all.gnrh.recluster.clustree, 
#         reduction = "pca", 
#         label = TRUE,
#         repel = TRUE)
# ggsave('gnrh/gnrh.pca.PC1.PC2.png',
#        width = 10,
#        height = 10)
# #PC 2 and 3
# DimPlot(burtoni.snseq.combined.sct.all.gnrh.recluster.clustree, 
#         reduction = "pca", 
#         label = TRUE,
#         repel = TRUE,
#         dims = c(3,2))
# ggsave('gnrh/gnrh.pca.PC2.PC3.png',
#        width = 10,
#        height = 10)
# 
# 
# #heatmap PCA
# DimHeatmap(burtoni.snseq.combined.sct.all.gnrh.recluster.clustree, 
#            dims = 1:6, 
#            cells = 500, 
#            balanced = TRUE)
# ggsave('gnrh/gnrh.pca.heatmap.png',
#        width = 10,
#        height = 10)
# #loadings
# VizDimLoadings(burtoni.snseq.combined.sct.all.gnrh.recluster.clustree, 
#                dims = 1:6, 
#                reduction = "pca")
# ggsave('gnrh/gnrh.pca.loadings.png',
#        width = 10,
#        height = 10)
# 
# 
# ### identify marker genes
# ## all
# burtoni.snseq.combined.sct.all.avp.oxt.recluster.clustree.markers <- FindAllMarkers(burtoni.snseq.combined.sct.all.avp.oxt.recluster.clustree, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# 
# #### sst ####
# ## sst
# #avp vs oxt
# ## umap plots
# FeaturePlot(object = burtoni.snseq.combined.sct.all, 
#             features = c("ENSONIG00000033642"),
#             pt.size = 0.5,
#             min.cutoff = 4,
#             max.cutoff = 5,
#             order = TRUE)
# ggsave('sst/sst.umap.png',
#        width = 20,
#        height = 20)
# #across samples
# FeaturePlot(object = burtoni.snseq.combined.sct.all, 
#             features = c("ENSONIG00000033642"),
#             pt.size = 0.5,
#             min.cutoff = 4,
#             max.cutoff = 6,
#             order = TRUE,
#             split.by = "orig.ident")
# ggsave('sst/sst.subvsdom.umap.png',
#        width = 10,
#        height = 10)
# 
# ##vinplot
# #sst
# VlnPlot(burtoni.snseq.combined.sct.all, 
#         features = c("ENSONIG00000033642"),
#         split.by = "orig.ident") +
#   ylim(4,25)
# ggsave('sst/sst.subvsdom.vlnplot.png',
#        width = 10,
#        height = 10)
# 
# ##check number of cells expressing gene above threshold
# #sst
# #334
# sum(GetAssayData(object = burtoni.snseq.combined.sct.all, 
#                  slot = "data")["ENSONIG00000033642",]>4)
# 
# ### subset data to relevant cells
# burtoni.snseq.combined.sct.all.sst = subset(burtoni.snseq.combined.sct.all,
#                                             subset = ENSONIG00000033642 > 4 )
# 
# ## run PCA, UMAP, and cluster 
# #use 0.4 resolution
# burtoni.snseq.combined.sct.all.sst.recluster = burtoni.snseq.combined.sct.all.sst %>% 
#   RunPCA() %>%
#   FindNeighbors(dims = 1:15) %>%
#   RunUMAP(dims = 1:15) %>%
#   FindClusters(resolution = 0.4)
# 
# ## check cells per group per cluster
# #get table
# sst.cells.per.cluster.table = table(burtoni.snseq.combined.sct.all.sst.recluster@active.ident, 
#                                     burtoni.snseq.combined.sct.all.sst.recluster@meta.data$orig.ident) %>% 
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
# sst.cells.per.cluster.table %>% 
#   as.data.frame() %>% 
#   mutate(cluster.id = as.numeric(cluster.id)) %>% 
#   ggplot(aes(x = cluster.id,
#              y = count,
#              group = orig.ident,
#              color = orig.ident)) +
#   geom_point() +
#   theme_bw()
# ggsave('sst/Cells.per.cluster.DomvsSub.sst.png',
#        width = 10,
#        height = 10)
# 
# #percent
# sst.cells.per.cluster.table %>% 
#   as.data.frame() %>% 
#   mutate(cluster.id = as.numeric(cluster.id)) %>% 
#   ggplot(aes(x = cluster.id,
#              y = percentage,
#              group = orig.ident,
#              color = orig.ident)) +
#   geom_point() +
#   theme_bw()
# ggsave('sst/Cells.per.cluster.DomvsSub.sst.percentage.png',
#        width = 10,
#        height = 10)
# 
# ### graph 
# ##dom vs sub
# DimPlot(burtoni.snseq.combined.sct.all.sst.recluster, 
#         reduction = "umap", 
#         group.by = "orig.ident")
# ggsave('sst/DomVsSub.dimplot.sst.all.png')
# 
# ## clusters
# DimPlot(burtoni.snseq.combined.sct.all.sst.recluster, 
#         reduction = "umap", 
#         label = TRUE,
#         repel = TRUE)
# ggsave('sst/Clusters.dimplot.sst.all.png',
#        width = 10,
#        height = 10)
# 
# ###dotplot
# ## cell type
# DotPlot(burtoni.snseq.combined.sct.all.sst.recluster, 
#         features = burtoni.snseq.combined.sct.all.sst.recluster@meta.data %>% 
#           select(ends_with(".score1")) %>% 
#           colnames(), 
#         cols = c("grey", 
#                  "red"), 
#         dot.scale = 8,
#         col.min = 1,
#         dot.min = .4) + 
#   RotatedAxis()
# ggsave('sst/sst.celltype.dotplot.png',
#        width = 10,
#        height = 10)
# 
# 
# 
# ### clustree
# # cluster across resolutions
# burtoni.snseq.combined.sct.all.sst.recluster.clustree <- Seurat::FindClusters(object = burtoni.snseq.combined.sct.all.sst.recluster, 
#                                                                               resolution = resolution.range.reduced)
# #check data
# head(burtoni.snseq.combined.sct.all.sst.recluster.clustree[[]])
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
# clustree(burtoni.snseq.combined.sct.all.sst.recluster.clustree, 
#          prefix = "integrated_snn_res.",
#          node_colour = 'cluster',
#          node_size_range = c(10,20),
#          scale_node_text = TRUE) +
#   scale_edge_color_continuous(low = "black", 
#                               high = "black") +
#   scale_color_manual(values = presentation.color)+
#   theme(legend.position = "bottom")
# ggsave('sst/sst.clustree.png',
#        width = 10,
#        height = 10)
# 
# ###pca
# # clusters
# DimPlot(burtoni.snseq.combined.sct.all.sst.recluster.clustree, 
#         reduction = "pca", 
#         label = TRUE,
#         repel = TRUE)
# ggsave('sst/sst.pca.PC1.PC2.png',
#        width = 10,
#        height = 10)
# #PC 2 and 3
# DimPlot(burtoni.snseq.combined.sct.all.sst.recluster.clustree, 
#         reduction = "pca", 
#         label = TRUE,
#         repel = TRUE,
#         dims = c(3,2))
# ggsave('sst/sst.pca.PC2.PC3.png',
#        width = 10,
#        height = 10)
# 
# 
# #heatmap PCA
# DimHeatmap(burtoni.snseq.combined.sct.all.sst.recluster.clustree, 
#            dims = 1:6, 
#            cells = 500, 
#            balanced = TRUE)
# ggsave('sst/sst.pca.heatmap.png',
#        width = 10,
#        height = 10)
# #loadings
# VizDimLoadings(burtoni.snseq.combined.sct.all.sst.recluster.clustree, 
#                dims = 1:6, 
#                reduction = "pca")
# ggsave('sst/sst.pca.loadings.png',
#        width = 10,
#        height = 10)
# 
# 
# 
# #### Esr1 neurons ####
# ### esr1 neurons
# #esr1
# VlnPlot(burtoni.snseq.combined.sct.all.neurons, 
#         features = c("esr1"),
#         group.by = "orig.ident")
# ggsave('neuropeptides/esr1/esr1.neurons.DomvsSub.png',
#        width = 10,
#        height = 10)
# 
# ### subset data to relevant cells
# burtoni.snseq.combined.sct.all.esr1.neurons = burtoni.snseq.combined.sct.all.neurons
# 
# #set idents
# Idents(object = burtoni.snseq.combined.sct.all.esr1.neurons) <- "Cell.type"
# 
# #subset
# burtoni.snseq.combined.sct.all.esr1.neurons = subset(burtoni.snseq.combined.sct.all.esr1.neurons,
#                                                      subset = esr1 > 2)
# 
# 
# ## run PCA, UMAP, and cluster 
# #use 0.4 resolution
# burtoni.snseq.combined.sct.all.esr1.neurons.recluster = burtoni.snseq.combined.sct.all.esr1.neurons %>% 
#   RunPCA() %>%
#   FindNeighbors(dims = 1:15) %>%
#   RunUMAP(dims = 1:15) %>%
#   FindClusters(resolution = 0.8)
# 
# ### clustree
# # cluster across resolutions
# burtoni.snseq.combined.sct.all.esr1.neurons.clustree <- Seurat::FindClusters(object = burtoni.snseq.combined.sct.all.esr1.neurons.recluster, 
#                                                                              resolution = resolution.range.reduced)
# 
# 
# #clustree
# clustree(burtoni.snseq.combined.sct.all.esr1.neurons.clustree, 
#          prefix = "integrated_snn_res.",
#          node_colour = 'cluster',
#          node_size_range = c(10,20),
#          scale_node_text = TRUE) +
#   scale_edge_color_continuous(low = "black", 
#                               high = "black") +
#   theme(legend.position = "bottom")
# ggsave('neuropeptides/esr1/esr1.neurons.clustree.png',
#        width = 10,
#        height = 10)
# 
# 
# ### graph 
# ##dom vs sub
# DimPlot(burtoni.snseq.combined.sct.all.esr1.neurons.recluster, 
#         reduction = "umap", 
#         group.by = "orig.ident")
# ggsave('neuropeptides/esr1/DomVsSub.dimplot.esr1.neurons.all.png')
# 
# ## clusters
# DimPlot(burtoni.snseq.combined.sct.all.esr1.neurons.recluster, 
#         reduction = "umap", 
#         label = TRUE,
#         repel = TRUE,
#         pt.size = 5)
# ggsave('neuropeptides/esr1/Clusters.dimplot.esr1.neurons.all.png',
#        width = 10,
#        height = 10)
# 
# ## cell type
# DimPlot(burtoni.snseq.combined.sct.all.esr1.neurons.recluster, 
#         reduction = "umap", 
#         group.by = "Cell.type")
# ggsave('neuropeptides/esr1/Celltype.dimplot.esr1.neurons.all.png')
# 
# 
# 
# ##expression level
# burtoni.snseq.combined.sct.all.esr1.neurons.recluster.expression = full_join(full_join(burtoni.snseq.combined.sct.all.esr1.neurons.recluster@reductions$umap@cell.embeddings %>% 
#                                                                                          as.data.frame() %>% 
#                                                                                          rownames_to_column("Cell.id"),
#                                                                                        burtoni.snseq.combined.sct.all.esr1.neurons.recluster@meta.data %>% 
#                                                                                          rownames_to_column("Cell.id")),
#                                                                              burtoni.snseq.combined.sct.all.esr1.neurons.recluster@assays$integrated@scale.data %>% 
#                                                                                as.data.frame() %>% 
#                                                                                filter(rownames(burtoni.snseq.combined.sct.all.esr1.neurons.recluster@assays$integrated@scale.data) %in% c('avp',
#                                                                                                                                                                                           'oxt')) %>% 
#                                                                                t() %>% as.data.frame() %>% 
#                                                                                rownames_to_column('Cell.id')) %>% 
#   full_join(burtoni.snseq.combined.sct.all.esr1.neurons.recluster@reductions$pca@cell.embeddings %>% 
#               as.data.frame() %>% 
#               select(PC_1,
#                      PC_2,
#                      PC_3,
#                      PC_4,
#                      PC_5) %>% 
#               rownames_to_column("Cell.id")) %>% 
#   mutate(Oxt.cell = ifelse(oxt > 2 & avp > 2,
#                            'both',
#                            ifelse(oxt > 2,
#                                   'oxt',
#                                   'avp')))
# 
# # seperation by cell type
# table(burtoni.snseq.combined.sct.all.esr1.neurons.recluster.expression %>%
#         select(Cell.type,
#                integrated_snn_res.0.8))
# 
# # difference social status
# table(burtoni.snseq.combined.sct.all.esr1.neurons.recluster.expression %>%
#         select(integrated_snn_res.0.8,
#                orig.ident))
# 
# # difference social status
# # percent
# # 167 dom to 106 sub
# table(burtoni.snseq.combined.sct.all.esr1.neurons.recluster.expression %>%
#         select(integrated_snn_res.0.8,
#                orig.ident)) %>% 
#   as.data.frame() %>% 
#   pivot_wider(names_from = orig.ident,
#               values_from = Freq) %>% 
#   mutate(DomSubPerct = 100 * dom_burtoni_snseq/(sub_burtoni_snseq + dom_burtoni_snseq),
#          TotalCount = sub_burtoni_snseq + dom_burtoni_snseq) %>% 
#   ggplot(aes(x= DomSubPerct,
#              y = TotalCount,
#              color = integrated_snn_res.0.8)) +
#   geom_vline(xintercept = 61.2) +
#   geom_vline(xintercept = 50, 
#              linetype="dotted") +
#   geom_vline(xintercept = 72.4, 
#              linetype="dotted") +
#   geom_point(size = 10) +
#   theme_classic()+ 
#   theme(text = element_text(size = 20)) +
#   xlim(40,80)
# ggsave('neuropeptides/esr1/Clusters.DomvsSub.ratio.esr1.neurons.all.png',
#        width = 10,
#        height = 10)
# 
# 
# #### nucleobindin2 ####
# 
# ### NUCB2 cell counts
# # a paralog
# NUBC2.cell.counts = data.frame(Gene.names = c("nucb2a",
#                                               "nucb2b"),
#                                Total.nuclei.count = c(sum(GetAssayData(object = burtoni.snseq.combined.sct.all, 
#                                                                        slot = "data")["nucb2a",]>2),
#                                                       sum(GetAssayData(object = burtoni.snseq.combined.sct.all, 
#                                                                        slot = "data")["nucb2b",]>2)),
#                                Neuron.nuclei.count = c(sum(GetAssayData(object = burtoni.snseq.combined.sct.all.neurons, 
#                                                                         slot = "data")["nucb2a",]>2),
#                                                        sum(GetAssayData(object = burtoni.snseq.combined.sct.all.neurons, 
#                                                                         slot = "data")["nucb2b",]>2))
# ) %>% 
#   mutate(Percentage.neuron = 100*Neuron.nuclei.count/Total.nuclei.count)
# 
# ### graph
# ##vinplot
# ##dom vs sub
# #nucb2a
# VlnPlot(burtoni.snseq.combined.sct.all, 
#         features = c("nucb2a"),
#         group.by = "Cell.type",
#         split.by = "orig.ident") +
#   ylim(2,15)
# ggsave('nucb2/vinplot.nucb2a.cell.types.and.status.png',
#        width = 10,
#        height = 10)
# 
# #nucb2b
# VlnPlot(burtoni.snseq.combined.sct.all, 
#         features = c("nucb2b"),
#         group.by = "Cell.type",
#         split.by = "orig.ident") +
#   ylim(2,15)
# ggsave('nucb2/vinplot.nucb2b.cell.types.and.status.png',
#        width = 10,
#        height = 10)
# 
# 
# ### nucb2b neurons
# ### subset data to relevant cells
# burtoni.snseq.combined.sct.all.nucb2.neurons = burtoni.snseq.combined.sct.all
# 
# #set idents
# Idents(object = burtoni.snseq.combined.sct.all.nucb2.neurons) <- "Cell.type"
# 
# #subset
# burtoni.snseq.combined.sct.all.nucb2.neurons = subset(burtoni.snseq.combined.sct.all.nucb2.neurons,
#                                                       subset = nucb2b > 4 | nucb2a > 4, idents = c("neurons",
#                                                                                                    "excitatory",
#                                                                                                    "inhibitory"))
# 
# 
# ## run PCA, UMAP, and cluster 
# #use 0.4 resolution
# burtoni.snseq.combined.sct.all.nucb2.neurons.recluster = burtoni.snseq.combined.sct.all.nucb2.neurons %>% 
#   RunPCA() %>%
#   FindNeighbors(dims = 1:15) %>%
#   RunUMAP(dims = 1:15) %>%
#   FindClusters(resolution = 0.2)
# 
# ### clustree
# # cluster across resolutions
# burtoni.snseq.combined.sct.all.nucb2.neurons.clustree <- Seurat::FindClusters(object = burtoni.snseq.combined.sct.all.nucb2.neurons.recluster, 
#                                                                               resolution = resolution.range.reduced)
# #check data
# head(burtoni.snseq.combined.sct.all.nucb2.neurons.clustree[[]])
# #clustree
# clustree(burtoni.snseq.combined.sct.all.nucb2.neurons.clustree, 
#          prefix = "integrated_snn_res.",
#          node_colour = 'cluster',
#          node_size_range = c(10,20),
#          scale_node_text = TRUE) +
#   scale_edge_color_continuous(low = "black", 
#                               high = "black") +
#   theme(legend.position = "bottom")
# ggsave('nucb2/nucb2.neurons.clustree.png',
#        width = 10,
#        height = 10)
# 
# 
# ### graph 
# ##dom vs sub
# DimPlot(burtoni.snseq.combined.sct.all.nucb2.neurons.recluster, 
#         reduction = "umap", 
#         group.by = "orig.ident") +
#   ggtitle('Dom vs Sub nucb2 neurons UMAP')
# ggsave('nucb2/DomVsSub.dimplot.nucb2.neurons.all.png')
# 
# ## clusters
# DimPlot(burtoni.snseq.combined.sct.all.nucb2.neurons.recluster, 
#         reduction = "umap", 
#         label = TRUE,
#         repel = TRUE) +
#   ggtitle('nucb2 neurons UMAP clusters')
# ggsave('nucb2/Clusters.dimplot.nucb2.neurons.all.png',
#        width = 10,
#        height = 10)
# 
# 
# ###expression level
# burtoni.snseq.combined.sct.all.nucb2.neurons.recluster.expression = full_join(full_join(burtoni.snseq.combined.sct.all.nucb2.neurons.recluster@reductions$umap@cell.embeddings %>% 
#                                                                                           as.data.frame() %>% 
#                                                                                           rownames_to_column("Cell.id"),
#                                                                                         burtoni.snseq.combined.sct.all.nucb2.neurons.recluster@meta.data %>% 
#                                                                                           rownames_to_column("Cell.id")),
#                                                                               burtoni.snseq.combined.sct.all.nucb2.neurons.recluster@assays$integrated@scale.data %>% 
#                                                                                 as.data.frame() %>% 
#                                                                                 filter(rownames(burtoni.snseq.combined.sct.all.nucb2.neurons.recluster@assays$integrated@scale.data) %in% c('nucb2a',
#                                                                                                                                                                                             'nucb2b')) %>% 
#                                                                                 t() %>% as.data.frame() %>% 
#                                                                                 rownames_to_column('Cell.id')) %>% 
#   full_join(burtoni.snseq.combined.sct.all.nucb2.neurons.recluster@reductions$pca@cell.embeddings %>% 
#               as.data.frame() %>% 
#               select(PC_1,
#                      PC_2,
#                      PC_3,
#                      PC_4,
#                      PC_5) %>% 
#               rownames_to_column("Cell.id"))
# 
# #no seperation by cell type
# table(burtoni.snseq.combined.sct.all.nucb2.neurons.recluster.expression %>% 
#         select(Cell.type,
#                integrated_snn_res.0.2))
# 
# # difference social status
# table(burtoni.snseq.combined.sct.all.nucb2.neurons.recluster.expression %>% 
#         select(integrated_snn_res.0.2,
#                orig.ident))
# 
# ## graph nucb2 
# burtoni.snseq.combined.sct.all.nucb2.neurons.recluster.expression %>% 
#   mutate(Gene.expressed = ifelse(nucb2a >= 3,
#                                  'nucb2a',
#                                  ifelse(nucb2b >= 3,
#                                         'nucb2b',
#                                         NA))) %>% 
#   droplevels() %>% 
#   ggplot(aes(x = UMAP_1,
#              y = UMAP_2,
#              color = Gene.expressed)) +
#   geom_point() +
#   theme_classic() +
#   ggtitle('nucb2a vs nucb2b neurons UMAP')
# ggsave('nucb2/UMAP paralogs across neuron clusters.png',
#        width = 10,
#        height = 10)
# 
# 
# 
# 
# 
# #### AVP  ####
# ### avp 
# ### subset data to relevant cells
# burtoni.snseq.combined.sct.all.avp = burtoni.snseq.combined.sct.all
# 
# #subset
# burtoni.snseq.combined.sct.all.avp = subset(burtoni.snseq.combined.sct.all.avp,
#                                             subset = avp > 4)
# 
# ## run PCA, UMAP, and cluster 
# #use 0.4 resolution
# burtoni.snseq.combined.sct.all.avp.recluster = burtoni.snseq.combined.sct.all.avp %>% 
#   RunPCA() %>%
#   FindNeighbors(dims = 1:15) %>%
#   RunUMAP(dims = 1:15) %>%
#   FindClusters(resolution = 0.6)
# 
# ### clustree
# # cluster across resolutions
# burtoni.snseq.combined.sct.all.avp.clustree <- Seurat::FindClusters(object = burtoni.snseq.combined.sct.all.avp.recluster, 
#                                                                     resolution = resolution.range.reduced)
# #check data
# head(burtoni.snseq.combined.sct.all.avp.clustree[[]])
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
# clustree(burtoni.snseq.combined.sct.all.avp.clustree, 
#          prefix = "integrated_snn_res.",
#          node_colour = 'cluster',
#          node_size_range = c(10,20),
#          scale_node_text = TRUE) +
#   scale_edge_color_continuous(low = "black", 
#                               high = "black") +
#   scale_color_manual(values = presentation.color)+
#   theme(legend.position = "bottom")
# ggsave('avp.oxt/avp/avp.clustree.png',
#        width = 10,
#        height = 10)
# 
# 
# ### graph 
# ##dom vs sub
# DimPlot(burtoni.snseq.combined.sct.all.avp.recluster, 
#         reduction = "umap", 
#         group.by = "orig.ident")
# ggsave('avp.oxt/avp/DomVsSub.dimplot.avp.all.png')
# 
# ## clusters
# DimPlot(burtoni.snseq.combined.sct.all.avp.recluster, 
#         reduction = "umap", 
#         label = TRUE,
#         repel = TRUE)
# ggsave('avp.oxt/avp/Clusters.dimplot.avp.all.png',
#        width = 10,
#        height = 10)
# 
# #### AVP neurons 4 ####
# ### avp neurons
# ### subset data to relevant cells
# burtoni.snseq.combined.sct.all.avp.neurons.4 = burtoni.snseq.combined.sct.all
# 
# #set idents
# Idents(object = burtoni.snseq.combined.sct.all.avp.neurons.4) <- "Cell.type"
# 
# #subset
# burtoni.snseq.combined.sct.all.avp.neurons.4 = subset(burtoni.snseq.combined.sct.all.avp.neurons.4,
#                                                       subset = avp > 4, idents = c("neurons",
#                                                                                    "excitatory",
#                                                                                    "inhibitory"))
# 
# 
# ## run PCA, UMAP, and cluster 
# #use 0.4 resolution
# burtoni.snseq.combined.sct.all.avp.neurons.4.recluster = burtoni.snseq.combined.sct.all.avp.neurons.4 %>% 
#   RunPCA() %>%
#   FindNeighbors(dims = 1:15) %>%
#   RunUMAP(dims = 1:15) %>%
#   FindClusters(resolution = 1)
# 
# ### clustree
# # cluster across resolutions
# burtoni.snseq.combined.sct.all.avp.neurons.4.clustree <- Seurat::FindClusters(object = burtoni.snseq.combined.sct.all.avp.neurons.4.recluster, 
#                                                                               resolution = resolution.range.reduced)
# #clustree
# clustree(burtoni.snseq.combined.sct.all.avp.neurons.4.clustree, 
#          prefix = "integrated_snn_res.",
#          node_colour = 'cluster',
#          node_size_range = c(10,20),
#          scale_node_text = TRUE) +
#   scale_edge_color_continuous(low = "black", 
#                               high = "black")
# ggsave('avp.oxt/avp.4.threshold/avp.neurons.clustree.png',
#        width = 10,
#        height = 10)
# 
# 
# ### graph 
# ##dom vs sub
# DimPlot(burtoni.snseq.combined.sct.all.avp.neurons.4.recluster, 
#         reduction = "umap", 
#         group.by = "orig.ident")
# ggsave('avp.oxt/avp.4.threshold/DomVsSub.dimplot.avp.neurons.all.png')
# 
# ## clusters
# DimPlot(burtoni.snseq.combined.sct.all.avp.neurons.4.recluster, 
#         reduction = "umap", 
#         label = TRUE,
#         repel = TRUE)
# ggsave('avp.oxt/avp.4.threshold/Clusters.dimplot.avp.neurons.all.png',
#        width = 10,
#        height = 10)
# 
# 
# ## clusters
# # avp expression
# FeaturePlot(burtoni.snseq.combined.sct.all.avp.neurons.4.recluster, 
#             reduction = "umap",
#             features = c('avp'),
#             min.cutoff = 2,
#             max.cutoff = 5,
#             label = TRUE,
#             repel = TRUE)
# ggsave('avp.oxt/avp.4.threshold/Clusters.dimplot.avp.neurons.all.expression.png',
#        width = 10,
#        height = 10)
# 
# ##expression level
# burtoni.snseq.combined.sct.all.avp.neurons.4.recluster.expression = full_join(full_join(burtoni.snseq.combined.sct.all.avp.neurons.4.recluster@reductions$umap@cell.embeddings %>% 
#                                                                                           as.data.frame() %>% 
#                                                                                           rownames_to_column("Cell.id"),
#                                                                                         burtoni.snseq.combined.sct.all.avp.neurons.4.recluster@meta.data %>% 
#                                                                                           rownames_to_column("Cell.id")),
#                                                                               burtoni.snseq.combined.sct.all.avp.neurons.4.recluster@assays$integrated@scale.data %>% 
#                                                                                 as.data.frame() %>% 
#                                                                                 filter(rownames(burtoni.snseq.combined.sct.all.avp.neurons.4.recluster@assays$integrated@scale.data) %in% c('avp',
#                                                                                                                                                                                             'oxt')) %>% 
#                                                                                 t() %>% as.data.frame() %>% 
#                                                                                 rownames_to_column('Cell.id')) %>% 
#   full_join(burtoni.snseq.combined.sct.all.avp.neurons.4.recluster@reductions$pca@cell.embeddings %>% 
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
# # seperation by cell type
# table(burtoni.snseq.combined.sct.all.avp.neurons.4.recluster.expression %>% 
#         filter(Oxt.cell != 'both') %>% 
#         select(Cell.type,
#                integrated_snn_res.1))
# 
# # difference social status
# table(burtoni.snseq.combined.sct.all.avp.neurons.4.recluster.expression %>% 
#         filter(Oxt.cell != 'both') %>% 
#         select(integrated_snn_res.1,
#                orig.ident))
# 

#### AVP neurons ####
### avp neurons
### subset data to relevant cells
burtoni.snseq.combined.sct.all.avp.neurons = burtoni.snseq.combined.sct.all

#set idents
Idents(object = burtoni.snseq.combined.sct.neurons) <- "sctypemarkers.hypo"

#subset to neurons
burtoni.snseq.combined.sct.neurons = subset(burtoni.snseq.combined.sct.neurons,
                                            idents = c("C7-1: GLU",
                                                       "C7-2: GABA"))

#avp
# VlnPlot(burtoni.snseq.combined.sct.all, 
#         features = c("avp"),
#         split.by = "orig.ident")
# 
# VlnPlot(burtoni.snseq.combined.sct.all, 
#         features = c("avp"),
#         group.by = "orig.ident") +
#   ylim(0,5)


#set idents
Idents(object = burtoni.snseq.combined.sct.all.avp.neurons) <- "Cell.type"

# subset with SCT data
DefaultAssay(burtoni.snseq.combined.sct.all.avp.neurons) = 'SCT'

#subset
burtoni.snseq.combined.sct.all.avp.neurons = subset(burtoni.snseq.combined.sct.all.avp.neurons,
                                                    subset = avp >= 1, idents = c("neurons",
                                                                                  "excitatory",
                                                                                  "inhibitory"))


# cluster with integrated
DefaultAssay(burtoni.snseq.combined.sct.all.avp.neurons) = 'integrated'

## run PCA, UMAP, and cluster 
#use 0.8 resolution
burtoni.snseq.combined.sct.all.avp.neurons.recluster = burtoni.snseq.combined.sct.all.avp.neurons %>% 
  RunPCA() %>%
  FindNeighbors(dims = 1:15) %>%
  RunUMAP(dims = 1:15) %>%
  FindClusters(resolution = 0.8)

#set idents
Idents(object = burtoni.snseq.combined.sct.all.avp.neurons.recluster) <- "integrated_snn_res.0.8"

#subset
burtoni.snseq.combined.sct.all.avp.neurons.recluster = subset(burtoni.snseq.combined.sct.all.avp.neurons.recluster, 
                                                                 idents = c("0",
                                                                            "1",
                                                                            "2",
                                                                            "3",
                                                                            "4"))

# ### clustree
# # cluster across resolutions
# burtoni.snseq.combined.sct.all.avp.neurons.clustree <- Seurat::FindClusters(object = burtoni.snseq.combined.sct.all.avp.neurons.recluster, 
#                                                                             resolution = resolution.range.reduced)
# #check data
# #head(burtoni.snseq.combined.sct.all.avp.neurons.clustree[[]])
# 
# #clustree
# clustree(burtoni.snseq.combined.sct.all.avp.neurons.clustree, 
#          prefix = "integrated_snn_res.",
#          node_colour = 'cluster',
#          node_size_range = c(10,20),
#          scale_node_text = TRUE) +
#   scale_edge_color_continuous(low = "black", 
#                               high = "black") +
#   theme(legend.position = "bottom")
# ggsave('avp.oxt/avp/avp.neurons.clustree.png',
#        width = 10,
#        height = 10)

# 
# ### graph 
# ##dom vs sub
# DimPlot(burtoni.snseq.combined.sct.all.avp.neurons.recluster, 
#         reduction = "umap", 
#         group.by = "orig.ident")
# ggsave('avp.oxt/avp/DomVsSub.dimplot.avp.neurons.all.png')
# 
# ## clusters
# DimPlot(burtoni.snseq.combined.sct.all.avp.neurons.recluster, 
#         reduction = "umap", 
#         label = TRUE,
#         repel = TRUE)
# ggsave('avp.oxt/avp/Clusters.dimplot.avp.neurons.all.png',
#        width = 10,
#        height = 10)
# 
# 
# ## clusters
# # avp expression
# FeaturePlot(burtoni.snseq.combined.sct.all.avp.neurons.recluster, 
#             reduction = "umap",
#             features = c('avp'),
#             min.cutoff = 2,
#             max.cutoff = 5,
#             label = TRUE,
#             repel = TRUE)
# ggsave('avp.oxt/avp/Clusters.dimplot.avp.neurons.all.expression.png',
#        width = 10,
#        height = 10)

##expression level
burtoni.snseq.combined.sct.all.avp.neurons.recluster.expression = full_join(full_join(burtoni.snseq.combined.sct.all.avp.neurons.recluster@reductions$umap@cell.embeddings %>% 
                                                                                        as.data.frame() %>% 
                                                                                        rownames_to_column("Cell.id"),
                                                                                      burtoni.snseq.combined.sct.all.avp.neurons.recluster@meta.data %>% 
                                                                                        rownames_to_column("Cell.id")),
                                                                            burtoni.snseq.combined.sct.all.avp.neurons.recluster@assays$SCT@data %>% 
                                                                              as.data.frame() %>% 
                                                                              filter(rownames(burtoni.snseq.combined.sct.all.avp.neurons.recluster@assays$SCT@data) %in% c('avp',
                                                                                                                                                                                        'oxt')) %>% 
                                                                              t() %>% as.data.frame() %>% 
                                                                              rownames_to_column('Cell.id')) %>% 
  full_join(burtoni.snseq.combined.sct.all.avp.neurons.recluster@reductions$pca@cell.embeddings %>% 
              as.data.frame() %>% 
              select(PC_1,
                     PC_2,
                     PC_3,
                     PC_4,
                     PC_5) %>% 
              rownames_to_column("Cell.id")) %>% 
  mutate(Oxt.cell = ifelse(oxt >= 1 & avp >= 1,
                           'both',
                           ifelse(oxt >= 1,
                                  'oxt',
                                  'avp')))

# # seperation by cell type
# table(burtoni.snseq.combined.sct.all.avp.neurons.recluster.expression %>% 
#         filter(Oxt.cell != 'both') %>% 
#         select(Cell.type,
#                integrated_snn_res.0.8))
# 
# # difference social status
# table(burtoni.snseq.combined.sct.all.avp.neurons.recluster.expression %>% 
#         filter(Oxt.cell != 'both') %>% 
#         select(integrated_snn_res.0.8,
#                orig.ident))
# 
# # difference social status
# # percent
# # 347 dom to 214 sub
# table(burtoni.snseq.combined.sct.all.avp.neurons.recluster.expression %>% 
#         filter(Oxt.cell != 'both') %>% 
#         select(integrated_snn_res.0.8,
#                orig.ident)) %>% 
#   as.data.frame() %>% 
#   pivot_wider(names_from = orig.ident,
#               values_from = Freq) %>% 
#   mutate(DomSubPerct = 100 * dom_burtoni_snseq/(sub_burtoni_snseq + dom_burtoni_snseq),
#          TotalCount = sub_burtoni_snseq + dom_burtoni_snseq) %>% 
#   ggplot(aes(x= DomSubPerct,
#              y = TotalCount,
#              color = integrated_snn_res.0.8)) +
#   geom_vline(xintercept = 61.9) +
#   geom_point(size = 5) +
#   theme_classic() + 
#   theme(text = element_text(size = 20)) +
#   ggtitle('Clusters.DomvsSub.ratio.avp.neurons.all')
# ggsave('avp.oxt/avp/Clusters.DomvsSub.ratio.avp.neurons.all.png',
#        width = 10,
#        height = 10)

# burtoni.snseq.combined.sct.all.avp.neurons.recluster.expression %>%
#   ggplot(aes(x = integrated_snn_res.0.8,
#          y = avp,
#          fill = orig.ident)) +
#   geom_boxplot()










# #### DEsingle approach ####
# ###DEsingle
# ##https://miaozhun.github.io/DEsingle/
# #### Use DEsingle to calculate DEG between cell profile groups
# 
# ### calculate variable genes
# ## identify top 500 variable genes
# burtoni.snseq.combined.sct.all.avp.neurons.recluster.variable <- FindVariableFeatures(burtoni.snseq.combined.sct.all.avp.neurons.recluster, 
#                                                                                       selection.method = "vst", 
#                                                                                       nfeatures = 500, 
#                                                                                       verbose = F)
# # identify top 500 variable genes
# avp.neuron.group.topgenes <- head(VariableFeatures(burtoni.snseq.combined.sct.all.avp.neurons.recluster.variable), 
#                                   500)
# 
# ### dom vs sub
# ### cluster 0 
# ### group vector
# ## vector of factor specifies groups
# ## corresponds to columns in counts
# ## subset to groups of interest 
# avp.neuron.group.vector.0 = burtoni.snseq.combined.sct.all.avp.neurons.recluster.expression %>% 
#   filter(integrated_snn_res.0.8 == 0)
# 
# ## create vector of factor
# avp.neuron.0.vector.list = avp.neuron.group.vector.0 %>% 
#   mutate(avp.neuron.0.orig.ident = orig.ident %>% 
#            as.factor()) %>% 
#   pull(avp.neuron.0.orig.ident) %>% 
#   droplevels()
# 
# ### counts matrix
# ## raw read count matrix
# ## rows = genes, columns = cells
# # subset to groups of interest (0 and 1 of AVP neurons)
# # only keep 500 variable genes
# # set negative values to 0
# avp.neuron.0.vector.count = GetAssayData(burtoni.snseq.combined.sct.all.avp.neurons.recluster) %>% 
#   as_tibble(rownames = NA) %>% 
#   rownames_to_column('gene') %>% 
#   select(c(gene,
#            avp.neuron.group.vector.0 %>% 
#              pull(Cell.id))) %>% 
#   filter(gene %in% avp.neuron.group.topgenes) %>% 
#   column_to_rownames('gene') %>% 
#   as.matrix() %>% 
#   pmax(0)
# 
# ### Detecting the DE genes
# avp.neuron.0.results = DEsingle(counts = avp.neuron.0.vector.count, 
#                                 group = avp.neuron.0.vector.list,
#                                 parallel = TRUE)
# 
# # Dividing the DE genes into 3 categories at threshold of FDR < 0.1
# avp.neuron.0.results.classified <- DEtype(results = avp.neuron.0.results, 
#                                           threshold = 0.1)
# 
# 
# # create color groups
# avp.neuron.0.results.classified = avp.neuron.0.results.classified %>% 
#   mutate(Sig = ifelse(pvalue <= 0.01,
#                       'Sig',
#                       'Not Sig'),
#          Direction = total_mean_1 - total_mean_2,
#          Direction.type = ifelse(Direction > 0,
#                                  'up',
#                                  'down'),
#          Sig.direction = ifelse(Sig == 'Sig',
#                                 Direction.type,
#                                 Sig))
# 
# 
# ##graph volcano plot
# # create log fold change and -log pvalue
# avp.neuron.0.results.classified %>% 
#   rownames_to_column('gene') %>% 
#   mutate(logFoldChange = log(foldChange)) %>% 
#   mutate(logFoldChange = ifelse(logFoldChange < -3.5,
#                                 -3.3,
#                                 logFoldChange),
#          logFoldChange = ifelse(logFoldChange > 3,
#                                 3.3,
#                                 logFoldChange)) %>% 
#   mutate(sig.label = ifelse(Sig == 'Sig',
#                             gene,
#                             '')) %>% 
#   ggplot(aes(x = logFoldChange,
#              y = -log10(pvalue),
#              color = Sig.direction))+
#   geom_vline(xintercept = 3.2, 
#              linetype="dotted")+
#   geom_vline(xintercept = -3.2,
#              linetype="dotted") +
#   geom_point(size = 5) +
#   theme_classic() +
#   xlim(c(-3.5,3.5)) + 
#   scale_color_manual(values=c("21B9CA", 
#                               "grey", 
#                               "orange3")) +
#   geom_text(aes(label = sig.label),
#             vjust = 0, 
#             nudge_y = 0.10,
#             size = 5) +
#   theme(text = element_text(size = 20),
#         legend.position = 'none')  
# ggsave('neuropeptides/avp.oxt/avp/DEsingle.avp.neurons.0.DomvsSub.volcano.png',
#        width = 5,
#        height = 5)
# 
# ### dom vs sub
# ### cluster 2
# ### group vector
# ## vector of factor specifies groups
# ## corresponds to columns in counts
# ## subset to groups of interest 
# avp.neuron.group.vector.2 = burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.expression %>% 
#   filter(integrated_snn_res.0.2 == 2)
# 
# ## create vector of factor
# avp.neuron.2.vector.list = avp.neuron.group.vector.2 %>% 
#   mutate(avp.neuron.2.orig.ident = orig.ident %>% 
#            as.factor()) %>% 
#   pull(avp.neuron.2.orig.ident) %>% 
#   droplevels()
# 
# ### counts matrix
# ## raw read count matrix
# ## rows = genes, columns = cells
# # subset to groups of interest (0 and 1 of AVP neurons)
# # only keep 500 variable genes
# # set negative values to 0
# avp.neuron.2.vector.count = GetAssayData(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster) %>% 
#   as_tibble(rownames = NA) %>% 
#   rownames_to_column('gene') %>% 
#   select(c(gene,
#            avp.neuron.group.vector.2 %>% 
#              pull(Cell.id))) %>% 
#   filter(gene %in% avp.neuron.reduce.group.topgenes) %>% 
#   column_to_rownames('gene') %>% 
#   as.matrix() %>% 
#   pmax(0)
# 
# ### Detecting the DE genes
# avp.neuron.2.results = DEsingle(counts = avp.neuron.2.vector.count, 
#                                 group = avp.neuron.2.vector.list,
#                                 parallel = TRUE)
# 
# # Dividing the DE genes into 3 categories at threshold of FDR < 0.1
# avp.neuron.2.results.classified <- DEtype(results = avp.neuron.2.results, 
#                                           threshold = 0.1)
# 
# 
# # create color groups
# avp.neuron.2.results.classified = avp.neuron.2.results.classified %>% 
#   mutate(Sig = ifelse(pvalue <= 0.01,
#                       'Sig',
#                       'Not Sig'),
#          Direction = total_mean_1 - total_mean_2,
#          Direction.type = ifelse(Direction > 0,
#                                  'up',
#                                  'down'),
#          Sig.direction = ifelse(Sig == 'Sig',
#                                 Direction.type,
#                                 Sig))
# 
# 
# ##graph volcano plot
# # create log fold change and -log pvalue
# avp.neuron.2.results.classified %>% 
#   rownames_to_column('gene') %>% 
#   mutate(logFoldChange = log(foldChange)) %>% 
#   mutate(logFoldChange = ifelse(logFoldChange < -3.5,
#                                 -3.3,
#                                 logFoldChange),
#          logFoldChange = ifelse(logFoldChange > 3,
#                                 3.3,
#                                 logFoldChange)) %>% 
#   mutate(sig.label = ifelse(Sig == 'Sig',
#                             gene,
#                             '')) %>% 
#   ggplot(aes(x = logFoldChange,
#              y = -log10(pvalue),
#              color = Sig.direction))+
#   geom_vline(xintercept = 3.2, 
#              linetype="dotted")+
#   geom_vline(xintercept = -3.2,
#              linetype="dotted") +
#   geom_point(size = 5) +
#   theme_classic() +
#   xlim(c(-3.5,3.5)) + 
#   scale_color_manual(values=c("21B9CA", 
#                               "grey", 
#                               "orange3")) +
#   geom_text(aes(label = sig.label),
#             vjust = 0, 
#             nudge_y = 0.10,
#             size = 5) +
#   theme(text = element_text(size = 20),
#         legend.position = 'none')  
# ggsave('neuropeptides/avp.oxt/avp/DEsingle.avp.neurons.2.DomvsSub.volcano.png',
#        width = 5,
#        height = 5)
# 
# ### dom vs sub
# ### group vector
# ## vector of factor specifies groups
# ## corresponds to columns in counts
# ## subset to groups of interest 
# avp.neuron.DomvsSub.vector = burtoni.snseq.combined.sct.all.avp.neurons.recluster.expression %>% 
#   filter(integrated_snn_res.0.8 != 7,
#          integrated_snn_res.0.8 != 6,
#          integrated_snn_res.0.8 != 5) 
# 
# ## create vector of factor
# avp.neuron.DomvsSub.vector.list = avp.neuron.DomvsSub.vector %>% 
#   mutate(avp.neuron.orig.ident = orig.ident %>% 
#            as.factor()) %>% 
#   pull(avp.neuron.orig.ident) %>% 
#   droplevels()
# 
# ### counts matrix
# ## raw read count matrix
# ## rows = genes, columns = cells
# # subset to groups of interest (0 and 1 of AVP neurons)
# # only keep 500 variable genes
# # set negative values to 0
# avp.neuron.DomvsSub.vector.count = GetAssayData(burtoni.snseq.combined.sct.all.avp.neurons.recluster) %>% 
#   as_tibble(rownames = NA) %>% 
#   rownames_to_column('gene') %>% 
#   select(c(gene,
#            avp.neuron.DomvsSub.vector %>% 
#              pull(Cell.id))) %>% 
#   filter(gene %in% avp.neuron.group.topgenes) %>% 
#   column_to_rownames('gene') %>% 
#   as.matrix() %>% 
#   pmax(0)
# 
# ### Detecting the DE genes
# avp.neuron.DomvsSub.results = DEsingle(counts = avp.neuron.DomvsSub.vector.count, 
#                                        group = avp.neuron.DomvsSub.vector.list,
#                                        parallel = TRUE)
# 
# # Dividing the DE genes into 3 categories at threshold of FDR < 0.1
# avp.neuron.DomvsSub.results.classified <- DEtype(results = avp.neuron.DomvsSub.results, 
#                                                  threshold = 0.1)
# 
# 
# # create color groups
# avp.neuron.DomvsSub.results.classified = avp.neuron.DomvsSub.results.classified %>% 
#   mutate(Sig = ifelse(pvalue <= 0.01,
#                       'Sig',
#                       'Not Sig'),
#          Direction = total_mean_1 - total_mean_2,
#          Direction.type = ifelse(Direction > 0,
#                                  'up',
#                                  'down'),
#          Sig.direction = ifelse(Sig == 'Sig',
#                                 Direction.type,
#                                 Sig))
# 
# 
# ##graph volcano plot
# # create log fold change and -log pvalue
# avp.neuron.DomvsSub.results.classified %>% 
#   rownames_to_column('gene') %>% 
#   mutate(logFoldChange = log(foldChange)) %>% 
#   mutate(logFoldChange = ifelse(logFoldChange < -3.5,
#                                 -3.3,
#                                 logFoldChange),
#          logFoldChange = ifelse(logFoldChange > 3,
#                                 3.3,
#                                 logFoldChange)) %>% 
#   mutate(sig.label = ifelse(Sig == 'Sig',
#                             gene,
#                             '')) %>% 
#   ggplot(aes(x = logFoldChange,
#              y = -log10(pvalue),
#              color = Sig.direction))+
#   geom_vline(xintercept = 3.2, 
#              linetype="dotted")+
#   geom_vline(xintercept = -3.2,
#              linetype="dotted") +
#   geom_point(size = 5) +
#   theme_classic() +
#   xlim(c(-3.5,3.5)) + 
#   scale_color_manual(values=c("21B9CA", 
#                               "grey", 
#                               "orange3")) +
#   geom_text(aes(label = sig.label),
#             vjust = 0, 
#             nudge_y = 0.10,
#             size = 5) +
#   theme(text = element_text(size = 20),
#         legend.position = 'none')  
# ggsave('neuropeptides/avp.oxt/avp/DEsingle.avp.neurons.DomvsSub.volcano.png',
#        width = 5,
#        height = 5)
# 
# #TCF4, mpped1, pcdh1b, gh1, prl, COX3, ND4
# #### edgeRQLFDetRate approach ####
# #https://github.com/csoneson/conquer_comparison/blob/master/scripts/apply_edgeRQLFDetRate.R
# 
# #load edgeR
# suppressPackageStartupMessages(library(edgeR))
# 
# #create function
# run_edgeRQLFDetRate <- function(L) {
#   message("edgeRQLFDetRate")
#   session_info <- sessionInfo()
#   tryCatch({
#     timing <- system.time({
#       dge <- DGEList(L$count, group = L$condt)
#       dge <- calcNormFactors(dge)
#       cdr <- scale(colMeans(L$count > 0))
#       design <- model.matrix(~ cdr + L$condt)
#       dge <- estimateDisp(dge, design = design)
#       fit <- glmQLFit(dge, design = design)
#       qlf <- glmQLFTest(fit)
#       tt <- topTags(qlf, n = Inf)
#     })
#     
#     plotBCV(dge)
#     plotQLDisp(fit)
#     hist(tt$table$PValue, 50)
#     hist(tt$table$FDR, 50)
#     limma::plotMDS(dge, col = as.numeric(as.factor(L$condt)), pch = 19)
#     plotSmear(qlf)
#     
#     list(session_info = session_info,
#          timing = timing,
#          tt = tt,
#          df = data.frame(pval = tt$table$PValue,
#                          padj = tt$table$FDR,
#                          row.names = rownames(tt$table)))
#   }, error = function(e) {
#     "edgeRQLFDetRate results could not be calculated"
#     list(session_info = session_info)
#   })
# }
# 
# 
# ## run with all avp neurons
# #need list with L with  count and condt
# avp.neuron.DomvsSub.vector.edgeR = list(count = avp.neuron.DomvsSub.vector.count,
#                                         condt = avp.neuron.DomvsSub.vector.list)
# 
# 
# #run function  
# avp.neuron.DomvsSub.vector.edgeR.results = run_edgeRQLFDetRate(avp.neuron.DomvsSub.vector.edgeR)
# 
# # save results to dataframe
# avp.neuron.DomvsSub.vector.edgeR.results.df = avp.neuron.DomvsSub.vector.edgeR.results$tt@.Data[[1]]
# 
# #add color for significance  
# avp.neuron.DomvsSub.vector.edgeR.results.df = avp.neuron.DomvsSub.vector.edgeR.results.df %>% 
#   mutate(Sig = ifelse(FDR <= 0.05,
#                       'Sig',
#                       'Not Sig'),
#          Direction.type = ifelse(logFC > 0,
#                                  'up',
#                                  'down'),
#          Sig.direction = ifelse(Sig == 'Sig',
#                                 Direction.type,
#                                 Sig))
# 
# 
# ##graph volcano plot
# # create log fold change and -log pvalue
# avp.neuron.DomvsSub.vector.edgeR.results.df %>% 
#   rownames_to_column('gene') %>% 
#   mutate(sig.label = ifelse(Sig == 'Sig',
#                             gene,
#                             '')) %>% 
#   ggplot(aes(x = logFC,
#              y = -log10(FDR),
#              color = Sig.direction))+
#   geom_point(size = 5) +
#   theme_classic() + 
#   scale_color_manual(values=c("21B9CA", 
#                               "grey", 
#                               "orange3")) +
#   geom_text(aes(label = sig.label),
#             vjust = 0, 
#             nudge_y = 0.10,
#             size = 5) +
#   theme(text = element_text(size = 20),
#         legend.position = 'none')
# ggsave('neuropeptides/avp.oxt/avp/EdgeR.avp.neurons.DomvsSub.volcano.png',
#        width = 5,
#        height = 5)
# 
# ## run with avp cluster 0
# #need list with L with  count and condt
# avp.neuron.0.vector.edgeR = list(count = avp.neuron.0.vector.count,
#                                  condt = avp.neuron.0.vector.list)
# 
# 
# #run function  
# avp.neuron.0.vector.edgeR.results = run_edgeRQLFDetRate(avp.neuron.0.vector.edgeR)
# 
# # save results to dataframe
# avp.neuron.0.vector.edgeR.results.df = avp.neuron.0.vector.edgeR.results$tt@.Data[[1]]
# 
# #add color for significance  
# avp.neuron.0.vector.edgeR.results.df = avp.neuron.0.vector.edgeR.results.df %>% 
#   mutate(Sig = ifelse(FDR <= 0.05,
#                       'Sig',
#                       'Not Sig'),
#          Direction.type = ifelse(logFC > 0,
#                                  'up',
#                                  'down'),
#          Sig.direction = ifelse(Sig == 'Sig',
#                                 Direction.type,
#                                 Sig))
# 
# 
# ##graph volcano plot
# # create log fold change and -log pvalue
# avp.neuron.0.vector.edgeR.results.df %>% 
#   rownames_to_column('gene') %>% 
#   mutate(sig.label = ifelse(Sig == 'Sig',
#                             gene,
#                             '')) %>% 
#   ggplot(aes(x = logFC,
#              y = -log10(FDR),
#              color = Sig.direction))+
#   geom_point(size = 5) +
#   theme_classic() + 
#   scale_color_manual(values=c("21B9CA", 
#                               "grey", 
#                               "orange3")) +
#   geom_text(aes(label = sig.label),
#             vjust = 0, 
#             nudge_y = 0.10,
#             size = 5) +
#   theme(text = element_text(size = 20),
#         legend.position = 'none')
# ggsave('neuropeptides/avp.oxt/avp/EdgeR.avp.neurons.0.DomvsSub.volcano.png',
#        width = 5,
#        height = 5)
# 
# 
# #### MASTcpmDetRate approach ####
# #https://github.com/csoneson/conquer_comparison/blob/master/scripts/apply_MASTcpmDetRate.R
# #load libraries
# suppressPackageStartupMessages(library(MAST))
# suppressPackageStartupMessages(library(edgeR))
# 
# #create function
# run_MASTcpmDetRate <- function(L) {
#   message("MAST, CPM (including detection rate)")
#   session_info <- sessionInfo()
#   tryCatch({
#     timing <- system.time({
#       stopifnot(all(names(L$condt) == colnames(L$count)))
#       grp <- L$condt
#       cdr <- scale(colMeans(L$count > 0))
#       dge <- DGEList(counts = L$count)
#       dge <- edgeR::calcNormFactors(dge)
#       cpms <- cpm(dge)
#       sca <- FromMatrix(exprsArray = log2(cpms + 1), 
#                         cData = data.frame(wellKey = names(grp), 
#                                            grp = grp, cdr = cdr))
#       zlmdata <- zlm.SingleCellAssay(~cdr + grp, sca)
#       mast <- lrTest(zlmdata, "grp")
#     })
#     
#     hist(mast[, "hurdle", "Pr(>Chisq)"], 50)
#     
#     list(session_info = session_info,
#          timing = timing,
#          res = mast,
#          df = data.frame(pval = mast[, "hurdle", "Pr(>Chisq)"],
#                          row.names = names(mast[, "hurdle", "Pr(>Chisq)"])))
#   }, error = function(e) {
#     "MASTcpmDetRate results could not be calculated"
#     list(session_info = session_info)
#   })
# }
# 
# 
# #something wrong with cpm function...
# 
# 
# ## run with all avp neurons
# #need list with L with  count and condt
# avp.neuron.DomvsSub.vector.MAST = list(count = avp.neuron.DomvsSub.vector.count,
#                                         condt = avp.neuron.DomvsSub.vector.list)
# 
# #run function  
# avp.neuron.DomvsSub.vector.MAST.results = run_MASTcpmDetRate(avp.neuron.DomvsSub.vector.MAST)
# 
# # save results to dataframe
# avp.neuron.DomvsSub.vector.MAST.results.df = avp.neuron.DomvsSub.vector.MAST.results$tt@.Data[[1]]
# 
# #add color for significance  
# avp.neuron.DomvsSub.vector.MAST.results.df = avp.neuron.DomvsSub.vector.MAST.results.df %>% 
#   mutate(Sig = ifelse(FDR <= 0.05,
#                       'Sig',
#                       'Not Sig'),
#          Direction.type = ifelse(logFC > 0,
#                                  'up',
#                                  'down'),
#          Sig.direction = ifelse(Sig == 'Sig',
#                                 Direction.type,
#                                 Sig))
# 
# 
# ##graph volcano plot
# # create log fold change and -log pvalue
# avp.neuron.DomvsSub.vector.MAST.results.df %>% 
#   rownames_to_column('gene') %>% 
#   mutate(sig.label = ifelse(Sig == 'Sig',
#                             gene,
#                             '')) %>% 
#   ggplot(aes(x = logFC,
#              y = -log10(FDR),
#              color = Sig.direction))+
#   geom_point(size = 5) +
#   theme_classic() + 
#   scale_color_manual(values=c("21B9CA", 
#                               "grey", 
#                               "orange3")) +
#   geom_text(aes(label = sig.label),
#             vjust = 0, 
#             nudge_y = 0.10,
#             size = 5) +
#   theme(text = element_text(size = 20),
#         legend.position = 'none')
# ggsave('neuropeptides/avp.oxt/avp/MAST.avp.neurons.DomvsSub.volcano.png',
#        width = 5,
#        height = 5)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# #### limmatrend approach ####
# #https://github.com/csoneson/conquer_comparison/blob/master/scripts/apply_limmatrend.R
# #load libraries
# suppressPackageStartupMessages(library(limma))
# suppressPackageStartupMessages(library(edgeR))
# 
# #create function
# run_limmatrend <- function(L) {
#   message("limmatrend")
#   session_info <- sessionInfo()
#   timing <- system.time({
#     dge <- DGEList(L$count, group = L$condt)
#     dge <- calcNormFactors(dge)
#     design <- model.matrix(~L$condt)
#     y <- new("EList")
#     y$E <- edgeR::cpm(dge, log = TRUE, prior.count = 3)
#     fit <- lmFit(y, design = design)
#     fit <- eBayes(fit, trend = TRUE, robust = TRUE)
#     tt <- topTable(fit, n = Inf, adjust.method = "BH")
#   })
#   
#   hist(tt$P.Value, 50)
#   hist(tt$adj.P.Val, 50)
#   limma::plotMDS(dge, col = as.numeric(as.factor(L$condt)), pch = 19)
#   plotMD(fit)
#   
#   list(session_info = session_info,
#        timing = timing,
#        tt = tt,
#        df = data.frame(pval = tt$P.Value,
#                        padj = tt$adj.P.Val,
#                        row.names = rownames(tt)))
# }
# 
# ## run with all avp neurons
# #need list with L with  count and condt
# avp.neuron.DomvsSub.vector.limma = list(count = avp.neuron.DomvsSub.vector.count,
#                                        condt = avp.neuron.DomvsSub.vector.list)
# 
# #run function  
# avp.neuron.DomvsSub.vector.limma.results = run_limmatrend(avp.neuron.DomvsSub.vector.limma)
# 
# # save results to dataframe
# avp.neuron.DomvsSub.vector.limma.results.df = avp.neuron.DomvsSub.vector.limma.results$tt
# 
# #add color for significance  
# avp.neuron.DomvsSub.vector.limma.results.df = avp.neuron.DomvsSub.vector.limma.results.df %>% 
#   mutate(Sig = ifelse(adj.P.Val <= 0.05,
#                       'Sig',
#                       'Not Sig'),
#          Direction.type = ifelse(logFC > 0,
#                                  'up',
#                                  'down'),
#          Sig.direction = ifelse(Sig == 'Sig',
#                                 Direction.type,
#                                 Sig))
# 
# 
# ##graph volcano plot
# # create log fold change and -log pvalue
# avp.neuron.DomvsSub.vector.limma.results.df %>% 
#   rownames_to_column('gene') %>% 
#   mutate(sig.label = ifelse(Sig == 'Sig',
#                             gene,
#                             '')) %>% 
#   ggplot(aes(x = logFC,
#              y = -log10(adj.P.Val),
#              color = Sig.direction))+
#   geom_point(size = 5) +
#   theme_classic() + 
#   scale_color_manual(values=c("21B9CA", 
#                               "grey", 
#                               "orange3")) +
#   geom_text(aes(label = sig.label),
#             vjust = 0, 
#             nudge_y = 0.10,
#             size = 5) +
#   theme(text = element_text(size = 20),
#         legend.position = 'none')
# ggsave('neuropeptides/avp.oxt/avp/limma.avp.neurons.DomvsSub.volcano.png',
#        width = 5,
#        height = 5)
# 
# ## run with avp cluster 0
# #need list with L with  count and condt
# avp.neuron.0.vector.limma = list(count = avp.neuron.0.vector.count,
#                                  condt = avp.neuron.0.vector.list)
# 
# 
# #run function  
# avp.neuron.0.vector.limma.results = run_limmatrend(avp.neuron.0.vector.limma)
# 
# # save results to dataframe
# avp.neuron.0.vector.limma.results.df = avp.neuron.0.vector.limma.results$tt
# 
# #add color for significance  
# avp.neuron.0.vector.limma.results.df = avp.neuron.0.vector.limma.results.df %>% 
#   mutate(Sig = ifelse(adj.P.Val <= 0.05,
#                       'Sig',
#                       'Not Sig'),
#          Direction.type = ifelse(logFC > 0,
#                                  'up',
#                                  'down'),
#          Sig.direction = ifelse(Sig == 'Sig',
#                                 Direction.type,
#                                 Sig))
# 
# 
# ##graph volcano plot
# # create log fold change and -log pvalue
# avp.neuron.0.vector.limma.results.df %>% 
#   rownames_to_column('gene') %>% 
#   mutate(sig.label = ifelse(Sig == 'Sig',
#                             gene,
#                             '')) %>% 
#   ggplot(aes(x = logFC,
#              y = -log10(adj.P.Val),
#              color = Sig.direction))+
#   geom_point(size = 5) +
#   theme_classic() + 
#   scale_color_manual(values=c("21B9CA", 
#                               "grey", 
#                               "orange3")) +
#   geom_text(aes(label = sig.label),
#             vjust = 0, 
#             nudge_y = 0.10,
#             size = 5) +
#   theme(text = element_text(size = 20),
#         legend.position = 'none')
# ggsave('neuropeptides/avp.oxt/avp/limma.avp.neurons.0.DomvsSub.volcano.png',
#        width = 5,
#        height = 5)
# 
# ## run with avp cluster 2
# #need list with L with  count and condt
# avp.neuron.2.vector.limma = list(count = avp.neuron.2.vector.count,
#                                  condt = avp.neuron.2.vector.list)
# 
# 
# #run function  
# avp.neuron.2.vector.limma.results = run_limmatrend(avp.neuron.2.vector.limma)
# 
# # save results to dataframe
# avp.neuron.2.vector.limma.results.df = avp.neuron.2.vector.limma.results$tt
# 
# #add color for significance  
# avp.neuron.2.vector.limma.results.df = avp.neuron.2.vector.limma.results.df %>% 
#   mutate(Sig = ifelse(adj.P.Val <= 0.05,
#                       'Sig',
#                       'Not Sig'),
#          Direction.type = ifelse(logFC > 0,
#                                  'up',
#                                  'down'),
#          Sig.direction = ifelse(Sig == 'Sig',
#                                 Direction.type,
#                                 Sig))
# 
# 
# ##graph volcano plot
# # create log fold change and -log pvalue
# avp.neuron.2.vector.limma.results.df %>% 
#   rownames_to_column('gene') %>% 
#   mutate(sig.label = ifelse(Sig == 'Sig',
#                             gene,
#                             '')) %>% 
#   ggplot(aes(x = logFC,
#              y = -log10(adj.P.Val),
#              color = Sig.direction))+
#   geom_hline(yintercept = -log10(0.05),
#              linetype = 'dotted')+
#   geom_point(size = 5) +
#   theme_classic() + 
#   scale_color_manual(values=c(#"21B9CA", 
#                               "grey", 
#                               "orange3")) +
#   geom_text(aes(label = sig.label),
#             vjust = 0, 
#             nudge_y = 0.10,
#             size = 5) +
#   theme(text = element_text(size = 20),
#         legend.position = 'none')
# ggsave('neuropeptides/avp.oxt/avp/limma.avp.neurons.2.DomvsSub.volcano.png',
#        width = 5,
#        height = 5)
# 
# #### compare methods #### 
# ### compare to desingle
# 
# ### cluster 0
# #rename columns
# #edgeR
# avp.neuron.0.vector.edgeR.results.df.rename = avp.neuron.0.vector.edgeR.results.df
# colnames(avp.neuron.0.vector.edgeR.results.df.rename) = paste(colnames(avp.neuron.0.vector.edgeR.results.df.rename),
#                                                               "edgeR",
#                                                               sep=".")
# #DEsingle
# avp.neuron.0.results.classified.rename = avp.neuron.0.results.classified
# colnames(avp.neuron.0.results.classified.rename) = paste(colnames(avp.neuron.0.results.classified.rename),
#                                                               "DEsingle",
#                                                               sep=".")
# #limma
# avp.neuron.0.vector.limma.results.df.rename = avp.neuron.0.vector.limma.results.df
# colnames(avp.neuron.0.vector.limma.results.df.rename) = paste(colnames(avp.neuron.0.vector.limma.results.df.rename),
#                                                               "limma",
#                                                               sep=".")
# 
# #join data
# avp.neuron.0.combined.results = full_join(avp.neuron.0.vector.edgeR.results.df.rename %>% 
#                                             rownames_to_column('gene'),
#                                           avp.neuron.0.results.classified.rename %>% 
#                                             rownames_to_column('gene'),
#                                           by = c('gene')) %>% 
#   full_join(avp.neuron.0.vector.limma.results.df.rename %>% 
#               rownames_to_column('gene'))
# 
# 
# 
# #logFC
# png('neuropeptides/avp.oxt/avp/comparison/cluster0/logFC comparison DEG methods.png')
# chart.Correlation(avp.neuron.0.combined.results %>%
#                     mutate(logFC.DEsingle = log10(foldChange.DEsingle)) %>% 
#                     select(c(logFC.edgeR,
#                              logFC.limma,
#                              logFC.DEsingle)),
#                   pch=19)
# dev.off()
# 
# #logFC
# png('neuropeptides/avp.oxt/avp/comparison/cluster0/pvalue comparison DEG methods.png')
# chart.Correlation(avp.neuron.0.combined.results %>%
#                     select(c(pvalue.DEsingle,
#                              PValue.edgeR,
#                              P.Value.limma)),
#                   pch=19)
# dev.off()
# 
# ### DomvsSub
# #rename columns
# #edgeR
# avp.neuron.DomvsSub.vector.edgeR.results.df.rename = avp.neuron.DomvsSub.vector.edgeR.results.df
# colnames(avp.neuron.DomvsSub.vector.edgeR.results.df.rename) = paste(colnames(avp.neuron.DomvsSub.vector.edgeR.results.df.rename),
#                                                               "edgeR",
#                                                               sep=".")
# #DEsingle
# avp.neuron.DomvsSub.results.classified.rename = avp.neuron.DomvsSub.results.classified
# colnames(avp.neuron.DomvsSub.results.classified.rename) = paste(colnames(avp.neuron.DomvsSub.results.classified.rename),
#                                                          "DEsingle",
#                                                          sep=".")
# #limma
# avp.neuron.DomvsSub.vector.limma.results.df.rename = avp.neuron.DomvsSub.vector.limma.results.df
# colnames(avp.neuron.DomvsSub.vector.limma.results.df.rename) = paste(colnames(avp.neuron.DomvsSub.vector.limma.results.df.rename),
#                                                               "limma",
#                                                               sep=".")
# 
# #join data
# avp.neuron.DomvsSub.combined.results = full_join(avp.neuron.DomvsSub.vector.edgeR.results.df.rename %>% 
#                                             rownames_to_column('gene'),
#                                           avp.neuron.DomvsSub.results.classified.rename %>% 
#                                             rownames_to_column('gene'),
#                                           by = c('gene')) %>% 
#   full_join(avp.neuron.DomvsSub.vector.limma.results.df.rename %>% 
#               rownames_to_column('gene'))
# 
# 
# 
# #logFC
# png('neuropeptides/avp.oxt/avp/comparison/logFC comparison DEG methods.png')
# chart.Correlation(avp.neuron.DomvsSub.combined.results %>%
#                     mutate(logFC.DEsingle = log10(foldChange.DEsingle)) %>% 
#                     select(c(logFC.edgeR,
#                              logFC.limma,
#                              logFC.DEsingle)),
#                   pch=19)
# dev.off()
# 
# #logFC
# png('neuropeptides/avp.oxt/avp/comparison/pvalue comparison DEG methods.png')
# chart.Correlation(avp.neuron.DomvsSub.combined.results %>%
#                     select(c(pvalue.DEsingle,
#                              PValue.edgeR,
#                              P.Value.limma)),
#                   pch=19)
# dev.off()
# 

# #### AVP neurons reduce ####
# ### subset data to relevant cells
# burtoni.snseq.combined.sct.reduce.avp.neurons.recluster = burtoni.snseq.combined.sct.all.avp.neurons.recluster
# 
# #set idents
# Idents(object = burtoni.snseq.combined.sct.reduce.avp.neurons.recluster) <- "integrated_snn_res.0.8"
# 
# #subset
# burtoni.snseq.combined.sct.reduce.avp.neurons.recluster = subset(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster, 
#                                                                  idents = c("0",
#                                                                             "1",
#                                                                             "2",
#                                                                             "3",
#                                                                             "4"))
# 
# 
# ## run PCA, UMAP, and cluster 
# #use 0.2 resolution
# burtoni.snseq.combined.sct.reduce.avp.neurons.recluster = burtoni.snseq.combined.sct.reduce.avp.neurons.recluster %>% 
#   RunPCA() %>%
#   FindNeighbors(dims = 1:15) %>%
#   RunUMAP(dims = 1:15) %>%
#   FindClusters(resolution = 0.2)
# 
# # ### clustree
# # # cluster across resolutions
# # burtoni.snseq.combined.sct.reduce.avp.neurons.clustree <- Seurat::FindClusters(object = burtoni.snseq.combined.sct.reduce.avp.neurons.recluster, 
# #                                                                             resolution = resolution.range.reduced)
# # #check data
# # #head(burtoni.snseq.combined.sct.all.avp.neurons.clustree[[]])
# # 
# # #clustree
# # clustree(burtoni.snseq.combined.sct.reduce.avp.neurons.clustree, 
# #          prefix = "integrated_snn_res.",
# #          node_colour = 'cluster',
# #          node_size_range = c(10,20),
# #          scale_node_text = TRUE) +
# #   scale_edge_color_continuous(low = "black", 
# #                               high = "black") +
# #   theme(legend.position = "bottom")
# # ggsave('neuropeptides/avp.oxt/avp.reduce/avp.neurons.clustree.png',
# #        width = 10,
# #        height = 10)
# 
# # 
# # ### graph 
# # ##dom vs sub
# # DimPlot(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster, 
# #         reduction = "umap", 
# #         group.by = "orig.ident")
# # ggsave('neuropeptides/avp.oxt/avp.reduce/DomVsSub.dimplot.avp.neurons.all.png')
# # 
# # ## clusters
# # DimPlot(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster,
# #         reduction = "umap",
# #         label = TRUE,
# #         repel = TRUE)
# # ggsave('neuropeptides/avp.oxt/avp.reduce/Clusters.dimplot.avp.neurons.all.png',
# #        width = 10,
# #        height = 10)
# # 
# # # PCA heatmap
# # DimHeatmap(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster, 
# #            dims = 1:6, 
# #            cells = 500, 
# #            balanced = TRUE)
# 
# ##expression level
# burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.expression = full_join(full_join(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster@reductions$umap@cell.embeddings %>% 
#                                                                                         as.data.frame() %>% 
#                                                                                         rownames_to_column("Cell.id"),
#                                                                                       burtoni.snseq.combined.sct.reduce.avp.neurons.recluster@meta.data %>% 
#                                                                                         rownames_to_column("Cell.id")),
#                                                                             burtoni.snseq.combined.sct.reduce.avp.neurons.recluster@assays$SCT@data %>% 
#                                                                               as.data.frame() %>% 
#                                                                               filter(rownames(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster@assays$SCT@data) %in% c('avp',
#                                                                                                                                                                                         'oxt')) %>% 
#                                                                               t() %>% as.data.frame() %>% 
#                                                                               rownames_to_column('Cell.id')) %>% 
#   full_join(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster@reductions$pca@cell.embeddings %>% 
#               as.data.frame() %>% 
#               select(PC_1,
#                      PC_2,
#                      PC_3,
#                      PC_4,
#                      PC_5) %>% 
#               rownames_to_column("Cell.id")) %>% 
#   mutate(Oxt.cell = ifelse(oxt >= 1 & avp >= 1,
#                            'both',
#                            ifelse(oxt >= 1,
#                                   'oxt',
#                                   'avp')))
# 
# # # seperation by cell type
# # table(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.expression %>% 
# #         filter(Oxt.cell != 'both') %>% 
# #         select(Cell.type,
# #                integrated_snn_res.0.2))
# # 
# # # difference social status
# # table(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.expression %>% 
# #         filter(Oxt.cell != 'both') %>% 
# #         select(integrated_snn_res.0.2,
# #                orig.ident))
# # 
# # # difference social status
# # # percent
# # # 313 dom to 194 sub
# # table(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.expression %>% 
# #         filter(Oxt.cell != 'both') %>% 
# #         select(integrated_snn_res.0.2,
# #                orig.ident)) %>% 
# #   as.data.frame() %>% 
# #   pivot_wider(names_from = orig.ident,
# #               values_from = Freq) %>% 
# #   mutate(DomSubPerct = 100 * dom_burtoni_snseq/(sub_burtoni_snseq + dom_burtoni_snseq),
# #          TotalCount = sub_burtoni_snseq + dom_burtoni_snseq) %>% 
# #   ggplot(aes(x= DomSubPerct,
# #              y = TotalCount,
# #              color = integrated_snn_res.0.2)) +
# #   geom_vline(xintercept = 61.9) +
# #   geom_point(size = 5) +
# #   theme_classic() + 
# #   theme(text = element_text(size = 20)) +
# #   xlim(c(49,75))+
# #   ylim(c(0,125))+
# #   ggtitle('Clusters.DomvsSub.ratio.avp.neurons.reduce')
# # ggsave('neuropeptides/avp.oxt/avp.reduce/Clusters.DomvsSub.ratio.avp.neurons.reduce.png',
# #        width = 10,
# #        height = 10)
# # 
# # ## Cell type by cluster
# # table(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.expression %>% 
# #         select(integrated_snn_res.0.2,
# #                Cell.type)) %>% 
# #   as.data.frame()  %>% 
# #   ggplot(aes(x= integrated_snn_res.0.2,
# #              y = Cell.type)) +
# #   geom_tile(aes(fill = Freq)) +
# #   geom_label(aes(label = Freq)) +
# #   scale_fill_gradientn(colours=c("blue","red"))
# # ggsave('neuropeptides/avp.oxt/avp.reduce/Cell.type vs cluster.png',
# #        width = 10,
# #        height = 10)
# 
# ## compare cluster
# burtoni.snseq.combined.sct.reduce.avp.neurons.recluster@meta.data %>% 
#   select(c(Genotype.id,
#            integrated_snn_res.0.2,
#            orig.ident)) %>% 
#   table() %>% 
#   as.data.frame() %>% 
#   filter(Freq > 0) %>% 
#   ggplot(aes(x = integrated_snn_res.0.2,
#              y = Freq,
#              color = orig.ident,
#              group = Genotype.id)) +
#   geom_jitter(height = 0,
#               width = 0.1) +
#   geom_line()+
#   theme_classic()
# ggsave('neuropeptides/avp.oxt/avp.reduce/Genotype.id across clusters avp neuron reduce.png',
#        width = 10,
#        height = 10)

# #### AVP neurons PCA ####
# ## graph PCA
# ##clusters
# #PC 1 vs PC 2
# burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.expression %>%
#   ggplot(aes(x = PC_1,
#              y = PC_2,
#              color = integrated_snn_res.0.2)) +
#   geom_point() +
#   theme_classic()
# ggsave('neuropeptides/avp.oxt/avp.reduce/PCA/PC1.PC2.png',
#        width = 10,
#        height = 10)
# 
# #PC 1 vs PC 3
# burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.expression %>%
#   ggplot(aes(x = PC_1,
#              y = PC_3,
#              color = integrated_snn_res.0.2)) +
#   geom_point() +
#   theme_classic()
# ggsave('neuropeptides/avp.oxt/avp.reduce/PCA/PC1.PC3.png',
#        width = 10,
#        height = 10)
# 
# #PC 2 vs PC 3
# burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.expression %>%
#   ggplot(aes(x = PC_2,
#              y = PC_3,
#              color = integrated_snn_res.0.2)) +
#   geom_point() +
#   theme_classic()
# ggsave('neuropeptides/avp.oxt/avp.reduce/PCA/PC2.PC3.png',
#        width = 10,
#        height = 10)
# 
# ## status
# #PC 1 vs PC 2
# burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.expression %>%
#   ggplot(aes(x = PC_1,
#              y = PC_2,
#              color = orig.ident)) +
#   geom_point() +
#   theme_classic()
# ggsave('neuropeptides/avp.oxt/avp.reduce/PCA/PC1.PC2.status.png',
#        width = 10,
#        height = 10)
# 
# #PC 1 vs PC 3
# burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.expression %>%
#   ggplot(aes(x = PC_1,
#              y = PC_3,
#              color = orig.ident)) +
#   geom_point() +
#   theme_classic()
# ggsave('neuropeptides/avp.oxt/avp.reduce/PCA/PC1.PC3.status.png',
#        width = 10,
#        height = 10)
# 
# #PC 2 vs PC 3
# burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.expression %>%
#   ggplot(aes(x = PC_2,
#              y = PC_3,
#              color = orig.ident)) +
#   geom_point() +
#   theme_classic()
# ggsave('neuropeptides/avp.oxt/avp.reduce/PCA/PC2.PC3.status.png',
#        width = 10,
#        height = 10)
# 
# 
# 
# 
# 
# 
# 

#### limmatrend clusters and status approach ####
### calculate variable genes
## identify top 500 variable genes
burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.variable <- FindVariableFeatures(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster, 
                                                                                      selection.method = "vst", 
                                                                                      nfeatures = 500, 
                                                                                      verbose = F)
# identify top 500 variable genes
avp.neuron.reduce.group.topgenes <- head(VariableFeatures(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.variable), 
                                  500)

## create vector of factor
avp.neuron.reduce.DomvsSub.vector.list = burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.expression %>% 
  mutate(avp.neuron.orig.ident = orig.ident %>% 
           as.factor()) %>% 
  pull(avp.neuron.orig.ident) %>% 
  droplevels()

### counts matrix
## raw read count matrix
## rows = genes, columns = cells
# only keep 500 variable genes
# set negative values to 0
avp.neuron.reduce.DomvsSub.vector.count = GetAssayData(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster) %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column('gene') %>% 
  select(c(gene,
           burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.expression %>% 
             pull(Cell.id))) %>% 
  filter(gene %in% avp.neuron.reduce.group.topgenes) %>% 
  column_to_rownames('gene') %>% 
  as.matrix() %>% 
  pmax(0)


## run with reduce avp neurons
#need list with L with  count and condt
# add in cell cluster with condt.2
## create vector of factor
avp.neuron.reduce.DomvsSub.vector.list.cluster = burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.expression %>% 
  mutate(avp.neuron.cluster = integrated_snn_res.0.2 %>% 
           as.factor()) %>% 
  pull(avp.neuron.cluster) %>% 
  droplevels

#create list
avp.neuron.reduce.vector.limma = list(count = avp.neuron.reduce.DomvsSub.vector.count,
                                   condt = avp.neuron.reduce.DomvsSub.vector.list,
                                   condt.2 = avp.neuron.reduce.DomvsSub.vector.list.cluster)
# add genotype.ID?



#load libraries
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(edgeR))

#create function
run_limmatrend_complex <- function(L) {
  message("limmatrend")
  session_info <- sessionInfo()
  timing <- system.time({
    treat <- paste(L$condt, 
                   L$condt.2,
                   sep = '.')
    design <- model.matrix(~0+treat) # genotype.id?
    contrasts <- makeContrasts(DvsS_0 = treatdom_burtoni_snseq.0 - treatsub_burtoni_snseq.0, 
                               DvsS_1 = treatdom_burtoni_snseq.1 - treatsub_burtoni_snseq.1,
                               DvsS_2 = treatdom_burtoni_snseq.2 - treatsub_burtoni_snseq.2,
                               DvsS_3 = treatdom_burtoni_snseq.3 - treatsub_burtoni_snseq.3, 
                               DvsS_4 = treatdom_burtoni_snseq.4 - treatsub_burtoni_snseq.4,
                               levels = design)
    dge <- DGEList(L$count, 
                   group = treat)
    dge <- calcNormFactors(dge)
    
    y <- new("EList")
    y$E <- edgeR::cpm(dge, 
                      log = TRUE, 
                      prior.count = 3)
    fit <- lmFit(y, 
                 design = design)
    fit <- contrasts.fit(fit , 
                         contrasts)
    fit <- eBayes(fit, 
                  trend = TRUE,
                  robust = TRUE)
    tt0 <- topTable(fit, 
                   n = Inf,
                   coef = "DvsS_0",
                   adjust.method = "BH")
    tt1 <- topTable(fit, 
                    n = Inf,
                    coef = "DvsS_1",
                    adjust.method = "BH")
    tt2 <- topTable(fit, 
                    n = Inf,
                    coef = "DvsS_2",
                    adjust.method = "BH")
    tt3 <- topTable(fit, 
                    n = Inf,
                    coef = "DvsS_3",
                    adjust.method = "BH")
    tt4 <- topTable(fit, 
                    n = Inf,
                    coef = "DvsS_4",
                    adjust.method = "BH")
  })
  
  hist(tt0$P.Value, 50)
  hist(tt0$adj.P.Val, 50)
  hist(tt1$P.Value, 50)
  hist(tt1$adj.P.Val, 50)
  hist(tt2$P.Value, 50)
  hist(tt2$adj.P.Val, 50)
  hist(tt3$P.Value, 50)
  hist(tt3$adj.P.Val, 50)
  hist(tt4$P.Value, 50)
  hist(tt4$adj.P.Val, 50)
  limma::plotMDS(dge, 
                 col = as.numeric(as.factor(L$condt)), 
                 pch = 19)
  plotMD(fit)
  
  list(session_info = session_info,
       timing = timing,
       tt0 = tt0,
       tt1 = tt1,
       tt2 = tt2,
       tt3 = tt3,
       tt4 = tt4)
}


###run function  
avp.neuron.reduce.vector.limma.results = run_limmatrend_complex(avp.neuron.reduce.vector.limma)

# save results to dataframe
avp.neuron.reduce.vector.limma.results.df = full_join(avp.neuron.reduce.vector.limma.results$tt0 %>% 
                                                        rename_with(~paste0(.,"_DvsS_0")) %>% 
                                                        rownames_to_column("Gene"),
                                                      avp.neuron.reduce.vector.limma.results$tt1 %>% 
                                                        rename_with(~paste0(.,"_DvsS_1")) %>% 
                                                        rownames_to_column("Gene")) %>% 
  full_join(avp.neuron.reduce.vector.limma.results$tt2 %>% 
              rename_with(~paste0(.,"_DvsS_2")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(avp.neuron.reduce.vector.limma.results$tt3 %>% 
              rename_with(~paste0(.,"_DvsS_3")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(avp.neuron.reduce.vector.limma.results$tt4 %>% 
              rename_with(~paste0(.,"_DvsS_4")) %>% 
              rownames_to_column("Gene")) 



#add color for significance  
avp.neuron.reduce.vector.limma.results.df = avp.neuron.reduce.vector.limma.results.df %>% 
  mutate(Sig_DvsS_0 = ifelse(adj.P.Val_DvsS_0 <= 0.05,
                      'Sig',
                      'Not Sig'),
         Direction.type_DvsS_0 = ifelse(logFC_DvsS_0 > 0,
                                 'up',
                                 'down'),
         Sig.direction_DvsS_0 = ifelse(Sig_DvsS_0 == 'Sig',
                                Direction.type_DvsS_0,
                                Sig_DvsS_0)) %>% 
  mutate(Sig_DvsS_1 = ifelse(adj.P.Val_DvsS_1 <= 0.05,
                             'Sig',
                             'Not Sig'),
         Direction.type_DvsS_1 = ifelse(logFC_DvsS_1 > 0,
                                        'up',
                                        'down'),
         Sig.direction_DvsS_1 = ifelse(Sig_DvsS_1 == 'Sig',
                                       Direction.type_DvsS_1,
                                       Sig_DvsS_1)) %>% 
  mutate(Sig_DvsS_2 = ifelse(adj.P.Val_DvsS_2 <= 0.05,
                             'Sig',
                             'Not Sig'),
         Direction.type_DvsS_2 = ifelse(logFC_DvsS_2 > 0,
                                        'up',
                                        'down'),
         Sig.direction_DvsS_2 = ifelse(Sig_DvsS_2 == 'Sig',
                                       Direction.type_DvsS_2,
                                       Sig_DvsS_2)) %>% 
  mutate(Sig_DvsS_3 = ifelse(adj.P.Val_DvsS_3 <= 0.05,
                             'Sig',
                             'Not Sig'),
         Direction.type_DvsS_3 = ifelse(logFC_DvsS_3 > 0,
                                        'up',
                                        'down'),
         Sig.direction_DvsS_3 = ifelse(Sig_DvsS_3 == 'Sig',
                                       Direction.type_DvsS_3,
                                       Sig_DvsS_3)) %>% 
  mutate(Sig_DvsS_4 = ifelse(adj.P.Val_DvsS_4 <= 0.05,
                             'Sig',
                             'Not Sig'),
         Direction.type_DvsS_4 = ifelse(logFC_DvsS_4 > 0,
                                        'up',
                                        'down'),
         Sig.direction_DvsS_4 = ifelse(Sig_DvsS_4 == 'Sig',
                                       Direction.type_DvsS_4,
                                       Sig_DvsS_4))


##graph volcano plot
# create log fold change and -log pvalue
# cluter 0
avp.neuron.reduce.vector.limma.results.df %>% 
  mutate(sig.label_DvsS_0 = ifelse(Sig_DvsS_0 == 'Sig',
                                   Gene,
                                   '')) %>% 
  ggplot(aes(x = logFC_DvsS_0,
             y = -log10(adj.P.Val_DvsS_0),
             color = Sig.direction_DvsS_0))+
  geom_hline(yintercept = -log10(0.05),
             linetype = 'dotted') +
  geom_point(size = 5) +
  theme_classic() + 
  scale_color_manual(values=c("21B9CA", 
                              "grey", 
                              "orange3")) +
  geom_text(aes(label = sig.label_DvsS_0),
            vjust = 0, 
            nudge_y = 0.10,
            size = 5) +
  theme(text = element_text(size = 20),
        legend.position = 'none')
ggsave('neuropeptides/avp.oxt/avp.reduce/limmatrend/limma.avp.neurons.reduce.volcano.cluster0.png',
       width = 5,
       height = 5)

# cluter 1
avp.neuron.reduce.vector.limma.results.df %>% 
  mutate(sig.label_DvsS_1 = ifelse(Sig_DvsS_1 == 'Sig',
                                   Gene,
                                   '')) %>% 
  ggplot(aes(x = logFC_DvsS_1,
             y = -log10(adj.P.Val_DvsS_1),
             color = Sig.direction_DvsS_1))+
  geom_hline(yintercept = -log10(0.05),
             linetype = 'dotted') +
  geom_point(size = 5) +
  theme_classic() + 
  scale_color_manual(values=c("21B9CA", 
                              "grey", 
                              "orange3")) +
  geom_text(aes(label = sig.label_DvsS_1),
            vjust = 0, 
            nudge_y = 0.10,
            size = 5) +
  theme(text = element_text(size = 20),
        legend.position = 'none')
ggsave('neuropeptides/avp.oxt/avp.reduce/limmatrend/limma.avp.neurons.reduce.volcano.cluster1.png',
       width = 5,
       height = 5)

# cluter 2
avp.neuron.reduce.vector.limma.results.df %>% 
  mutate(sig.label_DvsS_2 = ifelse(Sig_DvsS_2 == 'Sig',
                                   Gene,
                                   '')) %>% 
  ggplot(aes(x = logFC_DvsS_2,
             y = -log10(adj.P.Val_DvsS_2),
             color = Sig.direction_DvsS_2))+
  geom_hline(yintercept = -log10(0.05),
             linetype = 'dotted') +
  geom_point(size = 5) +
  theme_classic() + 
  scale_color_manual(values=c(#"21B9CA", 
                              "grey", 
                              "orange3")) +
  geom_text(aes(label = sig.label_DvsS_2),
            vjust = 0, 
            nudge_y = 0.10,
            size = 5) +
  theme(text = element_text(size = 20),
        legend.position = 'none')
ggsave('neuropeptides/avp.oxt/avp.reduce/limmatrend/limma.avp.neurons.reduce.volcano.cluster2.png',
       width = 5,
       height = 5)

# cluter 3
avp.neuron.reduce.vector.limma.results.df %>% 
  mutate(sig.label_DvsS_3 = ifelse(Sig_DvsS_3 == 'Sig',
                                   Gene,
                                   '')) %>% 
  ggplot(aes(x = logFC_DvsS_3,
             y = -log10(adj.P.Val_DvsS_3),
             color = Sig.direction_DvsS_3)) +
  geom_hline(yintercept = -log10(0.05),
             linetype = 'dotted') +
  geom_point(size = 5) +
  theme_classic() + 
  scale_color_manual(values=c("21B9CA", 
                              "grey", 
                              "orange3")) +
  geom_text(aes(label = sig.label_DvsS_3),
            vjust = 0, 
            nudge_y = 0.10,
            size = 5) +
  theme(text = element_text(size = 20),
        legend.position = 'none')
ggsave('neuropeptides/avp.oxt/avp.reduce/limmatrend/limma.avp.neurons.reduce.volcano.cluster3.png',
       width = 5,
       height = 5)

# cluter 4
  avp.neuron.reduce.vector.limma.results.df %>% 
  mutate(sig.label_DvsS_4 = ifelse(Sig_DvsS_4 == 'Sig',
                            Gene,
                            '')) %>% 
  ggplot(aes(x = logFC_DvsS_4,
             y = -log10(adj.P.Val_DvsS_4),
             color = Sig.direction_DvsS_4))+
    geom_hline(yintercept = -log10(0.05),
               linetype = 'dotted') +
  geom_point(size = 5) +
  theme_classic() + 
  scale_color_manual(values=c("21B9CA", 
                              "grey", 
                              "orange3")) +
  geom_text(aes(label = sig.label_DvsS_4),
            vjust = 0, 
            nudge_y = 0.10,
            size = 5) +
  theme(text = element_text(size = 20),
        legend.position = 'none')
ggsave('neuropeptides/avp.oxt/avp.reduce/limmatrend/limma.avp.neurons.reduce.volcano.cluster4.png',
       width = 5,
       height = 5)

#### limmatrend genotype, clusters combine, and status approach ####
### calculate variable genes
## identify top 500 variable genes
burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.variable <- FindVariableFeatures(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster, 
                                                                                         selection.method = "vst", 
                                                                                         nfeatures = 500, 
                                                                                         verbose = F)
# identify top 500 variable genes
avp.neuron.reduce.group.topgenes <- head(VariableFeatures(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.variable), 
                                         500)

## create vector of factor
avp.neuron.reduce.DomvsSub.vector.list = burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.expression %>% 
  mutate(avp.neuron.orig.ident = orig.ident %>% 
           as.factor()) %>% 
  pull(avp.neuron.orig.ident) %>% 
  droplevels()

### counts matrix
## raw read count matrix
## rows = genes, columns = cells
# only keep 500 variable genes
# set negative values to 0
avp.neuron.reduce.DomvsSub.vector.count = GetAssayData(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster) %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column('gene') %>% 
  select(c(gene,
           burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.expression %>% 
             pull(Cell.id))) %>% 
  filter(gene %in% avp.neuron.reduce.group.topgenes) %>% 
  column_to_rownames('gene') %>% 
  as.matrix() %>% 
  pmax(0)


## run with reduce avp neurons
#need list with L with  count and condt
# add in cell cluster with condt.2
## create vector of factor
avp.neuron.reduce.DomvsSub.vector.list.cluster = burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.expression %>% 
  mutate(avp.neuron.cluster = as.character(integrated_snn_res.0.2),
         avp.neuron.cluster.replace = ifelse(avp.neuron.cluster == '4',
                                             '2',
                                             avp.neuron.cluster)) %>% 
  mutate(avp.neuron.cluster.replace = avp.neuron.cluster.replace %>% 
           as.factor()) %>% 
  pull(avp.neuron.cluster.replace) %>% 
  droplevels()


# add in genotype
## create vector of factor
avp.neuron.reduce.DomvsSub.vector.list.genotype = burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.expression %>% 
  mutate(genotype = Genotype.id %>% 
           as.factor()) %>% 
  pull(genotype) %>% 
  droplevels()

#create list
avp.neuron.reduce.vector.limma = list(count = avp.neuron.reduce.DomvsSub.vector.count,
                                      condt = avp.neuron.reduce.DomvsSub.vector.list,
                                      condt.2 = avp.neuron.reduce.DomvsSub.vector.list.cluster,
                                      genotype = avp.neuron.reduce.DomvsSub.vector.list.genotype)
# add genotype.ID?



#load libraries
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(edgeR))

#create function
run_limmatrend_complex_genotype <- function(L) {
  message("limmatrend")
  session_info <- sessionInfo()
  timing <- system.time({
    treat <- paste(L$condt, 
                   L$condt.2,
                   sep = '.')
    genotype <- L$genotype
    design <- model.matrix(~0+treat+genotype) 
    contrasts <- makeContrasts(DvsS_0 = treatdom_burtoni_snseq.0 - treatsub_burtoni_snseq.0, 
                               DvsS_1 = treatdom_burtoni_snseq.1 - treatsub_burtoni_snseq.1,
                               DvsS_2 = treatdom_burtoni_snseq.2 - treatsub_burtoni_snseq.2,
                               DvsS_3 = treatdom_burtoni_snseq.3 - treatsub_burtoni_snseq.3,
                               DvsS=(treatdom_burtoni_snseq.0+treatdom_burtoni_snseq.1+treatdom_burtoni_snseq.2+treatdom_burtoni_snseq.3)/4-(treatsub_burtoni_snseq.0+treatsub_burtoni_snseq.1+treatsub_burtoni_snseq.2+treatsub_burtoni_snseq.3)/4,
                               levels = design)
    dge <- DGEList(L$count, 
                   group = treat)
    dge <- calcNormFactors(dge)
    
    y <- new("EList")
    y$E <- edgeR::cpm(dge, 
                      log = TRUE, 
                      prior.count = 3)
    fit <- lmFit(y, 
                 design = design)
    fit <- contrasts.fit(fit , 
                         contrasts)
    fit <- eBayes(fit, 
                  trend = TRUE,
                  robust = TRUE)
    tt0 <- topTable(fit, 
                    n = Inf,
                    coef = "DvsS_0",
                    adjust.method = "BH")
    tt1 <- topTable(fit, 
                    n = Inf,
                    coef = "DvsS_1",
                    adjust.method = "BH")
    tt2 <- topTable(fit, 
                    n = Inf,
                    coef = "DvsS_2",
                    adjust.method = "BH")
    tt3 <- topTable(fit, 
                    n = Inf,
                    coef = "DvsS_3",
                    adjust.method = "BH")
    ttDS <- topTable(fit, 
                    n = Inf,
                    coef = "DvsS",
                    adjust.method = "BH")
  })
  # Open pdf file
  pdf(file= "./neuropeptides/avp.oxt/avp/limmatrend/genotype.combine.cluster/limmatrend.histograms.pdf" )
  # create a 2X2 grid
  par( mfrow= c(2,2) )
  #graph
  hist(tt0$P.Value, 50)
  hist(tt0$adj.P.Val, 50)
  hist(tt1$P.Value, 50)
  hist(tt1$adj.P.Val, 50)
  hist(tt2$P.Value, 50)
  hist(tt2$adj.P.Val, 50)
  hist(tt3$P.Value, 50)
  hist(tt3$adj.P.Val, 50)
  hist(ttDS$P.Value, 50)
  hist(ttDS$adj.P.Val, 50)
  dev.off()
  
  # Open pdf file
  pdf(file= "./neuropeptides/avp.oxt/avp/limmatrend/genotype.combine.cluster/limmatrend.MDS.pdf" )
  # create a 2X1 grid
  par( mfrow= c(2,1) )
  limma::plotMDS(dge, 
                 col = as.numeric(as.factor(L$condt)), 
                 pch = 19)
  plotMD(fit)
  limma::plotMDS(dge, 
                 col = as.numeric(as.factor(L$condt.2)), 
                 pch = 19)
  plotMD(fit)
  dev.off()
  
  #print results
  list(session_info = session_info,
       timing = timing,
       tt0 = tt0,
       tt1 = tt1,
       tt2 = tt2,
       tt3 = tt3,
       ttDS = ttDS)
}


###run function  
avp.neuron.reduce.genotype.vector.limma.results = run_limmatrend_complex_genotype(avp.neuron.reduce.vector.limma)

# save results to dataframe
avp.neuron.reduce.genotype.vector.limma.results.df = full_join(avp.neuron.reduce.genotype.vector.limma.results$tt0 %>% 
                                                        rename_with(~paste0(.,"_DvsS_0")) %>% 
                                                        rownames_to_column("Gene"),
                                                        avp.neuron.reduce.genotype.vector.limma.results$tt1 %>% 
                                                        rename_with(~paste0(.,"_DvsS_1")) %>% 
                                                        rownames_to_column("Gene")) %>% 
  full_join(avp.neuron.reduce.genotype.vector.limma.results$tt2 %>% 
              rename_with(~paste0(.,"_DvsS_2")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(avp.neuron.reduce.genotype.vector.limma.results$tt3 %>% 
              rename_with(~paste0(.,"_DvsS_3")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(avp.neuron.reduce.genotype.vector.limma.results$ttDS %>% 
              rename_with(~paste0(.,"_DvsS")) %>% 
              rownames_to_column("Gene")) 



#add color for significance  
avp.neuron.reduce.genotype.vector.limma.results.df = avp.neuron.reduce.genotype.vector.limma.results.df %>% 
  mutate(Sig_DvsS_0 = ifelse(adj.P.Val_DvsS_0 <= 0.05,
                             'Sig',
                             'Not Sig'),
         Direction.type_DvsS_0 = ifelse(logFC_DvsS_0 > 0,
                                        'up',
                                        'down'),
         Sig.direction_DvsS_0 = ifelse(Sig_DvsS_0 == 'Sig',
                                       Direction.type_DvsS_0,
                                       Sig_DvsS_0)) %>% 
  mutate(Sig_DvsS_1 = ifelse(adj.P.Val_DvsS_1 <= 0.05,
                             'Sig',
                             'Not Sig'),
         Direction.type_DvsS_1 = ifelse(logFC_DvsS_1 > 0,
                                        'up',
                                        'down'),
         Sig.direction_DvsS_1 = ifelse(Sig_DvsS_1 == 'Sig',
                                       Direction.type_DvsS_1,
                                       Sig_DvsS_1)) %>% 
  mutate(Sig_DvsS_2 = ifelse(adj.P.Val_DvsS_2 <= 0.05,
                             'Sig',
                             'Not Sig'),
         Direction.type_DvsS_2 = ifelse(logFC_DvsS_2 > 0,
                                        'up',
                                        'down'),
         Sig.direction_DvsS_2 = ifelse(Sig_DvsS_2 == 'Sig',
                                       Direction.type_DvsS_2,
                                       Sig_DvsS_2)) %>% 
  mutate(Sig_DvsS_3 = ifelse(adj.P.Val_DvsS_3 <= 0.05,
                             'Sig',
                             'Not Sig'),
         Direction.type_DvsS_3 = ifelse(logFC_DvsS_3 > 0,
                                        'up',
                                        'down'),
         Sig.direction_DvsS_3 = ifelse(Sig_DvsS_3 == 'Sig',
                                       Direction.type_DvsS_3,
                                       Sig_DvsS_3))  %>% 
  mutate(Sig_DvsS = ifelse(adj.P.Val_DvsS <= 0.05,
                             'Sig',
                             'Not Sig'),
         Direction.type_DvsS = ifelse(logFC_DvsS > 0,
                                        'up',
                                        'down'),
         Sig.direction_DvsS = ifelse(Sig_DvsS == 'Sig',
                                       Direction.type_DvsS,
                                       Sig_DvsS))


##graph volcano plot
# create comparison list
limma.genotpye.vector = c("DvsS_0",
                          "DvsS_1",
                          "DvsS_2",
                          "DvsS_3",
                          "DvsS")

# 
# 
# avp.neuron.reduce.genotype.vector.limma.results.df %>% 
#   mutate(Sig_DvsS = ifelse(adj.P.Val_DvsS <= 0.1,
#                            'Sig',
#                            'Not Sig'),
#          Direction.type_DvsS = ifelse(logFC_DvsS > 0,
#                                       'up',
#                                       'down'),
#          Sig.direction_DvsS = ifelse(Sig_DvsS == 'Sig',
#                                      Direction.type_DvsS,
#                                      Sig_DvsS)) %>% 
#   mutate(sig.label = ifelse(get(paste0("Sig_", i)) == 'Sig',
#                             Gene,
#                             '')) %>% 
#   ggplot(aes(x = get(paste0("logFC_", i)),
#              y = -log10(get(paste0("adj.P.Val_", i))),
#              color = get(paste0("Sig.direction_", i))))+
#   geom_hline(yintercept = -log10(0.1),
#              linetype = 'dotted') +
#   geom_point(size = 5) +
#   theme_classic() + 
#   scale_color_manual(values=c("21B9CA", 
#                               "grey", 
#                               "orange3")) +
#   geom_text(aes(label = sig.label),
#             vjust = 0, 
#             nudge_y = 0.10,
#             size = 5) +
#   theme(text = element_text(size = 20),
#         legend.position = 'none') +
#   xlab(paste0("logFC_", i)) +
#   ylab( paste0("-log10(adj.P.Val_", i,")")) +
#   ggtitle(paste0(i, " volcano plot"))
# ggsave(paste0('neuropeptides/avp.oxt/avp/limmatrend/limma.avp.neurons.reduce.volcano.', i, '.png'),
#        width = 5,
#        height = 5)


for (i in limma.genotpye.vector) {
  # graph volcano plot
  avp.neuron.reduce.genotype.vector.limma.results.df %>% 
    mutate(sig.label = ifelse(get(paste0("Sig_", i)) == 'Sig',
                                     Gene,
                                     '')) %>% 
    ggplot(aes(x = get(paste0("logFC_", i)),
               y = -log10(get(paste0("adj.P.Val_", i))),
               color = get(paste0("Sig.direction_", i))))+
    geom_hline(yintercept = -log10(0.05),
               linetype = 'dotted') +
    geom_point(size = 5) +
    theme_classic() + 
    scale_color_manual(values=c("21B9CA", 
                                "grey", 
                                "orange3")) +
    geom_text(aes(label = sig.label),
              vjust = 0, 
              nudge_y = 0.10,
              size = 5) +
    theme(text = element_text(size = 20),
          legend.position = 'none') +
    xlab(paste0("logFC_", i)) +
    ylab( paste0("-log10(adj.P.Val_", i,")")) +
    ggtitle(paste0(i, " volcano plot"))
  ggsave(paste0('neuropeptides/avp.oxt/avp/limmatrend/genotype.combine.cluster/limma.avp.neurons.reduce.volcano.', i, '.png'),
         width = 5,
         height = 5)
}

## get list of significant genes
DEGs.cluster.genotype = avp.neuron.reduce.genotype.vector.limma.results.df %>% 
  mutate(Keep = 0,
         Keep = case_when(Sig.direction_DvsS_0 != "Not Sig" ~ Keep + 1,
                          Sig.direction_DvsS_1 != "Not Sig" ~ Keep + 1,
                          Sig.direction_DvsS_2 != "Not Sig" ~ Keep + 1,
                          Sig.direction_DvsS_3 != "Not Sig" ~ Keep + 1,
                          Sig.direction_DvsS_4 != "Not Sig" ~ Keep + 1,
                          TRUE ~ Keep)) %>% 
  filter(Keep != 0) %>%  
  select(Gene,
         Sig.direction_DvsS_0,
         Sig.direction_DvsS_1,
         Sig.direction_DvsS_2,
         Sig.direction_DvsS_3,
         Sig.direction_DvsS_4,
         adj.P.Val_DvsS_0,
         adj.P.Val_DvsS_1,
         adj.P.Val_DvsS_2,
         adj.P.Val_DvsS_3,
         adj.P.Val_DvsS_4) %>% 
  mutate(Keep = case_when(Sig.direction_DvsS_0 == "up" | Sig.direction_DvsS_1 == "up" | Sig.direction_DvsS_2 == "up" | Sig.direction_DvsS_3 == "up" | Sig.direction_DvsS_4 == "up" ~ "Dom",
                          TRUE ~ "Sub")) %>% 
  select(-c(Sig.direction_DvsS_0,
            Sig.direction_DvsS_1,
            Sig.direction_DvsS_2,
            Sig.direction_DvsS_3,
            Sig.direction_DvsS_4)) %>% 
  group_by(Gene) %>% 
  mutate(adj.P.Val = min(c(adj.P.Val_DvsS_0,
                           adj.P.Val_DvsS_1,
                           adj.P.Val_DvsS_2,
                           adj.P.Val_DvsS_3,
                           adj.P.Val_DvsS_4))) %>% 
  select(-c(adj.P.Val_DvsS_0,
            adj.P.Val_DvsS_1,
            adj.P.Val_DvsS_2,
            adj.P.Val_DvsS_3,
            adj.P.Val_DvsS_4)) 
write_csv(DEGs.cluster.genotype,
          'neuropeptides/avp.oxt/avp/limmatrend/genotype/DEGs.cluster.csv')

# save hub genes for GO analysis with:
# to get GO term accession use biomart: http://useast.ensembl.org/biomart/martview/

## combine go terms
DEGs.cluster.genotype.GO.terms = read.csv('neuropeptides/avp.oxt/avp/limmatrend/genotype/DEGs.cluster.GO.terms.csv')
  
# join with pvalue 
DEGs.cluster.genotype.GO.terms.pvalue = DEGs.cluster.genotype.GO.terms %>% 
  left_join(DEGs.cluster.genotype ,
            by = c("Gene.name" = "Gene")) %>% 
  filter(!is.na(adj.P.Val)) %>% 
  rbind(DEGs.cluster.genotype.GO.terms %>% 
          left_join(DEGs.cluster.genotype ,
                    by = c("Gene.stable.ID" = "Gene")) %>% 
          filter(!is.na(adj.P.Val)))  
  
# save
DEGs.cluster.genotype.GO.terms.pvalue %>% 
  select(c(Gene.stable.ID,
           GO.term.accession,
           Gene.name,
           Keep,
           adj.P.Val)) %>% 
  distinct() %>% 
  relocate(GO.term.accession, 
           .after = last_col())%>% 
  relocate(adj.P.Val, 
           .after = last_col()) %>% 
  arrange(Keep) %>% 
write_csv('neuropeptides/avp.oxt/avp/limmatrend/genotype/DEGs.cluster.GO.terms.pvalue.csv')
  
# to look for enrichment and create figures use GO term accession and kME (higher is better) for each module: http://revigo.irb.hr/
# remember to sort by kME and filter out blank GO terms before copying to revigo
# all module genes




# 




#### limmatrend genotype, clusters, and status approach ####
### calculate variable genes
## identify top 500 variable genes
burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.variable <- FindVariableFeatures(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster, 
                                                                                         selection.method = "vst", 
                                                                                         nfeatures = 500, 
                                                                                         verbose = F)
# identify top 500 variable genes
avp.neuron.reduce.group.topgenes <- head(VariableFeatures(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.variable), 
                                         500)

## create vector of factor
avp.neuron.reduce.DomvsSub.vector.list = burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.expression %>% 
  mutate(avp.neuron.orig.ident = orig.ident %>% 
           as.factor()) %>% 
  pull(avp.neuron.orig.ident) %>% 
  droplevels()

### counts matrix
## raw read count matrix
## rows = genes, columns = cells
# only keep 500 variable genes
# set negative values to 0
avp.neuron.reduce.DomvsSub.vector.count = GetAssayData(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster) %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column('gene') %>% 
  select(c(gene,
           burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.expression %>% 
             pull(Cell.id))) %>% 
  filter(gene %in% avp.neuron.reduce.group.topgenes) %>% 
  column_to_rownames('gene') %>% 
  as.matrix() %>% 
  pmax(0)


## run with reduce avp neurons
#need list with L with  count and condt
# add in cell cluster with condt.2
## create vector of factor
avp.neuron.reduce.DomvsSub.vector.list.cluster = burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.expression %>% 
  mutate(avp.neuron.cluster = integrated_snn_res.0.2 %>% 
           as.factor()) %>% 
  pull(avp.neuron.cluster) %>% 
  droplevels


# add in genotype
## create vector of factor
avp.neuron.reduce.DomvsSub.vector.list.genotype = burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.expression %>% 
  mutate(genotype = Genotype.id %>% 
           as.factor()) %>% 
  pull(genotype) %>% 
  droplevels()

#create list
avp.neuron.reduce.vector.limma = list(count = avp.neuron.reduce.DomvsSub.vector.count,
                                      condt = avp.neuron.reduce.DomvsSub.vector.list,
                                      condt.2 = avp.neuron.reduce.DomvsSub.vector.list.cluster,
                                      genotype = avp.neuron.reduce.DomvsSub.vector.list.genotype)
# add genotype.ID?



#load libraries
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(edgeR))

#create function
run_limmatrend_complex_genotype <- function(L) {
  message("limmatrend")
  session_info <- sessionInfo()
  timing <- system.time({
    treat <- paste(L$condt, 
                   L$condt.2,
                   sep = '.')
    genotype <- L$genotype
    design <- model.matrix(~0+treat+genotype) 
    contrasts <- makeContrasts(DvsS_0 = treatdom_burtoni_snseq.0 - treatsub_burtoni_snseq.0, 
                               DvsS_1 = treatdom_burtoni_snseq.1 - treatsub_burtoni_snseq.1,
                               DvsS_2 = treatdom_burtoni_snseq.2 - treatsub_burtoni_snseq.2,
                               DvsS_3 = treatdom_burtoni_snseq.3 - treatsub_burtoni_snseq.3, 
                               DvsS_4 = treatdom_burtoni_snseq.4 - treatsub_burtoni_snseq.4,
                               DvsS=(treatdom_burtoni_snseq.0+treatdom_burtoni_snseq.1+treatdom_burtoni_snseq.2+treatdom_burtoni_snseq.3+treatdom_burtoni_snseq.4)/5-(treatsub_burtoni_snseq.0+treatsub_burtoni_snseq.1+treatsub_burtoni_snseq.2+treatsub_burtoni_snseq.3+treatsub_burtoni_snseq.4)/5,
                               levels = design)
    dge <- DGEList(L$count, 
                   group = treat)
    dge <- calcNormFactors(dge)
    
    y <- new("EList")
    y$E <- edgeR::cpm(dge, 
                      log = TRUE, 
                      prior.count = 3)
    fit <- lmFit(y, 
                 design = design)
    fit <- contrasts.fit(fit , 
                         contrasts)
    fit <- eBayes(fit, 
                  trend = TRUE,
                  robust = TRUE)
    tt0 <- topTable(fit, 
                    n = Inf,
                    coef = "DvsS_0",
                    adjust.method = "BH")
    tt1 <- topTable(fit, 
                    n = Inf,
                    coef = "DvsS_1",
                    adjust.method = "BH")
    tt2 <- topTable(fit, 
                    n = Inf,
                    coef = "DvsS_2",
                    adjust.method = "BH")
    tt3 <- topTable(fit, 
                    n = Inf,
                    coef = "DvsS_3",
                    adjust.method = "BH")
    tt4 <- topTable(fit, 
                    n = Inf,
                    coef = "DvsS_4",
                    adjust.method = "BH")
    ttDS <- topTable(fit, 
                     n = Inf,
                     coef = "DvsS",
                     adjust.method = "BH")
  })
  # Open pdf file
  pdf(file= "./neuropeptides/avp.oxt/avp/limmatrend/genotype/limmatrend.histograms.pdf" )
  # create a 2X2 grid
  par( mfrow= c(2,2) )
  #graph
  hist(tt0$P.Value, 50)
  hist(tt0$adj.P.Val, 50)
  hist(tt1$P.Value, 50)
  hist(tt1$adj.P.Val, 50)
  hist(tt2$P.Value, 50)
  hist(tt2$adj.P.Val, 50)
  hist(tt3$P.Value, 50)
  hist(tt3$adj.P.Val, 50)
  hist(tt4$P.Value, 50)
  hist(tt4$adj.P.Val, 50)
  hist(ttDS$P.Value, 50)
  hist(ttDS$adj.P.Val, 50)
  dev.off()
  
  # Open pdf file
  pdf(file= "./neuropeptides/avp.oxt/avp/limmatrend/genotype/limmatrend.MDS.pdf" )
  # create a 2X1 grid
  par( mfrow= c(2,1) )
  limma::plotMDS(dge, 
                 col = as.numeric(as.factor(L$condt)), 
                 pch = 19)
  plotMD(fit)
  limma::plotMDS(dge, 
                 col = as.numeric(as.factor(L$condt.2)), 
                 pch = 19)
  plotMD(fit)
  dev.off()
  
  #print results
  list(session_info = session_info,
       timing = timing,
       tt0 = tt0,
       tt1 = tt1,
       tt2 = tt2,
       tt3 = tt3,
       tt4 = tt4,
       ttDS = ttDS)
}


###run function  
avp.neuron.reduce.genotype.vector.limma.results = run_limmatrend_complex_genotype(avp.neuron.reduce.vector.limma)

# save results to dataframe
avp.neuron.reduce.genotype.vector.limma.results.df = full_join(avp.neuron.reduce.genotype.vector.limma.results$tt0 %>% 
                                                                 rename_with(~paste0(.,"_DvsS_0")) %>% 
                                                                 rownames_to_column("Gene"),
                                                               avp.neuron.reduce.genotype.vector.limma.results$tt1 %>% 
                                                                 rename_with(~paste0(.,"_DvsS_1")) %>% 
                                                                 rownames_to_column("Gene")) %>% 
  full_join(avp.neuron.reduce.genotype.vector.limma.results$tt2 %>% 
              rename_with(~paste0(.,"_DvsS_2")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(avp.neuron.reduce.genotype.vector.limma.results$tt3 %>% 
              rename_with(~paste0(.,"_DvsS_3")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(avp.neuron.reduce.genotype.vector.limma.results$tt4 %>% 
              rename_with(~paste0(.,"_DvsS_4")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(avp.neuron.reduce.genotype.vector.limma.results$ttDS %>% 
              rename_with(~paste0(.,"_DvsS")) %>% 
              rownames_to_column("Gene")) 



#add color for significance  
avp.neuron.reduce.genotype.vector.limma.results.df = avp.neuron.reduce.genotype.vector.limma.results.df %>% 
  mutate(Sig_DvsS_0 = ifelse(adj.P.Val_DvsS_0 <= 0.05,
                             'Sig',
                             'Not Sig'),
         Direction.type_DvsS_0 = ifelse(logFC_DvsS_0 > 0,
                                        'up',
                                        'down'),
         Sig.direction_DvsS_0 = ifelse(Sig_DvsS_0 == 'Sig',
                                       Direction.type_DvsS_0,
                                       Sig_DvsS_0)) %>% 
  mutate(Sig_DvsS_1 = ifelse(adj.P.Val_DvsS_1 <= 0.05,
                             'Sig',
                             'Not Sig'),
         Direction.type_DvsS_1 = ifelse(logFC_DvsS_1 > 0,
                                        'up',
                                        'down'),
         Sig.direction_DvsS_1 = ifelse(Sig_DvsS_1 == 'Sig',
                                       Direction.type_DvsS_1,
                                       Sig_DvsS_1)) %>% 
  mutate(Sig_DvsS_2 = ifelse(adj.P.Val_DvsS_2 <= 0.05,
                             'Sig',
                             'Not Sig'),
         Direction.type_DvsS_2 = ifelse(logFC_DvsS_2 > 0,
                                        'up',
                                        'down'),
         Sig.direction_DvsS_2 = ifelse(Sig_DvsS_2 == 'Sig',
                                       Direction.type_DvsS_2,
                                       Sig_DvsS_2)) %>% 
  mutate(Sig_DvsS_3 = ifelse(adj.P.Val_DvsS_3 <= 0.05,
                             'Sig',
                             'Not Sig'),
         Direction.type_DvsS_3 = ifelse(logFC_DvsS_3 > 0,
                                        'up',
                                        'down'),
         Sig.direction_DvsS_3 = ifelse(Sig_DvsS_3 == 'Sig',
                                       Direction.type_DvsS_3,
                                       Sig_DvsS_3)) %>% 
  mutate(Sig_DvsS_4 = ifelse(adj.P.Val_DvsS_4 <= 0.05,
                             'Sig',
                             'Not Sig'),
         Direction.type_DvsS_4 = ifelse(logFC_DvsS_4 > 0,
                                        'up',
                                        'down'),
         Sig.direction_DvsS_4 = ifelse(Sig_DvsS_4 == 'Sig',
                                       Direction.type_DvsS_4,
                                       Sig_DvsS_4)) %>% 
  mutate(Sig_DvsS = ifelse(adj.P.Val_DvsS <= 0.05,
                           'Sig',
                           'Not Sig'),
         Direction.type_DvsS = ifelse(logFC_DvsS > 0,
                                      'up',
                                      'down'),
         Sig.direction_DvsS = ifelse(Sig_DvsS == 'Sig',
                                     Direction.type_DvsS,
                                     Sig_DvsS))


##graph volcano plot
# create comparison list
limma.genotpye.vector = c("DvsS_0",
                          "DvsS_1",
                          "DvsS_2",
                          "DvsS_3",
                          "DvsS_4",
                          "DvsS")

for (i in limma.genotpye.vector) {
  # graph volcano plot
  avp.neuron.reduce.genotype.vector.limma.results.df %>% 
    mutate(sig.label = ifelse(get(paste0("Sig_", i)) == 'Sig',
                              Gene,
                              '')) %>% 
    ggplot(aes(x = get(paste0("logFC_", i)),
               y = -log10(get(paste0("adj.P.Val_", i))),
               color = get(paste0("Sig.direction_", i))))+
    geom_hline(yintercept = -log10(0.05),
               linetype = 'dotted') +
    geom_point(size = 5) +
    theme_classic() + 
    scale_color_manual(values=c("21B9CA", 
                                "grey", 
                                "orange3")) +
    geom_text(aes(label = sig.label),
              vjust = 0, 
              nudge_y = 0.10,
              size = 5) +
    theme(text = element_text(size = 20),
          legend.position = 'none') +
    xlab(paste0("logFC_", i)) +
    ylab( paste0("-log10(adj.P.Val_", i,")")) +
    ggtitle(paste0(i, " volcano plot"))
  ggsave(paste0('neuropeptides/avp.oxt/avp/limmatrend/genotype/limma.avp.neurons.reduce.volcano.', i, '.png'),
         width = 5,
         height = 5)
}

## get list of significant genes
DEGs.cluster.genotype = avp.neuron.reduce.genotype.vector.limma.results.df %>% View()
  mutate(Keep = 0,
         Keep.2 = case_when(Sig.direction_DvsS_0 != "Not Sig" ~ Keep + 1,
                          Sig.direction_DvsS_1 != "Not Sig" ~ Keep + 1,
                          Sig.direction_DvsS_2 != "Not Sig" ~ Keep + 1,
                          Sig.direction_DvsS_3 != "Not Sig" ~ Keep + 1,
                          Sig.direction_DvsS_4 != "Not Sig" ~ Keep + 1,
                          Sig.direction_DvsS != "Not Sig" ~ Keep + 1,
                          TRUE ~ Keep)) %>% 
  filter(Keep.2 != 0) %>% 
  select(Gene,
         Sig.direction_DvsS_0,
         Sig.direction_DvsS_1,
         Sig.direction_DvsS_2,
         Sig.direction_DvsS_3,
         Sig.direction_DvsS_4,
         Sig.direction_DvsS,
         adj.P.Val_DvsS_0,
         adj.P.Val_DvsS_1,
         adj.P.Val_DvsS_2,
         adj.P.Val_DvsS_3,
         adj.P.Val_DvsS_4,
         adj.P.Val_DvsS) %>% 
  mutate(Keep = case_when(Sig.direction_DvsS_0 == "up" | Sig.direction_DvsS_1 == "up" | Sig.direction_DvsS_2 == "up" | Sig.direction_DvsS_3 == "up" | Sig.direction_DvsS_4 == "up" ~ "Dom"| Sig.direction_DvsS == "up" ~ "Dom",
                          TRUE ~ "Sub")) %>% 
  group_by(Gene) %>% 
  mutate(adj.P.Val = min(c(adj.P.Val_DvsS_0,
                           adj.P.Val_DvsS_1,
                           adj.P.Val_DvsS_2,
                           adj.P.Val_DvsS_3,
                           adj.P.Val_DvsS_4,
                           adj.P.Val_DvsS))) %>% 
  select(-c(adj.P.Val_DvsS_0,
            adj.P.Val_DvsS_1,
            adj.P.Val_DvsS_2,
            adj.P.Val_DvsS_3,
            adj.P.Val_DvsS_4,
            adj.P.Val_DvsS)) %>% 
write_csv(DEGs.cluster.genotype,
          'neuropeptides/avp.oxt/avp/limmatrend/genotype/DEGs.cluster.csv')

# save hub genes for GO analysis with:
# to get GO term accession use biomart: http://useast.ensembl.org/biomart/martview/

## combine go terms
DEGs.cluster.genotype.GO.terms = read.csv('neuropeptides/avp.oxt/avp/limmatrend/genotype/DEGs.cluster.GO.terms.csv')

# join with pvalue 
DEGs.cluster.genotype.GO.terms.pvalue = DEGs.cluster.genotype.GO.terms %>% 
  left_join(DEGs.cluster.genotype ,
            by = c("Gene.name" = "Gene")) %>% 
  filter(!is.na(adj.P.Val)) %>% 
  rbind(DEGs.cluster.genotype.GO.terms %>% 
          left_join(DEGs.cluster.genotype ,
                    by = c("Gene.stable.ID" = "Gene")) %>% 
          filter(!is.na(adj.P.Val)))  

# save
DEGs.cluster.genotype.GO.terms.pvalue %>% 
  select(c(Gene.stable.ID,
           GO.term.accession,
           Gene.name,
           Keep,
           adj.P.Val)) %>% 
  distinct() %>% 
  relocate(GO.term.accession, 
           .after = last_col())%>% 
  relocate(adj.P.Val, 
           .after = last_col()) %>% 
  arrange(Keep) %>% 
  write_csv('neuropeptides/avp.oxt/avp/limmatrend/genotype/DEGs.cluster.GO.terms.pvalue.csv')

# to look for enrichment and create figures use GO term accession and kME (higher is better) for each module: http://revigo.irb.hr/
# remember to sort by kME and filter out blank GO terms before copying to revigo
# all module genes












#### limmatrend genotype, and status approach ####
### calculate variable genes
## identify top 500 variable genes
burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.variable <- FindVariableFeatures(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster, 
                                                                                         selection.method = "vst", 
                                                                                         nfeatures = 500, 
                                                                                         verbose = F)
# identify top 500 variable genes
avp.neuron.reduce.group.topgenes <- head(VariableFeatures(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.variable), 
                                         500)

## create vector of factor
avp.neuron.reduce.DomvsSub.vector.list = burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.expression %>% 
  mutate(avp.neuron.orig.ident = orig.ident %>% 
           as.factor()) %>% 
  pull(avp.neuron.orig.ident) %>% 
  droplevels()

### counts matrix
## raw read count matrix
## rows = genes, columns = cells
# only keep 500 variable genes
# set negative values to 0
avp.neuron.reduce.DomvsSub.vector.count = GetAssayData(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster) %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column('gene') %>% 
  select(c(gene,
           burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.expression %>% 
             pull(Cell.id))) %>% 
  filter(gene %in% avp.neuron.reduce.group.topgenes) %>% 
  column_to_rownames('gene') %>% 
  as.matrix() %>% 
  pmax(0)


## run with reduce avp neurons
#need list with L with  count and condt
# add in cell cluster with condt.2
## create vector of factor
avp.neuron.reduce.DomvsSub.vector.list.cluster = burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.expression %>% 
  mutate(avp.neuron.cluster = integrated_snn_res.0.2 %>% 
           as.factor()) %>% 
  pull(avp.neuron.cluster) %>% 
  droplevels


# add in genotype
## create vector of factor
avp.neuron.reduce.DomvsSub.vector.list.genotype = burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.expression %>% 
  mutate(genotype = Genotype.id %>% 
           as.factor()) %>% 
  pull(genotype) %>% 
  droplevels()

#create list
avp.neuron.reduce.vector.limma = list(count = avp.neuron.reduce.DomvsSub.vector.count,
                                      condt = avp.neuron.reduce.DomvsSub.vector.list,
                                      condt.2 = avp.neuron.reduce.DomvsSub.vector.list.cluster,
                                      genotype = avp.neuron.reduce.DomvsSub.vector.list.genotype)
# add genotype.ID?



#load libraries
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(edgeR))

#create function
run_limmatrend_complex_genotype <- function(L) {
  message("limmatrend")
  session_info <- sessionInfo()
  timing <- system.time({
    treat <- L$condt
    genotype <- L$genotype
    cluster <- L$condt.2
    design <- model.matrix(~0+treat+genotype+cluster) 
    contrasts <- makeContrasts(DvsS = treatdom_burtoni_snseq - treatsub_burtoni_snseq,
                               levels = design)
    dge <- DGEList(L$count, 
                   group = treat)
    dge <- calcNormFactors(dge)
    
    y <- new("EList")
    y$E <- edgeR::cpm(dge, 
                      log = TRUE, 
                      prior.count = 3)
    fit <- lmFit(y, 
                 design = design)
    fit <- contrasts.fit(fit , 
                         contrasts)
    fit <- eBayes(fit, 
                  trend = TRUE,
                  robust = TRUE)
    ttDS <- topTable(fit, 
                     n = Inf,
                     coef = "DvsS",
                     adjust.method = "BH")
  })
  # Open pdf file
  pdf(file= "./neuropeptides/avp.oxt/avp/limmatrend/genotype.DvS/limmatrend.histograms.pdf" )
  # create a 2X2 grid
  par( mfrow= c(2,2) )
  #graph
  hist(ttDS$P.Value, 50)
  hist(ttDS$adj.P.Val, 50)
  dev.off()
  
  # Open pdf file
  pdf(file= "./neuropeptides/avp.oxt/avp/limmatrend/genotype.DvS/limmatrend.MDS.pdf" )
  # create a 2X1 grid
  par( mfrow= c(2,1) )
  limma::plotMDS(dge, 
                 col = as.numeric(as.factor(L$condt)), 
                 pch = 19)
  plotMD(fit)
  limma::plotMDS(dge, 
                 col = as.numeric(as.factor(L$condt.2)), 
                 pch = 19)
  plotMD(fit)
  limma::plotMDS(dge, 
                 col = as.numeric(as.factor(L$genotype)), 
                 pch = 19)
  plotMD(fit)
  dev.off()
  
  #print results
  list(session_info = session_info,
       timing = timing,
       ttDS = ttDS)
}


###run function  
avp.neuron.reduce.genotype.vector.limma.results = run_limmatrend_complex_genotype(avp.neuron.reduce.vector.limma)

# save results to dataframe
avp.neuron.reduce.genotype.vector.limma.results.df = avp.neuron.reduce.genotype.vector.limma.results$ttDS %>% 
  rename_with(~paste0(.,"_DvsS")) %>% 
  rownames_to_column("Gene")



#add color for significance  
avp.neuron.reduce.genotype.vector.limma.results.df = avp.neuron.reduce.genotype.vector.limma.results.df %>% 
  mutate(Sig_DvsS = ifelse(adj.P.Val_DvsS <= 0.05,
                           'Sig',
                           'Not Sig'),
         Direction.type_DvsS = ifelse(logFC_DvsS > 0,
                                      'up',
                                      'down'),
         Sig.direction_DvsS = ifelse(Sig_DvsS == 'Sig',
                                     Direction.type_DvsS,
                                     Sig_DvsS))


##graph volcano plot
# create comparison list
limma.genotpye.vector = c("DvsS")

for (i in limma.genotpye.vector) {
  # graph volcano plot
  avp.neuron.reduce.genotype.vector.limma.results.df %>% 
    mutate(sig.label = ifelse(get(paste0("Sig_", i)) == 'Sig',
                              Gene,
                              '')) %>% 
    ggplot(aes(x = get(paste0("logFC_", i)),
               y = -log10(get(paste0("adj.P.Val_", i))),
               color = get(paste0("Sig.direction_", i))))+
    geom_hline(yintercept = -log10(0.05),
               linetype = 'dotted') +
    geom_point(size = 5) +
    theme_classic() + 
    scale_color_manual(values=c("21B9CA", 
                                "grey", 
                                "orange3")) +
    geom_text(aes(label = sig.label),
              vjust = 0, 
              nudge_y = 0.10,
              size = 5) +
    theme(text = element_text(size = 20),
          legend.position = 'none') +
    xlab(paste0("logFC_", i)) +
    ylab( paste0("-log10(adj.P.Val_", i,")")) +
    ggtitle(paste0(i, " volcano plot"))
  ggsave(paste0('neuropeptides/avp.oxt/avp/limmatrend/genotype.DvS//limma.avp.neurons.reduce.volcano.', i, '.png'),
         width = 5,
         height = 5)
}

## get list of significant genes
DEGs.cluster.genotype = avp.neuron.reduce.genotype.vector.limma.results.df %>% 
  mutate(Keep = 0,
         Keep = case_when(Sig.direction_DvsS != "Not Sig" ~ Keep + 1,
                          TRUE ~ Keep)) %>% 
  filter(Keep != 0) %>%  
  select(Gene,
         Sig.direction_DvsS,
         adj.P.Val_DvsS) %>% 
  mutate(Keep = case_when(Sig.direction_DvsS == "up"  ~ "Dom",
                          TRUE ~ "Sub")) %>% 
  select(-c(Sig.direction_DvsS)) %>% 
  group_by(Gene) %>% 
  mutate(adj.P.Val = min(c(adj.P.Val_DvsS))) %>% 
  select(-c(adj.P.Val_DvsS)) 
write_csv(DEGs.cluster.genotype,
          'neuropeptides/avp.oxt/avp/limmatrend/genotype.DvS/DEGs.cluster.csv')

# save hub genes for GO analysis with:
# to get GO term accession use biomart: http://useast.ensembl.org/biomart/martview/

## combine go terms
DEGs.cluster.genotype.GO.terms = read.csv('neuropeptides/avp.oxt/avp/limmatrend/genotype.DvS//DEGs.cluster.GO.terms.csv')

# join with pvalue 
DEGs.cluster.genotype.GO.terms.pvalue = DEGs.cluster.genotype.GO.terms %>% 
  left_join(DEGs.cluster.genotype ,
            by = c("Gene.name" = "Gene")) %>% 
  filter(!is.na(adj.P.Val)) %>% 
  rbind(DEGs.cluster.genotype.GO.terms %>% 
          left_join(DEGs.cluster.genotype ,
                    by = c("Gene.stable.ID" = "Gene")) %>% 
          filter(!is.na(adj.P.Val)))  

# save
DEGs.cluster.genotype.GO.terms.pvalue %>% 
  select(c(Gene.stable.ID,
           GO.term.accession,
           Gene.name,
           Keep,
           adj.P.Val)) %>% 
  distinct() %>% 
  filter(GO.term.accession != '') %>% 
  relocate(GO.term.accession, 
           .after = last_col())%>% 
  relocate(adj.P.Val, 
           .after = last_col()) %>% 
  arrange(Keep) %>%
  write_csv('neuropeptides/avp.oxt/avp/limmatrend/genotype.DvS//DEGs.cluster.GO.terms.pvalue.csv')

# to look for enrichment and create figures use GO term accession and kME (higher is better) for each module: http://revigo.irb.hr/
# remember to sort by kME and filter out blank GO terms before copying to revigo
# all module genes

#### avp analysis #### 

## violin plot
#cluster and genotype
VlnPlot(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster,
        features = c("avp"),
        split.by = 'Genotype.id',
        group.by = 'integrated_snn_res.0.2',
        slot = 'scale.data')
ggsave('neuropeptides/avp.oxt/avp/AVP expression per cluster and genotype.png',
       height = 10,
       width = 10)

#cluster and status
VlnPlot(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster,
        features = c("avp"),
        split.by = 'orig.ident',
        group.by = 'integrated_snn_res.0.2',
        slot = 'scale.data')
ggsave('neuropeptides/avp.oxt/avp/AVP expression per cluster and social status.png',
       height = 10,
       width = 10)

#for poster
# VlnPlot(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster,
#         features = c("avp"),
#         split.by = 'orig.ident',
#         group.by = 'integrated_snn_res.0.2',
#         slot = 'scale.data',
#         pt.size = 0.5)+
#   scale_fill_manual(values = c("#4e499e",
#                                "#60bb46"))+
#   theme(legend.position = "none") +
#   theme(panel.border = element_rect(color = "black",
#                                     fill = NA,
#                                     size = 1))+
#   theme(axis.text = element_text(size = 15))  +
#   theme(axis.title = element_text(size = 20)) +
#   theme(plot.title = element_blank())+
#   ylab('AVP expression level') +
#   xlab('Cluster')
# ggsave('neuropeptides/avp.oxt/avp/AVP expression per cluster and social status.pdf',
#        height = 5,
#        width = 5,
#        units = "in",
#        dpi = 320)

#for poster
# update
VlnPlot(burtoni.snseq.combined.sct.all.avp.neurons.recluster,
        features = c("avp"),
        split.by = 'orig.ident',
        group.by = 'integrated_snn_res.0.8',
        slot = 'data',
        pt.size = 0.5,
        assay = 'SCT')+
  scale_fill_manual(values = c("#4e499e",
                               "#60bb46"))+
  theme(legend.position = "none") +
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1))+
  theme(axis.text = element_text(size = 15))  +
  theme(axis.title = element_text(size = 20)) +
  theme(plot.title = element_blank())+
  ylab('AVP expression level') +
  xlab('AVP neuron custer')
ggsave('neuropeptides/avp.oxt/avp/AVP expression per cluster and social status update.pdf',
       height = 5,
       width = 5,
       units = "in",
       dpi = 320)



#get average expression data
avp.avg.exp.cluster.genotype = AverageExpression(AddMetaData(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster,
                              metadata = paste(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster@meta.data$Genotype.id,
                                               burtoni.snseq.combined.sct.reduce.avp.neurons.recluster@meta.data$integrated_snn_res.0.2,
                                               sep = '_'),
                              col.name = 'Genotype.cluster'),
  assays = 'integrated',
  features = 'avp',
  return.seurat = FALSE,
  group.by = "Genotype.cluster",
  slot = "scale.data",
  verbose = TRUE) %>% 
  as.data.frame() %>% 
  pivot_longer(cols = everything(),
               names_to = 'ID',
               values_to = 'Avg.expression') %>% 
  separate(ID,
           sep = '_',
           into = c('Test',
                    'Cluster')) %>% 
  separate(Test,
           sep = "\\.",
           into = c('Data.type',
                    'orig.ident',
                    'Name')) %>%  
  mutate(Genotype.id = paste(orig.ident,
                              Name,
                             sep = '.')) %>% 
  select(-c(Data.type)) 



##graph AVP per cluster across social status
avp.avg.exp.cluster.genotype %>%
  ggplot(aes(x = Cluster,
                    y = Avg.expression,
                    color = orig.ident,
             shape = Name,
             group = orig.ident)) +
  geom_point(position = position_dodge(width = 0.3)) +
  stat_summary(
    fun = mean, 
    geom = "errorbar", 
    aes(ymax = ..y.., ymin = ..y..), 
    position = position_dodge(width = 0.3), 
    width = 0.25) +
  theme_classic() +
  ylim(c(0,max(avp.avg.exp.cluster.genotype$Avg.expression)))
ggsave('neuropeptides/avp.oxt/avp/Average AVP expression per cluster and genotype.png',
       height = 10,
       width = 10)

#for poster
##graph AVP per cluster across social status
# avp.avg.exp.cluster.genotype %>%
#   ggplot(aes(x = Cluster,
#              y = Avg.expression,
#              color = orig.ident,
#              group = orig.ident)) +
#   stat_summary(
#     fun = mean,
#     geom = "errorbar",
#     aes(ymax = ..y.., ymin = ..y..),
#     position = position_dodge(width = 0.75),
#     width = 0.5,
#     size = 5,
#     color = 'black') +
#   geom_point(position = position_dodge(width = 0.75),
#              shape = 21,
#              color = 'black',
#              aes(fill = orig.ident),
#              size = 8) +
#   ylim(c(0,max(avp.avg.exp.cluster.genotype$Avg.expression)))+
#   theme_classic() +
#   theme(legend.position = 'none') +
#   scale_fill_manual(values = c("#4e499e",
#                                "#60bb46"))+
#   scale_color_manual(values = c("#4e499e",
#                                "#60bb46"))+
#   ylab("AVP average expression") +
#   theme(panel.border = element_rect(color = "black",
#                                     fill = NA,
#                                     size = 1))+
#   theme(axis.text = element_text(size = 15))  +
#   theme(axis.title = element_text(size = 20))
# ggsave('neuropeptides/avp.oxt/avp/Average AVP expression per cluster and genotype poster.pdf',
#        height = 5,
#        width = 5,
#        units = "in",
#        dpi = 320)



#get average expression data
# update
avp.avg.exp.cluster.genotype = AverageExpression(AddMetaData(burtoni.snseq.combined.sct.all.avp.neurons.recluster,
                                                             metadata = paste(burtoni.snseq.combined.sct.all.avp.neurons.recluster@meta.data$Genotype.id,
                                                                              burtoni.snseq.combined.sct.all.avp.neurons.recluster@meta.data$integrated_snn_res.0.8,
                                                                              sep = '_'),
                                                             col.name = 'Genotype.cluster'),
                                                 assays = 'SCT',
                                                 features = 'avp',
                                                 return.seurat = FALSE,
                                                 group.by = "Genotype.cluster",
                                                 slot = "data",
                                                 verbose = TRUE) %>% 
  as.data.frame() %>% 
  pivot_longer(cols = everything(),
               names_to = 'ID',
               values_to = 'Avg.expression') %>% 
  separate(ID,
           sep = '_',
           into = c('Test',
                    'Cluster')) %>% 
  separate(Test,
           sep = "\\.",
           into = c('Data.type',
                    'orig.ident',
                    'Name')) %>%  
  mutate(Genotype.id = paste(orig.ident,
                             Name,
                             sep = '.')) %>% 
  select(-c(Data.type)) 


#for poster
##graph AVP per cluster across social status
avp.avg.exp.cluster.genotype %>%
  ggplot(aes(x = Cluster,
             y = Avg.expression,
             color = orig.ident,
             group = orig.ident)) +
  stat_summary(
    fun = mean,
    geom = "errorbar",
    aes(ymax = ..y.., ymin = ..y..),
    position = position_dodge(width = 0.75),
    width = 0.5,
    size = 5,
    color = 'black') +
  geom_point(position = position_dodge(width = 0.75),
             shape = 21,
             color = 'black',
             aes(fill = orig.ident),
             size = 8) +
  ylim(c(0,max(avp.avg.exp.cluster.genotype$Avg.expression)))+
  theme_classic() +
  theme(legend.position = 'none') +
  scale_fill_manual(values = c("#4e499e",
                               "#60bb46"))+
  scale_color_manual(values = c("#4e499e",
                               "#60bb46"))+
  ylab("AVP average expression") +
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1))+
  theme(axis.text = element_text(size = 15))  +
  theme(axis.title = element_text(size = 20)) +
  xlab('AVP neuron custer')
ggsave('neuropeptides/avp.oxt/avp/Average AVP expression per cluster and genotype poster update.pdf',
       height = 5,
       width = 5,
       units = "in",
       dpi = 320)

#get average expression data
# update
avp.avg.exp.cluster.genotype.read = AverageExpression(AddMetaData(burtoni.snseq.combined.sct.all.avp.neurons.recluster,
                                                             metadata = paste(burtoni.snseq.combined.sct.all.avp.neurons.recluster@meta.data$Genotype.id,
                                                                              burtoni.snseq.combined.sct.all.avp.neurons.recluster@meta.data$integrated_snn_res.0.8,
                                                                              sep = '_'),
                                                             col.name = 'Genotype.cluster'),
                                                 assays = 'RNA',
                                                 features = 'avp',
                                                 return.seurat = FALSE,
                                                 group.by = "Genotype.cluster",
                                                 slot = "counts",
                                                 verbose = TRUE) %>% 
  as.data.frame() %>% 
  pivot_longer(cols = everything(),
               names_to = 'ID',
               values_to = 'Avg.expression') %>% 
  separate(ID,
           sep = '_',
           into = c('Test',
                    'Cluster')) %>% 
  separate(Test,
           sep = "\\.",
           into = c('Data.type',
                    'orig.ident',
                    'Name')) %>%  
  mutate(Genotype.id = paste(orig.ident,
                             Name,
                             sep = '.')) %>% 
  select(-c(Data.type)) 


#for poster
##graph AVP per cluster across social status
avp.avg.exp.cluster.genotype.read %>%
  ggplot(aes(x = Cluster,
             y = Avg.expression,
             color = orig.ident,
             group = orig.ident)) +
  stat_summary(
    fun = mean,
    geom = "errorbar",
    aes(ymax = ..y.., ymin = ..y..),
    position = position_dodge(width = 0.75),
    width = 0.5,
    size = 5,
    color = 'black') +
  geom_point(position = position_dodge(width = 0.75),
             shape = 21,
             color = 'black',
             aes(fill = orig.ident),
             size = 8) +
  ylim(c(0,max(avp.avg.exp.cluster.genotype.read$Avg.expression)))+
  theme_classic() +
  theme(legend.position = 'none') +
  scale_fill_manual(values = c("#4e499e",
                               "#60bb46"))+
  scale_color_manual(values = c("#4e499e",
                                "#60bb46"))+
  ylab("AVP average expression") +
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1))+
  theme(axis.text = element_text(size = 15))  +
  theme(axis.title = element_text(size = 20)) +
  xlab('AVP neuron custer')
ggsave('neuropeptides/avp.oxt/avp/Average AVP expression per cluster and genotype poster update.pdf',
       height = 5,
       width = 5,
       units = "in",
       dpi = 320)

#### alluvial ####
### avp neurons
## status, cell type, cluster
# color cluster
burtoni.snseq.combined.sct.reduce.avp.neurons.recluster@meta.data %>% 
  select(c(orig.ident,
           Cell.type,
           integrated_snn_res.0.2)) %>% 
  table() %>% 
  as.data.frame() %>% 
  ggplot(aes(axis1 = orig.ident, 
             axis2 = Cell.type,
             axis3 = integrated_snn_res.0.2,
             y = Freq)) +
  geom_alluvium(aes(fill = integrated_snn_res.0.2)) +
  geom_stratum() +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("orig.ident", "integrated_snn_res.0.2"),
                   expand = c(0.15, 0.05)) +
  scale_fill_viridis_d() +
  theme_classic() +
  theme(axis.text.x=element_blank())
ggsave('./neuropeptides/avp.oxt/avp/alluvial/AVP neuron status celltype by cluster.png',
       height = 10,
       width = 10)

# color status
burtoni.snseq.combined.sct.reduce.avp.neurons.recluster@meta.data %>% 
  select(c(orig.ident,
           Cell.type,
           integrated_snn_res.0.2)) %>% 
  table() %>% 
  as.data.frame() %>% 
  group_by(orig.ident) %>% 
  mutate(Count = sum(Freq)) %>% 
  ungroup() %>% 
  mutate(Freq.scale = 100*Freq/Count) %>% 
  ggplot(aes(axis1 = orig.ident, 
             axis2 = Cell.type,
             axis3 = integrated_snn_res.0.2,
             y = Freq.scale)) +
  geom_alluvium(aes(fill = orig.ident)) +
  geom_stratum() +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("orig.ident", "integrated_snn_res.0.2"),
                   expand = c(0.15, 0.05)) +
  scale_fill_viridis_d() +
  theme_classic() +
  theme(axis.text.x=element_blank())
ggsave('./neuropeptides/avp.oxt/avp/alluvial/AVP neuron cluster celltype by status.png',
       height = 10,
       width = 10)

## status cluster
# color status
burtoni.snseq.combined.sct.reduce.avp.neurons.recluster@meta.data %>% 
  select(c(orig.ident,
           Cell.type,
           integrated_snn_res.0.2)) %>% 
  table() %>% 
  as.data.frame() %>% 
  ggplot(aes(axis1 = orig.ident, 
             axis2 = integrated_snn_res.0.2,
             y = Freq)) +
  geom_alluvium(aes(fill = orig.ident)) +
  geom_stratum() +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("orig.ident", "integrated_snn_res.0.2"),
                   expand = c(0.15, 0.05)) +
  scale_fill_viridis_d() +
  theme_classic() +
  theme(axis.text.x=element_blank())
ggsave('./neuropeptides/avp.oxt/avp/alluvial/AVP neuron cluster by status.png',
       height = 10,
       width = 10)

# color cluster
burtoni.snseq.combined.sct.reduce.avp.neurons.recluster@meta.data %>% 
  select(c(orig.ident,
           Cell.type,
           integrated_snn_res.0.2)) %>% 
  table() %>% 
  as.data.frame() %>% 
  ggplot(aes(axis1 = orig.ident, 
             axis2 = integrated_snn_res.0.2,
             y = Freq)) +
  geom_alluvium(aes(fill = integrated_snn_res.0.2)) +
  geom_stratum() +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("orig.ident", "integrated_snn_res.0.2"),
                   expand = c(0.15, 0.05)) +
  scale_fill_viridis_d() +
  theme_classic() +
  theme(axis.text.x=element_blank())
ggsave('./neuropeptides/avp.oxt/avp/alluvial/AVP neuron status by cluster.png',
       height = 10,
       width = 10)



#### limmatrend genotype, clusters, and status approach with SCT ####
#### Variable have same name as 'limmatrend genotype, clusters, and status approach'
### create new expression matrix with SCT
##expression level
burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.expression.sct = full_join(full_join(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster@reductions$umap@cell.embeddings %>% 
                                                                                           as.data.frame() %>% 
                                                                                           rownames_to_column("Cell.id"),
                                                                                         burtoni.snseq.combined.sct.reduce.avp.neurons.recluster@meta.data %>% 
                                                                                           rownames_to_column("Cell.id")),
                                                                               burtoni.snseq.combined.sct.reduce.avp.neurons.recluster@assays$SCT@scale.data %>% 
                                                                                 as.data.frame() %>% 
                                                                                 filter(rownames(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster@assays$SCT@scale.data) %in% c('avp',
                                                                                                                                                                                              'oxt')) %>% 
                                                                                 t() %>% as.data.frame() %>% 
                                                                                 rownames_to_column('Cell.id')) %>% 
  full_join(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster@reductions$pca@cell.embeddings %>% 
              as.data.frame() %>% 
              select(PC_1,
                     PC_2,
                     PC_3,
                     PC_4,
                     PC_5) %>% 
              rownames_to_column("Cell.id")) %>% 
  mutate(Oxt.cell = ifelse(oxt > 2 & avp > 2,
                           'both',
                           ifelse(oxt > 2,
                                  'oxt',
                                  'avp')))

### calculate variable genes
## identify top 500 variable genes
burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.variable <- FindVariableFeatures(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster, 
                                                                                         selection.method = "vst", 
                                                                                         nfeatures = 500, 
                                                                                         verbose = F)
# identify top 500 variable genes
avp.neuron.reduce.group.topgenes <- head(VariableFeatures(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.variable), 
                                         500)

## create vector of factor
avp.neuron.reduce.DomvsSub.vector.list.sct = burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.expression.sct %>% 
  mutate(avp.neuron.orig.ident = orig.ident %>% 
           as.factor()) %>% 
  pull(avp.neuron.orig.ident) %>% 
  droplevels()

### counts matrix
## raw read count matrix
## rows = genes, columns = cells
# only keep 500 variable genes
# set negative values to 0
avp.neuron.reduce.DomvsSub.vector.count.sct = GetAssayData(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster,
                                                           assay = 'SCT') %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column('gene') %>% 
  select(c(gene,
           burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.expression.sct %>% 
             pull(Cell.id))) %>% 
  filter(gene %in% avp.neuron.reduce.group.topgenes) %>% 
  column_to_rownames('gene') %>% 
  as.matrix() %>% 
  pmax(0)


## run with reduce avp neurons
#need list with L with  count and condt
# add in cell cluster with condt.2
## create vector of factor
avp.neuron.reduce.DomvsSub.vector.list.cluster.sct = burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.expression.sct %>% 
  mutate(avp.neuron.cluster = integrated_snn_res.0.2 %>% 
           as.factor()) %>% 
  pull(avp.neuron.cluster) %>% 
  droplevels


# add in genotype
## create vector of factor
avp.neuron.reduce.DomvsSub.vector.list.genotype.sct = burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.expression.sct %>% 
  mutate(genotype = Genotype.id %>% 
           as.factor()) %>% 
  pull(genotype) %>% 
  droplevels()

#create list
avp.neuron.reduce.vector.limma.sct = list(count = avp.neuron.reduce.DomvsSub.vector.count.sct,
                                      condt = avp.neuron.reduce.DomvsSub.vector.list.sct,
                                      condt.2 = avp.neuron.reduce.DomvsSub.vector.list.cluster.sct,
                                      genotype = avp.neuron.reduce.DomvsSub.vector.list.genotype.sct)
# add genotype.ID?



#load libraries
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(edgeR))

#create function
run_limmatrend_complex_genotype.sct <- function(L) {
  message("limmatrend")
  session_info <- sessionInfo()
  timing <- system.time({
    treat <- paste(L$condt, 
                   L$condt.2,
                   sep = '.')
    genotype <- L$genotype
    design <- model.matrix(~0+treat+genotype) 
    contrasts <- makeContrasts(DvsS_0 = treatdom_burtoni_snseq.0 - treatsub_burtoni_snseq.0, 
                               DvsS_1 = treatdom_burtoni_snseq.1 - treatsub_burtoni_snseq.1,
                               DvsS_2 = treatdom_burtoni_snseq.2 - treatsub_burtoni_snseq.2,
                               DvsS_3 = treatdom_burtoni_snseq.3 - treatsub_burtoni_snseq.3, 
                               DvsS_4 = treatdom_burtoni_snseq.4 - treatsub_burtoni_snseq.4,
                               DvsS=(treatdom_burtoni_snseq.0+treatdom_burtoni_snseq.1+treatdom_burtoni_snseq.2+treatdom_burtoni_snseq.3+treatdom_burtoni_snseq.4)/5-(treatsub_burtoni_snseq.0+treatsub_burtoni_snseq.1+treatsub_burtoni_snseq.2+treatsub_burtoni_snseq.3+treatsub_burtoni_snseq.4)/5,
                               levels = design)
    dge <- DGEList(L$count, 
                   group = treat)
    dge <- calcNormFactors(dge)
    
    y <- new("EList")
    y$E <- edgeR::cpm(dge, 
                      log = TRUE, 
                      prior.count = 3)
    fit <- lmFit(y, 
                 design = design)
    fit <- contrasts.fit(fit , 
                         contrasts)
    fit <- eBayes(fit, 
                  trend = TRUE,
                  robust = TRUE)
    tt0 <- topTable(fit, 
                    n = Inf,
                    coef = "DvsS_0",
                    adjust.method = "BH")
    tt1 <- topTable(fit, 
                    n = Inf,
                    coef = "DvsS_1",
                    adjust.method = "BH")
    tt2 <- topTable(fit, 
                    n = Inf,
                    coef = "DvsS_2",
                    adjust.method = "BH")
    tt3 <- topTable(fit, 
                    n = Inf,
                    coef = "DvsS_3",
                    adjust.method = "BH")
    tt4 <- topTable(fit, 
                    n = Inf,
                    coef = "DvsS_4",
                    adjust.method = "BH")
    ttDS <- topTable(fit, 
                     n = Inf,
                     coef = "DvsS",
                     adjust.method = "BH")
  })
  # Open pdf file
  pdf(file= "./neuropeptides/avp.oxt/avp/limmatrend/genotype.sct/limmatrend.histograms.pdf" )
  # create a 2X2 grid
  par( mfrow= c(2,2) )
  #graph
  hist(tt0$P.Value, 50)
  hist(tt0$adj.P.Val, 50)
  hist(tt1$P.Value, 50)
  hist(tt1$adj.P.Val, 50)
  hist(tt2$P.Value, 50)
  hist(tt2$adj.P.Val, 50)
  hist(tt3$P.Value, 50)
  hist(tt3$adj.P.Val, 50)
  hist(tt4$P.Value, 50)
  hist(tt4$adj.P.Val, 50)
  hist(ttDS$P.Value, 50)
  hist(ttDS$adj.P.Val, 50)
  dev.off()
  
  # Open pdf file
  pdf(file= "./neuropeptides/avp.oxt/avp/limmatrend/genotype.sct/limmatrend.MDS.pdf" )
  # create a 2X1 grid
  par( mfrow= c(2,1) )
  limma::plotMDS(dge, 
                 col = as.numeric(as.factor(L$condt)), 
                 pch = 19)
  plotMD(fit)
  limma::plotMDS(dge, 
                 col = as.numeric(as.factor(L$condt.2)), 
                 pch = 19)
  plotMD(fit)
  dev.off()
  
  #print results
  list(session_info = session_info,
       timing = timing,
       tt0 = tt0,
       tt1 = tt1,
       tt2 = tt2,
       tt3 = tt3,
       tt4 = tt4,
       ttDS = ttDS)
}


###run function  
avp.neuron.reduce.genotype.vector.limma.results.sct = run_limmatrend_complex_genotype.sct(avp.neuron.reduce.vector.limma.sct)

# save results to dataframe
avp.neuron.reduce.genotype.vector.limma.results.df.sct = full_join(avp.neuron.reduce.genotype.vector.limma.results.sct$tt0 %>% 
                                                                 rename_with(~paste0(.,"_DvsS_0")) %>% 
                                                                 rownames_to_column("Gene"),
                                                               avp.neuron.reduce.genotype.vector.limma.results.sct$tt1 %>% 
                                                                 rename_with(~paste0(.,"_DvsS_1")) %>% 
                                                                 rownames_to_column("Gene")) %>% 
  full_join(avp.neuron.reduce.genotype.vector.limma.results.sct$tt2 %>% 
              rename_with(~paste0(.,"_DvsS_2")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(avp.neuron.reduce.genotype.vector.limma.results.sct$tt3 %>% 
              rename_with(~paste0(.,"_DvsS_3")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(avp.neuron.reduce.genotype.vector.limma.results.sct$tt4 %>% 
              rename_with(~paste0(.,"_DvsS_4")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(avp.neuron.reduce.genotype.vector.limma.results.sct$ttDS %>% 
              rename_with(~paste0(.,"_DvsS")) %>% 
              rownames_to_column("Gene")) 



#add color for significance  
avp.neuron.reduce.genotype.vector.limma.results.df.sct = avp.neuron.reduce.genotype.vector.limma.results.df.sct %>% 
  mutate(Sig_DvsS_0 = ifelse(adj.P.Val_DvsS_0 <= 0.05,
                             'Sig',
                             'Not Sig'),
         Direction.type_DvsS_0 = ifelse(logFC_DvsS_0 > 0,
                                        'up',
                                        'down'),
         Sig.direction_DvsS_0 = ifelse(Sig_DvsS_0 == 'Sig',
                                       Direction.type_DvsS_0,
                                       Sig_DvsS_0)) %>% 
  mutate(Sig_DvsS_1 = ifelse(adj.P.Val_DvsS_1 <= 0.05,
                             'Sig',
                             'Not Sig'),
         Direction.type_DvsS_1 = ifelse(logFC_DvsS_1 > 0,
                                        'up',
                                        'down'),
         Sig.direction_DvsS_1 = ifelse(Sig_DvsS_1 == 'Sig',
                                       Direction.type_DvsS_1,
                                       Sig_DvsS_1)) %>% 
  mutate(Sig_DvsS_2 = ifelse(adj.P.Val_DvsS_2 <= 0.05,
                             'Sig',
                             'Not Sig'),
         Direction.type_DvsS_2 = ifelse(logFC_DvsS_2 > 0,
                                        'up',
                                        'down'),
         Sig.direction_DvsS_2 = ifelse(Sig_DvsS_2 == 'Sig',
                                       Direction.type_DvsS_2,
                                       Sig_DvsS_2)) %>% 
  mutate(Sig_DvsS_3 = ifelse(adj.P.Val_DvsS_3 <= 0.05,
                             'Sig',
                             'Not Sig'),
         Direction.type_DvsS_3 = ifelse(logFC_DvsS_3 > 0,
                                        'up',
                                        'down'),
         Sig.direction_DvsS_3 = ifelse(Sig_DvsS_3 == 'Sig',
                                       Direction.type_DvsS_3,
                                       Sig_DvsS_3)) %>% 
  mutate(Sig_DvsS_4 = ifelse(adj.P.Val_DvsS_4 <= 0.05,
                             'Sig',
                             'Not Sig'),
         Direction.type_DvsS_4 = ifelse(logFC_DvsS_4 > 0,
                                        'up',
                                        'down'),
         Sig.direction_DvsS_4 = ifelse(Sig_DvsS_4 == 'Sig',
                                       Direction.type_DvsS_4,
                                       Sig_DvsS_4)) %>% 
  mutate(Sig_DvsS = ifelse(adj.P.Val_DvsS <= 0.05,
                           'Sig',
                           'Not Sig'),
         Direction.type_DvsS = ifelse(logFC_DvsS > 0,
                                      'up',
                                      'down'),
         Sig.direction_DvsS = ifelse(Sig_DvsS == 'Sig',
                                     Direction.type_DvsS,
                                     Sig_DvsS))


##graph volcano plot
# create comparison list
limma.genotpye.vector = c("DvsS_0",
                          "DvsS_1",
                          "DvsS_2",
                          "DvsS_3",
                          "DvsS_4",
                          "DvsS")

for (i in limma.genotpye.vector) {
  # graph volcano plot
  avp.neuron.reduce.genotype.vector.limma.results.df.sct %>% 
    mutate(sig.label = ifelse(get(paste0("Sig_", i)) == 'Sig',
                              Gene,
                              '')) %>% 
    ggplot(aes(x = get(paste0("logFC_", i)),
               y = -log10(get(paste0("adj.P.Val_", i))),
               color = get(paste0("Sig.direction_", i))))+
    geom_hline(yintercept = -log10(0.05),
               linetype = 'dotted') +
    geom_point(size = 5) +
    theme_classic() + 
    scale_color_manual(values=c("21B9CA", 
                                "grey", 
                                "orange3")) +
    geom_text(aes(label = sig.label),
              vjust = 0, 
              nudge_y = 0.10,
              size = 5) +
    theme(text = element_text(size = 20),
          legend.position = 'none') +
    xlab(paste0("logFC_", i)) +
    ylab( paste0("-log10(adj.P.Val_", i,")")) +
    ggtitle(paste0(i, " volcano plot"))
  ggsave(paste0('neuropeptides/avp.oxt/avp/limmatrend/genotype.sct/limma.avp.neurons.reduce.volcano.', i, '.png'),
         width = 5,
         height = 5)
}

## get list of significant genes
DEGs.cluster.genotype.sct = avp.neuron.reduce.genotype.vector.limma.results.df.sct %>% 
mutate(Keep = 0,
       Keep.2 = case_when(Sig.direction_DvsS_0 != "Not Sig" ~ Keep + 1,
                          Sig.direction_DvsS_1 != "Not Sig" ~ Keep + 1,
                          Sig.direction_DvsS_2 != "Not Sig" ~ Keep + 1,
                          Sig.direction_DvsS_3 != "Not Sig" ~ Keep + 1,
                          Sig.direction_DvsS_4 != "Not Sig" ~ Keep + 1,
                          Sig.direction_DvsS != "Not Sig" ~ Keep + 1,
                          TRUE ~ Keep)) %>% 
  filter(Keep.2 != 0) %>% 
  select(Gene,
         Sig.direction_DvsS_0,
         Sig.direction_DvsS_1,
         Sig.direction_DvsS_2,
         Sig.direction_DvsS_3,
         Sig.direction_DvsS_4,
         Sig.direction_DvsS,
         adj.P.Val_DvsS_0,
         adj.P.Val_DvsS_1,
         adj.P.Val_DvsS_2,
         adj.P.Val_DvsS_3,
         adj.P.Val_DvsS_4,
         adj.P.Val_DvsS) %>% 
  mutate(Keep = case_when(Sig.direction_DvsS_0 == "up" | Sig.direction_DvsS_1 == "up" | Sig.direction_DvsS_2 == "up" | Sig.direction_DvsS_3 == "up" | Sig.direction_DvsS_4 == "up" ~ "Dom"| Sig.direction_DvsS == "up" ~ "Dom",
                          TRUE ~ "Sub")) %>% 
  group_by(Gene) %>% 
  mutate(adj.P.Val = min(c(adj.P.Val_DvsS_0,
                           adj.P.Val_DvsS_1,
                           adj.P.Val_DvsS_2,
                           adj.P.Val_DvsS_3,
                           adj.P.Val_DvsS_4,
                           adj.P.Val_DvsS))) %>% 
  select(-c(adj.P.Val_DvsS_0,
            adj.P.Val_DvsS_1,
            adj.P.Val_DvsS_2,
            adj.P.Val_DvsS_3,
            adj.P.Val_DvsS_4,
            adj.P.Val_DvsS)) %>% 
  write_csv(DEGs.cluster.genotype,
            'neuropeptides/avp.oxt/avp/limmatrend/genotype.sct/DEGs.cluster.csv')

# save hub genes for GO analysis with:
# to get GO term accession use biomart: http://useast.ensembl.org/biomart/martview/

## combine go terms
DEGs.cluster.genotype.GO.terms = read.csv('neuropeptides/avp.oxt/avp/limmatrend/genotype.sct/DEGs.cluster.GO.terms.csv')

# join with pvalue 
DEGs.cluster.genotype.GO.terms.pvalue = DEGs.cluster.genotype.GO.terms %>% 
  left_join(DEGs.cluster.genotype ,
            by = c("Gene.name" = "Gene")) %>% 
  filter(!is.na(adj.P.Val)) %>% 
  rbind(DEGs.cluster.genotype.GO.terms %>% 
          left_join(DEGs.cluster.genotype ,
                    by = c("Gene.stable.ID" = "Gene")) %>% 
          filter(!is.na(adj.P.Val)))  

# save
DEGs.cluster.genotype.GO.terms.pvalue %>% 
  select(c(Gene.stable.ID,
           GO.term.accession,
           Gene.name,
           Keep,
           adj.P.Val)) %>% 
  distinct() %>% 
  relocate(GO.term.accession, 
           .after = last_col())%>% 
  relocate(adj.P.Val, 
           .after = last_col()) %>% 
  arrange(Keep) %>% 
  write_csv('neuropeptides/avp.oxt/avp/limmatrend/genotype.sct/DEGs.cluster.GO.terms.pvalue.csv')

# to look for enrichment and create figures use GO term accession and kME (higher is better) for each module: http://revigo.irb.hr/
# remember to sort by kME and filter out blank GO terms before copying to revigo
# all module genes




#### compare integrated vs sct limmatrend ####
### graph scatterplot of logfc for each comparison
## combine data
# select logfc columns
# pivot longer to join by gene
avp.neuron.reduce.genotype.vector.limma.results.df.sctvsintegrated = avp.neuron.reduce.genotype.vector.limma.results.df.sct %>% 
  select(Gene,
         logFC_DvsS,
         logFC_DvsS_0,
         logFC_DvsS_1,
         logFC_DvsS_2,
         logFC_DvsS_3,
         logFC_DvsS_4) %>% 
  pivot_longer(names_to = 'contrast.logFC',
               values_to = 'logFC.sct',
               cols = c(logFC_DvsS,
                        logFC_DvsS_0,
                        logFC_DvsS_1,
                        logFC_DvsS_2,
                        logFC_DvsS_3,
                        logFC_DvsS_4)) %>% 
  full_join(avp.neuron.reduce.genotype.vector.limma.results.df %>% 
              select(Gene,
                     logFC_DvsS,
                     logFC_DvsS_0,
                     logFC_DvsS_1,
                     logFC_DvsS_2,
                     logFC_DvsS_3,
                     logFC_DvsS_4) %>% 
              pivot_longer(names_to = 'contrast.logFC',
                           values_to = 'logFC',
                           cols = c(logFC_DvsS,
                                    logFC_DvsS_0,
                                    logFC_DvsS_1,
                                    logFC_DvsS_2,
                                    logFC_DvsS_3,
                                    logFC_DvsS_4)))

# create pvalue dataframe
tmp.df = avp.neuron.reduce.genotype.vector.limma.results.df.sct %>% 
  select(Gene,
         P.Value_DvsS,
         P.Value_DvsS_0,
         P.Value_DvsS_1,
         P.Value_DvsS_2,
         P.Value_DvsS_3,
         P.Value_DvsS_4) %>% 
  pivot_longer(names_to = 'contrast.P.Value',
               values_to = 'P.Value.sct',
               cols = c(P.Value_DvsS,
                        P.Value_DvsS_0,
                        P.Value_DvsS_1,
                        P.Value_DvsS_2,
                        P.Value_DvsS_3,
                        P.Value_DvsS_4)) %>% 
  full_join(avp.neuron.reduce.genotype.vector.limma.results.df %>% 
              select(Gene,
                     P.Value_DvsS,
                     P.Value_DvsS_0,
                     P.Value_DvsS_1,
                     P.Value_DvsS_2,
                     P.Value_DvsS_3,
                     P.Value_DvsS_4) %>% 
              pivot_longer(names_to = 'contrast.P.Value',
                           values_to = 'P.Value',
                           cols = c(P.Value_DvsS,
                                    P.Value_DvsS_0,
                                    P.Value_DvsS_1,
                                    P.Value_DvsS_2,
                                    P.Value_DvsS_3,
                                    P.Value_DvsS_4)))

# combine pvalue with logfc data
avp.neuron.reduce.genotype.vector.limma.results.df.sctvsintegrated = avp.neuron.reduce.genotype.vector.limma.results.df.sctvsintegrated %>% 
  cbind(tmp.df%>% 
          select(-c(Gene))) 

# remove temp dataframe
rm(tmp.df)

# create significance and direction
avp.neuron.reduce.genotype.vector.limma.results.df.sctvsintegrated = avp.neuron.reduce.genotype.vector.limma.results.df.sctvsintegrated %>% 
  mutate(direction.sct = ifelse(logFC.sct > 0,
                                1,
                                -1),
         direction = ifelse(logFC > 0,
                                1,
                                -1),
         sig = case_when(P.Value <= 0.05 | P.Value.sct <= 0.05 ~ '0.05',
                         P.Value <= 0.1 | P.Value.sct <= 0.1 ~ '0.1',
                         TRUE ~ 'not.sig'))

### graph
## create scatterplot
avp.neuron.reduce.genotype.vector.limma.results.df.sctvsintegrated %>% 
  ggplot(aes(x = logFC,
             y= logFC.sct,
             color = sig)) +
  geom_point() +
  geom_abline(slope = 1,
              intercept = 0) +
  facet_wrap(.~contrast.logFC) +
  theme_bw() +
  xlim(-max(abs(avp.neuron.reduce.genotype.vector.limma.results.df.sctvsintegrated$logFC)), 
       max(abs(avp.neuron.reduce.genotype.vector.limma.results.df.sctvsintegrated$logFC)))+
  ylim(-max(abs(avp.neuron.reduce.genotype.vector.limma.results.df.sctvsintegrated$logFC.sct)), 
       max(abs(avp.neuron.reduce.genotype.vector.limma.results.df.sctvsintegrated$logFC.sct)))  +
  scale_color_manual(values = c('red',
                                'pink',
                                'black'))
ggsave('./neuropeptides/avp.oxt/avp/limmatrend/genotype.sct/comparison.logfc.png',
       width = 10,
       height = 10)

## create scatterplot
avp.neuron.reduce.genotype.vector.limma.results.df.sctvsintegrated %>% 
  mutate(P.Value = P.Value * direction,
         P.Value.sct = P.Value.sct * direction.sct) %>% 
  ggplot(aes(x = P.Value,
             y= P.Value.sct,
             color = sig)) +
  geom_point() +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  facet_wrap(.~contrast.P.Value) +
  theme_classic() +
  scale_color_manual(values = c('red',
                                'pink',
                                'black'))
ggsave('./neuropeptides/avp.oxt/avp/limmatrend/genotype.sct/comparison.pvalue.png',
       width = 10,
       height = 10)














#### limmatrend genotype, clusters, and status approach treat ####
### calculate variable genes
## identify top 500 variable genes
burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.variable <- FindVariableFeatures(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster, 
                                                                                         selection.method = "vst", 
                                                                                         nfeatures = 500, 
                                                                                         verbose = F)
# identify top 500 variable genes
avp.neuron.reduce.group.topgenes <- head(VariableFeatures(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.variable), 
                                         500)

## create vector of factor
avp.neuron.reduce.DomvsSub.vector.list = burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.expression %>% 
  mutate(avp.neuron.orig.ident = orig.ident %>% 
           as.factor()) %>% 
  pull(avp.neuron.orig.ident) %>% 
  droplevels()

### counts matrix
## raw read count matrix
## rows = genes, columns = cells
# only keep 500 variable genes
# set negative values to 0
avp.neuron.reduce.DomvsSub.vector.count = GetAssayData(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster) %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column('gene') %>% 
  select(c(gene,
           burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.expression %>% 
             pull(Cell.id))) %>% 
  filter(gene %in% avp.neuron.reduce.group.topgenes) %>% 
  column_to_rownames('gene') %>% 
  as.matrix() %>% 
  pmax(0)


## run with reduce avp neurons
#need list with L with  count and condt
# add in cell cluster with condt.2
## create vector of factor
avp.neuron.reduce.DomvsSub.vector.list.cluster = burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.expression %>% 
  mutate(avp.neuron.cluster = integrated_snn_res.0.2 %>% 
           as.factor()) %>% 
  pull(avp.neuron.cluster) %>% 
  droplevels


# add in genotype
## create vector of factor
avp.neuron.reduce.DomvsSub.vector.list.genotype = burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.expression %>% 
  mutate(genotype = Genotype.id %>% 
           as.factor()) %>% 
  pull(genotype) %>% 
  droplevels()

#create list
avp.neuron.reduce.vector.limma = list(count = avp.neuron.reduce.DomvsSub.vector.count,
                                      condt.2 = avp.neuron.reduce.DomvsSub.vector.list.cluster,
                                      genotype = avp.neuron.reduce.DomvsSub.vector.list.genotype)
# add genotype.ID?



#load libraries
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(edgeR))

#create function
run_limmatrend_complex_genotype <- function(L) {
  message("limmatrend")
  session_info <- sessionInfo()
  timing <- system.time({
    treat <- paste(L$genotype, 
                   L$condt.2,
                   sep = '.')
    design <- model.matrix(~0+treat) 
    contrasts <- makeContrasts(DvsS_0 = (treatDom.A.0 + treatDom.B.0 + treatDom.C.0)/3 - (treatSub.A.0 + treatSub.B.0 + treatSub.C.0)/3, 
                               DvsS_1 = (treatDom.A.1 + treatDom.B.1 + treatDom.C.1)/3 - (treatSub.A.1 + treatSub.B.1 + treatSub.C.1)/3,
                               DvsS_2 = (treatDom.A.2 + treatDom.B.2 + treatDom.C.2)/3 - (treatSub.A.2 + treatSub.B.2 + treatSub.C.2)/3,
                               DvsS_3 = (treatDom.A.3 + treatDom.B.3 + treatDom.C.3)/3 - (treatSub.A.3 + treatSub.B.3 + treatSub.C.3)/3, 
                               DvsS_4 = (treatDom.A.4 + treatDom.B.4 + treatDom.C.4)/3 - (treatSub.A.4 + treatSub.B.4 + treatSub.C.4)/3,
                               DvsS=(treatDom.A.0 + treatDom.B.0 + treatDom.C.0+treatDom.A.1 + treatDom.B.1 + treatDom.C.1+treatDom.A.2 + treatDom.B.2 + treatDom.C.2+treatDom.A.3 + treatDom.B.3 + treatDom.C.3+treatDom.A.4 + treatDom.B.4 + treatDom.C.4)/15-(treatSub.A.0 + treatSub.B.0 + treatSub.C.0+treatSub.A.1 + treatSub.B.1 + treatSub.C.1+treatSub.A.2 + treatSub.B.2 + treatSub.C.2+treatSub.A.3 + treatSub.B.3 + treatSub.C.3+treatSub.A.4 + treatSub.B.4 + treatSub.C.4)/15,
                               levels = design)
    dge <- DGEList(L$count, 
                   group = treat)
    dge <- calcNormFactors(dge)
    
    y <- new("EList")
    y$E <- edgeR::cpm(dge, 
                      log = TRUE, 
                      prior.count = 3)
    fit <- lmFit(y, 
                 design = design)
    fit <- contrasts.fit(fit , 
                         contrasts)
    fit <- eBayes(fit, 
                  trend = TRUE,
                  robust = TRUE)
    tt0 <- topTable(fit, 
                    n = Inf,
                    coef = "DvsS_0",
                    adjust.method = "BH")
    tt1 <- topTable(fit, 
                    n = Inf,
                    coef = "DvsS_1",
                    adjust.method = "BH")
    tt2 <- topTable(fit, 
                    n = Inf,
                    coef = "DvsS_2",
                    adjust.method = "BH")
    tt3 <- topTable(fit, 
                    n = Inf,
                    coef = "DvsS_3",
                    adjust.method = "BH")
    tt4 <- topTable(fit, 
                    n = Inf,
                    coef = "DvsS_4",
                    adjust.method = "BH")
    ttDS <- topTable(fit, 
                     n = Inf,
                     coef = "DvsS",
                     adjust.method = "BH")
  })
  # Open pdf file
  pdf(file= "./neuropeptides/avp.oxt/avp/limmatrend/genotype.treat/limmatrend.histograms.pdf" )
  # create a 2X2 grid
  par( mfrow= c(2,2) )
  #graph
  hist(tt0$P.Value, 50)
  hist(tt0$adj.P.Val, 50)
  hist(tt1$P.Value, 50)
  hist(tt1$adj.P.Val, 50)
  hist(tt2$P.Value, 50)
  hist(tt2$adj.P.Val, 50)
  hist(tt3$P.Value, 50)
  hist(tt3$adj.P.Val, 50)
  hist(tt4$P.Value, 50)
  hist(tt4$adj.P.Val, 50)
  hist(ttDS$P.Value, 50)
  hist(ttDS$adj.P.Val, 50)
  dev.off()
  
  # Open pdf file
  pdf(file= "./neuropeptides/avp.oxt/avp/limmatrend/genotype.treat/limmatrend.MDS.pdf" )
  # create a 2X1 grid
  par( mfrow= c(2,1) )
  limma::plotMDS(dge, 
                 col = as.numeric(as.factor(L$condt)), 
                 pch = 19)
  plotMD(fit)
  limma::plotMDS(dge, 
                 col = as.numeric(as.factor(L$condt.2)), 
                 pch = 19)
  plotMD(fit)
  dev.off()
  
  #print results
  list(session_info = session_info,
       timing = timing,
       tt0 = tt0,
       tt1 = tt1,
       tt2 = tt2,
       tt3 = tt3,
       tt4 = tt4,
       ttDS = ttDS)
}


###run function  
avp.neuron.reduce.genotype.vector.limma.results = run_limmatrend_complex_genotype(avp.neuron.reduce.vector.limma)

# save results to dataframe
avp.neuron.reduce.genotype.vector.limma.results.df.treat = full_join(avp.neuron.reduce.genotype.vector.limma.results$tt0 %>% 
                                                                 rename_with(~paste0(.,"_DvsS_0")) %>% 
                                                                 rownames_to_column("Gene"),
                                                               avp.neuron.reduce.genotype.vector.limma.results$tt1 %>% 
                                                                 rename_with(~paste0(.,"_DvsS_1")) %>% 
                                                                 rownames_to_column("Gene")) %>% 
  full_join(avp.neuron.reduce.genotype.vector.limma.results$tt2 %>% 
              rename_with(~paste0(.,"_DvsS_2")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(avp.neuron.reduce.genotype.vector.limma.results$tt3 %>% 
              rename_with(~paste0(.,"_DvsS_3")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(avp.neuron.reduce.genotype.vector.limma.results$tt4 %>% 
              rename_with(~paste0(.,"_DvsS_4")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(avp.neuron.reduce.genotype.vector.limma.results$ttDS %>% 
              rename_with(~paste0(.,"_DvsS")) %>% 
              rownames_to_column("Gene")) 



#add color for significance  
avp.neuron.reduce.genotype.vector.limma.results.df.treat = avp.neuron.reduce.genotype.vector.limma.results.df.treat %>% 
  mutate(Sig_DvsS_0 = ifelse(adj.P.Val_DvsS_0 <= 0.05,
                             'Sig',
                             'Not Sig'),
         Direction.type_DvsS_0 = ifelse(logFC_DvsS_0 > 0,
                                        'up',
                                        'down'),
         Sig.direction_DvsS_0 = ifelse(Sig_DvsS_0 == 'Sig',
                                       Direction.type_DvsS_0,
                                       Sig_DvsS_0)) %>% 
  mutate(Sig_DvsS_1 = ifelse(adj.P.Val_DvsS_1 <= 0.05,
                             'Sig',
                             'Not Sig'),
         Direction.type_DvsS_1 = ifelse(logFC_DvsS_1 > 0,
                                        'up',
                                        'down'),
         Sig.direction_DvsS_1 = ifelse(Sig_DvsS_1 == 'Sig',
                                       Direction.type_DvsS_1,
                                       Sig_DvsS_1)) %>% 
  mutate(Sig_DvsS_2 = ifelse(adj.P.Val_DvsS_2 <= 0.05,
                             'Sig',
                             'Not Sig'),
         Direction.type_DvsS_2 = ifelse(logFC_DvsS_2 > 0,
                                        'up',
                                        'down'),
         Sig.direction_DvsS_2 = ifelse(Sig_DvsS_2 == 'Sig',
                                       Direction.type_DvsS_2,
                                       Sig_DvsS_2)) %>% 
  mutate(Sig_DvsS_3 = ifelse(adj.P.Val_DvsS_3 <= 0.05,
                             'Sig',
                             'Not Sig'),
         Direction.type_DvsS_3 = ifelse(logFC_DvsS_3 > 0,
                                        'up',
                                        'down'),
         Sig.direction_DvsS_3 = ifelse(Sig_DvsS_3 == 'Sig',
                                       Direction.type_DvsS_3,
                                       Sig_DvsS_3)) %>% 
  mutate(Sig_DvsS_4 = ifelse(adj.P.Val_DvsS_4 <= 0.05,
                             'Sig',
                             'Not Sig'),
         Direction.type_DvsS_4 = ifelse(logFC_DvsS_4 > 0,
                                        'up',
                                        'down'),
         Sig.direction_DvsS_4 = ifelse(Sig_DvsS_4 == 'Sig',
                                       Direction.type_DvsS_4,
                                       Sig_DvsS_4)) %>% 
  mutate(Sig_DvsS = ifelse(adj.P.Val_DvsS <= 0.05,
                           'Sig',
                           'Not Sig'),
         Direction.type_DvsS = ifelse(logFC_DvsS > 0,
                                      'up',
                                      'down'),
         Sig.direction_DvsS = ifelse(Sig_DvsS == 'Sig',
                                     Direction.type_DvsS,
                                     Sig_DvsS))


##graph volcano plot
# create comparison list
limma.genotpye.vector = c("DvsS_0",
                          "DvsS_1",
                          "DvsS_2",
                          "DvsS_3",
                          "DvsS_4",
                          "DvsS")

for (i in limma.genotpye.vector) {
  # graph volcano plot
  avp.neuron.reduce.genotype.vector.limma.results.df.treat %>% 
    mutate(sig.label = ifelse(get(paste0("Sig_", i)) == 'Sig',
                              Gene,
                              '')) %>% 
    ggplot(aes(x = get(paste0("logFC_", i)),
               y = -log10(get(paste0("adj.P.Val_", i))),
               color = get(paste0("Sig.direction_", i))))+
    geom_hline(yintercept = -log10(0.05),
               linetype = 'dotted') +
    geom_point(size = 5) +
    theme_classic() + 
    scale_color_manual(values=c("21B9CA", 
                                "grey", 
                                "orange3")) +
    geom_text(aes(label = sig.label),
              vjust = 0, 
              nudge_y = 0.10,
              size = 5) +
    theme(text = element_text(size = 20),
          legend.position = 'none') +
    xlab(paste0("logFC_", i)) +
    ylab( paste0("-log10(adj.P.Val_", i,")")) +
    ggtitle(paste0(i, " volcano plot"))
  ggsave(paste0('neuropeptides/avp.oxt/avp/limmatrend/genotype.treat/limma.avp.neurons.reduce.volcano.', i, '.png'),
         width = 5,
         height = 5)
}





#### limmatrend genotype, clusters, and status approach with SCT treat ####
#### Variable have same name as 'limmatrend genotype, clusters, and status approach'
### create new expression matrix with SCT
##expression level
burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.expression.sct = full_join(full_join(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster@reductions$umap@cell.embeddings %>% 
                                                                                               as.data.frame() %>% 
                                                                                               rownames_to_column("Cell.id"),
                                                                                             burtoni.snseq.combined.sct.reduce.avp.neurons.recluster@meta.data %>% 
                                                                                               rownames_to_column("Cell.id")),
                                                                                   burtoni.snseq.combined.sct.reduce.avp.neurons.recluster@assays$SCT@scale.data %>% 
                                                                                     as.data.frame() %>% 
                                                                                     filter(rownames(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster@assays$SCT@scale.data) %in% c('avp',
                                                                                                                                                                                           'oxt')) %>% 
                                                                                     t() %>% as.data.frame() %>% 
                                                                                     rownames_to_column('Cell.id')) %>% 
  full_join(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster@reductions$pca@cell.embeddings %>% 
              as.data.frame() %>% 
              select(PC_1,
                     PC_2,
                     PC_3,
                     PC_4,
                     PC_5) %>% 
              rownames_to_column("Cell.id")) %>% 
  mutate(Oxt.cell = ifelse(oxt > 2 & avp > 2,
                           'both',
                           ifelse(oxt > 2,
                                  'oxt',
                                  'avp')))

### calculate variable genes
## identify top 500 variable genes
burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.variable <- FindVariableFeatures(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster, 
                                                                                         selection.method = "vst", 
                                                                                         nfeatures = 500, 
                                                                                         verbose = F)
# identify top 500 variable genes
avp.neuron.reduce.group.topgenes <- head(VariableFeatures(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.variable), 
                                         500)

## create vector of factor
avp.neuron.reduce.DomvsSub.vector.list.sct = burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.expression.sct %>% 
  mutate(avp.neuron.orig.ident = orig.ident %>% 
           as.factor()) %>% 
  pull(avp.neuron.orig.ident) %>% 
  droplevels()

### counts matrix
## raw read count matrix
## rows = genes, columns = cells
# only keep 500 variable genes
# set negative values to 0
avp.neuron.reduce.DomvsSub.vector.count.sct = GetAssayData(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster,
                                                           assay = 'SCT') %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column('gene') %>% 
  select(c(gene,
           burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.expression.sct %>% 
             pull(Cell.id))) %>% 
  filter(gene %in% avp.neuron.reduce.group.topgenes) %>% 
  column_to_rownames('gene') %>% 
  as.matrix() %>% 
  pmax(0)


## run with reduce avp neurons
#need list with L with  count and condt
# add in cell cluster with condt.2
## create vector of factor
avp.neuron.reduce.DomvsSub.vector.list.cluster.sct = burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.expression.sct %>% 
  mutate(avp.neuron.cluster = integrated_snn_res.0.2 %>% 
           as.factor()) %>% 
  pull(avp.neuron.cluster) %>% 
  droplevels


# add in genotype
## create vector of factor
avp.neuron.reduce.DomvsSub.vector.list.genotype.sct = burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.expression.sct %>% 
  mutate(genotype = Genotype.id %>% 
           as.factor()) %>% 
  pull(genotype) %>% 
  droplevels()

#create list
avp.neuron.reduce.vector.limma.sct = list(count = avp.neuron.reduce.DomvsSub.vector.count.sct,
                                          condt.2 = avp.neuron.reduce.DomvsSub.vector.list.cluster.sct,
                                          genotype = avp.neuron.reduce.DomvsSub.vector.list.genotype.sct)
# add genotype.ID?



#load libraries
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(edgeR))

#create function
run_limmatrend_complex_genotype.sct <- function(L) {
  message("limmatrend")
  session_info <- sessionInfo()
  timing <- system.time({
    treat <- paste(L$genotype, 
                   L$condt.2,
                   sep = '.')
    design <- model.matrix(~0+treat) 
    contrasts <- makeContrasts(DvsS_0 = (treatDom.A.0 + treatDom.B.0 + treatDom.C.0)/3 - (treatSub.A.0 + treatSub.B.0 + treatSub.C.0)/3, 
                               DvsS_1 = (treatDom.A.1 + treatDom.B.1 + treatDom.C.1)/3 - (treatSub.A.1 + treatSub.B.1 + treatSub.C.1)/3,
                               DvsS_2 = (treatDom.A.2 + treatDom.B.2 + treatDom.C.2)/3 - (treatSub.A.2 + treatSub.B.2 + treatSub.C.2)/3,
                               DvsS_3 = (treatDom.A.3 + treatDom.B.3 + treatDom.C.3)/3 - (treatSub.A.3 + treatSub.B.3 + treatSub.C.3)/3, 
                               DvsS_4 = (treatDom.A.4 + treatDom.B.4 + treatDom.C.4)/3 - (treatSub.A.4 + treatSub.B.4 + treatSub.C.4)/3,
                               DvsS=(treatDom.A.0 + treatDom.B.0 + treatDom.C.0+treatDom.A.1 + treatDom.B.1 + treatDom.C.1+treatDom.A.2 + treatDom.B.2 + treatDom.C.2+treatDom.A.3 + treatDom.B.3 + treatDom.C.3+treatDom.A.4 + treatDom.B.4 + treatDom.C.4)/15-(treatSub.A.0 + treatSub.B.0 + treatSub.C.0+treatSub.A.1 + treatSub.B.1 + treatSub.C.1+treatSub.A.2 + treatSub.B.2 + treatSub.C.2+treatSub.A.3 + treatSub.B.3 + treatSub.C.3+treatSub.A.4 + treatSub.B.4 + treatSub.C.4)/15,
                               levels = design)
    dge <- DGEList(L$count, 
                   group = treat)
    dge <- calcNormFactors(dge)
    
    y <- new("EList")
    y$E <- edgeR::cpm(dge, 
                      log = TRUE, 
                      prior.count = 3)
    fit <- lmFit(y, 
                 design = design)
    fit <- contrasts.fit(fit , 
                         contrasts)
    fit <- eBayes(fit, 
                  trend = TRUE,
                  robust = TRUE)
    tt0 <- topTable(fit, 
                    n = Inf,
                    coef = "DvsS_0",
                    adjust.method = "BH")
    tt1 <- topTable(fit, 
                    n = Inf,
                    coef = "DvsS_1",
                    adjust.method = "BH")
    tt2 <- topTable(fit, 
                    n = Inf,
                    coef = "DvsS_2",
                    adjust.method = "BH")
    tt3 <- topTable(fit, 
                    n = Inf,
                    coef = "DvsS_3",
                    adjust.method = "BH")
    tt4 <- topTable(fit, 
                    n = Inf,
                    coef = "DvsS_4",
                    adjust.method = "BH")
    ttDS <- topTable(fit, 
                     n = Inf,
                     coef = "DvsS",
                     adjust.method = "BH")
  })
  # Open pdf file
  pdf(file= "./neuropeptides/avp.oxt/avp/limmatrend/genotype.sct.treat/limmatrend.histograms.pdf" )
  # create a 2X2 grid
  par( mfrow= c(2,2) )
  #graph
  hist(tt0$P.Value, 50)
  hist(tt0$adj.P.Val, 50)
  hist(tt1$P.Value, 50)
  hist(tt1$adj.P.Val, 50)
  hist(tt2$P.Value, 50)
  hist(tt2$adj.P.Val, 50)
  hist(tt3$P.Value, 50)
  hist(tt3$adj.P.Val, 50)
  hist(tt4$P.Value, 50)
  hist(tt4$adj.P.Val, 50)
  hist(ttDS$P.Value, 50)
  hist(ttDS$adj.P.Val, 50)
  dev.off()
  
  # Open pdf file
  pdf(file= "./neuropeptides/avp.oxt/avp/limmatrend/genotype.sct.treat/limmatrend.MDS.pdf" )
  # create a 2X1 grid
  par( mfrow= c(2,1) )
  limma::plotMDS(dge, 
                 col = as.numeric(as.factor(L$condt)), 
                 pch = 19)
  plotMD(fit)
  limma::plotMDS(dge, 
                 col = as.numeric(as.factor(L$condt.2)), 
                 pch = 19)
  plotMD(fit)
  dev.off()
  
  #print results
  list(session_info = session_info,
       timing = timing,
       tt0 = tt0,
       tt1 = tt1,
       tt2 = tt2,
       tt3 = tt3,
       tt4 = tt4,
       ttDS = ttDS)
}


###run function  
avp.neuron.reduce.genotype.vector.limma.results.sct = run_limmatrend_complex_genotype.sct(avp.neuron.reduce.vector.limma.sct)

# save results to dataframe
avp.neuron.reduce.genotype.vector.limma.results.df.sct.treat = full_join(avp.neuron.reduce.genotype.vector.limma.results.sct$tt0 %>% 
                                                                           rename_with(~paste0(.,"_DvsS_0")) %>% 
                                                                           rownames_to_column("Gene"),
                                                                         avp.neuron.reduce.genotype.vector.limma.results.sct$tt1 %>% 
                                                                           rename_with(~paste0(.,"_DvsS_1")) %>% 
                                                                           rownames_to_column("Gene")) %>% 
  full_join(avp.neuron.reduce.genotype.vector.limma.results.sct$tt2 %>% 
              rename_with(~paste0(.,"_DvsS_2")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(avp.neuron.reduce.genotype.vector.limma.results.sct$tt3 %>% 
              rename_with(~paste0(.,"_DvsS_3")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(avp.neuron.reduce.genotype.vector.limma.results.sct$tt4 %>% 
              rename_with(~paste0(.,"_DvsS_4")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(avp.neuron.reduce.genotype.vector.limma.results.sct$ttDS %>% 
              rename_with(~paste0(.,"_DvsS")) %>% 
              rownames_to_column("Gene")) 



#add color for significance  
avp.neuron.reduce.genotype.vector.limma.results.df.sct.treat = avp.neuron.reduce.genotype.vector.limma.results.df.sct.treat %>% 
  mutate(Sig_DvsS_0 = ifelse(adj.P.Val_DvsS_0 <= 0.05,
                             'Sig',
                             'Not Sig'),
         Direction.type_DvsS_0 = ifelse(logFC_DvsS_0 > 0,
                                        'up',
                                        'down'),
         Sig.direction_DvsS_0 = ifelse(Sig_DvsS_0 == 'Sig',
                                       Direction.type_DvsS_0,
                                       Sig_DvsS_0)) %>% 
  mutate(Sig_DvsS_1 = ifelse(adj.P.Val_DvsS_1 <= 0.05,
                             'Sig',
                             'Not Sig'),
         Direction.type_DvsS_1 = ifelse(logFC_DvsS_1 > 0,
                                        'up',
                                        'down'),
         Sig.direction_DvsS_1 = ifelse(Sig_DvsS_1 == 'Sig',
                                       Direction.type_DvsS_1,
                                       Sig_DvsS_1)) %>% 
  mutate(Sig_DvsS_2 = ifelse(adj.P.Val_DvsS_2 <= 0.05,
                             'Sig',
                             'Not Sig'),
         Direction.type_DvsS_2 = ifelse(logFC_DvsS_2 > 0,
                                        'up',
                                        'down'),
         Sig.direction_DvsS_2 = ifelse(Sig_DvsS_2 == 'Sig',
                                       Direction.type_DvsS_2,
                                       Sig_DvsS_2)) %>% 
  mutate(Sig_DvsS_3 = ifelse(adj.P.Val_DvsS_3 <= 0.05,
                             'Sig',
                             'Not Sig'),
         Direction.type_DvsS_3 = ifelse(logFC_DvsS_3 > 0,
                                        'up',
                                        'down'),
         Sig.direction_DvsS_3 = ifelse(Sig_DvsS_3 == 'Sig',
                                       Direction.type_DvsS_3,
                                       Sig_DvsS_3)) %>% 
  mutate(Sig_DvsS_4 = ifelse(adj.P.Val_DvsS_4 <= 0.05,
                             'Sig',
                             'Not Sig'),
         Direction.type_DvsS_4 = ifelse(logFC_DvsS_4 > 0,
                                        'up',
                                        'down'),
         Sig.direction_DvsS_4 = ifelse(Sig_DvsS_4 == 'Sig',
                                       Direction.type_DvsS_4,
                                       Sig_DvsS_4)) %>% 
  mutate(Sig_DvsS = ifelse(adj.P.Val_DvsS <= 0.05,
                           'Sig',
                           'Not Sig'),
         Direction.type_DvsS = ifelse(logFC_DvsS > 0,
                                      'up',
                                      'down'),
         Sig.direction_DvsS = ifelse(Sig_DvsS == 'Sig',
                                     Direction.type_DvsS,
                                     Sig_DvsS))


##graph volcano plot
# create comparison list
limma.genotpye.vector = c("DvsS_0",
                          "DvsS_1",
                          "DvsS_2",
                          "DvsS_3",
                          "DvsS_4",
                          "DvsS")

for (i in limma.genotpye.vector) {
  # graph volcano plot
  avp.neuron.reduce.genotype.vector.limma.results.df.sct.treat %>% 
    mutate(sig.label = ifelse(get(paste0("Sig_", i)) == 'Sig',
                              Gene,
                              '')) %>% 
    ggplot(aes(x = get(paste0("logFC_", i)),
               y = -log10(get(paste0("adj.P.Val_", i))),
               color = get(paste0("Sig.direction_", i))))+
    geom_hline(yintercept = -log10(0.05),
               linetype = 'dotted') +
    geom_point(size = 5) +
    theme_classic() + 
    scale_color_manual(values=c("21B9CA", 
                                "grey", 
                                "orange3")) +
    geom_text(aes(label = sig.label),
              vjust = 0, 
              nudge_y = 0.10,
              size = 5) +
    theme(text = element_text(size = 20),
          legend.position = 'none') +
    xlab(paste0("logFC_", i)) +
    ylab( paste0("-log10(adj.P.Val_", i,")")) +
    ggtitle(paste0(i, " volcano plot"))
  ggsave(paste0('neuropeptides/avp.oxt/avp/limmatrend/genotype.sct.treat/limma.avp.neurons.reduce.volcano.', i, '.png'),
         width = 5,
         height = 5)
}




#### compare integrated limmatrend models ####
### graph scatterplot of logfc for each comparison
## combine data
# select logfc columns
# pivot longer to join by gene
avp.neuron.reduce.genotype.vector.limma.results.df.integrated.models = avp.neuron.reduce.genotype.vector.limma.results.df.treat %>% 
  select(Gene,
         logFC_DvsS,
         logFC_DvsS_0,
         logFC_DvsS_1,
         logFC_DvsS_2,
         logFC_DvsS_3,
         logFC_DvsS_4) %>% 
  pivot_longer(names_to = 'contrast.logFC',
               values_to = 'logFC.sct',
               cols = c(logFC_DvsS,
                        logFC_DvsS_0,
                        logFC_DvsS_1,
                        logFC_DvsS_2,
                        logFC_DvsS_3,
                        logFC_DvsS_4)) %>% 
  full_join(avp.neuron.reduce.genotype.vector.limma.results.df %>% 
              select(Gene,
                     logFC_DvsS,
                     logFC_DvsS_0,
                     logFC_DvsS_1,
                     logFC_DvsS_2,
                     logFC_DvsS_3,
                     logFC_DvsS_4) %>% 
              pivot_longer(names_to = 'contrast.logFC',
                           values_to = 'logFC',
                           cols = c(logFC_DvsS,
                                    logFC_DvsS_0,
                                    logFC_DvsS_1,
                                    logFC_DvsS_2,
                                    logFC_DvsS_3,
                                    logFC_DvsS_4)))

# create pvalue dataframe
tmp.df = avp.neuron.reduce.genotype.vector.limma.results.df.treat %>% 
  select(Gene,
         P.Value_DvsS,
         P.Value_DvsS_0,
         P.Value_DvsS_1,
         P.Value_DvsS_2,
         P.Value_DvsS_3,
         P.Value_DvsS_4) %>% 
  pivot_longer(names_to = 'contrast.P.Value',
               values_to = 'P.Value.sct',
               cols = c(P.Value_DvsS,
                        P.Value_DvsS_0,
                        P.Value_DvsS_1,
                        P.Value_DvsS_2,
                        P.Value_DvsS_3,
                        P.Value_DvsS_4)) %>% 
  full_join(avp.neuron.reduce.genotype.vector.limma.results.df %>% 
              select(Gene,
                     P.Value_DvsS,
                     P.Value_DvsS_0,
                     P.Value_DvsS_1,
                     P.Value_DvsS_2,
                     P.Value_DvsS_3,
                     P.Value_DvsS_4) %>% 
              pivot_longer(names_to = 'contrast.P.Value',
                           values_to = 'P.Value',
                           cols = c(P.Value_DvsS,
                                    P.Value_DvsS_0,
                                    P.Value_DvsS_1,
                                    P.Value_DvsS_2,
                                    P.Value_DvsS_3,
                                    P.Value_DvsS_4)))

# combine pvalue with logfc data
avp.neuron.reduce.genotype.vector.limma.results.df.integrated.models = avp.neuron.reduce.genotype.vector.limma.results.df.integrated.models %>% 
  cbind(tmp.df%>% 
          select(-c(Gene))) 

# remove temp dataframe
rm(tmp.df)

# create significance and direction
avp.neuron.reduce.genotype.vector.limma.results.df.integrated.models = avp.neuron.reduce.genotype.vector.limma.results.df.integrated.models %>% 
  mutate(direction.sct = ifelse(logFC.sct > 0,
                                1,
                                -1),
         direction = ifelse(logFC > 0,
                            1,
                            -1),
         sig = case_when(P.Value <= 0.05 | P.Value.sct <= 0.05 ~ '0.05',
                         P.Value <= 0.1 | P.Value.sct <= 0.1 ~ '0.1',
                         TRUE ~ 'not.sig'))

### graph
## create scatterplot
avp.neuron.reduce.genotype.vector.limma.results.df.integrated.models %>% 
  ggplot(aes(x = logFC,
             y= logFC.sct,
             color = sig)) +
  geom_point() +
  geom_abline(slope = 1,
              intercept = 0) +
  facet_wrap(.~contrast.logFC) +
  theme_bw() +
  xlim(-max(abs(avp.neuron.reduce.genotype.vector.limma.results.df.integrated.models$logFC)), 
       max(abs(avp.neuron.reduce.genotype.vector.limma.results.df.integrated.models$logFC)))+
  ylim(-max(abs(avp.neuron.reduce.genotype.vector.limma.results.df.integrated.models$logFC.sct)), 
       max(abs(avp.neuron.reduce.genotype.vector.limma.results.df.integrated.models$logFC.sct)))  +
  scale_color_manual(values = c('red',
                                'pink',
                                'black'))
ggsave('./neuropeptides/avp.oxt/avp/limmatrend/genotype.treat/comparison.logfc.png',
       width = 10,
       height = 10)

## create scatterplot
avp.neuron.reduce.genotype.vector.limma.results.df.integrated.models %>% 
  mutate(P.Value = P.Value * direction,
         P.Value.sct = P.Value.sct * direction.sct) %>% 
  ggplot(aes(x = P.Value,
             y= P.Value.sct,
             color = sig)) +
  geom_point() +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  facet_wrap(.~contrast.P.Value) +
  theme_classic() +
  scale_color_manual(values = c('red',
                                'pink',
                                'black'))
ggsave('./neuropeptides/avp.oxt/avp/limmatrend/genotype.treat/comparison.pvalue.png',
       width = 10,
       height = 10)














#### compare integrated vs sct limmatrend treat ####
### graph scatterplot of logfc for each comparison
## combine data
# select logfc columns
# pivot longer to join by gene
avp.neuron.reduce.genotype.vector.limma.results.df.sctvsintegrated.treat = avp.neuron.reduce.genotype.vector.limma.results.df.sct.treat %>% 
  select(Gene,
         logFC_DvsS,
         logFC_DvsS_0,
         logFC_DvsS_1,
         logFC_DvsS_2,
         logFC_DvsS_3,
         logFC_DvsS_4) %>% 
  pivot_longer(names_to = 'contrast.logFC',
               values_to = 'logFC.sct',
               cols = c(logFC_DvsS,
                        logFC_DvsS_0,
                        logFC_DvsS_1,
                        logFC_DvsS_2,
                        logFC_DvsS_3,
                        logFC_DvsS_4)) %>% 
  full_join(avp.neuron.reduce.genotype.vector.limma.results.df.treat %>% 
              select(Gene,
                     logFC_DvsS,
                     logFC_DvsS_0,
                     logFC_DvsS_1,
                     logFC_DvsS_2,
                     logFC_DvsS_3,
                     logFC_DvsS_4) %>% 
              pivot_longer(names_to = 'contrast.logFC',
                           values_to = 'logFC',
                           cols = c(logFC_DvsS,
                                    logFC_DvsS_0,
                                    logFC_DvsS_1,
                                    logFC_DvsS_2,
                                    logFC_DvsS_3,
                                    logFC_DvsS_4)))

# create pvalue dataframe
tmp.df = avp.neuron.reduce.genotype.vector.limma.results.df.sct.treat %>% 
  select(Gene,
         P.Value_DvsS,
         P.Value_DvsS_0,
         P.Value_DvsS_1,
         P.Value_DvsS_2,
         P.Value_DvsS_3,
         P.Value_DvsS_4) %>% 
  pivot_longer(names_to = 'contrast.P.Value',
               values_to = 'P.Value.sct',
               cols = c(P.Value_DvsS,
                        P.Value_DvsS_0,
                        P.Value_DvsS_1,
                        P.Value_DvsS_2,
                        P.Value_DvsS_3,
                        P.Value_DvsS_4)) %>% 
  full_join(avp.neuron.reduce.genotype.vector.limma.results.df.treat %>% 
              select(Gene,
                     P.Value_DvsS,
                     P.Value_DvsS_0,
                     P.Value_DvsS_1,
                     P.Value_DvsS_2,
                     P.Value_DvsS_3,
                     P.Value_DvsS_4) %>% 
              pivot_longer(names_to = 'contrast.P.Value',
                           values_to = 'P.Value',
                           cols = c(P.Value_DvsS,
                                    P.Value_DvsS_0,
                                    P.Value_DvsS_1,
                                    P.Value_DvsS_2,
                                    P.Value_DvsS_3,
                                    P.Value_DvsS_4)))

# combine pvalue with logfc data
avp.neuron.reduce.genotype.vector.limma.results.df.sctvsintegrated.treat = avp.neuron.reduce.genotype.vector.limma.results.df.sctvsintegrated.treat %>% 
  cbind(tmp.df%>% 
          select(-c(Gene))) 

# remove temp dataframe
rm(tmp.df)

# create significance and direction
avp.neuron.reduce.genotype.vector.limma.results.df.sctvsintegrated.treat = avp.neuron.reduce.genotype.vector.limma.results.df.sctvsintegrated.treat %>% 
  mutate(direction.sct = ifelse(logFC.sct > 0,
                                1,
                                -1),
         direction = ifelse(logFC > 0,
                            1,
                            -1),
         sig = case_when(P.Value <= 0.05 | P.Value.sct <= 0.05 ~ '0.05',
                         P.Value <= 0.1 | P.Value.sct <= 0.1 ~ '0.1',
                         TRUE ~ 'not.sig'))

### graph
## create scatterplot
avp.neuron.reduce.genotype.vector.limma.results.df.sctvsintegrated.treat %>% 
  ggplot(aes(x = logFC,
             y= logFC.sct,
             color = sig)) +
  geom_point() +
  geom_abline(slope = 1,
              intercept = 0) +
  facet_wrap(.~contrast.logFC) +
  theme_bw() +
  xlim(-max(abs(avp.neuron.reduce.genotype.vector.limma.results.df.sctvsintegrated.treat$logFC)), 
       max(abs(avp.neuron.reduce.genotype.vector.limma.results.df.sctvsintegrated.treat$logFC)))+
  ylim(-max(abs(avp.neuron.reduce.genotype.vector.limma.results.df.sctvsintegrated.treat$logFC.sct)), 
       max(abs(avp.neuron.reduce.genotype.vector.limma.results.df.sctvsintegrated.treat$logFC.sct)))  +
  scale_color_manual(values = c('red',
                                'pink',
                                'black'))
ggsave('./neuropeptides/avp.oxt/avp/limmatrend/genotype.sct.treat/comparison.logfc.png',
       width = 10,
       height = 10)

## create scatterplot
avp.neuron.reduce.genotype.vector.limma.results.df.sctvsintegrated.treat %>% 
  mutate(P.Value = P.Value * direction,
         P.Value.sct = P.Value.sct * direction.sct) %>% 
  ggplot(aes(x = P.Value,
             y= P.Value.sct,
             color = sig)) +
  geom_point() +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  facet_wrap(.~contrast.P.Value) +
  theme_classic() +
  scale_color_manual(values = c('red',
                                'pink',
                                'black'))
ggsave('./neuropeptides/avp.oxt/avp/limmatrend/genotype.sct.treat/comparison.pvalue.png',
       width = 10,
       height = 10)














#### limmatrend genotype, clusters, and status approach with SCT treat prep ####
#### Variable have same name as 'limmatrend genotype, clusters, and status approach'
### calculate variable genes
## identify top 500 variable genes
burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.variable.prep <- FindVariableFeatures(burtoni.snseq.combined.sct.all.avp.neurons.recluster, 
                                                                                         selection.method = "vst", 
                                                                                         nfeatures = 500, 
                                                                                         verbose = F)
# identify top 500 variable genes
avp.neuron.reduce.group.topgenes.prep <- head(VariableFeatures(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.variable.prep), 
                                         500)

# create dummy
burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.expression.sct.prep = burtoni.snseq.combined.sct.all.avp.neurons.recluster.expression

## create vector of factor
avp.neuron.reduce.DomvsSub.vector.list.sct.prep = burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.expression.sct.prep %>% 
  mutate(avp.neuron.orig.ident = orig.ident %>% 
           as.factor()) %>% 
  pull(avp.neuron.orig.ident) %>% 
  droplevels()

# create dummy
burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.prep = burtoni.snseq.combined.sct.all.avp.neurons.recluster

### counts matrix
## raw read count matrix
## rows = genes, columns = cells
# only keep 500 variable genes
# set negative values to 0
avp.neuron.reduce.DomvsSub.vector.count.sct.prep = GetAssayData(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.prep,
                                                           assay = 'SCT') %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column('gene') %>% 
  select(c(gene,
           burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.expression.sct.prep %>% 
             pull(Cell.id))) %>% 
  filter(gene %in% avp.neuron.reduce.group.topgenes.prep) %>% 
  column_to_rownames('gene') %>% 
  as.matrix() %>% 
  pmax(0)


## run with reduce avp neurons
#need list with L with  count and condt
# add in cell cluster with condt.2
## create vector of factor
avp.neuron.reduce.DomvsSub.vector.list.cluster.sct.prep = burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.expression.sct.prep %>% 
  mutate(avp.neuron.cluster = integrated_snn_res.0.8 %>% 
           as.factor()) %>% 
  pull(avp.neuron.cluster) %>% 
  droplevels


# add in genotype
## create vector of factor
avp.neuron.reduce.DomvsSub.vector.list.genotype.sct.prep = burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.expression.sct.prep %>% 
  mutate(genotype = Genotype.id %>% 
           as.factor()) %>% 
  pull(genotype) %>% 
  droplevels()

#create list
avp.neuron.reduce.vector.limma.sct.prep = list(count = avp.neuron.reduce.DomvsSub.vector.count.sct.prep,
                                          condt.2 = avp.neuron.reduce.DomvsSub.vector.list.cluster.sct.prep,
                                          genotype = avp.neuron.reduce.DomvsSub.vector.list.genotype.sct.prep)
# add genotype.ID?



#load libraries
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(edgeR))

#create function
run_limmatrend_complex_genotype.sct <- function(L) {
  message("limmatrend")
  session_info <- sessionInfo()
  timing <- system.time({
    treat <- paste(L$genotype, 
                   L$condt.2,
                   sep = '.')
    design <- model.matrix(~0+treat) 
    contrasts <- makeContrasts(DvsS_0 = (treatDom.A.0 + treatDom.B.0 + treatDom.C.0)/3 - (treatSub.A.0 + treatSub.B.0 + treatSub.C.0)/3, 
                               DvsS_1 = (treatDom.A.1 + treatDom.B.1 + treatDom.C.1)/3 - (treatSub.A.1 + treatSub.B.1 + treatSub.C.1)/3,
                               DvsS_2 = (treatDom.A.2 + treatDom.B.2 + treatDom.C.2)/3 - (treatSub.A.2 + treatSub.B.2 + treatSub.C.2)/3,
                               DvsS_3 = (treatDom.A.3 + treatDom.B.3 + treatDom.C.3)/3 - (treatSub.A.3 + treatSub.B.3 + treatSub.C.3)/3, 
                               DvsS_4 = (treatDom.A.4 + treatDom.B.4 + treatDom.C.4)/3 - (treatSub.A.4 + treatSub.B.4 + treatSub.C.4)/3,
                               DvsS=(treatDom.A.0 + treatDom.B.0 + treatDom.C.0+treatDom.A.1 + treatDom.B.1 + treatDom.C.1+treatDom.A.2 + treatDom.B.2 + treatDom.C.2+treatDom.A.3 + treatDom.B.3 + treatDom.C.3+treatDom.A.4 + treatDom.B.4 + treatDom.C.4)/15-(treatSub.A.0 + treatSub.B.0 + treatSub.C.0+treatSub.A.1 + treatSub.B.1 + treatSub.C.1+treatSub.A.2 + treatSub.B.2 + treatSub.C.2+treatSub.A.3 + treatSub.B.3 + treatSub.C.3+treatSub.A.4 + treatSub.B.4 + treatSub.C.4)/15,
                               levels = design)
    dge <- DGEList(L$count, 
                   group = treat)
    dge <- calcNormFactors(dge)
    
    y <- new("EList")
    y$E <- edgeR::cpm(dge, 
                      log = TRUE, 
                      prior.count = 3)
    fit <- lmFit(y, 
                 design = design)
    fit <- contrasts.fit(fit , 
                         contrasts)
    fit <- eBayes(fit, 
                  trend = TRUE,
                  robust = TRUE)
    tt0 <- topTable(fit, 
                    n = Inf,
                    coef = "DvsS_0",
                    adjust.method = "BH")
    tt1 <- topTable(fit, 
                    n = Inf,
                    coef = "DvsS_1",
                    adjust.method = "BH")
    tt2 <- topTable(fit, 
                    n = Inf,
                    coef = "DvsS_2",
                    adjust.method = "BH")
    tt3 <- topTable(fit, 
                    n = Inf,
                    coef = "DvsS_3",
                    adjust.method = "BH")
    tt4 <- topTable(fit, 
                    n = Inf,
                    coef = "DvsS_4",
                    adjust.method = "BH")
    ttDS <- topTable(fit, 
                     n = Inf,
                     coef = "DvsS",
                     adjust.method = "BH")
  })
  # Open pdf file
  pdf(file= "./neuropeptides/avp.oxt/avp/limmatrend/genotype.sct.treat.prep/limmatrend.histograms.pdf" )
  # create a 2X2 grid
  par( mfrow= c(2,2) )
  #graph
  hist(tt0$P.Value, 50)
  hist(tt0$adj.P.Val, 50)
  hist(tt1$P.Value, 50)
  hist(tt1$adj.P.Val, 50)
  hist(tt2$P.Value, 50)
  hist(tt2$adj.P.Val, 50)
  hist(tt3$P.Value, 50)
  hist(tt3$adj.P.Val, 50)
  hist(tt4$P.Value, 50)
  hist(tt4$adj.P.Val, 50)
  hist(ttDS$P.Value, 50)
  hist(ttDS$adj.P.Val, 50)
  dev.off()
  
  # Open pdf file
  pdf(file= "./neuropeptides/avp.oxt/avp/limmatrend/genotype.sct.treat.prep/limmatrend.MDS.pdf" )
  # create a 2X1 grid
  par( mfrow= c(2,1) )
  limma::plotMDS(dge, 
                 col = as.numeric(as.factor(L$condt)), 
                 pch = 19)
  plotMD(fit)
  limma::plotMDS(dge, 
                 col = as.numeric(as.factor(L$condt.2)), 
                 pch = 19)
  plotMD(fit)
  dev.off()
  
  #print results
  list(session_info = session_info,
       timing = timing,
       tt0 = tt0,
       tt1 = tt1,
       tt2 = tt2,
       tt3 = tt3,
       tt4 = tt4,
       ttDS = ttDS)
}


###run function  
avp.neuron.reduce.genotype.vector.limma.results.sct = run_limmatrend_complex_genotype.sct(avp.neuron.reduce.vector.limma.sct.prep)

# save results to dataframe
avp.neuron.reduce.genotype.vector.limma.results.df.sct.treat.prep = full_join(avp.neuron.reduce.genotype.vector.limma.results.sct$tt0 %>% 
                                                                           rename_with(~paste0(.,"_DvsS_0")) %>% 
                                                                           rownames_to_column("Gene"),
                                                                         avp.neuron.reduce.genotype.vector.limma.results.sct$tt1 %>% 
                                                                           rename_with(~paste0(.,"_DvsS_1")) %>% 
                                                                           rownames_to_column("Gene")) %>% 
  full_join(avp.neuron.reduce.genotype.vector.limma.results.sct$tt2 %>% 
              rename_with(~paste0(.,"_DvsS_2")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(avp.neuron.reduce.genotype.vector.limma.results.sct$tt3 %>% 
              rename_with(~paste0(.,"_DvsS_3")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(avp.neuron.reduce.genotype.vector.limma.results.sct$tt4 %>% 
              rename_with(~paste0(.,"_DvsS_4")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(avp.neuron.reduce.genotype.vector.limma.results.sct$ttDS %>% 
              rename_with(~paste0(.,"_DvsS")) %>% 
              rownames_to_column("Gene")) 



#add color for significance  
avp.neuron.reduce.genotype.vector.limma.results.df.sct.treat.prep = avp.neuron.reduce.genotype.vector.limma.results.df.sct.treat.prep %>% 
  mutate(Sig_DvsS_0 = ifelse(adj.P.Val_DvsS_0 <= 0.05,
                             'Sig',
                             'Not Sig'),
         Direction.type_DvsS_0 = ifelse(logFC_DvsS_0 > 0,
                                        'up',
                                        'down'),
         Sig.direction_DvsS_0 = ifelse(Sig_DvsS_0 == 'Sig',
                                       Direction.type_DvsS_0,
                                       Sig_DvsS_0)) %>% 
  mutate(Sig_DvsS_1 = ifelse(adj.P.Val_DvsS_1 <= 0.05,
                             'Sig',
                             'Not Sig'),
         Direction.type_DvsS_1 = ifelse(logFC_DvsS_1 > 0,
                                        'up',
                                        'down'),
         Sig.direction_DvsS_1 = ifelse(Sig_DvsS_1 == 'Sig',
                                       Direction.type_DvsS_1,
                                       Sig_DvsS_1)) %>% 
  mutate(Sig_DvsS_2 = ifelse(adj.P.Val_DvsS_2 <= 0.05,
                             'Sig',
                             'Not Sig'),
         Direction.type_DvsS_2 = ifelse(logFC_DvsS_2 > 0,
                                        'up',
                                        'down'),
         Sig.direction_DvsS_2 = ifelse(Sig_DvsS_2 == 'Sig',
                                       Direction.type_DvsS_2,
                                       Sig_DvsS_2)) %>% 
  mutate(Sig_DvsS_3 = ifelse(adj.P.Val_DvsS_3 <= 0.05,
                             'Sig',
                             'Not Sig'),
         Direction.type_DvsS_3 = ifelse(logFC_DvsS_3 > 0,
                                        'up',
                                        'down'),
         Sig.direction_DvsS_3 = ifelse(Sig_DvsS_3 == 'Sig',
                                       Direction.type_DvsS_3,
                                       Sig_DvsS_3)) %>% 
  mutate(Sig_DvsS_4 = ifelse(adj.P.Val_DvsS_4 <= 0.05,
                             'Sig',
                             'Not Sig'),
         Direction.type_DvsS_4 = ifelse(logFC_DvsS_4 > 0,
                                        'up',
                                        'down'),
         Sig.direction_DvsS_4 = ifelse(Sig_DvsS_4 == 'Sig',
                                       Direction.type_DvsS_4,
                                       Sig_DvsS_4)) %>% 
  mutate(Sig_DvsS = ifelse(adj.P.Val_DvsS <= 0.05,
                           'Sig',
                           'Not Sig'),
         Direction.type_DvsS = ifelse(logFC_DvsS > 0,
                                      'up',
                                      'down'),
         Sig.direction_DvsS = ifelse(Sig_DvsS == 'Sig',
                                     Direction.type_DvsS,
                                     Sig_DvsS))


##graph volcano plot
# create comparison list
limma.genotpye.vector = c("DvsS_0",
                          "DvsS_1",
                          "DvsS_2",
                          "DvsS_3",
                          "DvsS_4",
                          "DvsS")

for (i in limma.genotpye.vector) {
  # graph volcano plot
  avp.neuron.reduce.genotype.vector.limma.results.df.sct.treat.prep %>% 
    mutate(sig.label = ifelse(get(paste0("Sig_", i)) == 'Sig',
                              Gene,
                              '')) %>% 
    ggplot(aes(x = get(paste0("logFC_", i)),
               y = -log10(get(paste0("adj.P.Val_", i))),
               color = get(paste0("Sig.direction_", i))))+
    geom_hline(yintercept = -log10(0.05),
               linetype = 'dotted') +
    geom_point(size = 5) +
    theme_classic() + 
    scale_color_manual(values=c("21B9CA", 
                                "grey", 
                                "orange3")) +
    geom_text(aes(label = sig.label),
              vjust = 0, 
              nudge_y = 0.10,
              size = 5) +
    theme(text = element_text(size = 20),
          legend.position = 'none') +
    xlab(paste0("logFC_", i)) +
    ylab( paste0("-log10(adj.P.Val_", i,")")) +
    ggtitle(paste0(i, " volcano plot"))
  ggsave(paste0('neuropeptides/avp.oxt/avp/limmatrend/genotype.sct.treat.prep/limma.avp.neurons.reduce.volcano.', i, '.png'),
         width = 5,
         height = 5)
}





##### poster heatmap ####
# pheatmap
# avp.neuron.reduce.genotype.vector.limma.results.df.sct.treat.prep  %>%
#   mutate(Keep = case_when(Sig_DvsS_0 != "Not Sig" ~ "Keep",
#                           Sig_DvsS_1 != "Not Sig" ~ "Keep",
#                           Sig_DvsS_2 != "Not Sig" ~ "Keep",
#                           Sig_DvsS_3 != "Not Sig" ~ "Keep",
#                           Sig_DvsS_4 != "Not Sig" ~ "Keep",
#                           Sig_DvsS != "Not Sig" ~ "Keep",
#                           TRUE ~ "Remove")) %>%
#   filter(Keep == "Keep") %>%
#   mutate(Gene = case_when(Gene == "ENSONIG00000040658" ~ "HBE1",
#                           Gene == "ENSONIG00000004376" ~ "EBF4",
#                           Gene == "ENSONIG00000036497" ~ "FBXL19",
#                           Gene == "ENSONIG00000010896" ~ "PBX3",
#                           Gene == "si:dkey-22o22.2" ~ "CDH2",
#                           TRUE ~ Gene)) %>%
#   column_to_rownames('Gene') %>%
#   mutate(All = case_when(Direction.type_DvsS == 'down' ~ log10(adj.P.Val_DvsS),
#                          Direction.type_DvsS == 'up' ~ -log10(adj.P.Val_DvsS),
#                          TRUE ~ 0),
#          '0' = case_when(Direction.type_DvsS_0 == 'down' ~ log10(adj.P.Val_DvsS_0),
#                          Direction.type_DvsS_0 == 'up' ~ -log10(adj.P.Val_DvsS_0),
#                          TRUE ~ 0),
#          '1' = case_when(Direction.type_DvsS_1 == 'down' ~ log10(adj.P.Val_DvsS_1),
#                          Direction.type_DvsS_1 == 'up' ~ -log10(adj.P.Val_DvsS_1),
#                          TRUE ~ 0),
#          '2' = case_when(Direction.type_DvsS_2 == 'down' ~ log10(adj.P.Val_DvsS_2),
#                          Direction.type_DvsS_2 == 'up' ~ -log10(adj.P.Val_DvsS_2),
#                          TRUE ~ 0),
#          '3' = case_when(Direction.type_DvsS_3 == 'down' ~ log10(adj.P.Val_DvsS_3),
#                          Direction.type_DvsS_3 == 'up' ~ -log10(adj.P.Val_DvsS_3),
#                          TRUE ~ 0),
#          '4' = case_when(Direction.type_DvsS_4 == 'down' ~ log10(adj.P.Val_DvsS_4),
#                          Direction.type_DvsS_4 == 'up' ~ -log10(adj.P.Val_DvsS_4),
#                          TRUE ~ 0)) %>%
#   select(c(All,
#            '0',
#            '1',
#            '2',
#            '3',
#            '4')) %>%
#   arrange('4',
#           '3',
#           '2',
#           '1',
#           '0',
#           All) %>%
#   as.matrix() %>%
#   t() %>%
#   pheatmap(cluster_rows = F,
#            cluster_cols = T,
#            scale = 'none',
#            border_color = 'black',
#            color = colorRampPalette(c("#60bb46", "white","#4e499e"))(100),
#            legend = F,
#            treeheight_col = 0,
#            treeheight_row = 0,
#            angle_col = 315,
#            fontsize = 15,
#            breaks = seq(-4, 4, length.out = 100),
#            filename = "./neuropeptides/avp.oxt/avp/limmatrend/Heatmap DEGs across clusters update.pdf",
#            width = 9.71,
#            height = 6.5
#   )
# dev.off()

# create matrix
avp.neuron.deg.matrix = avp.neuron.reduce.genotype.vector.limma.results.df.sct.treat.prep %>% 
  mutate(Keep = case_when(Sig_DvsS_0 != "Not Sig" ~ "Keep",
                          Sig_DvsS_1 != "Not Sig" ~ "Keep",
                          Sig_DvsS_2 != "Not Sig" ~ "Keep",
                          Sig_DvsS_3 != "Not Sig" ~ "Keep",
                          Sig_DvsS_4 != "Not Sig" ~ "Keep",
                          Sig_DvsS != "Not Sig" ~ "Keep",
                          TRUE ~ "Remove")) %>% 
  filter(Keep == "Keep") %>% 
  mutate(Gene = case_when(Gene == "ENSONIG00000040658" ~ "HBE1",
                          Gene == "ENSONIG00000004376" ~ "EBF4",
                          Gene == "ENSONIG00000036497" ~ "FBXL19",
                          Gene == "ENSONIG00000010896" ~ "PBX3",
                          Gene == "si:dkey-22o22.2" ~ "CDH2",
                          Gene == "ENSONIG00000034988" ~ 'MT: 1,086',
                          Gene == "ENSONIG00000010896" ~ 'Pbx3',
                          Gene == "ENSONIG00000002603" ~ 'GnRHR-II',
                          Gene == "ENSONIG00000031366" ~ 'MT: 70',
                          Gene == "ENSONIG00000004376" ~ "Ebf4",
                          Gene == "ENSONIG00000003414" ~ 'slc47a2.1',
                          Gene == "ENSONIG00000003310" ~ 'Epha5',
                          Gene == "im:7142702"        ~ 'Bcl11b', 
                          Gene == "ENSONIG00000040640"  ~ 'ZNF804A',      
                          Gene == "ENSONIG00000040652" ~ 'or106-8',
                          Gene == "ENSONIG00000041695" ~ 'CCBE1',
                          Gene == "ENSONIG00000036306" ~ '36306',
                          TRUE ~ Gene)) %>% 
  column_to_rownames('Gene') %>% 
  mutate(All = case_when(Direction.type_DvsS == 'down' ~ log10(adj.P.Val_DvsS),
                         Direction.type_DvsS == 'up' ~ -log10(adj.P.Val_DvsS),
                         TRUE ~ 0),
         '0' = case_when(Direction.type_DvsS_0 == 'down' ~ log10(adj.P.Val_DvsS_0),
                         Direction.type_DvsS_0 == 'up' ~ -log10(adj.P.Val_DvsS_0),
                         TRUE ~ 0),
         '1' = case_when(Direction.type_DvsS_1 == 'down' ~ log10(adj.P.Val_DvsS_1),
                         Direction.type_DvsS_1 == 'up' ~ -log10(adj.P.Val_DvsS_1),
                         TRUE ~ 0),
         '2' = case_when(Direction.type_DvsS_2 == 'down' ~ log10(adj.P.Val_DvsS_2),
                         Direction.type_DvsS_2 == 'up' ~ -log10(adj.P.Val_DvsS_2),
                         TRUE ~ 0),
         '3' = case_when(Direction.type_DvsS_3 == 'down' ~ log10(adj.P.Val_DvsS_3),
                         Direction.type_DvsS_3 == 'up' ~ -log10(adj.P.Val_DvsS_3),
                         TRUE ~ 0),
         '4' = case_when(Direction.type_DvsS_4 == 'down' ~ log10(adj.P.Val_DvsS_4),
                         Direction.type_DvsS_4 == 'up' ~ -log10(adj.P.Val_DvsS_4),
                         TRUE ~ 0)) %>% 
  select(c(All,
           '0',
           '1',
           '2',
           '3',
           '4')) %>% 
  arrange(All,
          '0',
          '1',
          '2',
          '3',
          '4') %>% 
  as.matrix() %>% 
  t() 

# create heatmap
avp.neuron.deg.phtmap <- pheatmap::pheatmap(avp.neuron.deg.matrix,
                                            cluster_rows = F,
                                            cluster_cols = T,
                                            scale = 'none',
                                            silent = T)

# get order list
avp.neuron.deg.order.list = c(avp.neuron.deg.phtmap$tree_col$order)
# create label list
avp.neuron.deg.order = data.frame(label = avp.neuron.deg.phtmap$tree_col$labels) %>% 
  rownames_to_column('position')
# combine labels and position
avp.neuron.deg.order = avp.neuron.deg.order[match(avp.neuron.deg.order.list,
                                                  avp.neuron.deg.order$position),] %>% 
  mutate(avp.neuron.deg.order.list = 1:length(avp.neuron.deg.order.list))
# create new order list
# move Dom bias genes to end
avp.neuron.deg.col_dend <- avp.neuron.deg.phtmap[[2]]
avp.neuron.deg.col_dend <- dendextend::rotate(avp.neuron.deg.col_dend, 
                                              order = rbind(avp.neuron.deg.order[c(7:9,13:30),],avp.neuron.deg.order[c(1:6,10:12,31:38),]) %>% pull(label) )
# graph pheat map
# for poster
pheatmap(avp.neuron.deg.matrix, 
         cluster_cols=as.hclust(avp.neuron.deg.col_dend),
         cluster_rows = F,
         gaps_row = 1,
         scale = 'none',
         border_color = 'black',
         color = colorRampPalette(c("#60bb46", "white","#4e499e"))(100),
         legend = T,
         treeheight_col = 0,
         treeheight_row = 0,
         angle_col = 315,
         fontsize = 15,
         breaks = seq(-2, 2, length.out = 100),
         filename = "./neuropeptides/avp.oxt/avp/limmatrend/Heatmap DEGs across clusters reorder poster update.png",
         width = 12,
         height = 5
)
dev.off()  






#### graph avp and OXT ####
## graph proportion of AVP cells that express Oxt
burtoni.snseq.combined.sct.all.avp.neurons.recluster.expression %>% 
  dplyr::select(Oxt.cell,
                integrated_snn_res.0.8,
                orig.ident) %>% 
  table() %>% 
  as.data.frame() %>% 
  filter(Freq != 0) %>% 
  ggplot(aes(x = integrated_snn_res.0.8,
             y = Freq,
             color = orig.ident,
             shape = Oxt.cell,
             group = paste0(orig.ident,Oxt.cell))) + 
  geom_line()+
  geom_point(size = 5) +
  theme_classic()+ 
  theme(text = element_text(size = 20)) +
  xlab('AVP neuron custer') +
  scale_color_manual(values = c("#4e499e",
                                "#60bb46"))
ggsave('./neuropeptides/avp.oxt/avp/AVP expressing oxt neuron by status and cluster.png',
       height = 10,
       width = 15)

## graph proportion of AVP cells that express Oxt
burtoni.snseq.combined.sct.all.avp.neurons.recluster.expression %>% 
  dplyr::select(Oxt.cell,
                integrated_snn_res.0.8,
                orig.ident) %>% 
  table() %>% 
  as.data.frame() %>% 
  filter(Freq != 0) %>% 
  pivot_wider(names_from = Oxt.cell,
              values_from = Freq) %>% 
  mutate(Percentage = 100*both/(avp+both)) %>% 
  ggplot(aes(x = integrated_snn_res.0.8,
             y = Percentage,
             color = orig.ident,
             group = orig.ident)) + 
  geom_line()+
  geom_point(size = 5) +
  theme_classic()+ 
  theme(text = element_text(size = 20)) +
  xlab('AVP neuron custer') +
  ylab('Percentage co-express Oxt') +
  ylim(c(0,100)) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#4e499e",
                                "#60bb46"))
  ggsave('./neuropeptides/avp.oxt/avp/AVP expressing oxt neuron by status and cluster percentage.pdf',
         height = 5,
         width = 5)


  
## graph mean UMI and read count per sample
AVP.neuron.read.UMI.stats  = burtoni.snseq.combined.sct.all.avp.neurons.recluster@meta.data %>%
    dplyr::select(Genotype.id,
                  integrated_snn_res.0.8,
                  orig.ident,
                  nCount_RNA,
                  nFeature_RNA,
                  nCount_SCT,
                  nFeature_SCT) %>% 
  group_by(Genotype.id,
           integrated_snn_res.0.8,
           orig.ident) %>% 
  summarise(nCount_RNA.mean = mean(nCount_RNA),
            nFeature_RNA.mean = mean(nFeature_RNA),
            nCount_SCT.mean = mean(nCount_SCT),
            nFeature_SCT.mean = mean(nFeature_SCT)) %>% 
  ungroup()
    
    
Features.list.avp.stats = c('nCount_RNA.mean',
                            'nCount_SCT.mean',
                            'nFeature_RNA.mean',
                            'nFeature_SCT.mean')

## graph
for (i in Features.list.avp.stats) {
  #genotype
  AVP.neuron.read.UMI.stats %>% 
    ggplot(aes(x = integrated_snn_res.0.8,
               y = get(i),
               color = orig.ident,
               group = Genotype.id)) + 
    geom_line()+
    geom_point(size = 5) +
    theme_classic()+ 
    theme(text = element_text(size = 20)) +
    xlab('AVP neuron custer') +
    ylab(paste('Mean',
               i,
               sep = ' ')) +
    scale_color_manual(values = c("#4e499e",
                                  "#60bb46")) 
  ggsave(paste('./neuropeptides/avp.oxt/avp/stats/',
               'Mean ',
               i,
               ' across genotype.png',
               sep = ''),
         height = 5,
         width = 10)
}    


  
  











