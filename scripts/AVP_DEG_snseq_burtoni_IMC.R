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
library(sp)
library(SeuratObject)
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
    dplyr::select(Cell.type,
           Cell.id) %>% 
    column_to_rownames(var = "Cell.id"),
  col.name = 'Cell.type'
)

## genotype
burtoni.snseq.combined.sct = AddMetaData(
  object = burtoni.snseq.combined.sct,
  metadata = burtoni.souporcell.filtered %>% 
    dplyr::select(Genotype.id,
           Cell.id) %>% 
    column_to_rownames(var = "Cell.id"),
  col.name = 'Genotype.id'
)

## sctype.hypo
burtoni.snseq.combined.sct = AddMetaData(
  object = burtoni.snseq.combined.sct,
  metadata = burtoni.sctypemarkers.hypo %>% 
    dplyr::select(sctypemarkers.hypo,
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

#### neurons ####
burtoni.snseq.combined.sct.neurons = burtoni.snseq.combined.sct

## keep cluster ids names
burtoni.snseq.combined.sct.neurons = AddMetaData(burtoni.snseq.combined.sct.neurons,
                                                 burtoni.snseq.combined.sct.neurons@meta.data$integrated_snn_res.0.8,
                                                        col.name = 'broad_integrated_snn_res.0.8')

#set idents
Idents(object = burtoni.snseq.combined.sct.neurons) <- "sctypemarkers.hypo"

#subset to neurons
burtoni.snseq.combined.sct.neurons = subset(burtoni.snseq.combined.sct.neurons,
                                            idents = c("C7-1: GLU",
                                                       "C7-2: GABA"))

# need to set to integrated for clustering
DefaultAssay(burtoni.snseq.combined.sct.neurons) = 'integrated'

#check data loaded correctly
## run PCA, UMAP, and cluster 
#use 0.8 resolution
burtoni.snseq.combined.sct.neurons = burtoni.snseq.combined.sct.neurons %>% 
  RunPCA() %>%
  FindNeighbors(dims = 1:15) %>%
  RunUMAP(dims = 1:15) %>%
  FindClusters(resolution = 0.8)

# need to set to integrated for clustering
DefaultAssay(burtoni.snseq.combined.sct.neurons) = 'SCT'

## graph 
# idents to new clusters
Idents(object = burtoni.snseq.combined.sct.neurons) <- "integrated_snn_res.0.8"

DimPlot(burtoni.snseq.combined.sct.neurons,
        group.by='sctypemarkers.hypo',
        label=TRUE) +
  ggtitle('Neurons') 

##expression level
# avp, with UMAP data, with meta data
burtoni.snseq.combined.sct.neurons.expression = full_join(full_join(burtoni.snseq.combined.sct.neurons@reductions$umap@cell.embeddings %>% 
                                                                                           as.data.frame() %>% 
                                                                                           rownames_to_column("Cell.id"),
                                                         burtoni.snseq.combined.sct.neurons@meta.data %>%
                                                                                           rownames_to_column("Cell.id")),
                                               burtoni.snseq.combined.sct.neurons@assays$SCT@data %>% 
                                                                                 as.data.frame() %>% 
                                                                                 filter(rownames(burtoni.snseq.combined.sct.neurons@assays$SCT@data) %in% c('avp')) %>% 
                                                                                 t() %>% as.data.frame() %>% 
                                                                                 rownames_to_column('Cell.id'))

# ## hoverlocator
# 
# HoverLocator(plot = DimPlot(burtoni.snseq.combined.sct.neurons), 
#              information = FetchData(object = burtoni.snseq.combined.sct.neurons,
#                                      vars = 'integrated_snn_res.0.8'))


#### AVP neurons ####
### avp neurons
### subset data to relevant cells
burtoni.snseq.combined.sct.all.avp.neurons = burtoni.snseq.combined.sct.neurons

## keep cluster ids names
burtoni.snseq.combined.sct.all.avp.neurons = AddMetaData(burtoni.snseq.combined.sct.all.avp.neurons,
                                                         burtoni.snseq.combined.sct.all.avp.neurons@meta.data$integrated_snn_res.0.8,
                                                 col.name = 'neurons_integrated_snn_res.0.8')

# subset with SCT data
DefaultAssay(burtoni.snseq.combined.sct.all.avp.neurons) = 'SCT'

#subset
burtoni.snseq.combined.sct.all.avp.neurons = subset(burtoni.snseq.combined.sct.all.avp.neurons,
                                                    subset = avp >= 1)


# cluster with integrated
DefaultAssay(burtoni.snseq.combined.sct.all.avp.neurons) = 'integrated'

## run PCA, UMAP, and cluster 
#use 0.8 resolution
burtoni.snseq.combined.sct.all.avp.neurons.recluster = burtoni.snseq.combined.sct.all.avp.neurons %>% 
  RunPCA() %>%
  FindNeighbors(dims = 1:15) %>%
  RunUMAP(dims = 1:15) %>%
  FindClusters(resolution = 0.8)

## keep cluster ids names
burtoni.snseq.combined.sct.all.avp.neurons.recluster = AddMetaData(burtoni.snseq.combined.sct.all.avp.neurons.recluster,
                                                                   burtoni.snseq.combined.sct.all.avp.neurons.recluster@meta.data$integrated_snn_res.0.8,
                                                         col.name = 'neurons.avp_integrated_snn_res.0.8')

#set idents
Idents(object = burtoni.snseq.combined.sct.all.avp.neurons.recluster) <- "integrated_snn_res.0.8"

# use SCT
DefaultAssay(burtoni.snseq.combined.sct.all.avp.neurons) = 'SCT'

# # graph
# # cell types
# DimPlot(burtoni.snseq.combined.sct.all.avp.neurons.recluster,
#         group.by='sctypemarkers.hypo',
#         label=TRUE) +
#   ggtitle('Neurons AVP')
# ggsave('neuropeptides/avp.oxt/avp/UMAP.sctype.hypo.avp.neurons.all.png',
#        width = 10,
#        height = 10)
# # clusters
# DimPlot(burtoni.snseq.combined.sct.all.avp.neurons.recluster,
#         label=TRUE) +
#   ggtitle('Neurons AVP')
# ggsave('neuropeptides/avp.oxt/avp/UMAP.clusters.avp.neurons.all.png',
#        width = 10,
#        height = 10)
# ## clusters
# # avp expression
# FeaturePlot(burtoni.snseq.combined.sct.all.avp.neurons.recluster,
#             reduction = "umap",
#             features = c('avp'),
#             min.cutoff = 1,
#             max.cutoff = 5,
#             label = TRUE,
#             repel = TRUE)
# ggsave('neuropeptides/avp.oxt/avp/Clusters.dimplot.avp.neurons.all.expression.png',
#        width = 10,
#        height = 10)

## clustree
# cluster across resolutions
burtoni.snseq.combined.sct.all.avp.neurons.clustree <- Seurat::FindClusters(object = burtoni.snseq.combined.sct.all.avp.neurons.recluster,
                                                                            resolution = resolution.range.reduced)
#check data
head(burtoni.snseq.combined.sct.all.avp.neurons.clustree[[]])

# #clustree
# clustree(burtoni.snseq.combined.sct.all.avp.neurons.clustree,
#          prefix = "integrated_snn_res.",
#          node_colour = 'cluster',
#          node_size_range = c(10,20),
#          scale_node_text = TRUE) +
#   scale_edge_color_continuous(low = "black",
#                               high = "black") +
#   theme(legend.position = "bottom")
# ggsave('neuropeptides/avp.oxt/avp/Clustree.avp.neurons.png',
#        width = 10,
#        height = 10)

# remove data
rm(burtoni.snseq.combined.sct.all.avp.neurons.clustree)


## count cells per cluster
burtoni.snseq.combined.sct.all.avp.neurons.recluster@meta.data %>% 
  dplyr::select(orig.ident,
                integrated_snn_res.0.8) %>% 
  table()

## combine clusters 4 and 5
# get metadata
burtoni.snseq.combined.sct.all.avp.neurons.recluster@meta.data = burtoni.snseq.combined.sct.all.avp.neurons.recluster@meta.data %>% 
  mutate(integrated_snn_res.0.8 = as.character(integrated_snn_res.0.8)) %>% 
  mutate(integrated_snn_res.0.8 = ifelse(integrated_snn_res.0.8 == 5,
                                         4,
                                         integrated_snn_res.0.8))


#subset to remove small noise clusters
Idents(burtoni.snseq.combined.sct.all.avp.neurons.recluster) <- "integrated_snn_res.0.8" 
cluster_values <- c(0,
                    1,
                    2,
                    3,
                    4)
burtoni.snseq.combined.sct.reduce.avp.neurons.recluster = subset(burtoni.snseq.combined.sct.all.avp.neurons.recluster, 
                                                                 idents = cluster_values,
                                                                 invert = FALSE)

## count cells per cluster
burtoni.snseq.combined.sct.reduce.avp.neurons.recluster@meta.data %>% 
  dplyr::select(orig.ident,
                integrated_snn_res.0.8) %>% 
  table()


# use SCT
DefaultAssay(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster) = 'SCT'

# ## graph
# # cell types
# DimPlot(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster,
#         group.by='sctypemarkers.hypo',
#         label=TRUE) +
#   ggtitle('Neurons AVP')
# ggsave('neuropeptides/avp.oxt/avp/UMAP.sctype.hypo.avp.neurons.png',
#        width = 10,
#        height = 10)
# # clusters
# DimPlot(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster,
#         label=TRUE) +
#   ggtitle('Neurons AVP')
# ggsave('neuropeptides/avp.oxt/avp/UMAP.clusters.avp.neurons.png',
#        width = 10,
#        height = 10)

##expression level
# avp, with UMAP data, with meta data
burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.expression = full_join(full_join(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster@reductions$umap@cell.embeddings %>% 
                                                                                           as.data.frame() %>% 
                                                                                           rownames_to_column("Cell.id"),
                                                                                         burtoni.snseq.combined.sct.reduce.avp.neurons.recluster@meta.data %>%
                                                                                           rownames_to_column("Cell.id")),
                                                                               burtoni.snseq.combined.sct.reduce.avp.neurons.recluster@assays$SCT@data %>% 
                                                                                 as.data.frame() %>% 
                                                                                 filter(rownames(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster@assays$SCT@data) %in% c('avp')) %>% 
                                                                                 t() %>% as.data.frame() %>% 
                                                                                 rownames_to_column('Cell.id'))

## project avp neuron clusters onto neuron umap
# add avp clusters
burtoni.snseq.combined.sct.neurons.expression = burtoni.snseq.combined.sct.neurons.expression %>% full_join(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.expression %>% 
                                                                                                              dplyr::select(Cell.id,
                                                                                                                            orig.ident,
                                                                                                                            neurons.avp_integrated_snn_res.0.8))

# #get average expression data
# avp.avg.exp.cluster.genotype = AverageExpression(AddMetaData(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster,
#                                                              metadata = paste(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster@meta.data$Genotype.id,
#                                                                               burtoni.snseq.combined.sct.reduce.avp.neurons.recluster@meta.data$integrated_snn_res.0.8,
#                                                                               sep = '_'),
#                                                              col.name = 'Genotype.cluster'),
#                                                  assays = 'SCT',
#                                                  features = c('avp'),
#                                                  return.seurat = FALSE,
#                                                  group.by = "Genotype.cluster",
#                                                  slot = "data",
#                                                  verbose = TRUE) %>% 
#   as.data.frame() %>% 
#   pivot_longer(cols = everything(),
#                names_to = 'ID',
#                values_to = 'Avg.expression') %>% 
#   separate(ID,
#            sep = '_',
#            into = c('Test',
#                     'Cluster')) %>% 
#   separate(Test,
#            sep = "\\.",
#            into = c('Data.type',
#                     'orig.ident',
#                     'Name')) %>%  
#   mutate(Genotype.id = paste(orig.ident,
#                              Name,
#                              sep = '.')) %>% 
#   dplyr::select(-c(Data.type)) 

#get average expression data
avp.avg.exp.cluster.genotype = AverageExpression(AddMetaData(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster,
                                                             metadata = paste(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster@meta.data$Genotype.id,
                                                                              burtoni.snseq.combined.sct.reduce.avp.neurons.recluster@meta.data$integrated_snn_res.0.8,
                                                                              sep = '_'),
                                                             col.name = 'Genotype.cluster'),
                                                 assays = 'SCT',
                                                 features = c('avp'),
                                                 return.seurat = FALSE,
                                                 group.by = "Genotype.cluster",
                                                 slot = "data",
                                                 verbose = TRUE) %>% 
  as.data.frame() %>% 
  pivot_longer(cols = everything(),
               names_to = 'ID',
               values_to = 'Avg.expression') %>% 
  separate_wider_delim(ID,
           delim = '.',
           names = c('Data.type',
                    'orig.ident',
                    'Genotype.id',
                    'Cluster')) %>% 
  mutate(Genotype.id = paste(orig.ident,
                             Genotype.id,
                             sep = '.')) %>% 
  dplyr::select(-c(Data.type)) 

## add avg expression
# avg expression per social status
avp.avg.exp.cluster.genotype.all = avp.avg.exp.cluster.genotype %>% 
  group_by(orig.ident,
           Cluster) %>% 
  mutate(Avg.avp.expression = mean(Avg.expression)) %>% 
  ungroup() %>% 
  dplyr::select(Cluster,
                Avg.avp.expression,
                Genotype.id) %>% 
  dplyr::rename(neurons.avp_integrated_snn_res.0.8 = Cluster)
# add to expression metadata frame
burtoni.snseq.combined.sct.neurons.expression.avg = avp.avg.exp.cluster.genotype.all %>% 
  right_join(burtoni.snseq.combined.sct.neurons.expression)

# change level orders
levels(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster) = c("0","1","2","3","4")


#### avp analysis #### 
## difference social status
# percent
# 254 dom to 236 sub
table(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.expression %>%
        select(integrated_snn_res.0.8,
               orig.ident)) %>%
  as.data.frame() %>%
  pivot_wider(names_from = orig.ident,
              values_from = Freq) %>%
  mutate(DomSubPerct = 100 * dom_burtoni_snseq/(sub_burtoni_snseq + dom_burtoni_snseq),
         TotalCount = sub_burtoni_snseq + dom_burtoni_snseq) %>%
  ggplot(aes(x= DomSubPerct,
             y = TotalCount,
             color = integrated_snn_res.0.8)) +
  geom_point(size = 5) +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  ggtitle('Clusters.DomvsSub.ratio.avp.neurons.reduce')
ggsave('neuropeptides/avp.oxt/avp/Clusters.DomvsSub.ratio.avp.neurons.reduce.png',
       width = 10,
       height = 10)

## violin plot
#for poster
VlnPlot(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster,
        features = c("avp"),
        split.by = 'orig.ident',
        group.by = 'integrated_snn_res.0.8',
        slot = "scale.data", #data and scale.data are flipped?
        pt.size = 0.5,
        assay = "SCT")+
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
ggsave('neuropeptides/avp.oxt/avp/AVP expression per cluster and social status.pdf',
       height = 5,
       width = 5,
       units = "in",
       dpi = 320)

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
  # theme(panel.border = element_rect(color = "black",
                                    # fill = NA,
                                    # size = 1))+
  theme(axis.text = element_text(size = 15))  +
  theme(axis.title = element_text(size = 20)) +
  xlab('AVP neuron cluster')
ggsave('neuropeptides/avp.oxt/avp/Average AVP expression per cluster and genotype poster update.pdf',
       height = 5,
       width = 5,
       units = "in",
       dpi = 320)

## project avp neuron clusters onto neuron umap
# graph
burtoni.snseq.combined.sct.neurons.expression %>% 
  ggplot(aes(x = umap_1,
             y = umap_2)) +
  geom_point(color = "grey") +
  geom_point(data = burtoni.snseq.combined.sct.neurons.expression %>% 
               filter(!is.na(neurons.avp_integrated_snn_res.0.8)),
             aes(color = neurons.avp_integrated_snn_res.0.8)) +
  theme_classic()
ggsave('neuropeptides/avp.oxt/avp/UMAP avp clusters projected neurons.png',
       height = 10,
       width = 15)

# # graph presentation
# burtoni.snseq.combined.sct.neurons.expression %>% 
#   ggplot(aes(x = umap_1,
#              y = umap_2)) +
#   geom_point(color = "grey") +
#   geom_point(data = burtoni.snseq.combined.sct.neurons.expression %>% 
#                filter(!is.na(neurons.avp_integrated_snn_res.0.8)),
#              aes(color = neurons.avp_integrated_snn_res.0.8)) +
#   theme_classic() +
#   NoLegend()
# ggsave('neuropeptides/avp.oxt/avp/UMAP avp clusters projected neurons presentation.png',
#        height = 7,
#        width = 7)
# 
# # graph presentation
# burtoni.snseq.combined.sct.neurons.expression %>% 
#   ggplot(aes(x = umap_1,
#              y = umap_2)) +
#   geom_point(color = "grey") +
#   geom_point(data = burtoni.snseq.combined.sct.neurons.expression %>% 
#                filter(!is.na(neurons.avp_integrated_snn_res.0.8)),
#              aes(color = neurons.avp_integrated_snn_res.0.8)) +
#   theme_classic() +
#   labs(color = 'AVP neuron cluster')
# ggsave('neuropeptides/avp.oxt/avp/UMAP avp clusters projected neurons presentation legend.png',
#        height = 10,
       # width = 15)

## graph paper
# add meta data to object
burtoni.snseq.combined.sct.neurons.add.avp = burtoni.snseq.combined.sct.neurons %>% 
  AddMetaData(burtoni.snseq.combined.sct.neurons.expression %>% 
                dplyr::select(Cell.id,
                              neurons.avp_integrated_snn_res.0.8) %>% 
                column_to_rownames('Cell.id') %>% 
                mutate(neurons.avp_integrated_snn_res.0.8 = as.character(neurons.avp_integrated_snn_res.0.8)),
              col.name = 'AVP.neuron.cluster')

# make list of AVP cells
AVP.neurons.list <- list("0" = burtoni.snseq.combined.sct.neurons.expression %>% 
                           filter((neurons.avp_integrated_snn_res.0.8== 0)) %>% 
                           pull(Cell.id),
                         "1" = burtoni.snseq.combined.sct.neurons.expression %>% 
                           filter((neurons.avp_integrated_snn_res.0.8== 1)) %>% 
                           pull(Cell.id) ,
                         "2" = burtoni.snseq.combined.sct.neurons.expression %>% 
                           filter((neurons.avp_integrated_snn_res.0.8== 2)) %>% 
                           pull(Cell.id),
                         "3" = burtoni.snseq.combined.sct.neurons.expression %>% 
                           filter((neurons.avp_integrated_snn_res.0.8== 3)) %>% 
                           pull(Cell.id) ,
                         "4" = burtoni.snseq.combined.sct.neurons.expression %>% 
                           filter((neurons.avp_integrated_snn_res.0.8== 4)) %>% 
                           pull(Cell.id) )
# graph cells
DimPlot(burtoni.snseq.combined.sct.neurons.add.avp,
        group.by = 'AVP.neuron.cluster',
        # alpha = 0,
        # label = T,
        sizes.highlight = 2,
        cells.highlight = AVP.neurons.list,
        # label.box = T,
        cols.highlight = c('#F8766D',
                           '#B79F00',
                           '#00BA38',
                           '#00BFC4',
                           '#619CFF',
                           '#F564E3')
        ) +
  NoLegend() +
  theme(
    plot.title = element_blank()
  )
ggsave('neuropeptides/avp.oxt/avp/UMAP avp clusters projected neurons paper.png',
       height = 7,
       width = 7)

# graph cells
DimPlot(burtoni.snseq.combined.sct.neurons.add.avp,
        group.by = 'AVP.neuron.cluster',
        alpha = 0,
        label = T,
        sizes.highlight = 2,
        cells.highlight = AVP.neurons.list,
        # label.box = T,
        cols.highlight = c('#F8766D',
                           '#B79F00',
                           '#00BA38',
                           '#00BFC4',
                           '#619CFF',
                           '#F564E3')
) +
  NoLegend() +
  theme(
    plot.title = element_blank()
  )
ggsave('neuropeptides/avp.oxt/avp/UMAP avp clusters projected neurons paper label.png',
       height = 7,
       width = 7)
# 
# # graph
# burtoni.snseq.combined.sct.neurons.expression.avg %>% 
#   ggplot(aes(x = umap_1,
#              y = umap_2)) +
#   geom_point(color = "grey") +
#   geom_point(data = burtoni.snseq.combined.sct.neurons.expression.avg %>% 
#                filter(!is.na(neurons.avp_integrated_snn_res.0.8)),
#              aes(color = Avg.avp.expression)) +
#   theme_classic() +
#   labs(color = 'AVP neuron cluster') +
#   scale_color_gradientn(colors = c('lightpink', 'red4'))
# ggsave('neuropeptides/avp.oxt/avp/UMAP avp avg expression projected neurons.png',
#        height = 10,
#        width = 15)

### alluvial 
## avp neuron cluster by neuron cluster
# color cluster
burtoni.snseq.combined.sct.reduce.avp.neurons.recluster@meta.data %>% 
  droplevels() %>% 
  select(c(neurons.avp_integrated_snn_res.0.8,
           neurons_integrated_snn_res.0.8)) %>% 
  table() %>% 
  as.data.frame() %>% 
  ggplot(aes(axis1 = neurons_integrated_snn_res.0.8,
             axis2 = neurons.avp_integrated_snn_res.0.8,
             y = Freq)) +
  geom_alluvium(aes(fill = neurons.avp_integrated_snn_res.0.8)) +
  geom_stratum() +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("neurons_integrated_snn_res.0.8", "neurons.avp_integrated_snn_res.0.8"),
                   expand = c(0.15, 0.05)) +
  scale_fill_viridis_d() +
  theme_classic() +
  theme(axis.text.x=element_blank())
ggsave('./neuropeptides/avp.oxt/avp/Alluvial AVP neuron cluster by neuron cluster.png',
       height = 10,
       width = 10)

## graph mean UMI and read count per sample
AVP.neuron.read.UMI.stats  = burtoni.snseq.combined.sct.reduce.avp.neurons.recluster@meta.data %>%
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
  ggsave(paste('neuropeptides/avp.oxt/avp/stats/',
               'Mean ',
               i,
               ' across genotype.png',
               sep = ''),
         height = 5,
         width = 10)
}    


## dotplot
# KCNMB4
DotPlot(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster,
        features = c('KCNMB4'
                     # 'avp'
                     ),
        assay = 'SCT',
        # group.by = 'neurons.avp_integrated_snn_res.0.8'
        )
ggsave('./neuropeptides/avp.oxt/avp/Dotplot avp gigantocellular cluster.png',
       height = 10,
       width = 10)

# otpa
DotPlot(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster,
        features = c('otpa'),
        assay = 'SCT',
        col.min = 0.5
        # group.by = 'neurons.avp_integrated_snn_res.0.8'
        )
ggsave('./neuropeptides/avp.oxt/avp/Dotplot avp otpa cluster.png',
       height = 10,
       width = 10)

# otpa by status
DotPlot(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster,
        features = c('otpa'),
        assay = 'SCT',
        # group.by = 'neurons.avp_integrated_snn_res.0.8',
        split.by = 'orig.ident',
        col.min = 0.5,
        cols = c('blue',
                 'blue')
)
ggsave('./neuropeptides/avp.oxt/avp/Dotplot avp otpa cluster status.png',
       height = 10,
       width = 10)

# otpa by genotype
DotPlot(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster,
        features = c('otpa'),
        assay = 'SCT',
        # group.by = 'neurons.avp_integrated_snn_res.0.8',
        split.by = 'Genotype.id',
        cols = c('blue',
                 'blue',
                 'blue',
                 'blue',
                 'blue',
                 'blue'),
        col.min = 0.5
)
ggsave('./neuropeptides/avp.oxt/avp/Dotplot avp otpa cluster genotype.png',
       height = 10,
       width = 10)

# IEG
DotPlot(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster,
        features = c('egr1',
                     'npas4',
                     'fosb',
                     'fosab',
                     'fosaa'),
        assay = 'SCT',
        col.min = 0.5,
        # group.by = 'neurons.avp_integrated_snn_res.0.8'
        )
ggsave('./neuropeptides/avp.oxt/avp/Dotplot avp IEG cluster.png',
       height = 10,
       width = 10)

# IEG by status
DotPlot(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster,
        features = c('egr1',
                     'npas4',
                     'fosb',
                     'fosab',
                     'fosaa'),
        assay = 'SCT',
        # split.by = 'neurons.avp_integrated_snn_res.0.8',
        split.by = 'orig.ident',
        col.min = 0.5,
        cols = c('blue',
                 'blue')
        )
ggsave('./neuropeptides/avp.oxt/avp/Dotplot avp IEG cluster status.png',
       height = 10,
       width = 10)


# IEG by status
DotPlot(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster,
        features = c('egr1',
                     'npas4',
                     'fosb',
                     'fosab',
                     'fosaa'),
        assay = 'SCT',
        # split.by = 'neurons.avp_integrated_snn_res.0.8',
        split.by = 'orig.ident',
        col.min = 0.5,
        cols = c('blue',
                 'blue')
)
ggsave('./neuropeptides/avp.oxt/avp/Dotplot avp IEG cluster status.png',
       height = 10,
       width = 10)

# IEG by status
DotPlot(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster,
        features = c('egr1',
                     'npas4',
                     'fosb'),
        assay = 'SCT',
        # split.by = 'neurons.avp_integrated_snn_res.0.8',
        split.by = 'orig.ident',
        col.min = 0,
        cols = 'Reds',
        dot.min = .05
        )+ 
  scale_color_distiller(palette = 'Reds',
                        direction = 0) +
  scale_y_discrete(labels = c('Dom 0',
                              'Sub 0',
                              'Dom 1',
                              'Sub 1',
                              'Dom 2',
                              'Sub 2',
                              'Dom 3',
                              'Sub 3',
                              'Dom 4',
                              'Sub 4')) +
  ylab('AVP neuron cluster by status') +
  xlab('IEGs')
ggsave('./neuropeptides/avp.oxt/avp/Dotplot avp IEG cluster status paper.png',
       height = 5,
       width = 5)

# check neurons
DotPlot(burtoni.snseq.combined.sct.neurons,
        features = c('egr1',
                     'npas4',
                     'fosb',
                     'fosab',
                     'fosaa'),
        assay = 'SCT',
        col.min = 0.5,
        # group.by = 'integrated_snn_res.0.8'
        )
ggsave('./neurons/Dotplot avp IEG cluster.png',
       height = 10,
       width = 10)

#### AVP and oxt analysis ####
### graph avp and oxt
# VlnPlot(burtoni.snseq.combined.sct.neurons,
#         features = c('avp',
#                      'oxt'),
#         assay = 'SCT',
#         slot = 'data',
#         group.by = 'orig.ident')

### get meta data and expression data for oxt/avp cells
burtoni.snseq.combined.sct.neurons.expression.avp.oxt = full_join(full_join(burtoni.snseq.combined.sct.neurons@reductions$umap@cell.embeddings %>% 
                                                                              as.data.frame() %>% 
                                                                              rownames_to_column("Cell.id"),
                                                                    burtoni.snseq.combined.sct.neurons@meta.data %>%
                                                                              rownames_to_column("Cell.id")),
                                                          burtoni.snseq.combined.sct.neurons@assays$SCT@counts %>% 
                                                                    as.data.frame() %>% 
                                                                    filter(rownames(burtoni.snseq.combined.sct.neurons@assays$SCT@counts) %in% c('avp',
                                                                                                                                                       'oxt')) %>% 
                                                                    t() %>% as.data.frame() %>% 
                                                                    rownames_to_column('Cell.id'))

## assign AVP and oxt cell type cell type
burtoni.snseq.combined.sct.neurons.expression.avp.oxt = burtoni.snseq.combined.sct.neurons.expression.avp.oxt %>% 
  mutate(AVP.oxt.cell = case_when(avp >= 2 & oxt >= 2 ~ "Both",
                                  avp >= 2 & oxt < 2 ~ "Avp",
                                  avp < 2 & oxt >= 2 ~ "Oxt",
                                  avp < 2 & oxt < 2 ~ "None",
                                  TRUE ~ 'NA'))

## calculate average separately 
burtoni.snseq.combined.sct.neurons.expression.avp.oxt.avg = burtoni.snseq.combined.sct.neurons.expression.avp.oxt %>% 
  group_by(integrated_snn_res.0.8,
           orig.ident) %>% 
  mutate(Avg.oxt.expression = mean(oxt, na.rm = T),
         Avg.avp.expression = mean(avp, na.rm = T)) %>% 
  ungroup() %>% 
  mutate(Avg.avp.expression.scale = 100*Avg.avp.expression/max(Avg.avp.expression, na.rm = T),
         Avg.oxt.expression.scale = 100*Avg.oxt.expression/max(Avg.oxt.expression, na.rm = T))
  


#scale avp and oxt levels
burtoni.snseq.combined.sct.neurons.expression.avp.oxt = burtoni.snseq.combined.sct.neurons.expression.avp.oxt %>% 
  mutate(avp.scale = 100*avp/max(avp),
         oxt.scale = 100*oxt/max(oxt),
         both.avg = (avp+oxt)/2,
         both.avg.scale = (avp.scale+oxt.scale)/2)


# add avp averages
burtoni.snseq.combined.sct.neurons.expression.avp.oxt = burtoni.snseq.combined.sct.neurons.expression.avp.oxt %>% 
  full_join(burtoni.snseq.combined.sct.neurons.expression.avg %>% 
              dplyr::select(-c(avp)))

# calculate oxt averages
# use neuron clusters
#get average expression data
oxt.avg.exp.cluster.genotype = AverageExpression(AddMetaData(burtoni.snseq.combined.sct.neurons,
                                                             metadata = paste(burtoni.snseq.combined.sct.neurons@meta.data$Genotype.id,
                                                                              burtoni.snseq.combined.sct.neurons@meta.data$integrated_snn_res.0.8,
                                                                              sep = '_'),
                                                             col.name = 'Genotype.cluster'),
                                                 assays = 'SCT',
                                                 features = c('oxt'),
                                                 return.seurat = FALSE,
                                                 group.by = "Genotype.cluster",
                                                 slot = "counts",
                                                 verbose = TRUE) %>% 
  as.data.frame() %>% 
  pivot_longer(cols = everything(),
               names_to = 'ID',
               values_to = 'Avg.oxt.expression') %>% 
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
  dplyr::select(-c(Data.type)) 

## add avg expression
# filter low expressing values 
# avg expression per social status
oxt.avg.exp.cluster.genotype.all = oxt.avg.exp.cluster.genotype %>% 
  group_by(orig.ident,
           Cluster) %>% 
  mutate(Avg.oxt.expression = mean(Avg.oxt.expression)) %>% 
  ungroup() %>% 
  dplyr::select(Cluster,
                Avg.oxt.expression,
                Genotype.id) %>% 
  dplyr::rename(integrated_snn_res.0.8 = Cluster)

# add to expression metadata frame
burtoni.snseq.combined.sct.neurons.expression.avp.oxt = burtoni.snseq.combined.sct.neurons.expression.avp.oxt %>% 
  full_join(oxt.avg.exp.cluster.genotype.all)

#scale average avp and oxt levels
burtoni.snseq.combined.sct.neurons.expression.avp.oxt = burtoni.snseq.combined.sct.neurons.expression.avp.oxt %>% 
  mutate(Avg.avp.expression.scale = 100*Avg.avp.expression/max(Avg.avp.expression, na.rm = T),
         Avg.oxt.expression.scale = 100*Avg.oxt.expression/max(Avg.oxt.expression, na.rm = T),
         Avg.both.avg = (Avg.avp.expression+Avg.oxt.expression)/2,
         Avg.both.avg.scale = (Avg.avp.expression.scale+Avg.oxt.expression.scale)/2)

#### AVP and oxt graph ####
##graph oxt per cluster across social status
oxt.avg.exp.cluster.genotype %>%
  ggplot(aes(x = reorder(Cluster,
                         -Avg.oxt.expression),
             y = Avg.oxt.expression,
             color = orig.ident,
             group = orig.ident)) +
  geom_hline(yintercept = 2)+
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
  ylim(c(0,max(oxt.avg.exp.cluster.genotype$Avg.oxt.expression)))+
  theme_classic() +
  # theme(legend.position = 'none') +
  scale_fill_manual(values = c("#4e499e",
                               "#60bb46"))+
  scale_color_manual(values = c("#4e499e",
                                "#60bb46"))+
  ylab("OXT average expression") +
  # theme(panel.border = element_rect(color = "black",
  # fill = NA,
  # size = 1))+
  theme(axis.text = element_text(size = 15))  +
  theme(axis.title = element_text(size = 20)) +
  xlab('Neuron cluster')
ggsave('neuropeptides/avp.oxt/oxt/Average oxt expression per cluster and genotype.png',
       height = 5,
       width = 10)

#scaled
oxt.avg.exp.cluster.genotype %>%
  mutate(Avg.oxt.expression = log1p(Avg.oxt.expression)) %>% 
  ggplot(aes(x = reorder(Cluster,
                         -Avg.oxt.expression),
             y = Avg.oxt.expression,
             color = orig.ident,
             group = orig.ident)) +
  geom_hline(yintercept = log1p(2))+
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
  ylim(c(0,max(log1p(oxt.avg.exp.cluster.genotype$Avg.oxt.expression))))+
  theme_classic() +
  # theme(legend.position = 'none') +
  scale_fill_manual(values = c("#4e499e",
                               "#60bb46"))+
  scale_color_manual(values = c("#4e499e",
                                "#60bb46"))+
  ylab("OXT average expression") +
  # theme(panel.border = element_rect(color = "black",
  # fill = NA,
  # size = 1))+
  theme(axis.text = element_text(size = 15))  +
  theme(axis.title = element_text(size = 20)) +
  xlab('Neuron cluster') +
  ylab('OXT average expression scaled')
ggsave('neuropeptides/avp.oxt/oxt/Average oxt expression per cluster and genotype scaled.png',
       height = 5,
       width = 10)

# projection 
# cluster avg
# avp: avp neuron clusters
# oxt: neuron clusters
# need new color scale
library(ggnewscale)
burtoni.snseq.combined.sct.neurons.expression.avp.oxt %>% 
  ggplot(aes(x= umap_1,
             umap_2)) +
  geom_point(data = burtoni.snseq.combined.sct.neurons.expression.avp.oxt %>% 
               filter(AVP.oxt.cell == 'None'),
             color = 'grey',
             size = 2) +
  geom_point(data = burtoni.snseq.combined.sct.neurons.expression.avp.oxt %>% 
               filter(AVP.oxt.cell == 'Avp' | AVP.oxt.cell == 'Both'),
             aes(color = Avg.avp.expression.scale,
                 size = Avg.avp.expression.scale)) +
  scale_colour_gradientn(colours = c("lightpink", "red4"),
                         limits = c(0,100)) +
  scale_size(range = c(1,2))+
  new_scale_colour() +
  new_scale('size') +
  theme_classic() +
burtoni.snseq.combined.sct.neurons.expression.avp.oxt %>% 
  ggplot(aes(x= umap_1,
             umap_2)) +
  geom_point(data = burtoni.snseq.combined.sct.neurons.expression.avp.oxt %>% 
               filter(AVP.oxt.cell == 'None'),
             color = 'grey',
             size = 2) +
  geom_point(data = burtoni.snseq.combined.sct.neurons.expression.avp.oxt %>% 
               filter(AVP.oxt.cell == 'Oxt'| AVP.oxt.cell == 'Both'),
             aes(color = Avg.oxt.expression.scale,
                 size = Avg.oxt.expression.scale,
                 alpha = Avg.oxt.expression.scale)) +
  scale_colour_gradientn(colours = c("dodgerblue", "blue4"),
                         limits = c(0,100)) +
  scale_size(range = c(1,2))+
  new_scale_colour() +
  new_scale('size') +
  theme_classic()
ggsave('neuropeptides/avp.oxt/UMAP projection avp vs oxt average expression.png',
       height = 10,
       width = 20)

## for paper
# labeled by average expression
burtoni.snseq.combined.sct.neurons.expression.avp.oxt.avg %>% 
  ggplot(aes(x= umap_1,
             umap_2)) +
  geom_point(data = burtoni.snseq.combined.sct.neurons.expression.avp.oxt.avg %>% 
               filter(AVP.oxt.cell == 'None'),
             color = 'grey',
             size = 2) +
  geom_point(data = burtoni.snseq.combined.sct.neurons.expression.avp.oxt.avg %>% 
               filter(AVP.oxt.cell == 'Avp' | AVP.oxt.cell == 'Both'),
             aes(color = Avg.avp.expression.scale,
                 size = Avg.avp.expression.scale)) +
  scale_colour_gradientn(colours = c("lightpink", "red4"),
                         limits = c(0,100)
                         ) +
  scale_size(range = c(1,2))+
  new_scale_colour() +
  new_scale('size') +
  theme_classic() +
  theme(legend.position = 'none') +
  burtoni.snseq.combined.sct.neurons.expression.avp.oxt.avg %>% 
  ggplot(aes(x= umap_1,
             umap_2)) +
  geom_point(data = burtoni.snseq.combined.sct.neurons.expression.avp.oxt.avg %>% 
               filter(AVP.oxt.cell == 'None'),
             color = 'grey',
             size = 2) +
  geom_point(data = burtoni.snseq.combined.sct.neurons.expression.avp.oxt.avg %>% 
               filter(AVP.oxt.cell == 'Oxt'| AVP.oxt.cell == 'Both'),
             aes(color = Avg.oxt.expression.scale,
                 size = Avg.oxt.expression.scale,
                 alpha = Avg.oxt.expression.scale)) +
  scale_colour_gradientn(colours = c("dodgerblue", "blue4"),
                         limits = c(0,100)
                         ) +
  scale_size(range = c(1,2))+
  new_scale_colour() +
  new_scale('size') +
  theme_classic() +
  theme(legend.position = 'none') 
ggsave('neuropeptides/avp.oxt/UMAP projection avp vs oxt average expression paper.png',
       height = 6.5,
       width = 13)

# with legends
burtoni.snseq.combined.sct.neurons.expression.avp.oxt.avg %>% 
  ggplot(aes(x= umap_1,
             umap_2)) +
  geom_point(data = burtoni.snseq.combined.sct.neurons.expression.avp.oxt.avg %>% 
               filter(AVP.oxt.cell == 'None'),
             color = 'grey',
             size = 2) +
  geom_point(data = burtoni.snseq.combined.sct.neurons.expression.avp.oxt.avg %>% 
               filter(AVP.oxt.cell == 'Avp' | AVP.oxt.cell == 'Both'),
             aes(color = Avg.avp.expression.scale,
                 size = Avg.avp.expression.scale)) +
  scale_colour_gradientn(colours = c("lightpink", "red4"),
                         limits = c(0,100)
  ) +
  scale_size(range = c(1,2))+
  new_scale_colour() +
  new_scale('size') +
  theme_classic() +
  burtoni.snseq.combined.sct.neurons.expression.avp.oxt.avg %>% 
  ggplot(aes(x= umap_1,
             umap_2)) +
  geom_point(data = burtoni.snseq.combined.sct.neurons.expression.avp.oxt.avg %>% 
               filter(AVP.oxt.cell == 'None'),
             color = 'grey',
             size = 2) +
  geom_point(data = burtoni.snseq.combined.sct.neurons.expression.avp.oxt.avg %>% 
               filter(AVP.oxt.cell == 'Oxt'| AVP.oxt.cell == 'Both'),
             aes(color = Avg.oxt.expression.scale,
                 size = Avg.oxt.expression.scale,
                 alpha = Avg.oxt.expression.scale)) +
  scale_colour_gradientn(colours = c("dodgerblue", "blue4"),
                         limits = c(0,100)
  ) +
  scale_size(range = c(1,2))+
  new_scale_colour() +
  new_scale('size') +
  theme_classic() 
ggsave('neuropeptides/avp.oxt/UMAP projection avp vs oxt average expression paper legend.png',
       height = 6.5,
       width = 13)


# compare overlap
burtoni.snseq.combined.sct.neurons.expression.avp.oxt %>% 
  ggplot(aes(x= umap_1,
             umap_2)) +
  geom_point(data = burtoni.snseq.combined.sct.neurons.expression.avp.oxt %>% 
               filter(AVP.oxt.cell == 'None'),
             color = 'grey',
             size = 2) +
  geom_point(data = burtoni.snseq.combined.sct.neurons.expression.avp.oxt %>% 
               filter(AVP.oxt.cell == 'Both'),
             aes(color = Avg.both.avg.scale,
                 size = Avg.both.avg.scale)) +
  scale_colour_gradientn(colours = c("orchid", "purple4"),
                         limits = c(0,100)) +
  scale_size(range = c(1,2))+
  new_scale_colour() +
  new_scale('size') +
  theme_classic() +
burtoni.snseq.combined.sct.neurons.expression.avp.oxt %>% 
  ggplot(aes(x= umap_1,
             umap_2)) +
  geom_point(data = burtoni.snseq.combined.sct.neurons.expression.avp.oxt %>% 
               filter(AVP.oxt.cell == 'None'),
             color = 'grey',
             size = 2) +
  geom_point(data = burtoni.snseq.combined.sct.neurons.expression.avp.oxt %>% 
               filter(AVP.oxt.cell == 'Oxt'| AVP.oxt.cell == 'Both'),
             aes(color = Avg.oxt.expression.scale,
                 size = Avg.oxt.expression.scale,
                 alpha = Avg.oxt.expression.scale)) +
  scale_colour_gradientn(colours = c("dodgerblue", "blue4"),
                         limits = c(0,100)) +
  scale_size(range = c(1,2))+
  new_scale_colour() +
  new_scale('size') +
  new_scale('alpha') +
  geom_point(data = burtoni.snseq.combined.sct.neurons.expression.avp.oxt %>% 
               filter(AVP.oxt.cell == 'Avp' | AVP.oxt.cell == 'Both'),
             aes(color = Avg.avp.expression.scale,
                 size = Avg.avp.expression.scale)) +
  scale_colour_gradientn(colours = c("lightpink", "red4"),
                         limits = c(0,100)) +
  scale_size(range = c(1,2))+
  new_scale_colour() +
  new_scale('size') +
  geom_point(data = burtoni.snseq.combined.sct.neurons.expression.avp.oxt %>% 
               filter(AVP.oxt.cell == 'Both'),
             aes(color = Avg.both.avg.scale,
                 size = Avg.both.avg.scale)) +
  scale_colour_gradientn(colours = c("orchid", "purple4"),
                         limits = c(0,100)) +
  scale_size(range = c(1,2)) +
  theme_classic()
ggsave('neuropeptides/avp.oxt/UMAP projection avp vs oxt average expression overlap.png',
       height = 10,
       width = 20)



## graph results
burtoni.snseq.combined.sct.neurons.expression.avp.oxt %>% 
  dplyr::select(c(AVP.oxt.cell,
                  orig.ident)) %>% 
  table() %>% 
  as.data.frame() %>% 
  group_by(orig.ident) %>% 
  mutate(Total = sum(Freq)) %>% 
  ungroup() %>% 
  mutate(Percent = 100*Freq/Total) %>% 
  ggplot(aes(x = orig.ident,
             y = Percent,
             fill = AVP.oxt.cell,
             group = AVP.oxt.cell)) +
  geom_bar(position="dodge",
           stat = 'identity') +
  labs(fill = 'Neuron type')  +
  theme_classic()
ggsave('neuropeptides/avp.oxt/avp/Percentage neurons expressing AVP and OXT.png',
       height = 10,
       width = 10)


# just avp cells
## graph results
burtoni.snseq.combined.sct.neurons.expression.avp.oxt %>% 
  dplyr::select(c(AVP.oxt.cell,
                  Genotype.id,
                  orig.ident)) %>% 
  table() %>% 
  as.data.frame() %>%
  filter(AVP.oxt.cell == c('Avp', 
                           'Both')) %>%  
  group_by(orig.ident) %>% 
  mutate(Total = sum(Freq)) %>% 
  ungroup() %>% 
  mutate(Percent = 100*Freq/Total) %>% 
  filter(Percent != 0) %>% 
  ggplot(aes(x = AVP.oxt.cell,
             y = Percent,
             color = orig.ident,
             group = AVP.oxt.cell)) +
  geom_point() +
  labs(color = 'AVP neuron type') +
  theme_classic()
ggsave('neuropeptides/avp.oxt/avp/Percentage AVP neurons expressing OXT.png',
       height = 10,
       width = 10)

# count
burtoni.snseq.combined.sct.neurons.expression.avp.oxt %>% 
  dplyr::select(c(AVP.oxt.cell,
                  orig.ident)) %>% 
  table() %>% 
  as.data.frame() %>%
  filter(AVP.oxt.cell == c('Avp', 
                           'Both')) %>%  
  group_by(orig.ident) %>% 
  mutate(Total = sum(Freq)) %>% 
  ungroup() %>% 
  mutate(Percent = 100*Freq/Total) %>% 
  filter(Percent != 0) %>% 
  ggplot(aes(x = orig.ident,
             y = Freq,
             fill = AVP.oxt.cell,
             group = AVP.oxt.cell)) +
  geom_bar(position="dodge",
           stat = 'identity') +
  labs(fill = 'AVP neuron type') +
  theme_classic()
ggsave('neuropeptides/avp.oxt/avp/Count AVP neurons expressing OXT.png',
       height = 10,
       width = 10)






### get meta data and expression data for avp cells
burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.avp.oxt.expression = full_join(full_join(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster@reductions$umap@cell.embeddings %>% 
                                                                              as.data.frame() %>% 
                                                                              rownames_to_column("Cell.id"),
                                                                              burtoni.snseq.combined.sct.reduce.avp.neurons.recluster@meta.data %>%
                                                                              rownames_to_column("Cell.id")),
                                                                            burtoni.snseq.combined.sct.reduce.avp.neurons.recluster@assays$SCT@data %>% 
                                                                    as.data.frame() %>% 
                                                                    filter(rownames(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster@assays$SCT@data) %in% c('avp',
                                                                                                                                               'oxt')) %>% 
                                                                    t() %>% as.data.frame() %>% 
                                                                    rownames_to_column('Cell.id'))


## assign AVP and oxt cell type cell type
burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.avp.oxt.expression = burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.avp.oxt.expression %>% 
  mutate(AVP.oxt.cell = case_when(avp > 1 & oxt > 1 ~ "Both",
                                  avp > 1 & oxt < 1 ~ "Avp",
                                  avp < 1 & oxt > 1 ~ "Oxt",
                                  avp < 1 & oxt < 1 ~ "None",
                                  TRUE ~ 'NA'))






## bar chart of percent expressing both
# create genotype.id and orig.ident matrix
Genotype.id.orig.ident = burtoni.snseq.combined.sct.reduce.avp.neurons.recluster@meta.data %>% 
  dplyr::select(orig.ident,
                Genotype.id) %>% 
  distinct()

# graph
burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.avp.oxt.expression %>% 
  dplyr::select(c(AVP.oxt.cell,
                  Genotype.id,
                  neurons.avp_integrated_snn_res.0.8)) %>% 
  table() %>% 
  as.data.frame() %>%
  filter(AVP.oxt.cell == c('Avp', 
                           'Both')) %>%  
  group_by(Genotype.id,
           neurons.avp_integrated_snn_res.0.8) %>% 
  mutate(Total = sum(Freq)) %>% 
  ungroup() %>% 
  mutate(Percent = 100*Freq/Total) %>% 
  filter(AVP.oxt.cell == 'Both') %>% 
  filter(!is.na(Percent)) %>%  
  left_join(Genotype.id.orig.ident) %>%
  group_by(orig.ident,
           neurons.avp_integrated_snn_res.0.8) %>% 
  mutate(mean = mean(Percent),
         se = sd(Percent) / sqrt(length(.))) %>% 
  ungroup() %>% 
  ggplot(aes(x = neurons.avp_integrated_snn_res.0.8,
             group = orig.ident)) +
  geom_bar(aes(y=mean,
               fill = orig.ident), 
            stat="identity", 
            position="dodge") +
  geom_errorbar( aes(ymin=mean-se, 
                     ymax=mean+se), 
                 width=0.4, 
                 colour="black", 
                 alpha=0.9, 
                 size=1.3,
                 position = position_dodge(width = 0.9))+
  theme_classic() +
  ylab('Percentage AVP neurons OXT+') +
  xlab('AVP neuron clusters') +
  scale_fill_manual(values = c("#4e499e",
                               "#60bb46")) + 
  theme(text = element_text(size = 20))
ggsave('neuropeptides/avp.oxt/avp/Bar chart overlap AVP neurons expressing OXT.png',
       height = 10,
       width = 15)




#### Cluster 17 identification ####
### get marker orthologs from hypomap C66-22 cluster
# use Hypomap supplemental material "5_hypomap_marker_genes"
# filter out cluster C66-22 from 'cluster_id' and make new excel csv
C66.22.data = read.csv('../Gene.lists/Cluster C66 22 hypomap.csv')

#### use ensembl to get orthologs
## load ensembl
library(biomaRt)
#selecting biomart database
ensembl = useMart('ensembl')
# #list datasets
# dataset = listDatasets(ensembl)
##select dataset
#mouse
ensembl.mouse = useDataset('mmusculus_gene_ensembl',
                             mart=ensembl)

# # get list of all attributes
# listAttributes(ensembl.mouse) %>% View


#create attributes lists
tilapia.attributes = c('oniloticus_homolog_ensembl_gene',
                       'oniloticus_homolog_associated_gene_name',
                       'external_gene_name')

# use gene names
tilapia.C66.22.gene = getBM(attributes = tilapia.attributes,
                                    mart = ensembl.mouse,
                                    values = C66.22.data$gene,
                                    filter = 'external_gene_name',
                                    useCache = FALSE) # useCache has to do with version of R not being up to date?

# combine gene and ensembl gene IDs
tilapia.C66.22.gene.list = data.frame(Gene.orthologs = c(tilapia.C66.22.gene$oniloticus_homolog_ensembl_gene,
                                                         tilapia.C66.22.gene$oniloticus_homolog_associated_gene_name))

# remove empty rows
tilapia.C66.22.gene.list = tilapia.C66.22.gene.list %>% 
  mutate(Gene.orthologs = na_if(Gene.orthologs, "")) %>% 
  na.omit() %>% 
  distinct()

## compare to gene list 
tilapia.C66.22.gene.list = tilapia.C66.22.gene.list %>% 
  filter(Gene.orthologs %in% rownames(burtoni.snseq.combined.sct.neurons@assays$SCT@counts))

## marker genes from ensembl for hypomap C66-22 cluster
geneSet.all = tilapia.C66.22.gene.list %>% 
  pull(Gene.orthologs)

## reduce to genes in data
geneSet.all.expression= burtoni.snseq.combined.sct.neurons@assays$SCT@counts %>%
  as.data.frame() %>% 
    rownames_to_column('Gene') %>% 
  filter(Gene %in% c(geneSet.all,
                     'avp',
                     'oxt')) %>% 
  column_to_rownames('Gene') 

# set all values below 2 to 0
geneSet.all.expression[geneSet.all.expression < 2] <- 0

# set all values above 1 to 1
geneSet.all.expression[geneSet.all.expression > 1] <- 1

# 
geneSet.all.expression = t(geneSet.all.expression)

# check reduced gene list
geneSet.all.expression.count = geneSet.all.expression %>%  
  colSums() %>%
  as.data.frame() %>%
  dplyr::rename("Cell.count" = ".") %>% 
  rownames_to_column('Gene') 
  
# graph
geneSet.all.expression.count %>% 
  ggplot() +
  geom_point(aes(x = reorder(Gene, -Cell.count),
               y = Cell.count)) +
  theme_classic()  +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1)) +
  xlab('') +
  geom_hline(yintercept = 50) +
  geom_hline(yintercept = 500) 
ggsave('neuropeptides/avp.oxt/avp.vs.oxt/Hypomap C66-22 cell count.png',
       height = 5,
       width = 5)

# remove low expressing genes
# remove high expressing genes
geneSet.reduce = geneSet.all.expression.count %>% 
  filter(Cell.count > 50) %>% 
  filter(Cell.count < 500) %>% 
  pull(Gene)


# cluster 17 markers
tmp = PrepSCTFindMarkers(burtoni.snseq.combined.sct.neurons)
cluster17.markers.df <- FindMarkers(tmp, ident.1 = 17, min.pct = 0.25, assay = 'SCT')

cluster17.markers = cluster17.markers.df %>% 
  rownames_to_column('Gene')

cluster17.markers = cluster17.markers %>% 
  mutate(specificity = avg_log2FC*(pct.1/pct.2))

cluster17.markers = cluster17.markers %>% 
  filter(specificity > 0) %>% 
  filter(p_val_adj < 0.1) 

# get overlap of markers
geneSet.overlap = cluster17.markers %>% 
  filter(Gene %in%
           intersect(cluster17.markers %>% 
                       pull(Gene),
                     geneSet.all)) %>% 
  pull(Gene)

## remove avp and oxt
geneSet.all.remove = geneSet.all[geneSet.all != c('avp', 'oxt')]
geneSet.reduce.remove  = geneSet.reduce[geneSet.reduce != c('avp', 'oxt')]
geneSet.17 = cluster17.markers$Gene
geneSet.17.remove = geneSet.17[geneSet.17 != c('avp', 'oxt')]

#make dummy dataframe
tmp = burtoni.snseq.combined.sct.neurons

## all marker genes
tmp <- AddModuleScore(object = burtoni.snseq.combined.sct.neurons, 
                      features = list(geneSet.all.remove),
                      name = "geneSet.all.remove",
                      assay = 'SCT')

# #graph umap
FeaturePlot(object = tmp, features = "geneSet.all.remove1",
            max.cutoff = 1) +
  theme_classic()
ggsave('neuropeptides/avp.oxt/avp.vs.oxt/Hypomap C66-22 all UMAP.png',
       width = 5,
       height = 5)

# graph dotplot
DotPlot(object = tmp, 
        features = "geneSet.all.remove1")+
  theme_classic()
ggsave('neuropeptides/avp.oxt/avp.vs.oxt/Hypomap C66-22 all dotplot.png',
       width = 5,
       height = 5)

DotPlot(object = tmp, features = c(geneSet.all.remove)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1))
ggsave('neuropeptides/avp.oxt/avp.vs.oxt/Hypomap C66-22 all genes dotplot.png',
       width = 10,
       height = 10)

## all marker genes
# add oxt and avp
tmp <- AddModuleScore(object = burtoni.snseq.combined.sct.neurons, 
                      features = list(geneSet.all),
                      name = "geneSet.all",
                      assay = 'SCT')

# #graph umap
FeaturePlot(object = tmp, features = "geneSet.all1")+
  theme_classic()
ggsave('neuropeptides/avp.oxt/avp.vs.oxt/Hypomap C66-22 all UMAP avp and oxt.png',
       width = 5,
       height = 5)

# graph dotplot
DotPlot(object = tmp, 
        features = "geneSet.all1")+
  theme_classic()
ggsave('neuropeptides/avp.oxt/avp.vs.oxt/Hypomap C66-22 all dotplot avp and oxt.png',
       width = 5,
       height = 5)

DotPlot(object = tmp, features = c(geneSet.all)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1))
ggsave('neuropeptides/avp.oxt/avp.vs.oxt/Hypomap C66-22 all genes dotplot avp and oxt.png',
       width = 5,
       height = 5)

## reduce marker genes
tmp <- AddModuleScore(object = burtoni.snseq.combined.sct.neurons, 
                      features = list(geneSet.reduce.remove), 
                      name = "geneSet.reduce.remove",
                      assay = 'SCT')

# #graph umap
FeaturePlot(object = tmp, features = "geneSet.reduce.remove1")+
  theme_classic()
ggsave('neuropeptides/avp.oxt/avp.vs.oxt/Hypomap C66-22 reduce UMAP.png',
       width = 5,
       height = 5)

# graph dotplot
DotPlot(object = tmp, 
        features = "geneSet.reduce.remove1")+
  theme_classic()
ggsave('neuropeptides/avp.oxt/avp.vs.oxt/Hypomap C66-22 reduce dotplot.png',
       width = 5,
       height = 5)

DotPlot(object = tmp, features = c(geneSet.reduce.remove)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1))
ggsave('neuropeptides/avp.oxt/avp.vs.oxt/Hypomap C66-22 reduce genes dotplot.png',
       width = 5,
       height = 5)


## reduce marker genes
# add oxt and avp
tmp <- AddModuleScore(object = burtoni.snseq.combined.sct.neurons, 
                      features = list(geneSet.reduce), 
                      name = "geneSet.reduce",
                      assay = 'SCT')

# #graph umap
FeaturePlot(object = tmp, features = "geneSet.reduce1")+
  theme_classic()
ggsave('neuropeptides/avp.oxt/avp.vs.oxt/Hypomap C66-22 reduce UMAP avp and oxt.png',
       width = 5,
       height = 5)

# graph dotplot
DotPlot(object = tmp, 
        features = "geneSet.reduce1")+
  theme_classic()
ggsave('neuropeptides/avp.oxt/avp.vs.oxt/Hypomap C66-22 reduce dotplot avp and oxt.png',
       width = 5,
       height = 5)

DotPlot(object = tmp, features = geneSet.reduce)+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1))
ggsave('neuropeptides/avp.oxt/avp.vs.oxt/Hypomap C66-22 reduce genes dotplot avp and oxt.png',
       width = 5,
       height = 5)



# compare overlap genes
tmp <- AddModuleScore(object = burtoni.snseq.combined.sct.neurons, 
                      features = list(geneSet.overlap), 
                      name = "geneSet",
                      assay = 'SCT')
 
# #graph umap
FeaturePlot(object = tmp, features = "geneSet1")+
  theme_classic()
ggsave('neuropeptides/avp.oxt/avp.vs.oxt/Cluster 17 and Hypomap C66-22 UMAP dotplot.png',
       width = 5,
       height = 5)

# graph dotplot
DotPlot(object = tmp, features = "geneSet1")+
  theme_classic()
ggsave('neuropeptides/avp.oxt/avp.vs.oxt/Cluster 17 and Hypomap C66-22 overlap dotplot.png',
       width = 5,
       height = 5)

DotPlot(object = tmp, features = geneSet.overlap)+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1))
ggsave('neuropeptides/avp.oxt/avp.vs.oxt/Cluster 17 and Hypomap C66-22 overlap genes dotplot.png',
       width = 5,
       height = 5)

## add avp and oxt
tmp <- AddModuleScore(object = burtoni.snseq.combined.sct.neurons, 
                      features = list(c(geneSet.overlap,
                                   'avp',
                                   'oxt')),
                      name = "geneSet.avp.oxt",
                      assay = 'SCT')

# #graph umap
FeaturePlot(object = tmp, features = "geneSet.avp.oxt1")+
  theme_classic()
ggsave('neuropeptides/avp.oxt/avp.vs.oxt/Cluster 17 and Hypomap C66-22 overlap UMAP avp and oxt.png',
       width = 5,
       height = 5)

# graph dotplot
DotPlot(object = tmp, features = "geneSet.avp.oxt1")+
  theme_classic()
ggsave('neuropeptides/avp.oxt/avp.vs.oxt/Cluster 17 and Hypomap C66-22 overlap dotplot avp and oxt.png',
       width = 5,
       height = 5)

DotPlot(object = tmp, features = c(geneSet.overlap,'avp','oxt'))+
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1))
ggsave('neuropeptides/avp.oxt/avp.vs.oxt/Cluster 17 and Hypomap C66-22 overlap genes dotplot avp and oxt.png',
       width = 5,
       height = 5)

## 17 marker genes
tmp <- AddModuleScore(object = burtoni.snseq.combined.sct.neurons, 
                      features = list(geneSet.17),
                      name = "geneSet.17.",
                      assay = 'SCT')

#graph umap
FeaturePlot(object = tmp, features = "geneSet.17.1")+
  theme_classic()
ggsave('neuropeptides/avp.oxt/avp.vs.oxt/Cluster 17 all UMAP.png',
       width = 5,
       height = 5)

# graph dotplot
DotPlot(object = tmp, 
        features = "geneSet.17.1")+
  theme_classic()
ggsave('neuropeptides/avp.oxt/avp.vs.oxt/Cluster 17 all dotplot.png',
       width = 5,
       height = 5)

DotPlot(object = tmp, features = c(geneSet.17)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1))
ggsave('neuropeptides/avp.oxt/avp.vs.oxt/Cluster 17 all genes dotplot.png',
       width = 5,
       height = 5)











#### All cluster identification ####
### get marker orthologs from hypomap C66-22 cluster
# use Hypomap supplemental material "5_hypomap_marker_genes"
# filter out cluster C66-22 from 'cluster_id' and make new excel csv
C66.data = read.csv('../Gene.lists/Cluster C66 1 to 50 hypomap.csv')

C66.data = C66.data %>% 
  dplyr::rename(Mouse_gene = gene)

## get brain region information
hypomap.region.data = read.csv('../Gene.lists/Clustser C286 brain annotation hypomap.csv')
# rename column
hypomap.region.data = hypomap.region.data %>% 
  dplyr::rename(cluster_id_286 = cluster_id)

## get cluster hierarchy information
C286.data = read.csv('../Gene.lists/Clustser C286 to C185 hypomap.csv')
C185.data = read.csv('../Gene.lists/Clustser C185 to C66 hypomap.csv')

## combine data
C286.data.C185.data = full_join(C286.data,
                                C185.data)
# add to brain data
hypomap.region.data = hypomap.region.data %>% 
  full_join(C286.data.C185.data %>% 
              dplyr::select(cluster_id_66,
                            cluster_id_286))

# graph overlap of clusters and brain regions
hypomap.region.data %>% 
  dplyr::select(cluster_id_66,
                Region_summarized,
                ncells) %>% 
  na.omit() %>%
  group_by(cluster_id_66,
           Region_summarized) %>% 
  summarise(ncells = sum(ncells)) %>% 
  ggplot(aes(x = cluster_id_66,
             y = Region_summarized,
             fill = log(ncells),
             label = ncells)) +
  geom_tile() +
  geom_text() +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1))
ggsave('neuropeptides/avp.oxt/avp.vs.oxt/HypomapMarkers/Hypomap C66 clusters overlap brain regions.png',
       height = 15,
       width = 30)


# reduce data frame and combine brain regions
hypomap.region.data.reduce = hypomap.region.data %>% 
  dplyr::select(cluster_id_66,
                Region_summarized) %>% 
  distinct() %>% 
  na.omit() %>% 
  group_by(cluster_id_66) %>% 
  summarise(Region_summarized=toString(unique(Region_summarized), .groups='drop')) %>% 
  dplyr::rename(cluster_id = cluster_id_66)

#### use ensembl to get orthologs
## load ensembl
library(biomaRt)
#selecting biomart database
ensembl = useMart('ensembl')
# #list datasets
# dataset = listDatasets(ensembl)
##select dataset
#mouse
ensembl.mouse = useDataset('mmusculus_gene_ensembl',
                           mart=ensembl)

# # get list of all attributes
# listAttributes(ensembl.mouse) %>% View


#create attributes lists
tilapia.attributes.gene = c('oniloticus_homolog_associated_gene_name',
                       'external_gene_name')

tilapia.attributes.id = c('oniloticus_homolog_ensembl_gene',
                       'external_gene_name')

# use gene names
tilapia.C66.gene = getBM(attributes = tilapia.attributes.gene,
                            mart = ensembl.mouse,
                            values = C66.data$Mouse_gene,
                            filter = 'external_gene_name',
                            useCache = FALSE) # useCache has to do with version of R not being up to date?

# use gene ID
tilapia.C66.ID = getBM(attributes = tilapia.attributes.id,
                            mart = ensembl.mouse,
                            values = C66.data$Mouse_gene,
                            filter = 'external_gene_name',
                            useCache = FALSE) # useCache has to do with version of R not being up to date?

# combine gene and ensembl gene IDs
tilapia.C66.gene.list = tilapia.C66.gene %>% 
  filter(oniloticus_homolog_associated_gene_name != "") %>%
  dplyr::rename(gene = oniloticus_homolog_associated_gene_name) %>% 
  rbind(tilapia.C66.ID %>% 
          filter(oniloticus_homolog_ensembl_gene != "")%>%
          dplyr::rename(gene = oniloticus_homolog_ensembl_gene) ) %>% 
  distinct() %>% 
  dplyr::rename(Mouse_gene = external_gene_name)

# combine with specificity data
tilapia.C66.gene.list = tilapia.C66.gene.list %>% 
  full_join(C66.data) %>% 
  distinct()

# 4108 genes
nrow(tilapia.C66.gene.list)

# remove NA
tilapia.C66.gene.list = tilapia.C66.gene.list %>% 
  na.omit()

# 3755 genes
nrow(tilapia.C66.gene.list)

# remove genes not in dataset
tilapia.C66.gene.list = tilapia.C66.gene.list %>% 
  filter(gene %in% rownames(burtoni.snseq.combined.sct.neurons@assays$SCT@counts))

# 1380 genes
nrow(tilapia.C66.gene.list)

# graph number of genes per cluster
tilapia.C66.gene.list %>% 
  dplyr::select(cluster_id) %>% 
  table() %>% 
  data.frame() %>% 
ggplot(aes(x = reorder(cluster_id, 
                       -Freq),
           y = Freq)) +
  geom_bar(stat = 'identity') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1)) +
  xlab('Hypomap Cluster_ID') +
  ylab('Tilapia gene count')
ggsave('neuropeptides/avp.oxt/avp.vs.oxt/HypomapMarkers/Hypomap C66 clusters gene count.png',
       height = 5,
       width = 10)

### create list of marker genes for each AVP cell type
AVP.clusters.to.neuron.clusters = burtoni.snseq.combined.sct.reduce.avp.neurons.recluster@meta.data %>% 
  dplyr::select(neurons_integrated_snn_res.0.8,
                neurons.avp_integrated_snn_res.0.8)

#create matrix
AVP.clusters.to.neuron.clusters.mat = AVP.clusters.to.neuron.clusters %>% 
  droplevels() %>% 
  table() %>% 
  as.matrix() %>% 
  t()

# convert to percentage
AVP.clusters.to.neuron.clusters.mat.percent = 
  round((AVP.clusters.to.neuron.clusters.mat/rowSums(AVP.clusters.to.neuron.clusters.mat))*100,2)

# graph count heatmap
AVP.clusters.to.neuron.clusters.mat %>%  
  pheatmap(display_numbers = T,
           number_format = "%.f",
           color = colorRampPalette(c("white","red"))(100),
           cluster_cols = F,
           cluster_rows = F)
ggsave('neuropeptides/avp.oxt/avp.vs.oxt/AVP neuron clusters to neuron clusters.png',
       height = 10,
       width = 10)

# graph percentage heatmap
AVP.clusters.to.neuron.clusters.mat.percent %>%  
  pheatmap(display_numbers = T,
           number_format = "%.f",
           color = colorRampPalette(c("white","red"))(100),
           cluster_cols = F,
           cluster_rows = F) 
ggsave('neuropeptides/avp.oxt/avp.vs.oxt/AVP neuron clusters to neuron clusters percent.png',
       height = 10,
       width = 10)

## AVP neuron cluster to neuron cluster
# filter out to keep only above 10%
AVP.clusters.to.neuron.clusters.df = AVP.clusters.to.neuron.clusters.mat.percent %>% 
  as.data.frame() %>% 
  filter(Freq > 10)

# graph
AVP.clusters.to.neuron.clusters.mat.percent %>% 
  as.data.frame() %>% 
  ggplot(aes(x = neurons.avp_integrated_snn_res.0.8,
             label = neurons_integrated_snn_res.0.8,
             y = Freq)) +
  geom_label() +
  theme_classic() +
  geom_hline(yintercept = 25)
ggsave('neuropeptides/avp.oxt/avp.vs.oxt/AVP neuron clusters to neuron clusters.png',
                                    height = 10,
                                    width = 10)

### get cluster markers
## keep positive markers for cluster (specificity above 0)
## use weighted 2x FC and 1.5 ratio of in cluster/out cluster (specificity > 0.87) [ 0.87 = 0.58 * 1.5]
## only significant markers (adjusted pvalue < 0.05)

# prepare data
burtoni.snseq.combined.sct.neurons = PrepSCTFindMarkers(burtoni.snseq.combined.sct.neurons)

## AVP cluster 0
avp.cluster0.markers.df <- FindMarkers(burtoni.snseq.combined.sct.neurons, 
                                       ident.1 = 0, 
                                       min.pct = 0.25, 
                                       assay = 'SCT')
#create gene column
avp.cluster0.markers.df = avp.cluster0.markers.df %>% 
  rownames_to_column('gene')
# create specificity score
avp.cluster0.markers.df = avp.cluster0.markers.df %>% 
  mutate(specificity = avg_log2FC*(pct.1/pct.2))
# filter out poor markers
avp.cluster0.markers.df = avp.cluster0.markers.df %>% 
  filter(specificity > 0.87) %>% 
  filter(p_val_adj < 0.05) 
# add cluster id
avp.cluster0.markers.df = avp.cluster0.markers.df %>% 
  mutate(Avp.cluster.id = "0")

## AVP cluster 1
avp.cluster1.markers.df <- FindMarkers(burtoni.snseq.combined.sct.neurons, 
                                       ident.1 = 1, 
                                       min.pct = 0.25, 
                                       assay = 'SCT')
#create gene column
avp.cluster1.markers.df = avp.cluster1.markers.df %>% 
  rownames_to_column('gene')
# create specificity score
avp.cluster1.markers.df = avp.cluster1.markers.df %>% 
  mutate(specificity = avg_log2FC*(pct.1/pct.2))
# filter out poor markers
avp.cluster1.markers.df = avp.cluster1.markers.df %>% 
  filter(specificity > 0.87) %>% 
  filter(p_val_adj < 0.05) 
# add cluster id
avp.cluster1.markers.df = avp.cluster1.markers.df %>% 
  mutate(Avp.cluster.id = "1")

## AVP cluster 2
avp.cluster2.markers.df <- FindMarkers(burtoni.snseq.combined.sct.neurons, 
                                       ident.1 = c(3,5), 
                                       min.pct = 0.25, 
                                       assay = 'SCT')
#create gene column
avp.cluster2.markers.df = avp.cluster2.markers.df %>% 
  rownames_to_column('gene')
# create specificity score
avp.cluster2.markers.df = avp.cluster2.markers.df %>% 
  mutate(specificity = avg_log2FC*(pct.1/pct.2))
# filter out poor markers
avp.cluster2.markers.df = avp.cluster2.markers.df %>% 
  filter(specificity > 0.87) %>% 
  filter(p_val_adj < 0.05) 
# add cluster id
avp.cluster2.markers.df = avp.cluster2.markers.df %>% 
  mutate(Avp.cluster.id = "2")

## AVP cluster 3
avp.cluster3.markers.df <- FindMarkers(burtoni.snseq.combined.sct.neurons, 
                                       ident.1 = 2, 
                                       min.pct = 0.25, 
                                       assay = 'SCT')
#create gene column
avp.cluster3.markers.df = avp.cluster3.markers.df %>% 
  rownames_to_column('gene')
# create specificity score
avp.cluster3.markers.df = avp.cluster3.markers.df %>% 
  mutate(specificity = avg_log2FC*(pct.1/pct.2))
# filter out poor markers
avp.cluster3.markers.df = avp.cluster3.markers.df %>% 
  filter(specificity > 0.87) %>% 
  filter(p_val_adj < 0.05) 
# add cluster id
avp.cluster3.markers.df = avp.cluster3.markers.df %>% 
  mutate(Avp.cluster.id = "3")

## AVP cluster 4
avp.cluster4.markers.df <- FindMarkers(burtoni.snseq.combined.sct.neurons, 
                                       ident.1 = 11, 
                                       min.pct = 0.25, 
                                       assay = 'SCT')
#create gene column
avp.cluster4.markers.df = avp.cluster4.markers.df %>% 
  rownames_to_column('gene')
# create specificity score
avp.cluster4.markers.df = avp.cluster4.markers.df %>% 
  mutate(specificity = avg_log2FC*(pct.1/pct.2))
# filter out poor markers
avp.cluster4.markers.df = avp.cluster4.markers.df %>% 
  filter(specificity > 0.87) %>% 
  filter(p_val_adj < 0.05) 
# add cluster id
avp.cluster4.markers.df = avp.cluster4.markers.df %>% 
  mutate(Avp.cluster.id = "4")

### combine markers together
avp.cluster.markers.df = avp.cluster0.markers.df %>% 
  rbind(avp.cluster1.markers.df) %>% 
  rbind(avp.cluster2.markers.df) %>% 
  rbind(avp.cluster3.markers.df) %>% 
  rbind(avp.cluster4.markers.df) 

# graph marker counts
avp.cluster.markers.df %>% 
  dplyr::select(Avp.cluster.id) %>% 
  table() %>% 
  data.frame() %>%
  # add_row(Avp.cluster.id = "0", 
  #         Freq = 0) %>% 
  ggplot(aes(x = Avp.cluster.id,
             y = Freq)) +
  geom_bar(stat = 'identity') +
  geom_label(aes(label = Freq)) +
  theme_classic() +
  ylab('Number of marker genes')
ggsave('neuropeptides/avp.oxt/avp.vs.oxt/HypomapMarkers/AVP neuron clusters marker gene count.png',
       height = 5,
       width = 5)

# compare count to ensembl list
avp.cluster.markers.df %>% 
  filter(gene %in% tilapia.C66.gene.list$gene) %>% 
  dplyr::select(Avp.cluster.id) %>% 
  table() %>% 
  data.frame() %>%
  # add_row(Avp.cluster.id = "0", 
  #         Freq = 0) %>% 
  mutate(Mouse_gene_data = "Present") %>% 
  rbind(
avp.cluster.markers.df %>% 
  dplyr::select(Avp.cluster.id) %>% 
  table() %>% 
  data.frame() %>%
  # add_row(Avp.cluster.id = "0", 
  #         Freq = 0) %>% 
  mutate(Mouse_gene_data = "Absent")) %>% 
  ggplot(aes(x = Avp.cluster.id,
             y = Freq)) +
  geom_bar(stat = 'identity',
           aes(fill = Mouse_gene_data)) +
  geom_label(aes(label = Freq)) +
  theme_classic() +
  ylab('Number of marker genes')
ggsave('neuropeptides/avp.oxt/avp.vs.oxt/HypomapMarkers/AVP neuron clusters marker gene count mouse data.png',
       height = 5,
       width = 5)

## for paper
avp.cluster.markers.df %>% 
  filter(gene %in% tilapia.C66.gene.list$gene) %>% 
  dplyr::select(Avp.cluster.id) %>% 
  table() %>% 
  data.frame() %>%
  # add_row(Avp.cluster.id = "0", 
  #         Freq = 0) %>% 
  mutate(Mouse_gene_data = "Present") %>% 
  rbind(
    avp.cluster.markers.df %>% 
      dplyr::select(Avp.cluster.id) %>% 
      table() %>% 
      data.frame() %>%
      # add_row(Avp.cluster.id = "0", 
      #         Freq = 0) %>% 
      mutate(Mouse_gene_data = "Absent")) %>% 
  ggplot(aes(x = Avp.cluster.id,
             y = Freq)) +
  geom_bar(stat = 'identity',
           aes(fill = Mouse_gene_data)) +
  geom_label(aes(label = Freq)) +
  scale_fill_manual(values = c('grey75',
                               'grey25')) +
  theme_classic() +
  ylab('Number of marker genes')+
  xlab('AVP neuron cluster') +
  theme(legend.position = 'none')
ggsave('neuropeptides/avp.oxt/avp.vs.oxt/HypomapMarkers/AVP neuron clusters marker gene count mouse data paper.png',
       height = 3.5,
       width = 3.5)

# filter out low expressing genes
avp.cluster.markers.df.reduce = avp.cluster.markers.df 

# identify duplicated genes in AVP markers present in mouse data
avp.cluster.markers.df.reduce %>% 
  filter(n() > 1, .by = gene) %>% 
  filter(gene %in% tilapia.C66.gene.list$gene)

# identify duplicated genes in mouse data
tilapia.C66.gene.list %>% 
  filter(n() > 1, .by = gene) %>% 
  filter(gene %in% tilapia.C66.gene.list$gene)

## combine marker genes with mouse data
avp.cluster.markers.df.reduce.C66 = tilapia.C66.gene.list %>% 
  left_join(avp.cluster.markers.df.reduce,
            by = 'gene',
            suffix = c("", ".avp")) %>% 
  mutate(Avp.cluster.id = ifelse(is.na(Avp.cluster.id),
                                 'None',
                                 Avp.cluster.id))

## sum number of genes present for each cluster
avp.cluster.markers.df.reduce.C66.sum = avp.cluster.markers.df.reduce.C66 %>% 
  mutate(Total = 1,
         Avp.count = ifelse(Avp.cluster.id == 'None',
                            0,
                            1)) %>% 
  group_by(cluster_id) %>% 
  mutate(Total.sum = sum(Total)) %>% 
  ungroup() %>% 
  group_by(cluster_id,
           Avp.cluster.id,
           Total.sum) %>% 
  summarise(Avp.count.sum = sum(Avp.count),
            mouse.specificty.sum = sum(specificity),
            Avp.specificty.sum = sum(specificity.avp)) %>% 
  ungroup() %>% 
  mutate(specificity.prod.avg = mouse.specificty.sum * Avp.specificty.sum)

# add brain data
avp.cluster.markers.df.reduce.C66.sum = avp.cluster.markers.df.reduce.C66.sum %>% 
  full_join(hypomap.region.data.reduce)


#### save all cluster results ####
write.csv(avp.cluster.markers.df.reduce.C66.sum %>% 
            filter(Avp.cluster.id != 'None'),
          file = 'neuropeptides/avp.oxt/avp.vs.oxt/HypomapMarkers/avp.cluster.markers.df.reduce.C66.sum.reduce.csv',
          row.names = F)

write.csv(avp.cluster.markers.df.reduce.C66.sum,
          file = 'neuropeptides/avp.oxt/avp.vs.oxt/HypomapMarkers/avp.cluster.markers.df.reduce.C66.sum.csv',
          row.names = F)

# avp.cluster.markers.df.reduce.C66.sum = read.csv('neuropeptides/avp.oxt/avp.cluster.markers.df.reduce.C66.sum.csv')

#### graph Hypomap AVP cluster identification continue ####
## Make summary table 
avp.cluster.markers.df.reduce.C66.sum.all = hypomap.region.data %>% 
  dplyr::select(cluster_id_66,
                Region_summarized,
                ncells) %>% 
  na.omit() %>%
  group_by(cluster_id_66) %>% 
  mutate(ncells.prop = ncells/sum(ncells)) %>% 
  ungroup() %>% 
  dplyr::rename(cluster_id = cluster_id_66) %>% 
  dplyr::rename(Region = Region_summarized) %>% 
  right_join(avp.cluster.markers.df.reduce.C66.sum) %>% 
  mutate(specificity.prod.avg.weight = specificity.prod.avg * ncells.prop) %>% 
  droplevels() %>% 
  na.omit() 

# region
avp.cluster.markers.df.reduce.C66.sum.all.region = avp.cluster.markers.df.reduce.C66.sum.all %>% 
  group_by(Avp.cluster.id,
           Region) %>% 
  summarise(specificity.prod.avg.weight.sum = sum(specificity.prod.avg.weight)) %>% 
  na.omit() 

# region
avp.cluster.markers.df.reduce.C66.sum.all.cluster = avp.cluster.markers.df.reduce.C66.sum.all %>% 
  group_by(Avp.cluster.id,
           cluster_id) %>% 
  summarise(specificity.prod.avg.weight.sum = sum(specificity.prod.avg.weight)) %>% 
  na.omit() 
  
## graph
# scaled by proportion of cells
# region
for (i in unique(avp.cluster.markers.df.reduce.C66.sum.all$Avp.cluster.id)) {
  avp.cluster.markers.df.reduce.C66.sum.all.region %>% 
    filter(Avp.cluster.id == i) %>% 
  ggplot(aes(x = reorder(Region,
                         -specificity.prod.avg.weight.sum),
             y = specificity.prod.avg.weight.sum,
             label = round(specificity.prod.avg.weight.sum))) +
  geom_bar(stat = 'identity') +
  geom_label() +
  theme_classic() +
    xlab('') +
    ggtitle(paste('AVP neuron cluster: ',
                  i,
                  sep = '')) +
    theme(axis.text.x = element_text(angle = 45, 
                                     hjust = 1))
  ggsave(paste('neuropeptides/avp.oxt/avp.vs.oxt/HypomapMarkers/Brain region Avp cluster ',
               i,
               '.png',
               sep = ''),
         height = 10,
         width = 10)
}

# cluster
for (i in unique(avp.cluster.markers.df.reduce.C66.sum.all$Avp.cluster.id)) {
  avp.cluster.markers.df.reduce.C66.sum.all.cluster %>% 
    filter(Avp.cluster.id == i) %>% 
    ggplot(aes(x = reorder(cluster_id,
                           -specificity.prod.avg.weight.sum),
               y = specificity.prod.avg.weight.sum,
               label = round(specificity.prod.avg.weight.sum))) +
    geom_bar(stat = 'identity') +
    geom_label() +
    theme_classic() +
    xlab('') +
    ggtitle(paste('AVP neuron cluster: ',
                  i,
                  sep = '')) +
    theme(axis.text.x = element_text(angle = 45, 
                                     hjust = 1))
  ggsave(paste('neuropeptides/avp.oxt/avp.vs.oxt/HypomapMarkers/Brain cluster Avp cluster ',
               i,
               '.png',
               sep = ''),
         height = 10,
         width = 10)
}

## graph heatmap
# brain regions
avp.cluster.markers.df.reduce.C66.sum.all.region %>% 
  filter(specificity.prod.avg.weight.sum >= 1) %>% 
  ggplot(aes(x = reorder(Region,
                         -specificity.prod.avg.weight.sum),
             y = Avp.cluster.id,
             label = round(specificity.prod.avg.weight.sum),
             fill = log(specificity.prod.avg.weight.sum))) +
  geom_tile() +
  geom_text(color='black')+
  scale_fill_gradient2(low = 'white',
                       # mid = 'red1',
                       high = 'red4') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1),
        legend.position = 'none')
ggsave('neuropeptides/avp.oxt/avp.vs.oxt/HypomapMarkers/Heatmap Brain region Avp cluster.png',
       height = 6.5,
       width = 13)

# for paper
avp.cluster.markers.df.reduce.C66.sum.all.region %>% 
  filter(specificity.prod.avg.weight.sum >= 1) %>% 
  ggplot(aes(x = reorder(Region,
                         -specificity.prod.avg.weight.sum),
             y = Avp.cluster.id,
             label = round(specificity.prod.avg.weight.sum),
             fill = log(specificity.prod.avg.weight.sum))) +
  geom_tile() +
  geom_text(color='black')+
  scale_fill_gradient2(low = 'white',
                       # mid = 'red1',
                       high = 'red4') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1),
        legend.position = 'none',
        text = element_text(size = 8)) +
  xlab('') +
  ylab('AVP neuron cluster')
ggsave('neuropeptides/avp.oxt/avp.vs.oxt/HypomapMarkers/Heatmap Brain region Avp cluster paper.png',
       height = 6.5,
       width = 13)

# brain clusters
avp.cluster.markers.df.reduce.C66.sum.all.cluster %>% 
  filter(specificity.prod.avg.weight.sum >= 1) %>% 
  ggplot(aes(x = cluster_id,
             y = Avp.cluster.id,
             label = round(specificity.prod.avg.weight.sum),
             fill = log(specificity.prod.avg.weight.sum))) +
  geom_tile() +
  geom_text(color='black')+
  scale_fill_gradient2(low = 'white',
                       # mid = 'red1',
                       high = 'red4') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1),
        legend.position = 'none')
ggsave('neuropeptides/avp.oxt/avp.vs.oxt/HypomapMarkers/Heatmap Brain cluster Avp cluster.png',
       height = 6.5,
       width = 13)















#### limmatrend genotype, clusters, and status approach with SCT treat prep ####
#### Variable have same name as 'limmatrend genotype, clusters, and status approach'
### calculate variable genes
## identify top 500 variable genes
# use integrated assay for variable features
burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.variable.prep <- FindVariableFeatures(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster, 
                                                                                              assay = 'integrated',
                                                                                         selection.method = "vst", 
                                                                                         nfeatures = 500, 
                                                                                         verbose = F)
# identify top 500 variable genes
avp.neuron.reduce.group.topgenes.prep <- head(VariableFeatures(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.variable.prep,
                                                               assay = 'integrated'), 
                                         500)

# create dummy
burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.expression.sct.prep = burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.expression

## create vector of factor
avp.neuron.reduce.DomvsSub.vector.list.sct.prep = burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.expression.sct.prep %>% 
  mutate(avp.neuron.orig.ident = orig.ident %>% 
           as.factor()) %>% 
  pull(avp.neuron.orig.ident) %>% 
  droplevels()

# create dummy
burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.prep = burtoni.snseq.combined.sct.reduce.avp.neurons.recluster

### counts matrix
## raw read count matrix
## rows = genes, columns = cells
# only keep 500 variable genes
# set negative values to 0
avp.neuron.reduce.DomvsSub.vector.count.sct.prep = GetAssayData(burtoni.snseq.combined.sct.reduce.avp.neurons.recluster.prep,
                                                           assay = 'SCT') %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column('gene') %>% 
  dplyr::select(c(gene,
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





#### poster heatmap ####
# pheatmap
avp.neuron.reduce.genotype.vector.limma.results.df.sct.treat.prep  %>%
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
  arrange('4',
          '3',
          '2',
          '1',
          '0',
          All) %>%
  as.matrix() %>%
  t() %>%
  pheatmap(cluster_rows = F,
           cluster_cols = T,
           scale = 'none',
           border_color = 'black',
           color = colorRampPalette(c("#60bb46", "white","#4e499e"))(100),
           legend = F,
           treeheight_col = 0,
           treeheight_row = 0,
           angle_col = 315,
           fontsize = 15,
           breaks = seq(-4, 4, length.out = 100),
           filename = "neuropeptides/avp.oxt/avp/limmatrend/Heatmap DEGs across clusters update.pdf",
           width = 9.71,
           height = 6.5
  )
dev.off()

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
                          Gene == "ENSONIG00000029183" ~ 'pro-MCH 2',
                          Gene == "ENSONIG00000036400" ~ '36400',
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
                                              order = rbind(avp.neuron.deg.order[c(24:25,16:23),],avp.neuron.deg.order[c(8:15,3:7,26:29,1:2),]) %>% pull(label) )
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







  
  












