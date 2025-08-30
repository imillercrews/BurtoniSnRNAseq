#### Burtoni snseq seurat analysis
### hdWGCNA neurons
###Note: Seurat requires R version > 4
## use lambcomp1 to run R with command 
# > R-4.0.3
###hdWGCNA
##https://smorabit.github.io/hdWGCNA/


### set working directory
setwd("/stor/work/Hofmann/All_projects/A_burtoni_snseq/seurat/")

#### installation ####

### hdwgcna
# # OLD
# # create new conda environment for R
# conda create -n hdWGCNA -c conda-forge r-base r-essentials
# 
# # activate conda environment
# conda activate hdWGCNA

# install BiocManager
# install.packages("BiocManager")

# install Bioconductor core packages
# BiocManager::install()

# # install additional packages:
# install.packages(c("Seurat", "WGCNA", "igraph", "devtools"))

#Now you can install the hdWGCNA package using devtools.
# devtools::install_github('smorabit/hdWGCNA', ref='dev')

#load gene overlap
# BiocManager::install("GeneOverlap")

#### load libraries ####

#load libraries
library(Seurat)
library(GeneOverlap)
library(igraph)
library(WGCNA)
library(hdWGCNA)
library(cowplot)
library(patchwork)
library(corrplot)
library(ggrepel)
library(tidyverse)
library(emmeans)
library(multcomp)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

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

# ### Add metadata old
# ## cell type
# burtoni.snseq.combined.sct.all = AddMetaData(
#   object = burtoni.snseq.combined.sct.all,
#   metadata = burtoni.scsorter.data.scsort.output %>% 
#     select(Cell.type,
#            Cell.id) %>% 
#     column_to_rownames(var = "Cell.id"),
#   col.name = 'Cell.type'
# )
# 
# ## genotype
# burtoni.snseq.combined.sct.all = AddMetaData(
#   object = burtoni.snseq.combined.sct.all,
#   metadata = burtoni.souporcell.filtered %>% 
#     select(Genotype.id,
#            Cell.id) %>% 
#     column_to_rownames(var = "Cell.id"),
#   col.name = 'Genotype.id'
# )
# 
# ## sctype.hypo
# burtoni.snseq.combined.sct.all = AddMetaData(
#   object = burtoni.snseq.combined.sct.all,
#   metadata = burtoni.sctypemarkers.hypo %>% 
#     select(sctypemarkers.hypo,
#            Cell.id) %>% 
#     column_to_rownames(var = "Cell.id"),
#   col.name = 'sctypemarkers.hypo'
# )


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

### neuropeptide list
neuropeptides.list = read_csv("../Gene.lists/neuropeptides.list_orthologs.csv")

#reduced range 2
resolution.range.reduced.2 <- seq(from = 0.2, to = 1, by = 0.1)

#### functions ####
### create function to get dotplot data
DotPlot.data = function (object, assay = NULL, features, cols = c("lightgrey", 
                                                                  "blue"), col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, 
                         idents = NULL, group.by = NULL, split.by = NULL, cluster.idents = FALSE, 
                         scale = TRUE, scale.by = "radius", scale.min = NA, scale.max = NA) 
{
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  split.colors <- !is.null(x = split.by) && !any(cols %in% 
                                                   rownames(x = brewer.pal.info))
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
    if (split.colors) {
      if (length(x = unique(x = splits)) > length(x = cols)) {
        stop("Not enough colors for the number of groups")
      }
      cols <- cols[1:length(x = unique(x = splits))]
      names(x = cols) <- unique(x = splits)
    }
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
  if (split.colors) {
    avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, 
                                         breaks = 20))
  }
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(x = data.plot$features.plot, 
                                    levels = features)
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  if (split.colors) {
    splits.use <- vapply(X = as.character(x = data.plot$id), 
                         FUN = gsub, FUN.VALUE = character(length = 1L), pattern = paste0("^((", 
                                                                                          paste(sort(x = levels(x = object), decreasing = TRUE), 
                                                                                                collapse = "|"), ")_)"), replacement = "", 
                         USE.NAMES = FALSE)
    data.plot$colors <- mapply(FUN = function(color, value) {
      return(colorRampPalette(colors = c("grey", color))(20)[value])
    }, color = cols[splits.use], value = avg.exp.scaled)
  }
  color.by <- ifelse(test = split.colors, yes = "colors", no = "avg.exp.scaled")
  if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
  }
  if (!is.null(x = feature.groups)) {
    data.plot$feature.groups <- factor(x = feature.groups[data.plot$features.plot], 
                                       levels = unique(x = feature.groups))
  }
  return(data.plot)
}

### create function for bigger network graphs
ModuleUMAPPlot.size = function (seurat_obj, sample_edges = TRUE, edge_prop = 0.2, 
                                label_hubs = 5, edge.alpha = 0.25, vertex.label.cex = 0.5, 
                                label_genes = NULL, return_graph = FALSE, keep_grey_edges = TRUE, dot.size = 3, edge.size = 0.5,
                                wgcna_name = NULL, ...) 
{
  if (is.null(wgcna_name)) {
    wgcna_name <- seurat_obj@misc$active_wgcna
  }
  TOM <- GetTOM(seurat_obj, wgcna_name)
  modules <- GetModules(seurat_obj, wgcna_name)
  umap_df <- GetModuleUMAP(seurat_obj, wgcna_name)
  mods <- levels(umap_df$module)
  mods <- mods[mods != "grey"]
  subset_TOM <- TOM[umap_df$gene, umap_df$gene[umap_df$hub == 
                                                 "hub"]]
  hub_list <- lapply(mods, function(cur_mod) {
    cur <- subset(modules, module == cur_mod)
    cur[, c("gene_name", paste0("kME_", cur_mod))] %>% top_n(label_hubs) %>% 
      .$gene_name
  })
  names(hub_list) <- mods
  hub_labels <- as.character(unlist(hub_list))
  print("hub labels")
  print(hub_labels)
  print(label_genes)
  if (is.null(label_genes)) {
    label_genes <- hub_labels
  }
  else {
    if (!any(label_genes %in% umap_df$gene)) {
      stop("Some genes in label_genes not found in the UMAP.")
    }
    label_genes <- unique(c(label_genes, hub_labels))
  }
  print(label_genes)
  selected_modules <- modules[umap_df$gene, ]
  selected_modules <- cbind(selected_modules, umap_df[, c("UMAP1", 
                                                          "UMAP2", "hub", "kME")])
  selected_modules$label <- ifelse(selected_modules$gene_name %in% 
                                     label_genes, selected_modules$gene_name, "")
  selected_modules$fontcolor <- ifelse(selected_modules$color == 
                                         "black", "gray50", "black")
  selected_modules$framecolor <- ifelse(selected_modules$gene_name %in% 
                                          label_genes, "black", selected_modules$color)
  edge_df <- subset_TOM %>% reshape2::melt()
  print(dim(edge_df))
  edge_df$color <- future.apply::future_sapply(1:nrow(edge_df), 
                                               function(i) {
                                                 gene1 = as.character(edge_df[i, "Var1"])
                                                 gene2 = as.character(edge_df[i, "Var2"])
                                                 col1 <- selected_modules[selected_modules$gene_name == 
                                                                            gene1, "color"]
                                                 col2 <- selected_modules[selected_modules$gene_name == 
                                                                            gene2, "color"]
                                                 if (col1 == col2) {
                                                   col = col1
                                                 }
                                                 else {
                                                   col = "grey90"
                                                 }
                                                 col
                                               })
  if (!keep_grey_edges) {
    edge_df <- edge_df %>% subset(color != "grey90")
  }
  groups <- unique(edge_df$color)
  if (sample_edges) {
    temp <- do.call(rbind, lapply(groups, function(cur_group) {
      cur_df <- edge_df %>% subset(color == cur_group)
      n_edges <- nrow(cur_df)
      cur_sample <- sample(1:n_edges, round(n_edges * 
                                              edge_prop))
      cur_df[cur_sample, ]
    }))
  }
  else {
    temp <- do.call(rbind, lapply(groups, function(cur_group) {
      cur_df <- edge_df %>% subset(color == cur_group)
      n_edges <- nrow(cur_df)
      cur_df %>% dplyr::top_n(round(n_edges * edge_prop), 
                              wt = value)
    }))
  }
  edge_df <- temp
  print(dim(edge_df))
  edge_df <- edge_df %>% group_by(color) %>% mutate(value = scale01(value))
  edge_df <- edge_df %>% arrange(value)
  edge_df <- rbind(subset(edge_df, color == "grey90"), subset(edge_df, 
                                                              color != "grey90"))
  edge_df$color_alpha <- ifelse(edge_df$color == "grey90", 
                                alpha(edge_df$color, alpha = edge_df$value/2), alpha(edge_df$color, 
                                                                                     alpha = edge_df$value))
  selected_modules <- rbind(subset(selected_modules, hub == 
                                     "other"), subset(selected_modules, hub != "other"))
  selected_modules <- rbind(subset(selected_modules, label == 
                                     ""), subset(selected_modules, label != ""))
  g <- igraph::graph_from_data_frame(edge_df, directed = FALSE, 
                                     vertices = selected_modules)
  print("making net")
  print(head(edge_df))
  print(head(selected_modules))
  if (return_graph) {
    return(g)
  }
  plot(g, layout = as.matrix(selected_modules[, c("UMAP1", 
                                                  "UMAP2")]), edge.color = adjustcolor(E(g)$color_alpha, 
                                                                                       alpha.f = edge.alpha), vertex.size = V(g)$kME * dot.size, edge.curved = 0, 
       edge.width = edge.size, vertex.color = V(g)$color, vertex.label = V(g)$label, 
       vertex.label.dist = 1.1, vertex.label.degree = -pi/4, 
       vertex.label.family = "Helvetica", vertex.label.font = 3, 
       vertex.label.color = V(g)$fontcolor, vertex.label.cex = 0, 
       vertex.frame.color = V(g)$framecolor, margin = 0)
}

#### subset neuron data ####
burtoni.snseq.combined.sct.neurons = burtoni.snseq.combined.sct

# keep original clusters
burtoni.snseq.combined.sct.neurons@meta.data$seurat_clusters_all = burtoni.snseq.combined.sct.neurons@meta.data$seurat_clusters

#set idents
Idents(object = burtoni.snseq.combined.sct.neurons) <- "sctypemarkers.hypo"

#subset to neurons
burtoni.snseq.combined.sct.neurons = subset(burtoni.snseq.combined.sct.neurons,
                                                idents = c("C7-1: GLU",
                                                           "C7-2: GABA"))

#remove large file
rm(burtoni.snseq.combined.sct)

# need to set to integrated for clustering
DefaultAssay(burtoni.snseq.combined.sct.neurons) = 'integrated'

#check data loaded correctly
## run PCA, UMAP, and cluster 
#use 0.8 resolution
burtoni.snseq.combined.sct.neurons.recluster = burtoni.snseq.combined.sct.neurons %>% 
  RunPCA() %>%
  FindNeighbors(dims = 1:15) %>%
  RunUMAP(dims = 1:15) %>%
  FindClusters(resolution = 0.8)

#remove large file
rm(burtoni.snseq.combined.sct.neurons)

## graph 
# idents to new clusters
Idents(object = burtoni.snseq.combined.sct.neurons.recluster) <- "integrated_snn_res.0.8"

# DimPlot(burtoni.snseq.combined.sct.neurons.recluster,
#         group.by='sctypemarkers.hypo',
#         label=TRUE) +
#   umap_theme() +
#   ggtitle('Neurons') +
#   NoLegend()

# #check data
# #use all neurons
# #with variable features (top 2000 variable genes)
# png('neurons/Clusters.dimplot.neurons.all.png',
#     width = 10,
#     height = 10,
#     units = 'in',
#     res = 300)
# DimPlot(burtoni.snseq.combined.sct.neurons.recluster,
#         group.by='integrated_snn_res.0.8',
#         label=TRUE) +
#   umap_theme() +
#   ggtitle('Neurons') +
#   NoLegend()
# dev.off()
# # type
# png('neurons/Celltype.dimplot.neurons.all.png',
#     width = 10,
#     height = 10,
#     units = 'in',
#     res = 300)
# DimPlot(burtoni.snseq.combined.sct.neurons.recluster,
#         group.by='sctypemarkers.hypo',
#         label=TRUE) +
#   umap_theme() +
#   ggtitle('Neurons') +
#   NoLegend()
# dev.off()


## for paper
# just cluster labels
DimPlot(burtoni.snseq.combined.sct.neurons.recluster,
        group.by='integrated_snn_res.0.8',
        alpha = 0,
        label=TRUE) +
  ggtitle('Neurons') +
  NoLegend()+
  theme(
    plot.title = element_blank()
  )
ggsave('neurons/Clusters.labels.dimplot.neurons.all.paper.png',
       height = 7,
       width = 7)

# just cell type label
DimPlot(burtoni.snseq.combined.sct.neurons.recluster,
        group.by='sctypemarkers.hypo',
        alpha = 0,
        label=T,
        label.box = T,
        cols = c('#00BA38',
                 '#B79F00')) +
  ggtitle('Neurons') +
  NoLegend()+
  theme(
    plot.title = element_blank()
  )
ggsave('neurons/Celltype.labels.dimplot.neurons.all.paper.png',
       height = 7,
       width = 7)

# just color cells
DimPlot(burtoni.snseq.combined.sct.neurons.recluster,
        group.by='sctypemarkers.hypo',
        label=F,
        cols = c('#00BA38',
                 '#B79F00')) +
  ggtitle('Neurons') +
  NoLegend()+
  theme(
    plot.title = element_blank()
  )
ggsave('neurons/Celltype.dimplot.neurons.all.paper.png',
       height = 7,
       width = 7)
# 
# # for poster
# #check data
# #use all neurons
# #with variable features (top 2000 variable genes)
# pdf('neurons/Clusters.dimplot.neurons.all.poster.pdf',
#     width = 5,
#     height = 5)
# DimPlot(burtoni.snseq.combined.sct.neurons.recluster,
#         group.by='integrated_snn_res.0.8',
#         label=TRUE) +
#   umap_theme() +
#   ggtitle('Neurons UMAP') +
#   NoLegend()
# dev.off()
# 
# #check per genotype
# #use all neurons 
# #with variable features (top 2000 variable genes)
# png('neurons/Clusters.dimplot.neurons.all.genotype.png',
#     width = 10,
#     height = 10,
#     units = 'in',
#     res = 300)
# DimPlot(burtoni.snseq.combined.sct.neurons.recluster,
#         group.by='Genotype.id') +
#   umap_theme() +
#   ggtitle('Neurons')
# dev.off()
# 
# ## create count across genotype per cluster
# burtoni.snseq.combined.sct.neurons.recluster@meta.data %>% 
#   as.data.frame() %>% 
#   select(c(Genotype.id,
#            integrated_snn_res.0.8)) %>% 
#   table() %>% 
#   as.data.frame() %>% 
#   mutate(orig.ident = ifelse(grepl("Dom", 
#                                    Genotype.id),
#                              "dom",
#                              "sub")) %>% 
#   ggplot(aes(x = integrated_snn_res.0.8,
#              y = Freq,
#              color = orig.ident)) +
#   geom_point() +
#   theme_classic()
# ggsave("neurons/Genotype.cluster.neurons.all.count.png",
#        width = 10,
#        height = 10)

# # select 3000 variable genes
# burtoni.snseq.combined.sct.neurons.recluster <- FindVariableFeatures(burtoni.snseq.combined.sct.neurons.recluster,
#                                                                          selection.method = "vst",
#                                                                          nfeatures = 3000,
#                                                                          verbose = F)

# # Identify the 10 most highly variable genes
# top10.neurons <- head(VariableFeatures(burtoni.snseq.combined.sct.neurons.recluster), 
#                       10)

# # plot variable features with labels
# png('neurons/variable_genes/variance by expression.png',
#     width = 10,
#     height = 10,
#     units = 'in',
#     res = 300)
# LabelPoints(plot = VariableFeaturePlot(burtoni.snseq.combined.sct.neurons.recluster), 
#             points = top10.neurons, 
#             repel = TRUE)
# dev.off()
# 
# #graph variable features
# HVFInfo(object = burtoni.snseq.combined.sct.neurons.recluster,
#         status = TRUE) %>% 
#   ggplot(aes(x = residual_variance,
#              fill = variable)) + 
#   geom_histogram() +
#   theme_classic() +
#   ggtitle('Top 2000 variable genes')
# ggsave('neurons/variable_genes/histogram variable genes.png',
#        width = 10,
#        height = 10)

# HVFInfo(object = burtoni.snseq.combined.sct.neurons.recluster,
#         status = TRUE) %>% 
#   ggplot(aes(x = residual_variance,
#              fill = variable)) + 
#   geom_density() +
#   theme_classic()+
#   ggtitle('Top 2000 variable genes')
# ggsave('neurons/variable_genes/density variable genes.png',
#        width = 10,
#        height = 10)

#gene_select parameter:

#   variable: use the genes stored in the Seurat object’s VariableFeatures.
# fraction: use genes that are expressed in a certain fraction of cells for in the whole dataset or in each group of cells, specified by group.by.
# custom: use genes that are specified in a custom list.

# burtoni.snseq.combined.sct.neurons.recluster <- SetupForWGCNA(
#   burtoni.snseq.combined.sct.neurons.recluster,
#   gene_select = "variable", # the gene selection approach
#   wgcna_name = "neurons" # the name of the hdWGCNA experiment
# )

#### neuron tree and percentage figures ####
# https://romanhaa.github.io/projects/scrnaseq_workflow/#clustering

### create tree of cell types
burtoni.snseq.combined.sct.neurons.recluster.tree <- BuildClusterTree(
  burtoni.snseq.combined.sct.neurons.recluster,
  dims = 1:15,
  reorder = FALSE,
  reorder.numeric = FALSE,
  slot = 'data',
  assay = "SCT"
)

## create tree
tree <- burtoni.snseq.combined.sct.neurons.recluster.tree@tools$BuildClusterTree
# tree$tip.label <- paste0("Cluster ", tree$tip.label)

### graph tree of clusters
ggtree::ggtree(tree, aes(x, y)) +
  scale_y_reverse() +
  ggtree::geom_tree() +
  ggtree::theme_tree() +
  ggtree::geom_tiplab(offset = 1) +
  ggtree::geom_tippoint(
    shape = 16,
    size = 5) +
  coord_cartesian(clip = 'off') +
  theme(plot.margin = unit(c(0,2.5,0,0),
                           'cm'))
ggsave('neurons/Cluster tree.png',
       height = 6,
       width = 6)

# # convert tree to tibble
# tree = as_tibble(tree)
# 
# ## add cell type data
# tree = full_join(tree,
#                  burtoni.snseq.combined.sct.neurons.recluster.tree@meta.data %>%
#                    dplyr::select(seurat_clusters,
#                                  sctypemarkers.hypo) %>%
#                    distinct(),
#                  by = c("label" = "seurat_clusters"))
# 
# # convert back to tree
# tree = treeio::as.treedata(tree)
# 
# ### graph tree of celltypes
# ggtree::ggtree(tree, aes(x, y)) +
#   scale_y_reverse() +
#   ggtree::geom_tree() +
#   ggtree::theme_tree() +
#   ggtree::geom_tiplab(offset = 1) +
#   ggtree::geom_tippoint(
#     aes(color = sctypemarkers.hypo),
#                         shape = 16,
#                         size = 5) +
#   coord_cartesian(clip = 'off') +
#   theme(plot.margin = unit(c(0,2.5,0,0),
#                            'cm'))
# ggsave('neurons/Cluster cell type tree.png',
#        height = 6,
#        width = 6)

### percentages by cluster
## get counts of nuclei per cluster per genotype
# create table of samples by cluster
table_samples_by_clusters <- burtoni.snseq.combined.sct.neurons.recluster@meta.data %>% 
  group_by(Genotype.id, 
           seurat_clusters) %>%
  summarize(count = n()) %>%
  spread(seurat_clusters, 
         count, 
         fill = 0) %>%
  ungroup() %>%
  mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c('Genotype.id', 
                  'total_cell_count', 
                  everything())) %>%
  arrange(factor(Genotype.id, levels = levels(burtoni.snseq.combined.sct.neurons.recluster@meta.data$Genotype.id)))

# create table of clusters by sample
table_clusters_by_samples <- burtoni.snseq.combined.sct.neurons.recluster@meta.data %>%
  group_by(Genotype.id, 
           seurat_clusters) %>%
  summarize(count = n()) %>%
  spread(Genotype.id,
         count, 
         fill = 0) %>%
  ungroup() %>%
  mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  select(c('seurat_clusters', 
           'total_cell_count', 
           everything())) %>%
  arrange(factor(seurat_clusters, levels = levels(burtoni.snseq.combined.sct.neurons.recluster@meta.data$seurat_clusters)))

## graph percentages by cluster
# get count per sample
temp_labels <- burtoni.snseq.combined.sct.neurons.recluster@meta.data %>%
  group_by(Genotype.id) %>%
  tally()

# create graph of sample by cluster
p1 <- table_samples_by_clusters %>%
  dplyr::select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'Genotype.id') %>%
  # mutate(Genotype.id = factor(Genotype.id, levels = levels(burtoni.snseq.combined.sct.neurons.recluster@meta.data$Genotype.id))) %>%
  ggplot(aes(Genotype.id, value)) +
  geom_bar(aes(fill = variable), position = 'fill', stat = 'identity') +
  geom_text(
    data = temp_labels,
    aes(x = Genotype.id, y = Inf, label = paste0(format(n, big.mark = ',', trim = TRUE)), vjust = -1),
    color = 'black', size = 2.8
  ) +
  # scale_fill_manual(name = 'Cluster') +
  scale_y_continuous(name = 'Percentage [%]', labels = scales::percent_format(), expand = c(0.01,0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  theme(
    legend.position = 'left',
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    plot.margin = margin(t = 20, r = 0, b = 0, l = 0, unit = 'pt')
  )

# get count per cluster
temp_labels <- burtoni.snseq.combined.sct.neurons.recluster@meta.data %>%
  group_by(seurat_clusters) %>%
  tally()

# create graph of cluster by sample
p2 <- table_clusters_by_samples %>%
  select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'seurat_clusters') %>%
  mutate(seurat_clusters = factor(seurat_clusters, levels = levels(burtoni.snseq.combined.sct.neurons.recluster@meta.data$seurat_clusters))) %>%
  ggplot(aes(seurat_clusters, value)) +
  geom_bar(aes(fill = variable), position = 'fill', stat = 'identity', color = 'black') +
  geom_text(
    data = temp_labels, aes(x = seurat_clusters, y = Inf, label = paste0(format(n, big.mark = ',', trim = TRUE)), vjust = -1),
    color = 'black', size = 2.8
  ) +
  scale_fill_manual(name = 'Sample', values = c("#60bb46", 
                                                "#60bb46", 
                                                "#60bb46", 
                                                "#4e499e",
                                                "#4e499e",
                                                "#4e499e"
  )) +
  scale_y_continuous(name = 'Percentage [%]', labels = scales::percent_format(), expand = c(0.01,0)) +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  theme(
    legend.position = 'right',
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 16),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_blank(),
    plot.margin = margin(t = 20, r = 0, b = 0, l = 10, unit = 'pt')
  )

ggsave(
  'neurons/Percentage genotype by clusters.png',
  p1 + p2 +
    patchwork::plot_layout(ncol = 2, widths = c(
      burtoni.snseq.combined.sct.neurons.recluster@meta.data$Genotype.id %>% unique() %>% length(),
      burtoni.snseq.combined.sct.neurons.recluster@meta.data$seurat_clusters %>% unique() %>% length()
    )),
  width = 18, height = 8
)

## graph neuron bias
# count number of neurons for each sample type
table_clusters_by_samples %>%
  reshape2::melt(id.vars = c('seurat_clusters',
                             'total_cell_count')) %>% 
  mutate(percentage = 100*value/total_cell_count,
         orig.ident = ifelse(grepl("Dom",
                                   variable),
                             "Dom",
                             "Sub")) %>% 
  group_by(orig.ident) %>% 
  summarise(count = sum(value))
# orig.ident count
# <chr>      <dbl>
#   1 Dom         6454
# 2 Sub         6187
# 6456/12641
#51%

table_clusters_by_samples %>%
  reshape2::melt(id.vars = c('seurat_clusters',
                             'total_cell_count')) %>% 
  mutate(percentage = 100*value/total_cell_count,
         orig.ident = ifelse(grepl("Dom",
                                   variable),
                             "Dom",
                             "Sub")) %>% 
  group_by(orig.ident,
           seurat_clusters) %>% 
  summarise(total.percentage = sum(percentage)) %>% 
  ggplot(aes(x = seurat_clusters,
             y = total.percentage)) +
  geom_bar(aes(fill = orig.ident), 
           position = 'fill', 
           stat = 'identity') +
  geom_hline(yintercept = 0.51,
             linetype = 'solid')+
  geom_hline(yintercept = 0.66,
             linetype = 'dashed') +
  geom_hline(yintercept = 0.33,
             linetype = 'dashed') +
  theme_classic() +
  scale_fill_manual(values = c("#60bb46",
                               "#4e499e")) +
  ggtitle("Neuron cluster bias") +
  xlab('neuron cluster') 
ggsave('neurons/Percentage bias by clusters.png',
       height = 8,
       width = 18)
  
  
  
#### WGCNA for Neurons ####
burtoni.snseq.combined.sct.neurons.recluster = SetupForWGCNA(burtoni.snseq.combined.sct.neurons.recluster,
                                                                 features = VariableFeatures(burtoni.snseq.combined.sct.neurons.recluster),
                                                                 wgcna_name = "neurons"
)

# construct metacells  in each group
# burtoni.snseq.combined.sct.neurons.recluster <- MetacellsByGroups(
#   seurat_obj = burtoni.snseq.combined.sct.neurons.recluster,
#   group.by = c("Genotype.id"), # specify the columns in seurat_obj@meta.data to group by
#   k = 25, # nearest-neighbors parameter
#   max_shared = 5, # maximum number of shared cells between two metacells
#   ident.group = 'Genotype.id', # set the Idents of the metacell seurat object
#   slot = "data", #set for counts
#   assay = "integrated", #set assay type
#   mode = "average", #metacell expression profile determined by average of constituent single cells
#   min_cells = 50 #set minimum number of cells per metacell
# )

burtoni.snseq.combined.sct.neurons.recluster <- MetacellsByGroups(
  seurat_obj = burtoni.snseq.combined.sct.neurons.recluster,
  group.by = c("Genotype.id",
               "orig.ident"), # specify the columns in seurat_obj@meta.data to group by
  k = 25, # nearest-neighbors parameter
  max_shared = 5, # maximum number of shared cells between two metacells
  ident.group = 'orig.ident', # set the Idents of the metacell seurat object
  slot = "data", #set for counts
  assay = "SCT", #set assay type
  mode = "average", #metacell expression profile determined by average of constituent single cells
  min_cells = 50 #set minimum number of cells per metacell
)

# ## check metacell clusters
# #cluster metacells
# seurat_obj = burtoni.snseq.combined.sct.neurons.recluster
# # seurat_obj <- NormalizeMetacells(seurat_obj)
# seurat_obj <- ScaleMetacells(seurat_obj, 
#                              features=VariableFeatures(seurat_obj))
# seurat_obj <- RunPCAMetacells(seurat_obj, 
#                               features=VariableFeatures(seurat_obj))
# seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_metacell_obj <- FindNeighbors(seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_metacell_obj, 
#                                                                                     dims = 1:15)
# # seurat_obj <- RunHarmonyMetacells(seurat_obj, 
#                                   # group.by.vars='Genotype.id'
#                                   # )
# seurat_obj <- RunUMAPMetacells(seurat_obj,  
#                                dims=1:15)
# seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_metacell_obj <- Seurat::FindClusters(seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_metacell_obj, 
#                                                                                            resolution = 0.6)
# #graph 
# # metacell by genotype
# png('neurons/wgcna/metacell/Metacell UMAP genotype.png',
#     width = 10,
#     height = 10,
#     units = 'in',
#     res = 300)
# DimPlotMetacells(seurat_obj,
#                  group.by = c('Genotype.id')) +
#   umap_theme() 
# dev.off()
# # metacell by status
# png('neurons/wgcna/metacell/Metacell UMAP social status.png',
#     width = 10,
#     height = 10,
#     units = 'in',
#     res = 300)
# DimPlotMetacells(seurat_obj,
#                  group.by = c('orig.ident')) +
#   umap_theme() 
# dev.off()
# # metacell by cluster
# png('neurons/wgcna/metacell/Metacell UMAP clusters.png',
#     width = 10,
#     height = 10,
#     units = 'in',
#     res = 300)
# DimPlotMetacells(seurat_obj,
#                  group.by = c('seurat_clusters')) +
#   umap_theme() 
# dev.off()
# # metacell proportion by genotype
# seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_metacell_obj@meta.data %>% 
#   select(c(orig.ident,
#            seurat_clusters)) %>% 
#   table() %>% 
#   as.data.frame() %>% 
#   ggplot(aes(y = orig.ident,
#              x = seurat_clusters,
#              fill = Freq,
#              label = Freq)) +
#   geom_tile(color = 'black') +
#   geom_text() + 
#   scale_fill_gradient(low = "white", high = "red") +
#   coord_fixed()
# ggsave('neurons/wgcna/metacell/Metacell heatmap clusters by status.png')
# # metacell proportion by genotype
# seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_metacell_obj@meta.data %>% 
#   select(c(Genotype.id,
#            seurat_clusters)) %>% 
#   table() %>% 
#   as.data.frame() %>% 
#   ggplot(aes(y = Genotype.id,
#              x = seurat_clusters,
#              fill = Freq,
#              label = Freq)) +
#   geom_tile(color = 'black') +
#   geom_text() + 
#   scale_fill_gradient(low = "white", high = "red") +
#   coord_fixed()
# ggsave('neurons/wgcna/metacell/Metacell heatmap clusters by genotype.png')

# ### clustree
# # cluster across resolutions
# seurat_obj@misc[[seurat_9obj@misc$active_wgcna]]$wgcna_metacell_obj <- Seurat::FindClusters(seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_metacell_obj, 
#                                                                               resolution = resolution.range.reduced.2)
# #set presentation colors
# presentation.color <- c('#66c2a5',
#                         '#fc8d62',
#                         '#8da0cb',
#                         '#e78ac3',
#                         '#a6d854',
#                         '#ffd92f',
#                         '#e5c494',
#                         '#b3b3b3')
# # #clustree
# DimPlotMetacells(seurat_obj,
#                  group.by = c('integrated_snn_res.0.8')) +
#   umap_theme()

#remove object
# rm(seurat_obj)

## setup expression matrix
#use metacell expression data
burtoni.snseq.combined.sct.neurons.recluster <- SetDatExpr(
  burtoni.snseq.combined.sct.neurons.recluster,
  assay = 'SCT', # using SCT assay
  slot = 'data' # using normalized data
)

##select soft-power threshold
burtoni.snseq.combined.sct.neurons.recluster <- TestSoftPowers(
  burtoni.snseq.combined.sct.neurons.recluster,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list.neuron <- PlotSoftPowers(burtoni.snseq.combined.sct.neurons.recluster)

# assemble with patchwork
png('neurons/wgcna/Soft power threshold.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)
wrap_plots(plot_list.neuron,
           ncol=2)
dev.off()

#use soft power threshold 9?


### construct co-expression network
# construct co-expression network:
burtoni.snseq.combined.sct.neurons.recluster <- ConstructNetwork(burtoni.snseq.combined.sct.neurons.recluster, 
                                                                     soft_power=10,
                                                                     # setDatExpr=FALSE,
                                                                     tom_name = 'neuron.unknown', # name of the topoligical overlap matrix written to disk
                                                                     overwrite_tom = TRUE,
                                                                 
)

#graph dendrogram
png('neurons/wgcna/WGCNA dendrogram.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)
PlotDendrogram(burtoni.snseq.combined.sct.neurons.recluster,
               main='neurons hdWGCNA Dendrogram')
dev.off()

#graph dendrogram for poster
# pdf('neurons/wgcna/WGCNA dendrogram poster.pdf',
#     width = 6.18,
#     height = 5.65)
# PlotDendrogram(burtoni.snseq.combined.sct.neurons.recluster,
#                main='Neurons hdWGCNA Dendrogram')
# dev.off()
# ended up exporting the right size with plot viewer

##compute module eigengenes
# need to run ScaleData first or else harmony throws an error:
# burtoni.snseq.combined.sct.neurons.recluster <- ScaleData(burtoni.snseq.combined.sct.neurons.recluster, 
#                                                               features=VariableFeatures(burtoni.snseq.combined.sct.neurons.recluster))

# compute all MEs in the full single-cell dataset

burtoni.snseq.combined.sct.neurons.recluster <- ModuleEigengenes(
  burtoni.snseq.combined.sct.neurons.recluster,
  exclude_grey = TRUE
  # , group.by.vars = 'Genotype.id'
)

# module eigengenes:
MEs.neurons <- GetMEs(burtoni.snseq.combined.sct.neurons.recluster,
                      harmonized=FALSE)


# compute eigengene-based connectivity (kME):
burtoni.snseq.combined.sct.neurons.recluster <- ModuleConnectivity(
  burtoni.snseq.combined.sct.neurons.recluster)


#### Graph WGCNNA for neurons ####

# plot genes ranked by kME for each module
png('neurons/wgcna/Plot kMEs neurons.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)
PlotKMEs(burtoni.snseq.combined.sct.neurons.recluster, 
         ncol=7)
dev.off()

# correlation between modules
png('neurons/wgcna/Correlogram modules neurons.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)
ModuleCorrelogram(burtoni.snseq.combined.sct.neurons.recluster,
                  features = "MEs",
                  col = rev(COL2('RdBu', 50)),
                  addCoef.col = 'black',
                  sig.level = c(0.001), 
                  pch.cex = 0.9,
                  insig = 'blank',
                  diag = FALSE,
                  order = 'hclust')
dev.off()

# make a featureplot of hMEs for each module
plot_list.neurons.me <- ModuleFeaturePlot(
  burtoni.snseq.combined.sct.neurons.recluster,
  features='MEs', # plot the MEs
  order=TRUE # order so the points with highest MEs are on top
)

# stitch together with patchwork
png('neurons/wgcna/UMAP MEs neurons.png',
    width = 6.5,
    height = 6.5,
    units = 'in',
    res = 300)
wrap_plots(plot_list.neurons.me, ncol=4)
dev.off()

# for poster
# stitch together with patchwork
# pdf('neurons/wgcna/UMAP MEs neurons poster.pdf',
#     width = 10,
#     height = 5)
# wrap_plots(plot_list.neurons.me, ncol=4)
# dev.off()

## create one umap plot
# combine umap data and meta data
# neuron.umap = burtoni.snseq.combined.sct.neurons.recluster@reductions$umap@cell.embeddings %>% 
#   as.data.frame() %>% 
#   rownames_to_column('Cell.id') %>% 
#   full_join(burtoni.snseq.combined.sct.neurons.recluster@meta.data %>% 
#               as.data.frame() %>% 
#               rownames_to_column('Cell.id'))

# add cluster_color for each cell
# neuron.umap.long = neuron.umap %>% 
#   select(c(Cell.id,
#            UMAP_1,
#            UMAP_2,
#            orig.ident,
#            Genotype.id,
#            Cell.type,
#            integrated_snn_res.0.8,
#            turquoise,
#            green,
#            yellow,
#            red,
#            pink,
#            blue,
#            brown,
#            black,
#            grey)) %>% 
#   pivot_longer(cols = c(turquoise,
#                         green,
#                         yellow,
#                         red,
#                         pink,
#                         blue,
#                         brown,
#                         black,
#                         grey),
#                names_to = 'module',
#                values_to = 'ME_score') %>% 
#   group_by(Cell.id) %>% 
#   mutate(max.ME_score = max(ME_score),
#          second.max.ME_score = max(ME_score[ME_score != max(ME_score)]),
#          max.module = ifelse(max.ME_score == ME_score,
#                              module,
#                              NA),
#          second.max.module = ifelse(second.max.ME_score == ME_score,
#                              module,
#                              NA),
#          ratio.max.mes = second.max.ME_score/max.ME_score,
#          keep = case_when(
#            ratio.max.mes < 0 ~ max.module,
#            ratio.max.mes < 0.75 ~ max.module,
#            ratio.max.mes > 1  ~ max.module,
#            TRUE  ~ paste(max.module,
#                          second.max.module,
#                          sep = ':'),
#          )) %>% 
#   filter(!is.na(max.module)) 
# 
# # graph umap
# neuron.umap.long %>% 
#   ggplot(aes(x = UMAP_1,
#              y = UMAP_2,
#              color = max.module)) +
#   geom_point() +
#   theme_classic()


# get mods from object
mods.neurons <- colnames(MEs.neurons); mods.neurons <- mods.neurons[mods.neurons != 'grey']

# add MEs to Seurat meta-data:
burtoni.snseq.combined.sct.neurons.recluster@meta.data <- cbind(burtoni.snseq.combined.sct.neurons.recluster@meta.data,
                                                                    MEs.neurons)

##neuropeptides
modules.neurons <- GetModules(burtoni.snseq.combined.sct.neurons.recluster)
# 
# #heatmap of kme and neuropeptide
# neuropeptides.df = data.frame(gene_name = neuropeptides.list %>%
#                                 filter(!is.na(Gene.name.nile.tilapia)) %>%
#                                 pull(Gene.name.nile.tilapia) %>%
#                                 unique(),
#                               neuropeptide = TRUE)
# 
# #combine modules with neuropeptide
# modules.neurons = modules.neurons %>%
#   full_join(neuropeptides.df)

# #convert to long
# modules.long.neurons = modules.neurons %>% 
#   pivot_longer(kME_turquoise:kME_red ,
#                names_to = 'module.kME',
#                values_to = 'kME')
# 
# 
# #graph
# modules.long.neurons %>% 
#   filter(neuropeptide == TRUE) %>% 
#   drop_na() %>%
#   mutate(module.fct=as.integer(module)) %>% 
#   ggplot(aes(x = module.kME,
#              y = reorder(gene_name,
#                          module.fct),
#              label = round(kME,
#                            2))) +
#   geom_tile(aes(fill = kME)) +
#   geom_text() +
#   geom_point(aes(x=-Inf,
#                  color = module),
#              size = 5) +
#   scale_fill_gradientn(colours=c("blue",
#                                  "white",
#                                  "red"), 
#                        limits = c(-0.6, 0.6)) +
#   theme_classic() +
#   ylab('Neuropeptides') +
#   scale_color_manual(values = c("blue"="blue",
#                                 "grey"="grey",
#                                 "turquoise"="turquoise",
#                                 "brown"="brown",
#                                 "red"="red",
#                                 "green"="green",
#                                 "black"="black",
#                                 "yellow"="yellow")) +
#   coord_cartesian(clip = 'off')+ 
#   theme(axis.ticks.y = element_blank())
# ggsave('neurons/wgcna/Module and neuropeptide heatmap.png',
#        width = 10,
#        height = 10)
# 
# ## create dotplot of modules ME and clusters
# DotPlot(burtoni.snseq.combined.sct.neurons.recluster,
#         features = c("blue", 
#                      "red", 
#                      "yellow", 
#                      "green",
#                      "turquoise",
#                      "brown"),
#         cols = c('grey',
#                  'red'),
#         group.by = "integrated_snn_res.0.8",
#         col.min = 2,
#         dot.min = .5,
#         scale = FALSE) +
#   RotatedAxis() +
#   theme_bw()
# ggsave('neurons/wgcna/ME by cluster dotplot.png',
#        width = 10,
#        height = 10)
# 
# ## get dot plot data
# # relabel column 
burtoni.neuron.wgcna.module.cluster = DotPlot.data(burtoni.snseq.combined.sct.neurons.recluster,
                                                   features = c("blue",
                                                                "red",
                                                                "yellow",
                                                                "green",
                                                                "turquoise",
                                                                "brown"),
                                                   cols = c('grey',
                                                            'red'),
                                                   group.by = "integrated_snn_res.0.8",
                                                   col.min = 2,
                                                   dot.min = .5,
                                                   scale = FALSE)
# 
# 
# ## trait correlation
# Idents(burtoni.snseq.combined.sct.neurons.recluster) <- burtoni.snseq.combined.sct.neurons.recluster$orig.ident
# 
# # convert orig.ident to factor
# burtoni.snseq.combined.sct.neurons.recluster$orig.ident.fct <- as.factor(burtoni.snseq.combined.sct.neurons.recluster$orig.ident)
# 
# # list of traits to correlate
# cur_traits.neurons <- c('orig.ident.fct', 
#                         'nCount_SCT', 
#                         'nFeature_SCT')
# 
# burtoni.snseq.combined.sct.neurons.recluster <- ModuleTraitCorrelation(
#   burtoni.snseq.combined.sct.neurons.recluster,
#   traits = cur_traits.neurons,
#   features = 'MEs'
# )
# 
# #plot
# png('neurons/wgcna/Module and trait correlation heatmap.png',
#     width = 10,
#     height = 10,
#     units = 'in',
#     res = 300)
# PlotModuleTraitCorrelation(
#   burtoni.snseq.combined.sct.neurons.recluster,
#   label = 'fdr',
#   label_symbol = 'stars',
#   text_size = 5,
#   text_digits = 5,
#   text_color = 'black',
#   high_color = 'red',
#   mid_color = 'white',
#   low_color = 'blue',
#   plot_max = 0.5,
#   combine=TRUE
# )
# dev.off()

#### Network WGCNA for neurons ####
# visualize network with UMAP
burtoni.snseq.combined.sct.neurons.recluster <- RunModuleUMAP(
  burtoni.snseq.combined.sct.neurons.recluster,
  n_hubs = 5, #number of hub genes to include for the umap embedding
  n_neighbors=15, #neighbors parameter for umap
  min_dist=0.3, #min distance between points in umap space
  spread=5
)


# compute cell-type marker genes with Seurat:
#set idents to cluster
Idents(burtoni.snseq.combined.sct.neurons.recluster) <- burtoni.snseq.combined.sct.neurons.recluster$integrated_snn_res.0.8

#calculate marker genes per cluster
# with 3000 variable genes
markers.neuron <- Seurat::FindAllMarkers(
  burtoni.snseq.combined.sct.neurons.recluster,
  only.pos = TRUE,
  logfc.threshold=1,
  features = VariableFeatures(burtoni.snseq.combined.sct.neurons.recluster)
)

# compute marker gene overlaps
burtoni.snseq.combined.sct.neurons.recluster.overlap_df <- OverlapModulesDEGs(
  burtoni.snseq.combined.sct.neurons.recluster,
  deg_df = markers.neuron,
  fc_cutoff = 1 # log fold change cutoff for overlap analysis
)

# overlap barplot, produces a plot for each cluster
plot_list.neuron.overlap <- OverlapBarPlot(burtoni.snseq.combined.sct.neurons.recluster.overlap_df)

# stitch plots with patchwork
png('neurons/wgcna/Vlnplot modules vs cluster odds ratio.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)
wrap_plots(plot_list.neuron.overlap, 
           ncol=4)
dev.off()

# plot odds ratio of the overlap as a dot plot
png('neurons/wgcna/Dotplot modules vs cluster odds ratio.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)
OverlapDotPlot(burtoni.snseq.combined.sct.neurons.recluster.overlap_df,
               plot_var = 'odds_ratio') +
  ggtitle('Overlap of modules & cluster markers')
dev.off()

##graph network
#umap
png('neurons/wgcna/gene_network/Gene network UMAP.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)
ModuleUMAPPlot(burtoni.snseq.combined.sct.neurons.recluster,
               edge.alpha=0.5,
               sample_edges=TRUE,
               keep_grey_edges=FALSE,
               edge_prop=0.075, # taking the top 20% strongest edges in each module
               label_hubs=10 # how many hub genes to plot per module?
)
dev.off()

# for poster
##graph network
#umap
# pdf('neurons/wgcna/gene_network/Gene network UMAP poster.pdf',
#     width = 10,
#     height = 10)
# ModuleUMAPPlot.size(burtoni.snseq.combined.sct.neurons.recluster,
#                edge.alpha= 0.75,
#                sample_edges=TRUE,
#                keep_grey_edges=FALSE,
#                edge_prop=0.075, # taking the top 20% strongest edges in each module
#                label_hubs=0, # how many hub genes to plot per module?
#                vertex.label.cex = 1,
#                dot.size = 8,
#                edge.size = 4
# )
# dev.off()

# module specific network
ModuleNetworkPlot(burtoni.snseq.combined.sct.neurons.recluster,
                  outdir = 'neurons/wgcna/gene_network/')

# hubgene network
png('neurons/wgcna/gene_network/Hubgene network.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)
HubGeneNetworkPlot(burtoni.snseq.combined.sct.neurons.recluster,
                   n_hubs = 5, 
                   n_other=50,
                   edge_prop = 0.75,
                   mods = 'all'
)
dev.off()

# 
# # get list of marker genes
# write_csv(markers.neuron,
#           './neurons/markers.neuron.csv')



#### DMEs WGCNA for neurons ####
## Differential module
# get list of dom and sub cells
group1.neurons.dom <- burtoni.snseq.combined.sct.neurons.recluster@meta.data %>% 
  subset(orig.ident == 'dom_burtoni_snseq') %>%
  rownames
group2.neurons.sub <- burtoni.snseq.combined.sct.neurons.recluster@meta.data %>% 
  subset(orig.ident == 'sub_burtoni_snseq') %>% 
  rownames

# calculate differential module eigengene
DMEs.neurons <- FindDMEs(
  burtoni.snseq.combined.sct.neurons.recluster,
  barcodes1 = group1.neurons.dom,
  barcodes2 = group2.neurons.sub,
  test.use='wilcox',
  harmonized = FALSE
)

# graph DME
PlotDMEsVolcano(
  burtoni.snseq.combined.sct.neurons.recluster,
  DMEs.neurons) +
  theme_classic() +
  xlim(-0.5, 0.5)
ggsave('neurons/wgcna/ME by cluster volcanoplot.png',
       width = 10,
       height = 10)

# for poster
# create mod color list
mod_colors <- DMEs.neurons$module
names(mod_colors) <- as.character(DMEs.neurons$module)
# graph DME
DMEs.neurons %>%
  mutate(anno = ifelse(p_val_adj < 0.05,
                       module,
                       "")) %>%
  ggplot(aes(x = avg_log2FC,
             y = -log10(p_val_adj)))  +
  geom_rect(aes(xmin = -Inf,
                xmax = Inf,
                ymin = -Inf,
                ymax = -log10(0.05)),
            fill = "grey75",
            alpha = 0.8,
            color = NA) +
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "grey75",
             alpha = 0.8) +
  geom_point(aes(fill = module),
             size = 8,
             pch = 21,
             color = "black") +
  geom_text_repel(aes(label = anno),
                  color = "black",
                  min.segment.length = 0,
                  max.overlaps = Inf,
                  size = 8,
                  point.padding = 4,
                  nudge_x = 0.06) +
  theme_classic() +
  theme(legend.position = 'none') +
  xlim(c(-0.55,
         0.55)) +
  scale_fill_manual(values = mod_colors)+
  xlab(bquote("Average log"[2] ~ "(Fold Change)")) +
  ylab(bquote("-log"[10] ~ "(Adj. P-value)")) +
  ggtitle('Differential module eigengene analysis') +
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1))+
  theme(axis.text = element_text(size = 15))  +
  theme(axis.title = element_text(size = 20))+
  theme(plot.title = element_text(size=20))
ggsave('neurons/wgcna/ME by cluster volcanoplot poster.pdf',
       width = 10,
       height = 5,
       units = "in",
       dpi = 320)

# for presentation
DMEs.neurons %>%
  mutate(anno = ifelse(p_val_adj < 0.05,
                       module,
                       "")) %>%
  ggplot(aes(x = avg_log2FC,
             y = -log10(p_val_adj)))  +
  geom_rect(aes(xmin = -Inf,
                xmax = Inf,
                ymin = -Inf,
                ymax = -log10(0.01)),
            fill = "grey75",
            alpha = 0.8,
            color = NA) +
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "grey75",
             alpha = 0.8) +
  geom_point(aes(fill = module),
             size = 8,
             pch = 21,
             color = "black")  +
  theme_classic() +
  theme(legend.position = 'none') +
  xlim(c(-0.55,
         0.55)) +
  scale_fill_manual(values = mod_colors)+
  xlab(bquote("Average log"[2] ~ "(Fold Change)")) +
  ylab(bquote("-log"[10] ~ "(Adj. P-value)")) +
  ggtitle('Differential module eigengene analysis') +
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1))+
  theme(axis.text = element_text(size = 15))  +
  theme(axis.title = element_text(size = 20))+
  theme(plot.title = element_text(size=20))
ggsave('neurons/wgcna/ME by cluster volcanoplot presentation.pdf',
       width = 10,
       height = 5,
       units = "in",
       dpi = 320)

# ### graph MEs for relevant clusters
# for (i in mods.neurons) {
#   #create list of clusters with module eigengene above:
#   # average scaled expression of 1 and 50% percent expression
#   cluster.list = burtoni.neuron.wgcna.module.cluster %>% 
#     filter(features.plot == i) %>% 
#     mutate(keep = ifelse(avg.exp.scaled >= 1 & pct.exp >= 50,
#                          "keep",
#                          "remove")) %>% 
#     filter(keep == "keep") %>%
#     mutate(id = as.character(id)) %>% 
#     pull(id) 
#   
#   ##graph ME per cluster across social status
#   burtoni.snseq.combined.sct.neurons.recluster@meta.data %>% 
#     mutate(cluster = as.character(integrated_snn_res.0.8)) %>% 
#     filter(cluster %in% cluster.list) %>% 
#     droplevels() %>%
#     ggplot(aes_string(x = 'integrated_snn_res.0.8',
#                       y = i,
#                       color = 'orig.ident')) +
#     geom_boxplot(outlier.shape = NA) +
#     geom_point(position=position_jitterdodge(jitter.width = 0.25,
#                                              jitter.height = 0),
#                alpha = 0.2)+
#     theme_classic() +
#     ggtitle(i)
#   ggsave(paste('neurons/wgcna/modules/',
#                i,
#                ' ME by cluster boxplot.png',
#                sep = ''),
#          width = 10,
#          height = 10)
#   
#   ## ME across social status
#   burtoni.snseq.combined.sct.neurons.recluster@meta.data %>% 
#     mutate(cluster = as.character(integrated_snn_res.0.8)) %>% 
#     filter(cluster %in% cluster.list) %>% 
#     droplevels() %>%
#     ggplot(aes_string(x = 'orig.ident',
#                       y = i)) +
#     geom_violin(draw_quantiles = c(0.5)) +
#     geom_jitter(width = 0.25,
#                 height = 0)+
#     theme_classic() +
#     ggtitle(i)
#   ggsave(paste('neurons/wgcna/modules/',
#                i,
#                ' ME violinplot.png',
#                sep = ''),
#          width = 10,
#          height = 10)
#   
#   # cluster bias 
#   burtoni.snseq.combined.sct.neurons.recluster@meta.data %>% 
#     mutate(cluster = as.character(integrated_snn_res.0.8)) %>% 
#     filter(cluster %in% cluster.list) %>% 
#     mutate(keep = ifelse(get(i) >= 1,
#                          "present",
#                          "absence")) %>% 
#     droplevels() %>%
#     select(c(cluster,
#              keep,
#              orig.ident)) %>% 
#     table() %>% 
#     as.data.frame() %>% 
#     pivot_wider(id_cols = c(cluster,
#                             keep),
#                 names_from = orig.ident,
#                 values_from = Freq) %>% 
#     filter(keep == "present") %>% 
#     mutate(dom_burtoni_snseq.scaled = 0.783645 * dom_burtoni_snseq) %>% 
#     ggplot(aes(x = dom_burtoni_snseq.scaled,
#                y = sub_burtoni_snseq,
#                label = cluster)) +
#     geom_abline(slope = 1,
#                 intercept = 0) +
#     geom_text() +
#     theme_classic() +
#     ggtitle(paste(i,
#                   'sub vs dom scaled counts per cluster'))
#   ggsave(paste('neurons/wgcna/modules/',
#                i,
#                ' ME count bias.png',
#                sep = ''),
#          width = 10,
#          height = 10)
#   
# }

#### save point ####
# save hub genes for GO analysis with:
hub.genes.neurons = GetModules(burtoni.snseq.combined.sct.neurons.recluster)

#save hub genes
write_csv(hub.genes.neurons,
          './neurons/wgcna/hub.genes.neurons.csv')
# load hubgenes
hub.genes.neurons = read.csv('./neurons/wgcna/hub.genes.neurons.csv')

# get list of marker genes
write_csv(markers.neuron,
          './neurons/markers.neuron.csv')
# load marker genes neurons
markers.neuron = read.csv('./neurons/markers.neuron.csv')

## just single cell object
save(burtoni.snseq.combined.sct.neurons.recluster,
     file = "./neurons/wgcna/burtoni.snseq.combined.sct.neurons.recluster.RData")

load('./neurons/wgcna/burtoni.snseq.combined.sct.neurons.recluster.RData')
#### Create GO database ####
# making GO library
library(AnnotationForge)

## make library from NCBI
makeOrgPackageFromNCBI(version="0.1",
                       maintainer="Isaac Miller-Crews <imillerc@iu.edu>",
                       author="Isaac Miller-Crews <imillerc@iu.edu>",
                       outputDir = "./",
                       tax_id = "8128",
                       genus = "Oreochromis",
                       species = "niloticus")


# then you can call install.packages based on the return value
install.packages("./org.Oniloticus.eg.db", 
                 repos=NULL)


#### Create GO figures per module ####
# load libraries
library(enrichplot)
library(clusterProfiler)
library(org.Oniloticus.eg.db) # need to make database in code above

# load all genes for GO comparison
all.genes = burtoni.snseq.combined.sct.neurons.recluster@misc[["neurons"]][["wgcna_genes"]]

## get list of module colors
module.colors = module.genes %>% 
  filter(module != 'grey') %>% 
  pull(color) %>% 
  unique()

# convert ensembl gene name to gene name using enrichgo
tmp.enrichgo.ensembl <- enrichGO(gene = all.genes, 
                                 universe = all.genes,
                                 OrgDb = "org.Oniloticus.eg.db",
                                 keyType = "ENSEMBL",
                                 ont = "BP",
                                 pvalueCutoff = 0.15,
                                 pAdjustMethod = 'fdr',
                                 minGSSize = 10,
                                 maxGSSize = 500,
                                 readable = T)
# get gene names
tmp.enrichgo.ensembl.genes = tmp.enrichgo.ensembl@gene2Symbol %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  na.omit() %>% 
  dplyr::rename(ENSEMBL = V1) %>% 
  rownames_to_column('gene_name') %>% 
  separate_wider_delim(gene_name,
                       delim = '.',
                       names = c('gene_name',
                                 NA),
                       too_few = 'align_start') #some ID have duplicate genes
  

# get dataframe of genes and modules
module.genes = burtoni.snseq.combined.sct.neurons.recluster@misc[["neurons"]][["wgcna_modules"]] %>% 
  dplyr::select(gene_name,
                module)

# combine with ensembl genes
module.genes.ensembl = module.genes %>% 
  left_join(tmp.enrichgo.ensembl.genes) %>% 
  mutate(gene_name = ifelse(is.na(ENSEMBL),
                            gene_name,
                            ENSEMBL)) %>% 
  dplyr::select(-c(ENSEMBL)) %>% 
  distinct()

# add ensembl genes to all genes
all.genes.ensembl = c(all.genes,
                      tmp.enrichgo.ensembl.genes$ENSEMBL) %>% 
  unique()

### loop through modules for enriched genes
# create empty data frame to save results
enrichGO.results.modules = data.frame()

for (i in module.colors) {
  
  # run enrichgo  
  # gene symbol
  tmp.enrichgo <- enrichGO(gene = module.genes.ensembl %>% 
                             filter(module == i) %>% 
                             pull(gene_name) %>% 
                             unique(), 
                           universe = all.genes.ensembl,
                           OrgDb = "org.Oniloticus.eg.db",
                           keyType = "SYMBOL",
                           ont = "BP",
                           pvalueCutoff = 0.1,
                           pAdjustMethod = 'fdr',
                           minGSSize = 10,
                           maxGSSize = 500,
                           readable = T)
  
  # save results
  enrichGO.results.modules = enrichGO.results.modules %>% 
    rbind(tmp.enrichgo@result %>% 
                mutate(module = i,
                       ont = 'BP'))
  
  write.csv(GO.results.modules,
            'neurons/wgcna/go_terms/enrichGO.results.modules.csv',
            row.names = F)
  
  # check number of sig GO
  tmp.num = tmp.enrichgo@result %>% 
    filter(p.adjust < 0.1) %>% 
    nrow()
  
  if (tmp.num != 0) {
  # graph GO results
  # create plot
  tmp.enrichgo.fit <- plot(barplot(tmp.enrichgo,
                                   showCategory = 15))+
    ggtitle(paste0(i,
                   ': module enriched BP GO'))
  
  # save plot
  png(paste0('neurons/wgcna/go_terms/enrichplot/',
             i,
             ' module enriched BP GO.png'),
      res = 300,
      width = 10,
      height = 10,
      units = 'in')
  print(tmp.enrichgo.fit)
  dev.off()
  
  # create GO term tree
  # https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html#tree-plot
  # download enrichplot
  # remotes::install_github("GuangchuangYu/enrichplot")
  
  # graph tree
  png(paste0('neurons/wgcna/go_terms/enrichplot/',
             i,
             ' module enriched BP GO treeplot.png'), 
      res = 300, 
      width = 10, 
      height = 10,
      units = 'in')
  tmp.enrichgo %>% 
    pairwise_termsim() %>% 
    treeplot() +
    ggtitle(paste0(i,
                   ': module enriched BP GO tree'))
  dev.off()
  }
  if (tmp.num == 0) {
    # save plot
    png(paste0('neurons/wgcna/go_terms/enrichplot/',
               i,
               ' module enriched BP GO.png'),
        res = 300,
        width = 10,
        height = 10,
        units = 'in')
    plot.new()
    dev.off()
    
    # create GO term tree
    # https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html#tree-plot
    # download enrichplot
    # remotes::install_github("GuangchuangYu/enrichplot")
    
    # graph tree
    png(paste0('neurons/wgcna/go_terms/enrichplot/',
               i,
               ' module enriched BP GO treeplot.png'), 
        res = 300, 
        width = 10, 
        height = 10,
        units = 'in')
    plot.new()
    dev.off()
  }
}

# #### Get GO terms for MWU GO ####
# ## load ensembl
# library(biomaRt)
# ###selecting biomart database
# ensembl = useMart('ensembl')
# # #list datasets
# # dataset = listDatasets(ensembl)
# ##select dataset
# #nile tilapia
# ensembl.tilapia = useDataset('oniloticus_gene_ensembl',
#                            mart=ensembl)
# 
# # # get list of all attributes
# # listAttributes(ensembl.tilapia) %>% View
# 
# 
# #create attributes lists
# tilapia.attributes = c('external_gene_name',
#                                 'ensembl_gene_id',
#                                 'go_id')
# 
# tilapia.attributes.2 = c('ensembl_gene_id',
#                        'go_id')
# 
# ##identify GO terms for WGCNA tilapia genes
# # use gene names
# tilapia.wgcna.go.terms.gene = getBM(attributes = tilapia.attributes,
#                                   mart = ensembl.tilapia,
#                                   values = hub.genes.neurons$gene_name,
#                                   filter = 'external_gene_name',
#                                   useCache = FALSE) # useCache has to do with version of R not being up to date?
# # use ensembl gene IDs
# tilapia.wgcna.go.terms.ensemblID = getBM(attributes = tilapia.attributes.2,
#                                     mart = ensembl.tilapia,
#                                     values = hub.genes.neurons$gene_name,
#                                     filter = 'ensembl_gene_id',
#                                     useCache = FALSE) # useCache has to do with version of R not being up to date?
# # combine gene and ensembl gene IDs
# tilapia.wgcna.go.terms = full_join(tilapia.wgcna.go.terms.gene,
#                                tilapia.wgcna.go.terms.ensemblID)
# 
# # single gene column 
# tilapia.wgcna.go.terms = tilapia.wgcna.go.terms %>% 
#   mutate(gene_name = ifelse(is.na(external_gene_name),
#                             ensembl_gene_id,
#                             external_gene_name)) %>% 
#   dplyr::select(c(gene_name,
#                   go_id))
# 
# #check length
# #12162
# #14276
# tilapia.wgcna.go.terms %>%
#   nrow()
# 
# #check number of tilapia genes?
# #2429
# #2939
# tilapia.wgcna.go.terms %>%
#   pull(gene_name) %>%
#   unique() %>%
#   length()
# 
# #check number of GO IDs?
# #1906
# #2019
# tilapia.wgcna.go.terms %>%
#   pull(go_id) %>%
#   unique() %>%
#   length()
# 
# ## check duplicates
# # 11644
# # 13990
# tilapia.wgcna.go.terms %>%
#   distinct() %>% 
#   nrow()
# # 12162 - 11644 = 518 duplicates
# # 14276 - 13990 = 286 duplicates
# 
# # remove duplicates
# tilapia.wgcna.go.terms = tilapia.wgcna.go.terms %>% 
#   distinct()
# 
# ## collapse go_id into gene list
# tilapia.wgcna.go.terms = tilapia.wgcna.go.terms %>%
#   group_by(gene_name) %>%
#   summarize(go_id = str_c(go_id, collapse = ";"))
# 
# # add 'unknown' 
# tilapia.wgcna.go.terms = tilapia.wgcna.go.terms %>% 
#   mutate(go_id = ifelse(go_id == '',
#                         'unknown',
#                         go_id))
# 
# # save to csv 
# write.csv(tilapia.wgcna.go.terms,
#           file = './neurons/wgcna/tilapia.wgcna.go.terms.csv')
# # save to tab delimited with no column names
# write_tsv(tilapia.wgcna.go.terms,
#           file = './neurons/wgcna/go_terms/tilapia.wgcna.go.terms.tsv',
#           col_names = FALSE)

# #### prepare for GO_MWU ####
# ### create dataframe of WGCNA genes with average module log fold change
# wgcna.DMEs.neurons.log2FC = full_join(hub.genes.neurons %>% 
#                                         dplyr::select(c(gene_name,
#                                                        module)),
#                                       DMEs.neurons %>% 
#   dplyr::select(c(avg_log2FC, 
#                   module))) %>% 
#   filter(module != 'grey') %>% 
#   dplyr::select(-c(module))
# # save
# write.csv(wgcna.DMEs.neurons.log2FC,
#           file = './neurons/wgcna/go_terms/wgcna.DMEs.neurons.log2FC.csv',
#           row.names = FALSE,
#           quote = FALSE)
# 
# ### create data input list for each module (remove grey) 
# # want a list of genes, the kme if they are in that module, and a 0 kME if they are not in that module
# module.name.list = hub.genes.neurons %>% 
#   pull(module) %>% 
#   unique() 
# # remove grey module
# module.name.list = module.name.list[! module.name.list %in% c("grey")]
# 
# ## loop through to create data frame for each module
# for (i in module.name.list) {
#   tmp.name = paste("kME_",
#                    i,
#                    sep='')
#   # create temporary dataframe with module kme
#   tmp = hub.genes.neurons %>% 
#     dplyr::select(c(gene_name,
#                     module,
#                     tmp.name)) %>% 
#     mutate(kme = ifelse(module != i,
#                              0,
#                              .[[tmp.name]])) %>% 
#     dplyr::select(c(gene_name,
#                     kme))
#   # save file
#   write.csv(tmp,
#             file = paste('./neurons/wgcna/go_terms/wgcna.DMEs.neurons.kme.',
#             i,
#             '.csv',
#             sep = ''),
#             row.names = FALSE,
#             quote = FALSE)
# }
# 
# 
# 
# ## need wgcna.DMEs.neurons.log2FC.csv and tilapia.wgcna.go.terms.tsv for MWU_GO
# # need to run on desktop 
# # https://github.com/z0on/GO_MWU
# # https://github.com/schmidte10/GO-MWU-automation-


#### GLM binomial neuron ####
## scale counts
# celltype
# genotype
# dataframe of scale
neuron.genotype.scale = burtoni.snseq.combined.sct.neurons.recluster@meta.data %>%
  dplyr::select(c(orig.ident,
                  Genotype.id)) %>% 
  table() %>%
  as.data.frame() %>% 
  dplyr::rename(total.neuron.count = Freq) %>% 
  filter(total.neuron.count != 0)

#dataframe of cluster counts
neuron.cluster.genotype.count = burtoni.snseq.combined.sct.neurons.recluster@meta.data %>%
  dplyr::select(c(orig.ident,
                  seurat_clusters,
                  Genotype.id)) %>% 
  table() %>%
  as.data.frame() %>% 
  filter(Freq != 0) 

#combine counts and scale
neuron.cluster.genotype.count = neuron.cluster.genotype.count %>% 
  full_join(neuron.genotype.scale) %>% 
  mutate(Percent = 100*Freq/total.neuron.count)
### add variable for status
neuron.cluster.genotype.count = neuron.cluster.genotype.count %>% 
  mutate(Status = ifelse(orig.ident == "dom_burtoni_snseq",
                         "dom",
                         "sub")) 

### run glm on every cluster
## use with proportion data
## create empty data frame
# p-value
neuron.cluster.glm = data.frame(matrix(ncol = 9, nrow = 0))

#provide column names
colnames(neuron.cluster.glm) <- c(
  "seurat_clusters",
  "contrast",
  "Estimate",
  "Std.error",
  "z.value",
  "adj.pvalue",
  "z.ratio",
  'lwr',
  'upr')


## loop through each cluster
for (i in unique(neuron.cluster.genotype.count$seurat_clusters)) {
  ## run interaction binomial GLM on status
  tmp = glm(cbind(Freq, Others) ~ Status, 
            data = neuron.cluster.genotype.count %>% 
              filter(seurat_clusters == i) %>% 
              mutate(Others = total.neuron.count - Freq,
                     Status = as.factor(Status)), 
            family = binomial)
  # multivariate t distribution to correct for multiple testing
  tmp.res = summary(glht(tmp))
  
  # save results
  tmp.res.df = data.frame(Estimate = tmp.res$test$coefficients,
                          Std.error = tmp.res$test$sigma,
                          z.value = tmp.res$test$tstat,
                          adj.pvalue = tmp.res$test$pvalues
                          # z.ratio = NA
  ) %>% 
    rownames_to_column('contrast') %>% 
    mutate(seurat_clusters = i)
  
  # graph confidence intervals
  # get confidence intervals
  ci_out <- confint(tmp.res)
  ci_df <- as.data.frame(ci_out$confint)
  ci_df$contrast <- row.names(ci_df)
  ci_df$seurat_clusters = i
  
  # ## run all pairwise comparison
  # # Tukey needed for all pairwise comparisons
  # tmp2 = glm(cbind(Freq, Others) ~ orig.ident, 
  #            data = neuron.cluster.genotype.count %>% 
  #              filter(seurat_clusters == i) %>% 
  #              mutate(Others = total.neuron.count - Freq,
  #                     Status = as.factor(Status)), 
  #            family = binomial)
  # 
  # # use emmeans to get pairwise differences
  # emm = emmeans(tmp2,
  #               ~ "orig.ident")
  # emm.res = pairs(emm)
  # 
  # 
  # # save results
  # emm.res.df = data.frame(Estimate = summary(emm.res)$estimate,
  #                         Std.error = summary(emm.res)$SE,
  #                         z.ratio = summary(emm.res)$z.ratio,
  #                         adj.pvalue = summary(emm.res)$p.value,
  #                         contrast = summary(emm.res)$contrast,
  #                         z.value = NA) %>%
  #   mutate(seurat_clusters = i)
  
  ## combine in data frame
  neuron.cluster.glm = neuron.cluster.glm %>% 
    # rbind(emm.res.df) %>% 
    rbind(tmp.res.df %>% 
            full_join(ci_df)) 
  
}


## FDR correction for all comparisons
neuron.cluster.glm.fdr = neuron.cluster.glm %>% 
  filter(contrast != "(Intercept)") %>% 
  mutate(p.adjust.fdr = p.adjust(adj.pvalue,
                                 method = 'fdr'))

# combine with data frame
neuron.cluster.glm = neuron.cluster.glm %>% 
  left_join(neuron.cluster.glm.fdr)



## create rounded p.value to make it easier to read
neuron.cluster.glm = neuron.cluster.glm %>% 
  mutate(adj.pvalue.round = ifelse(is.na(p.adjust.fdr),
                                   round(adj.pvalue, digits = 4),
                                   round(p.adjust.fdr, digits = 4)))

# ## create p value dataframe
# neuron.cluster.glm.p = neuron.cluster.glm %>% 
#   dplyr::select(seurat_clusters,
#                 adj.pvalue.round,
#                 contrast) %>% 
#   pivot_wider(names_from = contrast,
#               values_from = adj.pvalue.round) 

## graph results
# graph clusters per status
neuron.cluster.genotype.count %>% 
  ggplot() +
  geom_point(aes(x = seurat_clusters,
                 y = Percent,
                 group = Status,
                 color = Status),
             position = position_dodge(width = 0.5),
             size = 3) +
  geom_text(data = neuron.cluster.glm %>% 
              filter(contrast != "(Intercept)") %>% 
              dplyr::select(contrast,
                            seurat_clusters,
                            adj.pvalue.round) %>% 
              mutate(sig = case_when(adj.pvalue.round < 0.01  ~ '*',
                                     TRUE ~ NA)) %>% 
              na.omit(),
            aes(x = seurat_clusters,
                y = max(neuron.cluster.genotype.count$Percent),
                label = sig),
            size = 5) +
  theme_classic() +
  geom_vline(xintercept = seq(0.5,
                              length(unique(neuron.cluster.genotype.count$seurat_clusters)),
                              by = 1),
             color = 'gray75',
             size = 0.5,
             linetype = 2
  ) +
  scale_color_manual(values = c("#4e499e",
                                "#60bb46")) +
  ggtitle("Neuron clusters: binomial GLM cell count") +
  xlab('neuron clusters')
ggsave('neurons/neuron_cluster/Neuron clusters binomial GLM cell percent.png',
       height = 4,
       width = 8,
       units = 'in')

## forest plot
# get max difference value for x axis limits on graph
max.diff.value = neuron.cluster.glm %>% 
  filter(contrast != "(Intercept)") %>% 
  pull(Estimate) %>% 
  abs() %>% 
  max()
# add pvalue
neuron.cluster.glm %>% 
  filter(contrast != "(Intercept)") %>% 
  mutate(effect = case_when(adj.pvalue.round < 0.01 & Estimate > 0 ~ 'sub_bias',
                            adj.pvalue.round < 0.01 & Estimate < 0 ~ 'dom_bias',
                            TRUE ~ 'none')) %>% 
  ggplot(aes(x = Estimate, 
             y = as.factor(as.numeric(seurat_clusters)))) +
  geom_errorbar(aes(xmin = lwr, 
                    xmax = upr), 
                linewidth = 1) +
  geom_point(aes(color = effect),
             size = 3) +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = 0.5,
             linetype = 2) +
  geom_vline(xintercept = -0.5,
             linetype = 2) +
  geom_text(aes(label = seurat_clusters),
            color = 'white',
            size = 2) +
  labs(x = "Effect size", 
       y = "neuron clusters") +
  theme_classic() +
  scale_color_manual(breaks = c('dom_bias',
                                'sub_bias',
                                'none'),
                     values = c(
                       "#4e499e",
                       "#60bb46",
                       "grey25")) +
  # ggtitle("Neuron clusters: binomial GLM cell count") +
  theme(legend.position = 'none') +
  xlim(-max.diff.value-0.25,
       max.diff.value+0.25) 
ggsave('neurons/neuron_cluster/Neuron clusters binomial GLM forest plot.png',
       height = 3.76,
       width = 2.76,
       units = 'in')

# ## check which neuron clusters are GABA or GLU
# burtoni.snseq.combined.sct.neurons.recluster@meta.data %>%
#   dplyr::select(c(seurat_clusters,
#                   sctypemarkers.hypo)) %>% 
#   table() %>% 
#   data.frame() %>% 
#   pivot_wider(values_from = 'Freq',
#               names_from = 'sctypemarkers.hypo') %>% 
#   mutate(ratio = `C7-1: GLU`/`C7-2: GABA`,
#          log.ratio = log(ratio))

### compare to average ME 
neuron.cluster.glm = burtoni.snseq.combined.sct.neurons.recluster@meta.data %>% 
  group_by(seurat_clusters) %>% 
  summarise(brown = mean(brown),
            green = mean(green),
            red = mean(red),
            blue = mean(blue),
            black = mean(black),
            turquoise = mean(turquoise),
            yellow = mean(yellow)) %>% 
  pivot_longer(cols = -c(seurat_clusters),
               values_to = 'average.ME',
               names_to = 'module') %>% 
  group_by(seurat_clusters) %>% 
  mutate(max.average.ME = max(average.ME),
         max = ifelse(max.average.ME == average.ME,
                      1,
                      0),
         ratio = 100*average.ME/max.average.ME) %>% 
  filter(max == 1) %>% 
  dplyr::select(seurat_clusters,
                module) %>% 
  distinct() %>% 
  right_join(neuron.cluster.glm)

## save glm table
write_csv(neuron.cluster.glm,
          file = 'neurons/neuron_cluster/neuron_seurat_cluster_glm.csv')


## check module eigengene vs significant clusters
neuron.cluster.glm %>% 
  filter(contrast != "(Intercept)") %>% 
  filter(p.adjust.fdr < 0.01) %>% 
  view()


#### GLM binomial all ####
## scale counts
# celltype
# genotype
# dataframe of scale
genotype.scale = burtoni.snseq.combined.sct@meta.data %>%
  dplyr::select(c(orig.ident,
                  Genotype.id)) %>% 
  table() %>%
  as.data.frame() %>% 
  dplyr::rename(total.count = Freq) %>% 
  filter(total.count != 0)

#dataframe of cluster counts
cluster.genotype.count = burtoni.snseq.combined.sct@meta.data %>%
  dplyr::select(c(orig.ident,
                  seurat_clusters,
                  Genotype.id)) %>% 
  table() %>%
  as.data.frame() %>% 
  filter(Freq != 0) 

#combine counts and scale
cluster.genotype.count = cluster.genotype.count %>% 
  full_join(genotype.scale) %>% 
  mutate(Percent = 100*Freq/total.count)
### add variable for status
cluster.genotype.count = cluster.genotype.count %>% 
  mutate(Status = ifelse(orig.ident == "dom_burtoni_snseq",
                         "dom",
                         "sub")) 

### run glm on every cluster
## use with proportion data
## create empty data frame
# p-value
cluster.glm = data.frame(matrix(ncol = 9, nrow = 0))

#provide column names
colnames(cluster.glm) <- c(
  "seurat_clusters",
  "contrast",
  "Estimate",
  "Std.error",
  "z.value",
  "adj.pvalue",
  "z.ratio",
  'lwr',
  'upr')


## loop through each cluster
for (i in unique(cluster.genotype.count$seurat_clusters)) {
  ## run interaction binomial GLM on status
  tmp = glm(cbind(Freq, Others) ~ Status, 
            data = cluster.genotype.count %>% 
              filter(seurat_clusters == i) %>% 
              mutate(Others = total.count - Freq,
                     Status = as.factor(Status)), 
            family = binomial)
  # multivariate t distribution to correct for multiple testing
  tmp.res = summary(glht(tmp))
  
  # save results
  tmp.res.df = data.frame(Estimate = tmp.res$test$coefficients,
                          Std.error = tmp.res$test$sigma,
                          z.value = tmp.res$test$tstat,
                          adj.pvalue = tmp.res$test$pvalues
                          # z.ratio = NA
  ) %>% 
    rownames_to_column('contrast') %>% 
    mutate(seurat_clusters = i)
  
  # graph confidence intervals
  # get confidence intervals
  ci_out <- confint(tmp.res)
  ci_df <- as.data.frame(ci_out$confint)
  ci_df$contrast <- row.names(ci_df)
  ci_df$seurat_clusters = i
  
  # ## run all pairwise comparison
  # # Tukey needed for all pairwise comparisons
  # tmp2 = glm(cbind(Freq, Others) ~ orig.ident, 
  #            data = neuron.cluster.genotype.count %>% 
  #              filter(seurat_clusters == i) %>% 
  #              mutate(Others = total.neuron.count - Freq,
  #                     Status = as.factor(Status)), 
  #            family = binomial)
  # 
  # # use emmeans to get pairwise differences
  # emm = emmeans(tmp2,
  #               ~ "orig.ident")
  # emm.res = pairs(emm)
  # 
  # 
  # # save results
  # emm.res.df = data.frame(Estimate = summary(emm.res)$estimate,
  #                         Std.error = summary(emm.res)$SE,
  #                         z.ratio = summary(emm.res)$z.ratio,
  #                         adj.pvalue = summary(emm.res)$p.value,
  #                         contrast = summary(emm.res)$contrast,
  #                         z.value = NA) %>%
  #   mutate(seurat_clusters = i)
  
  ## combine in data frame
  cluster.glm = cluster.glm %>% 
    # rbind(emm.res.df) %>% 
    rbind(tmp.res.df %>% 
            full_join(ci_df)) 
  
}


## FDR correction for all comparisons
cluster.glm.fdr = cluster.glm %>% 
  filter(contrast != "(Intercept)") %>% 
  mutate(p.adjust.fdr = p.adjust(adj.pvalue,
                                 method = 'fdr'))

# combine with data frame
cluster.glm = cluster.glm %>% 
  left_join(cluster.glm.fdr)



## create rounded p.value to make it easier to read
cluster.glm = cluster.glm %>% 
  mutate(adj.pvalue.round = ifelse(is.na(p.adjust.fdr),
                                   round(adj.pvalue, digits = 4),
                                   round(p.adjust.fdr, digits = 4)))

## save glm table
write_csv(cluster.glm,
          file = 'neurons/neuron_cluster/seurat_cluster_glm.csv')

# ## create p value dataframe
# neuron.cluster.glm.p = neuron.cluster.glm %>% 
#   dplyr::select(seurat_clusters,
#                 adj.pvalue.round,
#                 contrast) %>% 
#   pivot_wider(names_from = contrast,
#               values_from = adj.pvalue.round) 

## graph results
# graph clusters per status
cluster.genotype.count %>% 
  ggplot() +
  geom_point(aes(x = seurat_clusters,
                 y = Percent,
                 group = Status,
                 color = Status),
             position = position_dodge(width = 0.5),
             size = 3) +
  # geom_text(data = cluster.glm %>% 
  #             filter(contrast != "(Intercept)") %>% 
  #             dplyr::select(contrast,
  #                           seurat_clusters,
  #                           adj.pvalue.round) %>% 
  #             mutate(sig = case_when(adj.pvalue.round < 0.05 & adj.pvalue.round > 0.001 ~ '*',
  #                                    adj.pvalue.round < 0.001 & adj.pvalue.round > 0.0001 ~ '**',
  #                                    adj.pvalue.round < 0.0001 ~ '***',
  #                                    TRUE ~ NA)) %>% 
  #             na.omit(),
  #           aes(x = seurat_clusters,
  #               y = max(cluster.genotype.count$Percent),
  #               label = sig),
  #           size = 5) +
  geom_text(data = cluster.glm %>% 
              filter(contrast != "(Intercept)") %>% 
              dplyr::select(contrast,
                            seurat_clusters,
                            adj.pvalue.round,
                            Estimate) %>% 
              mutate(sig = case_when(adj.pvalue.round < 0.01  ~ '*',
                                     TRUE ~ NA)) %>% 
              na.omit(),
            aes(x = seurat_clusters,
                y = max(cluster.genotype.count$Percent),
                label = sig),
            size = 5) +
  theme_classic() +
  geom_vline(xintercept = seq(0.5,
                              length(unique(cluster.genotype.count$seurat_clusters)),
                              by = 1),
             color = 'gray75',
             size = 0.5,
             linetype = 2
  ) +
  scale_color_manual(values = c("#4e499e",
                                "#60bb46")) +
  ggtitle("All clusters: binomial GLM cell count") +
  xlab('clusters') 
ggsave('neurons/neuron_cluster/clusters binomial GLM cell percent.png',
       height = 4,
       width = 8,
       units = 'in')

## forest plot
# get max difference value for x axis limits on graph
max.diff.value = cluster.glm %>% 
  filter(contrast != "(Intercept)") %>% 
  pull(Estimate) %>% 
  abs() %>% 
  max()
# add pvalue
cluster.glm %>% 
  filter(contrast != "(Intercept)") %>% 
  mutate(effect = case_when(adj.pvalue.round < 0.01 & Estimate > 0 ~ 'sub_bias',
                            adj.pvalue.round < 0.01 & Estimate < 0 ~ 'dom_bias',
                            TRUE ~ 'none')) %>% 
  ggplot(aes(x = Estimate, 
             y = as.factor(as.numeric(seurat_clusters)))) +
  geom_errorbar(aes(xmin = lwr, 
                    xmax = upr), 
                linewidth = 1) +
  geom_point(aes(color = effect),
             size = 3) +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = 0.5,
             linetype = 2) +
  geom_vline(xintercept = -0.5,
             linetype = 2) +
  geom_text(aes(label = seurat_clusters),
            color = 'white',
            size = 2) +
  labs(x = "Effect size", 
       y = "all clusters") +
  theme_classic() +
  scale_color_manual(breaks = c('dom_bias',
                                'sub_bias',
                                'none'),
                     values = c(
                       "#4e499e",
                       "#60bb46",
                       "grey25")) +
  xlim(-max.diff.value-0.25,
       max.diff.value+0.25) +
  # ggtitle("All clusters: binomial GLM cell count") +
  theme(legend.position = 'none')
ggsave('neurons/neuron_cluster/clusters binomial GLM forest plot.png',
       height = 6,
       width = 3,
       units = 'in')

# ## check which neuron clusters are GABA or GLU
# burtoni.snseq.combined.sct.neurons.recluster@meta.data %>%
#   dplyr::select(c(seurat_clusters,
#                   sctypemarkers.hypo)) %>% 
#   table() %>% 
#   data.frame() %>% 
#   pivot_wider(values_from = 'Freq',
#               names_from = 'sctypemarkers.hypo') %>% 
#   mutate(ratio = `C7-1: GLU`/`C7-2: GABA`,
#          log.ratio = log(ratio))



#### create figures for poster old ####
### dotplot of broad cell types across UMAP

#assign broad cell class
burtoni.scsorter.data.scsort.output.broad = burtoni.scsorter.data.scsort.output %>% 
  select(Cell.type,
         Cell.id) %>%
  mutate(Cell.class.broad = case_when(Cell.type == 'neurons' ~ 'neurons',
                                      Cell.type == 'inhibitory' ~ 'neurons',
                                      Cell.type == 'excitatory' ~ 'neurons',
                                      Cell.type == 'mural' ~ 'vascular',
                                      Cell.type == 'VSM' ~ 'vascular',
                                      Cell.type == 'endothelial' ~ 'vascular',
                                      Cell.type == 'Unknown' ~ 'unknown',
                                      TRUE ~ 'glia'))

##combine umap with metadata
burtoni.scsorter.data.scsort.output.umap =  full_join(burtoni.scsorter.data.scsort.output.broad,
                                                      burtoni.scsorter.data.scsort.output.umap)

##graph umap with cell markers
burtoni.scsorter.data.scsort.output.umap %>% 
  ggplot(aes(x=UMAP_1,
             y= UMAP_2,
             color = Cell.type)) +
  geom_point() +
  theme_bw()
ggsave('scsorter/UMAP scsorter cell types.png',
       width = 10,
       height = 10)


#graph umap with cell markers
#broad cell class
burtoni.scsorter.data.scsort.output.umap %>% 
  ggplot(aes(x=UMAP_1,
             y= UMAP_2,
             color = Cell.class.broad)) +
  geom_point() +
  theme_classic() +
  scale_color_manual(values = c('3C51A3',
                                'orange3',
                                'grey',
                                '4AB85C'))
ggsave('scsorter/UMAP scsorter cell cell types broad.pdf',
       width = 10,
       height = 10)

#graph umap with cell markers
#broad cell class
# seperate doms and subs
burtoni.scsorter.data.scsort.output.umap %>% 
  ggplot(aes(x=UMAP_1,
             y= UMAP_2,
             color = Cell.class.broad)) +
  geom_point() +
  theme_classic() +
  scale_color_manual(values = c('3C51A3',
                                'orange3',
                                'grey',
                                '4AB85C')) +
  facet_grid( ~ orig.ident) +
  labs(color='Broad cell classes') +
  theme(legend.text=element_text(size=20),
        legend.title=element_text(size=20))+
  guides(colour = guide_legend(override.aes = list(size=10)))
ggsave('scsorter/UMAP scsorter cell cell types broad sub vs dom.pdf',
       width = 12,
       height = 10)





#### WGCNA for neurons old ####
burtoni.snseq.combined.sct.all.neurons = burtoni.snseq.combined.sct.all

#set idents
Idents(object = burtoni.snseq.combined.sct.all.neurons) <- "Cell.type"

#subset to neurons
burtoni.snseq.combined.sct.all.neurons = subset(burtoni.snseq.combined.sct.all.neurons,
                                                idents = c("neurons",
                                                           "excitatory",
                                                           "inhibitory"))

#remove large file
rm(burtoni.snseq.combined.sct.all)

# need to set to integrated for clustering
DefaultAssay(burtoni.snseq.combined.sct.all.neurons) = 'integrated'

#check data loaded correctly
## run PCA, UMAP, and cluster 
#use 0.8 resolution
burtoni.snseq.combined.sct.all.neurons.recluster = burtoni.snseq.combined.sct.all.neurons %>% 
  RunPCA() %>%
  FindNeighbors(dims = 1:15) %>%
  RunUMAP(dims = 1:15) %>%
  FindClusters(resolution = 0.8)

#remove large file
rm(burtoni.snseq.combined.sct.all.neurons)

## graph 
# idents to new clusters
Idents(object = burtoni.snseq.combined.sct.all.neurons.recluster) <- "integrated_snn_res.0.8"

# #check data
# #use all neurons 
# #with variable features (top 2000 variable genes)
# png('neurons/Clusters.dimplot.neurons.all.png',
#     width = 10,
#     height = 10,
#     units = 'in',
#     res = 300)
# DimPlot(burtoni.snseq.combined.sct.all.neurons.recluster,
#         group.by='integrated_snn_res.0.8',
#         label=TRUE) +
#   umap_theme() +
#   ggtitle('Neurons') +
#   NoLegend()
# dev.off()
# 
# # for poster
# #check data
# #use all neurons
# #with variable features (top 2000 variable genes)
# pdf('neurons/Clusters.dimplot.neurons.all.poster.pdf',
#     width = 5,
#     height = 5)
# DimPlot(burtoni.snseq.combined.sct.all.neurons.recluster,
#         group.by='integrated_snn_res.0.8',
#         label=TRUE) +
#   umap_theme() +
#   ggtitle('Neurons UMAP') +
#   NoLegend()
# dev.off()
# 
# #check per genotype
# #use all neurons 
# #with variable features (top 2000 variable genes)
# png('neurons/Clusters.dimplot.neurons.all.genotype.png',
#     width = 10,
#     height = 10,
#     units = 'in',
#     res = 300)
# DimPlot(burtoni.snseq.combined.sct.all.neurons.recluster,
#         group.by='Genotype.id') +
#   umap_theme() +
#   ggtitle('Neurons')
# dev.off()
# 
# ## create count across genotype per cluster
# burtoni.snseq.combined.sct.all.neurons.recluster@meta.data %>% 
#   as.data.frame() %>% 
#   select(c(Genotype.id,
#            integrated_snn_res.0.8)) %>% 
#   table() %>% 
#   as.data.frame() %>% 
#   mutate(orig.ident = ifelse(grepl("Dom", 
#                                    Genotype.id),
#                              "dom",
#                              "sub")) %>% 
#   ggplot(aes(x = integrated_snn_res.0.8,
#              y = Freq,
#              color = orig.ident)) +
#   geom_point() +
#   theme_classic()
# ggsave("neurons/Genotype.cluster.neurons.all.count.png",
#        width = 10,
#        height = 10)

# # select 3000 variable genes
# burtoni.snseq.combined.sct.all.neurons.recluster <- FindVariableFeatures(burtoni.snseq.combined.sct.all.neurons.recluster,
#                                                                          selection.method = "vst",
#                                                                          nfeatures = 3000,
#                                                                          verbose = F)

# # Identify the 10 most highly variable genes
# top10.neurons <- head(VariableFeatures(burtoni.snseq.combined.sct.all.neurons.recluster), 
#                       10)

# # plot variable features with labels
# png('neurons/variable_genes/variance by expression.png',
#     width = 10,
#     height = 10,
#     units = 'in',
#     res = 300)
# LabelPoints(plot = VariableFeaturePlot(burtoni.snseq.combined.sct.all.neurons.recluster), 
#             points = top10.neurons, 
#             repel = TRUE)
# dev.off()
# 
# #graph variable features
# HVFInfo(object = burtoni.snseq.combined.sct.all.neurons.recluster,
#         status = TRUE) %>% 
#   ggplot(aes(x = residual_variance,
#              fill = variable)) + 
#   geom_histogram() +
#   theme_classic() +
#   ggtitle('Top 2000 variable genes')
# ggsave('neurons/variable_genes/histogram variable genes.png',
#        width = 10,
#        height = 10)

# HVFInfo(object = burtoni.snseq.combined.sct.all.neurons.recluster,
#         status = TRUE) %>% 
#   ggplot(aes(x = residual_variance,
#              fill = variable)) + 
#   geom_density() +
#   theme_classic()+
#   ggtitle('Top 2000 variable genes')
# ggsave('neurons/variable_genes/density variable genes.png',
#        width = 10,
#        height = 10)

#gene_select parameter:

#   variable: use the genes stored in the Seurat object’s VariableFeatures.
# fraction: use genes that are expressed in a certain fraction of cells for in the whole dataset or in each group of cells, specified by group.by.
# custom: use genes that are specified in a custom list.

# burtoni.snseq.combined.sct.all.neurons.recluster <- SetupForWGCNA(
#   burtoni.snseq.combined.sct.all.neurons.recluster,
#   gene_select = "variable", # the gene selection approach
#   wgcna_name = "neurons" # the name of the hdWGCNA experiment
# )


burtoni.snseq.combined.sct.all.neurons.recluster = SetupForWGCNA(burtoni.snseq.combined.sct.all.neurons.recluster,
                                                                 features = VariableFeatures(burtoni.snseq.combined.sct.all.neurons.recluster),
                                                                 wgcna_name = "neurons"
)

# construct metacells  in each group
# burtoni.snseq.combined.sct.all.neurons.recluster <- MetacellsByGroups(
#   seurat_obj = burtoni.snseq.combined.sct.all.neurons.recluster,
#   group.by = c("Genotype.id"), # specify the columns in seurat_obj@meta.data to group by
#   k = 25, # nearest-neighbors parameter
#   max_shared = 5, # maximum number of shared cells between two metacells
#   ident.group = 'Genotype.id', # set the Idents of the metacell seurat object
#   slot = "data", #set for counts
#   assay = "integrated", #set assay type
#   mode = "average", #metacell expression profile determined by average of constituent single cells
#   min_cells = 50 #set minimum number of cells per metacell
# )

burtoni.snseq.combined.sct.all.neurons.recluster <- MetacellsByGroups(
  seurat_obj = burtoni.snseq.combined.sct.all.neurons.recluster,
  group.by = c("Genotype.id",
               "orig.ident"), # specify the columns in seurat_obj@meta.data to group by
  k = 25, # nearest-neighbors parameter
  max_shared = 5, # maximum number of shared cells between two metacells
  ident.group = 'orig.ident', # set the Idents of the metacell seurat object
  slot = "data", #set for counts
  assay = "SCT", #set assay type
  mode = "average", #metacell expression profile determined by average of constituent single cells
  min_cells = 50 #set minimum number of cells per metacell
)

# ## check metacell clusters
# #cluster metacells
# seurat_obj = burtoni.snseq.combined.sct.all.neurons.recluster
# # seurat_obj <- NormalizeMetacells(seurat_obj)
# seurat_obj <- ScaleMetacells(seurat_obj, 
#                              features=VariableFeatures(seurat_obj))
# seurat_obj <- RunPCAMetacells(seurat_obj, 
#                               features=VariableFeatures(seurat_obj))
# seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_metacell_obj <- FindNeighbors(seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_metacell_obj, 
#                                                                                     dims = 1:15)
# # seurat_obj <- RunHarmonyMetacells(seurat_obj, 
#                                   # group.by.vars='Genotype.id'
#                                   # )
# seurat_obj <- RunUMAPMetacells(seurat_obj,  
#                                dims=1:15)
# seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_metacell_obj <- Seurat::FindClusters(seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_metacell_obj, 
#                                                                                            resolution = 0.6)
# #graph 
# # metacell by genotype
# png('neurons/wgcna/metacell/Metacell UMAP genotype.png',
#     width = 10,
#     height = 10,
#     units = 'in',
#     res = 300)
# DimPlotMetacells(seurat_obj,
#                  group.by = c('Genotype.id')) +
#   umap_theme() 
# dev.off()
# # metacell by status
# png('neurons/wgcna/metacell/Metacell UMAP social status.png',
#     width = 10,
#     height = 10,
#     units = 'in',
#     res = 300)
# DimPlotMetacells(seurat_obj,
#                  group.by = c('orig.ident')) +
#   umap_theme() 
# dev.off()
# # metacell by cluster
# png('neurons/wgcna/metacell/Metacell UMAP clusters.png',
#     width = 10,
#     height = 10,
#     units = 'in',
#     res = 300)
# DimPlotMetacells(seurat_obj,
#                  group.by = c('seurat_clusters')) +
#   umap_theme() 
# dev.off()
# # metacell proportion by genotype
# seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_metacell_obj@meta.data %>% 
#   select(c(orig.ident,
#            seurat_clusters)) %>% 
#   table() %>% 
#   as.data.frame() %>% 
#   ggplot(aes(y = orig.ident,
#              x = seurat_clusters,
#              fill = Freq,
#              label = Freq)) +
#   geom_tile(color = 'black') +
#   geom_text() + 
#   scale_fill_gradient(low = "white", high = "red") +
#   coord_fixed()
# ggsave('neurons/wgcna/metacell/Metacell heatmap clusters by status.png')
# # metacell proportion by genotype
# seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_metacell_obj@meta.data %>% 
#   select(c(Genotype.id,
#            seurat_clusters)) %>% 
#   table() %>% 
#   as.data.frame() %>% 
#   ggplot(aes(y = Genotype.id,
#              x = seurat_clusters,
#              fill = Freq,
#              label = Freq)) +
#   geom_tile(color = 'black') +
#   geom_text() + 
#   scale_fill_gradient(low = "white", high = "red") +
#   coord_fixed()
# ggsave('neurons/wgcna/metacell/Metacell heatmap clusters by genotype.png')

# ### clustree
# # cluster across resolutions
# seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_metacell_obj <- Seurat::FindClusters(seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_metacell_obj, 
#                                                                               resolution = resolution.range.reduced.2)
# #set presentation colors
# presentation.color <- c('#66c2a5',
#                         '#fc8d62',
#                         '#8da0cb',
#                         '#e78ac3',
#                         '#a6d854',
#                         '#ffd92f',
#                         '#e5c494',
#                         '#b3b3b3')
# # #clustree
# DimPlotMetacells(seurat_obj,
#                  group.by = c('integrated_snn_res.0.8')) +
#   umap_theme()

#remove object
# rm(seurat_obj)

## setup expression matrix
#use metacell expression data

##select soft-power threshold
burtoni.snseq.combined.sct.all.neurons.recluster <- TestSoftPowers(
  burtoni.snseq.combined.sct.all.neurons.recluster,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list.neuron <- PlotSoftPowers(burtoni.snseq.combined.sct.all.neurons.recluster)

# assemble with patchwork
png('neurons/wgcna/Soft power threshold.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)
wrap_plots(plot_list.neuron, 
           ncol=2)
dev.off()

#use soft power threshold 10?


### construct co-expression network
# construct co-expression network:
burtoni.snseq.combined.sct.all.neurons.recluster <- ConstructNetwork(burtoni.snseq.combined.sct.all.neurons.recluster, 
                                                                     soft_power=8,
                                                                     # setDatExpr=FALSE,
                                                                     tom_name = 'neuron', # name of the topoligical overlap matrix written to disk
                                                                     overwrite_tom = TRUE
)

#graph dendrogram
png('neurons/wgcna/WGCNA dendrogram.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)
PlotDendrogram(burtoni.snseq.combined.sct.all.neurons.recluster, 
               main='neurons hdWGCNA Dendrogram')
dev.off()

#graph dendrogram for poster
# pdf('neurons/wgcna/WGCNA dendrogram poster.pdf',
#     width = 10,
#     height = 2.5)
# PlotDendrogram(burtoni.snseq.combined.sct.all.neurons.recluster,
#                main='Neurons hdWGCNA Dendrogram')
# dev.off()
# ended up exporting the right size with plot viewer

##compute module eigengenes
# need to run ScaleData first or else harmony throws an error:
# burtoni.snseq.combined.sct.all.neurons.recluster <- ScaleData(burtoni.snseq.combined.sct.all.neurons.recluster, 
#                                                               features=VariableFeatures(burtoni.snseq.combined.sct.all.neurons.recluster))

# compute all MEs in the full single-cell dataset

burtoni.snseq.combined.sct.all.neurons.recluster <- ModuleEigengenes(
  burtoni.snseq.combined.sct.all.neurons.recluster,
  exclude_grey = TRUE
  # , group.by.vars = 'Genotype.id'
)

# module eigengenes:
MEs.neurons <- GetMEs(burtoni.snseq.combined.sct.all.neurons.recluster,
                      harmonized=FALSE)


# compute eigengene-based connectivity (kME):
burtoni.snseq.combined.sct.all.neurons.recluster <- ModuleConnectivity(
  burtoni.snseq.combined.sct.all.neurons.recluster)


#### Graph WGCNNA for neurons old ####

# plot genes ranked by kME for each module
png('neurons/wgcna/Plot kMEs neurons.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)
PlotKMEs(burtoni.snseq.combined.sct.all.neurons.recluster, 
         ncol=7)
dev.off()

# correlation between modules
png('neurons/wgcna/Correlogram modules neurons.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)
ModuleCorrelogram(burtoni.snseq.combined.sct.all.neurons.recluster,
                  features = "MEs",
                  col = rev(COL2('RdBu', 50)),
                  addCoef.col = 'black',
                  sig.level = c(0.001), 
                  pch.cex = 0.9,
                  insig = 'blank',
                  diag = FALSE,
                  order = 'hclust')
dev.off()

# make a featureplot of hMEs for each module
plot_list.neurons.me <- ModuleFeaturePlot(
  burtoni.snseq.combined.sct.all.neurons.recluster,
  features='MEs', # plot the MEs
  order=TRUE # order so the points with highest MEs are on top
)

# stitch together with patchwork
png('neurons/wgcna/UMAP MEs neurons.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)
wrap_plots(plot_list.neurons.me, ncol=4)
dev.off()

# for poster
# stitch together with patchwork
# pdf('neurons/wgcna/UMAP MEs neurons poster.pdf',
#     width = 10,
#     height = 5)
# wrap_plots(plot_list.neurons.me, ncol=4)
# dev.off()

## create one umap plot
# combine umap data and meta data
# neuron.umap = burtoni.snseq.combined.sct.all.neurons.recluster@reductions$umap@cell.embeddings %>% 
#   as.data.frame() %>% 
#   rownames_to_column('Cell.id') %>% 
#   full_join(burtoni.snseq.combined.sct.all.neurons.recluster@meta.data %>% 
#               as.data.frame() %>% 
#               rownames_to_column('Cell.id'))

# add cluster_color for each cell
# neuron.umap.long = neuron.umap %>% 
#   select(c(Cell.id,
#            UMAP_1,
#            UMAP_2,
#            orig.ident,
#            Genotype.id,
#            Cell.type,
#            integrated_snn_res.0.8,
#            turquoise,
#            green,
#            yellow,
#            red,
#            pink,
#            blue,
#            brown,
#            black,
#            grey)) %>% 
#   pivot_longer(cols = c(turquoise,
#                         green,
#                         yellow,
#                         red,
#                         pink,
#                         blue,
#                         brown,
#                         black,
#                         grey),
#                names_to = 'module',
#                values_to = 'ME_score') %>% 
#   group_by(Cell.id) %>% 
#   mutate(max.ME_score = max(ME_score),
#          second.max.ME_score = max(ME_score[ME_score != max(ME_score)]),
#          max.module = ifelse(max.ME_score == ME_score,
#                              module,
#                              NA),
#          second.max.module = ifelse(second.max.ME_score == ME_score,
#                              module,
#                              NA),
#          ratio.max.mes = second.max.ME_score/max.ME_score,
#          keep = case_when(
#            ratio.max.mes < 0 ~ max.module,
#            ratio.max.mes < 0.75 ~ max.module,
#            ratio.max.mes > 1  ~ max.module,
#            TRUE  ~ paste(max.module,
#                          second.max.module,
#                          sep = ':'),
#          )) %>% 
#   filter(!is.na(max.module)) 
# 
# # graph umap
# neuron.umap.long %>% 
#   ggplot(aes(x = UMAP_1,
#              y = UMAP_2,
#              color = max.module)) +
#   geom_point() +
#   theme_classic()


# get mods from object
mods.neurons <- colnames(MEs.neurons); mods.neurons <- mods.neurons[mods.neurons != 'grey']

# add MEs to Seurat meta-data:
burtoni.snseq.combined.sct.all.neurons.recluster@meta.data <- cbind(burtoni.snseq.combined.sct.all.neurons.recluster@meta.data, 
                                                                    MEs.neurons)

##neuropeptides
modules.neurons <- GetModules(burtoni.snseq.combined.sct.all.neurons.recluster)

#heatmap of kme and neuropeptide
neuropeptides.df = data.frame(gene_name = neuropeptides.list %>% 
                                filter(!is.na(Gene.name.nile.tilapia)) %>% 
                                pull(Gene.name.nile.tilapia) %>% 
                                unique(),
                              neuropeptide = TRUE)

#combine modules with neuropeptide
modules.neurons = modules.neurons %>% 
  full_join(neuropeptides.df)

# #convert to long
# modules.long.neurons = modules.neurons %>% 
#   pivot_longer(kME_turquoise:kME_red ,
#                names_to = 'module.kME',
#                values_to = 'kME')
# 
# 
# #graph
# modules.long.neurons %>% 
#   filter(neuropeptide == TRUE) %>% 
#   drop_na() %>%
#   mutate(module.fct=as.integer(module)) %>% 
#   ggplot(aes(x = module.kME,
#              y = reorder(gene_name,
#                          module.fct),
#              label = round(kME,
#                            2))) +
#   geom_tile(aes(fill = kME)) +
#   geom_text() +
#   geom_point(aes(x=-Inf,
#                  color = module),
#              size = 5) +
#   scale_fill_gradientn(colours=c("blue",
#                                  "white",
#                                  "red"), 
#                        limits = c(-0.6, 0.6)) +
#   theme_classic() +
#   ylab('Neuropeptides') +
#   scale_color_manual(values = c("blue"="blue",
#                                 "grey"="grey",
#                                 "turquoise"="turquoise",
#                                 "brown"="brown",
#                                 "red"="red",
#                                 "green"="green",
#                                 "black"="black",
#                                 "yellow"="yellow")) +
#   coord_cartesian(clip = 'off')+ 
#   theme(axis.ticks.y = element_blank())
# ggsave('neurons/wgcna/Module and neuropeptide heatmap.png',
#        width = 10,
#        height = 10)

## create dotplot of modules ME and clusters
DotPlot(burtoni.snseq.combined.sct.all.neurons.recluster,
        features = c("blue", 
                     "red", 
                     "yellow", 
                     "green", 
                     "pink", 
                     "turquoise",
                     "brown",
                     "black"),
        cols = c('grey',
                 'red'),
        group.by = "integrated_snn_res.0.8",
        col.min = 2,
        dot.min = .5,
        scale = FALSE) +
  RotatedAxis() +
  theme_bw()
ggsave('neurons/wgcna/ME by cluster dotplot.png',
       width = 10,
       height = 10)

## get dot plot data
# relabel column 
burtoni.neuron.wgcna.module.cluster = DotPlot.data(burtoni.snseq.combined.sct.all.neurons.recluster,
                                                   features = c("blue", 
                                                                "red", 
                                                                "yellow", 
                                                                "green", 
                                                                "pink", 
                                                                "turquoise",
                                                                "brown",
                                                                "black"),
                                                   cols = c('grey',
                                                            'red'),
                                                   group.by = "integrated_snn_res.0.8",
                                                   col.min = 2,
                                                   dot.min = .5,
                                                   scale = FALSE)


## trait correlation
Idents(burtoni.snseq.combined.sct.all.neurons.recluster) <- burtoni.snseq.combined.sct.all.neurons.recluster$orig.ident

# convert orig.ident to factor
burtoni.snseq.combined.sct.all.neurons.recluster$orig.ident.fct <- as.factor(burtoni.snseq.combined.sct.all.neurons.recluster$orig.ident)

# list of traits to correlate
cur_traits.neurons <- c('orig.ident.fct', 
                        'nCount_SCT', 
                        'nFeature_SCT')

burtoni.snseq.combined.sct.all.neurons.recluster <- ModuleTraitCorrelation(
  burtoni.snseq.combined.sct.all.neurons.recluster,
  traits = cur_traits.neurons,
  features = 'MEs'
)

#plot
png('neurons/wgcna/Module and trait correlation heatmap.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)
PlotModuleTraitCorrelation(
  burtoni.snseq.combined.sct.all.neurons.recluster,
  label = 'fdr',
  label_symbol = 'stars',
  text_size = 5,
  text_digits = 5,
  text_color = 'black',
  high_color = 'red',
  mid_color = 'white',
  low_color = 'blue',
  plot_max = 0.5,
  combine=TRUE
)
dev.off()

#### Network WGCNA for neurons old ####
# visualize network with UMAP
burtoni.snseq.combined.sct.all.neurons.recluster <- RunModuleUMAP(
  burtoni.snseq.combined.sct.all.neurons.recluster,
  n_hubs = 5, #number of hub genes to include for the umap embedding
  n_neighbors=15, #neighbors parameter for umap
  min_dist=0.3, #min distance between points in umap space
  spread=5
)


# compute cell-type marker genes with Seurat:
#set idents to cluster
Idents(burtoni.snseq.combined.sct.all.neurons.recluster) <- burtoni.snseq.combined.sct.all.neurons.recluster$integrated_snn_res.0.8

#calculate marker genes per cluster
# with 3000 variable genes
markers.neuron <- Seurat::FindAllMarkers(
  burtoni.snseq.combined.sct.all.neurons.recluster,
  only.pos = TRUE,
  logfc.threshold=1,
  features = VariableFeatures(burtoni.snseq.combined.sct.all.neurons.recluster)
)

# compute marker gene overlaps
burtoni.snseq.combined.sct.all.neurons.recluster.overlap_df <- OverlapModulesDEGs(
  burtoni.snseq.combined.sct.all.neurons.recluster,
  deg_df = markers.neuron,
  fc_cutoff = 1 # log fold change cutoff for overlap analysis
)

# overlap barplot, produces a plot for each cluster
plot_list.neuron.overlap <- OverlapBarPlot(burtoni.snseq.combined.sct.all.neurons.recluster.overlap_df)

# stitch plots with patchwork
png('neurons/wgcna/Vlnplot modules vs cluster odds ratio.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)
wrap_plots(plot_list.neuron.overlap, 
           ncol=4)
dev.off()

# plot odds ratio of the overlap as a dot plot
png('neurons/wgcna/Dotplot modules vs cluster odds ratio.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)
OverlapDotPlot(burtoni.snseq.combined.sct.all.neurons.recluster.overlap_df,
               plot_var = 'odds_ratio') +
  ggtitle('Overlap of modules & cluster markers')
dev.off()

##graph network
#umap
png('neurons/wgcna/gene_network/Gene network UMAP.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)
ModuleUMAPPlot(burtoni.snseq.combined.sct.all.neurons.recluster,
               edge.alpha=0.5,
               sample_edges=TRUE,
               keep_grey_edges=FALSE,
               edge_prop=0.075, # taking the top 20% strongest edges in each module
               label_hubs=10 # how many hub genes to plot per module?
)
dev.off()

# for poster
##graph network
#umap
# pdf('neurons/wgcna/gene_network/Gene network UMAP poster.pdf',
#     width = 10,
#     height = 10)
# ModuleUMAPPlot.size(burtoni.snseq.combined.sct.all.neurons.recluster,
#                edge.alpha= 0.75,
#                sample_edges=TRUE,
#                keep_grey_edges=FALSE,
#                edge_prop=0.075, # taking the top 20% strongest edges in each module
#                label_hubs=0, # how many hub genes to plot per module?
#                vertex.label.cex = 1,
#                dot.size = 8,
#                edge.size = 4
# )
# dev.off()

# module specific network
ModuleNetworkPlot(burtoni.snseq.combined.sct.all.neurons.recluster,
                  outdir = 'neurons/wgcna/gene_network/')

# hubgene network
png('neurons/wgcna/gene_network/Hubgene network.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)
HubGeneNetworkPlot(burtoni.snseq.combined.sct.all.neurons.recluster,
                   n_hubs = 5, 
                   n_other=50,
                   edge_prop = 0.75,
                   mods = 'all'
)
dev.off()






#### DMEs WGCNA for neurons old ####
## Differential module
# get list of dom and sub cells
group1.neurons.dom <- burtoni.snseq.combined.sct.all.neurons.recluster@meta.data %>% 
  subset(orig.ident == 'dom_burtoni_snseq') %>%
  rownames
group2.neurons.sub <- burtoni.snseq.combined.sct.all.neurons.recluster@meta.data %>% 
  subset(orig.ident == 'sub_burtoni_snseq') %>% 
  rownames

# calculate differential module eigengene
DMEs.neurons <- FindDMEs(
  burtoni.snseq.combined.sct.all.neurons.recluster,
  barcodes1 = group1.neurons.dom,
  barcodes2 = group2.neurons.sub,
  test.use='wilcox',
  harmonized = FALSE
)

# graph DME
PlotDMEsVolcano(
  burtoni.snseq.combined.sct.all.neurons.recluster,
  DMEs.neurons) +
  theme_classic() +
  xlim(-0.5, 0.5)
ggsave('neurons/wgcna/ME by cluster volcanoplot.png',
       width = 10,
       height = 10)

# for poster
# create mod color list
mod_colors <- DMEs.neurons$module
names(mod_colors) <- as.character(DMEs.neurons$module)
# graph DME
DMEs.neurons %>%
  mutate(anno = ifelse(p_val_adj < 0.05,
                       module,
                       "")) %>%
  ggplot(aes(x = avg_log2FC,
             y = -log10(p_val_adj)))  +
  geom_rect(aes(xmin = -Inf,
                xmax = Inf,
                ymin = -Inf,
                ymax = -log10(0.05)),
            fill = "grey75",
            alpha = 0.8,
            color = NA) +
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "grey75",
             alpha = 0.8) +
  geom_point(aes(fill = module),
             size = 8,
             pch = 21,
             color = "black") +
  geom_text_repel(aes(label = anno),
                  color = "black",
                  min.segment.length = 0,
                  max.overlaps = Inf,
                  size = 8,
                  point.padding = 4,
                  nudge_x = 0.06) +
  theme_classic() +
  theme(legend.position = 'none') +
  xlim(c(-0.45,
         0.45)) +
  scale_fill_manual(values = mod_colors)+
  xlab(bquote("Average log"[2] ~ "(Fold Change)")) +
  ylab(bquote("-log"[10] ~ "(Adj. P-value)")) +
  ggtitle('Differential module eigengene analysis') +
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1))+
  theme(axis.text = element_text(size = 15))  +
  theme(axis.title = element_text(size = 20))+
  theme(plot.title = element_text(size=20))
ggsave('neurons/wgcna/ME by cluster volcanoplot poster.pdf',
       width = 10,
       height = 5,
       units = "in",
       dpi = 320)

### graph MEs for relevant clusters
for (i in mods.neurons) {
  #create list of clusters with module eigengene above:
  # average scaled expression of 1 and 50% percent expression
  cluster.list = burtoni.neuron.wgcna.module.cluster %>% 
    filter(features.plot == i) %>% 
    mutate(keep = ifelse(avg.exp.scaled >= 1 & pct.exp >= 50,
                         "keep",
                         "remove")) %>% 
    filter(keep == "keep") %>%
    mutate(id = as.character(id)) %>% 
    pull(id) 
  
  ##graph ME per cluster across social status
  burtoni.snseq.combined.sct.all.neurons.recluster@meta.data %>% 
    mutate(cluster = as.character(integrated_snn_res.0.8)) %>% 
    filter(cluster %in% cluster.list) %>% 
    droplevels() %>%
    ggplot(aes_string(x = 'integrated_snn_res.0.8',
                      y = i,
                      color = 'orig.ident')) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(position=position_jitterdodge(jitter.width = 0.25,
                                             jitter.height = 0),
               alpha = 0.2)+
    theme_classic() +
    ggtitle(i)
  ggsave(paste('neurons/wgcna/modules/',
               i,
               ' ME by cluster boxplot.png',
               sep = ''),
         width = 10,
         height = 10)
  
  ## ME across social status
  burtoni.snseq.combined.sct.all.neurons.recluster@meta.data %>% 
    mutate(cluster = as.character(integrated_snn_res.0.8)) %>% 
    filter(cluster %in% cluster.list) %>% 
    droplevels() %>%
    ggplot(aes_string(x = 'orig.ident',
                      y = i)) +
    geom_violin(draw_quantiles = c(0.5)) +
    geom_jitter(width = 0.25,
                height = 0)+
    theme_classic() +
    ggtitle(i)
  ggsave(paste('neurons/wgcna/modules/',
               i,
               ' ME violinplot.png',
               sep = ''),
         width = 10,
         height = 10)
  
  # cluster bias 
  burtoni.snseq.combined.sct.all.neurons.recluster@meta.data %>% 
    mutate(cluster = as.character(integrated_snn_res.0.8)) %>% 
    filter(cluster %in% cluster.list) %>% 
    mutate(keep = ifelse(get(i) >= 1,
                         "present",
                         "absence")) %>% 
    droplevels() %>%
    select(c(cluster,
             keep,
             orig.ident)) %>% 
    table() %>% 
    as.data.frame() %>% 
    pivot_wider(id_cols = c(cluster,
                            keep),
                names_from = orig.ident,
                values_from = Freq) %>% 
    filter(keep == "present") %>% 
    mutate(dom_burtoni_snseq.scaled = 0.783645 * dom_burtoni_snseq) %>% 
    ggplot(aes(x = dom_burtoni_snseq.scaled,
               y = sub_burtoni_snseq,
               label = cluster)) +
    geom_abline(slope = 1,
                intercept = 0) +
    geom_text() +
    theme_classic() +
    ggtitle(paste(i,
                  'sub vs dom scaled counts per cluster'))
  ggsave(paste('neurons/wgcna/modules/',
               i,
               ' ME count bias.png',
               sep = ''),
         width = 10,
         height = 10)
  
}

#### marker genes, hub genes, and go terms WGCNA for neurons old ####
hub.genes.neurons = GetModules(burtoni.snseq.combined.sct.all.neurons.recluster)

write_csv(hub.genes.neurons,
          './neurons/wgcna/hub.genes.neurons.csv')

## get list of marker genes
write_csv(markers.neuron,
          './neurons/markers.neuron.csv')

# get go terms for cluster 7

# save hub genes for GO analysis with:
# to get GO term accession use biomart: http://useast.ensembl.org/biomart/martview/

## load module go term data
# brown module
brown_go_terms = read.csv('./neurons/wgcna/brown_go_terms.csv')
# black module
black_go_terms = read.csv('./neurons/wgcna/black_go_terms.csv')

# marker neurons
markers.neuron_go_terms = read.csv('./neurons/markers.neuron_go_terms.csv')

## combine go terms with kME scores
#brown
brown_go_terms.kme = brown_go_terms %>% 
  left_join(hub.genes.neurons %>% 
              filter(module == 'brown') %>% 
              select(c(gene_name,
                       module,
                       kME_brown)),
            by = c("Gene.name" = "gene_name")) %>% 
  filter(!is.na(kME_brown)) %>% 
  rbind(brown_go_terms %>% 
          left_join(hub.genes.neurons %>% 
                      filter(module == 'brown') %>% 
                      select(c(gene_name,
                               module,
                               kME_brown)),
                    by = c("Gene.stable.ID" = "gene_name")) %>% 
          filter(!is.na(kME_brown)))

#black
black_go_terms.kme = black_go_terms %>% 
  left_join(hub.genes.neurons %>% 
              filter(module == 'black') %>% 
              select(c(gene_name,
                       module,
                       kME_black)),
            by = c("Gene.name" = "gene_name")) %>% 
  filter(!is.na(kME_black)) %>% 
  rbind(black_go_terms %>% 
          left_join(hub.genes.neurons %>% 
                      filter(module == 'black') %>% 
                      select(c(gene_name,
                               module,
                               kME_black)),
                    by = c("Gene.stable.ID" = "gene_name")) %>% 
          filter(!is.na(kME_black)))

## combine marker go terms with pval
#marker
markers.neuron_go_terms.sig = markers.neuron_go_terms %>% 
  left_join(markers.neuron,
            by = c("Gene.name" = "gene"))  %>% 
  filter(!is.na(p_val)) %>% 
  rbind(markers.neuron_go_terms %>% 
          left_join(markers.neuron,
                    by = c("Gene.stable.ID" = "gene"))) %>% 
  filter(!is.na(p_val)) %>% 
  relocate(GO.term.accession, 
           .after = last_col())%>% 
  relocate(p_val_adj, 
           .after = last_col())

## increase kME to positive values
# move GO terms and kME to end
#brown
brown_go_terms.kme = brown_go_terms.kme %>% 
  mutate(kME_brown.scaled = kME_brown + 
           abs(min(kME_brown)) + 
           0.0001) %>% 
  relocate(GO.term.accession, 
           .after = last_col())%>% 
  relocate(kME_brown.scaled, 
           .after = last_col())

#black
black_go_terms.kme = black_go_terms.kme %>% 
  mutate(kME_black.scaled = kME_black + 
           abs(min(kME_black)) + 
           0.0001) %>% 
  relocate(GO.term.accession, 
           .after = last_col())%>% 
  relocate(kME_black.scaled, 
           .after = last_col())

## create list of top 25 genes
#brown
brown_go_terms.kme.top25 = brown_go_terms %>% 
  left_join(hub.genes.neurons %>% 
              filter(module == 'brown') %>% 
              top_n(25,
                    kME_brown) %>% 
              select(c(gene_name,
                       module,
                       kME_brown)),
            by = c("Gene.name" = "gene_name")) %>% 
  filter(!is.na(kME_brown)) %>% 
  rbind(brown_go_terms %>% 
          left_join(hub.genes.neurons %>% 
                      filter(module == 'brown') %>% 
                      top_n(25,
                            kME_brown) %>% 
                      select(c(gene_name,
                               module,
                               kME_brown)),
                    by = c("Gene.stable.ID" = "gene_name")) %>% 
          filter(!is.na(kME_brown)))  %>% 
  relocate(GO.term.accession, 
           .after = last_col())%>% 
  relocate(kME_brown, 
           .after = last_col())

#black
black_go_terms.kme.top25 = black_go_terms %>% 
  left_join(hub.genes.neurons %>% 
              filter(module == 'black') %>% 
              top_n(25,
                    kME_black) %>% 
              select(c(gene_name,
                       module,
                       kME_black)),
            by = c("Gene.name" = "gene_name")) %>% 
  filter(!is.na(kME_black)) %>% 
  rbind(black_go_terms %>% 
          left_join(hub.genes.neurons %>% 
                      filter(module == 'black') %>% 
                      top_n(25,
                            kME_black) %>% 
                      select(c(gene_name,
                               module,
                               kME_black)),
                    by = c("Gene.stable.ID" = "gene_name")) %>% 
          filter(!is.na(kME_black))) %>% 
  relocate(GO.term.accession, 
           .after = last_col())%>% 
  relocate(kME_black, 
           .after = last_col())

## create list for revigo
#brown
brown_go_terms.kme.top25.filter = brown_go_terms.kme.top25 %>% 
  select(c(GO.term.accession,
           kME_brown)) %>% 
  filter(GO.term.accession != "") %>% 
  distinct() 

#black
black_go_terms.kme.top25.filter = black_go_terms.kme.top25 %>% 
  select(c(GO.term.accession,
           kME_black)) %>% 
  filter(GO.term.accession != "") %>% 
  distinct() 



# to look for enrichment and create figures use GO term accession and kME (higher is better) for each module: http://revigo.irb.hr/
# remember to sort by kME and filter out blank GO terms before copying to revigo
# all module genes
#brown
write_csv(brown_go_terms.kme,
          './neurons/wgcna/go_terms/brown_go_terms.kme.csv')
#black
write_csv(black_go_terms.kme,
          './neurons/wgcna/go_terms/black_go_terms.kme.csv')

#top 25 hub genes filtered
#brown
write_csv(brown_go_terms.kme.top25.filter,
          './neurons/wgcna/go_terms/brown_go_terms.kme.top25.filter.csv')
#black
write_csv(black_go_terms.kme.top25.filter,
          './neurons/wgcna/go_terms/black_go_terms.kme.top25.filter.csv')

#marker
write_csv(markers.neuron_go_terms.sig,
          './neurons/markers.neuron_go_terms.sig.csv')
#use 0.05 p_val_adj cutoff










#### create figures for poster old ####
### dotplot of broad cell types across UMAP

#assign broad cell class
burtoni.scsorter.data.scsort.output.broad = burtoni.scsorter.data.scsort.output %>% 
  select(Cell.type,
         Cell.id) %>%
  mutate(Cell.class.broad = case_when(Cell.type == 'neurons' ~ 'neurons',
                                      Cell.type == 'inhibitory' ~ 'neurons',
                                      Cell.type == 'excitatory' ~ 'neurons',
                                      Cell.type == 'mural' ~ 'vascular',
                                      Cell.type == 'VSM' ~ 'vascular',
                                      Cell.type == 'endothelial' ~ 'vascular',
                                      Cell.type == 'Unknown' ~ 'unknown',
                                      TRUE ~ 'glia'))

##combine umap with metadata
burtoni.scsorter.data.scsort.output.umap =  full_join(burtoni.scsorter.data.scsort.output.broad,
                                                      burtoni.scsorter.data.scsort.output.umap)

##graph umap with cell markers
burtoni.scsorter.data.scsort.output.umap %>% 
  ggplot(aes(x=UMAP_1,
             y= UMAP_2,
             color = Cell.type)) +
  geom_point() +
  theme_bw()
ggsave('scsorter/UMAP scsorter cell types.png',
       width = 10,
       height = 10)


#graph umap with cell markers
#broad cell class
burtoni.scsorter.data.scsort.output.umap %>% 
  ggplot(aes(x=UMAP_1,
             y= UMAP_2,
             color = Cell.class.broad)) +
  geom_point() +
  theme_classic() +
  scale_color_manual(values = c('3C51A3',
                                'orange3',
                                'grey',
                                '4AB85C'))
ggsave('scsorter/UMAP scsorter cell cell types broad.pdf',
       width = 10,
       height = 10)

#graph umap with cell markers
#broad cell class
# seperate doms and subs
burtoni.scsorter.data.scsort.output.umap %>% 
  ggplot(aes(x=UMAP_1,
             y= UMAP_2,
             color = Cell.class.broad)) +
  geom_point() +
  theme_classic() +
  scale_color_manual(values = c('3C51A3',
                                'orange3',
                                'grey',
                                '4AB85C')) +
  facet_grid( ~ orig.ident) +
  labs(color='Broad cell classes') +
  theme(legend.text=element_text(size=20),
        legend.title=element_text(size=20))+
  guides(colour = guide_legend(override.aes = list(size=10)))
ggsave('scsorter/UMAP scsorter cell cell types broad sub vs dom.pdf',
       width = 12,
       height = 10)




