#### Burtoni snseq seurat analysis
### limmatrend vascular
###Note: Seurat requires R version > 4
## use lambcomp1 to run R with command 
# > R-4.0.3



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
library(cowplot)
library(patchwork)
library(corrplot)
library(ggrepel)
library(tidyverse)

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
burtoni.sctypemarkers.hypo = read.csv('./sctype.hypo/sctypemarkers.hypo.csv')

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

### neuropeptide list
neuropeptides.list = read_csv("../Gene.lists/neuropeptides.list_orthologs.csv")

#reduced range 2
resolution.range.reduced.2 <- seq(from = 0.2, to = 1, by = 0.1)

#reduced range 2
resolution.range.reduced.2 <- seq(from = 0, to = 1, by = 0.1)

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

#### clustering for vascular ####
burtoni.snseq.combined.sct.vascular = burtoni.snseq.combined.sct

#set idents
Idents(object = burtoni.snseq.combined.sct.vascular) <- "sctypemarkers.hypo"

# check celltype names
burtoni.snseq.combined.sct.vascular@meta.data$sctypemarkers.hypo %>% 
  unique()

#subset to vascular
burtoni.snseq.combined.sct.vascular = subset(burtoni.snseq.combined.sct.vascular,
                                            idents = c("C7-7: Vascular"))

#remove large file
# rm(burtoni.snseq.combined.sct)

# need to set to integrated for clustering
DefaultAssay(burtoni.snseq.combined.sct.vascular) = 'integrated'

#check data loaded correctly
## run PCA, UMAP, and cluster 
#use 0.1 resolution
burtoni.snseq.combined.sct.vascular = burtoni.snseq.combined.sct.vascular %>% 
  RunPCA() %>%
  FindNeighbors(dims = 1:15) %>%
  RunUMAP(dims = 1:15) %>%
  FindClusters(resolution = 0.1)

## graph 
# idents to new clusters
Idents(object = burtoni.snseq.combined.sct.vascular) <- "integrated_snn_res.0.1"

DimPlot(burtoni.snseq.combined.sct.vascular,
        group.by='integrated_snn_res.0.1',
        label=TRUE) +
  ggtitle('vascular') +
  NoLegend() +
  theme_classic()
ggsave('vascular/UMAP vascular clusters.png',
       height = 5.5,
       width = 5.5)

# #check per genotype
# #use all vascular
DimPlot(burtoni.snseq.combined.sct.vascular,
        group.by='Genotype.id') +
  ggtitle('vascular')+
  theme_classic()
ggsave('vascular/UMAP vascular genotypes.png',
       height = 5.5,
       width = 5.5)

## create count across genotype per cluster
burtoni.snseq.combined.sct.vascular@meta.data %>%
  as.data.frame() %>%
  select(c(Genotype.id,
           integrated_snn_res.0.1)) %>%
  table() %>%
  as.data.frame() %>%
  mutate(orig.ident = ifelse(grepl("Dom",
                                   Genotype.id),
                             "dom",
                             "sub")) %>%
  ggplot(aes(x = integrated_snn_res.0.1,
             y = Freq,
             color = orig.ident)) +
  geom_jitter(height = 0,
              width = 0.1) +
  theme_classic()
ggsave('vascular/Genotype per clusters count.png',
       height = 5.5,
       width = 5.5)

## project clusters back onto umap
# get original UMAP data
UMAP.orig = burtoni.snseq.combined.sct[["umap"]]@cell.embeddings %>% 
  as.data.frame() %>% 
  rownames_to_column('Cell.id')

  
# add vascular metadata
UMAP.orig = UMAP.orig %>% 
  full_join(burtoni.snseq.combined.sct.vascular@meta.data %>%
                        as.data.frame() %>% 
                        rownames_to_column('Cell.id') %>% 
                        select(c(Cell.id,
                                 integrated_snn_res.0.1))) 

# graph
UMAP.orig %>% 
  ggplot(aes(x = UMAP_1,
             y= UMAP_2)) +
  geom_point(color = 'grey') +
  geom_point(aes(color = integrated_snn_res.0.1)) +
  theme_classic() +
  labs(color = 'vascular clusters')
ggsave('vascular/UMAP projection vascular clusters.png',
       height = 10,
       width = 12.5)

### clustree
library(clustree)

# cluster across resolutions
burtoni.snseq.combined.sct.vascular.clustree <- Seurat::FindClusters(object = burtoni.snseq.combined.sct.vascular,
                                                                               resolution = resolution.range.reduced.2)
#check data
# head(burtoni.snseq.combined.sct.vascular.recluster.clustree[[]])
# #set presentation colors
# presentation.color <- c('#66c2a5',
#                         '#fc8d62',
#                         '#8da0cb',
#                         '#e78ac3',
#                         '#a6d854',
#                         '#ffd92f',
#                         '#e5c494',
#                         '#b3b3b3')
#clustree
clustree(burtoni.snseq.combined.sct.vascular.clustree,
         prefix = "integrated_snn_res.",
         node_colour = 'cluster',
         node_size_range = c(10,20),
         scale_node_text = TRUE) +
  scale_edge_color_continuous(low = "black",
                              high = "black") +
  # scale_color_manual(values = presentation.color)+
  theme(legend.position = "bottom")
ggsave('vascular/Clustree vascular.png',
       height = 10,
       width = 10)

# remove large file
rm(burtoni.snseq.combined.sct.vascular.clustree)

#### limmatrend genotype, clusters, and status approach with SCT treat prep ####
### calculate variable genes
## identify top 1500 variable genes
# use integrated assay for variable features
burtoni.snseq.combined.sct.vascular.variable.prep <- FindVariableFeatures(burtoni.snseq.combined.sct.vascular, 
                                                                                              assay = 'integrated',
                                                                                              selection.method = "vst", 
                                                                                              nfeatures = 1500, 
                                                                                              verbose = F)
# identify top 1500 variable genes
vascular.topgenes.prep <- head(VariableFeatures(burtoni.snseq.combined.sct.vascular.variable.prep,
                                                               assay = 'integrated'), 
                                              1500)

# create dummy
burtoni.snseq.combined.sct.vascular.expression.sct.prep = full_join(full_join(burtoni.snseq.combined.sct.vascular@reductions$umap@cell.embeddings %>% 
                                                                                                   as.data.frame() %>% 
                                                                                                   rownames_to_column("Cell.id"),
                                                                                burtoni.snseq.combined.sct.vascular@meta.data %>%
                                                                                                   rownames_to_column("Cell.id")),
                                                                      burtoni.snseq.combined.sct.vascular@assays$SCT@data %>% 
                                                                                         as.data.frame()  %>% 
                                                                                         t() %>% as.data.frame() %>% 
                                                                                         rownames_to_column('Cell.id'))

## create vector of factor
vascular.DomvsSub.vector.list.sct.prep = burtoni.snseq.combined.sct.vascular.expression.sct.prep %>% 
  mutate(vascular.orig.ident = orig.ident %>% 
           as.factor()) %>% 
  pull(vascular.orig.ident) %>% 
  droplevels()

### counts matrix
## raw read count matrix
## rows = genes, columns = cells
# only keep 1500 variable genes
# set negative values to 0
vascular.DomvsSub.vector.count.sct.prep = GetAssayData(burtoni.snseq.combined.sct.vascular,
                                                                assay = 'SCT') %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column('gene') %>% 
  dplyr::select(c(gene,
                  burtoni.snseq.combined.sct.vascular.expression.sct.prep %>% 
                    pull(Cell.id))) %>% 
  filter(gene %in% vascular.topgenes.prep) %>% 
  column_to_rownames('gene') %>% 
  as.matrix() %>% 
  pmax(0)


## run with vascular clusters
#need list with L with  count and condt
# add in cell cluster with condt.2
## create vector of factor
vascular.DomvsSub.vector.list.cluster.sct.prep = burtoni.snseq.combined.sct.vascular.expression.sct.prep %>% 
  mutate(vascular.cluster = integrated_snn_res.0.1 %>% 
           as.factor()) %>% 
  pull(vascular.cluster) %>% 
  droplevels


# add in genotype
## create vector of factor
vascular.DomvsSub.vector.list.genotype.sct.prep = burtoni.snseq.combined.sct.vascular.expression.sct.prep %>% 
  mutate(genotype = Genotype.id %>% 
           as.factor()) %>% 
  pull(genotype) %>% 
  droplevels()

#create list
vascular.vector.limma.sct.prep = list(count = vascular.DomvsSub.vector.count.sct.prep,
                                               condt.2 = vascular.DomvsSub.vector.list.cluster.sct.prep,
                                               genotype = vascular.DomvsSub.vector.list.genotype.sct.prep)
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
  # # Open pdf file
  # pdf(file= "./vascular/limmatrend/limmatrend.histograms.pdf" )
  # # create a 2X2 grid
  # par( mfrow= c(2,2) )
  # #graph
  # hist(tt0$P.Value, 50)
  # hist(tt0$adj.P.Val, 50)
  # hist(tt1$P.Value, 50)
  # hist(tt1$adj.P.Val, 50)
  # hist(tt2$P.Value, 50)
  # hist(tt2$adj.P.Val, 50)
  # hist(tt3$P.Value, 50)
  # hist(tt3$adj.P.Val, 50)
  # hist(tt4$P.Value, 50)
  # hist(tt4$adj.P.Val, 50)
  # hist(ttDS$P.Value, 50)
  # hist(ttDS$adj.P.Val, 50)
  # dev.off()
  # 
  # # Open pdf file
  # pdf(file= "./vascular/limmatrend/limmatrend.MDS.pdf" )
  # # create a 2X1 grid
  # par( mfrow= c(2,1) )
  # limma::plotMDS(dge, 
  #                col = as.numeric(as.factor(L$condt)), 
  #                pch = 19)
  # plotMD(fit)
  # limma::plotMDS(dge, 
  #                col = as.numeric(as.factor(L$condt.2)), 
  #                pch = 19)
  # plotMD(fit)
  # dev.off()
  # 
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
vascular.genotype.vector.limma.results.sct = run_limmatrend_complex_genotype.sct(vascular.vector.limma.sct.prep)

# save results to dataframe
vascular.genotype.vector.limma.results.df.sct.treat.prep = full_join(vascular.genotype.vector.limma.results.sct$tt0 %>% 
                                                                                rename_with(~paste0(.,"_DvsS_0")) %>% 
                                                                                rownames_to_column("Gene"),
                                                                              vascular.genotype.vector.limma.results.sct$tt1 %>% 
                                                                                rename_with(~paste0(.,"_DvsS_1")) %>% 
                                                                                rownames_to_column("Gene")) %>% 
  full_join(vascular.genotype.vector.limma.results.sct$tt2 %>% 
              rename_with(~paste0(.,"_DvsS_2")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(vascular.genotype.vector.limma.results.sct$tt3 %>% 
              rename_with(~paste0(.,"_DvsS_3")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(vascular.genotype.vector.limma.results.sct$tt4 %>% 
              rename_with(~paste0(.,"_DvsS_4")) %>% 
              rownames_to_column("Gene")) %>% 
  full_join(vascular.genotype.vector.limma.results.sct$ttDS %>% 
              rename_with(~paste0(.,"_DvsS")) %>% 
              rownames_to_column("Gene")) 



#add color for significance  
vascular.genotype.vector.limma.results.df.sct.treat.prep = vascular.genotype.vector.limma.results.df.sct.treat.prep %>% 
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
  vascular.genotype.vector.limma.results.df.sct.treat.prep %>% 
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
  ggsave(paste0('vascular/limmatrend/limma.vascular.volcano.', i, '.png'),
         width = 5,
         height = 5)
}

vascular.genotype.vector.limma.results.df.sct.treat.prep = vascular.genotype.vector.limma.results.df.sct.treat.prep %>% 
  mutate(Keep = case_when(Sig_DvsS_0 != "Not Sig" ~ "Keep",
                          Sig_DvsS_1 != "Not Sig" ~ "Keep",
                          Sig_DvsS_2 != "Not Sig" ~ "Keep",
                          Sig_DvsS_3 != "Not Sig" ~ "Keep",
                          Sig_DvsS_4 != "Not Sig" ~ "Keep",
                          Sig_DvsS != "Not Sig" ~ "Keep",
                          TRUE ~ "Remove")) %>% 
  # mutate(Gene = case_when(Gene == "ENSONIG00000040658" ~ "HBE1",
  #                         Gene == "ENSONIG00000004376" ~ "EBF4",
  #                         Gene == "ENSONIG00000036497" ~ "FBXL19",
  #                         Gene == "ENSONIG00000010896" ~ "PBX3",
  #                         Gene == "si:dkey-22o22.2" ~ "CDH2",
  #                         Gene == "ENSONIG00000034988" ~ 'MT: 1,086',
  #                         Gene == "ENSONIG00000010896" ~ 'Pbx3',
  #                         Gene == "ENSONIG00000002603" ~ 'GnRHR-II',
  #                         Gene == "ENSONIG00000031366" ~ 'MT: 70',
  #                         Gene == "ENSONIG00000004376" ~ "Ebf4",
  #                         Gene == "ENSONIG00000003414" ~ 'slc47a2.1',
#                         Gene == "ENSONIG00000003310" ~ 'Epha5',
#                         Gene == "im:7142702"        ~ 'Bcl11b', 
#                         Gene == "ENSONIG00000040640"  ~ 'ZNF804A',      
#                         Gene == "ENSONIG00000040652" ~ 'or106-8',
#                         Gene == "ENSONIG00000041695" ~ 'CCBE1',
#                         Gene == "ENSONIG00000036306" ~ '36306',
#                         TRUE ~ Gene)) %>% 
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
                       TRUE ~ 0))
#### save results ####
# save to csv
write_csv(vascular.genotype.vector.limma.results.df.sct.treat.prep,
          file = 'vascular/limmatrend/vascular.limmatrend.results.csv')

##### poster heatmap ####
library(pheatmap)

# create matrix
vascular.deg.matrix = vascular.genotype.vector.limma.results.df.sct.treat.prep %>%
  filter(Keep == "Keep") %>% 
  column_to_rownames('Gene') %>% 
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

## replace 'inf' values with max 
vascular.deg.matrix[!is.finite(vascular.deg.matrix)] <- vascular.deg.matrix %>% 
  as.data.frame() %>% 
  filter_all(all_vars(is.finite(.))) %>% 
  max()


# graph pheatmap
#pdf
vascular.deg.matrix %>% 
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
           filename = "vascular/limmatrend/Heatmap DEGs across clusters.pdf",
           width = 12,
           height = 5
  )
#png
vascular.deg.matrix %>% 
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
           filename = "vascular/limmatrend/Heatmap DEGs across clusters.png",
           width = 12,
           height = 5
  )

# create heatmap
vascular.deg.phtmap <- pheatmap::pheatmap(vascular.deg.matrix,
                                            cluster_rows = F,
                                            cluster_cols = T,
                                            scale = 'none',
                                            silent = T)

# get order list
vascular.deg.order.list = c(vascular.deg.phtmap$tree_col$order)
# create label list
vascular.deg.order = data.frame(label = vascular.deg.phtmap$tree_col$labels) %>% 
  rownames_to_column('position')
# combine labels and position
vascular.deg.order = vascular.deg.order[match(vascular.deg.order.list,
                                                  vascular.deg.order$position),] %>% 
  mutate(vascular.deg.order.list = 1:length(vascular.deg.order.list)) %>% 
  mutate(position = as.numeric(position))
# create new order list
vascular.deg.order.list.new = vascular.deg.order %>%
  arrange(position) %>% 
  pull(label)

# move Dom bias genes to end
vascular.deg.col_dend <- vascular.deg.phtmap[[2]]
vascular.deg.col_dend <- dendextend::rotate(vascular.deg.col_dend, 
                                              order = vascular.deg.order.list.new)
# graph pheat map
# pdf
pheatmap(vascular.deg.matrix, 
         cluster_cols=as.hclust(vascular.deg.col_dend),
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
         filename = "vascular/limmatrend/Heatmap DEGs across clusters reorder.pdf",
         width = 12,
         height = 5
)

# png
pheatmap(vascular.deg.matrix, 
         cluster_cols=as.hclust(vascular.deg.col_dend),
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
         filename = "vascular/limmatrend/Heatmap DEGs across clusters reorder.png",
         width = 12,
         height = 5
)
