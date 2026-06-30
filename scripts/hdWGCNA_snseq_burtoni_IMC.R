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

# # create new conda environment for R
# conda create -n hdWGCNA -c conda-forge -c bioconda r-base=4.4 mamba
# 
# # activate conda environment
# conda activate hdWGCNA
# 
# # install critical R packages
# mamba install -c conda-forge -c bioconda r-seurat r-hdf5r r-wgcna r-igraph r-tidyverse r-ggraph r-harmony r-enrichr r-devtools
# 
# # install Bioconductor packages
# mamba install -c conda-forge -c bioconda bioconductor-ucell bioconductor-genomicranges bioconductor-geneoverlap 
# Next open R and install hdWGCNA.
# 
# # install Bioconductor
# install.packages("BiocManager")
# BiocManager::install()
# 
# # install hdWGCNA from GitHub
# devtools::install_github('smorabit/hdWGCNA', ref='dev')


# # need to reinstall hdWGCNA
# install_github("wjawaid/enrichR")
# library(remotes)
# install.packages("spam")
# remotes::install_github("carmonalab/UCell", ref="v2.2")
# devtools::install_github("smorabit/hdWGCNA")

# devtools::install_github("neurorestore/Libra")

#### load libraries ####

#load libraries
library(gtable, lib.loc = "/usr/lib/R/site-library")
library(cowplot)
library(sp)
library(SeuratObject, lib.loc = "/usr/local/lib/R/site-library")
library(Seurat)
library(GeneOverlap)
library(igraph)
library(WGCNA)
library(hdWGCNA)
library(patchwork)
library(corrplot)
library(ggrepel)
library(tidyverse)
library(emmeans)
library(multcomp)
library(ggtree)

# DEG
library(Libra)
library(UpSetR)
library(ComplexHeatmap)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

#### load data ####
### check mito counts
gtf.dataframe = rtracklayer::import("/stor/work/Hofmann/All_projects/A_burtoni_snseq/cellranger/ref/NileTilapia.reference/genes/genes.gtf.gz") %>% 
  as.data.frame()

# # check for mito genes
# gtf.dataframe %>% 
#   filter(seqnames == "MT") %>%
#   filter(type == 'gene') %>% 
#   View()
# # 37 MT genes

# get list of mito genes 
mito.genes.list = gtf.dataframe %>%
  filter(seqnames == "MT") %>%
  filter(type == 'gene') %>%
  mutate(gene_ID_name = ifelse(is.na(gene_name),
                               gene_id,
                               gene_name)) %>% 
  pull(gene_ID_name) 

# remove gtf dataframe
rm(gtf.dataframe)

### single cell data
# load('burtoni.snseq.combined.sct.all.RData')
load('burtoni.snseq.combined.sct.RData')

# ### scsorter data
# load("burtoni.scsorter.data.scsort.output.RData")

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
# ## cell type
# burtoni.snseq.combined.sct = AddMetaData(
#   object = burtoni.snseq.combined.sct,
#   metadata = burtoni.scsorter.data.scsort.output %>% 
#     dplyr::select(Cell.type,
#            Cell.id) %>% 
#     column_to_rownames(var = "Cell.id"),
#   col.name = 'Cell.type'
# )

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

#### GLM binomial all ####
### get summary stats across genotypes
burtoni.snseq.combined.sct@meta.data %>% 
  dplyr::select(orig.ident,
                Genotype.id) %>% 
  table() %>% 
  as.data.frame() %>% 
  filter(Freq > 0) %>% 
  dplyr::select(-c(Genotype.id)) %>% 
  psych::describeBy(group = 'orig.ident')

### create tree of cell types
burtoni.snseq.combined.sct.tree <- BuildClusterTree(
  burtoni.snseq.combined.sct,
  dims = 1:30,
  reorder = FALSE,
  reorder.numeric = FALSE,
  slot = 'data',
  assay = "SCT"
)

## create tree
tree <- burtoni.snseq.combined.sct.tree@tools$BuildClusterTree
# tree$tip.label <- paste0("Cluster ", tree$tip.label)

# convert tree to tibble
tree = as_tibble(tree)

## add cell type data
tree = full_join(tree,
                 burtoni.snseq.combined.sct@meta.data %>% 
                   dplyr::select(seurat_clusters,
                                 sctypemarkers.hypo) %>% 
                   distinct(),
                 by = c("label" = "seurat_clusters"))

# convert back to tree
# paper
tree = treeio::as.treedata(tree)
ggtree(tree) +
  scale_y_reverse() +
  geom_tree() +
  theme_tree() +
  geom_tiplab(offset = 1.5,
              angle = 90,
              vjust = 0.5,
              hjust =0.5,
              size = 3) +
  geom_tippoint(aes(color = sctypemarkers.hypo),
                        shape = 16,
                        size = 5) +
  coord_cartesian(clip = 'off') +
  theme(plot.margin = unit(c(0,0.5,0,0),
                           'in')) +
  scale_color_manual(values = c('#F8766D',
                                '#B79F00',
                                '#00BA38',
                                '#00BFC4',
                                '#619CFF',
                                '#F564E3',
                                'grey')) +
  theme(legend.position = 'none') 
ggsave('sctype.hypo/unknown/Cluster cell type tree.pdf',
       height = 6.5,
       width = 3,
       units = 'in',
       dpi = 720,
       limitsize = FALSE)

# remove extra 
rm(tree)
rm(burtoni.snseq.combined.sct.tree)


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

# paper
cluster.genotype.count %>% 
  ggplot() +
  geom_point(aes(x = seurat_clusters,
                 y = Percent,
                 group = Status,
                 color = Status),
             position = position_dodge(width = 0.5),
             size = 5) +
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
          size = 8) +
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
  xlab('clusters') +
theme(legend.position = 'inside',
      legend.position.inside = c(0.9,
                                 0.75),
      text = element_text(size = 16))
ggsave('neurons/neuron_cluster/clusters binomial GLM cell percent.pdf',
       height = 6,
       width = 13,
       units = 'in',
       dpi = 720)

# paper
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
             size = 5) +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = 0.5,
             linetype = 2) +
  geom_vline(xintercept = -0.5,
             linetype = 2) +
  # geom_text(aes(label = seurat_clusters),
  #           color = 'white',
  #           size = 3) +
  labs(x = "Effect size", 
       y = "all clusters") +
  theme_classic() +
  scale_color_manual(breaks = c('dom_bias',
                                'sub_bias',
                                'none'),
                     values = c(
                       "#4e499e",
                       "#60bb46",
                       "grey")) +
  xlim(-max.diff.value-0.25,
       max.diff.value+0.25) +
  # ggtitle("All clusters: binomial GLM cell count") +
  theme(legend.position = 'none')+
  theme(axis.text = element_text(size = 16),
        axis.title = element_text(size = 16))+
  scale_x_reverse()
ggsave('neurons/neuron_cluster/clusters binomial GLM forest plot.pdf',
       height = 12,
       width = 6,
       units = 'in',
       dpi = 720)


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



#### DEG analysis all ####
# create dummy dataset
all.obj = burtoni.snseq.combined.sct

# use counts assay for limmatrend
# needed for pseudobulking to renormalize
DefaultAssay(all.obj) = "RNA"

### pre-filter genes
## remove genes that are not in 10% of cells within any cluster and treatment
all.obj$cluster_orig.ident = paste0(all.obj$seurat_clusters, 
                                       "_", 
                                       all.obj$orig.ident)

# get count data
all.obj.counts = GetAssayData(all.obj, 
                                 assay = "RNA", 
                                 slot = "counts")

# get list of cluster by treatment
all.obj.groups = all.obj$cluster_orig.ident

# get percent of cells with > 0 reads per group
all.obj.pct = sapply(unique(all.obj.groups), function(g) {
  cells = which(all.obj.groups == g)
  rowMeans(all.obj.counts[, cells, drop = FALSE] > 0) * 100
})

# get max percent value per gene
all.obj.pct.max = apply(all.obj.pct, 1, max)

# keep genes that have a max of at least 10% 
all.obj.pct.keep = names(all.obj.pct.max[all.obj.pct.max >= 10])
# down to 7208 genes

# remove mito genes
all.obj.pct.keep = all.obj.pct.keep[!(all.obj.pct.keep %in% mito.genes.list)]
# down to 7194 genes

# subset seurat object
all.obj = subset(all.obj, 
                    features = all.obj.pct.keep)
# get an error but appears to have worked
Features(all.obj) %>%
  length()
# [1] 7194

### use run_de to get DEG across clusters using pseudobulk
## get cell information
# set cell ID info that run_de will recognize
all.obj$cell_type = all.obj$seurat_clusters
all.obj$replicate = all.obj$Genotype.id
# make sure order is sub vs dom 
all.obj$label = factor(all.obj$orig.ident, levels = c("sub_burtoni_snseq", 
                                                            "dom_burtoni_snseq"))

## run limma trend pseudobulk
all.deg.results = run_de(
  all.obj,
  de_family = "pseudobulk",
  de_method = "limma",
  de_type   = "trend",
  min_reps = 3)

# # remove large dummy seurat object
# rm(all.obj)

## filter 
# needs to be in at least 10% of cluster in at least one group
all.deg.results = all.deg.results %>%
  filter(sub_burtoni_snseq.pct >= 0.1 | dom_burtoni_snseq.pct >= 0.1) %>%
  group_by(cell_type) %>%
  mutate(p_val_adj = p.adjust(p_val, method = "fdr")) %>%
  ungroup()

# save results
write.csv(all.deg.results,
          'DEG/all.deg.results.csv',
          row.names = F)


### compare DEG to proportion shift
## get data
# get counts data
all.deg.counts = GetAssayData(all.obj, 
                                 assay = "RNA",
                                 slot = "counts")
# get metadata
all.deg.meta   <- all.obj@meta.data

# create genotype by cluster ID
all.deg.group.id = paste(all.deg.meta$seurat_clusters,
                            all.deg.meta$Genotype.id, 
                            sep = "_") 


# get percentage of cells expressing gene 
# threshold RNA count of > 2
all.deg.pct.list <- lapply(unique(all.deg.group.id), function(g) {
  cells <- which(all.deg.group.id == g)
  
  data.frame(
    gene     = rownames(all.deg.counts),
    group    = g,
    pct_expr = rowMeans(all.deg.counts[, cells, drop = FALSE] > 2)
  )
})

# convert to data frame
all.deg.pct.list.df = do.call(rbind, 
                                 all.deg.pct.list)

# add meta data
all.deg.pct.list.df = all.deg.pct.list.df %>%
  separate(group,
           into = c("seurat_clusters", 
                    "Genotype.id"), 
           sep = "_") %>%
  left_join(
    all.deg.meta %>% 
      dplyr::select(seurat_clusters,
                    Genotype.id, 
                    orig.ident) %>%
      distinct(),
    by = c("seurat_clusters", 
           "Genotype.id"))


# get summary 
all.deg.pct.list.df.summary = all.deg.pct.list.df %>%
  group_by(gene,
           seurat_clusters,
           orig.ident) %>%
  summarise(mean_pct = mean(pct_expr), 
            .groups = "drop")

# get difference between dom and sub percentage
all.deg.pct.list.df.summary.wide = all.deg.pct.list.df.summary %>%
  tidyr::pivot_wider(names_from = orig.ident, 
                     values_from = mean_pct) %>%
  mutate(delta_pct = dom_burtoni_snseq - sub_burtoni_snseq)


# add to DEG results
all.deg.results = all.deg.results %>%
  left_join(all.deg.pct.list.df.summary.wide,
            by = c("gene", 
                   "cell_type" = "seurat_clusters"))

# look for concordance
# use 10% percentage difference as threshold
all.deg.results = all.deg.results %>%
  mutate(concordance = ifelse(sign(avg_logFC) == sign(delta_pct),
                              'concordant',
                              'discordant'),
         Likely_DEG.type = case_when(
           abs(delta_pct) > 0.2 & concordance == 'concordant' ~ "Fraction-driven",
           abs(delta_pct) <= 0.2 & concordance == 'concordant' ~ "Intensity-driven",
           abs(delta_pct) <= 0.2 & concordance == 'discordant' ~ "Intensity-driven",
           TRUE ~ "Mixed"
         ),
         DEG.type = ifelse(avg_logFC > 0,
                           'Dom',
                           'Sub'))

# add cell type
all.deg.results = all.deg.results %>% 
  left_join(all.deg.meta %>% 
              distinct(cell_type,
                       sctypemarkers.hypo))
# save results
write.csv(all.deg.results,
          'DEG/all.deg.results.csv',
          row.names = F)

### graph
## compare logFC to DEG
# scatterplot
all.deg.results %>% 
  filter(p_val_adj < 0.05) %>% 
  ggplot(aes(x = delta_pct, 
             y = avg_logFC, 
             shape = Likely_DEG.type,
             color = DEG.type)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0.2,
             linetype = "dashed") +
  geom_vline(xintercept = -0.2,
             linetype = "dashed") +
  geom_vline(xintercept = 0) +
  geom_point(size = 0.8) +
  facet_wrap(~ cell_type, scales = "free") +
  theme_classic() +
  scale_color_manual(values = c('Sub' = '#60bb46',
                                'Dom' = '#4e499e'))+
  labs(
    x = "Percent expression difference (Dom - Sub)",
    y = "Log2 Fold Change (Pseudobulk)",
    shape = "Likely DEG type",
    color = 'Social status'
  )
ggsave('DEG/DEG logfc vs percentage all.png')

# barplot
all.deg.results %>% 
  filter(p_val_adj < 0.05) %>% 
  group_by(cell_type) %>% 
  dplyr::count(Likely_DEG.type,
               DEG.type) %>% 
  mutate(n = ifelse(DEG.type == 'Sub',
                    -n,
                    n)) %>%
  ggplot(aes(x = n, 
             y = as.factor(cell_type), 
             fill = DEG.type,
             color = Likely_DEG.type,
             alpha = Likely_DEG.type)) +
  geom_bar(stat = 'identity') +
  theme_classic() +
  scale_fill_manual(values = c('Sub' = '#60bb46',
                               'Dom' = '#4e499e')) +
  scale_color_manual(values = c('Fraction-driven' = 'black',
                                'Intensity-driven' = 'black')) +
  scale_alpha_manual(values = c('Fraction-driven' = 0.6,
                                'Intensity-driven' = 1)) +
  labs(
    y = "Cluster",
    x = "DEG count",
    color = "Likely DEG type",
    fill = 'Social status'
  ) +
  theme(legend.position = 'inside',
        legend.position.inside = c(0.85,0.5))
ggsave('DEG/DEG counts per cluster all.png')

# paper
all.deg.results %>% 
  filter(!(sctypemarkers.hypo %in% c("C7-1: GLU",
                                     "C7-2: GABA"))) %>% 
  filter(p_val_adj < 0.05) %>% 
  group_by(cell_type,
           sctypemarkers.hypo) %>% 
  dplyr::count(Likely_DEG.type,
               DEG.type) %>% 
  mutate(n = ifelse(DEG.type == 'Sub',
                    -n,
                    n)) %>% 
  mutate(cell_type = as.factor(cell_type)) %>% 
  ggplot(aes(x = n, 
             y = fct_reorder(cell_type, 
                             as.integer(factor(sctypemarkers.hypo))), 
             fill = DEG.type,
             # color = Likely_DEG.type,
             alpha = Likely_DEG.type)) +
  geom_bar(stat = 'identity') +
  theme_classic() +
  scale_fill_manual(values = c('Sub' = '#60bb46',
                               'Dom' = '#4e499e')) +
  # scale_color_manual(values = c('Fraction-driven' = 'black',
  #                               'Intensity-driven' = 'black')) +
  scale_alpha_manual(values = c('Fraction-driven' = 0.6,
                                'Intensity-driven' = 1)) +
  labs(
    y = "Non-neural clusters",
    x = "DEG count",
  ) +
  theme(legend.position = 'none',
        text = element_text(size = 12),
        plot.margin = margin(l = 0.25, unit = "in"))
ggsave('DEG/DEG counts per cluster all.pdf',
       height = 3.25,
       width = 3.25,
       units = 'in',
       dpi = 720)


## graph upset
# get list of DEGs
all.significant_degs = all.deg.results %>%
  filter(!(sctypemarkers.hypo %in% c("C7-1: GLU",
                                     "C7-2: GABA"))) %>% 
  filter(p_val_adj < 0.05) %>%
  distinct(cell_type,
           gene) %>% 
  dplyr::select(cell_type,
                gene)

# create name list
all.deg_list = split(all.significant_degs$gene, 
                        all.significant_degs$cell_type)

# just graph overlap
ComplexUpset::upset(fromList(all.deg_list), 
                    intersect = names(all.deg_list),
                    min_size = 3) 
ggsave('DEG/Upset all.png')
ggsave('DEG/Upset all.pdf',
       width = 6.5,
       height = 6.5,
       units = 'in',
       dpi = 720)

# overlap
# get list of DEGs
all.significant_degs = all.deg.results %>%
  filter(p_val_adj < 0.05) %>%
  distinct(cell_type,
           gene) %>% 
  dplyr::select(cell_type,
                gene)

# create name list
all.deg_list = split(all.significant_degs$gene, 
                        all.significant_degs$cell_type)

# just graph overlap
ComplexUpset::upset(fromList(all.deg_list), 
                    intersect = names(all.deg_list),
                    sort_intersections_by = c('degree', 
                                              'cardinality'),
                    # min_size = 2,
                    min_degree = 10) 
ggsave('DEG/Upset overlap all.png')

# paper
ComplexUpset::upset(fromList(all.deg_list), 
                    intersect = names(all.deg_list),
                    sort_intersections_by = c('degree', 
                                              'cardinality'),
                    # min_size = 2,
                    min_degree = 10) 
ggsave('DEG/Upset overlap all.pdf',
       width = 6.5,
       height = 3,
       units = 'in',
       dpi = 720)


# up vs down
# get list of DEGs
all.significant_degs = all.deg.results %>%
  filter(p_val_adj < 0.05) %>%
  distinct(cell_type,
           gene,
           DEG.type) %>% 
  mutate(cell_DEG_type = paste0(cell_type,
                                '_',
                                DEG.type)) %>% 
  dplyr::select(cell_DEG_type,
                gene,) 

# create name list
all.deg_list = split(all.significant_degs$gene, 
                        all.significant_degs$cell_DEG_type)

# just graph overlap
ComplexUpset::upset(fromList(all.deg_list), 
                    intersect = names(all.deg_list),
                    min_size = 3) 
ggsave('alls/Upset all status.png')

# check genes in all clusters (or most)
all.deg.results %>%
  filter(!(sctypemarkers.hypo %in% c("C7-1: GLU",
                                     "C7-2: GABA"))) %>% 
  filter(p_val_adj < 0.05) %>%
  dplyr::count(gene) %>%
  filter(n > 11)

all.deg.results %>% 
  filter(!(sctypemarkers.hypo %in% c("C7-1: GLU",
                                     "C7-2: GABA"))) %>% 
  filter(gene %in% c(all.deg.results %>%
                       filter(!(sctypemarkers.hypo %in% c("C7-1: GLU",
                                                          "C7-2: GABA"))) %>% 
                       filter(p_val_adj < 0.05) %>% 
                       dplyr::count(gene) %>% 
                       filter(n > 11) %>% 
                       pull(gene))) %>% 
  distinct(Likely_DEG.type,
           DEG.type,
           gene) 

# Likely_DEG.type  DEG.type gene              
# <chr>            <chr>    <chr>             
#   1 Intensity-driven Dom      ENSONIG00000011023
# 2 Intensity-driven Dom      cga               
# 3 Fraction-driven  Dom      gh1               
# 4 Intensity-driven Dom      lhb               
# 5 Fraction-driven  Dom      oxt               
# 6 Fraction-driven  Dom      prl   


# # check DEG numbers
# all.deg.results %>%
#   filter(p_val_adj < 0.15) %>%
#   dplyr::count(cell_type)
# # A tibble: 20 × 2
# cell_type     n
# <chr>     <int>
#   1 0           501
# 2 1           385
# 3 10           20
# 4 11           16
# 5 12           53
# 6 13           32
# 7 14           20
# 8 15           20
# 9 16           20
# 10 17           10
# 11 18            8
# 12 19           17
# 13 2           175
# 14 3           142
# 15 4            94
# 16 5            81
# 17 6            57
# 18 7            74
# 19 8            42
# 20 9            40



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

# 
# ## for paper
# # just cluster labels
# DimPlot(burtoni.snseq.combined.sct.neurons.recluster,
#         group.by='integrated_snn_res.0.8',
#         alpha = 0,
#         label=TRUE) +
#   ggtitle('Neurons') +
#   NoLegend()+
#   theme(
#     plot.title = element_blank()
#   )
# ggsave('neurons/Clusters.labels.dimplot.neurons.all.paper.pdf',
#        height = 7,
#        width = 7,
#        units = 'in',
#        dpi = 720)
# 
# # just cell type label
# DimPlot(burtoni.snseq.combined.sct.neurons.recluster,
#         group.by='sctypemarkers.hypo',
#         alpha = 0,
#         label=T,
#         label.box = T,
#         cols = c('#00BA38',
#                  '#B79F00')) +
#   ggtitle('Neurons') +
#   NoLegend()+
#   theme(
#     plot.title = element_blank()
#   )
# ggsave('neurons/Celltype.labels.dimplot.neurons.all.paper.pdf',
#        height = 7,
#        width = 7,
#        units = 'in',
#        dpi = 720)
# 
# # just color cells
# DimPlot(burtoni.snseq.combined.sct.neurons.recluster,
#         group.by='sctypemarkers.hypo',
#         label=F,
#         cols = c('#00BA38',
#                  '#B79F00')) +
#   ggtitle('Neurons') +
#   NoLegend()+
#   theme(
#     plot.title = element_blank()
#   ) +
#   xlab('Dim 1') +
#   ylab('Dim 2')
#  ggsave('neurons/Celltype.dimplot.neurons.all.paper.pdf',
#         height = 7,
#         width = 7,
#         units = 'in',
#         dpi = 720)

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

# get cell id and colors
mapping = burtoni.snseq.combined.sct.neurons.recluster@meta.data %>%
  dplyr::select(seurat_clusters,
                sctypemarkers.hypo) %>%
  distinct() %>%
  dplyr::rename(label = seurat_clusters, 
         sctypemarkers.hypo = sctypemarkers.hypo) %>%
  mutate(label = as.character(label))

# paper
ggtree(tree) %<+% mapping +
  scale_y_reverse() +
  geom_tree() +
  theme_tree() +
  geom_tiplab(offset = 1.5,
              angle = 90,
              vjust = 0.5,
              hjust =0.5,
              size = 3) +
  geom_tippoint(aes(color = sctypemarkers.hypo),
                shape = 16,
                size = 5) +
  coord_cartesian(clip = 'off') +
  theme(plot.margin = unit(c(0,0.5,0,0),
                           'in')) +
  scale_color_manual(values = c(
                                'C7-2: GABA' = '#B79F00',
                                'C7-1: GLU' = '#00BA38')) +
  theme(legend.position = 'none') 
ggsave('sctype.hypo/unknown/Neurons cluster cell type tree.pdf',
       height = 6.5,
       width = 3,
       units = 'in',
       dpi = 720,
       limitsize = FALSE)

# remove extra 
rm(tree)
rm(burtoni.snseq.combined.sct.neurons.recluster.tree)

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
  dplyr::select(c('seurat_clusters', 
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
  dplyr::select(-c('total_cell_count')) %>%
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
  
  rm(burtoni.snseq.combined.sct.neurons.recluster.tree)
  
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
         height = 3,
         width = 6.5,
         units = 'in')
  

  # paper
  # graph clusters per status
  neuron.cluster.genotype.count %>% 
    ggplot() +
    geom_point(aes(x = seurat_clusters,
                   y = Percent,
                   group = Status,
                   color = Status),
               position = position_dodge(width = 0.5),
               size = 5) +
    geom_text(data = neuron.cluster.glm %>% 
                filter(contrast != "(Intercept)") %>% 
                dplyr::select(contrast,
                              seurat_clusters,
                              adj.pvalue.round) %>% 
                mutate(sig = case_when(adj.pvalue.round < 0.01  ~ '*',
                                       TRUE ~ NA)) %>% 
                na.omit(),
              aes(x = seurat_clusters,
                  y = max(neuron.cluster.genotype.count$Percent)+2,
                  label = sig),
              size = 8) +
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
    xlab('neuron clusters') +
    theme(legend.position = 'inside',
          legend.position.inside = c(0.9,
                                     0.75),
          text = element_text(size = 16))
  ggsave('neurons/neuron_cluster/Neuron clusters binomial GLM cell percent.pdf',
         height = 6,
         width = 13,
         units = 'in',
         dpi = 720)
  
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
             size = 5) +
    geom_vline(xintercept = 0) +
    geom_vline(xintercept = 0.5,
               linetype = 2) +
    geom_vline(xintercept = -0.5,
               linetype = 2) +
    # geom_text(aes(label = seurat_clusters),
    #           color = 'white',
    #           size = 3) +
    labs(x = "Effect size", 
         y = "Neuron clusters") +
    theme_classic() +
    scale_color_manual(breaks = c('dom_bias',
                                  'sub_bias',
                                  'none'),
                       values = c(
                         "#4e499e",
                         "#60bb46",
                         "grey")) +
    xlim(-max.diff.value-0.25,
         max.diff.value+0.25) +
    # ggtitle("All clusters: binomial GLM cell count") +
    theme(legend.position = 'none')+
    theme(axis.text = element_text(size = 16),
          axis.title = element_text(size = 16))+
    scale_x_reverse()
  ggsave('neurons/neuron_cluster/Neuron clusters binomial GLM forest plot.pdf',
         height = 6,
         width = 6,
         units = 'in',
         dpi = 720)
  
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
  
  ## save glm table
  write_csv(neuron.cluster.glm,
            file = 'neurons/neuron_cluster/neuron_seurat_cluster_glm.csv')
  
  
#### IEG analysis Neuron ####
### compare proportions of IEG expressing neurons
IEG.list = c('egr1',
               'npas4',
               'fosb',
               'fosab')

## get IEG count data
IEG.neuron.counts = GetAssayData(burtoni.snseq.combined.sct.neurons.recluster, 
                                 assay = 'RNA',
                                 slot = "counts") 
  
# subset to IEG list
IEG.neuron.counts = IEG.neuron.counts[IEG.list,]

# get cell ID and clusters
IEG.neuron.counts.meta = burtoni.snseq.combined.sct.neurons.recluster@meta.data %>%
    rownames_to_column("cell") %>%
    dplyr::select(
      cell,
      seurat_clusters,
      orig.ident,
      replicate = Genotype.id
    )

# create 1 IEG combined score
IEG.neuron.total = colSums(IEG.neuron.counts)

# combine 
IEG.neuron.counts = rbind(IEG.neuron.counts,
                          IEG.neuron.total)

# if count is above 1 indicate as IEG expressed
IEG.neuron.detected = IEG.neuron.counts > 1
  
# ccombine with meta data
IEG.neuron = as.data.frame(as.matrix(IEG.neuron.detected)) %>%
    rownames_to_column("gene") %>%
    pivot_longer(
      cols = -gene,
      names_to = "cell",
      values_to = "expressed"
    ) %>%
    left_join(IEG.neuron.counts.meta, 
              by = "cell")
  
# filter out 
# needs to be in at least 5% of cells in a cluster
IEG.neuron = IEG.neuron %>%
    group_by(
      gene,
      seurat_clusters,
      orig.ident,
      replicate
    ) %>%
    summarize(
      pct = mean(expressed) * 100,
      .groups = "drop"
    ) %>%
    group_by(gene, 
             seurat_clusters) %>%
    filter(
      max(pct) >= 5 &          
        n_distinct(pct) > 1      
    ) %>%
    ungroup()
  
## run FDR corrected t.test of propotion expressing 
# make sure that it is not all 0 and have all genotypes
# significant if: 95% CI doesn't cross 0, OR that FDR < 0.15
IEG.neuron.prop  = IEG.neuron %>%
    group_by(gene, 
             seurat_clusters) %>%
    filter(n() == 6) %>%                   
    filter(sd(pct) > 0) %>%                 
    do(tidy(t.test(pct ~ orig.ident, 
                   data = .))) %>%
    ungroup() %>%
    group_by(seurat_clusters) %>% 
    mutate(fdr = p.adjust(p.value, method = "fdr")) %>% 
    ungroup()  %>% 
    mutate(
           ci_excludes_zero = conf.low > 0 | conf.high < 0,
           direction = case_when(
             conf.low > 0  ~ "Dom > Sub",
             conf.high < 0 ~ "Dom < Sub",
             fdr < 0.15 & estimate > 0  ~ "Dom > Sub",
             fdr < 0.15 & estimate < 0 ~ "Dom < Sub",
             TRUE ~ "not.sig"))
  
## graph IEG proportions across neuron clusters
IEG.neuron.prop %>% 
  filter(gene == 'IEG.neuron.total') %>% 
  ggplot(aes(x = estimate1,
             y = estimate2)) +
  geom_abline(slope = 1,
              intercept = 0) +
  geom_point(aes(color = direction),
             size = 5) +
  theme_bw() +
  xlab('DOM proportion IEG+ neurons')+
  ylab('SUB proportion IEG+ neurons') +
  coord_fixed(xlim = c(0,60),
              ylim = c(0,60)) +
  geom_text_repel(aes(label = seurat_clusters),
                  box.padding = 0.5, 
                  max.overlaps = Inf)+
  scale_color_manual(values = c("Dom < Sub" = '#60bb46',
                                "Dom > Sub" = '#4e499e',
                                "not.sig" = 'black')) +
  theme(legend.position = 'none')
ggsave('neurons/IEG proportion social status.png')  

# paper
IEG.neuron.prop %>% 
  filter(gene == 'IEG.neuron.total') %>% 
  mutate(label = ifelse(estimate > 0,
                        as.character(seurat_clusters),
                        NA),
         label2 = ifelse(estimate < 0,
                         as.character(seurat_clusters),
                        NA)) %>% 
  ggplot(aes(x = estimate1,
             y = estimate2)) +
  geom_abline(slope = 1,
              intercept = 0) +
  geom_point(aes(color = direction),
             size = 3) +
  theme_bw() +
  xlab('DOM IEG+ neurons (%)')+
  ylab('SUB IEG+ neurons (%)') +
  coord_fixed(xlim = c(0,60),
              ylim = c(0,60)) +
  geom_text_repel(aes(label = label),
                  box.padding = 1,
                  min.segment.length = 0, 
                  point.padding = 0.5,
                  max.overlaps = Inf,
                  # hjust = 'outward',
                  nudge_x = 1)+
  geom_text_repel(aes(label = label2),
                  box.padding = 1,
                  min.segment.length = 0, 
                  point.padding = 0.5,
                  max.overlaps = Inf,
                  # hjust = 'outward',
                  nudge_x = -1)+
  scale_color_manual(values = c("Dom < Sub" = '#60bb46',
                                "Dom > Sub" = '#4e499e',
                                "not.sig" = 'black')) +
  theme(legend.position = 'none') 
ggsave('neurons/IEG proportion social status paper.pdf',
       height = 3.3,
       width =3.3,
       units = 'in',
       dpi = 720)  

  
  
#### DEG analysis Neuron ####
# create dummy dataset
neuron.obj = burtoni.snseq.combined.sct.neurons.recluster

# use counts assay for limmatrend
# needed for pseudobulking to renormalize
DefaultAssay(neuron.obj) = "RNA"

### pre-filter genes
## remove genes that are not in 10% of cells within any cluster and treatment
neuron.obj$cluster_orig.ident = paste0(neuron.obj$seurat_clusters, 
                                       "_", 
                                       neuron.obj$orig.ident)

# get count data
neuron.obj.counts = GetAssayData(neuron.obj, 
                                 assay = "RNA", 
                                 slot = "counts")

# get list of cluster by treatment
neuron.obj.groups = neuron.obj$cluster_orig.ident

# get percent of cells with > 0 reads per group
neuron.obj.pct = sapply(unique(neuron.obj.groups), function(g) {
  cells = which(neuron.obj.groups == g)
  rowMeans(neuron.obj.counts[, cells, drop = FALSE] > 0) * 100
  })

# get max percent value per gene
neuron.obj.pct.max = apply(neuron.obj.pct, 1, max)

# keep genes that have a max of at least 10% 
neuron.obj.pct.keep = names(neuron.obj.pct.max[neuron.obj.pct.max >= 10])
# down to 6019 genes

# remove mito genes
neuron.obj.pct.keep = neuron.obj.pct.keep[!(neuron.obj.pct.keep %in% mito.genes.list)]
# down to 6006 genes

# subset seurat object
neuron.obj = subset(neuron.obj, 
                    features = neuron.obj.pct.keep)
# get an error but appears to have worked
# Features(neuron.obj) %>%
#   length()
# # [1] 6006

### use run_de to get DEG across clusters using pseudobulk
## get cell information
# set cell ID info that run_de will recognize
neuron.obj$cell_type = neuron.obj$seurat_clusters
neuron.obj$replicate = neuron.obj$Genotype.id
# make sure order is sub vs dom 
neuron.obj$label = factor(neuron.obj$orig.ident, levels = c("sub_burtoni_snseq", 
                                                             "dom_burtoni_snseq"))

## run limma trend pseudobulk
neuron.deg.results = run_de(
  neuron.obj,
  de_family = "pseudobulk",
  de_method = "limma",
  de_type   = "trend",
  min_reps = 3)

# # remove large dummy seurat object
# rm(neuron.obj)

## filter 
# needs to be in at least 10% of cluster in at least one group
neuron.deg.results = neuron.deg.results %>%
  filter(sub_burtoni_snseq.pct >= 0.1 | dom_burtoni_snseq.pct >= 0.1) %>%
  group_by(cell_type) %>%
  mutate(p_val_adj = p.adjust(p_val, method = "fdr")) %>%
  ungroup()

# save results
write.csv(neuron.deg.results,
          'neurons/neuron.deg.results.csv',
          row.names = F)


### compare DEG to proportion shift
## get data
# get counts data
neruon.deg.counts = GetAssayData(neuron.obj, 
                      assay = "RNA",
                      slot = "counts")
# get metadata
neuron.deg.meta   <- neuron.obj@meta.data

# create genotype by cluster ID
neruon.deg.group.id = paste(neuron.deg.meta$seurat_clusters,
                  neuron.deg.meta$Genotype.id, 
                  sep = "_") 
  

# get percentage of cells expressing gene 
# threshold RNA count of > 2
neruon.deg.pct.list <- lapply(unique(neruon.deg.group.id), function(g) {
  cells <- which(neruon.deg.group.id == g)
  
  data.frame(
    gene     = rownames(neruon.deg.counts),
    group    = g,
    pct_expr = rowMeans(neruon.deg.counts[, cells, drop = FALSE] > 2)
  )
})

# convert to data frame
neruon.deg.pct.list.df = do.call(rbind, 
                                 neruon.deg.pct.list)

# add meta data
neruon.deg.pct.list.df = neruon.deg.pct.list.df %>%
  separate(group,
           into = c("seurat_clusters", 
                    "Genotype.id"), 
           sep = "_") %>%
  left_join(
    neuron.deg.meta %>% 
      dplyr::select(seurat_clusters,
                    Genotype.id, 
                    orig.ident) %>%
      distinct(),
    by = c("seurat_clusters", 
           "Genotype.id"))


# get summary 
neruon.deg.pct.list.df.summary = neruon.deg.pct.list.df %>%
  group_by(gene,
           seurat_clusters,
           orig.ident) %>%
  summarise(mean_pct = mean(pct_expr), 
            .groups = "drop")

# get difference between dom and sub percentage
neruon.deg.pct.list.df.summary.wide = neruon.deg.pct.list.df.summary %>%
  tidyr::pivot_wider(names_from = orig.ident, 
                     values_from = mean_pct) %>%
  mutate(delta_pct = dom_burtoni_snseq - sub_burtoni_snseq)


# add to DEG results
neuron.deg.results = neuron.deg.results %>%
  left_join(neruon.deg.pct.list.df.summary.wide %>% 
              mutate(seurat_clusters = as.numeric(seurat_clusters)),
            by = c("gene", 
                   "cell_type" = "seurat_clusters"))

# look for concordance
# use 10% percentage difference as threshold
neuron.deg.results = neuron.deg.results %>%
  mutate(concordance = ifelse(sign(avg_logFC) == sign(delta_pct),
                              'concordant',
                              'discordant'),
         Likely_DEG.type = case_when(
           abs(delta_pct) > 0.2 & concordance == 'concordant' ~ "Fraction-driven",
           abs(delta_pct) <= 0.2 & concordance == 'concordant' ~ "Intensity-driven",
           abs(delta_pct) <= 0.2 & concordance == 'discordant' ~ "Intensity-driven",
           TRUE ~ "Mixed"
         ),
         DEG.type = ifelse(avg_logFC > 0,
                           'Dom',
                           'Sub'))

# save results
write.csv(neuron.deg.results,
          'neurons/neuron.deg.results.csv',
          row.names = F)


### graph
## compare logFC to DEG
# scatterplot
neuron.deg.results %>% 
  filter(p_val_adj < 0.05) %>% 
ggplot(aes(x = delta_pct, 
           y = avg_logFC, 
           shape = Likely_DEG.type,
           color = DEG.type)) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0.2,
             linetype = "dashed") +
  geom_vline(xintercept = -0.2,
             linetype = "dashed") +
  geom_vline(xintercept = 0) +
  geom_point(size = 0.8) +
  facet_wrap(~ cell_type, scales = "free") +
  theme_classic() +
  scale_color_manual(values = c('Sub' = '#60bb46',
                                'Dom' = '#4e499e'))+
  labs(
    x = "Percent expression difference (Dom - Sub)",
    y = "Log2 Fold Change (Pseudobulk)",
    shape = "Likely DEG type",
    color = 'Social status'
  )
ggsave('neurons/DEG logfc vs percentage neruons.png')

# barplot
neuron.deg.results %>% 
  filter(p_val_adj < 0.05) %>% 
  group_by(cell_type) %>% 
  dplyr::count(Likely_DEG.type,
        DEG.type) %>% 
  mutate(n = ifelse(DEG.type == 'Sub',
                    -n,
                    n)) %>%
  ggplot(aes(x = n, 
             y = as.factor(cell_type), 
             fill = DEG.type,
             color = Likely_DEG.type,
             alpha = Likely_DEG.type)) +
  geom_bar(stat = 'identity') +
  theme_classic() +
  scale_fill_manual(values = c('Sub' = '#60bb46',
                                'Dom' = '#4e499e')) +
  scale_color_manual(values = c('Fraction-driven' = 'black',
                               'Intensity-driven' = 'black')) +
  scale_alpha_manual(values = c('Fraction-driven' = 0.6,
                                'Intensity-driven' = 1)) +
  labs(
    y = "Cluster",
    x = "DEG count",
    color = "Likely DEG type",
    fill = 'Social status'
  ) +
  theme(legend.position = 'inside',
        legend.position.inside = c(0.85,0.85))
ggsave('neurons/DEG counts per cluster neuron.png')

# paper
neuron.deg.results %>% 
  filter(p_val_adj < 0.05) %>% 
  group_by(cell_type) %>% 
  dplyr::count(Likely_DEG.type,
               DEG.type) %>% 
  mutate(n = ifelse(DEG.type == 'Sub',
                    -n,
                    n)) %>%
  ggplot(aes(x = n, 
             y = as.factor(cell_type), 
             fill = DEG.type,
             # color = Likely_DEG.type,
             alpha = Likely_DEG.type)) +
  geom_bar(stat = 'identity') +
  theme_classic() +
  scale_fill_manual(values = c('Sub' = '#60bb46',
                               'Dom' = '#4e499e')) +
  # scale_color_manual(values = c('Fraction-driven' = 'black',
  #                               'Intensity-driven' = 'black')) +
  scale_alpha_manual(values = c('Fraction-driven' = 0.6,
                                'Intensity-driven' = 1)) +
  labs(
    y = "neuron cluster",
    x = "DEG count",
  ) +
  theme(legend.position = 'none',
        text = element_text(size = 12))
ggsave('neurons/DEG counts per cluster neuron.pdf',
       height = 3.25,
       width = 3.25,
       units = 'in',
       dpi = 720)


## graph upset
# get list of DEGs
neuron.significant_degs = neuron.deg.results %>%
  filter(p_val_adj < 0.05) %>%
  distinct(cell_type,
           gene) %>% 
  dplyr::select(cell_type,
                gene)

# create name list
neuron.deg_list = split(neuron.significant_degs$gene, 
                        neuron.significant_degs$cell_type)

# just graph overlap
ComplexUpset::upset(fromList(neuron.deg_list), 
                    intersect = names(neuron.deg_list),
                    min_size = 2) 
ggsave('neurons/Upset neuron.png')


# overlap
# get list of DEGs
neuron.significant_degs = neuron.deg.results %>%
  filter(p_val_adj < 0.05) %>%
  distinct(cell_type,
           gene) %>% 
  dplyr::select(cell_type,
                gene)

# create name list
neuron.deg_list = split(neuron.significant_degs$gene, 
                        neuron.significant_degs$cell_type)

# just graph overlap
ComplexUpset::upset(fromList(neuron.deg_list), 
                    intersect = names(neuron.deg_list),
                    sort_intersections_by = c('degree', 
                                              'cardinality'),
                    # min_size = 2,
                    min_degree = 10) 
ggsave('neurons/Upset overlap neuron.png')

# paper
ComplexUpset::upset(fromList(neuron.deg_list), 
                    intersect = names(neuron.deg_list),
                    sort_intersections_by = c('degree', 
                                              'cardinality'),
                    # min_size = 2,
                    min_degree = 10) 
ggsave('neurons/Upset overlap neuron.pdf',
       width = 6.5,
       height = 3,
       units = 'in',
       dpi = 720)


# up vs down
# get list of DEGs
neuron.significant_degs = neuron.deg.results %>%
  filter(p_val_adj < 0.05) %>%
  distinct(cell_type,
           gene,
           DEG.type) %>% 
  mutate(cell_DEG_type = paste0(cell_type,
                                '_',
                                DEG.type)) %>% 
  dplyr::select(cell_DEG_type,
                gene,) 

# create name list
neuron.deg_list = split(neuron.significant_degs$gene, 
                        neuron.significant_degs$cell_DEG_type)

# just graph overlap
ComplexUpset::upset(fromList(neuron.deg_list), 
                    intersect = names(neuron.deg_list),
                    min_size = 3) 
ggsave('neurons/Upset neuron status.png')

# check genes in all clusters (or most)
neuron.deg.results %>%
     filter(p_val_adj < 0.05) %>%
     dplyr::count(gene) %>%
     filter(n > 15)

neuron.deg.results %>% 
  filter(gene %in% c(neuron.deg.results %>%
  filter(p_val_adj < 0.05) %>% 
  dplyr::count(gene) %>% 
  filter(n > 16) %>% 
    pull(gene))) %>% 
    distinct(Likely_DEG.type,
             DEG.type,
             gene) 

# Likely_DEG.type DEG.type               gene
# 1 Intensity-driven      Dom ENSONIG00000011023
# 2 Intensity-driven      Dom                cga
# 3  Fraction-driven      Dom                gh1
# 4 Intensity-driven      Dom                lhb
# 5  Fraction-driven      Dom                oxt
# 6  Fraction-driven      Dom                prl


# # check DEG numbers
# neuron.deg.results %>%
#   filter(p_val_adj < 0.15) %>%
#   dplyr::count(cell_type)
# # A tibble: 20 × 2
# cell_type     n
# <chr>     <int>
#   1 0           501
# 2 1           385
# 3 10           20
# 4 11           16
# 5 12           53
# 6 13           32
# 7 14           20
# 8 15           20
# 9 16           20
# 10 17           10
# 11 18            8
# 12 19           17
# 13 2           175
# 14 3           142
# 15 4            94
# 16 5            81
# 17 6            57
# 18 7            74
# 19 8            42
# 20 9            40



#### Get hypomap clusters data ####
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
  geom_text(color = 'white') +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1)) +
  scale_fill_gradient(low = 'grey',
                      high = 'grey1')
ggsave('neurons/Hypomap/Hypomap C66 clusters overlap brain regions.png',
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
  rbind(tilapia.C66.gene %>% 
          filter(oniloticus_homolog_associated_gene_name != "")%>%
          dplyr::rename(gene = oniloticus_homolog_associated_gene_name) %>% 
          mutate(gene = tolower(gene))) %>% 
  rbind(tilapia.C66.gene %>% 
          filter(oniloticus_homolog_associated_gene_name != "")%>%
          dplyr::rename(gene = oniloticus_homolog_associated_gene_name) %>% 
          mutate(gene = toupper(gene))) %>% 
  rbind(tilapia.C66.ID %>% 
          filter(oniloticus_homolog_ensembl_gene != "")%>%
          dplyr::rename(gene = oniloticus_homolog_ensembl_gene) ) %>% 
  distinct() %>% 
  dplyr::rename(Mouse_gene = external_gene_name)

# combine with specificity data and region
tilapia.C66.gene.list = tilapia.C66.gene.list %>% 
  left_join(C66.data) %>% 
  distinct() %>% 
  left_join(hypomap.region.data.reduce) %>% 
  distinct()

# 3845 genes
nrow(tilapia.C66.gene.list)

# remove NA
tilapia.C66.gene.list = tilapia.C66.gene.list %>% 
  filter(!(is.na(gene)))

# 3845 genes
nrow(tilapia.C66.gene.list)

# remove genes not in dataset
tilapia.C66.gene.list = tilapia.C66.gene.list %>% 
  filter(gene %in% rownames(burtoni.snseq.combined.sct.neurons.recluster@assays$SCT$counts))

# 1684 genes
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
  geom_point(data = C66.data %>% 
               dplyr::select(cluster_id) %>% 
               table() %>% 
               data.frame(),
             fill = 'black') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1)) +
  xlab('Hypomap Cluster_ID') +
  ylab('Tilapia gene count') +
  geom_hline(yintercept = 10) 
ggsave('neurons/Hypomap/Hypomap C66 clusters gene count.png',
       height = 5,
       width = 10)


# save data
write.csv(tilapia.C66.gene.list,
          'neurons/Hypomap/tilapia.C66.gene.list.csv',
          row.names = FALSE)

#### Neuron spatial prediction ####
## load hypomap marker
tilapia.C66.gene.list = read.csv('neurons/Hypomap/tilapia.C66.gene.list.csv')

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
# create region proportions 
hypomap.region.data.prop = hypomap.region.data %>% 
  dplyr::select(cluster_id_66,
                Region_summarized,
                ncells) %>% 
  na.omit() %>%
  group_by(cluster_id_66) %>% 
  mutate(ncells.prop = ncells/sum(ncells)) %>% 
  ungroup() %>% 
  dplyr::rename(cluster_id = cluster_id_66) %>% 
  dplyr::rename(Region = Region_summarized) 


### find markers for all neuron clusters
# prepare data
burtoni.snseq.combined.sct.neurons.recluster = PrepSCTFindMarkers(burtoni.snseq.combined.sct.neurons.recluster)

## Neuron markers
neuron.markers.df = FindAllMarkers(burtoni.snseq.combined.sct.neurons.recluster, 
                                       min.pct = 0.1, 
                                       assay = 'SCT',
                                   only.pos = T)

# create specficity score
neuron.markers.df = neuron.markers.df %>% 
  mutate(specificity = avg_log2FC*(pct.1/pct.2))

# 959 unique genes from hypomap
tilapia.C66.gene.list$gene %>% 
  unique() %>% 
  length()
# 4549 unique genes from findmarkers
neuron.markers.df$gene %>% 
  unique() %>% 
  length()
# 418 genes overlapping
intersect(unique(neuron.markers.df$gene),
      unique(tilapia.C66.gene.list$gene)) %>% 
  length()


# filter out clusters with less than 10 markers 
tilapia.C66.gene.list = tilapia.C66.gene.list %>% 
  filter(!(cluster_id %in% c(tilapia.C66.gene.list %>% 
                               filter(gene %in% unique(neuron.markers.df$gene)) %>% 
                               dplyr::select(cluster_id) %>% 
                               table() %>% 
                               data.frame() %>% 
                               filter(Freq < 10) %>% 
                               pull(cluster_id) %>% 
                               droplevels())))


# graph number of genes per cluster
tilapia.C66.gene.list %>% 
  dplyr::select(cluster_id) %>% 
  table() %>% 
  data.frame() %>% 
  ggplot(aes(x = reorder(cluster_id, 
                         -Freq),
             y = Freq)) +
  geom_bar(stat = 'identity') +
  geom_bar(data = tilapia.C66.gene.list %>% 
             filter(gene %in% unique(neuron.markers.df$gene)) %>% 
             dplyr::select(cluster_id) %>% 
             table() %>% 
             data.frame(),
           stat = 'identity',
           fill = 'black') +
  # geom_point(data = C66.data %>% 
  #              dplyr::select(cluster_id) %>% 
  #              filter(cluster_id %in% unique(tilapia.C66.gene.list$cluster_id)) %>% 
  #              table() %>% 
  #              data.frame(),
  #            fill = 'black',
  #            size = 3) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1)) +
  xlab('') +
  ylab('Tilapia gene count') +
  geom_hline(yintercept = 10,
             color = 'red',
             linetype = 'dashed') +
  ggtitle('Overlap of Hypomap markers and tilapia genes') +
  theme(text = element_text(size = 12))
ggsave('neurons/Hypomap/Hypomap C66 clusters gene count findmarkers filter.pdf',
       height = 3.25,
       width = 6.5,
       units = 'in',
       dpi = 320)


# graph overlap of clusters and brain regions
hypomap.region.data %>% 
  dplyr::select(cluster_id_66,
                Region_summarized,
                ncells) %>% 
  filter(cluster_id_66 %in% unique(tilapia.C66.gene.list$cluster_id)) %>% 
  na.omit() %>%
  group_by(cluster_id_66,
           Region_summarized) %>% 
  summarise(ncells = sum(ncells)) %>% 
mutate(Region.zone = case_when(Region_summarized %in% c('Medial preoptic area',
                                             'Lateral preoptic area',
                                             '(Anterior/Preoptic)Periventricular region') ~ 0,
                               Region_summarized %in% c('Anterior hypothalamic nucleus',
                                             'Paraventricular hypothalamic nucleus',
                                             'Suprachiasmatic nucleus') ~ 1,
                               Region_summarized %in% c('Arcuate hypothalamic nucleus',
                                             'Ventromedial hypothalamic nucleus') ~ 2,
                               Region_summarized %in% c('Dorsomedial nucleus of the hypothalamus') ~ 3,
                               Region_summarized %in% c('Lateral hypothalamic area') ~ 4,
                               Region_summarized %in% c('(Pre)Mammillary region',
                                             'Periventricular hypothalamic nucleus, posterior part') ~ 5,
                               TRUE ~ 6)) %>% 
ggplot(aes(x = cluster_id_66,
             y = fct_reorder(Region_summarized,
                         Region.zone),
             fill = log(ncells),
             label = ncells)) +
  geom_tile() +
  geom_segment(x=-6, xend = 0.4, y = 3.5, yend =3.5, color = 'black', linewidth = 0.6)+
  geom_segment(x=-6, xend = 0.4, y = 6.5, yend =6.5, color = 'black', linewidth = 0.6)+
  geom_segment(x=-6, xend = 0.4, y = 8.5, yend =8.5, color = 'black', linewidth = 0.6)+
  geom_segment(x=-6, xend = 0.4, y = 9.5, yend =9.5, color = 'black', linewidth = 0.6)+
  geom_segment(x=-6, xend = 0.4, y = 10.5, yend = 10.5, color = 'black', linewidth = 0.6)+
  geom_segment(x=-6, xend = 0.4, y = 12.5, yend = 12.5, color = 'black', linewidth = 0.6)+
  coord_cartesian(clip = 'off')+
  # geom_text(color = 'white',
  #           size = 2) +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1,
                                   size = 8),
        axis.text.y = element_text(size = 6),
        legend.position = 'none',
        plot.margin = margin(l = 0, t =0.5)) +
  scale_fill_gradient(low = 'grey',
                      high = 'grey1') +
  ylab('') +
  xlab('')
ggsave('neurons/Hypomap/Hypomap C66 clusters overlap brain regions.pdf',
       height = 3.5,
       width = 6.5,
       dpi = 720,
       units = 'in')


### create cluster specificity matrix
## neuron 
neuron.mat = neuron.markers.df %>% 
  filter(gene %in% unique(tilapia.C66.gene.list$gene)) %>% 
  dplyr::select(gene,
                cluster,
                specificity) %>% 
  pivot_wider(names_from = cluster,
              values_from = specificity,
              values_fill = 0) %>% 
  column_to_rownames('gene') %>% 
  as.matrix()

# log1p transform specificity 
neuron.mat = sign(neuron.mat) *log1p(abs(neuron.mat))

## hypomap 
tilapia.C66.mat = tilapia.C66.gene.list %>% 
  filter(gene %in% unique(neuron.markers.df$gene)) %>% 
  dplyr::select(gene,
                cluster_id,
                specificity) %>%
  group_by(cluster_id,
           gene) %>% 
  summarise(specificity = mean(specificity,
                               na.rm = TRUE), .groups = "drop") %>% 
  pivot_wider(names_from = cluster_id,
              values_from = specificity,
              values_fill = 0) %>% 
  column_to_rownames('gene') %>% 
  as.matrix()

# log1p transform specificity 
tilapia.C66.mat = sign(tilapia.C66.mat) *log1p(abs(tilapia.C66.mat))

### correlate across studies
## spearman correlation
neuron.corr.hypomap.mat =  cor(neuron.mat, 
               tilapia.C66.mat,
               method = "spearman")

# convert to z-score
neuron.corr.hypomap.mat.zscore = t(apply(neuron.corr.hypomap.mat, 1, function(x) {
  (x - mean(x)) / sd(x)
  }))


## permutation based z-score
# set number permutations
nperm = 1000

# create empty matrix
neuron.corr.hypomap.mat.zscore.null = matrix(NA, 
                                             nrow = ncol(neuron.mat), 
                                             ncol = ncol(tilapia.C66.mat))

rownames(neuron.corr.hypomap.mat.zscore.null) = colnames(neuron.mat)
colnames(neuron.corr.hypomap.mat.zscore.null) = colnames(tilapia.C66.mat)

# loop through each neuron cluster
for (i in seq_len(ncol(neuron.mat))) {
  
  neuron.mat_vec = neuron.mat[, i]
  
  # null distribution by permutation genes 
  # get correlation
  null_mat = replicate(nperm, {
    perm_vec = sample(neuron.mat_vec)
    apply(tilapia.C66.mat, 2, function(mv)
      cor(perm_vec, mv, method = "spearman"))
  })
  
  # get distributions 
  null_mean <- rowMeans(null_mat)
  null_sd   <- apply(null_mat, 1, sd)
  
  # get observation
  obs <- neuron.corr.hypomap.mat[i, ]
  
  # create z-score
  neuron.corr.hypomap.mat.zscore.null[i, ] <- (obs - null_mean) / null_sd
}


# best match by z score
neuron.corr.hypomap.mat.zscore.null.best = apply(neuron.corr.hypomap.mat.zscore.null, 1, function(x) {
  colnames(neuron.corr.hypomap.mat.zscore.null)[which.max(x)]
  }) %>% 
  as.data.frame() %>% 
  rownames_to_column('cluster')
colnames(neuron.corr.hypomap.mat.zscore.null.best)[2] = "cluster_id"

# compare top 2 best matches to see if one is much better
neuron.corr.hypomap.mat.zscore.null.gap = apply(neuron.corr.hypomap.mat.zscore.null, 1, function(x) {
  s <- sort(x, decreasing = TRUE)
  100*(s[1] - s[2])/s[1]
  }) %>% 
  as.data.frame()%>% 
  rownames_to_column('cluster')
colnames(neuron.corr.hypomap.mat.zscore.null.gap)[2] = "comparison.top.hits"

# combine
neuron.corr.hypomap.mat.zscore.null.best = neuron.corr.hypomap.mat.zscore.null.best %>% 
  left_join(neuron.corr.hypomap.mat.zscore.null.gap)

# filter out hypomap clusters that don't have a score above 2
neuron.corr.hypomap.mat.zscore.null.filter = neuron.corr.hypomap.mat.zscore.null
# neuron.corr.hypomap.mat.zscore.null.filter = neuron.corr.hypomap.mat.zscore.null[, colSums(neuron.corr.hypomap.mat.zscore.null > 3) > 0, drop = FALSE]

# remove low zscore values
neuron.corr.hypomap.mat.zscore.null.filter[neuron.corr.hypomap.mat.zscore.null.filter < 1] <- 0


## create region score
# convert to data frame
neuron.corr.hypomap.mat.zscore.null.filter.df = neuron.corr.hypomap.mat.zscore.null.filter %>% 
  as.data.frame() %>% 
  rownames_to_column('cluster') %>% 
  pivot_longer(cols = -c(cluster),
               names_to = 'cluster_id',
               values_to = 'z.score.null')

# add brain region data 
neuron.corr.hypomap.mat.zscore.null.filter.df.region = neuron.corr.hypomap.mat.zscore.null.filter.df %>% 
  left_join(hypomap.region.data.prop) %>% 
  mutate(z.score.null.weight = z.score.null*ncells.prop) %>% 
  group_by(cluster,
           Region) %>% 
  summarise(z.score.null.weight = sum(z.score.null.weight)) %>% 
  droplevels() %>% 
  na.omit() 

## collapse brain regions 
neuron.corr.hypomap.mat.zscore.null.filter.df.zone = neuron.corr.hypomap.mat.zscore.null.filter.df %>% 
  left_join(hypomap.region.data.prop) %>% 
  mutate(Region.zone = case_when(Region %in% c('Medial preoptic area',
                                               'Lateral preoptic area',
                                               '(Anterior/Preoptic)Periventricular region') ~ 'Preoptic',
                                 Region %in% c('Anterior hypothalamic nucleus',
                                               'Paraventricular hypothalamic nucleus',
                                               'Suprachiasmatic nucleus') ~ 'Anterior',
                                 Region %in% c('Arcuate hypothalamic nucleus',
                                               'Ventromedial hypothalamic nucleus') ~ 'Ventral',
                                 Region %in% c('Dorsomedial nucleus of the hypothalamus') ~ 'Dorsal',
                                 Region %in% c('Lateral hypothalamic area') ~ 'Lateral',
                                 Region %in% c('(Pre)Mammillary region',
                                               'Periventricular hypothalamic nucleus, posterior part') ~ 'Posterior',
                                 TRUE ~ Region)) %>% 
  mutate(z.score.null.weight = z.score.null*ncells.prop) %>% 
  group_by(cluster,
           Region.zone) %>% 
  summarise(z.score.null.weight = sum(z.score.null.weight)) %>% 
  droplevels() %>% 
  na.omit() 

# create cluster percentage
neuron.corr.hypomap.mat.zscore.null.filter.df.zone = neuron.corr.hypomap.mat.zscore.null.filter.df.zone %>% 
  group_by(cluster) %>% 
  mutate(Percent.cluster = 100*z.score.null.weight/sum(z.score.null.weight)) %>% 
  ungroup() %>% 
  group_by(Region.zone) %>% 
  mutate(Percent.zone = 100*z.score.null.weight/sum(z.score.null.weight)) %>% 
  ungroup()

### graph results
# ## heatmap
# # correlation
# pheatmap(neuron.corr.hypomap.mat,
#          annotation_col = tilapia.C66.gene.list %>% 
#            dplyr::select(cluster_id,
#                          parent_id) %>% 
#            distinct() %>% 
#            filter(cluster_id %in% colnames(neuron.corr.hypomap.mat)) %>% 
#            column_to_rownames('cluster_id'))
# 
# # z-score
# pheatmap(neuron.corr.hypomap.mat.zscore,
#          annotation_col = tilapia.C66.gene.list %>% 
#            dplyr::select(cluster_id,
#                          parent_id) %>% 
#            distinct() %>% 
#            filter(cluster_id %in% colnames(neuron.corr.hypomap.mat.zscore)) %>% 
#            column_to_rownames('cluster_id'))

# # permutation z-score
# pheatmap(neuron.corr.hypomap.mat.zscore.null,
#          annotation_col = tilapia.C66.gene.list %>%
#            dplyr::select(cluster_id,
#                          parent_id) %>%
#            distinct() %>%
#            filter(cluster_id %in% colnames(neuron.corr.hypomap.mat.zscore.null)) %>%
#            column_to_rownames('cluster_id'))


# filtered z-score
pdf('neurons/Hypomap/Hypomap clusters vs neuron clusters.pdf',
    height = 4,
    width = 6.5)
pheatmap(neuron.corr.hypomap.mat.zscore.null.filter,
         annotation_col = tilapia.C66.gene.list %>%
           dplyr::select(cluster_id,
                         parent_id) %>%
           distinct() %>%
           filter(cluster_id %in% colnames(neuron.corr.hypomap.mat.zscore.null.filter)) %>%
           column_to_rownames('cluster_id'),
         color = colorRampPalette(c("white", "darkred"))(50),
         legend = FALSE, 
         annotation_legend = FALSE,
         fontsize = 10)
dev.off()


# Broad regions
neuron.corr.hypomap.mat.zscore.null.filter.df.zone %>% 
  filter(Region.zone != 'Zona incerta') %>% 
  filter(z.score.null.weight >= 0.5) %>% 
  ggplot(aes(y = Region.zone,
             x = as.factor(as.numeric(cluster)),
             fill = log(z.score.null.weight),
             label = round(z.score.null.weight,
                           1))) +
  geom_tile() +
  geom_text(size = 2) +
  theme_classic() +
  scale_fill_gradient(low = 'white',
                      high = 'darkred',
                      na.value = 'white') +
  ylab('')+
  xlab('Neuron cluster') +
  theme(legend.position = 'none',
        text = element_text(size = 10))
ggsave('neurons/Hypomap/Hypomap broad zones vs neuron clusters.pdf',
       height = 2.5,
       width = 6.5,
       units = 'in',
       dpi = 720)


# AVP clusters
neuron.corr.hypomap.mat.zscore.null.filter.df.zone %>% 
  filter(Region.zone != 'Zona incerta') %>% 
  filter(z.score.null.weight >= 0.5) %>% 
  filter(cluster %in% c(0,
                        1,
                        2,
                        3,
                        4,
                        11,
                        6,
                        5)) %>% 
  ggplot(aes(x = Region.zone,
             y = as.factor(as.numeric(cluster)),
             fill = log(z.score.null.weight),
             label = round(z.score.null.weight,
                           1))) +
  geom_tile() +
  geom_text() +
  theme_classic() +
  scale_fill_gradient(low = 'white',
                      high = 'darkred',
                      na.value = 'white') +
  xlab('')+
  ylab('Neuron cluster') +
  theme(legend.position = 'none')
ggsave('neurons/Hypomap/Hypomap broad zones vs AVP neuron clusters.png',
       height = 6,
       width = 3.25,
       units = 'in')

# filter to regions of interest
neuron.corr.hypomap.mat.zscore.null.filter.df.zone %>% 
  filter(Region.zone != 'Zona incerta') %>% 
  filter(z.score.null.weight >= 0.5) %>% 
  filter(cluster %in% c(0,
                        1,
                        2,
                        3,
                        4,
                        11,
                        6,
                        5)) %>% 
  filter(Region.zone %in% c('Preoptic',
                            'Anterior',
                            'Ventral',
                            'Posterior')) %>% 
  ggplot(aes(x = Region.zone,
             y = as.factor(as.numeric(cluster)),
             fill = log(Percent.zone),
             label = round(z.score.null.weight,
                           1))) +
  geom_tile() +
  geom_text() +
  theme_classic() +
  scale_fill_gradient(low = 'white',
                      high = 'darkred',
                      na.value = 'white') +
  xlab('')+
  ylab('Neuron cluster') +
  theme(legend.position = 'none') +
  theme(text = element_text(size = 20))
ggsave('neurons/Hypomap/Hypomap broad zones vs AVP neuron clusters filter.pdf',
       height = 6.5,
       width = 6.5,
       units = 'in',
       dpi = 720)



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
  layer = 'counts', # using normalized data
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

#use soft power threshold 10?


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
         ncol=4)
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
# paper
png('neurons/wgcna/UMAP MEs neurons.png',
    width = 6.5,
    height = 6.5,
    units = 'in',
    res = 300)
wrap_plots(plot_list.neurons.me, 
           ncol=4)
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


# plot with Seurat's DotPlot function
p <- DotPlot(burtoni.snseq.combined.sct.neurons.recluster, 
             features=mods.neurons, 
             group.by = 'integrated_snn_res.0.8')

# flip the x/y axes, rotate the axis labels, and change color scheme:
p <- p +
  RotatedAxis() +
  scale_color_gradient2(high='red', 
                        mid='grey95', 
                        low='blue') +
  ylab('Neuron cluster') +
  xlab('Neuron hdWGCNA modules')

# plot output
p
ggsave('neurons/wgcna/Dotplot modules across clusters.pdf',
       height = 7,
       width = 7,
       units = 'in',
       dpi = 720)

# radar plot
# ModuleRadarPlot(
#   burtoni.snseq.combined.sct.neurons.recluster,
#   group.by = 'integrated_snn_res.0.8',
#   base.size = 4,
#   axis.label.size=4,
#   grid.label.size=0,
#   legend.text.size = 5,
#   combine = T,
#   ncol = 4,
#   background.circle.colour = 'white',
#   gridline.mid.colour = 'grey25') +
#   plot_layout(ncol = 4) 

# paper
# need to wrap to adjust titles
p = ModuleRadarPlot(
  burtoni.snseq.combined.sct.neurons.recluster,
  group.by = 'integrated_snn_res.0.8',
  base.size = 4,
  axis.label.size=4,
  grid.label.size=0,
  legend.text.size = 5,
  combine = F,
  ncol = 4,
  background.circle.colour = 'white',
  gridline.mid.colour = 'grey25') 

p = lapply(p, function(x)
  x + theme(plot.title = element_text(size = 16)))

wrap_plots(p,
           ncol = 4)
ggsave('neurons/wgcna/Radar plot modules across clusters.pdf',
       height = 3.5,
       width = 7,
       units = 'in',
       dpi = 720)

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
                                                   features = mods.neurons,
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

# # graph DME
# PlotDMEsVolcano(
#   burtoni.snseq.combined.sct.neurons.recluster,
#   DMEs.neurons) +
#   theme_classic() +
#   xlim(-0.5, 0.5)
# ggsave('neurons/wgcna/ME by cluster volcanoplot.png',
#        width = 10,
#        height = 10)

# for poster
# create mod color list
mod_colors <- DMEs.neurons$module
names(mod_colors) <- as.character(DMEs.neurons$module)
# graph DME
# DMEs.neurons %>%
#   mutate(anno = ifelse(p_val_adj < 0.05,
#                        module,
#                        "")) %>%
#   ggplot(aes(x = avg_log2FC,
#              y = -log10(p_val_adj)))  +
#   geom_rect(aes(xmin = -Inf,
#                 xmax = Inf,
#                 ymin = -Inf,
#                 ymax = -log10(0.05)),
#             fill = "grey75",
#             alpha = 0.8,
#             color = NA) +
#   geom_vline(xintercept = 0,
#              linetype = "dashed",
#              color = "grey75",
#              alpha = 0.8) +
#   geom_point(aes(fill = module),
#              size = 8,
#              pch = 21,
#              color = "black") +
#   geom_text_repel(aes(label = anno),
#                   color = "black",
#                   min.segment.length = 0,
#                   max.overlaps = Inf,
#                   size = 8,
#                   point.padding = 4,
#                   nudge_x = 0.06) +
#   theme_classic() +
#   theme(legend.position = 'none') +
#   # xlim(c(-0.55,
#   #        0.55)) +
#   scale_fill_manual(values = mod_colors)+
#   xlab(bquote("Average log"[2] ~ "(Fold Change)")) +
#   ylab(bquote("-log"[10] ~ "(Adj. P-value)")) +
#   ggtitle('Differential module eigengene analysis') +
#   theme(panel.border = element_rect(color = "black",
#                                     fill = NA,
#                                     size = 1))+
#   theme(axis.text = element_text(size = 15))  +
#   theme(axis.title = element_text(size = 20))+
#   theme(plot.title = element_text(size=20))
# ggsave('neurons/wgcna/ME by cluster volcanoplot poster.pdf',
#        width = 10,
#        height = 5,
#        units = "in",
#        dpi = 320)
# 
# # for presentation
# DMEs.neurons %>%
#   mutate(anno = ifelse(p_val_adj < 0.05,
#                        module,
#                        "")) %>%
#   ggplot(aes(x = avg_log2FC,
#              y = -log10(p_val_adj)))  +
#   geom_rect(aes(xmin = -Inf,
#                 xmax = Inf,
#                 ymin = -Inf,
#                 ymax = -log10(0.01)),
#             fill = "grey75",
#             alpha = 0.8,
#             color = NA) +
#   geom_vline(xintercept = 0,
#              linetype = "dashed",
#              color = "grey75",
#              alpha = 0.8) +
#   geom_point(aes(fill = module),
#              size = 8,
#              pch = 21,
#              color = "black")  +
#   theme_classic() +
#   theme(legend.position = 'none') +
#   # xlim(c(-0.55,
#   #        0.55)) +
#   scale_fill_manual(values = mod_colors)+
#   xlab(bquote("Average log"[2] ~ "(Fold Change)")) +
#   ylab(bquote("-log"[10] ~ "(Adj. P-value)")) +
#   ggtitle('Differential module eigengene analysis') +
#   theme(panel.border = element_rect(color = "black",
#                                     fill = NA,
#                                     size = 1))+
#   theme(axis.text = element_text(size = 15))  +
#   theme(axis.title = element_text(size = 20))+
#   theme(plot.title = element_text(size=20))
# ggsave('neurons/wgcna/ME by cluster volcanoplot presentation.pdf',
#        width = 10,
#        height = 5,
#        units = "in",
#        dpi = 720)
# 

# for paper
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
  # xlim(c(-0.55,
  #        0.55)) +
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
ggsave('neurons/wgcna/ME by cluster volcanoplot paper.pdf',
       width = 10,
       height = 5,
       units = "in",
       dpi = 720)



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

#### DMEs WGCNA for neuron clusters ####
### identify clusters specific to cell types
## pseudo bulk MEs across clusters
ME.by.neuron.cluster = burtoni.snseq.combined.sct.neurons.recluster@meta.data %>% 
  group_by(seurat_clusters) %>% 
  summarize(
    across(all_of(mods.neurons), mean, na.rm = TRUE),
    .groups = "drop"
  ) 

# get standard deviation for each module
SD.ME.by.neuron.cluster = ME.by.neuron.cluster %>%
  summarize(
    across(all_of(mods.neurons), sd)
  ) %>%
  pivot_longer(
    everything(),
    names_to = "module",
    values_to = "sd_across_clusters"
  ) %>% 
  left_join(ME.by.neuron.cluster %>%
              summarize(
                across(all_of(mods.neurons), mean)
              ) %>%
              pivot_longer(
                everything(),
                names_to = "module",
                values_to = "mean_across_clusters"
              )) %>% 
  mutate(relative_variability = abs(sd_across_clusters/mean_across_clusters)) %>% 
  filter(relative_variability > 1) %>% 
  pull(module)

# get zscore
ME.by.neuron.cluster.zscore = ME.by.neuron.cluster %>%
  dplyr::select(seurat_clusters, 
         all_of(SD.ME.by.neuron.cluster)) %>%
  pivot_longer(
    -seurat_clusters,
    names_to = "module",
    values_to = "mean_ME"
  ) %>%
  group_by(module) %>%
  mutate(
    z_ME = (mean_ME - mean(mean_ME)) / sd(mean_ME)
  ) %>%
  ungroup()
  
# select modules and cluster pairs to test
ME.by.neuron.cluster.zscore.modules = ME.by.neuron.cluster.zscore %>% 
  filter(z_ME > 2) %>% 
  dplyr::select(seurat_clusters,
                module)


### compare across social status
## pseudobulk ME across status, genotype, and cluster
# filter to clusters and modules of interest
ME.by.neuron.cluster.geno = burtoni.snseq.combined.sct.neurons.recluster@meta.data %>% 
  group_by(seurat_clusters,
           Genotype.id,
           orig.ident) %>% 
  summarize(
    n_cells = n(),
    across(all_of(mods.neurons), mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    mods.neurons,
    names_to = "module",
    values_to = "mean_ME"
  ) %>% 
  left_join(ME.by.neuron.cluster.zscore.modules %>% 
              mutate(keep = 1)) %>% 
  filter(keep == 1) %>% 
  mutate(module_cluster = paste0(module,
                                 '_',
                                 seurat_clusters)) %>%
  group_by(Genotype.id) %>%
  mutate(total_cells = sum(n_cells),
         rel_n_cells = n_cells / total_cells) %>%
  ungroup()


# graph results
# paper
ME.by.neuron.cluster.geno %>% 
  ggplot(aes(x = seurat_clusters,
             y = mean_ME,
             color = orig.ident,
             size = rel_n_cells)) +
  geom_point(position = position_dodge(width = 0.5),
             shape = 21,
             color = 'black',
             aes(fill = orig.ident)
             # size = 5
             ) +
  theme_classic() +
  theme(legend.position = 'none') +
  scale_fill_manual(values = c("#4e499e",
                               "#60bb46"))+
  scale_color_manual(values = c("#4e499e",
                                "#60bb46"))+
  ylab("Mean ME") +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        strip.text = element_text(size = 20)) +
  xlab('Neuron cluster') +  
  facet_wrap(~module,
             scales = 'free') +
  scale_size_continuous(name = "# nuclei",
                        range = c(2,6))
ggsave('neurons/wgcna/ME by cluster and genotype paper.pdf',
       height = 6,
       width = 6,
       units = "in",
       dpi = 320)

## run weighted lm
# get pvalue and estimate
ME.neuron.lmer.results <- ME.by.neuron.cluster.geno %>%
  group_by(seurat_clusters, module) %>%
  group_modify(~{
    df <- .x
    
    
    fit <- lm(
      mean_ME ~ orig.ident,
      weights = sqrt(rel_n_cells),
      data = df 
      
    )
    
    tidied <- broom.mixed::tidy(fit, effects = "fixed")
    
    tibble(
      estimate = -tidied$estimate[tidied$term == "orig.identsub_burtoni_snseq"],
      p.value = tidied$p.value[tidied$term == "orig.identsub_burtoni_snseq"]
    )
  }) %>%
  ungroup() %>%
  mutate(
    p.adj = p.adjust(p.value, 
                     method = "fdr")
  )

# save results
write_csv(ME.neuron.lmer.results,
          file = 'neurons/neuron_cluster/ME.neuron.lmer.results.csv')






# ## run wilcox test
# # get p-value 
# ME.neuron.wilcox.results <- ME.by.neuron.cluster.geno %>%
#   group_by(seurat_clusters, module) %>%
#   group_modify(~{
#     df <- .x
#     
#     wt <- wilcox.test(
#       mean_ME ~ orig.ident,
#       data = df,
#       exact = FALSE
#     )
#     
#     tibble(
#       p.value = wt$p.value,
#       statistic = wt$statistic
#     )
#   }) %>%
#   ungroup()
# 
# # get effect size
# effect_sizes <- ME.by.neuron.cluster.geno %>%
#   group_by(seurat_clusters, 
#            module, 
#            orig.ident) %>%
#   summarize(mean_ME = mean(mean_ME), 
#             .groups = "drop") %>%
#   pivot_wider(
#     names_from = orig.ident,
#     values_from = mean_ME
#   ) %>%
#   mutate(
#     delta_ME = dom_burtoni_snseq - sub_burtoni_snseq
#   ) %>%
#   dplyr::select(seurat_clusters,
#                 module,
#                 delta_ME)
# 
# # combine
# ME.neuron.wilcox.results = ME.neuron.wilcox.results %>%
#   left_join(effect_sizes,
#             by = c("seurat_clusters",
#                    "module"))
# 
# # adjust pvalue
# ME.neuron.wilcox.results <- ME.neuron.wilcox.results %>%
#   mutate(
#     p.adj = p.adjust(p.value, 
#                      method = "fdr")
#   )
# 
# write_csv(ME.neuron.wilcox.results,
#           file = 'neurons/neuron_cluster/ME.neuron.wilcox.results.csv')

#### save point ####
# save lmer results
write_csv(ME.neuron.lmer.results,
          file = 'neurons/neuron_cluster/ME.neuron.lmer.results.csv')

# # save wilcox results
# write_csv(ME.neuron.wilcox.results,
#           './neurons/wgcna/ME.neuron.wilcox.results')

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
#### compare to average ME ####
## load data
## save glm table
neuron.cluster.glm = read.csv('neurons/neuron_cluster/neuron_seurat_cluster_glm.csv')

# add ME data
neuron.cluster.glm = ME.by.neuron.cluster.zscore %>% 
  filter(z_ME > 2)  %>% 
  mutate(seurat_clusters = as.integer(as.character(seurat_clusters))) %>% 
  right_join(neuron.cluster.glm)

# ## check module eigengene vs significant clusters
# neuron.cluster.glm %>% 
#   filter(contrast != "(Intercept)") %>% 
#   filter(p.adjust.fdr < 0.01) %>% 
#   view()

## save glm table
write_csv(neuron.cluster.glm,
          file = 'neurons/neuron_cluster/neuron_seurat_cluster_glm.csv')

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
library(ggnewscale)

# load all genes for GO comparison
all.genes = burtoni.snseq.combined.sct.neurons.recluster@misc[["neurons"]][["wgcna_genes"]]

# get dataframe of genes and modules
module.genes = burtoni.snseq.combined.sct.neurons.recluster@misc[["neurons"]][["wgcna_modules"]] %>% 
  dplyr::select(gene_name,
                module)

# count number of genes
module.genes %>%
  dplyr::count(module)
# module    n
# 1      grey 1610
# 2     green   91
# 3     brown  205
# 4       red   83
# 5 turquoise  474
# 6    yellow  123
# 7      blue  291
# 8     black   73
# 9      pink   50

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


## get list of module colors
module.colors = module.genes %>% 
  filter(module != 'grey') %>% 
  pull(module) %>% 
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
                           pvalueCutoff = 0.15,
                           pAdjustMethod = 'fdr',
                           minGSSize = 10,
                           maxGSSize = 500,
                           readable = T)
  
  # save results
  enrichGO.results.modules = enrichGO.results.modules %>%
    rbind(tmp.enrichgo@result %>%
                mutate(module = i,
                       ont = 'BP'))

  write.csv(enrichGO.results.modules,
            'neurons/wgcna/go_terms/enrichGO.results.modules.csv',
            row.names = F)
  
  # check number of sig GO
  tmp.num = tmp.enrichgo@result %>% 
    filter(p.adjust < 0.15) %>% 
    nrow()
  
  if (tmp.num != 0) {
  # graph GO results
  # create plot
  tmp.enrichgo.fit <- barplot(tmp.enrichgo,
                                   showCategory = 15) +
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
  
  # create plot
  tmp.enrichgo.dot.fit <- dotplot(tmp.enrichgo,
                              showCategory = 15) +
    ggtitle(paste0(i,
                   ': module enriched BP GO'))
  
  # save plot
  png(paste0('neurons/wgcna/go_terms/enrichplot/',
             i,
             ' module enriched BP GO dotplot.png'),
      res = 300,
      width = 10,
      height = 10,
      units = 'in')
  print(tmp.enrichgo.dot.fit)
  dev.off()
  
  
  # graph simplified semantic networks
  ego <- simplify(
    tmp.enrichgo,
    cutoff = 0.7,
    by = "p.adjust",
    select_fun = min,
    measure = "Wang"
  )
  
  ego <- pairwise_termsim(ego)
  
  png(paste0('neurons/wgcna/go_terms/enrichplot/',
             i,
             ' module enriched BP GO summary network.png'), 
      res = 300, 
      width = 10, 
      height = 10,
      units = 'in')
  emapplot(
    ego,
    showCategory = 50,
    layout = "nicely"
  ) +
    ggtitle(paste0(i, ": module enriched BP GO (summarized)"))
  dev.off()
  
  # create GO term tree
  # https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html#tree-plot
  # download enrichplot
  # remotes::install_github("GuangchuangYu/enrichplot")

  

  
  
  # # graph tree
  # png(paste0('neurons/wgcna/go_terms/enrichplot/',
  #            i,
  #            ' module enriched BP GO treeplot.png'), 
  #     res = 300, 
  #     width = 10, 
  #     height = 10,
  #     units = 'in')
  # tmp.enrichgo %>% 
  #   pairwise_termsim() %>% 
  #   treeplot() +
  #   ggtitle(paste0(i,
  #                  ': module enriched BP GO tree'))
  # dev.off()
  
  
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
   
    # save plot
    png(paste0('neurons/wgcna/go_terms/enrichplot/',
               i,
               ' module enriched BP GO dotplot.png'),
        res = 300,
        width = 10,
        height = 10,
        units = 'in')
    plot.new()
    dev.off()

    # save plot
    png(paste0('neurons/wgcna/go_terms/enrichplot/',
               i,
               ' module enriched BP GO summary network.png'), 
        res = 300, 
        width = 10, 
        height = 10,
        units = 'in')
    plot.new()
    dev.off()
    
    # # create GO term tree
    # # https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html#tree-plot
    # # download enrichplot
    # # remotes::install_github("GuangchuangYu/enrichplot")
    # 
    # # graph tree
    # png(paste0('neurons/wgcna/go_terms/enrichplot/',
    #            i,
    #            ' module enriched BP GO treeplot.png'), 
    #     res = 300, 
    #     width = 10, 
    #     height = 10,
    #     units = 'in')
    # plot.new()
    # dev.off()
  }
}

## check enriched GO terms per module
enrichGO.results.modules %>% 
  filter(p.adjust < 0.15) %>% 
  dplyr::count(module)
# module  n
# 1     black  9
# 2      blue  4
# 3     brown  1
# 4     green 26
# 5      pink  3
# 6 turquoise 38

# paper
### combined figures for paper
# make dummy enrich object that has all genes
tmp = tmp.enrichgo.ensembl

# replace results with results from modules that have less than 5 enriched terms
# tmp@result = enrichGO.results.modules  %>% 
#   filter(module %in% c(enrichGO.results.modules %>% 
#                          filter(p.adjust < 0.15) %>% 
#                          dplyr::count(module) %>% 
#                          filter(n<5) %>% 
#                          pull(module)))

tmp@result = enrichGO.results.modules  %>% 
  group_by(module) %>% 
  top_n(n = 5,
        wt = -p.adjust)  %>% 
  full_join(enrichGO.results.modules  %>% 
              group_by(module) %>% 
              top_n(n = 5,
                    wt = Count)) %>% 
  full_join(
  enrichGO.results.modules %>% 
    filter(module %in% c(enrichGO.results.modules %>%
                           filter(p.adjust < 0.15) %>%
                           dplyr::count(module) %>%
                           filter(n<10) %>%
                           pull(module)))) %>% 
  distinct()

# create dotplot with colors as names
# p = tmp %>% 
#   dotplot(showCategory = 20) + 
#   aes(shape = I(21), 
#       stroke = 5) + 
#   aes(fill = p.adjust,
#       color = module) +
#   scale_color_manual(values = c("blue" = "blue",
#                                 "brown" = "brown",
#                                 "pink" = "pink")) +
#   scale_fill_continuous(low = "red",
#                         high = "blue") +
#   ggtitle(paste0(
#                  'Module enriched BP GO'))

p = tmp %>% 
  dotplot(showCategory = 50) + 
  aes(shape = I(21), 
      stroke = 5) + 
  aes(fill = -log10(p.adjust),
      color = module) +
  scale_color_manual(values = c("blue" = "blue",
                                "brown" = "brown",
                                "pink" = "pink",
                                "black" =  'grey25',
                                'green' = 'green',
                                'turquoise' = 'turquoise')) +
  scale_fill_continuous(low = "red",
                        high = "blue") +
  ggtitle(paste0(
    'Module enriched BP GO')) +
  theme(legend.position = 'inside',
        legend.position.inside = c(0.85,
                                   0.25))
  
# save plot
png(paste0('neurons/wgcna/go_terms/enrichplot/',
           'Low module enriched BP GO dotplot.png'),
    res = 720,
    width = 6.5,
    height = 13,
    units = 'in')
print(p)
dev.off()


# #### Get GO terms for MWU GO old ####
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

# #### prepare for GO_MWU old ####
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


# #### create figures for poster old ####
# ### dotplot of broad cell types across UMAP
# 
# #assign broad cell class
# burtoni.scsorter.data.scsort.output.broad = burtoni.scsorter.data.scsort.output %>% 
#   select(Cell.type,
#          Cell.id) %>%
#   mutate(Cell.class.broad = case_when(Cell.type == 'neurons' ~ 'neurons',
#                                       Cell.type == 'inhibitory' ~ 'neurons',
#                                       Cell.type == 'excitatory' ~ 'neurons',
#                                       Cell.type == 'mural' ~ 'vascular',
#                                       Cell.type == 'VSM' ~ 'vascular',
#                                       Cell.type == 'endothelial' ~ 'vascular',
#                                       Cell.type == 'Unknown' ~ 'unknown',
#                                       TRUE ~ 'glia'))
# 
# ##combine umap with metadata
# burtoni.scsorter.data.scsort.output.umap =  full_join(burtoni.scsorter.data.scsort.output.broad,
#                                                       burtoni.scsorter.data.scsort.output.umap)
# 
# ##graph umap with cell markers
# burtoni.scsorter.data.scsort.output.umap %>% 
#   ggplot(aes(x=UMAP_1,
#              y= UMAP_2,
#              color = Cell.type)) +
#   geom_point() +
#   theme_bw()
# ggsave('scsorter/UMAP scsorter cell types.png',
#        width = 10,
#        height = 10)
# 
# 
# #graph umap with cell markers
# #broad cell class
# burtoni.scsorter.data.scsort.output.umap %>% 
#   ggplot(aes(x=UMAP_1,
#              y= UMAP_2,
#              color = Cell.class.broad)) +
#   geom_point() +
#   theme_classic() +
#   scale_color_manual(values = c('3C51A3',
#                                 'orange3',
#                                 'grey',
#                                 '4AB85C'))
# ggsave('scsorter/UMAP scsorter cell cell types broad.pdf',
#        width = 10,
#        height = 10)
# 
# #graph umap with cell markers
# #broad cell class
# # seperate doms and subs
# burtoni.scsorter.data.scsort.output.umap %>% 
#   ggplot(aes(x=UMAP_1,
#              y= UMAP_2,
#              color = Cell.class.broad)) +
#   geom_point() +
#   theme_classic() +
#   scale_color_manual(values = c('3C51A3',
#                                 'orange3',
#                                 'grey',
#                                 '4AB85C')) +
#   facet_grid( ~ orig.ident) +
#   labs(color='Broad cell classes') +
#   theme(legend.text=element_text(size=20),
#         legend.title=element_text(size=20))+
#   guides(colour = guide_legend(override.aes = list(size=10)))
# ggsave('scsorter/UMAP scsorter cell cell types broad sub vs dom.pdf',
#        width = 12,
#        height = 10)
# 
# 
# 
# 
# 
# #### WGCNA for neurons old ####
# burtoni.snseq.combined.sct.all.neurons = burtoni.snseq.combined.sct.all
# 
# #set idents
# Idents(object = burtoni.snseq.combined.sct.all.neurons) <- "Cell.type"
# 
# #subset to neurons
# burtoni.snseq.combined.sct.all.neurons = subset(burtoni.snseq.combined.sct.all.neurons,
#                                                 idents = c("neurons",
#                                                            "excitatory",
#                                                            "inhibitory"))
# 
# #remove large file
# rm(burtoni.snseq.combined.sct.all)
# 
# # need to set to integrated for clustering
# DefaultAssay(burtoni.snseq.combined.sct.all.neurons) = 'integrated'
# 
# #check data loaded correctly
# ## run PCA, UMAP, and cluster 
# #use 0.8 resolution
# burtoni.snseq.combined.sct.all.neurons.recluster = burtoni.snseq.combined.sct.all.neurons %>% 
#   RunPCA() %>%
#   FindNeighbors(dims = 1:15) %>%
#   RunUMAP(dims = 1:15) %>%
#   FindClusters(resolution = 0.8)
# 
# #remove large file
# rm(burtoni.snseq.combined.sct.all.neurons)
# 
# ## graph 
# # idents to new clusters
# Idents(object = burtoni.snseq.combined.sct.all.neurons.recluster) <- "integrated_snn_res.0.8"
# 
# # #check data
# # #use all neurons 
# # #with variable features (top 2000 variable genes)
# # png('neurons/Clusters.dimplot.neurons.all.png',
# #     width = 10,
# #     height = 10,
# #     units = 'in',
# #     res = 300)
# # DimPlot(burtoni.snseq.combined.sct.all.neurons.recluster,
# #         group.by='integrated_snn_res.0.8',
# #         label=TRUE) +
# #   umap_theme() +
# #   ggtitle('Neurons') +
# #   NoLegend()
# # dev.off()
# # 
# # # for poster
# # #check data
# # #use all neurons
# # #with variable features (top 2000 variable genes)
# # pdf('neurons/Clusters.dimplot.neurons.all.poster.pdf',
# #     width = 5,
# #     height = 5)
# # DimPlot(burtoni.snseq.combined.sct.all.neurons.recluster,
# #         group.by='integrated_snn_res.0.8',
# #         label=TRUE) +
# #   umap_theme() +
# #   ggtitle('Neurons UMAP') +
# #   NoLegend()
# # dev.off()
# # 
# # #check per genotype
# # #use all neurons 
# # #with variable features (top 2000 variable genes)
# # png('neurons/Clusters.dimplot.neurons.all.genotype.png',
# #     width = 10,
# #     height = 10,
# #     units = 'in',
# #     res = 300)
# # DimPlot(burtoni.snseq.combined.sct.all.neurons.recluster,
# #         group.by='Genotype.id') +
# #   umap_theme() +
# #   ggtitle('Neurons')
# # dev.off()
# # 
# # ## create count across genotype per cluster
# # burtoni.snseq.combined.sct.all.neurons.recluster@meta.data %>% 
# #   as.data.frame() %>% 
# #   select(c(Genotype.id,
# #            integrated_snn_res.0.8)) %>% 
# #   table() %>% 
# #   as.data.frame() %>% 
# #   mutate(orig.ident = ifelse(grepl("Dom", 
# #                                    Genotype.id),
# #                              "dom",
# #                              "sub")) %>% 
# #   ggplot(aes(x = integrated_snn_res.0.8,
# #              y = Freq,
# #              color = orig.ident)) +
# #   geom_point() +
# #   theme_classic()
# # ggsave("neurons/Genotype.cluster.neurons.all.count.png",
# #        width = 10,
# #        height = 10)
# 
# # # select 3000 variable genes
# # burtoni.snseq.combined.sct.all.neurons.recluster <- FindVariableFeatures(burtoni.snseq.combined.sct.all.neurons.recluster,
# #                                                                          selection.method = "vst",
# #                                                                          nfeatures = 3000,
# #                                                                          verbose = F)
# 
# # # Identify the 10 most highly variable genes
# # top10.neurons <- head(VariableFeatures(burtoni.snseq.combined.sct.all.neurons.recluster), 
# #                       10)
# 
# # # plot variable features with labels
# # png('neurons/variable_genes/variance by expression.png',
# #     width = 10,
# #     height = 10,
# #     units = 'in',
# #     res = 300)
# # LabelPoints(plot = VariableFeaturePlot(burtoni.snseq.combined.sct.all.neurons.recluster), 
# #             points = top10.neurons, 
# #             repel = TRUE)
# # dev.off()
# # 
# # #graph variable features
# # HVFInfo(object = burtoni.snseq.combined.sct.all.neurons.recluster,
# #         status = TRUE) %>% 
# #   ggplot(aes(x = residual_variance,
# #              fill = variable)) + 
# #   geom_histogram() +
# #   theme_classic() +
# #   ggtitle('Top 2000 variable genes')
# # ggsave('neurons/variable_genes/histogram variable genes.png',
# #        width = 10,
# #        height = 10)
# 
# # HVFInfo(object = burtoni.snseq.combined.sct.all.neurons.recluster,
# #         status = TRUE) %>% 
# #   ggplot(aes(x = residual_variance,
# #              fill = variable)) + 
# #   geom_density() +
# #   theme_classic()+
# #   ggtitle('Top 2000 variable genes')
# # ggsave('neurons/variable_genes/density variable genes.png',
# #        width = 10,
# #        height = 10)
# 
# #gene_select parameter:
# 
# #   variable: use the genes stored in the Seurat object’s VariableFeatures.
# # fraction: use genes that are expressed in a certain fraction of cells for in the whole dataset or in each group of cells, specified by group.by.
# # custom: use genes that are specified in a custom list.
# 
# # burtoni.snseq.combined.sct.all.neurons.recluster <- SetupForWGCNA(
# #   burtoni.snseq.combined.sct.all.neurons.recluster,
# #   gene_select = "variable", # the gene selection approach
# #   wgcna_name = "neurons" # the name of the hdWGCNA experiment
# # )
# 
# 
# burtoni.snseq.combined.sct.all.neurons.recluster = SetupForWGCNA(burtoni.snseq.combined.sct.all.neurons.recluster,
#                                                                  features = VariableFeatures(burtoni.snseq.combined.sct.all.neurons.recluster),
#                                                                  wgcna_name = "neurons"
# )
# 
# # construct metacells  in each group
# # burtoni.snseq.combined.sct.all.neurons.recluster <- MetacellsByGroups(
# #   seurat_obj = burtoni.snseq.combined.sct.all.neurons.recluster,
# #   group.by = c("Genotype.id"), # specify the columns in seurat_obj@meta.data to group by
# #   k = 25, # nearest-neighbors parameter
# #   max_shared = 5, # maximum number of shared cells between two metacells
# #   ident.group = 'Genotype.id', # set the Idents of the metacell seurat object
# #   slot = "data", #set for counts
# #   assay = "integrated", #set assay type
# #   mode = "average", #metacell expression profile determined by average of constituent single cells
# #   min_cells = 50 #set minimum number of cells per metacell
# # )
# 
# burtoni.snseq.combined.sct.all.neurons.recluster <- MetacellsByGroups(
#   seurat_obj = burtoni.snseq.combined.sct.all.neurons.recluster,
#   group.by = c("Genotype.id",
#                "orig.ident"), # specify the columns in seurat_obj@meta.data to group by
#   k = 25, # nearest-neighbors parameter
#   max_shared = 5, # maximum number of shared cells between two metacells
#   ident.group = 'orig.ident', # set the Idents of the metacell seurat object
#   slot = "data", #set for counts
#   assay = "SCT", #set assay type
#   mode = "average", #metacell expression profile determined by average of constituent single cells
#   min_cells = 50 #set minimum number of cells per metacell
# )
# 
# # ## check metacell clusters
# # #cluster metacells
# # seurat_obj = burtoni.snseq.combined.sct.all.neurons.recluster
# # # seurat_obj <- NormalizeMetacells(seurat_obj)
# # seurat_obj <- ScaleMetacells(seurat_obj, 
# #                              features=VariableFeatures(seurat_obj))
# # seurat_obj <- RunPCAMetacells(seurat_obj, 
# #                               features=VariableFeatures(seurat_obj))
# # seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_metacell_obj <- FindNeighbors(seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_metacell_obj, 
# #                                                                                     dims = 1:15)
# # # seurat_obj <- RunHarmonyMetacells(seurat_obj, 
# #                                   # group.by.vars='Genotype.id'
# #                                   # )
# # seurat_obj <- RunUMAPMetacells(seurat_obj,  
# #                                dims=1:15)
# # seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_metacell_obj <- Seurat::FindClusters(seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_metacell_obj, 
# #                                                                                            resolution = 0.6)
# # #graph 
# # # metacell by genotype
# # png('neurons/wgcna/metacell/Metacell UMAP genotype.png',
# #     width = 10,
# #     height = 10,
# #     units = 'in',
# #     res = 300)
# # DimPlotMetacells(seurat_obj,
# #                  group.by = c('Genotype.id')) +
# #   umap_theme() 
# # dev.off()
# # # metacell by status
# # png('neurons/wgcna/metacell/Metacell UMAP social status.png',
# #     width = 10,
# #     height = 10,
# #     units = 'in',
# #     res = 300)
# # DimPlotMetacells(seurat_obj,
# #                  group.by = c('orig.ident')) +
# #   umap_theme() 
# # dev.off()
# # # metacell by cluster
# # png('neurons/wgcna/metacell/Metacell UMAP clusters.png',
# #     width = 10,
# #     height = 10,
# #     units = 'in',
# #     res = 300)
# # DimPlotMetacells(seurat_obj,
# #                  group.by = c('seurat_clusters')) +
# #   umap_theme() 
# # dev.off()
# # # metacell proportion by genotype
# # seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_metacell_obj@meta.data %>% 
# #   select(c(orig.ident,
# #            seurat_clusters)) %>% 
# #   table() %>% 
# #   as.data.frame() %>% 
# #   ggplot(aes(y = orig.ident,
# #              x = seurat_clusters,
# #              fill = Freq,
# #              label = Freq)) +
# #   geom_tile(color = 'black') +
# #   geom_text() + 
# #   scale_fill_gradient(low = "white", high = "red") +
# #   coord_fixed()
# # ggsave('neurons/wgcna/metacell/Metacell heatmap clusters by status.png')
# # # metacell proportion by genotype
# # seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_metacell_obj@meta.data %>% 
# #   select(c(Genotype.id,
# #            seurat_clusters)) %>% 
# #   table() %>% 
# #   as.data.frame() %>% 
# #   ggplot(aes(y = Genotype.id,
# #              x = seurat_clusters,
# #              fill = Freq,
# #              label = Freq)) +
# #   geom_tile(color = 'black') +
# #   geom_text() + 
# #   scale_fill_gradient(low = "white", high = "red") +
# #   coord_fixed()
# # ggsave('neurons/wgcna/metacell/Metacell heatmap clusters by genotype.png')
# 
# # ### clustree
# # # cluster across resolutions
# # seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_metacell_obj <- Seurat::FindClusters(seurat_obj@misc[[seurat_obj@misc$active_wgcna]]$wgcna_metacell_obj, 
# #                                                                               resolution = resolution.range.reduced.2)
# # #set presentation colors
# # presentation.color <- c('#66c2a5',
# #                         '#fc8d62',
# #                         '#8da0cb',
# #                         '#e78ac3',
# #                         '#a6d854',
# #                         '#ffd92f',
# #                         '#e5c494',
# #                         '#b3b3b3')
# # # #clustree
# # DimPlotMetacells(seurat_obj,
# #                  group.by = c('integrated_snn_res.0.8')) +
# #   umap_theme()
# 
# #remove object
# # rm(seurat_obj)
# 
# ## setup expression matrix
# #use metacell expression data
# 
# ##select soft-power threshold
# burtoni.snseq.combined.sct.all.neurons.recluster <- TestSoftPowers(
#   burtoni.snseq.combined.sct.all.neurons.recluster,
#   networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
# )
# 
# # plot the results:
# plot_list.neuron <- PlotSoftPowers(burtoni.snseq.combined.sct.all.neurons.recluster)
# 
# # assemble with patchwork
# png('neurons/wgcna/Soft power threshold.png',
#     width = 10,
#     height = 10,
#     units = 'in',
#     res = 300)
# wrap_plots(plot_list.neuron, 
#            ncol=2)
# dev.off()
# 
# #use soft power threshold 10?
# 
# 
# ### construct co-expression network
# # construct co-expression network:
# burtoni.snseq.combined.sct.all.neurons.recluster <- ConstructNetwork(burtoni.snseq.combined.sct.all.neurons.recluster, 
#                                                                      soft_power=8,
#                                                                      # setDatExpr=FALSE,
#                                                                      tom_name = 'neuron', # name of the topoligical overlap matrix written to disk
#                                                                      overwrite_tom = TRUE
# )
# 
# #graph dendrogram
# png('neurons/wgcna/WGCNA dendrogram.png',
#     width = 10,
#     height = 10,
#     units = 'in',
#     res = 300)
# PlotDendrogram(burtoni.snseq.combined.sct.all.neurons.recluster, 
#                main='neurons hdWGCNA Dendrogram')
# dev.off()
# 
# #graph dendrogram for poster
# # pdf('neurons/wgcna/WGCNA dendrogram poster.pdf',
# #     width = 10,
# #     height = 2.5)
# # PlotDendrogram(burtoni.snseq.combined.sct.all.neurons.recluster,
# #                main='Neurons hdWGCNA Dendrogram')
# # dev.off()
# # ended up exporting the right size with plot viewer
# 
# ##compute module eigengenes
# # need to run ScaleData first or else harmony throws an error:
# # burtoni.snseq.combined.sct.all.neurons.recluster <- ScaleData(burtoni.snseq.combined.sct.all.neurons.recluster, 
# #                                                               features=VariableFeatures(burtoni.snseq.combined.sct.all.neurons.recluster))
# 
# # compute all MEs in the full single-cell dataset
# 
# burtoni.snseq.combined.sct.all.neurons.recluster <- ModuleEigengenes(
#   burtoni.snseq.combined.sct.all.neurons.recluster,
#   exclude_grey = TRUE
#   # , group.by.vars = 'Genotype.id'
# )
# 
# # module eigengenes:
# MEs.neurons <- GetMEs(burtoni.snseq.combined.sct.all.neurons.recluster,
#                       harmonized=FALSE)
# 
# 
# # compute eigengene-based connectivity (kME):
# burtoni.snseq.combined.sct.all.neurons.recluster <- ModuleConnectivity(
#   burtoni.snseq.combined.sct.all.neurons.recluster)
# 
# 
# #### Graph WGCNNA for neurons old ####
# 
# # plot genes ranked by kME for each module
# png('neurons/wgcna/Plot kMEs neurons.png',
#     width = 10,
#     height = 10,
#     units = 'in',
#     res = 300)
# PlotKMEs(burtoni.snseq.combined.sct.all.neurons.recluster, 
#          ncol=7)
# dev.off()
# 
# # correlation between modules
# png('neurons/wgcna/Correlogram modules neurons.png',
#     width = 10,
#     height = 10,
#     units = 'in',
#     res = 300)
# ModuleCorrelogram(burtoni.snseq.combined.sct.all.neurons.recluster,
#                   features = "MEs",
#                   col = rev(COL2('RdBu', 50)),
#                   addCoef.col = 'black',
#                   sig.level = c(0.001), 
#                   pch.cex = 0.9,
#                   insig = 'blank',
#                   diag = FALSE,
#                   order = 'hclust')
# dev.off()
# 
# # make a featureplot of hMEs for each module
# plot_list.neurons.me <- ModuleFeaturePlot(
#   burtoni.snseq.combined.sct.all.neurons.recluster,
#   features='MEs', # plot the MEs
#   order=TRUE # order so the points with highest MEs are on top
# )
# 
# # stitch together with patchwork
# png('neurons/wgcna/UMAP MEs neurons.png',
#     width = 10,
#     height = 10,
#     units = 'in',
#     res = 300)
# wrap_plots(plot_list.neurons.me, ncol=4)
# dev.off()
# 
# # for poster
# # stitch together with patchwork
# # pdf('neurons/wgcna/UMAP MEs neurons poster.pdf',
# #     width = 10,
# #     height = 5)
# # wrap_plots(plot_list.neurons.me, ncol=4)
# # dev.off()
# 
# ## create one umap plot
# # combine umap data and meta data
# # neuron.umap = burtoni.snseq.combined.sct.all.neurons.recluster@reductions$umap@cell.embeddings %>% 
# #   as.data.frame() %>% 
# #   rownames_to_column('Cell.id') %>% 
# #   full_join(burtoni.snseq.combined.sct.all.neurons.recluster@meta.data %>% 
# #               as.data.frame() %>% 
# #               rownames_to_column('Cell.id'))
# 
# # add cluster_color for each cell
# # neuron.umap.long = neuron.umap %>% 
# #   select(c(Cell.id,
# #            UMAP_1,
# #            UMAP_2,
# #            orig.ident,
# #            Genotype.id,
# #            Cell.type,
# #            integrated_snn_res.0.8,
# #            turquoise,
# #            green,
# #            yellow,
# #            red,
# #            pink,
# #            blue,
# #            brown,
# #            black,
# #            grey)) %>% 
# #   pivot_longer(cols = c(turquoise,
# #                         green,
# #                         yellow,
# #                         red,
# #                         pink,
# #                         blue,
# #                         brown,
# #                         black,
# #                         grey),
# #                names_to = 'module',
# #                values_to = 'ME_score') %>% 
# #   group_by(Cell.id) %>% 
# #   mutate(max.ME_score = max(ME_score),
# #          second.max.ME_score = max(ME_score[ME_score != max(ME_score)]),
# #          max.module = ifelse(max.ME_score == ME_score,
# #                              module,
# #                              NA),
# #          second.max.module = ifelse(second.max.ME_score == ME_score,
# #                              module,
# #                              NA),
# #          ratio.max.mes = second.max.ME_score/max.ME_score,
# #          keep = case_when(
# #            ratio.max.mes < 0 ~ max.module,
# #            ratio.max.mes < 0.75 ~ max.module,
# #            ratio.max.mes > 1  ~ max.module,
# #            TRUE  ~ paste(max.module,
# #                          second.max.module,
# #                          sep = ':'),
# #          )) %>% 
# #   filter(!is.na(max.module)) 
# # 
# # # graph umap
# # neuron.umap.long %>% 
# #   ggplot(aes(x = UMAP_1,
# #              y = UMAP_2,
# #              color = max.module)) +
# #   geom_point() +
# #   theme_classic()
# 
# 
# # get mods from object
# mods.neurons <- colnames(MEs.neurons); mods.neurons <- mods.neurons[mods.neurons != 'grey']
# 
# # add MEs to Seurat meta-data:
# burtoni.snseq.combined.sct.all.neurons.recluster@meta.data <- cbind(burtoni.snseq.combined.sct.all.neurons.recluster@meta.data, 
#                                                                     MEs.neurons)
# 
# ##neuropeptides
# modules.neurons <- GetModules(burtoni.snseq.combined.sct.all.neurons.recluster)
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
# 
# # #convert to long
# # modules.long.neurons = modules.neurons %>% 
# #   pivot_longer(kME_turquoise:kME_red ,
# #                names_to = 'module.kME',
# #                values_to = 'kME')
# # 
# # 
# # #graph
# # modules.long.neurons %>% 
# #   filter(neuropeptide == TRUE) %>% 
# #   drop_na() %>%
# #   mutate(module.fct=as.integer(module)) %>% 
# #   ggplot(aes(x = module.kME,
# #              y = reorder(gene_name,
# #                          module.fct),
# #              label = round(kME,
# #                            2))) +
# #   geom_tile(aes(fill = kME)) +
# #   geom_text() +
# #   geom_point(aes(x=-Inf,
# #                  color = module),
# #              size = 5) +
# #   scale_fill_gradientn(colours=c("blue",
# #                                  "white",
# #                                  "red"), 
# #                        limits = c(-0.6, 0.6)) +
# #   theme_classic() +
# #   ylab('Neuropeptides') +
# #   scale_color_manual(values = c("blue"="blue",
# #                                 "grey"="grey",
# #                                 "turquoise"="turquoise",
# #                                 "brown"="brown",
# #                                 "red"="red",
# #                                 "green"="green",
# #                                 "black"="black",
# #                                 "yellow"="yellow")) +
# #   coord_cartesian(clip = 'off')+ 
# #   theme(axis.ticks.y = element_blank())
# # ggsave('neurons/wgcna/Module and neuropeptide heatmap.png',
# #        width = 10,
# #        height = 10)
# 
# ## create dotplot of modules ME and clusters
# DotPlot(burtoni.snseq.combined.sct.all.neurons.recluster,
#         features = c("blue", 
#                      "red", 
#                      "yellow", 
#                      "green", 
#                      "pink", 
#                      "turquoise",
#                      "brown",
#                      "black"),
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
# burtoni.neuron.wgcna.module.cluster = DotPlot.data(burtoni.snseq.combined.sct.all.neurons.recluster,
#                                                    features = c("blue", 
#                                                                 "red", 
#                                                                 "yellow", 
#                                                                 "green", 
#                                                                 "pink", 
#                                                                 "turquoise",
#                                                                 "brown",
#                                                                 "black"),
#                                                    cols = c('grey',
#                                                             'red'),
#                                                    group.by = "integrated_snn_res.0.8",
#                                                    col.min = 2,
#                                                    dot.min = .5,
#                                                    scale = FALSE)
# 
# 
# ## trait correlation
# Idents(burtoni.snseq.combined.sct.all.neurons.recluster) <- burtoni.snseq.combined.sct.all.neurons.recluster$orig.ident
# 
# # convert orig.ident to factor
# burtoni.snseq.combined.sct.all.neurons.recluster$orig.ident.fct <- as.factor(burtoni.snseq.combined.sct.all.neurons.recluster$orig.ident)
# 
# # list of traits to correlate
# cur_traits.neurons <- c('orig.ident.fct', 
#                         'nCount_SCT', 
#                         'nFeature_SCT')
# 
# burtoni.snseq.combined.sct.all.neurons.recluster <- ModuleTraitCorrelation(
#   burtoni.snseq.combined.sct.all.neurons.recluster,
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
#   burtoni.snseq.combined.sct.all.neurons.recluster,
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
# 
# #### Network WGCNA for neurons old ####
# # visualize network with UMAP
# burtoni.snseq.combined.sct.all.neurons.recluster <- RunModuleUMAP(
#   burtoni.snseq.combined.sct.all.neurons.recluster,
#   n_hubs = 5, #number of hub genes to include for the umap embedding
#   n_neighbors=15, #neighbors parameter for umap
#   min_dist=0.3, #min distance between points in umap space
#   spread=5
# )
# 
# 
# # compute cell-type marker genes with Seurat:
# #set idents to cluster
# Idents(burtoni.snseq.combined.sct.all.neurons.recluster) <- burtoni.snseq.combined.sct.all.neurons.recluster$integrated_snn_res.0.8
# 
# #calculate marker genes per cluster
# # with 3000 variable genes
# markers.neuron <- Seurat::FindAllMarkers(
#   burtoni.snseq.combined.sct.all.neurons.recluster,
#   only.pos = TRUE,
#   logfc.threshold=1,
#   features = VariableFeatures(burtoni.snseq.combined.sct.all.neurons.recluster)
# )
# 
# # compute marker gene overlaps
# burtoni.snseq.combined.sct.all.neurons.recluster.overlap_df <- OverlapModulesDEGs(
#   burtoni.snseq.combined.sct.all.neurons.recluster,
#   deg_df = markers.neuron,
#   fc_cutoff = 1 # log fold change cutoff for overlap analysis
# )
# 
# # overlap barplot, produces a plot for each cluster
# plot_list.neuron.overlap <- OverlapBarPlot(burtoni.snseq.combined.sct.all.neurons.recluster.overlap_df)
# 
# # stitch plots with patchwork
# png('neurons/wgcna/Vlnplot modules vs cluster odds ratio.png',
#     width = 10,
#     height = 10,
#     units = 'in',
#     res = 300)
# wrap_plots(plot_list.neuron.overlap, 
#            ncol=4)
# dev.off()
# 
# # plot odds ratio of the overlap as a dot plot
# png('neurons/wgcna/Dotplot modules vs cluster odds ratio.png',
#     width = 10,
#     height = 10,
#     units = 'in',
#     res = 300)
# OverlapDotPlot(burtoni.snseq.combined.sct.all.neurons.recluster.overlap_df,
#                plot_var = 'odds_ratio') +
#   ggtitle('Overlap of modules & cluster markers')
# dev.off()
# 
# ##graph network
# #umap
# png('neurons/wgcna/gene_network/Gene network UMAP.png',
#     width = 10,
#     height = 10,
#     units = 'in',
#     res = 300)
# ModuleUMAPPlot(burtoni.snseq.combined.sct.all.neurons.recluster,
#                edge.alpha=0.5,
#                sample_edges=TRUE,
#                keep_grey_edges=FALSE,
#                edge_prop=0.075, # taking the top 20% strongest edges in each module
#                label_hubs=10 # how many hub genes to plot per module?
# )
# dev.off()
# 
# # for poster
# ##graph network
# #umap
# # pdf('neurons/wgcna/gene_network/Gene network UMAP poster.pdf',
# #     width = 10,
# #     height = 10)
# # ModuleUMAPPlot.size(burtoni.snseq.combined.sct.all.neurons.recluster,
# #                edge.alpha= 0.75,
# #                sample_edges=TRUE,
# #                keep_grey_edges=FALSE,
# #                edge_prop=0.075, # taking the top 20% strongest edges in each module
# #                label_hubs=0, # how many hub genes to plot per module?
# #                vertex.label.cex = 1,
# #                dot.size = 8,
# #                edge.size = 4
# # )
# # dev.off()
# 
# # module specific network
# ModuleNetworkPlot(burtoni.snseq.combined.sct.all.neurons.recluster,
#                   outdir = 'neurons/wgcna/gene_network/')
# 
# # hubgene network
# png('neurons/wgcna/gene_network/Hubgene network.png',
#     width = 10,
#     height = 10,
#     units = 'in',
#     res = 300)
# HubGeneNetworkPlot(burtoni.snseq.combined.sct.all.neurons.recluster,
#                    n_hubs = 5, 
#                    n_other=50,
#                    edge_prop = 0.75,
#                    mods = 'all'
# )
# dev.off()
# 
# 
# 
# 
# 
# 
# #### DMEs WGCNA for neurons old ####
# ## Differential module
# # get list of dom and sub cells
# group1.neurons.dom <- burtoni.snseq.combined.sct.all.neurons.recluster@meta.data %>% 
#   subset(orig.ident == 'dom_burtoni_snseq') %>%
#   rownames
# group2.neurons.sub <- burtoni.snseq.combined.sct.all.neurons.recluster@meta.data %>% 
#   subset(orig.ident == 'sub_burtoni_snseq') %>% 
#   rownames
# 
# # calculate differential module eigengene
# DMEs.neurons <- FindDMEs(
#   burtoni.snseq.combined.sct.all.neurons.recluster,
#   barcodes1 = group1.neurons.dom,
#   barcodes2 = group2.neurons.sub,
#   test.use='wilcox',
#   harmonized = FALSE
# )
# 
# # graph DME
# PlotDMEsVolcano(
#   burtoni.snseq.combined.sct.all.neurons.recluster,
#   DMEs.neurons) +
#   theme_classic() +
#   xlim(-0.5, 0.5)
# ggsave('neurons/wgcna/ME by cluster volcanoplot.png',
#        width = 10,
#        height = 10)
# 
# # for poster
# # create mod color list
# mod_colors <- DMEs.neurons$module
# names(mod_colors) <- as.character(DMEs.neurons$module)
# # graph DME
# DMEs.neurons %>%
#   mutate(anno = ifelse(p_val_adj < 0.05,
#                        module,
#                        "")) %>%
#   ggplot(aes(x = avg_log2FC,
#              y = -log10(p_val_adj)))  +
#   geom_rect(aes(xmin = -Inf,
#                 xmax = Inf,
#                 ymin = -Inf,
#                 ymax = -log10(0.05)),
#             fill = "grey75",
#             alpha = 0.8,
#             color = NA) +
#   geom_vline(xintercept = 0,
#              linetype = "dashed",
#              color = "grey75",
#              alpha = 0.8) +
#   geom_point(aes(fill = module),
#              size = 8,
#              pch = 21,
#              color = "black") +
#   geom_text_repel(aes(label = anno),
#                   color = "black",
#                   min.segment.length = 0,
#                   max.overlaps = Inf,
#                   size = 8,
#                   point.padding = 4,
#                   nudge_x = 0.06) +
#   theme_classic() +
#   theme(legend.position = 'none') +
#   xlim(c(-0.45,
#          0.45)) +
#   scale_fill_manual(values = mod_colors)+
#   xlab(bquote("Average log"[2] ~ "(Fold Change)")) +
#   ylab(bquote("-log"[10] ~ "(Adj. P-value)")) +
#   ggtitle('Differential module eigengene analysis') +
#   theme(panel.border = element_rect(color = "black",
#                                     fill = NA,
#                                     size = 1))+
#   theme(axis.text = element_text(size = 15))  +
#   theme(axis.title = element_text(size = 20))+
#   theme(plot.title = element_text(size=20))
# ggsave('neurons/wgcna/ME by cluster volcanoplot poster.pdf',
#        width = 10,
#        height = 5,
#        units = "in",
#        dpi = 320)
# 
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
#   burtoni.snseq.combined.sct.all.neurons.recluster@meta.data %>% 
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
#   burtoni.snseq.combined.sct.all.neurons.recluster@meta.data %>% 
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
#   burtoni.snseq.combined.sct.all.neurons.recluster@meta.data %>% 
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
# 
# #### marker genes, hub genes, and go terms WGCNA for neurons old ####
# hub.genes.neurons = GetModules(burtoni.snseq.combined.sct.all.neurons.recluster)
# 
# write_csv(hub.genes.neurons,
#           './neurons/wgcna/hub.genes.neurons.csv')
# 
# ## get list of marker genes
# write_csv(markers.neuron,
#           './neurons/markers.neuron.csv')
# 
# # get go terms for cluster 7
# 
# # save hub genes for GO analysis with:
# # to get GO term accession use biomart: http://useast.ensembl.org/biomart/martview/
# 
# ## load module go term data
# # brown module
# brown_go_terms = read.csv('./neurons/wgcna/brown_go_terms.csv')
# # black module
# black_go_terms = read.csv('./neurons/wgcna/black_go_terms.csv')
# 
# # marker neurons
# markers.neuron_go_terms = read.csv('./neurons/markers.neuron_go_terms.csv')
# 
# ## combine go terms with kME scores
# #brown
# brown_go_terms.kme = brown_go_terms %>% 
#   left_join(hub.genes.neurons %>% 
#               filter(module == 'brown') %>% 
#               select(c(gene_name,
#                        module,
#                        kME_brown)),
#             by = c("Gene.name" = "gene_name")) %>% 
#   filter(!is.na(kME_brown)) %>% 
#   rbind(brown_go_terms %>% 
#           left_join(hub.genes.neurons %>% 
#                       filter(module == 'brown') %>% 
#                       select(c(gene_name,
#                                module,
#                                kME_brown)),
#                     by = c("Gene.stable.ID" = "gene_name")) %>% 
#           filter(!is.na(kME_brown)))
# 
# #black
# black_go_terms.kme = black_go_terms %>% 
#   left_join(hub.genes.neurons %>% 
#               filter(module == 'black') %>% 
#               select(c(gene_name,
#                        module,
#                        kME_black)),
#             by = c("Gene.name" = "gene_name")) %>% 
#   filter(!is.na(kME_black)) %>% 
#   rbind(black_go_terms %>% 
#           left_join(hub.genes.neurons %>% 
#                       filter(module == 'black') %>% 
#                       select(c(gene_name,
#                                module,
#                                kME_black)),
#                     by = c("Gene.stable.ID" = "gene_name")) %>% 
#           filter(!is.na(kME_black)))
# 
# ## combine marker go terms with pval
# #marker
# markers.neuron_go_terms.sig = markers.neuron_go_terms %>% 
#   left_join(markers.neuron,
#             by = c("Gene.name" = "gene"))  %>% 
#   filter(!is.na(p_val)) %>% 
#   rbind(markers.neuron_go_terms %>% 
#           left_join(markers.neuron,
#                     by = c("Gene.stable.ID" = "gene"))) %>% 
#   filter(!is.na(p_val)) %>% 
#   relocate(GO.term.accession, 
#            .after = last_col())%>% 
#   relocate(p_val_adj, 
#            .after = last_col())
# 
# ## increase kME to positive values
# # move GO terms and kME to end
# #brown
# brown_go_terms.kme = brown_go_terms.kme %>% 
#   mutate(kME_brown.scaled = kME_brown + 
#            abs(min(kME_brown)) + 
#            0.0001) %>% 
#   relocate(GO.term.accession, 
#            .after = last_col())%>% 
#   relocate(kME_brown.scaled, 
#            .after = last_col())
# 
# #black
# black_go_terms.kme = black_go_terms.kme %>% 
#   mutate(kME_black.scaled = kME_black + 
#            abs(min(kME_black)) + 
#            0.0001) %>% 
#   relocate(GO.term.accession, 
#            .after = last_col())%>% 
#   relocate(kME_black.scaled, 
#            .after = last_col())
# 
# ## create list of top 25 genes
# #brown
# brown_go_terms.kme.top25 = brown_go_terms %>% 
#   left_join(hub.genes.neurons %>% 
#               filter(module == 'brown') %>% 
#               top_n(25,
#                     kME_brown) %>% 
#               select(c(gene_name,
#                        module,
#                        kME_brown)),
#             by = c("Gene.name" = "gene_name")) %>% 
#   filter(!is.na(kME_brown)) %>% 
#   rbind(brown_go_terms %>% 
#           left_join(hub.genes.neurons %>% 
#                       filter(module == 'brown') %>% 
#                       top_n(25,
#                             kME_brown) %>% 
#                       select(c(gene_name,
#                                module,
#                                kME_brown)),
#                     by = c("Gene.stable.ID" = "gene_name")) %>% 
#           filter(!is.na(kME_brown)))  %>% 
#   relocate(GO.term.accession, 
#            .after = last_col())%>% 
#   relocate(kME_brown, 
#            .after = last_col())
# 
# #black
# black_go_terms.kme.top25 = black_go_terms %>% 
#   left_join(hub.genes.neurons %>% 
#               filter(module == 'black') %>% 
#               top_n(25,
#                     kME_black) %>% 
#               select(c(gene_name,
#                        module,
#                        kME_black)),
#             by = c("Gene.name" = "gene_name")) %>% 
#   filter(!is.na(kME_black)) %>% 
#   rbind(black_go_terms %>% 
#           left_join(hub.genes.neurons %>% 
#                       filter(module == 'black') %>% 
#                       top_n(25,
#                             kME_black) %>% 
#                       select(c(gene_name,
#                                module,
#                                kME_black)),
#                     by = c("Gene.stable.ID" = "gene_name")) %>% 
#           filter(!is.na(kME_black))) %>% 
#   relocate(GO.term.accession, 
#            .after = last_col())%>% 
#   relocate(kME_black, 
#            .after = last_col())
# 
# ## create list for revigo
# #brown
# brown_go_terms.kme.top25.filter = brown_go_terms.kme.top25 %>% 
#   select(c(GO.term.accession,
#            kME_brown)) %>% 
#   filter(GO.term.accession != "") %>% 
#   distinct() 
# 
# #black
# black_go_terms.kme.top25.filter = black_go_terms.kme.top25 %>% 
#   select(c(GO.term.accession,
#            kME_black)) %>% 
#   filter(GO.term.accession != "") %>% 
#   distinct() 
# 
# 
# 
# # to look for enrichment and create figures use GO term accession and kME (higher is better) for each module: http://revigo.irb.hr/
# # remember to sort by kME and filter out blank GO terms before copying to revigo
# # all module genes
# #brown
# write_csv(brown_go_terms.kme,
#           './neurons/wgcna/go_terms/brown_go_terms.kme.csv')
# #black
# write_csv(black_go_terms.kme,
#           './neurons/wgcna/go_terms/black_go_terms.kme.csv')
# 
# #top 25 hub genes filtered
# #brown
# write_csv(brown_go_terms.kme.top25.filter,
#           './neurons/wgcna/go_terms/brown_go_terms.kme.top25.filter.csv')
# #black
# write_csv(black_go_terms.kme.top25.filter,
#           './neurons/wgcna/go_terms/black_go_terms.kme.top25.filter.csv')
# 
# #marker
# write_csv(markers.neuron_go_terms.sig,
#           './neurons/markers.neuron_go_terms.sig.csv')
# #use 0.05 p_val_adj cutoff
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
# #### create figures for poster old ####
# ### dotplot of broad cell types across UMAP
# 
# #assign broad cell class
# burtoni.scsorter.data.scsort.output.broad = burtoni.scsorter.data.scsort.output %>% 
#   select(Cell.type,
#          Cell.id) %>%
#   mutate(Cell.class.broad = case_when(Cell.type == 'neurons' ~ 'neurons',
#                                       Cell.type == 'inhibitory' ~ 'neurons',
#                                       Cell.type == 'excitatory' ~ 'neurons',
#                                       Cell.type == 'mural' ~ 'vascular',
#                                       Cell.type == 'VSM' ~ 'vascular',
#                                       Cell.type == 'endothelial' ~ 'vascular',
#                                       Cell.type == 'Unknown' ~ 'unknown',
#                                       TRUE ~ 'glia'))
# 
# ##combine umap with metadata
# burtoni.scsorter.data.scsort.output.umap =  full_join(burtoni.scsorter.data.scsort.output.broad,
#                                                       burtoni.scsorter.data.scsort.output.umap)
# 
# ##graph umap with cell markers
# burtoni.scsorter.data.scsort.output.umap %>% 
#   ggplot(aes(x=UMAP_1,
#              y= UMAP_2,
#              color = Cell.type)) +
#   geom_point() +
#   theme_bw()
# ggsave('scsorter/UMAP scsorter cell types.png',
#        width = 10,
#        height = 10)
# 
# 
# #graph umap with cell markers
# #broad cell class
# burtoni.scsorter.data.scsort.output.umap %>% 
#   ggplot(aes(x=UMAP_1,
#              y= UMAP_2,
#              color = Cell.class.broad)) +
#   geom_point() +
#   theme_classic() +
#   scale_color_manual(values = c('3C51A3',
#                                 'orange3',
#                                 'grey',
#                                 '4AB85C'))
# ggsave('scsorter/UMAP scsorter cell cell types broad.pdf',
#        width = 10,
#        height = 10)
# 
# #graph umap with cell markers
# #broad cell class
# # seperate doms and subs
# burtoni.scsorter.data.scsort.output.umap %>% 
#   ggplot(aes(x=UMAP_1,
#              y= UMAP_2,
#              color = Cell.class.broad)) +
#   geom_point() +
#   theme_classic() +
#   scale_color_manual(values = c('3C51A3',
#                                 'orange3',
#                                 'grey',
#                                 '4AB85C')) +
#   facet_grid( ~ orig.ident) +
#   labs(color='Broad cell classes') +
#   theme(legend.text=element_text(size=20),
#         legend.title=element_text(size=20))+
#   guides(colour = guide_legend(override.aes = list(size=10)))
# ggsave('scsorter/UMAP scsorter cell cell types broad sub vs dom.pdf',
#        width = 12,
#        height = 10)
# 
# 
# 
# 
