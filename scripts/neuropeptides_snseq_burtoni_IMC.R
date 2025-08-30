#### Burtoni snseq seurat analysis
### neuropeptides
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
library(tidyverse)


#### load data ####
### single cell data
load('burtoni.snseq.combined.sct.all.RData')

### scsorter data
load("burtoni.scsorter.data.scsort.output.RData")

### Add metadata 
burtoni.snseq.combined.sct.all = AddMetaData(
  object = burtoni.snseq.combined.sct.all,
  metadata = burtoni.scsorter.data.scsort.output %>% 
    select(Cell.type,
           Cell.id) %>% 
    column_to_rownames(var = "Cell.id"),
  col.name = 'Cell.type'
)


# resolution range
resolution.range <- seq(from = 0, to = 2, by = 0.2)

#reduced range
resolution.range.reduced <- seq(from = 0, to = 1.2, by = 0.2)

#reduced range 2
resolution.range.reduced.2 <- seq(from = 0.2, to = 1, by = 0.1)

### neuropeptide list
neuropeptides.list = read_csv("../Gene.lists/neuropeptides.list_orthologs.csv")

#### recluster without vascular ####

burtoni.snseq.combined.sct.all.subset = burtoni.snseq.combined.sct.all

#set idents to broad cell class
Idents(object = burtoni.snseq.combined.sct.all.subset) <- "Cell.class.broad"

#subset to neurons
burtoni.snseq.combined.sct.all.subset = subset(burtoni.snseq.combined.sct.all.subset,
                                               idents = c("neurons",
                                                          "glia"))

#set idents to cell type
Idents(object = burtoni.snseq.combined.sct.all.subset) <- "Cell.type"

## run PCA, UMAP, and cluster 
#use 0.4 resolution
burtoni.snseq.combined.sct.all.subset.recluster = burtoni.snseq.combined.sct.all.subset %>% 
  RunPCA() %>%
  FindNeighbors(dims = 1:15) %>%
  RunUMAP(dims = 1:15) %>%
  FindClusters(resolution = 0.8)


### clustree
# cluster across resolutions
burtoni.snseq.combined.sct.all.subset.clustree <- Seurat::FindClusters(object = burtoni.snseq.combined.sct.all.subset.recluster, 
                                                                       resolution = resolution.range.reduced)

#clustree
clustree(burtoni.snseq.combined.sct.all.subset.clustree, 
         prefix = "integrated_snn_res.",
         scale_node_text = TRUE) +
  scale_edge_color_continuous(low = "black", 
                              high = "black") +
  theme(legend.position = "bottom")
ggsave('cell.subset/cell.subset.clustree.png',
       width = 10,
       height = 10)


### graph 
# # idents to new clusters
# Idents(object = burtoni.snseq.combined.sct.all.subset.recluster) <- "integrated_snn_res.0.8"
# 
# ## hover locator
# DimPlot(burtoni.snseq.combined.sct.all.subset.recluster, 
#         reduction = "umap", 
#         label = TRUE,
#         repel = TRUE, 
#         group.by = "SCT_snn_res.0.8") %>% 
#   HoverLocator(information = FetchData(object = burtoni.snseq.combined.sct.all.subset.recluster, 
#                                        vars = 'SCT_snn_res.0.8'))

##dom vs sub
DimPlot(burtoni.snseq.combined.sct.all.subset.recluster, 
        reduction = "umap", 
        group.by = "orig.ident")
ggsave('cell.subset/DomVsSub.dimplot.cell.subset.png')

## clusters
DimPlot(burtoni.snseq.combined.sct.all.subset.recluster, 
        reduction = "umap", 
        label = TRUE,
        repel = TRUE, 
        group.by = "integrated_snn_res.0.8")
ggsave('cell.subset/Clusters.dimplot.cell.subset.png',
       width = 10,
       height = 10)

## cell type
DimPlot(burtoni.snseq.combined.sct.all.subset.recluster, 
        reduction = "umap", 
        label = TRUE,
        repel = TRUE, 
        group.by = "Cell.type")
ggsave('cell.subset/Cell.type.dimplot.cell.subset.png',
       width = 10,
       height = 10)

## broad cell class
DimPlot(burtoni.snseq.combined.sct.all.subset.recluster, 
        reduction = "umap", 
        label = TRUE,
        repel = TRUE, 
        group.by = "Cell.class.broad")
ggsave('cell.subset/Cell.class.broad.dimplot.cell.subset.png',
       width = 10,
       height = 10)


##expression level
burtoni.snseq.combined.sct.all.subset.recluster.expression = full_join(full_join(burtoni.snseq.combined.sct.all.subset.recluster@reductions$umap@cell.embeddings %>% 
                                                                                   as.data.frame() %>% 
                                                                                   rownames_to_column("Cell.id"),
                                                                                 burtoni.snseq.combined.sct.all.subset.recluster@meta.data %>% 
                                                                                   rownames_to_column("Cell.id")),
                                                                       burtoni.snseq.combined.sct.all.subset.recluster@assays$integrated@scale.data %>% 
                                                                         as.data.frame()  %>% 
                                                                         filter(rownames(burtoni.snseq.combined.sct.all.nucb2.neurons.recluster@assays$integrated@scale.data) %in% c('avp')) %>% 
                                                                         t() %>% 
                                                                         as.data.frame() %>% 
                                                                         rownames_to_column('Cell.id')) %>% 
  full_join(burtoni.snseq.combined.sct.all.subset.recluster@reductions$pca@cell.embeddings %>% 
              as.data.frame() %>% 
              select(PC_1,
                     PC_2,
                     PC_3,
                     PC_4,
                     PC_5) %>% 
              rownames_to_column("Cell.id"))

# seperation by cell type
table(burtoni.snseq.combined.sct.all.subset.recluster.expression %>% 
        select(Cell.type,
               integrated_snn_res.0.8))

# seperation by cell class broad
table(burtoni.snseq.combined.sct.all.subset.recluster.expression %>% 
        select(Cell.class.broad,
               integrated_snn_res.0.8))

# difference social status
#bias
# 12500 to 10004
table(burtoni.snseq.combined.sct.all.subset.recluster.expression %>% 
        select(orig.ident))

#compare across clusters
table(burtoni.snseq.combined.sct.all.subset.recluster.expression %>% 
        select(integrated_snn_res.0.8,
               orig.ident)) %>% 
  as.data.frame() %>% 
  pivot_wider(names_from = orig.ident,
              values_from = Freq) %>% 
  mutate(Dom.sub.ratio = dom_burtoni_snseq/sub_burtoni_snseq) %>% 
  filter(Dom.sub.ratio > 1.25) 

#graph
table(burtoni.snseq.combined.sct.all.subset.recluster.expression %>% 
        select(integrated_snn_res.0.8,
               orig.ident)) %>% 
  as.data.frame() %>% 
  pivot_wider(names_from = orig.ident,
              values_from = Freq) %>% 
  mutate(Dom.sub.ratio = dom_burtoni_snseq/sub_burtoni_snseq,
         Total = dom_burtoni_snseq + sub_burtoni_snseq) %>% 
  ggplot(aes(x = Dom.sub.ratio,
             y=Total,
             label = integrated_snn_res.0.8)) +
  geom_label() +
  theme_classic() +
  geom_vline(xintercept = 1.25)
ggsave('cell.subset/Cluster.bias.DomvsSub.cell.subset.png',
       width = 10,
       height = 10)

## cell class across clusters
table(burtoni.snseq.combined.sct.all.subset.recluster.expression %>% 
        select(integrated_snn_res.0.8,
               Cell.class.broad)) %>% 
  as.data.frame() %>% 
  ggplot(aes(x = Freq,
             y = integrated_snn_res.0.8,
             fill = Cell.class.broad)) +
  geom_bar(stat = 'identity') +
  theme_classic() 
ggsave('cell.subset/Cluster.cell.class.broad.cell.subset.png',
       width = 10,
       height = 10)

# cell class across clusters
table(burtoni.snseq.combined.sct.all.subset.recluster.expression %>% 
        select(integrated_snn_res.0.8,
               Cell.class.broad)) %>% 
  as.data.frame() %>% 
  group_by(integrated_snn_res.0.8) %>% 
  mutate(Total.cluster = sum(Freq)) %>% 
  ungroup() %>% 
  mutate(Percent = 100*Freq/Total.cluster) %>% 
  ggplot(aes(x = Percent,
             y = integrated_snn_res.0.8,
             fill = Cell.class.broad)) +
  geom_bar(stat = 'identity') +
  theme_classic() 
ggsave('cell.subset/Cluster.cell.class.broad.cell.subset.percent.png',
       width = 10,
       height = 10)


#### neurons ####

burtoni.snseq.combined.sct.all.neurons = burtoni.snseq.combined.sct.all

#set idents
Idents(object = burtoni.snseq.combined.sct.all.neurons) <- "Cell.type"

#subset to neurons
burtoni.snseq.combined.sct.all.neurons = subset(burtoni.snseq.combined.sct.all.neurons,
                                                idents = c("neurons",
                                                           "excitatory",
                                                           "inhibitory"))

# clustering with integrated
DefaultAssay(burtoni.snseq.combined.sct.all.neurons) = 'integrated'

## run PCA, UMAP, and cluster 
#use 0.4 resolution
burtoni.snseq.combined.sct.all.neurons.recluster = burtoni.snseq.combined.sct.all.neurons %>% 
  RunPCA() %>%
  FindNeighbors(dims = 1:15) %>%
  RunUMAP(dims = 1:15) %>%
  FindClusters(resolution = 0.8)

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
burtoni.snseq.combined.sct.all.neurons.recluster.expression = full_join(full_join(burtoni.snseq.combined.sct.all.neurons.recluster@reductions$umap@cell.embeddings %>% 
                                                                                    as.data.frame() %>% 
                                                                                    rownames_to_column("Cell.id"),
                                                                                  burtoni.snseq.combined.sct.all.neurons.recluster@meta.data %>% 
                                                                                    rownames_to_column("Cell.id")),
                                                                        burtoni.snseq.combined.sct.all.neurons.recluster@assays$integrated@scale.data %>% 
                                                                          as.data.frame() %>% 
                                                                          filter(rownames(burtoni.snseq.combined.sct.all.neurons.recluster@assays$integrated@scale.data) %in% c('avp',
                                                                                                                                                                                'oxt')) %>% 
                                                                          t() %>% as.data.frame() %>% 
                                                                          rownames_to_column('Cell.id')) %>% 
  full_join(burtoni.snseq.combined.sct.all.neurons.recluster@reductions$pca@cell.embeddings %>% 
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
                                  ifelse(avp > 2,
                                         'avp',
                                         NA))))

# # seperation by cell type
# table(burtoni.snseq.combined.sct.all.neurons.recluster.expression %>% 
#         select(Cell.type,
#                integrated_snn_res.0.8))
# 
# # difference social status
# #bias
# # 8181 to 6411
# table(burtoni.snseq.combined.sct.all.neurons.recluster.expression %>% 
#         +         select(orig.ident))
# #compare across clusters
# table(burtoni.snseq.combined.sct.all.neurons.recluster.expression %>% 
#         select(integrated_snn_res.0.8,
#                orig.ident)) %>% 
#   as.data.frame() %>% 
#   pivot_wider(names_from = orig.ident,
#               values_from = Freq) %>% 
#   mutate(Dom.sub.ratio = dom_burtoni_snseq/sub_burtoni_snseq) %>% 
#   filter(Dom.sub.ratio > 1.33) 
# 
# 
# table(burtoni.snseq.combined.sct.all.neurons.recluster.expression %>% 
#         select(integrated_snn_res.0.8,
#                orig.ident)) %>% 
#   as.data.frame() %>% 
#   pivot_wider(names_from = orig.ident,
#               values_from = Freq) %>% 
#   mutate(Dom.sub.ratio = dom_burtoni_snseq/sub_burtoni_snseq,
#          Total = dom_burtoni_snseq + sub_burtoni_snseq) %>% 
#   ggplot(aes(x = Dom.sub.ratio,
#              y=Total)) +
#   geom_point() +
#   theme_classic() +
#   geom_vline(xintercept = 1.28)

# #avp and oxt
# table(burtoni.snseq.combined.sct.all.neurons.recluster.expression %>% 
#         select(Oxt.cell,
#                integrated_snn_res.0.8)) %>% View()





# 
# VlnPlot(burtoni.snseq.combined.sct.all, 
#         features = c("ENSONIG00000033642"),
#         group.by = "orig.ident") +
#   ylim(1,4)

#### neuropeptides ####
## neuropeptides
#avp vs oxt
## umap plots
FeaturePlot(object = burtoni.snseq.combined.sct.all, 
            features = c("oxt",
                         "avp"),
            pt.size = 0.5,
            min.cutoff = 4,
            max.cutoff = 5,
            order = TRUE,
            blend = T)
ggsave('avp.oxt/Avp.Oxt.umap.png',
       width = 20,
       height = 20)
#across samples
FeaturePlot(object = burtoni.snseq.combined.sct.all, 
            features = c("oxt",
                         "avp"),
            pt.size = 0.5,
            min.cutoff = 4,
            max.cutoff = 6,
            order = TRUE,
            blend = T,
            split.by = "orig.ident")
ggsave('avp.oxt/Avp.Oxt.subvsdom.umap.png',
       width = 10,
       height = 10)



##vinplot
#avp
VlnPlot(burtoni.snseq.combined.sct.all, 
        features = c("avp"),
        split.by = "orig.ident") +
  ylim(4,25)
ggsave('avp.oxt/Avp.subvsdom.vlnplot.png',
       width = 10,
       height = 10)
#oxt
VlnPlot(burtoni.snseq.combined.sct.all, 
        features = c("oxt"),
        split.by = "orig.ident") +
  ylim(4,25)
ggsave('avp.oxt/Oxt.subvsdom.vlnplot.png',
       width = 10,
       height = 10)


##check number of cells expressing gene above threshold
#oxt
#141
sum(GetAssayData(object = burtoni.snseq.combined.sct.all, 
                 slot = "data")["oxt",]>4)
#avp
#157
sum(GetAssayData(object = burtoni.snseq.combined.sct.all, 
                 slot = "data")["avp",]>4)


#create list of oxt cells?
# oxt.cell.ids = WhichCells(burtoni.snseq.combined.sct.all,
#                           subset.name = )

### subset data to relevant cells
burtoni.snseq.combined.sct.all.avp.oxt = subset(burtoni.snseq.combined.sct.all,
                                                subset = avp > 4 | oxt > 4)

## run PCA, UMAP, and cluster 
#use 0.6 resolution
burtoni.snseq.combined.sct.all.avp.oxt.recluster = burtoni.snseq.combined.sct.all.avp.oxt %>% 
  RunPCA() %>%
  FindNeighbors(dims = 1:15) %>%
  RunUMAP(dims = 1:15) %>%
  FindClusters(resolution = 1)

## check cells per group per cluster
#get table
Avp.cells.per.cluster.table = table(burtoni.snseq.combined.sct.all.avp.oxt.recluster@active.ident, 
                                    burtoni.snseq.combined.sct.all.avp.oxt.recluster@meta.data$orig.ident) %>% 
  as.data.frame.matrix() %>% 
  rownames_to_column("cluster.id") %>% 
  pivot_longer(cols = c("dom_burtoni_snseq",
                        "sub_burtoni_snseq"),
               names_to = "orig.ident",
               values_to = "count") %>% 
  group_by(orig.ident) %>% 
  mutate(total = sum(count)) %>% 
  ungroup() %>% 
  mutate(percentage = 100*count/total)

#graph
Avp.cells.per.cluster.table %>% 
  as.data.frame() %>% 
  mutate(cluster.id = as.numeric(cluster.id)) %>% 
  ggplot(aes(x = cluster.id,
             y = percentage,
             group = orig.ident,
             color = orig.ident)) +
  geom_point() +
  theme_bw()
ggsave('avp.oxt/Cells.per.cluster.DomvsSub.avp.oxt.png',
       width = 10,
       height = 10)

### graph 
##dom vs sub
DimPlot(burtoni.snseq.combined.sct.all.avp.oxt.recluster, 
        reduction = "umap", 
        group.by = "orig.ident")
# ggsave('avp.oxt/DomVsSub.dimplot.avp.oxt.all.png')

## clusters
DimPlot(burtoni.snseq.combined.sct.all.avp.oxt.recluster, 
        reduction = "umap", 
        label = TRUE,
        repel = TRUE)
ggsave('avp.oxt/Clusters.dimplot.avp.oxt.all.png',
       width = 10,
       height = 10)

#avp vs oxt
## umap plots
FeaturePlot(object = burtoni.snseq.combined.sct.all.avp.oxt.recluster, 
            features = c("oxt",
                         "avp"),
            pt.size = 0.5,
            max.cutoff = 6,
            order = TRUE,
            label = T,
            blend = T)
ggsave('avp.oxt/Avp.Oxt.clustering.umap.png',
       width = 15,
       height = 7)

#across samples
DimPlot(burtoni.snseq.combined.sct.all.avp.oxt.recluster, 
        reduction = "umap",
        group.by = 'orig.ident',
        pt.size = 1) +
  theme(legend.position = 'none') 
ggsave('avp.oxt/Avp.Oxt.clustering.umap.DomvsSub.png',
       width = 10,
       height = 10)


##pca
# clusters
DimPlot(burtoni.snseq.combined.sct.all.avp.oxt.recluster, 
        reduction = "pca", 
        label = TRUE,
        repel = TRUE)
ggsave('avp.oxt/Avp.Oxt.pca.PC1.PC2.png',
       width = 10,
       height = 10)
#PC 2 and 3
DimPlot(burtoni.snseq.combined.sct.all.avp.oxt.recluster, 
        reduction = "pca", 
        label = TRUE,
        repel = TRUE,
        dims = c(3,2))
ggsave('avp.oxt/Avp.Oxt.pca.PC2.PC3.png',
       width = 10,
       height = 10)
#avp and oxt
FeaturePlot(object = burtoni.snseq.combined.sct.all.avp.oxt.recluster, 
            features = c("oxt"),
            pt.size = 3,
            max.cutoff = 6,
            order = TRUE,
            label = T,
            reduction = 'pca',
            cols = c('blue',
                     'red'),
            dims = c(3,2))+
  theme(text = element_text(size = 30))
ggsave('avp.oxt/Avp.Oxt.pca.PC2.PC3.expression.png',
       width = 10,
       height = 10)

#avp and oxt
#no label
FeaturePlot(object = burtoni.snseq.combined.sct.all.avp.oxt.recluster, 
            features = c("avp"),
            pt.size = 3,
            max.cutoff = 6,
            order = TRUE,
            reduction = 'pca',
            cols = c('blue',
                     'red'),
            dims = c(3,2))+
  theme(text = element_text(size = 30))
ggsave('avp.oxt/Avp.Oxt.pca.PC2.PC3.expression.nolabel.png',
       width = 10,
       height = 10)

#heatmap PCA
DimHeatmap(burtoni.snseq.combined.sct.all.avp.oxt.recluster, 
           dims = 1:4, 
           cells = 500, 
           balanced = TRUE)
ggsave('avp.oxt/Avp.Oxt.pca.heatmap.png',
       width = 10,
       height = 10)
#loadings
VizDimLoadings(burtoni.snseq.combined.sct.all.avp.oxt.recluster, 
               dims = 1:4, 
               reduction = "pca")
ggsave('avp.oxt/Avp.Oxt.pca.loadings.png',
       width = 10,
       height = 10)

#calculate variance
pca <- burtoni.snseq.combined.sct.all.avp.oxt.recluster[["pca"]]
total_variance <- slot(burtoni.snseq.combined.sct.all.avp.oxt.recluster[["pca"]], "misc")$total.variance
eigValues = (pca@stdev)^2  ## EigenValues
varExplained = eigValues / total_variance



## marker genes clusters
burtoni.snseq.combined.sct.all.avp.oxt.recluster.markers <- FindAllMarkers(burtoni.snseq.combined.sct.all.avp.oxt.recluster, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

burtoni.snseq.combined.sct.all.avp.oxt.recluster.markers.top10 = burtoni.snseq.combined.sct.all.avp.oxt.recluster.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) %>% 
  arrange(avg_log2FC,
          .by_group = TRUE)

DoHeatmap(burtoni.snseq.combined.sct.all.avp.oxt.recluster, features = burtoni.snseq.combined.sct.all.avp.oxt.recluster.markers.top10$gene) + NoLegend()
ggsave('avp.oxt/Avp.Oxt.marker.heatmap.top10.png',
       width = 10,
       height = 10)


DoHeatmap(burtoni.snseq.combined.sct.all.avp.oxt.recluster, features = burtoni.snseq.combined.sct.all.avp.oxt.recluster.markers %>%
            group_by(cluster) %>%
            arrange(avg_log2FC,
                    .by_group = TRUE) %>% 
            pull(gene)) + NoLegend()
ggsave('avp.oxt/Avp.Oxt.marker.heatmap.png',
       width = 10,
       height = 10)


# #create dendrogram
# BuildClusterTree(burtoni.snseq.combined.sct.all.avp.oxt.recluster) %>% PlotClusterTree()
# ggsave('avp.oxt/Avp.Oxt.dendrogram.png',
#        width = 10,
#        height = 10)


##create cell type dendrogram
#remove rows with all zeros 
burtoni.snseq.combined.sct.all.avp.oxt.recluster.data <- as.data.frame(log1p(AverageExpression(burtoni.snseq.combined.sct.all.avp.oxt.recluster, verbose = FALSE)$RNA))
burtoni.snseq.combined.sct.all.avp.oxt.recluster.data$gene <- rownames(burtoni.snseq.combined.sct.all.avp.oxt.recluster.data)
#create dendrogram
burtoni.snseq.combined.sct.all.avp.oxt.recluster.tree = pvclust(burtoni.snseq.combined.sct.all.avp.oxt.recluster.data %>% 
                                                                  select(-c(gene)) %>%                                                            mutate(sum = rowSums(.)) %>% 
                                                                  filter(sum > 0) %>% 
                                                                  select(-c(sum)))

#graph
png(file = 'avp.oxt/Avp.Oxt.dendrogram.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 480)
plot(burtoni.snseq.combined.sct.all.avp.oxt.recluster.tree)
dev.off()

## gene expression between nodes
burtoni.snseq.combined.sct.all.avp.oxt.recluster.tree.domvssub = burtoni.snseq.combined.sct.all.avp.oxt.recluster

Idents(burtoni.snseq.combined.sct.all.avp.oxt.recluster.tree.domvssub) <- "orig.ident"

burtoni.snseq.combined.sct.all.avp.oxt.recluster.tree.domvssub <- as.data.frame(log1p(AverageExpression(burtoni.snseq.combined.sct.all.avp.oxt.recluster.tree.domvssub, verbose = FALSE)$RNA))
burtoni.snseq.combined.sct.all.avp.oxt.recluster.tree.domvssub$gene <- rownames(burtoni.snseq.combined.sct.all.avp.oxt.recluster.tree.domvssub)

ggplot(burtoni.snseq.combined.sct.all.avp.oxt.recluster.tree.domvssub, 
       aes(sub_burtoni_snseq, 
           dom_burtoni_snseq)) + 
  geom_point() +
  geom_abline(intercept = 0,
              slope = 1) +
  theme_bw() +
  xlim(0,200) +
  ylim(0,200)
# ggsave('avp.oxt/Avp.domvssub.expression.png',
#        width = 10,
#        height = 10)

## compare gene expression
#avp 1
burtoni.snseq.combined.sct.all.avp.subset = subset(burtoni.snseq.combined.sct.all.avp.oxt.recluster,
                                                   idents = c("1"))

Idents(burtoni.snseq.combined.sct.all.avp.subset) <- "orig.ident"

burtoni.snseq.combined.sct.all.avp.subset <- as.data.frame(log1p(AverageExpression(burtoni.snseq.combined.sct.all.avp.subset, verbose = FALSE)$RNA))
burtoni.snseq.combined.sct.all.avp.subset$gene <- rownames(burtoni.snseq.combined.sct.all.avp.subset)

ggplot(burtoni.snseq.combined.sct.all.avp.subset, 
       aes(sub_burtoni_snseq, 
           dom_burtoni_snseq)) + 
  geom_point() +
  geom_abline(intercept = 0,
              slope = 1) +
  theme_bw()
ggsave('avp.oxt/Avp.1.domvssub.expression.png',
       width = 10,
       height = 10)


#avp 3, 4, 5
burtoni.snseq.combined.sct.all.avp.subset.multi = subset(burtoni.snseq.combined.sct.all.avp.oxt.recluster,
                                                         idents = c("3",
                                                                    "4",
                                                                    "5"))

Idents(burtoni.snseq.combined.sct.all.avp.subset.multi) <- "orig.ident"

burtoni.snseq.combined.sct.all.avp.subset.multi <- as.data.frame(log1p(AverageExpression(burtoni.snseq.combined.sct.all.avp.subset.multi, verbose = FALSE)$RNA))
burtoni.snseq.combined.sct.all.avp.subset.multi$gene <- rownames(burtoni.snseq.combined.sct.all.avp.subset.multi)

ggplot(burtoni.snseq.combined.sct.all.avp.subset.multi, 
       aes(sub_burtoni_snseq, 
           dom_burtoni_snseq)) + 
  geom_point() +
  geom_abline(intercept = 0,
              slope = 1) +
  theme_bw() 
ggsave('avp.oxt/Avp.3.4.5.domvssub.expression.png',
       width = 10,
       height = 10)

# LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)


## create volcano plots
#across social status
#AVP 3 
burtoni.snseq.combined.sct.all.avp.oxt.recluster.3 = subset(burtoni.snseq.combined.sct.all.avp.oxt.recluster,
                                                            idents = 3)

burtoni.snseq.combined.sct.all.avp.oxt.recluster.3 = subset(burtoni.snseq.combined.sct.all.avp.oxt.recluster.3,
                                                            subset = avp > 4 )
#change ident
Idents(burtoni.snseq.combined.sct.all.avp.oxt.recluster.3) <- "orig.ident"
#find markers across samples
burtoni.snseq.combined.sct.all.avp.oxt.recluster.3.markers = FindMarkers(object = burtoni.snseq.combined.sct.all.avp.oxt.recluster.3,
                                                                         ident.1 = "dom_burtoni_snseq", 
                                                                         ident.2 = "sub_burtoni_snseq",
                                                                         test.use = "DESeq2",
                                                                         assay = 'RNA') 

#graph volcano plot
burtoni.snseq.combined.sct.all.avp.oxt.recluster.3.markers %>% 
  mutate(Sig = ifelse(p_val < 0.01 & avg_log2FC > 0.5,
                      "significant dom",
                      ifelse(p_val < 0.01 & avg_log2FC < -0.5,
                             "significant sub",
                             "not significant")))%>% 
  rownames_to_column(var = 'gene') %>% 
  mutate(gene.id = ifelse(Sig != "not significant",
                          gene,
                          NA)) %>%
  ggplot(aes(x = avg_log2FC,
             y = -log10(p_val),
             color = Sig)) +
  geom_point(size = 3) +
  theme_bw() +
  xlim(-2.5, 2.5) +
  ggtitle("Volcano plot AVP.3") +
  scale_color_manual(values = c('black',
                                "red",
                                "blue")) +
  geom_vline(xintercept = -0.5,
             linetype="dotted")+
  geom_vline(xintercept = 0.5,
             linetype="dotted") +
  geom_hline(yintercept = -log10(0.01),
             linetype="dotted")+ 
  geom_text_repel(aes(label = gene.id),
                  size = 5) +
  theme(text = element_text(size = 30),)   
ggsave('avp.oxt/Avp.3.domvssub.volcano.png',
       width = 10,
       height = 10)



#oxt 0 
burtoni.snseq.combined.sct.all.avp.oxt.recluster.0 = subset(burtoni.snseq.combined.sct.all.avp.oxt.recluster,
                                                            idents = 0)

burtoni.snseq.combined.sct.all.avp.oxt.recluster.0 = subset(burtoni.snseq.combined.sct.all.avp.oxt.recluster.0,
                                                            subset = oxt > 4 )
#change ident
Idents(burtoni.snseq.combined.sct.all.avp.oxt.recluster.0) <- "orig.ident"
#find markers across samples
burtoni.snseq.combined.sct.all.avp.oxt.recluster.0.markers = FindMarkers(object = burtoni.snseq.combined.sct.all.avp.oxt.recluster.0,
                                                                         ident.1 = "dom_burtoni_snseq", 
                                                                         ident.2 = "sub_burtoni_snseq") 

burtoni.snseq.combined.sct.all.avp.oxt.recluster.0.markers = FindMarkers(object = burtoni.snseq.combined.sct.all.avp.oxt.recluster.0,
                                                                         ident.1 = "dom_burtoni_snseq", 
                                                                         ident.2 = "sub_burtoni_snseq",
                                                                         test.use = "DESeq2",
                                                                         assay = 'RNA') 

#graph volcano plot
burtoni.snseq.combined.sct.all.avp.oxt.recluster.0.markers %>% 
  mutate(Sig = ifelse(p_val < 0.01 & avg_log2FC > 0.5,
                      "significant dom",
                      ifelse(p_val < 0.01 & avg_log2FC < -0.5,
                             "significant sub",
                             "not significant"))) %>% 
  rownames_to_column(var = 'gene') %>% 
  mutate(gene.id = ifelse(Sig != "not significant",
                          gene,
                          NA)) %>% 
  ggplot(aes(x = avg_log2FC,
             y = -log10(p_val),
             color = Sig)) +
  geom_point(size=3) +
  theme_bw() +
  xlim(-2.5, 2.5) +
  ggtitle("Volcano plot oxt.0") +
  scale_color_manual(values = c('black',
                                "red",
                                "blue")) +
  geom_vline(xintercept = -0.5,
             linetype="dotted")+
  geom_vline(xintercept = 0.5,
             linetype="dotted") +
  geom_hline(yintercept = -log10(0.01),
             linetype="dotted") + 
  geom_text_repel(aes(label = gene.id),
                  size = 5) +
  theme(text = element_text(size = 30))    
ggsave('avp.oxt/oxt.0.domvssub.volcano.png',
       width = 10,
       height = 10)



#AVP 1
burtoni.snseq.combined.sct.all.avp.oxt.recluster.1 = subset(burtoni.snseq.combined.sct.all.avp.oxt.recluster,
                                                            idents = 1)

burtoni.snseq.combined.sct.all.avp.oxt.recluster.1 = subset(burtoni.snseq.combined.sct.all.avp.oxt.recluster.1,
                                                            subset = avp > 4 )
#change ident
Idents(burtoni.snseq.combined.sct.all.avp.oxt.recluster.1) <- "orig.ident"
#find markers across samples
burtoni.snseq.combined.sct.all.avp.oxt.recluster.1.markers = FindMarkers(object = burtoni.snseq.combined.sct.all.avp.oxt.recluster.1,
                                                                         ident.1 = "dom_burtoni_snseq", 
                                                                         ident.2 = "sub_burtoni_snseq",
                                                                         test.use = "DESeq2",
                                                                         assay = 'RNA') 

#graph volcano plot
burtoni.snseq.combined.sct.all.avp.oxt.recluster.1.markers %>% 
  mutate(Sig = ifelse(p_val < 0.01 & avg_log2FC > 0.5,
                      "significant dom",
                      ifelse(p_val < 0.01 & avg_log2FC < -0.5,
                             "significant sub",
                             "not significant")))%>% 
  rownames_to_column(var = 'gene') %>% 
  mutate(gene.id = ifelse(Sig != "not significant",
                          gene,
                          NA)) %>%
  ggplot(aes(x = avg_log2FC,
             y = -log10(p_val),
             color = Sig)) +
  geom_point() +
  theme_bw() +
  xlim(-2.5, 2.5) +
  ggtitle("Volcano plot AVP.1") +
  scale_color_manual(values = c('black',
                                "red",
                                "blue")) +
  geom_vline(xintercept = -0.5,
             linetype="dotted")+
  geom_vline(xintercept = 0.5,
             linetype="dotted") +
  geom_hline(yintercept = -log10(0.01),
             linetype="dotted")+ 
  geom_text_repel(aes(label = gene.id)) 
ggsave('avp.oxt/Avp.1.domvssub.volcano.png',
       width = 10,
       height = 10)



#oxt 2 
burtoni.snseq.combined.sct.all.avp.oxt.recluster.2 = subset(burtoni.snseq.combined.sct.all.avp.oxt.recluster,
                                                            idents = 2)

burtoni.snseq.combined.sct.all.avp.oxt.recluster.2 = subset(burtoni.snseq.combined.sct.all.avp.oxt.recluster.2,
                                                            subset = oxt > 4 )
#change ident
Idents(burtoni.snseq.combined.sct.all.avp.oxt.recluster.2) <- "orig.ident"
#find markers across samples
burtoni.snseq.combined.sct.all.avp.oxt.recluster.2.markers = FindMarkers(object = burtoni.snseq.combined.sct.all.avp.oxt.recluster.2,
                                                                         ident.1 = "dom_burtoni_snseq", 
                                                                         ident.2 = "sub_burtoni_snseq") 

burtoni.snseq.combined.sct.all.avp.oxt.recluster.2.markers = FindMarkers(object = burtoni.snseq.combined.sct.all.avp.oxt.recluster.2,
                                                                         ident.1 = "dom_burtoni_snseq", 
                                                                         ident.2 = "sub_burtoni_snseq",
                                                                         test.use = "DESeq2",
                                                                         assay = 'RNA') 

#graph volcano plot
burtoni.snseq.combined.sct.all.avp.oxt.recluster.2.markers %>% 
  mutate(Sig = ifelse(p_val < 0.01 & avg_log2FC > 0.5,
                      "significant dom",
                      ifelse(p_val < 0.01 & avg_log2FC < -0.5,
                             "significant sub",
                             "not significant"))) %>% 
  rownames_to_column(var = 'gene') %>% 
  mutate(gene.id = ifelse(Sig != "not significant",
                          gene,
                          NA)) %>% 
  ggplot(aes(x = avg_log2FC,
             y = -log10(p_val),
             color = Sig)) +
  geom_point() +
  theme_bw() +
  xlim(-0.5, 0.5) +
  ggtitle("Volcano plot oxt.2") +
  scale_color_manual(values = c('black',
                                "red",
                                "blue")) +
  geom_vline(xintercept = -0.5,
             linetype="dotted")+
  geom_vline(xintercept = 0.5,
             linetype="dotted") +
  geom_hline(yintercept = -log10(0.01),
             linetype="dotted") + 
  geom_text_repel(aes(label = gene.id)) 
ggsave('avp.oxt/oxt.2.domvssub.volcano.png',
       width = 10,
       height = 10)


###dotplot
##oxt and avp
DotPlot(burtoni.snseq.combined.sct.all.avp.oxt.recluster, 
        features = c("avp",
                     "oxt"), 
        cols = c("grey", 
                 "red"), 
        dot.scale = 8,
        col.max = 1,
        dot.min = .5) + 
  RotatedAxis() +
  ggtitle("Avp.Oxt")
ggsave('avp.oxt/Avp.Oxt.clustering.dotplot.png',
       width = 10,
       height = 10)

## cell type
DotPlot(burtoni.snseq.combined.sct.all.avp.oxt.recluster, 
        features = burtoni.snseq.combined.sct.all.avp.oxt.recluster@meta.data %>% 
          select(ends_with(".score1")) %>% 
          colnames(), 
        cols = c("grey", 
                 "red"), 
        dot.scale = 8,
        col.min = 1,
        dot.min = .4) + 
  RotatedAxis()
ggsave('avp.oxt/Avp.Oxt.celltype.dotplot.png',
       width = 10,
       height = 10)

##markers parvo vs magno
DotPlot(burtoni.snseq.combined.sct.all.avp.oxt.recluster, 
        features = c('KCNMB4',
                     "RELN",
                     "cnr1",
                     "CNR1"), 
        cols = c("grey", 
                 "red"),
        col.min =  1,
        dot.min = .1) + 
  RotatedAxis() +
  ggtitle("Avp.Oxt")
ggsave('avp.oxt/Avp.Oxt.magno.parvo.markers.dotplot.png',
       width = 10,
       height = 10)

### clustree
# cluster across resolutions
burtoni.snseq.combined.sct.all.avp.oxt.recluster.clustree <- Seurat::FindClusters(object = burtoni.snseq.combined.sct.all.avp.oxt.recluster, 
                                                                                  resolution = resolution.range.reduced)
#check data
head(burtoni.snseq.combined.sct.all.avp.oxt.recluster.clustree[[]])
#set presentation colors
presentation.color <- c('#66c2a5',
                        '#fc8d62',
                        '#8da0cb',
                        '#e78ac3',
                        '#a6d854',
                        '#ffd92f',
                        '#e5c494',
                        '#b3b3b3')
#clustree
clustree(burtoni.snseq.combined.sct.all.avp.oxt.recluster.clustree, 
         prefix = "integrated_snn_res.",
         node_colour = 'cluster',
         node_size_range = c(10,20),
         scale_node_text = TRUE) +
  scale_edge_color_continuous(low = "black", 
                              high = "black") +
  scale_color_manual(values = presentation.color)+
  theme(legend.position = "bottom")
ggsave('avp.oxt/avp.oxt.clustree.png',
       width = 10,
       height = 10)

## loop through resolutions for presentation
for (x in resolution.range.reduced) {
  DimPlot(burtoni.snseq.combined.sct.all.avp.oxt.recluster.clustree, 
          reduction = "umap",
          group.by = paste('integrated_snn_res',
                           x,
                           sep = '.'),
          label = T,
          pt.size = 2,
          label.size = 10,
          repel = T) +
    theme(legend.position = 'none') +
    scale_color_manual(values = presentation.color)
  ggsave(paste('avp.oxt/resolution/Avp.oxt.resolution',
               x,
               'dimplot.avp.oxt.all.png',
               sep = '.'),
         width = 10,
         height = 10)
}

### identify marker genes
## all
burtoni.snseq.combined.sct.all.avp.oxt.recluster.clustree.markers <- FindAllMarkers(burtoni.snseq.combined.sct.all.avp.oxt.recluster.clustree, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

##compare clusters
#AVP
#1 vs 3
burtoni.snseq.combined.sct.all.avp.oxt.recluster.clustree.markers.1v3 <- FindMarkers(burtoni.snseq.combined.sct.all.avp.oxt.recluster.clustree, 
                                                                                     ident.1 = 1, 
                                                                                     ident.2 = 3, 
                                                                                     min.pct = 0.5, 
                                                                                     logfc.threshold = 1)

#OXT
#2 vs 0
burtoni.snseq.combined.sct.all.avp.oxt.recluster.clustree.markers.2v0 <- FindMarkers(burtoni.snseq.combined.sct.all.avp.oxt.recluster.clustree, 
                                                                                     ident.1 = 2, 
                                                                                     ident.2 = 0, 
                                                                                     min.pct = 0.55, 
                                                                                     logfc.threshold = 1)

## compare concordance across groups
burtoni.snseq.combined.sct.all.avp.oxt.recluster.clustree.markers.compare = full_join(burtoni.snseq.combined.sct.all.avp.oxt.recluster.clustree.markers.2v0 %>% 
                                                                                        rownames_to_column("gene"),
                                                                                      burtoni.snseq.combined.sct.all.avp.oxt.recluster.clustree.markers.1v3 %>% 
                                                                                        rownames_to_column("gene"),
                                                                                      by = "gene",
                                                                                      suffix = c(".2v0",
                                                                                                 ".1v3"))

## graph concordance
burtoni.snseq.combined.sct.all.avp.oxt.recluster.clustree.markers.compare %>% 
  ggplot(aes(x = avg_log2FC.2v0,
             y = avg_log2FC.1v3)) +
  geom_point() +
  theme_classic() +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0)


## 2 v 0 
burtoni.snseq.combined.sct.all.avp.oxt.recluster.clustree.markers.compare %>% 
  mutate(empty.1v3 = ifelse(is.na(avg_log2FC.1v3),
                            "NA",
                            "Overlap")) %>% 
  ggplot(aes(x = avg_log2FC.2v0,
             y = -log(p_val_adj.2v0),
             color = empty.1v3)) +
  geom_point() +
  theme_classic()

## 1 v 3 
burtoni.snseq.combined.sct.all.avp.oxt.recluster.clustree.markers.compare %>% 
  mutate(empty.2v0 = ifelse(is.na(avg_log2FC.2v0),
                            "NA",
                            "Overlap")) %>% 
  ggplot(aes(x = avg_log2FC.1v3,
             y = -log(p_val_adj.1v3),
             color = empty.2v0)) +
  geom_point() +
  theme_classic()

# remove na
burtoni.snseq.combined.sct.all.avp.oxt.recluster.clustree.markers.compare %>% 
  na.omit() %>% 
  View()

# remove pval adjust
# remove na
burtoni.snseq.combined.sct.all.avp.oxt.recluster.clustree.markers.compare %>% 
  na.omit() %>% 
  filter(p_val_adj.1v3<0.05) %>% 
  filter(p_val_adj.2v0<0.05) %>% 
  View()

### combine with scsorter cell type
## subset highly variable genes
## identify top variable genes
burtoni.snseq.combined.sct.all.avp.oxt.recluster.filter <- FindVariableFeatures(burtoni.snseq.combined.sct.all.avp.oxt.recluster, 
                                              selection.method = "vst", 
                                              nfeatures = 2000, 
                                              verbose = F)

# identify top 2000 variable genes
topgenes <- head(VariableFeatures(burtoni.snseq.combined.sct.all.avp.oxt.recluster.filter), 
                 2000)

#get assay data
#as matrix
burtoni.snseq.combined.sct.all.avp.oxt.recluster.filter.extract = GetAssayData(burtoni.snseq.combined.sct.all.avp.oxt.recluster.filter) %>% 
  as.matrix()

# filter genes
topgene_filter = rowSums(as.matrix(burtoni.snseq.combined.sct.all.avp.oxt.recluster.filter.extract)[topgenes, ]!=0) > ncol(burtoni.snseq.combined.sct.all.avp.oxt.recluster.filter.extract)*.1
#create filter gene list
topgenes = topgenes[topgene_filter]

# subset data by gene list
burtoni.snseq.combined.sct.all.avp.oxt.recluster.filter.extract = burtoni.snseq.combined.sct.all.avp.oxt.recluster.filter.extract[rownames(burtoni.snseq.combined.sct.all.avp.oxt.recluster.filter.extract) %in% topgenes, ]

#subset annotation file
burtoni.scsorter.data.scsort.output.avp.oxt = inner_join(burtoni.scsorter.data.scsort.output,
burtoni.snseq.combined.sct.all.avp.oxt.recluster.filter.extract %>% 
  colnames() %>% 
  as.data.frame() %>% 
  rename('Cell.id' = '.'))

##pheatmap
# create heatmap using pheatmap
pheatmap(burtoni.snseq.combined.sct.all.avp.oxt.recluster.filter.extract %>% 
           t(),
         annotation_row = burtoni.scsorter.data.scsort.output.avp.oxt %>% 
           select(Cell.type,
                  Cell.id) %>% 
           column_to_rownames("Cell.id"),
         scale = "column") +
  ggtitle('heatmap top 2000 variable genes')
ggsave('avp.oxt/heatmap cell types avp oxt top 2000 variable genes.png',
       width = 10,
       height = 10)

## create UMAP plot 
#extract umap data
## combine umap and metadata
burtoni.snseq.combined.sct.all.avp.oxt.recluster.umap = full_join(burtoni.snseq.combined.sct.all.avp.oxt.recluster@reductions$umap@cell.embeddings %>% 
                                                                    as.data.frame() %>% 
                                                                    rownames_to_column("Cell.id"),
                                                                  burtoni.snseq.combined.sct.all.avp.oxt.recluster@meta.data %>% 
                                                                    rownames_to_column("Cell.id"))

#graph umap
burtoni.snseq.combined.sct.all.avp.oxt.recluster.umap %>% 
  ggplot(aes(x=UMAP_1,
             y= UMAP_2,
             color = Cell.type)) +
  geom_point() +
  theme_bw()
ggsave('avp.oxt/scsorter/UMAP scsorter cell types avp oxt.png',
       width = 10,
       height = 10)

#check cell count
table(burtoni.snseq.combined.sct.all.avp.oxt.recluster.umap$Cell.type)

#check by orig.ident
table(burtoni.snseq.combined.sct.all.avp.oxt.recluster.umap %>% 
        select(Cell.type,
               orig.ident))


#neurons vs non-neurons
burtoni.snseq.combined.sct.all.avp.oxt.recluster.umap %>% 
  mutate(Cell.type.neuron = ifelse(Cell.type %in% c('neurons',
                                                     'excitatory',
                                                     'inhibitory'),
                                    "neurons",
                                    "other")) %>% 
  ggplot(aes(x=UMAP_1,
             y= UMAP_2,
             color = Cell.type.neuron,
             shape = integrated_snn_res.1)) +
  geom_point() +
  theme_bw()
ggsave('avp.oxt/scsorter/UMAP scsorter cell types neurons vs other avp oxt.png',
       width = 10,
       height = 10)





##expression level
burtoni.snseq.combined.sct.all.avp.oxt.recluster.expression = full_join(burtoni.snseq.combined.sct.all.avp.oxt.recluster.umap,
                                                                        burtoni.snseq.combined.sct.all.avp.oxt.recluster@assays$integrated@scale.data %>% 
                                                                          as.data.frame() %>% 
                                                                          filter(rownames(burtoni.snseq.combined.sct.all.avp.oxt.recluster@assays$integrated@scale.data) %in% c('avp',
                                                                                                                                                                                'oxt')) %>% 
                                                                          t() %>% as.data.frame() %>% 
                                                                          rownames_to_column('Cell.id'))

#graph
#oxt
burtoni.snseq.combined.sct.all.avp.oxt.recluster.expression %>%
  filter(oxt > 4) %>% 
  ggplot(aes(x=Cell.type ,
             y=avp)) + 
  geom_violin() +
  theme_bw()
ggsave('avp.oxt/scsorter/oxt scsorter cell types.png',
       width = 10,
       height = 10)

#avp
burtoni.snseq.combined.sct.all.avp.oxt.recluster.expression %>%
  filter(avp > 4) %>% 
  ggplot(aes(x=Cell.type ,
             y=avp)) + 
  geom_violin() +
  theme_bw()
ggsave('avp.oxt/scsorter/avp scsorter cell types.png',
       width = 10,
       height = 10)

#cell counts of oxt avp 
table(burtoni.snseq.combined.sct.all.avp.oxt.recluster.expression %>% 
        mutate(Oxt.cell = ifelse(oxt > 4 & avp > 4,
               'both',
               ifelse(oxt > 4,
                      'oxt',
                      'avp'))) %>% 
        select(Cell.type,
               Oxt.cell))

#cell counts of oxt avp cell neurons
table(burtoni.snseq.combined.sct.all.avp.oxt.recluster.expression %>% 
        filter(Cell.type %in% c("neurons", 
                                "excitatory",
                                "inhibitory")) %>% 
        mutate(Oxt.cell = ifelse(oxt > 4 & avp > 4,
                                 'both',
                                 ifelse(oxt > 4,
                                        'oxt',
                                        'avp'))) %>% 
        select(Cell.type,
               Oxt.cell))



#### AVP  ####
### avp 
### subset data to relevant cells
burtoni.snseq.combined.sct.all.avp = burtoni.snseq.combined.sct.all

#subset
burtoni.snseq.combined.sct.all.avp = subset(burtoni.snseq.combined.sct.all.avp,
                                                    subset = avp > 4)

## run PCA, UMAP, and cluster 
#use 0.4 resolution
burtoni.snseq.combined.sct.all.avp.recluster = burtoni.snseq.combined.sct.all.avp %>% 
  RunPCA() %>%
  FindNeighbors(dims = 1:15) %>%
  RunUMAP(dims = 1:15) %>%
  FindClusters(resolution = 0.6)

### clustree
# cluster across resolutions
burtoni.snseq.combined.sct.all.avp.clustree <- Seurat::FindClusters(object = burtoni.snseq.combined.sct.all.avp.recluster, 
                                                                            resolution = resolution.range.reduced)
#check data
head(burtoni.snseq.combined.sct.all.avp.clustree[[]])
#set presentation colors
presentation.color <- c('#66c2a5',
                        '#fc8d62',
                        '#8da0cb',
                        '#e78ac3',
                        '#a6d854',
                        '#ffd92f',
                        '#e5c494',
                        '#b3b3b3')
#clustree
clustree(burtoni.snseq.combined.sct.all.avp.clustree, 
         prefix = "integrated_snn_res.",
         node_colour = 'cluster',
         node_size_range = c(10,20),
         scale_node_text = TRUE) +
  scale_edge_color_continuous(low = "black", 
                              high = "black") +
  scale_color_manual(values = presentation.color)+
  theme(legend.position = "bottom")
ggsave('avp.oxt/avp/avp.clustree.png',
       width = 10,
       height = 10)


### graph 
##dom vs sub
DimPlot(burtoni.snseq.combined.sct.all.avp.recluster, 
        reduction = "umap", 
        group.by = "orig.ident")
ggsave('avp.oxt/avp/DomVsSub.dimplot.avp.all.png')

## clusters
DimPlot(burtoni.snseq.combined.sct.all.avp.recluster, 
        reduction = "umap", 
        label = TRUE,
        repel = TRUE)
ggsave('avp.oxt/avp/Clusters.dimplot.avp.all.png',
       width = 10,
       height = 10)

#### AVP neurons 4 ####
### avp neurons
### subset data to relevant cells
burtoni.snseq.combined.sct.all.avp.neurons.4 = burtoni.snseq.combined.sct.all

#set idents
Idents(object = burtoni.snseq.combined.sct.all.avp.neurons.4) <- "Cell.type"

#subset
burtoni.snseq.combined.sct.all.avp.neurons.4 = subset(burtoni.snseq.combined.sct.all.avp.neurons.4,
                                                    subset = avp > 4, idents = c("neurons",
                                                                                 "excitatory",
                                                                                 "inhibitory"))


## run PCA, UMAP, and cluster 
#use 0.4 resolution
burtoni.snseq.combined.sct.all.avp.neurons.4.recluster = burtoni.snseq.combined.sct.all.avp.neurons.4 %>% 
  RunPCA() %>%
  FindNeighbors(dims = 1:15) %>%
  RunUMAP(dims = 1:15) %>%
  FindClusters(resolution = 1)

### clustree
# cluster across resolutions
burtoni.snseq.combined.sct.all.avp.neurons.4.clustree <- Seurat::FindClusters(object = burtoni.snseq.combined.sct.all.avp.neurons.4.recluster, 
                                                                            resolution = resolution.range.reduced)
#clustree
clustree(burtoni.snseq.combined.sct.all.avp.neurons.4.clustree, 
         prefix = "integrated_snn_res.",
         node_colour = 'cluster',
         node_size_range = c(10,20),
         scale_node_text = TRUE) +
  scale_edge_color_continuous(low = "black", 
                              high = "black")
ggsave('avp.oxt/avp.4.threshold/avp.neurons.clustree.png',
       width = 10,
       height = 10)


### graph 
##dom vs sub
DimPlot(burtoni.snseq.combined.sct.all.avp.neurons.4.recluster, 
        reduction = "umap", 
        group.by = "orig.ident")
ggsave('avp.oxt/avp.4.threshold/DomVsSub.dimplot.avp.neurons.all.png')

## clusters
DimPlot(burtoni.snseq.combined.sct.all.avp.neurons.4.recluster, 
        reduction = "umap", 
        label = TRUE,
        repel = TRUE)
ggsave('avp.oxt/avp.4.threshold/Clusters.dimplot.avp.neurons.all.png',
       width = 10,
       height = 10)


## clusters
# avp expression
FeaturePlot(burtoni.snseq.combined.sct.all.avp.neurons.4.recluster, 
            reduction = "umap",
            features = c('avp'),
            min.cutoff = 2,
            max.cutoff = 5,
            label = TRUE,
            repel = TRUE)
ggsave('avp.oxt/avp.4.threshold/Clusters.dimplot.avp.neurons.all.expression.png',
       width = 10,
       height = 10)

##expression level
burtoni.snseq.combined.sct.all.avp.neurons.4.recluster.expression = full_join(full_join(burtoni.snseq.combined.sct.all.avp.neurons.4.recluster@reductions$umap@cell.embeddings %>% 
                                                                                        as.data.frame() %>% 
                                                                                        rownames_to_column("Cell.id"),
                                                                                        burtoni.snseq.combined.sct.all.avp.neurons.4.recluster@meta.data %>% 
                                                                                        rownames_to_column("Cell.id")),
                                                                              burtoni.snseq.combined.sct.all.avp.neurons.4.recluster@assays$integrated@scale.data %>% 
                                                                              as.data.frame() %>% 
                                                                              filter(rownames(burtoni.snseq.combined.sct.all.avp.neurons.4.recluster@assays$integrated@scale.data) %in% c('avp',
                                                                                                                                                                                        'oxt')) %>% 
                                                                              t() %>% as.data.frame() %>% 
                                                                              rownames_to_column('Cell.id')) %>% 
  full_join(burtoni.snseq.combined.sct.all.avp.neurons.4.recluster@reductions$pca@cell.embeddings %>% 
              as.data.frame() %>% 
              select(PC_1,
                     PC_2,
                     PC_3,
                     PC_4,
                     PC_5) %>% 
              rownames_to_column("Cell.id")) %>% 
  mutate(Oxt.cell = ifelse(oxt > 4 & avp > 4,
                           'both',
                           ifelse(oxt > 4,
                                  'oxt',
                                  'avp')))

# seperation by cell type
table(burtoni.snseq.combined.sct.all.avp.neurons.4.recluster.expression %>% 
        filter(Oxt.cell != 'both') %>% 
        select(Cell.type,
               integrated_snn_res.1))

# difference social status
table(burtoni.snseq.combined.sct.all.avp.neurons.4.recluster.expression %>% 
        filter(Oxt.cell != 'both') %>% 
        select(integrated_snn_res.1,
               orig.ident))


#### AVP neurons ####
### avp neurons
### subset data to relevant cells
burtoni.snseq.combined.sct.all.avp.neurons = burtoni.snseq.combined.sct.all

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

#subset
burtoni.snseq.combined.sct.all.avp.neurons = subset(burtoni.snseq.combined.sct.all.avp.neurons,
                                             subset = avp >= 2, idents = c("neurons",
                                                                             "excitatory",
                                                                             "inhibitory"))


## run PCA, UMAP, and cluster 
#use 0.8 resolution
burtoni.snseq.combined.sct.all.avp.neurons.recluster = burtoni.snseq.combined.sct.all.avp.neurons %>% 
  RunPCA() %>%
  FindNeighbors(dims = 1:15) %>%
  RunUMAP(dims = 1:15) %>%
  FindClusters(resolution = 0.8)

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
#         reduction = "umap",
#         features = c('avp'),
#         min.cutoff = 2,
#         max.cutoff = 5,
#         label = TRUE,
#         repel = TRUE)
# ggsave('avp.oxt/avp/Clusters.dimplot.avp.neurons.all.expression.png',
#        width = 10,
#        height = 10)

##expression level
burtoni.snseq.combined.sct.all.avp.neurons.recluster.expression = full_join(full_join(burtoni.snseq.combined.sct.all.avp.neurons.recluster@reductions$umap@cell.embeddings %>% 
                                                                                            as.data.frame() %>% 
                                                                                            rownames_to_column("Cell.id"),
                                                                                          burtoni.snseq.combined.sct.all.avp.neurons.recluster@meta.data %>% 
                                                                                            rownames_to_column("Cell.id")),
                                                                                burtoni.snseq.combined.sct.all.avp.neurons.recluster@assays$integrated@scale.data %>% 
                                                                                  as.data.frame() %>% 
                                                                                  filter(rownames(burtoni.snseq.combined.sct.all.avp.neurons.recluster@assays$integrated@scale.data) %in% c('avp',
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
  mutate(Oxt.cell = ifelse(oxt > 2 & avp > 2,
                           'both',
                           ifelse(oxt > 2,
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










#### neuropeptides neurons ####
### avp oxt neurons
### subset data to relevant cells
burtoni.snseq.combined.sct.all.avp.oxt.neurons = burtoni.snseq.combined.sct.all

#set idents
Idents(object = burtoni.snseq.combined.sct.all.avp.oxt.neurons) <- "Cell.type"

#subset neurons
burtoni.snseq.combined.sct.all.avp.oxt.neurons = subset(burtoni.snseq.combined.sct.all.avp.oxt.neurons,
                                                        idents = c("neurons",
                                                                   "excitatory",
                                                                   "inhibitory"))

#subset expression
burtoni.snseq.combined.sct.all.avp.oxt.neurons = subset(burtoni.snseq.combined.sct.all.avp.oxt.neurons,
                                                    subset = avp > 4 | oxt > 4)


## run PCA, UMAP, and cluster 
#use 0.4 resolution
burtoni.snseq.combined.sct.all.avp.oxt.neurons.recluster = burtoni.snseq.combined.sct.all.avp.oxt.neurons %>% 
  RunPCA() %>%
  FindNeighbors(dims = 1:15) %>%
  RunUMAP(dims = 1:15) %>%
  FindClusters(resolution = 0.8)

### clustree
# cluster across resolutions
burtoni.snseq.combined.sct.all.avp.oxt.neurons.clustree <- Seurat::FindClusters(object = burtoni.snseq.combined.sct.all.avp.oxt.neurons.recluster, 
                                                                            resolution = resolution.range)
#check data
head(burtoni.snseq.combined.sct.all.avp.oxt.neurons.clustree[[]])
#set presentation colors
presentation.color <- c('#66c2a5',
                        '#fc8d62',
                        '#8da0cb',
                        '#e78ac3',
                        '#a6d854',
                        '#ffd92f',
                        '#e5c494',
                        '#b3b3b3')
#clustree
clustree(burtoni.snseq.combined.sct.all.avp.oxt.neurons.clustree, 
         prefix = "integrated_snn_res.",
         node_colour = 'cluster',
         node_size_range = c(10,20),
         scale_node_text = TRUE) +
  scale_edge_color_continuous(low = "black", 
                              high = "black") +
  scale_color_manual(values = presentation.color)+
  theme(legend.position = "bottom")
ggsave('avp.oxt/avp/avp.oxt.neurons.clustree.png',
       width = 10,
       height = 10)


### graph 
##dom vs sub
DimPlot(burtoni.snseq.combined.sct.all.avp.oxt.neurons.recluster, 
        reduction = "umap", 
        group.by = "orig.ident")
ggsave('avp.oxt/avp/DomVsSub.dimplot.avp.oxt.neurons.all.png')

## clusters
DimPlot(burtoni.snseq.combined.sct.all.avp.oxt.neurons.recluster, 
        reduction = "umap", 
        label = TRUE,
        repel = TRUE)
ggsave('avp.oxt/avp/Clusters.dimplot.avp.oxt.neurons.all.png',
       width = 10,
       height = 10)

##expression level
burtoni.snseq.combined.sct.all.avp.oxt.neurons.recluster.expression = full_join(full_join(burtoni.snseq.combined.sct.all.avp.oxt.neurons.recluster@reductions$umap@cell.embeddings %>% 
                                                                                            as.data.frame() %>% 
                                                                                            rownames_to_column("Cell.id"),
                                                                                          burtoni.snseq.combined.sct.all.avp.oxt.neurons.recluster@meta.data %>% 
                                                                                            rownames_to_column("Cell.id")),
                                                                                burtoni.snseq.combined.sct.all.avp.oxt.neurons.recluster@assays$integrated@scale.data %>% 
                                                                          as.data.frame() %>% 
                                                                          filter(rownames(burtoni.snseq.combined.sct.all.avp.oxt.neurons.recluster@assays$integrated@scale.data) %in% c('avp',
                                                                                                                                                                                'oxt')) %>% 
                                                                          t() %>% as.data.frame() %>% 
                                                                          rownames_to_column('Cell.id')) %>% 
  full_join(burtoni.snseq.combined.sct.all.avp.oxt.neurons.recluster@reductions$pca@cell.embeddings %>% 
              as.data.frame() %>% 
              select(PC_1,
                     PC_2,
                     PC_3,
                     PC_4,
                     PC_5) %>% 
              rownames_to_column("Cell.id")) %>% 
  mutate(Oxt.cell = ifelse(oxt > 4 & avp > 4,
                           'both',
                           ifelse(oxt > 4,
                                  'oxt',
                                  'avp')))


##graph
burtoni.snseq.combined.sct.all.avp.oxt.neurons.recluster.expression %>% 
  filter(Oxt.cell != 'both') %>% 
  ggplot(aes(x = UMAP_1,
             y = UMAP_2,
             color = integrated_snn_res.0.8,
             shape = Oxt.cell)) +
  geom_point() +
  theme_bw()
ggsave('avp.oxt/avp/Umap.avp.oxt.neurons.all.png',
       width = 10,
       height = 10)

burtoni.snseq.combined.sct.all.avp.oxt.neurons.recluster.expression %>% 
  filter(Oxt.cell != 'both') %>% 
  ggplot(aes(x = UMAP_1,
             y = UMAP_2,
             shape = integrated_snn_res.0.8,
             color = Oxt.cell)) +
  geom_point() +
  theme_bw()
ggsave('avp.oxt/avp/Umap.expression.avp.oxt.neurons.all.png',
       width = 10,
       height = 10)

## cell count per cluster
table(burtoni.snseq.combined.sct.all.avp.oxt.neurons.recluster.expression %>% 
        filter(Oxt.cell != 'both') %>% 
        select(integrated_snn_res.0.8,
               Oxt.cell))

### PCA
#1 vs 2
burtoni.snseq.combined.sct.all.avp.oxt.neurons.recluster.expression %>% 
  filter(Oxt.cell != 'both') %>% 
  ggplot(aes(x = PC_1,
             y = PC_2,
             color = integrated_snn_res.0.8,
             shape = Oxt.cell)) +
  geom_point() +
  theme_bw()
ggsave('avp.oxt/avp/PCA.1.2.avp.oxt.neurons.all.png',
       width = 10,
       height = 10)
#2 vs 3
burtoni.snseq.combined.sct.all.avp.oxt.neurons.recluster.expression %>% 
  filter(Oxt.cell != 'both') %>% 
  ggplot(aes(x = PC_3,
             y = PC_2,
             shape = integrated_snn_res.0.8,
             color = Oxt.cell)) +
  geom_point() +
  theme_bw()
ggsave('avp.oxt/avp/PCA.2.3.expression.avp.oxt.neurons.all.png',
       width = 10,
       height = 10)

burtoni.snseq.combined.sct.all.avp.oxt.neurons.recluster.expression %>% 
  filter(Oxt.cell != 'both') %>% 
  ggplot(aes(x = PC_3,
             y = PC_2,
             color = integrated_snn_res.0.8,
             shape = Oxt.cell)) +
  geom_point() +
  theme_bw()
ggsave('avp.oxt/avp/PCA.2.3.avp.oxt.neurons.all.png',
       width = 10,
       height = 10)


#### AVP neurons differential expression ####
#### Use DEsingle to calculate DEG between cell profile groups

### calculate variable genes
## identify top 500 variable genes
burtoni.snseq.combined.sct.all.avp.neurons.recluster.variable <- FindVariableFeatures(burtoni.snseq.combined.sct.all.avp.neurons.recluster, 
                                              selection.method = "vst", 
                                              nfeatures = 500, 
                                              verbose = F)
# identify top 500 variable genes
avp.neuron.group.topgenes <- head(VariableFeatures(burtoni.snseq.combined.sct.all.avp.neurons.recluster.variable), 
                 500)


### group vector
## vector of factor specifies groups
## corresponds to columns in counts
## subset to groups of interest 

avp.neuron.group.vector.full = burtoni.snseq.combined.sct.all.avp.neurons.recluster.expression %>% 
  filter(integrated_snn_res.0.8 != 2)

## create vector of factor
avp.neuron.group.vector.list = avp.neuron.group.vector.full %>% 
  mutate(avp.neuron.group = integrated_snn_res.1 %>% 
           as.factor()) %>% 
  pull(avp.neuron.group) %>% 
  droplevels()

### counts matrix
## raw read count matrix
## rows = genes, columns = cells
# subset to groups of interest (0 and 1 of AVP neurons)
# only keep 500 variable genes
# set negative values to 0
avp.neuron.group.vector.count = GetAssayData(burtoni.snseq.combined.sct.all.avp.neurons.recluster) %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column('gene') %>% 
  select(c(gene,
           avp.neuron.group.vector.full %>% 
             pull(Cell.id))) %>% 
  filter(gene %in% avp.neuron.group.topgenes) %>% 
  column_to_rownames('gene') %>% 
  as.matrix() %>% 
  pmax(0)

### Detecting the DE genes
avp.neuron.group.results = DEsingle(counts = avp.neuron.group.vector.count, 
                    group = avp.neuron.group.vector.list,
                    parallel = TRUE)

# Dividing the DE genes into 3 categories at threshold of FDR < 0.1
avp.neuron.group.results.classified <- DEtype(results = avp.neuron.group.results, 
                             threshold = 0.1)


# create color groups
avp.neuron.group.results.classified = avp.neuron.group.results.classified %>% 
  mutate(Sig = ifelse(pvalue <= 0.01,
                      'Sig',
                      'Not Sig'),
         Direction = total_mean_1 - total_mean_2,
         Direction.type = ifelse(Direction > 0,
                                 'up',
                                 'down'),
         Sig.direction = ifelse(Sig == 'Sig',
                                Direction.type,
                                Sig))


##graph volcano plot
# create log fold change and -log pvalue
avp.neuron.group.results.classified %>% 
  rownames_to_column('gene') %>% 
  ggplot(aes(x = log(foldChange),
             y = -log10(pvalue),
             color = Sig.direction)) +
  geom_point(size = 5) +
  theme_classic() +
  xlim(c(-4.5,4.5)) + 
  scale_color_manual(values=c("21B9CA", 
                              "grey", 
                              "orange3")) +
  geom_text(data = avp.neuron.group.results.classified %>% 
              rownames_to_column('gene') %>% 
              filter(Sig == 'Sig'),
            aes(label = gene),
            vjust = 0, 
            nudge_y = 0.15)
ggsave('avp.oxt/avp/avp.neurons.groups.1vs0.volcano.pdf',
       width = 10,
       height = 10)

#TCF4, mpped1, pcdh1b, gh1, prl, COX3, ND4

### dom vs sub
### cluster 0 
### group vector
## vector of factor specifies groups
## corresponds to columns in counts
## subset to groups of interest 
avp.neuron.group.vector.0 = burtoni.snseq.combined.sct.all.avp.neurons.recluster.expression %>% 
  filter(integrated_snn_res.0.8 == 0)

## create vector of factor
avp.neuron.0.vector.list = avp.neuron.group.vector.0 %>% 
  mutate(avp.neuron.0.orig.ident = orig.ident %>% 
           as.factor()) %>% 
  pull(avp.neuron.0.orig.ident) %>% 
  droplevels()

### counts matrix
## raw read count matrix
## rows = genes, columns = cells
# subset to groups of interest (0 and 1 of AVP neurons)
# only keep 500 variable genes
# set negative values to 0
avp.neuron.0.vector.count = GetAssayData(burtoni.snseq.combined.sct.all.avp.neurons.recluster) %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column('gene') %>% 
  select(c(gene,
           avp.neuron.group.vector.0 %>% 
             pull(Cell.id))) %>% 
  filter(gene %in% avp.neuron.group.topgenes) %>% 
  column_to_rownames('gene') %>% 
  as.matrix() %>% 
  pmax(0)

### Detecting the DE genes
avp.neuron.0.results = DEsingle(counts = avp.neuron.0.vector.count, 
                                    group = avp.neuron.0.vector.list,
                                    parallel = TRUE)

# Dividing the DE genes into 3 categories at threshold of FDR < 0.1
avp.neuron.0.results.classified <- DEtype(results = avp.neuron.0.results, 
                                              threshold = 0.1)


# create color groups
avp.neuron.0.results.classified = avp.neuron.0.results.classified %>% 
  mutate(Sig = ifelse(pvalue <= 0.01,
                      'Sig',
                      'Not Sig'),
         Direction = total_mean_1 - total_mean_2,
         Direction.type = ifelse(Direction > 0,
                                 'up',
                                 'down'),
         Sig.direction = ifelse(Sig == 'Sig',
                                Direction.type,
                                Sig))


##graph volcano plot
# create log fold change and -log pvalue
avp.neuron.0.results.classified %>% 
  rownames_to_column('gene') %>% 
  mutate(logFoldChange = log(foldChange)) %>% 
  mutate(logFoldChange = ifelse(logFoldChange < -3.5,
                                -3.3,
                                logFoldChange),
         logFoldChange = ifelse(logFoldChange > 3,
                                3.3,
                                logFoldChange)) %>% 
  mutate(sig.label = ifelse(Sig == 'Sig',
                            gene,
                            '')) %>% 
  ggplot(aes(x = logFoldChange,
             y = -log10(pvalue),
             color = Sig.direction))+
  geom_vline(xintercept = 3.2, 
             linetype="dotted")+
  geom_vline(xintercept = -3.2,
             linetype="dotted") +
  geom_point(size = 5) +
  theme_classic() +
  xlim(c(-3.5,3.5)) + 
  scale_color_manual(values=c("21B9CA", 
                              "grey", 
                              "orange3")) +
  geom_text(aes(label = sig.label),
            vjust = 0, 
            nudge_y = 0.10,
            size = 5) +
  theme(text = element_text(size = 20),
        legend.position = 'none')  
ggsave('neuropeptides/avp.oxt/avp/avp.neurons.0.DomvsSub.volcano.pdf',
       width = 5,
       height = 5)

### dom vs sub
### group vector
## vector of factor specifies groups
## corresponds to columns in counts
## subset to groups of interest 
avp.neuron.DomvsSub.vector = burtoni.snseq.combined.sct.all.avp.neurons.recluster.expression 

## create vector of factor
avp.neuron.DomvsSub.vector.list = avp.neuron.DomvsSub.vector %>% 
  mutate(avp.neuron.orig.ident = orig.ident %>% 
           as.factor()) %>% 
  pull(avp.neuron.orig.ident) %>% 
  droplevels()

### counts matrix
## raw read count matrix
## rows = genes, columns = cells
# subset to groups of interest (0 and 1 of AVP neurons)
# only keep 500 variable genes
# set negative values to 0
avp.neuron.DomvsSub.vector.count = GetAssayData(burtoni.snseq.combined.sct.all.avp.neurons.recluster) %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column('gene') %>% 
  select(c(gene,
           avp.neuron.DomvsSub.vector %>% 
             pull(Cell.id))) %>% 
  filter(gene %in% avp.neuron.group.topgenes) %>% 
  column_to_rownames('gene') %>% 
  as.matrix() %>% 
  pmax(0)

### Detecting the DE genes
avp.neuron.DomvsSub.results = DEsingle(counts = avp.neuron.DomvsSub.vector.count, 
                                group = avp.neuron.DomvsSub.vector.list,
                                parallel = TRUE)

# Dividing the DE genes into 3 categories at threshold of FDR < 0.1
avp.neuron.DomvsSub.results.classified <- DEtype(results = avp.neuron.DomvsSub.results, 
                                          threshold = 0.1)


# create color groups
avp.neuron.DomvsSub.results.classified = avp.neuron.DomvsSub.results.classified %>% 
  mutate(Sig = ifelse(pvalue <= 0.01,
                      'Sig',
                      'Not Sig'),
         Direction = total_mean_1 - total_mean_2,
         Direction.type = ifelse(Direction > 0,
                                 'up',
                                 'down'),
         Sig.direction = ifelse(Sig == 'Sig',
                                Direction.type,
                                Sig))


##graph volcano plot
# create log fold change and -log pvalue
avp.neuron.DomvsSub.results.classified %>% 
  rownames_to_column('gene') %>% 
  mutate(logFoldChange = log(foldChange)) %>% 
  mutate(logFoldChange = ifelse(logFoldChange < -3.5,
                -3.3,
                logFoldChange),
         logFoldChange = ifelse(logFoldChange > 3,
                                3.3,
                                logFoldChange)) %>% 
  mutate(sig.label = ifelse(Sig == 'Sig',
                            gene,
                            '')) %>% 
  ggplot(aes(x = logFoldChange,
             y = -log10(pvalue),
             color = Sig.direction))+
  geom_vline(xintercept = 3.2, 
             linetype="dotted")+
  geom_vline(xintercept = -3.2,
             linetype="dotted") +
  geom_point(size = 5) +
  theme_classic() +
  xlim(c(-3.5,3.5)) + 
  scale_color_manual(values=c("21B9CA", 
                              "grey", 
                              "orange3")) +
  geom_text(aes(label = sig.label),
            vjust = 0, 
            nudge_y = 0.10,
            size = 5) +
  theme(text = element_text(size = 20),
        legend.position = 'none')  
ggsave('neuropeptides/avp.oxt/avp/avp.neurons.DomvsSub.volcano.pdf',
       width = 5,
       height = 5)

#TCF4, mpped1, pcdh1b, gh1, prl, COX3, ND4
#### try edgeR approach ####
#https://github.com/csoneson/conquer_comparison/blob/master/scripts/apply_edgeRQLFDetRate.R

#load edgeR
suppressPackageStartupMessages(library(edgeR))

#create function
run_edgeRQLFDetRate <- function(L) {
  message("edgeRQLFDetRate")
  session_info <- sessionInfo()
  tryCatch({
    timing <- system.time({
      dge <- DGEList(L$count, group = L$condt)
      dge <- calcNormFactors(dge)
      cdr <- scale(colMeans(L$count > 0))
      design <- model.matrix(~ cdr + L$condt)
      dge <- estimateDisp(dge, design = design)
      fit <- glmQLFit(dge, design = design)
      qlf <- glmQLFTest(fit)
      tt <- topTags(qlf, n = Inf)
    })
    
    plotBCV(dge)
    plotQLDisp(fit)
    hist(tt$table$PValue, 50)
    hist(tt$table$FDR, 50)
    limma::plotMDS(dge, col = as.numeric(as.factor(L$condt)), pch = 19)
    plotSmear(qlf)
    
    list(session_info = session_info,
         timing = timing,
         tt = tt,
         df = data.frame(pval = tt$table$PValue,
                         padj = tt$table$FDR,
                         row.names = rownames(tt$table)))
  }, error = function(e) {
    "edgeRQLFDetRate results could not be calculated"
    list(session_info = session_info)
  })
}


## run with all avp neurons
#need list with L with  count and condt
avp.neuron.DomvsSub.vector.edgeR = list(count = avp.neuron.DomvsSub.vector.count,
                                 condt = avp.neuron.DomvsSub.vector.list)


#run function  
avp.neuron.DomvsSub.vector.edgeR.results = run_edgeRQLFDetRate(avp.neuron.DomvsSub.vector.edgeR)

# save results to dataframe
avp.neuron.DomvsSub.vector.edgeR.results.df = avp.neuron.DomvsSub.vector.edgeR.results$tt@.Data[[1]]

#add color for significance  
avp.neuron.DomvsSub.vector.edgeR.results.df = avp.neuron.DomvsSub.vector.edgeR.results.df %>% 
  mutate(Sig = ifelse(FDR <= 0.05,
                      'Sig',
                      'Not Sig'),
         Direction.type = ifelse(logFC > 0,
                                 'up',
                                 'down'),
         Sig.direction = ifelse(Sig == 'Sig',
                                Direction.type,
                                Sig))


##graph volcano plot
# create log fold change and -log pvalue
avp.neuron.DomvsSub.vector.edgeR.results.df %>% 
  rownames_to_column('gene') %>% 
  mutate(sig.label = ifelse(Sig == 'Sig',
                            gene,
                            '')) %>% 
  ggplot(aes(x = logFC,
             y = -log10(FDR),
             color = Sig.direction))+
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
        legend.position = 'none')
ggsave('neuropeptides/avp.oxt/avp/avp.neurons.DomvsSub.volcano.edgeR.pdf',
       width = 5,
       height = 5)

## run with avp cluster 0
#need list with L with  count and condt
avp.neuron.0.vector.edgeR = list(count = avp.neuron.0.vector.count,
            condt = avp.neuron.0.vector.list)


#run function  
avp.neuron.0.vector.edgeR.results = run_edgeRQLFDetRate(avp.neuron.0.vector.edgeR)

# save results to dataframe
avp.neuron.0.vector.edgeR.results.df = avp.neuron.0.vector.edgeR.results$tt@.Data[[1]]

#add color for significance  
avp.neuron.0.vector.edgeR.results.df = avp.neuron.0.vector.edgeR.results.df %>% 
  mutate(Sig = ifelse(FDR <= 0.05,
                      'Sig',
                      'Not Sig'),
         Direction.type = ifelse(logFC > 0,
                                 'up',
                                 'down'),
         Sig.direction = ifelse(Sig == 'Sig',
                                Direction.type,
                                Sig))


##graph volcano plot
# create log fold change and -log pvalue
avp.neuron.0.vector.edgeR.results.df %>% 
  rownames_to_column('gene') %>% 
  mutate(sig.label = ifelse(Sig == 'Sig',
                            gene,
                            '')) %>% 
  ggplot(aes(x = logFC,
             y = -log10(FDR),
             color = Sig.direction))+
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
        legend.position = 'none')
ggsave('neuropeptides/avp.oxt/avp/avp.neurons.0.DomvsSub.volcano.edgeR.pdf',
       width = 5,
       height = 5)

## compare to desingle
avp.neuron.0.combined.results = full_join(avp.neuron.0.vector.edgeR.results.df %>% 
                                            rownames_to_column('gene'),
                                          avp.neuron.0.results.classified %>% 
                                            rownames_to_column('gene'),
                                          by = c('gene'))

# pvalue
avp.neuron.0.combined.results %>% 
  ggplot(aes(x = -log10(PValue),
             y= -log10(pvalue))) +
  geom_point() +
  ggtitle('EdgeR vs DEsingle AVP 0 pvalue') +
  theme_classic()
ggsave('neuropeptides/avp.oxt/avp/EdgeR vs DEsingle AVP 0 pvalue.pdf',
       width = 5,
       height = 5)


# FDR
avp.neuron.0.combined.results %>% 
  ggplot(aes(x = FDR,
             y= pvalue.adj.FDR)) +
  geom_point() +
  ggtitle('EdgeR vs DEsingle AVP 0 FDR') +
  theme_classic()
ggsave('neuropeptides/avp.oxt/avp/EdgeR vs DEsingle AVP 0 FDR.pdf',
       width = 5,
       height = 5)

#logFC
avp.neuron.0.combined.results %>% 
  ggplot(aes(x = logFC,
             y= log10(foldChange))) +
  geom_point() +
  geom_smooth(method = 'lm') +
  ggtitle('EdgeR vs DEsingle AVP 0 logFC') +
  theme_classic()
ggsave('neuropeptides/avp.oxt/avp/EdgeR vs DEsingle AVP 0 logFC.pdf',
       width = 5,
       height = 5)



#### gnrh ####
## gnrh
#avp vs oxt
## umap plots
FeaturePlot(object = burtoni.snseq.combined.sct.all, 
            features = c("ENSONIG00000011023"),
            pt.size = 0.5,
            min.cutoff = 4,
            max.cutoff = 5,
            order = TRUE)
ggsave('gnrh/gnrh.umap.png',
       width = 20,
       height = 20)
#across samples
FeaturePlot(object = burtoni.snseq.combined.sct.all, 
            features = c("ENSONIG00000011023"),
            pt.size = 0.5,
            min.cutoff = 4,
            max.cutoff = 6,
            order = TRUE,
            split.by = "orig.ident")
ggsave('gnrh/gnrh.subvsdom.umap.png',
       width = 10,
       height = 10)

##vinplot
#gnrh
VlnPlot(burtoni.snseq.combined.sct.all, 
        features = c("ENSONIG00000011023"),
        split.by = "orig.ident") +
  ylim(4,25)
ggsave('gnrh/gnrh.subvsdom.vlnplot.png',
       width = 10,
       height = 10)

##check number of cells expressing gene above threshold
#gnrh
#351
sum(GetAssayData(object = burtoni.snseq.combined.sct.all, 
                 slot = "data")["ENSONIG00000011023",]>4)

### subset data to relevant cells
burtoni.snseq.combined.sct.all.gnrh = subset(burtoni.snseq.combined.sct.all,
                                             subset = ENSONIG00000011023 > 4 )

## run PCA, UMAP, and cluster 
#use 0.4 resolution
burtoni.snseq.combined.sct.all.gnrh.recluster = burtoni.snseq.combined.sct.all.gnrh %>% 
  RunPCA() %>%
  FindNeighbors(dims = 1:15) %>%
  RunUMAP(dims = 1:15) %>%
  FindClusters(resolution = 0.6)

## check cells per group per cluster
#get table
Gnrh.cells.per.cluster.table = table(burtoni.snseq.combined.sct.all.gnrh.recluster@active.ident, 
                                     burtoni.snseq.combined.sct.all.gnrh.recluster@meta.data$orig.ident) %>% 
  as.data.frame.matrix() %>% 
  rownames_to_column("cluster.id") %>% 
  pivot_longer(cols = c("dom_burtoni_snseq",
                        "sub_burtoni_snseq"),
               names_to = "orig.ident",
               values_to = "count") %>% 
  group_by(orig.ident) %>% 
  mutate(total = sum(count)) %>% 
  ungroup() %>% 
  mutate(percentage = 100*count/total)

#graph
Gnrh.cells.per.cluster.table %>% 
  as.data.frame() %>% 
  mutate(cluster.id = as.numeric(cluster.id)) %>% 
  ggplot(aes(x = cluster.id,
             y = count,
             group = orig.ident,
             color = orig.ident)) +
  geom_point() +
  theme_bw()
ggsave('gnrh/Cells.per.cluster.DomvsSub.gnrh.png',
       width = 10,
       height = 10)

#percent
Gnrh.cells.per.cluster.table %>% 
  as.data.frame() %>% 
  mutate(cluster.id = as.numeric(cluster.id)) %>% 
  ggplot(aes(x = cluster.id,
             y = percentage,
             group = orig.ident,
             color = orig.ident)) +
  geom_point() +
  theme_bw()
ggsave('gnrh/Cells.per.cluster.DomvsSub.gnrh.percentage.png',
       width = 10,
       height = 10)

### graph 
##dom vs sub
DimPlot(burtoni.snseq.combined.sct.all.gnrh.recluster, 
        reduction = "umap", 
        group.by = "orig.ident")
ggsave('gnrh/DomVsSub.dimplot.gnrh.all.png')

## clusters
DimPlot(burtoni.snseq.combined.sct.all.gnrh.recluster, 
        reduction = "umap", 
        label = TRUE,
        repel = TRUE)
ggsave('gnrh/Clusters.dimplot.gnrh.all.png',
       width = 10,
       height = 10)

###dotplot
## cell type
DotPlot(burtoni.snseq.combined.sct.all.gnrh.recluster, 
        features = burtoni.snseq.combined.sct.all.gnrh.recluster@meta.data %>% 
          select(ends_with(".score1")) %>% 
          colnames(), 
        cols = c("grey", 
                 "red"), 
        dot.scale = 8,
        col.min = 1,
        dot.min = .4) + 
  RotatedAxis()
ggsave('gnrh/gnrh.celltype.dotplot.png',
       width = 10,
       height = 10)



### clustree
# cluster across resolutions
burtoni.snseq.combined.sct.all.gnrh.recluster.clustree <- Seurat::FindClusters(object = burtoni.snseq.combined.sct.all.gnrh.recluster, 
                                                                               resolution = resolution.range.reduced)
#check data
head(burtoni.snseq.combined.sct.all.gnrh.recluster.clustree[[]])
#set presentation colors
presentation.color <- c('#66c2a5',
                        '#fc8d62',
                        '#8da0cb',
                        '#e78ac3',
                        '#a6d854',
                        '#ffd92f',
                        '#e5c494',
                        '#b3b3b3')
#clustree
clustree(burtoni.snseq.combined.sct.all.gnrh.recluster.clustree, 
         prefix = "integrated_snn_res.",
         node_colour = 'cluster',
         node_size_range = c(10,20),
         scale_node_text = TRUE) +
  scale_edge_color_continuous(low = "black", 
                              high = "black") +
  scale_color_manual(values = presentation.color)+
  theme(legend.position = "bottom")
ggsave('gnrh/gnrh.clustree.png',
       width = 10,
       height = 10)

###pca
# clusters
DimPlot(burtoni.snseq.combined.sct.all.gnrh.recluster.clustree, 
        reduction = "pca", 
        label = TRUE,
        repel = TRUE)
ggsave('gnrh/gnrh.pca.PC1.PC2.png',
       width = 10,
       height = 10)
#PC 2 and 3
DimPlot(burtoni.snseq.combined.sct.all.gnrh.recluster.clustree, 
        reduction = "pca", 
        label = TRUE,
        repel = TRUE,
        dims = c(3,2))
ggsave('gnrh/gnrh.pca.PC2.PC3.png',
       width = 10,
       height = 10)


#heatmap PCA
DimHeatmap(burtoni.snseq.combined.sct.all.gnrh.recluster.clustree, 
           dims = 1:6, 
           cells = 500, 
           balanced = TRUE)
ggsave('gnrh/gnrh.pca.heatmap.png',
       width = 10,
       height = 10)
#loadings
VizDimLoadings(burtoni.snseq.combined.sct.all.gnrh.recluster.clustree, 
               dims = 1:6, 
               reduction = "pca")
ggsave('gnrh/gnrh.pca.loadings.png',
       width = 10,
       height = 10)


### identify marker genes
## all
burtoni.snseq.combined.sct.all.avp.oxt.recluster.clustree.markers <- FindAllMarkers(burtoni.snseq.combined.sct.all.avp.oxt.recluster.clustree, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#### sst ####
## sst
#avp vs oxt
## umap plots
FeaturePlot(object = burtoni.snseq.combined.sct.all, 
            features = c("ENSONIG00000033642"),
            pt.size = 0.5,
            min.cutoff = 4,
            max.cutoff = 5,
            order = TRUE)
ggsave('sst/sst.umap.png',
       width = 20,
       height = 20)
#across samples
FeaturePlot(object = burtoni.snseq.combined.sct.all, 
            features = c("ENSONIG00000033642"),
            pt.size = 0.5,
            min.cutoff = 4,
            max.cutoff = 6,
            order = TRUE,
            split.by = "orig.ident")
ggsave('sst/sst.subvsdom.umap.png',
       width = 10,
       height = 10)

##vinplot
#sst
VlnPlot(burtoni.snseq.combined.sct.all, 
        features = c("ENSONIG00000033642"),
        split.by = "orig.ident") +
  ylim(4,25)
ggsave('sst/sst.subvsdom.vlnplot.png',
       width = 10,
       height = 10)

##check number of cells expressing gene above threshold
#sst
#334
sum(GetAssayData(object = burtoni.snseq.combined.sct.all, 
                 slot = "data")["ENSONIG00000033642",]>4)

### subset data to relevant cells
burtoni.snseq.combined.sct.all.sst = subset(burtoni.snseq.combined.sct.all,
                                            subset = ENSONIG00000033642 > 4 )

## run PCA, UMAP, and cluster 
#use 0.4 resolution
burtoni.snseq.combined.sct.all.sst.recluster = burtoni.snseq.combined.sct.all.sst %>% 
  RunPCA() %>%
  FindNeighbors(dims = 1:15) %>%
  RunUMAP(dims = 1:15) %>%
  FindClusters(resolution = 0.4)

## check cells per group per cluster
#get table
sst.cells.per.cluster.table = table(burtoni.snseq.combined.sct.all.sst.recluster@active.ident, 
                                    burtoni.snseq.combined.sct.all.sst.recluster@meta.data$orig.ident) %>% 
  as.data.frame.matrix() %>% 
  rownames_to_column("cluster.id") %>% 
  pivot_longer(cols = c("dom_burtoni_snseq",
                        "sub_burtoni_snseq"),
               names_to = "orig.ident",
               values_to = "count") %>% 
  group_by(orig.ident) %>% 
  mutate(total = sum(count)) %>% 
  ungroup() %>% 
  mutate(percentage = 100*count/total)

#graph
sst.cells.per.cluster.table %>% 
  as.data.frame() %>% 
  mutate(cluster.id = as.numeric(cluster.id)) %>% 
  ggplot(aes(x = cluster.id,
             y = count,
             group = orig.ident,
             color = orig.ident)) +
  geom_point() +
  theme_bw()
ggsave('sst/Cells.per.cluster.DomvsSub.sst.png',
       width = 10,
       height = 10)

#percent
sst.cells.per.cluster.table %>% 
  as.data.frame() %>% 
  mutate(cluster.id = as.numeric(cluster.id)) %>% 
  ggplot(aes(x = cluster.id,
             y = percentage,
             group = orig.ident,
             color = orig.ident)) +
  geom_point() +
  theme_bw()
ggsave('sst/Cells.per.cluster.DomvsSub.sst.percentage.png',
       width = 10,
       height = 10)

### graph 
##dom vs sub
DimPlot(burtoni.snseq.combined.sct.all.sst.recluster, 
        reduction = "umap", 
        group.by = "orig.ident")
ggsave('sst/DomVsSub.dimplot.sst.all.png')

## clusters
DimPlot(burtoni.snseq.combined.sct.all.sst.recluster, 
        reduction = "umap", 
        label = TRUE,
        repel = TRUE)
ggsave('sst/Clusters.dimplot.sst.all.png',
       width = 10,
       height = 10)

###dotplot
## cell type
DotPlot(burtoni.snseq.combined.sct.all.sst.recluster, 
        features = burtoni.snseq.combined.sct.all.sst.recluster@meta.data %>% 
          select(ends_with(".score1")) %>% 
          colnames(), 
        cols = c("grey", 
                 "red"), 
        dot.scale = 8,
        col.min = 1,
        dot.min = .4) + 
  RotatedAxis()
ggsave('sst/sst.celltype.dotplot.png',
       width = 10,
       height = 10)



### clustree
# cluster across resolutions
burtoni.snseq.combined.sct.all.sst.recluster.clustree <- Seurat::FindClusters(object = burtoni.snseq.combined.sct.all.sst.recluster, 
                                                                              resolution = resolution.range.reduced)
#check data
head(burtoni.snseq.combined.sct.all.sst.recluster.clustree[[]])
#set presentation colors
presentation.color <- c('#66c2a5',
                        '#fc8d62',
                        '#8da0cb',
                        '#e78ac3',
                        '#a6d854',
                        '#ffd92f',
                        '#e5c494',
                        '#b3b3b3')
#clustree
clustree(burtoni.snseq.combined.sct.all.sst.recluster.clustree, 
         prefix = "integrated_snn_res.",
         node_colour = 'cluster',
         node_size_range = c(10,20),
         scale_node_text = TRUE) +
  scale_edge_color_continuous(low = "black", 
                              high = "black") +
  scale_color_manual(values = presentation.color)+
  theme(legend.position = "bottom")
ggsave('sst/sst.clustree.png',
       width = 10,
       height = 10)

###pca
# clusters
DimPlot(burtoni.snseq.combined.sct.all.sst.recluster.clustree, 
        reduction = "pca", 
        label = TRUE,
        repel = TRUE)
ggsave('sst/sst.pca.PC1.PC2.png',
       width = 10,
       height = 10)
#PC 2 and 3
DimPlot(burtoni.snseq.combined.sct.all.sst.recluster.clustree, 
        reduction = "pca", 
        label = TRUE,
        repel = TRUE,
        dims = c(3,2))
ggsave('sst/sst.pca.PC2.PC3.png',
       width = 10,
       height = 10)


#heatmap PCA
DimHeatmap(burtoni.snseq.combined.sct.all.sst.recluster.clustree, 
           dims = 1:6, 
           cells = 500, 
           balanced = TRUE)
ggsave('sst/sst.pca.heatmap.png',
       width = 10,
       height = 10)
#loadings
VizDimLoadings(burtoni.snseq.combined.sct.all.sst.recluster.clustree, 
               dims = 1:6, 
               reduction = "pca")
ggsave('sst/sst.pca.loadings.png',
       width = 10,
       height = 10)



#### Esr1 neurons ####
### esr1 neurons
#esr1
VlnPlot(burtoni.snseq.combined.sct.all.neurons, 
        features = c("esr1"),
        group.by = "orig.ident")
ggsave('neuropeptides/esr1/esr1.neurons.DomvsSub.png',
       width = 10,
       height = 10)

### subset data to relevant cells
burtoni.snseq.combined.sct.all.esr1.neurons = burtoni.snseq.combined.sct.all.neurons

#set idents
Idents(object = burtoni.snseq.combined.sct.all.esr1.neurons) <- "Cell.type"

#subset
burtoni.snseq.combined.sct.all.esr1.neurons = subset(burtoni.snseq.combined.sct.all.esr1.neurons,
                                                     subset = esr1 > 2)


## run PCA, UMAP, and cluster 
#use 0.4 resolution
burtoni.snseq.combined.sct.all.esr1.neurons.recluster = burtoni.snseq.combined.sct.all.esr1.neurons %>% 
  RunPCA() %>%
  FindNeighbors(dims = 1:15) %>%
  RunUMAP(dims = 1:15) %>%
  FindClusters(resolution = 0.8)

### clustree
# cluster across resolutions
burtoni.snseq.combined.sct.all.esr1.neurons.clustree <- Seurat::FindClusters(object = burtoni.snseq.combined.sct.all.esr1.neurons.recluster, 
                                                                             resolution = resolution.range.reduced)


#clustree
clustree(burtoni.snseq.combined.sct.all.esr1.neurons.clustree, 
         prefix = "integrated_snn_res.",
         node_colour = 'cluster',
         node_size_range = c(10,20),
         scale_node_text = TRUE) +
  scale_edge_color_continuous(low = "black", 
                              high = "black") +
  theme(legend.position = "bottom")
ggsave('neuropeptides/esr1/esr1.neurons.clustree.png',
       width = 10,
       height = 10)


### graph 
##dom vs sub
DimPlot(burtoni.snseq.combined.sct.all.esr1.neurons.recluster, 
        reduction = "umap", 
        group.by = "orig.ident")
ggsave('neuropeptides/esr1/DomVsSub.dimplot.esr1.neurons.all.png')

## clusters
DimPlot(burtoni.snseq.combined.sct.all.esr1.neurons.recluster, 
        reduction = "umap", 
        label = TRUE,
        repel = TRUE,
        pt.size = 5)
ggsave('neuropeptides/esr1/Clusters.dimplot.esr1.neurons.all.png',
       width = 10,
       height = 10)

## cell type
DimPlot(burtoni.snseq.combined.sct.all.esr1.neurons.recluster, 
        reduction = "umap", 
        group.by = "Cell.type")
ggsave('neuropeptides/esr1/Celltype.dimplot.esr1.neurons.all.png')



##expression level
burtoni.snseq.combined.sct.all.esr1.neurons.recluster.expression = full_join(full_join(burtoni.snseq.combined.sct.all.esr1.neurons.recluster@reductions$umap@cell.embeddings %>% 
                                                                                         as.data.frame() %>% 
                                                                                         rownames_to_column("Cell.id"),
                                                                                       burtoni.snseq.combined.sct.all.esr1.neurons.recluster@meta.data %>% 
                                                                                         rownames_to_column("Cell.id")),
                                                                             burtoni.snseq.combined.sct.all.esr1.neurons.recluster@assays$integrated@scale.data %>% 
                                                                               as.data.frame() %>% 
                                                                               filter(rownames(burtoni.snseq.combined.sct.all.esr1.neurons.recluster@assays$integrated@scale.data) %in% c('avp',
                                                                                                                                                                                          'oxt')) %>% 
                                                                               t() %>% as.data.frame() %>% 
                                                                               rownames_to_column('Cell.id')) %>% 
  full_join(burtoni.snseq.combined.sct.all.esr1.neurons.recluster@reductions$pca@cell.embeddings %>% 
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

# seperation by cell type
table(burtoni.snseq.combined.sct.all.esr1.neurons.recluster.expression %>%
        select(Cell.type,
               integrated_snn_res.0.8))

# difference social status
table(burtoni.snseq.combined.sct.all.esr1.neurons.recluster.expression %>%
        select(integrated_snn_res.0.8,
               orig.ident))

# difference social status
# percent
# 167 dom to 106 sub
table(burtoni.snseq.combined.sct.all.esr1.neurons.recluster.expression %>%
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
  geom_vline(xintercept = 61.2) +
  geom_vline(xintercept = 50, 
             linetype="dotted") +
  geom_vline(xintercept = 72.4, 
             linetype="dotted") +
  geom_point(size = 10) +
  theme_classic()+ 
  theme(text = element_text(size = 20)) +
  xlim(40,80)
ggsave('neuropeptides/esr1/Clusters.DomvsSub.ratio.esr1.neurons.all.png',
       width = 10,
       height = 10)


#### subset neuropeptides neurons ####
### avp, oxt, gnrh, sst, neurons
### subset data to relevant cells
burtoni.snseq.combined.sct.all.neuropeptides.neurons = burtoni.snseq.combined.sct.all

#set idents
Idents(object = burtoni.snseq.combined.sct.all.neuropeptides.neurons) <- "Cell.type"

#subset neurons
burtoni.snseq.combined.sct.all.neuropeptides.neurons = subset(burtoni.snseq.combined.sct.all.neuropeptides.neurons,
                                                        idents = c("neurons",
                                                                   "excitatory",
                                                                   "inhibitory"))

#subset expression
burtoni.snseq.combined.sct.all.neuropeptides.neurons = subset(burtoni.snseq.combined.sct.all.neuropeptides.neurons,
                                                        subset = avp > 2 | oxt > 2 | ENSONIG00000011023 > 2 | ENSONIG00000033642 > 2)


## run PCA, UMAP, and cluster 
#use 0.4 resolution
burtoni.snseq.combined.sct.all.neuropeptides.neurons.recluster = burtoni.snseq.combined.sct.all.neuropeptides.neurons %>% 
  RunPCA() %>%
  FindNeighbors(dims = 1:15) %>%
  RunUMAP(dims = 1:15) %>%
  FindClusters(resolution = 0.4)

### clustree
# cluster across resolutions
burtoni.snseq.combined.sct.all.neuropeptides.neurons.clustree <- Seurat::FindClusters(object = burtoni.snseq.combined.sct.all.neuropeptides.neurons.recluster, 
                                                                                resolution = resolution.range.reduced.2)
#check data
# head(burtoni.snseq.combined.sct.all.neuropeptides.neurons.clustree[[]])

#clustree
clustree(burtoni.snseq.combined.sct.all.neuropeptides.neurons.clustree, 
         prefix = "integrated_snn_res.",
         node_colour = 'cluster',
         node_size_range = c(10,20),
         scale_node_text = TRUE) +
  scale_edge_color_continuous(low = "black", 
                              high = "black") +
  theme(legend.position = "bottom")
ggsave('neuropeptides/neuropeptides.neurons.clustree.png',
       width = 10,
       height = 10)


### graph 
##dom vs sub
DimPlot(burtoni.snseq.combined.sct.all.neuropeptides.neurons.recluster, 
        reduction = "umap", 
        group.by = "orig.ident")
ggsave('neuropeptides/DomVsSub.dimplot.neuropeptides.neurons.all.png')

## clusters
DimPlot(burtoni.snseq.combined.sct.all.neuropeptides.neurons.recluster, 
        reduction = "umap", 
        label = TRUE,
        repel = TRUE)
ggsave('neuropeptides/Clusters.dimplot.neuropeptides.neurons.all.png',
       width = 10,
       height = 10)

##expression level
burtoni.snseq.combined.sct.all.neuropeptides.neurons.recluster.expression = full_join(full_join(burtoni.snseq.combined.sct.all.neuropeptides.neurons.recluster@reductions$umap@cell.embeddings %>% 
                                                                                            as.data.frame() %>% 
                                                                                            rownames_to_column("Cell.id"),
                                                                                          burtoni.snseq.combined.sct.all.neuropeptides.neurons.recluster@meta.data %>% 
                                                                                            rownames_to_column("Cell.id")),
                                                                                burtoni.snseq.combined.sct.all.neuropeptides.neurons.recluster@assays$integrated@scale.data %>% 
                                                                                  as.data.frame() %>% 
                                                                                  filter(rownames(burtoni.snseq.combined.sct.all.neuropeptides.neurons.recluster@assays$integrated@scale.data) %in% c('avp',
                                                                                                                                                                                                'oxt',
                                                                                                                                                                                                'ENSONIG00000011023',
                                                                                                                                                                                                'ENSONIG00000033642')) %>% 
                                                                                  t() %>% as.data.frame() %>% 
                                                                                  rownames_to_column('Cell.id')) %>% 
  full_join(burtoni.snseq.combined.sct.all.neuropeptides.neurons.recluster@reductions$pca@cell.embeddings %>% 
              as.data.frame() %>% 
              select(PC_1,
                     PC_2,
                     PC_3,
                     PC_4,
                     PC_5) %>% 
              rownames_to_column("Cell.id")) %>% 
  mutate(Cell.type.neuropeptides = case_when(avp > 2 ~ 'avp',
                                             oxt > 2 ~ 'oxt',
                                             ENSONIG00000011023 > 2 ~ 'gnrh',
                                             ENSONIG00000033642 > 2 ~ 'sst',
                                             TRUE ~ "all"))


##graph
burtoni.snseq.combined.sct.all.neuropeptides.neurons.recluster.expression %>% 
  ggplot(aes(x = UMAP_1,
             y = UMAP_2,
             color = integrated_snn_res.0.4,
             shape = Cell.type.neuropeptides)) +
  geom_point() +
  theme_bw()
ggsave('neuropeptides/Umap.neuropeptides.neurons.all.png',
       width = 10,
       height = 10)

#poster
burtoni.snseq.combined.sct.all.neuropeptides.neurons.recluster.expression %>% 
  ggplot(aes(x = UMAP_1,
             y = UMAP_2,
             color = Cell.type.neuropeptides)) +
  geom_point(size = 3) +
  theme_classic() +
  scale_color_manual(values = c('orange3',
                                "purple",
                                "21B9CA",
                                "yellow3"))+ 
  labs(color='Neuropeptides') +
  theme(legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        legend.position = c(0.8, 0.2)) +
  guides(colour = guide_legend(override.aes = list(size=10)))
ggsave('neuropeptides/Umap.neuropeptides.neurons.all.celltype.pdf',
       width = 12,
       height = 10)

## cell count per cluster
table(burtoni.snseq.combined.sct.all.neuropeptides.neurons.recluster.expression %>% 
        select(integrated_snn_res.0.4,
               Cell.type.neuropeptides))

table(burtoni.snseq.combined.sct.all.neuropeptides.neurons.recluster.expression %>% 
        select(integrated_snn_res.0.4,
               Cell.type.neuropeptides,
               orig.ident))

table(burtoni.snseq.combined.sct.all.neuropeptides.neurons.recluster.expression %>% 
        select(integrated_snn_res.0.4,
               orig.ident))

table(burtoni.snseq.combined.sct.all.neuropeptides.neurons.recluster.expression %>% 
        select(Cell.type.neuropeptides,
               orig.ident))

##chi squared test
chisq.test.nueropeptides.neurons = chisq.test(table(burtoni.snseq.combined.sct.all.neuropeptides.neurons.recluster.expression %>% 
                                                      select(Cell.type.neuropeptides,
                                                             orig.ident)) %>% 
                                                t())
#X-squared = 26.909, df = 3, p-value = 6.152e-06

## test if 50% for each
chisq.test(table(burtoni.snseq.combined.sct.all.neuropeptides.neurons.recluster.expression %>% 
                                                      select(Cell.type.neuropeptides,
                                                             orig.ident) %>% 
                   filter(Cell.type.neuropeptides == 'oxt')),
           p =c(0.5,0.5))



## graph

table(burtoni.snseq.combined.sct.all.neuropeptides.neurons.recluster.expression %>% 
        select(Cell.type.neuropeptides,
               orig.ident)) %>% 
        as_tibble() %>% 
  ggplot(aes(x = Cell.type.neuropeptides,
             y = n,
             group = orig.ident,
             fill = orig.ident)) +
  geom_bar(stat = 'identity',
           position=position_dodge(width=0.8),
           width = 0.7) +
  theme_classic() +
  scale_fill_manual(values = c('orange3',
                               "21B9CA"))
ggsave('neuropeptides/Bargraph.neuropeptides.neurons.all.celltype.count.pdf',
       width = 5,
       height = 5)

### PCA
#1 vs 2
burtoni.snseq.combined.sct.all.neuropeptides.neurons.recluster.expression %>% 
  ggplot(aes(x = PC_1,
             y = PC_2,
             color = integrated_snn_res.0.4,
             shape = Cell.type.neuropeptides)) +
  geom_point() +
  theme_bw()
ggsave('neuropeptides/PCA.1.2.neuropeptides.neurons.all.png',
       width = 10,
       height = 10)
#2 vs 3
burtoni.snseq.combined.sct.all.neuropeptides.neurons.recluster.expression %>% 
  ggplot(aes(x = PC_3,
             y = PC_2,
             shape = integrated_snn_res.0.4,
             color = Cell.type.neuropeptides)) +
  geom_point() +
  theme_bw()
ggsave('neuropeptides/PCA.2.3.expression.neuropeptides.neurons.all.png',
       width = 10,
       height = 10)



#### compare neuropeptide membership ####
### avp.oxt vs avp.oxt.neurons
png('avp.oxt/avp/heatmap avp.oxt vs avp.oxt.neurons.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)
table(left_join(burtoni.snseq.combined.sct.all.avp.oxt.neurons.recluster.expression %>% 
                  select(Cell.id,
                         integrated_snn_res.0.8,
                         Oxt.cell) %>% 
                  rename(Avp.oxt.neuron.cluster = integrated_snn_res.0.8),
                burtoni.snseq.combined.sct.all.avp.oxt.recluster.expression %>% 
                  select(Cell.id,
                         integrated_snn_res.1) %>% 
                  rename(Avp.oxt.cluster = integrated_snn_res.1)) %>% 
        filter(Oxt.cell != 'both') %>% 
        select(Avp.oxt.neuron.cluster,
               Avp.oxt.cluster,
               Oxt.cell))  %>% 
  unclass() %>% 
  as.data.frame() %>% 
  as.matrix() %>% 
  pheatmap(main = 'avp.oxt vs avp.oxt.neurons',
           color = colorRampPalette(c("white", "red"))(50))
dev.off()

### avp.neurons vs avp.oxt.neurons
png('avp.oxt/avp/heatmap avp.neurons vs avp.oxt.neurons.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)
table(right_join(burtoni.snseq.combined.sct.all.avp.oxt.neurons.recluster.expression %>% 
                   select(Cell.id,
                          integrated_snn_res.0.8,
                          Oxt.cell) %>% 
                   rename(Avp.oxt.neuron.cluster = integrated_snn_res.0.8),
                 burtoni.snseq.combined.sct.all.avp.neurons.recluster.expression %>% 
                   select(Cell.id,
                          integrated_snn_res.1) %>% 
                   rename(Avp.neuron.cluster = integrated_snn_res.1)) %>% 
        filter(Oxt.cell != 'both') %>% 
        select(Avp.oxt.neuron.cluster,
               Avp.neuron.cluster,
               Oxt.cell))  %>% 
  unclass() %>% 
  as.data.frame() %>% 
  as.matrix() %>% 
  pheatmap(main = 'avp.neurons vs avp.oxt.neurons',
           color = colorRampPalette(c("white", "red"))(50))
dev.off()

#better color pallette
png('avp.oxt/avp/heatmap avp.neurons vs avp.oxt.neurons orange.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)
table(right_join(burtoni.snseq.combined.sct.all.avp.oxt.neurons.recluster.expression %>% 
                   select(Cell.id,
                          integrated_snn_res.0.8,
                          Oxt.cell) %>% 
                   rename(Avp.oxt.neuron.cluster = integrated_snn_res.0.8),
                 burtoni.snseq.combined.sct.all.avp.neurons.recluster.expression %>% 
                   select(Cell.id,
                          integrated_snn_res.1) %>% 
                   rename(Avp.neuron.cluster = integrated_snn_res.1)) %>% 
        filter(Oxt.cell != 'both') %>% 
        select(Avp.oxt.neuron.cluster,
               Avp.neuron.cluster,
               Oxt.cell))  %>% 
  unclass() %>% 
  as.data.frame() %>% 
  as.matrix() %>% 
  pheatmap(main = 'avp.neurons vs avp.oxt.neurons',
           color = colorRampPalette(c("white", "orange3"))(30))
dev.off()


### avp.neurons vs avp.oxt
png('avp.oxt/avp/heatmap avp.oxt vs avp.neurons.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)
table(right_join(burtoni.snseq.combined.sct.all.avp.oxt.recluster.expression %>% 
                   select(Cell.id,
                          integrated_snn_res.1) %>% 
                   rename(Avp.oxt.cluster = integrated_snn_res.1),
                 burtoni.snseq.combined.sct.all.avp.neurons.recluster.expression %>% 
                   select(Cell.id,
                          integrated_snn_res.1,
                          Oxt.cell) %>% 
                   rename(Avp.cluster.neuron = integrated_snn_res.1)) %>% 
        filter(Oxt.cell != 'both') %>% 
        select(Avp.oxt.cluster,
               Avp.cluster.neuron,
               Oxt.cell)) %>% 
  unclass() %>% 
  as.data.frame() %>% 
  as.matrix() %>% 
  pheatmap(main = 'avp.oxt vs avp.neurons',
           color = colorRampPalette(c("white", "red"))(50))
dev.off()

### avp.neurons vs neuropeptides all
png('neuropeptides/avp.oxt/avp/heatmap all neuropeptides vs avp.neurons.png',
    width = 5,
    height = 5,
    units = 'in',
    res = 720)
table(right_join(burtoni.snseq.combined.sct.all.neuropeptides.neurons.recluster.expression %>% 
                   select(Cell.id,
                          integrated_snn_res.0.4) %>% 
                   rename(neuropeptides.cluster = integrated_snn_res.0.4),
                 burtoni.snseq.combined.sct.all.avp.neurons.recluster.expression %>% 
                   select(Cell.id,
                          integrated_snn_res.0.8,
                          Oxt.cell) %>% 
                   rename(Avp.cluster.neuron = integrated_snn_res.0.8)) %>% 
        select(neuropeptides.cluster,
               Avp.cluster.neuron)) %>% 
  unclass() %>% 
  as.data.frame() %>% 
  as.matrix() %>% 
  pheatmap(color = colorRampPalette(c("white", "orange3"))(50),
           fontsize = 20,
           cex = 1,
           legend = FALSE)
dev.off()

## updated
png('neuropeptides/avp.oxt/avp/heatmap all neuropeptides vs avp.neurons updated.png',
    width = 5,
    height = 5,
    units = 'in',
    res = 720)
table(right_join(burtoni.snseq.combined.sct.all.neuropeptides.neurons.recluster.expression %>% 
                   select(Cell.id,
                          integrated_snn_res.0.4) %>% 
                   rename(neuropeptides.cluster = integrated_snn_res.0.4),
                 burtoni.snseq.combined.sct.all.avp.neurons.recluster.expression %>% 
                   select(Cell.id,
                          integrated_snn_res.0.8,
                          Oxt.cell) %>% 
                   filter(integrated_snn_res.0.8 != 5,
                          integrated_snn_res.0.8 != 6,
                          integrated_snn_res.0.8 != 7) %>% 
                   droplevels() %>% 
                   rename(Avp.cluster.neuron = integrated_snn_res.0.8)) %>% 
        select(neuropeptides.cluster,
               Avp.cluster.neuron)) %>% 
  unclass() %>% 
  as.data.frame() %>% 
  as.matrix() %>% 
  pheatmap(color = colorRampPalette(c("white", "orange3"))(50),
           fontsize = 20,
           cex = 1,
           legend = FALSE)
dev.off()

### avp.neurons across threshold
png('avp.oxt/avp/heatmap avp.neurons across threshold.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 300)
table(left_join(burtoni.snseq.combined.sct.all.avp.neurons.4.recluster.expression %>% 
                   select(Cell.id,
                          integrated_snn_res.1) %>% 
                   rename(Avp.neuron.cluster.4 = integrated_snn_res.1),
                 burtoni.snseq.combined.sct.all.avp.neurons.recluster.expression %>% 
                   select(Cell.id,
                          integrated_snn_res.0.8) %>% 
                   rename(Avp.neuron.cluster.2 = integrated_snn_res.0.8)) %>%  
        select(Avp.neuron.cluster.4,
               Avp.neuron.cluster.2))  %>% 
  unclass() %>% 
  as.data.frame() %>% 
  as.matrix() %>% 
  pheatmap(main = 'avp.neurons across threshold',
           color = colorRampPalette(c("white", "red"))(50))
dev.off()

### all neurons vs neuropeptides all
png('avp.oxt/avp/heatmap all neurons vs avp.neurons.png',
    width = 5,
    height = 5,
    units = 'in',
    res = 720)
table(right_join(burtoni.snseq.combined.sct.all.neurons.recluster.expression %>% 
                   select(Cell.id,
                          integrated_snn_res.0.8) %>% 
                   rename(neuron.cluster = integrated_snn_res.0.8),
                 burtoni.snseq.combined.sct.all.avp.neurons.recluster.expression %>% 
                   select(Cell.id,
                          integrated_snn_res.0.8) %>% 
                   rename(Avp.cluster.neuron = integrated_snn_res.0.8)) %>% 
        select(neuron.cluster,
               Avp.cluster.neuron)) %>% 
  unclass() %>% 
  as.data.frame() %>% 
  as.matrix() %>% 
  pheatmap(color = colorRampPalette(c("white", "orange3"))(50),
           fontsize = 20,
           cex = 1,
           display_numbers = T,
           number_format = "%.0f") 
dev.off()

# for poster
### all neurons vs neuropeptides all
pdf('./neuropeptides/avp.oxt/avp/heatmap all neurons vs avp.neurons filter poster.pdf',
    width = 4.7,
    height = 4.7)
table(right_join(burtoni.snseq.combined.sct.all.neurons.recluster.expression %>% 
                   select(Cell.id,
                          integrated_snn_res.0.8) %>% 
                   rename(neuron.cluster = integrated_snn_res.0.8),
                 burtoni.snseq.combined.sct.all.avp.neurons.recluster.expression %>% 
                   mutate(Keep = as.numeric(as.character(integrated_snn_res.0.8))) %>% 
                   filter(Keep <= 4) %>% 
                   droplevels() %>% 
                   select(Cell.id,
                          integrated_snn_res.0.8) %>% 
                   rename(Avp.cluster.neuron = integrated_snn_res.0.8)) %>% 
        select(neuron.cluster,
               Avp.cluster.neuron)) %>% 
  unclass() %>% 
  as.data.frame() %>% 
  mutate(Total = rowSums(.)) %>% 
  filter(Total > 9) %>% 
  select(-c(Total)) %>% 
  as.matrix() %>% 
  pheatmap(color = colorRampPalette(c("white", "orange3"))(100),
           cex = 1,
           display_numbers = T,
           number_format = "%.0f", 
           border_color = 'black',
         treeheight_col = 0,
         treeheight_row = 0,
         angle_col = 0,
         fontsize = 15,
         legend = F) 
dev.off()  


# for poster
#update
### all neurons vs neuropeptides all
pdf('./neuropeptides/avp.oxt/avp/heatmap all neurons vs avp.neurons filter poster update.pdf',
    width = 4.7,
    height = 4.7)
table(right_join(burtoni.snseq.combined.sct.all.neurons.recluster.expression %>% 
                   select(Cell.id,
                          integrated_snn_res.0.8) %>% 
                   rename(neuron.cluster = integrated_snn_res.0.8),
                 burtoni.snseq.combined.sct.all.avp.neurons.recluster.expression %>% 
                   mutate(Keep = as.numeric(as.character(integrated_snn_res.0.8))) %>% 
                   filter(Keep <= 4) %>% 
                   droplevels() %>% 
                   select(Cell.id,
                          integrated_snn_res.0.8) %>% 
                   rename(Avp.cluster.neuron = integrated_snn_res.0.8)) %>% 
        select(neuron.cluster,
               Avp.cluster.neuron)) %>% 
  unclass() %>% 
  as.data.frame() %>% 
  mutate(Total = rowSums(.)) %>% 
  filter(Total > 9) %>% 
  select(-c(Total)) %>% 
  as.matrix() %>% 
  pheatmap(color = colorRampPalette(c("white", "orange3"))(100),
           cex = 1,
           display_numbers = T,
           number_format = "%.0f", 
           border_color = 'black',
           treeheight_col = 0,
           treeheight_row = 0,
           angle_col = 0,
           fontsize = 15,
           legend = F) 
dev.off()  









#### nucleobindin2 ####

### NUCB2 cell counts
# a paralog
NUBC2.cell.counts = data.frame(Gene.names = c("nucb2a",
                                              "nucb2b"),
                               Total.nuclei.count = c(sum(GetAssayData(object = burtoni.snseq.combined.sct.all, 
                                                                       slot = "data")["nucb2a",]>2),
                                                      sum(GetAssayData(object = burtoni.snseq.combined.sct.all, 
                                                                       slot = "data")["nucb2b",]>2)),
                               Neuron.nuclei.count = c(sum(GetAssayData(object = burtoni.snseq.combined.sct.all.neurons, 
                                                                        slot = "data")["nucb2a",]>2),
                                                       sum(GetAssayData(object = burtoni.snseq.combined.sct.all.neurons, 
                                                                        slot = "data")["nucb2b",]>2))
) %>% 
  mutate(Percentage.neuron = 100*Neuron.nuclei.count/Total.nuclei.count)

### graph
##vinplot
##dom vs sub
#nucb2a
VlnPlot(burtoni.snseq.combined.sct.all, 
        features = c("nucb2a"),
        group.by = "Cell.type",
        split.by = "orig.ident") +
  ylim(2,15)
ggsave('nucb2/vinplot.nucb2a.cell.types.and.status.png',
       width = 10,
       height = 10)

#nucb2b
VlnPlot(burtoni.snseq.combined.sct.all, 
        features = c("nucb2b"),
        group.by = "Cell.type",
        split.by = "orig.ident") +
  ylim(2,15)
ggsave('nucb2/vinplot.nucb2b.cell.types.and.status.png',
       width = 10,
       height = 10)


### nucb2b neurons
### subset data to relevant cells
burtoni.snseq.combined.sct.all.nucb2.neurons = burtoni.snseq.combined.sct.all

#set idents
Idents(object = burtoni.snseq.combined.sct.all.nucb2.neurons) <- "Cell.type"

#subset
burtoni.snseq.combined.sct.all.nucb2.neurons = subset(burtoni.snseq.combined.sct.all.nucb2.neurons,
                                                      subset = nucb2b > 2 | nucb2a > 2, idents = c("neurons",
                                                                                                   "excitatory",
                                                                                                   "inhibitory"))


## run PCA, UMAP, and cluster 
#use 0.4 resolution
burtoni.snseq.combined.sct.all.nucb2.neurons.recluster = burtoni.snseq.combined.sct.all.nucb2.neurons %>% 
  RunPCA() %>%
  FindNeighbors(dims = 1:15) %>%
  RunUMAP(dims = 1:15) %>%
  FindClusters(resolution = 0.2)

### clustree
# cluster across resolutions
burtoni.snseq.combined.sct.all.nucb2.neurons.clustree <- Seurat::FindClusters(object = burtoni.snseq.combined.sct.all.nucb2.neurons.recluster, 
                                                                              resolution = resolution.range.reduced)
#check data
head(burtoni.snseq.combined.sct.all.nucb2.neurons.clustree[[]])
#clustree
clustree(burtoni.snseq.combined.sct.all.nucb2.neurons.clustree, 
         prefix = "integrated_snn_res.",
         node_colour = 'cluster',
         node_size_range = c(10,20),
         scale_node_text = TRUE) +
  scale_edge_color_continuous(low = "black", 
                              high = "black") +
  theme(legend.position = "bottom")
ggsave('nucb2/nucb2.neurons.clustree.png',
       width = 10,
       height = 10)


### graph 
##dom vs sub
DimPlot(burtoni.snseq.combined.sct.all.nucb2.neurons.recluster, 
        reduction = "umap", 
        group.by = "orig.ident") +
  ggtitle('Dom vs Sub nucb2 neurons UMAP')
ggsave('nucb2/DomVsSub.dimplot.nucb2.neurons.all.png')

## clusters
DimPlot(burtoni.snseq.combined.sct.all.nucb2.neurons.recluster, 
        reduction = "umap", 
        label = TRUE,
        repel = TRUE) +
  ggtitle('nucb2 neurons UMAP clusters')
ggsave('nucb2/Clusters.dimplot.nucb2.neurons.all.png',
       width = 10,
       height = 10)


###expression level
burtoni.snseq.combined.sct.all.nucb2.neurons.recluster.expression = full_join(full_join(burtoni.snseq.combined.sct.all.nucb2.neurons.recluster@reductions$umap@cell.embeddings %>% 
                                                                                          as.data.frame() %>% 
                                                                                          rownames_to_column("Cell.id"),
                                                                                        burtoni.snseq.combined.sct.all.nucb2.neurons.recluster@meta.data %>% 
                                                                                          rownames_to_column("Cell.id")),
                                                                              burtoni.snseq.combined.sct.all.nucb2.neurons.recluster@assays$integrated@scale.data %>% 
                                                                                as.data.frame() %>% 
                                                                                filter(rownames(burtoni.snseq.combined.sct.all.nucb2.neurons.recluster@assays$integrated@scale.data) %in% c('avp',
                                                                                                                                                                                            'oxt', 
                                                                                                                                                                                            'nucb2a',
                                                                                                                                                                                            'nucb2b')) %>% 
                                                                                t() %>% as.data.frame() %>% 
                                                                                rownames_to_column('Cell.id')) %>% 
  full_join(burtoni.snseq.combined.sct.all.nucb2.neurons.recluster@reductions$pca@cell.embeddings %>% 
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
                                  'avp')))%>% 
  mutate(Nucb2.cell = ifelse(nucb2a > 2 & nucb2b > 2,
                           'both',
                           ifelse(nucb2a > 2,
                                  'nucb2a',
                                  'nucb2b')))

#no seperation by cell type
table(burtoni.snseq.combined.sct.all.nucb2.neurons.recluster.expression %>% 
        select(Cell.type,
               integrated_snn_res.0.2))

# difference social status
table(burtoni.snseq.combined.sct.all.nucb2.neurons.recluster.expression %>% 
        select(integrated_snn_res.0.2,
               orig.ident))

# difference social status
table(burtoni.snseq.combined.sct.all.nucb2.neurons.recluster.expression %>% 
        select(orig.ident,
               Oxt.cell))


# difference social status
table(burtoni.snseq.combined.sct.all.nucb2.neurons.recluster.expression %>% 
        select(orig.ident,
               Oxt.cell,
               Nucb2.cell))

## graph nucb2 
burtoni.snseq.combined.sct.all.nucb2.neurons.recluster.expression %>% 
  mutate(Gene.expressed = ifelse(nucb2a >= 3,
                                 'nucb2a',
                                 ifelse(nucb2b >= 3,
                                        'nucb2b',
                                        NA))) %>% 
  droplevels() %>% 
  ggplot(aes(x = UMAP_1,
             y = UMAP_2,
             color = Gene.expressed)) +
  geom_point() +
  theme_classic() +
  ggtitle('nucb2a vs nucb2b neurons UMAP')
ggsave('nucb2/UMAP paralogs across neuron clusters.png',
       width = 10,
       height = 10)





#### neuropeptides comparison ####
### create matrix of neuropeptide expression level 
## use neuropeptides.list
neuropeptides.list.gene.names =  neuropeptides.list %>% 
  filter(!is.na(Gene.name.nile.tilapia)) %>%
  pull(Gene.name.nile.tilapia) %>% 
  unique()
   
## use all neurons expression recluster
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides = full_join(full_join(burtoni.snseq.combined.sct.all.neurons.recluster@reductions$umap@cell.embeddings %>% 
                                                                                    as.data.frame() %>% 
                                                                                    rownames_to_column("Cell.id"),
                                                                                  burtoni.snseq.combined.sct.all.neurons.recluster@meta.data %>% 
                                                                                    rownames_to_column("Cell.id")),
                                                                        burtoni.snseq.combined.sct.all.neurons.recluster@assays$integrated@scale.data %>% 
                                                                          as.data.frame() %>% 
                                                                          filter(rownames(burtoni.snseq.combined.sct.all.neurons.recluster@assays$integrated@scale.data) %in% neuropeptides.list.gene.names) %>% t() %>% as.data.frame() %>% 
                                                                          rownames_to_column('Cell.id')) %>% 
  full_join(burtoni.snseq.combined.sct.all.neurons.recluster@reductions$pca@cell.embeddings %>% 
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
                                  ifelse(avp > 2,
                                         'avp',
                                         NA))))

## reduce columns to neuropeptides
# create presence absence matrix
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix = burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides %>% 
  select(any_of(neuropeptides.list.gene.names)) %>% 
  mutate(across(where(is.numeric), 
                function(x) ifelse(x < 2, 0, x))) %>% 
  mutate(across(where(is.numeric), 
                function(x) ifelse(x >= 2, 1, x))) 

###rename to actual gene name
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix = burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix %>% 
  rename(#ADM = ENSONIG00000036130,
         CHGA = ENSONIG00000000725,
         GHRH = ENSONIG00000020252,
         GNRH1 = ENSONIG00000011023,
         #KNG1 = ENSONIG00000008378,
         NMB = ENSONIG00000002537,
         NPFF = ENSONIG00000005229,
         SST = ENSONIG00000033642)


#how many neurons have neuropeptides
#total 14592
#~85% have at least one
#~56% have at least two

#graph number of neuropeptides per cell
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix %>%
  mutate(Total.neuropeptides = rowSums(.)) %>% 
  ggplot(aes(Total.neuropeptides)) +
  geom_histogram(binwidth = 0.3) +
  theme_classic() +
  scale_x_continuous(breaks=seq(0,11,1)) +
  ggtitle('Neuropeptides per neuron')
ggsave('neuropeptides/comparison/Neuropeptides per neuron.png',
       width = 10,
       height = 10)

# for poster
#graph number of neuropeptides per cell
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix %>%
  mutate(Total.neuropeptides = rowSums(.)) %>% 
  ggplot(aes(Total.neuropeptides)) +
  geom_histogram(binwidth = 0.3) +
  theme_classic() +
  scale_x_continuous(breaks=seq(0,11,1)) +
  xlab('Neuropeptides per neuron') +
  ylab('Number of neurons') +
  ggtitle('Neuropeptides per neuron') +
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1))+
  theme(axis.text = element_text(size = 15))  +
  theme(axis.title = element_text(size = 20))+
  theme(plot.title = element_text(size=20))
ggsave('neuropeptides/comparison/Neuropeptides per neuron.pdf',
       width = 10,
       height = 5)


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

##graph total count
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.count %>%
  ggplot(aes(Total.count)) +
  geom_histogram(binwidth = 100) +
  theme_classic()
ggsave('neuropeptides/comparison/Neuropeptides neuron counts.png',
       width = 10,
       height = 10)

# calculate the percentage
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent = t(burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap * 100 / diag(burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap))     
#convert data to long format
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent = burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent %>% 
  as.data.frame.table() 

##graph
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent %>% 
  filter(Freq < 100) %>% 
  ggplot(aes(x=Var1, 
             y=Var2)) +
  geom_tile(aes(fill=Freq)) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradientn(colours=c("blue","red")) 
ggsave('neuropeptides/comparison/Neuropeptides percentage overlap.png',
       width = 10,
       height = 10)

#add gene count data
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder = left_join(burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent,
                                                                                                             burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.count,
                                                                                                             by = c('Var1' = 'Gene.name')) %>% 
  rename(Var1.total.count = Total.count) %>% 
  left_join(burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.count,
            by = c('Var2' = 'Gene.name')) %>% 
  rename(Var2.total.count = Total.count)

#reorder 
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder <- burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder %>%
  arrange(Var1.total.count,
          Var2.total.count) %>%
  mutate(Var1 = fct_inorder(factor(Var1, ordered=TRUE)),
         Var2 = fct_inorder(factor(Var2, ordered=TRUE)))

#graph reorder
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder %>% 
  filter(Freq < 100) %>% 
  ggplot(aes(x=Var1, 
             y=Var2)) +
  geom_tile(aes(fill=Freq)) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradientn(colours=c("blue","red")) 
ggsave('neuropeptides/comparison/Neuropeptides percentage overlap reorder.png',
       width = 10,
       height = 10)

#graph reorder
#filter cell count
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder %>% 
  filter(Freq < 100) %>% 
  filter(Var2.total.count > 250) %>% 
  filter(Var1.total.count > 250) %>% 
  ggplot(aes(x=Var1, 
             y=Var2)) +
  geom_tile(aes(fill=Freq)) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradientn(colours=c("blue","red")) 
ggsave('neuropeptides/comparison/Neuropeptides percentage overlap reorder filter.png',
       width = 10,
       height = 10)



# calculate the percentage
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.df = burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap %>% 
  as.data.frame.table() 

##graph
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.df %>% 
  ggplot(aes(x=Var1, 
             y=Var2)) +
  geom_tile(aes(fill=Freq)) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradientn(colours=c("blue","red")) 
ggsave('neuropeptides/comparison/Neuropeptides overlap.png',
       width = 10,
       height = 10)

#add gene count data
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.df.reorder = left_join(burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.df,
                                                                                                                     burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.count,
                                                                                                                     by = c('Var1' = 'Gene.name')) %>% 
  rename(Var1.total.count = Total.count) %>% 
  left_join(burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.count,
            by = c('Var2' = 'Gene.name')) %>% 
  rename(Var2.total.count = Total.count)

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

#graph reorder
burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.df.reorder %>% 
  ggplot(aes(x=Var1, 
             y=Var2)) +
  geom_tile(aes(fill=Freq.na)) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradientn(colours=c("blue","red")) 
ggsave('neuropeptides/comparison/Neuropeptides overlap rescaled.png',
       width = 10,
       height = 10)





                                                                        
#### neuropeptides comparison across social status ####
### create matrix of neuropeptide expression level 
## use neuropeptides.list
## reduce columns to neuropeptides
## seperate into dom and sub
# create presence absence matrix
#dom
dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix = burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides %>% 
  filter(orig.ident == 'dom_burtoni_snseq') %>% 
  select(any_of(neuropeptides.list.gene.names)) %>% 
  mutate(across(where(is.numeric), 
                function(x) ifelse(x < 2, 0, x))) %>% 
  mutate(across(where(is.numeric), 
                function(x) ifelse(x >= 2, 1, x))) 
#sub
sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix = burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides %>% 
  filter(orig.ident == 'sub_burtoni_snseq') %>% 
  select(any_of(neuropeptides.list.gene.names)) %>% 
  mutate(across(where(is.numeric), 
                function(x) ifelse(x < 2, 0, x))) %>% 
  mutate(across(where(is.numeric), 
                function(x) ifelse(x >= 2, 1, x))) 

##rename to actual gene name
#dom
dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix = dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix %>% 
  rename(#ADM = ENSONIG00000036130,
    CHGA = ENSONIG00000000725,
    GHRH = ENSONIG00000020252,
    GNRH1 = ENSONIG00000011023,
    #KNG1 = ENSONIG00000008378,
    NMB = ENSONIG00000002537,
    NPFF = ENSONIG00000005229,
    SST = ENSONIG00000033642)

#sub
sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix = sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix %>% 
  rename(#ADM = ENSONIG00000036130,
    CHGA = ENSONIG00000000725,
    GHRH = ENSONIG00000020252,
    GNRH1 = ENSONIG00000011023,
    #KNG1 = ENSONIG00000008378,
    NMB = ENSONIG00000002537,
    NPFF = ENSONIG00000005229,
    SST = ENSONIG00000033642)

#how many neurons have neuropeptides
#dom = 8181
#sub = 6411

#graph number of neuropeptides per cell
#dom
dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix %>%
  mutate(Total.neuropeptides = rowSums(.)) %>% 
  ggplot(aes(Total.neuropeptides)) +
  geom_histogram(binwidth = 0.3) +
  theme_classic() +
  scale_x_continuous(breaks=seq(0,11,1)) +
  ggtitle('Neuropeptides per neuron dom')
ggsave('neuropeptides/comparison/social.status/Neuropeptides per neuron dom.png',
       width = 10,
       height = 10)

#sub
sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix %>%
  mutate(Total.neuropeptides = rowSums(.)) %>% 
  ggplot(aes(Total.neuropeptides)) +
  geom_histogram(binwidth = 0.3) +
  theme_classic() +
  scale_x_continuous(breaks=seq(0,11,1)) +
  ggtitle('Neuropeptides per neuron sub')
ggsave('neuropeptides/comparison/social.status/Neuropeptides per neuron sub.png',
       width = 10,
       height = 10)



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

##graph total count
#dom
dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.count %>%
  ggplot(aes(Total.count)) +
  geom_histogram(binwidth = 100) +
  theme_classic() +
  ggtitle('Neuropeptides neuron counts dom')
ggsave('neuropeptides/comparison/social.status/Neuropeptides neuron counts dom.png',
       width = 10,
       height = 10)
#sub
sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.count %>%
  ggplot(aes(Total.count)) +
  geom_histogram(binwidth = 100) +
  theme_classic() +
  ggtitle('Neuropeptides neuron counts sub')
ggsave('neuropeptides/comparison/social.status/Neuropeptides neuron counts sub.png',
       width = 10,
       height = 10)

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

##graph
#dom
dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent %>% 
  filter(Freq < 100) %>% 
  ggplot(aes(x=Var1, 
             y=Var2)) +
  geom_tile(aes(fill=Freq)) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradientn(colours=c("blue","red")) +
  ggtitle('Neuropeptides percentage overlap dom')
ggsave('neuropeptides/comparison/social.status/Neuropeptides percentage overlap dom.png',
       width = 10,
       height = 10)
#sub
sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent %>% 
  filter(Freq < 100) %>% 
  ggplot(aes(x=Var1, 
             y=Var2)) +
  geom_tile(aes(fill=Freq)) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradientn(colours=c("blue","red")) +
  ggtitle('Neuropeptides percentage overlap sub')
ggsave('neuropeptides/comparison/social.status/Neuropeptides percentage overlap sub.png',
       width = 10,
       height = 10)

#add gene count data
#dom
dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder = left_join(dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent,
                                                                                                                     dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.count,
                                                                                                                     by = c('Var1' = 'Gene.name')) %>% 
  rename(Var1.total.count = Total.count) %>% 
  left_join(dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.count,
            by = c('Var2' = 'Gene.name')) %>% 
  rename(Var2.total.count = Total.count)

#sub
sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder = left_join(sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent,
                                                                                                                         sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.count,
                                                                                                                         by = c('Var1' = 'Gene.name')) %>% 
  rename(Var1.total.count = Total.count) %>% 
  left_join(sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.count,
            by = c('Var2' = 'Gene.name')) %>% 
  rename(Var2.total.count = Total.count)

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

#graph reorder
#dom
dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder %>% 
  filter(Freq < 100) %>% 
  ggplot(aes(x=Var1, 
             y=Var2)) +
  geom_tile(aes(fill=Freq)) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradientn(colours=c("blue","red")) +
  ggtitle('Neuropeptides percentage overlap reorder dom')
ggsave('neuropeptides/comparison/social.status/Neuropeptides percentage overlap reorder dom.png',
       width = 10,
       height = 10)
#sub
sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder %>% 
  filter(Freq < 100) %>% 
  ggplot(aes(x=Var1, 
             y=Var2)) +
  geom_tile(aes(fill=Freq)) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradientn(colours=c("blue","red")) +
  ggtitle('Neuropeptides percentage overlap reorder sub')
ggsave('neuropeptides/comparison/social.status/Neuropeptides percentage overlap reorder sub.png',
       width = 10,
       height = 10)

#graph reorder
#filter cell count
#dom
dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder %>% 
  filter(Freq < 100) %>% 
  filter(Var2.total.count > 250) %>% 
  filter(Var1.total.count > 250) %>% 
  ggplot(aes(x=Var1, 
             y=Var2)) +
  geom_tile(aes(fill=Freq)) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradientn(colours=c("blue","red")) +
  ggtitle('Neuropeptides percentage overlap reorder filter dom')
ggsave('neuropeptides/comparison/social.status/Neuropeptides percentage overlap reorder filter dom.png',
       width = 10,
       height = 10)
#sub
sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder %>% 
  filter(Freq < 100) %>% 
  filter(Var2.total.count > 250) %>% 
  filter(Var1.total.count > 250) %>% 
  ggplot(aes(x=Var1, 
             y=Var2)) +
  geom_tile(aes(fill=Freq)) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradientn(colours=c("blue","red")) +
  ggtitle('Neuropeptides percentage overlap reorder filter sub')
ggsave('neuropeptides/comparison/social.status/Neuropeptides percentage overlap reorder filter sub.png',
       width = 10,
       height = 10)



### calculate  count
#dom
dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.df = dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap %>% 
  as.data.frame.table() 
#sub
sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.df = sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap %>% 
  as.data.frame.table() 

##graph
#dom
dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.df %>% 
  ggplot(aes(x=Var1, 
             y=Var2)) +
  geom_tile(aes(fill=Freq)) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradientn(colours=c("blue","red")) +
  ggtitle('Neuropeptides overlap dom')
ggsave('neuropeptides/comparison/social.status/Neuropeptides overlap.png',
       width = 10,
       height = 10)
#sub
sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.df %>% 
  ggplot(aes(x=Var1, 
             y=Var2)) +
  geom_tile(aes(fill=Freq)) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradientn(colours=c("blue","red")) +
  ggtitle('Neuropeptides overlap sub')
ggsave('neuropeptides/comparison/social.status/Neuropeptides overlap.png',
       width = 10,
       height = 10)

#add gene count data
#dom
dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.df.reorder = left_join(dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.df,
                                                                                                                    dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.count,
                                                                                                                by = c('Var1' = 'Gene.name')) %>% 
  rename(Var1.total.count = Total.count) %>% 
  left_join(dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.count,
            by = c('Var2' = 'Gene.name')) %>% 
  rename(Var2.total.count = Total.count)
#sub
sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.df.reorder = left_join(sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.df,
                                                                                                                    sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.count,
                                                                                                                    by = c('Var1' = 'Gene.name')) %>% 
  rename(Var1.total.count = Total.count) %>% 
  left_join(sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.count,
            by = c('Var2' = 'Gene.name')) %>% 
  rename(Var2.total.count = Total.count)

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

#graph reorder
#dom
dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.df.reorder %>% 
  ggplot(aes(x=Var1, 
             y=Var2)) +
  geom_tile(aes(fill=Freq.na)) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradientn(colours=c("blue","red")) +
  ggtitle('Neuropeptides overlap rescaled dom')
ggsave('neuropeptides/comparison/social.status/Neuropeptides overlap rescaled dom.png',
       width = 10,
       height = 10)
#sub
sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.df.reorder %>% 
  ggplot(aes(x=Var1, 
             y=Var2)) +
  geom_tile(aes(fill=Freq.na)) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradientn(colours=c("blue","red")) +
  ggtitle('Neuropeptides overlap rescaled sub')
ggsave('neuropeptides/comparison/social.status/Neuropeptides overlap rescaled sub.png',
       width = 10,
       height = 10)

### compare percentage overlap
## between doms and subs
#status
#set Var1 and Var2 to characters
#rename variables
status.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder = full_join(dom.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder %>% 
                                                                                                                              mutate(Var1 = as.character(Var1),
                                                                                                                                     Var2 = as.character(Var2)) %>% 
                                                                                                                              rename(Freq.dom = Freq,
                                                                                                                                     Var1.total.count.dom = Var1.total.count,
                                                                                                                                     Var2.total.count.dom = Var2.total.count),
                                                                                                                            sub.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder%>% 
                                                                                                                              mutate(Var1 = as.character(Var1),
                                                                                                                                     Var2 = as.character(Var2))%>% 
                                                                                                                              rename(Freq.sub = Freq,
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

#graph reorder
#status
status.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder %>% 
  ggplot(aes(x=Var1, 
             y=Var2)) +
  geom_tile(aes(fill=Freq.diff)) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradientn(colours=c('dark blue', "blue",'white',"red", 'dark red')) +
  ggtitle('Neuropeptides percentage overlap reorder difference')
ggsave('neuropeptides/comparison/social.status/Neuropeptides percentage difference overlap reorder status.png',
       width = 10,
       height = 10)

#graph reorder
#status
#filter
status.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder %>% 
  filter(Var2.total.count.dom > 200) %>% 
  filter(Var2.total.count.sub > 200) %>% view()
  ggplot(aes(x=Var1, 
             y=Var2)) +
  geom_tile(aes(fill=Freq.diff)) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradientn(colours=c('dark blue', "blue",'white',"red", 'dark red'), limits = c(-15, 15)) +
  ggtitle('Neuropeptides percentage overlap reorder difference filter')
ggsave('neuropeptides/comparison/social.status/Neuropeptides percentage difference overlap reorder status filter.png',
       width = 10,
       height = 10)

status.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder %>% 
  filter(Var2.total.count.dom > 200 & Var1.total.count.dom > 200) %>% 
  filter(Var2.total.count.sub > 200 & Var1.total.count.sub >200) %>%
ggplot(aes(x=Var1, 
           y=Var2)) +
  geom_tile(aes(fill=Freq.diff)) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradientn(colours=c('dark blue', "blue",'white',"red", 'dark red'), limits = c(-10, 10)) +
  ggtitle('Neuropeptides percentage overlap reorder difference filter')


##graph nueropeptide counts
#status

#count neuropeptides across social status
status.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder %>% select(Var2, Var2.total.count.dom, Var2.total.count.sub) %>% distinct() %>% select(Var2.total.count.dom, Var2.total.count.sub) %>% droplevels %>% colSums()

#Var2.total.count.dom=16527, Var2.total.count.sub =11248
#ratio = 0.6806
status.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder %>% 
  mutate(Ratio = Var2.total.count.dom/Var2.total.count.sub,
         Keep = ifelse(Ratio >= 2 | Ratio <= 0.66,
                       as.character(Var2),
                       NA),
         Keep =  case_when(Var2 == "avp"~ "avp",
                           Var2 == "oxt"~ "oxt",
                           Var2 == "SST"~ "SST",
                           Var2 == "GNRH1"~ "GNRH1",
           TRUE ~ Keep)) %>% 
  ggplot(aes(x=Var2.total.count.dom, 
             y=Var2.total.count.sub,
             label = Keep)) +
  geom_abline(slope = 1, 
              intercept = 0) +
  geom_abline(slope = 1.3194, 
              intercept = 0, 
              linetype = 'dashed') +
  geom_abline(slope = 0.6806,
              intercept = 0, 
              linetype = 'dashed') +
  geom_smooth(method = 'lm') +
  geom_point() +
  geom_text(hjust = 0, 
            nudge_x = 5,
            check_overlap = TRUE)+
  theme_classic() +
  lims(x = c(0,1000),
       y = c(0,1000)) +
  ggtitle('Neuropeptides count comparison across social status')
ggsave('neuropeptides/comparison/social.status/Neuropeptides count comparison across social status.png',
       width = 10,
       height = 10)
  


#Var2.total.count.dom=8181, Var2.total.count.sub =6411
#ratio = 0.783645
status.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder %>% 
  mutate(Dom.count.scaled_0.783645 = Var2.total.count.dom*0.783645,
         Sub.count = Var2.total.count.sub, 
         Ratio.scaled = Dom.count.scaled_0.783645/Sub.count,
         Keep = case_when(Ratio.scaled >= 1.5 & Var2.total.count.dom >= 100 ~ as.character(Var2),
                          Ratio.scaled <= 0.667 & Var2.total.count.sub >= 100 ~ as.character(Var2),
                          TRUE ~ as.character(NA)),
         Keep.name =  case_when(Var2 == "avp"~ "avp",
                                Var2 == "oxt"~ "oxt",
                                Var2 == "SST"~ "SST",
                                Var2 == "GNRH1"~ "GNRH1",
                                TRUE ~ Keep),
         Keep.color = case_when(Ratio.scaled >= 1.5 & Var2.total.count.dom >=100 ~ "Dominant bias",
                                Ratio.scaled <= 0.667 & Var2.total.count.sub >=100 ~ "Subordinate bias",
                                TRUE ~ "No bias")) %>% 
  select(c(Dom.count.scaled_0.783645,
           Sub.count,
           Keep.name,
           Keep.color)) %>% 
  distinct() %>% 
  ggplot(aes(x=Dom.count.scaled_0.783645, 
             y=Sub.count,
             label = Keep.name)) +
  geom_abline(slope = 1, 
              intercept = 0) +
  geom_abline(slope = 0.667, 
              intercept = 0, 
              linetype = 'dashed') +
  geom_abline(slope = 1.5,
              intercept = 0, 
              linetype = 'dashed') +
  geom_point(aes(fill = Keep.color),
             shape = 21,
             size = 4) +
  geom_label_repel()+
  theme_classic() +
  lims(x = c(0,750),
       y = c(0,750)) +
  ggtitle('Neuropeptide neuron count across social status') +
  xlab("Dominant neuron count scaled (78.36%)") +
  ylab("Subordinate neuron count") +
  scale_fill_manual(breaks = c("Dominant bias",
                               "Subordinate bias",
                               "No bias"),
                    values = c("#4e499e",
                               "#60bb46",
                               "black")) +
  theme(legend.position = c(0.8, 
                            0.25),
        legend.box.background = element_rect(colour = "black"),
        legend.background = element_blank(),
        legend.text=element_text(size=10),
        legend.title = element_blank())+
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 10),
        title =element_text(size=10, face='bold'))
ggsave('neuropeptides/comparison/social.status/Neuropeptides count comparison across social status scaled.png',
       width = 10,
       height = 10)

# for poster
#Var2.total.count.dom=8181, Var2.total.count.sub =6411
#ratio = 0.783645
status.burtoni.snseq.combined.sct.all.neurons.recluster.expression.neuropeptides.matrix.overlap.percent.reorder %>%
  mutate(Dom.count.scaled_0.783645 = Var2.total.count.dom*0.783645,
         Sub.count = Var2.total.count.sub,
         Ratio.scaled = Dom.count.scaled_0.783645/Sub.count,
         Keep = case_when(Ratio.scaled >= 1.5 & Var2.total.count.dom >= 100 ~ as.character(Var2),
                          Ratio.scaled <= 0.667 & Var2.total.count.sub >= 100 ~ as.character(Var2),
                          TRUE ~ as.character(NA)),
         Keep.name =  case_when(Var2 == "avp"~ "avp",
                                Var2 == "oxt"~ "oxt",
                                Var2 == "SST"~ "SST",
                                Var2 == "GNRH1"~ "GNRH1",
                                TRUE ~ Keep),
         Keep.color = case_when(Ratio.scaled >= 1.5 & Var2.total.count.dom >=100 ~ "Dominant bias",
                                Ratio.scaled <= 0.667 & Var2.total.count.sub >=100 ~ "Subordinate bias",
                                TRUE ~ "No bias")) %>%
  select(c(Dom.count.scaled_0.783645,
           Sub.count,
           Keep.name,
           Keep.color)) %>%
  distinct() %>%
  ggplot(aes(x=Dom.count.scaled_0.783645,
             y=Sub.count,
             label = Keep.name)) +
  geom_abline(slope = 1,
              intercept = 0) +
  geom_abline(slope = 0.667,
              intercept = 0,
              linetype = 'dashed') +
  geom_abline(slope = 1.5,
              intercept = 0,
              linetype = 'dashed') +
  geom_point(aes(fill = Keep.color),
             shape = 21,
             size = 6) +
  geom_label_repel(size = 8)+
  theme_classic() +
  lims(x = c(0,750),
       y = c(0,750)) +
  ggtitle('Neuropeptide neuron count across social status') +
  xlab("Dominant neuron count scaled (78.36%)") +
  ylab("Subordinate neuron count") +
  scale_fill_manual(breaks = c("Dominant bias",
                               "Subordinate bias",
                               "No bias"),
                    values = c("#4e499e",
                               "#60bb46",
                               "black")) +
  theme(legend.position = c(0.8,
                            0.25),
        legend.box.background = element_rect(colour = "black"),
        legend.background = element_blank(),
        legend.text=element_text(size=15),
        legend.title = element_blank())+
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1))+
  theme(axis.text = element_text(size = 15))  +
  theme(axis.title = element_text(size = 20))+
  theme(plot.title = element_text(size=20))
ggsave('neuropeptides/comparison/social.status/Neuropeptides count comparison across social status scaled poster.pdf',
       width = 10,
       height = 7.5)
