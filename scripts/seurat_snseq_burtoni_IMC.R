#### Burtoni snseq seurat analysis
### comparing dom and sub snseq data
###Note: Seurat requires R version > 4
## use lambcomp1 to run R with command 
# > R-4.0.3
## https://satijalab.org/seurat/articles/pbmc3k_tutorial.html

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
## update SCTtransform
# devtools::install_github("satijalab/seurat", ref = "develop")
# BiocManager::install("glmGamPoi")

#load libraries
library(tidyverse)
library(Seurat)
library(patchwork)
library(clustree)

#### load data ####
### Load data
## dom data
dom.data <- Read10X(data.dir = "../cellranger/count/NileTilapia.reference/dom/Burtoni-dom/outs/filtered_feature_bc_matrix/")

## sub data
sub.data <- Read10X(data.dir = "../cellranger/count/NileTilapia.reference/sub/Burtoni-sub/outs/filtered_feature_bc_matrix/")

### load marker genes
marker.genes = read.csv('../Gene.lists/Cell.type.gene.markers.csv')
## extract marker genes
marker.genes.list = marker.genes %>% 
  select(gene.id.tilapia,
         cell.type,
         cell.subtype.category,
         cell.subtype) %>% 
  filter(!is.na(gene.id.tilapia)) %>% 
  distinct()

## make list of broad cell types
cell.types = marker.genes.list %>% 
  filter(cell.subtype.category == 'broad') %>% 
  pull(cell.subtype) %>% 
  unique()

#### Clustering ####
###filter cells first
### try new normalization method
###dom
##create object
dom.object <- CreateSeuratObject(counts = dom.data,
                                 project = "dom_burtoni_snseq", 
                                 min.cells = 3,
                                 min.features = 200)
## label mitochondria genes
dom.object <- PercentageFeatureSet(dom.object, 
                                                      pattern = "^MT", 
                                                      col.name = "percent.mt")

## QC graph
VlnPlot(dom.object, 
        features = c("nFeature_RNA", 
                     "nCount_RNA",
                     "percent.mt"), 
        ncol = 3)
ggsave('dom/QCmetrics.dom.png')

### Scatterplot 
##features and counts
#dom
FeatureScatter(dom.object,
               feature1 = "nCount_RNA",
               feature2 = "nFeature_RNA") 
ggsave('dom/featurescatter.dom.png')
#filter
FeatureScatter(subset(dom.object, 
                      subset = nCount_RNA < 6000),
               feature1 = "nCount_RNA",
               feature2 = "nFeature_RNA")+
  geom_hline(yintercept = 1000)+
  geom_hline(yintercept = 2000)
ggsave('dom/featurescatter.dom.limits.png')
#histogram
dom.object@meta.data %>% 
  subset(nFeature_RNA < 3000) %>%
  mutate(Kept = ifelse(nFeature_RNA > 1000 & nFeature_RNA < 2000,
                       'Kept',
                       'Remove')) %>% 
  ggplot(aes(nFeature_RNA,
             fill = Kept))+
  geom_histogram(binwidth = 10) +
  theme_bw() +
  ggtitle('dom_burtoni_snseq UMI counts')  +
  scale_fill_manual(values = c('red',
                                 'black'))
ggsave('dom/dom_burtoni_snseq UMI counts filtered.png')


#filter
# FeatureScatter(subset(dom.object, subset = nFeature_RNA > 1000 & nFeature_RNA < 2000),
#                feature1 = "nCount_RNA",
#                feature2 = "nFeature_RNA")

#sub
FeatureScatter(sub.object,
               feature1 = "nCount_RNA",
               feature2 = "nFeature_RNA")
ggsave('sub/featurescatter.sub.png')
#filter
FeatureScatter(subset(sub.object, 
                      subset = nCount_RNA < 6000),
               feature1 = "nCount_RNA",
               feature2 = "nFeature_RNA")+
  geom_hline(yintercept = 1000)+
  geom_hline(yintercept = 2000)
ggsave('sub/featurescatter.sub.limits.png')
#histogram
sub.object@meta.data %>% 
  subset(nFeature_RNA < 3000) %>%
  mutate(Kept = ifelse(nFeature_RNA > 1000 & nFeature_RNA < 2000,
                       'Kept',
                       'Remove')) %>% 
  ggplot(aes(nFeature_RNA,
             fill = Kept))+
  geom_histogram(binwidth = 10) +
  theme_bw() +
  ggtitle('sub_burtoni_snseq UMI counts')  +
  scale_fill_manual(values = c('red',
                               'black'))
ggsave('sub/sub_burtoni_snseq UMI counts filtered.png')

## check filters
#count 15392 nuclei 
dom.object@meta.data %>% 
  subset(nFeature_RNA <= 2000) %>%
  subset(nFeature_RNA >= 1000) %>% 
  summarize(max = max (nFeature_RNA),
            min = min (nFeature_RNA),
            count = n())

## filter data
dom.object.SCTransform <- subset(dom.object, subset = nFeature_RNA > 1000 & nFeature_RNA < 2000)
#check filters worked
dom.object.SCTransform@meta.data %>% 
  summarize(max = max (nFeature_RNA),
            min = min (nFeature_RNA),
            count = n())

## run SCTransform, PCA, UMAP, and cluster 
#use 0.8 resolution from clustree
dom.object.SCTransform = dom.object.SCTransform %>% 
  SCTransform(vst.flavor = "v2") %>% #vars.to.regress = "percent.mt"?
  RunPCA() %>%
  FindNeighbors(dims = 1:15) %>%
  RunUMAP(dims = 1:15) %>%
  FindClusters(resolution = 0.8)

## check PC's in elbow plot
ElbowPlot(dom.object.SCTransform)
ggsave('dom/Elbowplot.dom.png')

###sub
##create object
sub.object <- CreateSeuratObject(counts = sub.data,
                                 project = "sub_burtoni_snseq", 
                                 min.cells = 3,
                                 min.features = 200)
## label mitochondria genes
sub.object <- PercentageFeatureSet(sub.object, 
                                   pattern = "^MT", 
                                   col.name = "percent.mt")

## QC graph
VlnPlot(sub.object, 
        features = c("nFeature_RNA", 
                     "nCount_RNA",
                     "percent.mt"), 
        ncol = 3)
ggsave('sub/QCmetrics.sub.png')

## check filters
#count 12514 nuclei 
sub.object@meta.data %>% 
  subset(nFeature_RNA <= 2000) %>%
  subset(nFeature_RNA >= 1000) %>% 
  summarize(max = max (nFeature_RNA),
            min = min (nFeature_RNA),
            count = n())

## filter data
sub.object.SCTransform <- subset(sub.object, subset = nFeature_RNA > 1000 & nFeature_RNA < 2000)
#check filters worked
sub.object.SCTransform@meta.data %>% 
  summarize(max = max (nFeature_RNA),
            min = min (nFeature_RNA),
            count = n())

## run SCTransform, PCA, UMAP, and cluster 
#use 0.8 resolution from clustree
sub.object.SCTransform = sub.object.SCTransform %>% 
  SCTransform(vst.flavor = "v2") %>% #vars.to.regress = "percent.mt"?
  RunPCA() %>%
  FindNeighbors(dims = 1:15) %>%
  RunUMAP(dims = 1:15) %>%
  FindClusters(resolution = 0.8)

## check PC's in elbow plot
ElbowPlot(sub.object.SCTransform)
ggsave('sub/Elbowplot.sub.png')

#### Clustree ####
# Select a range of resolutions
resolution.range <- seq(from = 0, to = 2, by = 0.2)
#reduced range
resolution.range.reduced <- seq(from = 0, to = 1.2, by = 0.2)
#reduced range
resolution.range.reduced.2 <- seq(from = 0, to = 1, by = 0.2)

### Find clusters using a range of resolutions
##dom
#create dummy dataset
dom.object.SCTransform.filter.2 = dom.object.SCTransform
# cluster across range
dom.object.SCTransform.filter.clustree <- Seurat::FindClusters(object = dom.object.SCTransform.filter.2, 
                                                               resolution = resolution.range)
#check data
head(dom.object.SCTransform.filter.clustree[[]])
#clustree
clustree(dom.object.SCTransform.filter.clustree, 
         prefix = "SCT_snn_res.")+
  ggtitle('Dom clusters across resolutions')+
  guides(edge_colour = FALSE, 
         edge_alpha = FALSE) +
  theme(legend.position = "bottom")
ggsave('comparison/dom.clustree.png')
#reduced range
dom.object.SCTransform.filter.clustree.reduced <- Seurat::FindClusters(object = dom.object.SCTransform.filter.2, 
                                                               resolution = resolution.range.reduced)
#clustree
clustree(dom.object.SCTransform.filter.clustree.reduced, 
         prefix = "SCT_snn_res.")+
  ggtitle('Dom clusters across resolutions reduced') +
  guides(edge_colour = FALSE, 
         edge_alpha = FALSE) +
  theme(legend.position = "bottom")+
  scale_edge_color_continuous(low = "black", high = "black")
ggsave('comparison/dom.clustree.reduced.png')

##sub
#create dummy dataset
sub.object.SCTransform.filter.2 = sub.object.SCTransform
# cluster across range
sub.object.SCTransform.filter.clustree <- Seurat::FindClusters(object = sub.object.SCTransform.filter.2, 
                                                               resolution = resolution.range)
#check data
head(sub.object.SCTransform.filter.clustree[[]])
#clustree
clustree(sub.object.SCTransform.filter.clustree, 
         prefix = "SCT_snn_res.")+
  ggtitle('sub clusters across resolutions')+
  guides(edge_colour = FALSE, 
         edge_alpha = FALSE) +
  theme(legend.position = "bottom")
ggsave('comparison/sub.clustree.png')
#reduced range
sub.object.SCTransform.filter.clustree.reduced <- Seurat::FindClusters(object = sub.object.SCTransform.filter.2, 
                                                                       resolution = resolution.range.reduced)
#clustree
clustree(sub.object.SCTransform.filter.clustree.reduced, 
         prefix = "SCT_snn_res.")+
  ggtitle('sub clusters across resolutions reduced') +
  guides(edge_colour = FALSE, 
         edge_alpha = FALSE) +
  theme(legend.position = "bottom")+
  scale_edge_color_continuous(low = "black", high = "black")
ggsave('comparison/sub.clustree.reduced.png')

#### integrate dom and sub data ####
### integrate data
## https://satijalab.org/seurat/articles/integration_introduction.html

## create list
burtoni.snseq.list <- list(dom_burtoni_snseq = dom.object.SCTransform, 
                           sub_burtoni_snseq = sub.object.SCTransform)

##select features
burtoni.snseq.features <- SelectIntegrationFeatures(object.list = burtoni.snseq.list,
                                                    nfeatures = 3000)
## identify anchors
burtoni.snseq.list <- PrepSCTIntegration(object.list = burtoni.snseq.list,
                                         anchor.features = burtoni.snseq.features)
## find anchors
burtoni.snseq.anchors <- FindIntegrationAnchors(object.list = burtoni.snseq.list,
                                         normalization.method = "SCT",
                                         anchor.features = burtoni.snseq.features)
## integrate data
burtoni.snseq.combined.sct <- IntegrateData(anchorset = burtoni.snseq.anchors,
                                     normalization.method = "SCT")
## Run PCA
burtoni.snseq.combined.sct <- RunPCA(burtoni.snseq.combined.sct, 
                                     verbose = FALSE)
#check elbowplot
ElbowPlot(burtoni.snseq.combined.sct)
ggsave('comparison/Elbowplot.comparison.png')

## run UMAP
burtoni.snseq.combined.sct <- RunUMAP(burtoni.snseq.combined.sct,
                                      reduction = "pca", 
                                      dims = 1:15)


# find neighbors
burtoni.snseq.combined.sct <- FindNeighbors(burtoni.snseq.combined.sct, reduction = "pca", dims = 1:30)
# # find clusters
burtoni.snseq.combined.sct <- FindClusters(burtoni.snseq.combined.sct, resolution = 0.8)

### graph 
##dom vs sub
DimPlot(burtoni.snseq.combined.sct, 
        reduction = "umap", 
        group.by = "orig.ident")
ggsave('comparison/DomVsSub.dimplot.comparison.png')
## clusters
DimPlot(burtoni.snseq.combined.sct, 
        reduction = "umap", 
        label = TRUE,
        repel = TRUE)
ggsave('comparison/Clusters.dimplot.comparison.png')


# ### oxt and avp
# ##across treatments
# FeaturePlot(burtoni.snseq.combined.sct, 
#             features = c("avp",
#                          "oxt"), 
#             split.by = "orig.ident", 
#             max.cutoff = 3,
#             cols = c("grey", "red"))
# ggsave('comparison/AVP.OXT.features.comparison.png')
# #per cluster
# VlnPlot(burtoni.snseq.combined.sct, 
#         features = c("avp", 
#                      "oxt"), 
#         split.by = "orig.ident",
#         pt.size = 0, 
#         combine = FALSE) %>% 
#   wrap_plots(ncol = 1)
# ggsave('comparison/AVP.OXT.clusters.comparison.png')


### clustree
# cluster across resolutions
burtoni.snseq.combined.sct.clustree <- Seurat::FindClusters(object = burtoni.snseq.combined.sct,                                                  resolution = resolution.range.reduced.2)
#check data
head(burtoni.snseq.combined.sct.clustree[[]])

#clustree
clustree(burtoni.snseq.combined.sct.clustree, 
         prefix = "integrated_snn_res.",
         node_size_range = c(10,20),
         node_colour = 'cluster') +
  scale_edge_color_continuous(low = "black", 
                              high = "black") +
  theme(legend.position = "bottom")
ggsave('comparison/Combined.clustree.png',
       width = 10,
       height = 10)

## check cells per group per cluster
#get table
DomvsSub.cells.per.cluster.table = table(burtoni.snseq.combined.sct@active.ident, 
                                         burtoni.snseq.combined.sct@meta.data$orig.ident) %>% 
  as.data.frame.matrix() %>% 
  rownames_to_column("cluster.id") %>% 
  pivot_longer(cols = c("dom_burtoni_snseq",
                        "sub_burtoni_snseq"),
               names_to = "orig.ident",
               values_to = "count") %>% 
  group_by(orig.ident) %>% 
  mutate(total = sum(count)) %>% 
  ungroup() %>% 
  mutate(percentage = 100*count/total) %>% 
  as.data.frame() %>% 
  mutate(cluster.id = as.numeric(cluster.id),
         count.2 = ifelse(orig.ident == 'sub_burtoni_snseq',
                          -count,
                          count),
         percentage.2 = ifelse(orig.ident == 'sub_burtoni_snseq',
                          -percentage,
                          percentage)) %>% 
  group_by(cluster.id) %>% 
  mutate(count.diff = sum(count.2),
         percentage.diff = sum(percentage.2)) %>% 
  ungroup() %>% 
  as.data.frame() %>% 
  mutate(cluster.id = as.numeric(cluster.id),
         count_bias = ifelse(count.diff > 0,
                             'Dom',
                             'Sub'),
         percentage_bias = ifelse(percentage.diff > 0,
                             'Dom',
                             'Sub')) 

#count
DomvsSub.cells.per.cluster.table %>% 
  ggplot(aes(x = cluster.id,
             y = count,
             group = orig.ident,
             color = orig.ident)) +
  geom_point() +
  theme_bw()
ggsave('comparison/Cells.per.cluster.DomvsSub.png',
       width = 10,
       height = 10)

#percentage
DomvsSub.cells.per.cluster.table  %>% 
  ggplot(aes(x = cluster.id,
             y = percentage,
             group = orig.ident,
             color = orig.ident)) +
  geom_point() +
  theme_bw()
ggsave('comparison/Cells.per.cluster.DomvsSub.percentage.png',
       width = 10,
       height = 10)

#count difference
DomvsSub.cells.per.cluster.table %>% 
  filter(orig.ident == 'dom_burtoni_snseq') %>% 
  ggplot(aes(y = cluster.id,
             x = count.diff,
             label = cluster.id,
             color = count_bias)) +
  geom_segment( aes(x=0, 
                    xend=count.diff, 
                    y=cluster.id, 
                    yend=cluster.id)) +
  geom_label() +
  theme_bw() +
  xlim(-max(abs(DomvsSub.cells.per.cluster.table$count.diff)),
       max(abs(DomvsSub.cells.per.cluster.table$count.diff)))
ggsave('comparison/Cells.per.cluster.DomvsSub.count.diff.png',
       width = 10,
       height = 10)

#percentage difference
DomvsSub.cells.per.cluster.table %>% 
  filter(orig.ident == 'dom_burtoni_snseq') %>% 
  ggplot(aes(y = cluster.id,
             x = percentage.diff,
             label = cluster.id,
             color = percentage_bias)) +
  geom_segment( aes(x=0, 
                    xend=percentage.diff, 
                    y=cluster.id, 
                    yend=cluster.id)) +
  geom_label() +
  theme_bw() +
  xlim(-max(abs(DomvsSub.cells.per.cluster.table$percentage.diff)),
       max(abs(DomvsSub.cells.per.cluster.table$percentage.diff)))
ggsave('comparison/Cells.per.cluster.DomvsSub.percentage.diff.png',
       width = 10,
       height = 10)




#### markers per clusters ####
### top markers per cluster
##find all markers
burtoni.snseq.combined.sct.markers <- FindAllMarkers(burtoni.snseq.combined.sct, 
                                                     only.pos = TRUE, 
                                                     min.pct = 0.25, 
                                                     logfc.threshold = 0.25)
# check top 2 for each cluster
# burtoni.snseq.combined.sct.markers %>%
#   group_by(cluster) %>%
#   slice_max(n = 2, 
#             order_by = avg_log2FC)

##graph
#identify top 10 per cluster
burtoni.snseq.combined.sct.markers.top10 = burtoni.snseq.combined.sct.markers %>% 
  group_by(cluster) %>%
  top_n(n = 10, 
        wt = avg_log2FC) %>% 
  mutate(cluster.position = as.numeric(as.character(cluster))) %>% 
  arrange(cluster.position)
#graph heatmap
DoHeatmap(burtoni.snseq.combined.sct,
          features = burtoni.snseq.combined.sct.markers.top10$gene) +
  NoLegend()
ggsave('comparison/Top 10 marker genes per cluster.png')

#ordered?
DoHeatmap(burtoni.snseq.combined.sct,
          features = burtoni.snseq.combined.sct.markers.top10 %>% 
            pull(gene)) +
  NoLegend()
ggsave('comparison/Top 10 marker genes per cluster order.png')


#### identify clusters ####

###View full gene list
# burtoni.snseq.all.genes %>%
#   View()

### check overlap of cluster markers and marker genes
## merge cell markers with cluster markers
#   burtoni.snseq.combined.sct.markers.celltype = full_join(burtoni.snseq.combined.sct.markers,
#                                                         marker.genes.list %>%
#                                                           rename(gene = gene.id.tilapia))
# 
# ## remove all NA
# burtoni.snseq.combined.sct.markers.celltype.subset = burtoni.snseq.combined.sct.markers.celltype %>%
#   na.omit()
#   
# 
# #check number of cells expressing gene above threshold
# sum(GetAssayData(object = burtoni.snseq.combined.sct.all, 
#                  slot = "data")["avp",]>5)
# 
# ###check top 10 markers per cluster
# #remove ensemble gene names
# burtoni.snseq.combined.sct.markers.top10 %>% 
#   filter(!str_detect(gene,
#                      "ENSONIG")) %>% 
#   View()
# 
# full_join(burtoni.snseq.combined.sct.markers.top10,
#           marker.genes.list %>% 
#             rename(gene = gene.id.tilapia)) %>% 
#   na.omit() %>% View()
# 
# ## heat map of markers per cluster
# #reorder levels
# levels(burtoni.snseq.combined.sct.all) <- factor(0:23)
# #reorder markers
# top20 <- burtoni.snseq.combined.sct.markers.celltype %>% 
#   group_by(cluster) %>% 
#   top_n(20, 
#         avg_logFC)
# 
# # across all clusters
# # DoHeatmap(burtoni.snseq.combined.sct.all, 
# #           features = burtoni.snseq.combined.sct.markers.celltype %>% 
# #             group_by(cluster) %>% 
# #             arrange(desc(avg_log2FC)) %>% 
# #             slice_head(n = 20) %>% 
# #             pull(gene), 
# #           size = 3)
# 
# #cluster specific
# DoHeatmap(burtoni.snseq.combined.sct.all, 
#           features = burtoni.snseq.combined.sct.markers.celltype.subset %>% 
#             pull(gene), 
#           size = 3)
# 
# #all markers
# DoHeatmap(burtoni.snseq.combined.sct.all, 
#           features = marker.genes.list %>% 
#             pull(gene.id.tilapia), 
#           size = 3)



## create marker score for each broad cell type
## Takes a while to fully run
for (x in cell.types) {
  #create marker score 
  burtoni.snseq.combined.sct <- AddModuleScore(object = burtoni.snseq.combined.sct,
                                               assay = 'SCT',
                                                   features = list(c(marker.genes.list %>% 
                                                                       filter(cell.subtype.category == 'broad') %>% 
                                                                       filter(cell.subtype == x) %>% 
                                                                       pull(gene.id.tilapia))), 
                                                   name = paste(x,
                                                                "score",
                                                                sep = '.'))
}



## graph cell type marker score
for (x in cell.types) {
  #graph cell type umap
  FeaturePlot(object = burtoni.snseq.combined.sct,
              features = paste(x,
                               ".score1",
                               sep =''),
              min.cutoff = 0.1,
              pt.size = 0.5,
              order = TRUE,
              cols = c("grey",
                       "red")) +
    ggtitle(paste(x,
                  "marker score",
                  sep = ' '),
            subtitle = paste("genes =",
                             marker.genes.list %>%
                               filter(cell.subtype.category == 'broad') %>%
                               filter(cell.subtype == x) %>%
                               pull(gene.id.tilapia) %>%
                               length))
  ggsave(paste('cell.markers/marker.umap/',
               x,
               ' marker score umap.png',
               sep = ''),
         width = 10,
         height = 10)
  
  #create dotplot for markers of each cell type 
  DotPlot(object = burtoni.snseq.combined.sct,
          assay = 'SCT',
          features = marker.genes.list %>% 
            filter(cell.subtype.category == 'broad') %>% 
            filter(cell.subtype == x) %>% 
            pull(gene.id.tilapia),
          col.min = 1,
          dot.min = 0.1) +
    ylab('Cluster ID') + 
    theme(axis.text.x = element_text(angle = 45,
                                     hjust = 1)) +
    theme_classic()+ 
    ggtitle(paste(x,
                  ' marker genes',
                  sep = ''))
  ggsave(paste('cell.markers/marker.dotplot/cell type marker enrichment ',
               x,
               ' dotplot.png'),
         width = 10,
         height = 10)
}

# ###hoverplot
# FeaturePlot(object = burtoni.snseq.combined.sct.all, 
#             features = "astrocytes.score1",
#             min.cutoff = 1,
#             max.cutoff = 2,
#             pt.size = 0.5,
#             order = TRUE,
#             cols = c("grey",
#                      "red")) %>% 
#   HoverLocator(information = FetchData(object = burtoni.snseq.combined.sct.all, 
#                                        vars = 'ident'))

### dotplot
## graph
DotPlot(object = burtoni.snseq.combined.sct, 
        features = burtoni.snseq.combined.sct@meta.data %>% 
          select(ends_with(".score1")) %>% 
          colnames(),
        col.min = 0,
        dot.min = .1) +
  ylab('Cluster ID') + 
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1))
ggsave('cell.markers/cell type marker enrichment per cluster dotplot.png',
       width = 10,
       height = 10)

##check data
##use marker score
cell.types.marker.per.cluster = DotPlot(object = burtoni.snseq.combined.sct, 
                                        features = burtoni.snseq.combined.sct@meta.data %>% 
                                          select(ends_with(".score1")) %>% 
                                          colnames())

#get counts of cell per cluster
cell.types.marker.per.cluster.count = table(burtoni.snseq.combined.sct@meta.data$seurat_clusters) %>% 
  as.data.frame() %>% 
  rename(id = Var1) %>% 
  rename(number.cells = Freq) %>% 
  full_join(cell.types.marker.per.cluster$data) %>% 
  mutate(number.exp = number.cells*pct.exp/100)

##use marker gene score
cell.types.gene.marker.per.cluster = DotPlot(object = burtoni.snseq.combined.sct,
                                             assay = 'SCT',
                                             features = marker.genes.list %>% 
                                               filter(cell.subtype.category == 'broad') %>% 
                                               pull(gene.id.tilapia))

#get counts of cell per cluster
cell.types.gene.marker.per.cluster.count = table(burtoni.snseq.combined.sct@meta.data$seurat_clusters) %>% 
  as.data.frame() %>% 
  rename(id = Var1) %>% 
  rename(number.cells = Freq) %>% 
  full_join(cell.types.gene.marker.per.cluster$data) %>% 
  mutate(number.exp = number.cells*pct.exp/100) %>% 
  rename(gene.id.tilapia = features.plot) %>% 
  left_join(marker.genes.list)

## combine gene marker and marker score
cell.types.per.cluster.count = full_join(cell.types.gene.marker.per.cluster.count, 
                                         cell.types.marker.per.cluster.count) %>% 
  mutate(feature.id = ifelse(is.na(gene.id.tilapia),
                             features.plot %>% 
                               as.character(),
                             cell.subtype),
         weight.pct.avg = avg.exp.scaled*pct.exp/100)

#graph
cell.types.per.cluster.count %>% 
  filter(weight.pct.avg > 0) %>% 
  filter(id == '5') %>% 
  ggplot(aes(x= avg.exp.scaled,
             y= pct.exp,
             label = feature.id,
             fill = weight.pct.avg)) +
  geom_label() +
  theme_bw()

#avg weight
cell.types.per.cluster.count %>% 
  filter(weight.pct.avg > 0) %>%
  ggplot(aes(x= id,
             y= weight.pct.avg,
             label = features.plot)) +
  geom_text(angle = 90) +
  theme_bw()



# ### remap data to clusters
# #make dummy dataset
# burtoni.snseq.combined.sct.all.label = burtoni.snseq.combined.sct.all
# #create cluster ids
# new.cluster.ids <- c("0", 
#                      '1',
#                      '2',
#                      'neuron.3',
#                      '4',
#                      'astrocyte',
#                      '6',
#                      '7',
#                      'oligodendrocyte',
#                      'neuron.9',
#                      '10',
#                      'inhibitory.11',
#                      '12',
#                      '13',
#                      '14',
#                      'macrophage',
#                      '16',
#                      'endothelial.17',
#                      'endothelial.18',
#                      '19',
#                      'inhibitory.20',
#                      'neuron.21',
#                      '22',
#                      '23')
# #add names as levels
# names(new.cluster.ids) <- levels(burtoni.snseq.combined.sct.all.label)
# #rename clusters
# burtoni.snseq.combined.sct.all.label <- RenameIdents(burtoni.snseq.combined.sct.all.label, 
#                      new.cluster.ids)
# #remap
# DimPlot(burtoni.snseq.combined.sct.all.label,
#         reduction = "umap",
#         label = T,
#         pt.size = 1,
#         label.size = 10,
#         repel = T) +
#   theme(legend.position = 'none')
# ggsave('comparison/Labelled.clusters.dimplot.comparison.all.png',
#        width = 10,
#        height = 10)
# #no label
# DimPlot(burtoni.snseq.combined.sct.all.label, 
#         reduction = "umap",
#         pt.size = 1) +
#   theme(legend.position = 'none') 
# ggsave('comparison/No.labelled.clusters.dimplot.comparison.all.png',
#        width = 10,
#        height = 10)

# ## hover plot
# DimPlot(burtoni.snseq.combined.sct.all,
#         reduction = "umap",
#         label = T,
#         pt.size = 1,
#         label.size = 10,
#         repel = T) %>%
#   HoverLocator(information = FetchData(object = burtoni.snseq.combined.sct.all,
#                                        vars = 'ident'))

#### integrate dom and sub data all genes ####
### integrate data
## https://satijalab.org/seurat/articles/integration_introduction.html

## create list
burtoni.snseq.list.all <- list(dom_burtoni_snseq = dom.object.SCTransform, 
                           sub_burtoni_snseq = sub.object.SCTransform)

##select features
burtoni.snseq.features.all <- SelectIntegrationFeatures(object.list = burtoni.snseq.list.all,
                                                        nfeatures = 3000)
## identify anchors
burtoni.snseq.list.all <- PrepSCTIntegration(object.list = burtoni.snseq.list.all,
                                         anchor.features = burtoni.snseq.features.all)
## find anchors
#20 min
burtoni.snseq.anchors.all <- FindIntegrationAnchors(object.list = burtoni.snseq.list.all,
                                                normalization.method = "SCT",
                                                anchor.features = burtoni.snseq.features.all)
##create all gene list
burtoni.snseq.all.genes <- lapply(burtoni.snseq.list.all, 
                       row.names) %>% 
  Reduce(intersect, .) 
## integrate data
#retain all genes
burtoni.snseq.combined.sct.all <- IntegrateData(anchorset = burtoni.snseq.anchors.all,
                                            normalization.method = "SCT",
                                            features.to.integrate = burtoni.snseq.all.genes)
## Run PCA
burtoni.snseq.combined.sct.all <- RunPCA(burtoni.snseq.combined.sct.all, 
                                     verbose = FALSE)
#check elbowplot
ElbowPlot(burtoni.snseq.combined.sct.all)
ggsave('comparison/Elbowplot.comparison.all.png')


## run find neighbors
burtoni.snseq.combined.sct.all <- FindNeighbors(burtoni.snseq.combined.sct.all)

## run UMAP
burtoni.snseq.combined.sct.all <- RunUMAP(burtoni.snseq.combined.sct.all,
                                      reduction = "pca", 
                                      dims = 1:15)

## 

### graph 
##dom vs sub
DimPlot(burtoni.snseq.combined.sct.all, 
        reduction = "umap", 
        group.by = "orig.ident")
ggsave('comparison/DomVsSub.dimplot.comparison.all.png',
       width = 10,
       height = 10)
## clusters
DimPlot(burtoni.snseq.combined.sct.all, 
        reduction = "umap", 
        label = TRUE,
        repel = TRUE)
# ggsave('comparison/Clusters.dimplot.comparison.all.png',
#        width = 10,
#        height = 10)

#hoverplot
DimPlot(burtoni.snseq.combined.sct.all, 
        reduction = "umap", 
        label = TRUE,
        repel = TRUE) %>% 
  HoverLocator(information = FetchData(object = burtoni.snseq.combined.sct.all, 
                                       vars = 'ident'))


### presentation 
#set presentation colors
presentation.color <- c('#66c2a5',
                        '#fc8d62',
                        '#8da0cb',
                        '#e78ac3',
                        '#a6d854',
                        '#ffd92f',
                        '#e5c494',
                        '#b3b3b3')
#empty cluster
DimPlot(burtoni.snseq.combined.sct.all, 
        reduction = "umap",
        cols = c(rep('grey', 
                     24)),
        pt.size = 1) +
  theme(legend.position = 'none')
ggsave('comparison/Empty.dimplot.comparison.all.png',
       width = 10,
       height = 10)
## presentation simple cluster
DimPlot(burtoni.snseq.combined.sct.all, 
        reduction = "umap",
        label = T,
        pt.size = 1,
        label.size = 10,
        repel = T) +
  theme(legend.position = 'none') 
ggsave('comparison/Clusters.dimplot.comparison.all.png',
       width = 10,
       height = 10)

## presentation across samples cluster
DimPlot(burtoni.snseq.combined.sct.all, 
        reduction = "umap",
        split.by = 'orig.ident',
        pt.size = 1) +
  theme(legend.position = 'none') 
ggsave('comparison/DomvsSub.clusters.dimplot.comparison.all.png',
       width = 20,
       height = 10)

DimPlot(burtoni.snseq.combined.sct.all, 
        split.by = 'orig.ident',
        pt.size = 1) +
  theme(legend.position = 'none') 

### clustree
# cluster across resolutions
burtoni.snseq.combined.sct.all.clustree <- Seurat::FindClusters(object = burtoni.snseq.combined.sct.all,                                                  resolution = resolution.range.reduced.2)
#check data
head(burtoni.snseq.combined.sct.all.clustree[[]])

#clustree
clustree(burtoni.snseq.combined.sct.all.clustree, 
         prefix = "integrated_snn_res.",
         node_size_range = c(10,20),
         node_colour = 'cluster') +
  scale_edge_color_continuous(low = "black", 
                              high = "black") +
  theme(legend.position = "bottom")
ggsave('comparison/Combined.clustree.png',
       width = 10,
       height = 10)

## check cells per group per cluster
#get table
DomvsSub.cells.per.cluster.table = table(burtoni.snseq.combined.sct.all@active.ident, 
                                    burtoni.snseq.combined.sct.all@meta.data$orig.ident) %>% 
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
DomvsSub.cells.per.cluster.table %>% 
  as.data.frame() %>% 
  mutate(cluster.id = as.numeric(cluster.id)) %>% 
  ggplot(aes(x = cluster.id,
             y = percentage,
             group = orig.ident,
             color = orig.ident)) +
  geom_point() +
  theme_bw()
ggsave('comparison/Cells.per.cluster.DomvsSub.png',
       width = 10,
       height = 10)

#### normalize SCT across samples ####
## rescale SCT assay for both samples
## subset genes
# use minimum count 1772
burtoni.snseq.combined.sct = PrepSCTFindMarkers(burtoni.snseq.combined.sct,
                                                    assay = 'SCT',
                                                    verbose = T)

## all genes
burtoni.snseq.combined.sct.all = PrepSCTFindMarkers(burtoni.snseq.combined.sct.all,
                                                    assay = 'SCT',
                                                    verbose = T)
#### identify clusters all genes ####

###View full gene list
# burtoni.snseq.all.genes %>%
#   View()

### check overlap of cluster markers and marker genes
## merge cell markers with cluster markers
#   burtoni.snseq.combined.sct.markers.celltype = full_join(burtoni.snseq.combined.sct.markers,
#                                                         marker.genes.list %>%
#                                                           rename(gene = gene.id.tilapia))
# 
# ## remove all NA
# burtoni.snseq.combined.sct.markers.celltype.subset = burtoni.snseq.combined.sct.markers.celltype %>%
#   na.omit()
#   
# 
# #check number of cells expressing gene above threshold
# sum(GetAssayData(object = burtoni.snseq.combined.sct.all, 
#                  slot = "data")["avp",]>5)
# 
# ###check top 10 markers per cluster
# #remove ensemble gene names
# burtoni.snseq.combined.sct.markers.top10 %>% 
#   filter(!str_detect(gene,
#                      "ENSONIG")) %>% 
#   View()
# 
# full_join(burtoni.snseq.combined.sct.markers.top10,
#           marker.genes.list %>% 
#             rename(gene = gene.id.tilapia)) %>% 
#   na.omit() %>% View()
# 
# ## heat map of markers per cluster
# #reorder levels
# levels(burtoni.snseq.combined.sct.all) <- factor(0:23)
# #reorder markers
# top20 <- burtoni.snseq.combined.sct.markers.celltype %>% 
#   group_by(cluster) %>% 
#   top_n(20, 
#         avg_logFC)
# 
# # across all clusters
# # DoHeatmap(burtoni.snseq.combined.sct.all, 
# #           features = burtoni.snseq.combined.sct.markers.celltype %>% 
# #             group_by(cluster) %>% 
# #             arrange(desc(avg_log2FC)) %>% 
# #             slice_head(n = 20) %>% 
# #             pull(gene), 
# #           size = 3)
# 
# #cluster specific
# DoHeatmap(burtoni.snseq.combined.sct.all, 
#           features = burtoni.snseq.combined.sct.markers.celltype.subset %>% 
#             pull(gene), 
#           size = 3)
# 
# #all markers
# DoHeatmap(burtoni.snseq.combined.sct.all, 
#           features = marker.genes.list %>% 
#             pull(gene.id.tilapia), 
#           size = 3)


## make list of broad cell types
cell.types = marker.genes.list %>% 
  filter(cell.subtype.category == 'broad') %>% 
  pull(cell.subtype) %>% 
  unique()

## create marker score for each broad cell type
## Takes a while to fully run
for (x in cell.types) {
  #create marker score 
  burtoni.snseq.combined.sct.all <- AddModuleScore(object = burtoni.snseq.combined.sct.all, 
                                                   features = list(c(marker.genes.list %>% 
                                                                       filter(cell.subtype.category == 'broad') %>% 
                                                                       filter(cell.subtype == x) %>% 
                                                                       pull(gene.id.tilapia))), 
                                                   name = paste(x,
                                                                "score",
                                                                sep = '.'))
}



## graph cell type marker score
for (x in cell.types) {
  #graph cell type umap
  FeaturePlot(object = burtoni.snseq.combined.sct.all, 
              features = paste(x,
                               ".score1",
                               sep =''),
              min.cutoff = 1,
              max.cutoff = 2,
              pt.size = 0.5,
              order = TRUE,
              cols = c("grey",
                       "red")) + 
    ggtitle(paste(x,
                  "marker score",
                  sep = ' '), 
            subtitle = paste("genes =", 
                             marker.genes.list %>% 
                                      filter(cell.subtype.category == 'broad') %>% 
                                      filter(cell.subtype == x) %>% 
                                      pull(gene.id.tilapia) %>% 
                               length))
  ggsave(paste('cell.markers/marker.umap/',
               x,
  ' marker score umap.png',
  sep = ''),
         width = 10,
         height = 10)
  
  #create dotplot for markers of each cell type 
  DotPlot(object = burtoni.snseq.combined.sct.all, 
          features = marker.genes.list %>% 
            filter(cell.subtype.category == 'broad') %>% 
            filter(cell.subtype == x) %>% 
            pull(gene.id.tilapia),
          col.min = 1,
          dot.min = 0.1) +
    ylab('Cluster ID') + 
    theme(axis.text.x = element_text(angle = 45,
                                     hjust = 1)) +
    ggtitle(paste(x,
                  ' marker genes',
                  sep = ''))
  ggsave(paste('cell.markers/marker.dotplot/cell type marker enrichment ',
               x,
               ' dotplot.png'),
         width = 10,
         height = 10)
}

# ###hoverplot
# FeaturePlot(object = burtoni.snseq.combined.sct.all, 
#             features = "astrocytes.score1",
#             min.cutoff = 1,
#             max.cutoff = 2,
#             pt.size = 0.5,
#             order = TRUE,
#             cols = c("grey",
#                      "red")) %>% 
#   HoverLocator(information = FetchData(object = burtoni.snseq.combined.sct.all, 
#                                        vars = 'ident'))

### dotplot
## graph
DotPlot(object = burtoni.snseq.combined.sct.all, 
        features = burtoni.snseq.combined.sct.all@meta.data %>% 
  select(ends_with(".score1")) %>% 
  colnames(),
  col.min = 0,
  dot.min = .1) +
  ylab('Cluster ID') + 
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1))
ggsave('cell.markers/cell type marker enrichment per cluster dotplot.png',
       width = 10,
       height = 10)

##check data
##use marker score
cell.types.marker.per.cluster = DotPlot(object = burtoni.snseq.combined.sct.all, 
        features = burtoni.snseq.combined.sct.all@meta.data %>% 
          select(ends_with(".score1")) %>% 
          colnames())

#get counts of cell per cluster
cell.types.marker.per.cluster.count = table(burtoni.snseq.combined.sct.all@meta.data$seurat_clusters) %>% 
  as.data.frame() %>% 
  rename(id = Var1) %>% 
  rename(number.cells = Freq) %>% 
  full_join(cell.types.marker.per.cluster$data) %>% 
  mutate(number.exp = number.cells*pct.exp/100)

##use marker gene score
cell.types.gene.marker.per.cluster = DotPlot(object = burtoni.snseq.combined.sct.all, 
                                        features = marker.genes.list %>% 
                                          filter(cell.subtype.category == 'broad') %>% 
                                          pull(gene.id.tilapia))

#get counts of cell per cluster
cell.types.gene.marker.per.cluster.count = table(burtoni.snseq.combined.sct.all@meta.data$seurat_clusters) %>% 
  as.data.frame() %>% 
  rename(id = Var1) %>% 
  rename(number.cells = Freq) %>% 
  full_join(cell.types.gene.marker.per.cluster$data) %>% 
  mutate(number.exp = number.cells*pct.exp/100) %>% 
  rename(gene.id.tilapia = features.plot) %>% 
  left_join(marker.genes.list)

## combine gene marker and marker score
cell.types.per.cluster.count = full_join(cell.types.gene.marker.per.cluster.count, 
                                         cell.types.marker.per.cluster.count) %>% 
  mutate(feature.id = ifelse(is.na(gene.id.tilapia),
                             features.plot %>% 
                               as.character(),
                             cell.subtype),
         weight.pct.avg = avg.exp.scaled*pct.exp/100)

#graph
cell.types.per.cluster.count %>% 
  filter(weight.pct.avg > 0) %>% 
  filter(id == '5') %>% 
ggplot(aes(x= avg.exp.scaled,
           y= pct.exp,
           label = feature.id,
           fill = weight.pct.avg)) +
  geom_label() +
  theme_bw()

#avg weight
cell.types.per.cluster.count %>% 
  filter(weight.pct.avg > 0) %>%
  ggplot(aes(x= id,
             y= weight.pct.avg,
             label = features.plot)) +
  geom_text(angle = 90) +
  theme_bw()



# ### remap data to clusters
# #make dummy dataset
# burtoni.snseq.combined.sct.all.label = burtoni.snseq.combined.sct.all
# #create cluster ids
# new.cluster.ids <- c("0", 
#                      '1',
#                      '2',
#                      'neuron.3',
#                      '4',
#                      'astrocyte',
#                      '6',
#                      '7',
#                      'oligodendrocyte',
#                      'neuron.9',
#                      '10',
#                      'inhibitory.11',
#                      '12',
#                      '13',
#                      '14',
#                      'macrophage',
#                      '16',
#                      'endothelial.17',
#                      'endothelial.18',
#                      '19',
#                      'inhibitory.20',
#                      'neuron.21',
#                      '22',
#                      '23')
# #add names as levels
# names(new.cluster.ids) <- levels(burtoni.snseq.combined.sct.all.label)
# #rename clusters
# burtoni.snseq.combined.sct.all.label <- RenameIdents(burtoni.snseq.combined.sct.all.label, 
#                      new.cluster.ids)
# #remap
# DimPlot(burtoni.snseq.combined.sct.all.label,
#         reduction = "umap",
#         label = T,
#         pt.size = 1,
#         label.size = 10,
#         repel = T) +
#   theme(legend.position = 'none')
# ggsave('comparison/Labelled.clusters.dimplot.comparison.all.png',
#        width = 10,
#        height = 10)
# #no label
# DimPlot(burtoni.snseq.combined.sct.all.label, 
#         reduction = "umap",
#         pt.size = 1) +
#   theme(legend.position = 'none') 
# ggsave('comparison/No.labelled.clusters.dimplot.comparison.all.png',
#        width = 10,
#        height = 10)

# ## hover plot
# DimPlot(burtoni.snseq.combined.sct.all,
#         reduction = "umap",
#         label = T,
#         pt.size = 1,
#         label.size = 10,
#         repel = T) %>%
#   HoverLocator(information = FetchData(object = burtoni.snseq.combined.sct.all,
#                                        vars = 'ident'))

#### save data point ####
##just single cell object
# variable genes for integration
save(burtoni.snseq.combined.sct,
     file = "burtoni.snseq.combined.sct.RData")
#all
save(burtoni.snseq.combined.sct.all,
file = "burtoni.snseq.combined.sct.all.RData")
# load('burtoni.snseq.combined.sct.all.RData')
# save.image('burtoni.snseq.SCTransform.RData')
# load('burtoni.snseq.SCTransform.RData')

## save out barcodes
burtoni.snseq.selected.barcodes = burtoni.snseq.combined.sct.all@assays$RNA@data@Dimnames[[2]] %>% 
  as.data.frame()
#rename column
colnames(burtoni.snseq.selected.barcodes) = 'Barcode.id'
#seperate barcode
burtoni.snseq.selected.barcodes = burtoni.snseq.selected.barcodes %>% 
  separate(Barcode.id, 
           into = c("Barcode","orig.ident"),
           sep = '_') %>% 
  mutate(orig.ident = ifelse(orig.ident == 1,
                             'dom',
                             'sub'))

#dom
dom.burtoni.snseq.selected.barcodes = burtoni.snseq.selected.barcodes %>% 
  filter(orig.ident == 'dom') %>% 
  select(Barcode)
#sub
sub.burtoni.snseq.selected.barcodes = burtoni.snseq.selected.barcodes %>% 
  filter(orig.ident == 'sub') %>% 
  select(Barcode)

##save tsv
#dom
write_tsv(dom.burtoni.snseq.selected.barcodes, 
          'dom.burtoni.snseq.selected.barcodes.tsv',
          col_names = FALSE)

#sub
write_tsv(dom.burtoni.snseq.selected.barcodes, 
          'sub.burtoni.snseq.selected.barcodes.tsv',
          col_names = FALSE)
