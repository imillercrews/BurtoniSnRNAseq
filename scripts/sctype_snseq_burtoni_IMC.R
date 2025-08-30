#### Burtoni snseq seurat analysis
### cell type assignment
###Note: Seurat requires R version > 4
## use lambcomp1 to run R with command 
# > R-4.0.3
###SCtype
##https://github.com/IanevskiAleksandr/sc-type/


### set working directory
setwd("/stor/work/Hofmann/All_projects/A_burtoni_snseq/seurat/")

#### load libraries ####
#load libraries
library(cowplot)
library(Seurat)
library(HGNChelper)
library(tidyverse)
library(ggalluvial)

# load gene set preparation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

#### load data ####
### single cell data
### single cell data
#load('burtoni.snseq.combined.sct.all.RData')
load('burtoni.snseq.combined.sct.RData')

### scsorter data
# load("burtoni.scsorter.data.scsort.output.RData")
load("burtoni.scsorter.data.scsort.output.new.RData")

### load souporcell data
burtoni.souporcell.filtered = read.csv('../souporcell/burtoni.souporcell.filtered.csv')


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

#### annotation data ####
## subset data from 'Cell.type.gene.markers.csv'
## load annotation data
# annotation.data = read.csv("../Gene.lists/Cell.type.annotation.genes.csv")
# ## remove NA
# annotation.data = annotation.data %>%
#   na.omit() %>%
#   distinct()
# 
# ##rename duplicate rownames
# ## can't handle upper/lowercase of same genes
# ## ex. ‘SLC32A1’, ‘SLCO1C1’
# annotation.data = annotation.data %>%
#   filter(Marker != "SLCO1C1") %>%
#   filter(Marker != "SLC32A1")
# 
# ## reduce to genes in data
# annotation.data.reduce= annotation.data %>% 
#   filter(Marker %in% rownames(burtoni.snseq.combined.sct@assays$SCT@scale.data))
# 
# # check reduced gene list
# annotation.data %>% 
#   dplyr::select(Type) %>% 
#   table() %>% 
#   as.data.frame() %>% 
#   dplyr::rename("Cell.type" = ".") %>% 
#   ggplot() +
#   geom_col(aes(x = reorder(Cell.type, -Freq),
#                y = Freq),
#            fill = 'black') + 
#   geom_col(data = annotation.data.reduce %>% 
#              dplyr::select(Type) %>% 
#              table() %>% 
#              as.data.frame() %>% 
#              dplyr::rename("Cell.type" = "."),
#            aes(x = reorder(Cell.type, -Freq),
#                y = Freq),
#            fill = 'red') +
#   theme_classic() +
#   ggtitle("Genes per cell type")
# ggsave('./sctype/Genes per cell type.png',
#        height = 10,
#        width = 10)

# ## collapse dataframe to match format needed
# annotation.data = annotation.data %>% 
#   group_by(Type) %>% 
#   summarise(geneSymbolmore1 = paste0(Marker, 
#                                      collapse = ",")) %>% 
#   rename(cellName = Type) %>% 
#   mutate(tissueType = 'Brain',
#          geneSymbolmore2 = 0) %>% 
#   relocate(cellName, 
#            .after = last_col()) %>% 
#   relocate(geneSymbolmore1, 
#            .after = last_col()) %>% 
#   relocate(geneSymbolmore2, 
#            .after = last_col()) %>% 
#   as.data.frame()

# annotation.data = annotation.data %>% 
#   filter(cellName != 'inhibitory') %>% 
#   filter(cellName != 'excitatory')

# # save as xlsx for sctype
# library(openxlsx)
# write.xlsx(annotation.data, 
#            "../Gene.lists/Cell.type.annotation.genes.sctype.xlsx",
#            sheetName = "Sheet1", 
#            colNames = TRUE, 
#            rowNames = FALSE, 
#            append = FALSE)


### load HypoMap annotation
annotation.data.hypomap = read.csv("../Gene.lists/Cell.type.gene.markers.HypoMap.annotation.ensemble.csv")

## remove NA
#get distinct 
annotation.data.hypomap = annotation.data.hypomap %>% 
  na.omit() %>% 
  dplyr::select(c(tilapia_gene_name,
                  cluster_name)) %>% 
  dplyr::rename('Type' = 'cluster_name') %>% 
  dplyr::rename('Marker' = 'tilapia_gene_name') %>% 
  distinct()

## reduce to genes in data
annotation.data.hypomap.reduce= annotation.data.hypomap %>% 
  filter(Marker %in% rownames(burtoni.snseq.combined.sct@assays$SCT@scale.data))

# # check reduced gene list
# annotation.data.hypomap %>% 
#   dplyr::select(Type) %>% 
#   table() %>% 
#   as.data.frame() %>% 
#   dplyr::rename("Cell.type" = ".") %>% 
#   ggplot() +
#   geom_col(aes(x = reorder(Cell.type, -Freq),
#                y = Freq),
#            fill = 'black') + 
#   geom_col(data = annotation.data.hypomap.reduce %>% 
#              dplyr::select(Type) %>% 
#              table() %>% 
#              as.data.frame() %>% 
#              dplyr::rename("Cell.type" = "."),
#            aes(x = reorder(Cell.type, -Freq),
#                y = Freq),
#            fill = 'red') +
#   theme_classic() +
#   ggtitle("HypoMap genes per cell type")
# ggsave('./sctype.hypo/HypoMap genes per cell type.png',
#        height = 10,
#        width = 10)

## collapse dataframe to match format needed
annotation.data.hypomap = annotation.data.hypomap %>% 
  group_by(Type) %>% 
  summarise(geneSymbolmore1 = paste0(Marker, 
                                     collapse = ",")) %>% 
  rename(cellName = Type) %>% 
  mutate(tissueType = 'Brain',
         geneSymbolmore2 = 0) %>% 
  relocate(cellName, 
           .after = last_col()) %>% 
  relocate(geneSymbolmore1, 
           .after = last_col()) %>% 
  relocate(geneSymbolmore2, 
           .after = last_col()) %>% 
  as.data.frame()



# # save as xlsx for sctype
# library(openxlsx)
# write.xlsx(annotation.data.hypomap,
#            "../Gene.lists/Cell.type.annotation.genes.sctype.hypomap.xlsx",
#            sheetName = "Sheet1",
#            colNames = TRUE,
#            rowNames = FALSE,
#            append = FALSE)

#### functions ####
### modify gene_sets_prepare function to integrate with user added data
## do not require xlsx for gene_sets_prepare
## remove checkGeneSymbols function 
# remove forcing genes to uppercase
gene_sets_prepare.df = function(db_file, cell_type){
  
  cell_markers = db_file #change to dataframe
  cell_markers = cell_markers[cell_markers$tissueType == cell_type,] 
  cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1); cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
  
  # correct gene symbols from the given DB (up-genes)
  cell_markers$geneSymbolmore1 = sapply(1:nrow(cell_markers), function(i){
    
    markers_all = gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore1[i],",")))
    # markers_all = toupper(markers_all[markers_all != "NA" & markers_all != ""])
    markers_all = sort(markers_all)
    
    if(length(markers_all) > 0){
      # suppressMessages({markers_all = unique(na.omit(checkGeneSymbols(markers_all)$Suggested.Symbol))})
      paste0(markers_all, collapse=",")
    } else {
      ""
    }
  })
  
  # correct gene symbols from the given DB (down-genes)
  cell_markers$geneSymbolmore2 = sapply(1:nrow(cell_markers), function(i){
    
    markers_all = gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore2[i],",")))
    # markers_all = toupper(markers_all[markers_all != "NA" & markers_all != ""])
    markers_all = sort(markers_all)
    
    if(length(markers_all) > 0){
      # suppressMessages({markers_all = unique(na.omit(checkGeneSymbols(markers_all)$Suggested.Symbol))})
      paste0(markers_all, collapse=",")
    } else {
      ""
    }
  })
  
  cell_markers$geneSymbolmore1 = gsub("///",",",cell_markers$geneSymbolmore1);cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1)
  cell_markers$geneSymbolmore2 = gsub("///",",",cell_markers$geneSymbolmore2);cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)
  
  gs = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore1[j]),",")))); names(gs) = cell_markers$cellName
  gs2 = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore2[j]),",")))); names(gs2) = cell_markers$cellName
  
  list(gs_positive = gs, gs_negative = gs2)
}

### marker specificity score
# marker_sensitivity_score = function(scRNAseqData, scaled = !0, gs, gs2 = NULL, gene_names_to_uppercase = !0, ...){
#   
#   # check input matrix
#   if(!is.matrix(scRNAseqData)){
#     warning("scRNAseqData doesn't seem to be a matrix")
#   } else {
#     if(sum(dim(scRNAseqData))==0){
#       warning("The dimension of input scRNAseqData matrix equals to 0, is it an empty matrix?")
#     }
#   }
#   
#   # marker sensitivity
#   marker_stat = sort(table(unlist(gs)), decreasing = T); 
#   marker_sensitivity = data.frame(score_marker_sensitivity = scales::rescale(as.numeric(marker_stat), to = c(0,1), from = c(length(gs),1)),
#                                   gene_ = names(marker_stat), stringsAsFactors = !1)
#   
#   # convert gene names to Uppercase
#   if(gene_names_to_uppercase){
#     rownames(scRNAseqData) = toupper(rownames(scRNAseqData));
#   }
#   
#   # subselect genes only found in data
#   names_gs_cp = names(gs); names_gs_2_cp = names(gs2);
#   gs = lapply(1:length(gs), function(d_){ 
#     GeneIndToKeep = rownames(scRNAseqData) %in% as.character(gs[[d_]]); rownames(scRNAseqData)[GeneIndToKeep]})
#   gs2 = lapply(1:length(gs2), function(d_){ 
#     GeneIndToKeep = rownames(scRNAseqData) %in% as.character(gs2[[d_]]); rownames(scRNAseqData)[GeneIndToKeep]})
#   names(gs) = names_gs_cp; names(gs2) = names_gs_2_cp;
#   cell_markers_genes_score = marker_sensitivity[marker_sensitivity$gene_ %in% unique(unlist(gs)),]
#   
#   gs
#   # marker_sensitivity
#   # cell_markers_genes_score
#   # GeneIndToKeep
#   # unique(unlist(gs))
#   # gs
#  }

### try with all marker genes
# sctype_score.all = function(scRNAseqData, scaled = !0, gs, gs2 = NULL, gene_names_to_uppercase = !0, ...){
#   
#   # check input matrix
#   if(!is.matrix(scRNAseqData)){
#     warning("scRNAseqData doesn't seem to be a matrix")
#   } else {
#     if(sum(dim(scRNAseqData))==0){
#       warning("The dimension of input scRNAseqData matrix equals to 0, is it an empty matrix?")
#     }
#   }
#   
#   # marker sensitivity
#   marker_stat = sort(table(unlist(gs)), decreasing = T); 
#   marker_sensitivity = data.frame(score_marker_sensitivity = scales::rescale(as.numeric(marker_stat), to = c(0,1), from = c(length(gs),1)),
#                                   gene_ = names(marker_stat), stringsAsFactors = !1)
#   
#   # convert gene names to Uppercase
#   if(gene_names_to_uppercase){
#     rownames(scRNAseqData) = toupper(rownames(scRNAseqData));
#   }
#   
#   # subselect genes only found in data
#   names_gs_cp = names(gs); names_gs_2_cp = names(gs2);
#   gs = lapply(1:length(gs), function(d_){ 
#     GeneIndToKeep = rownames(scRNAseqData) %in% as.character(gs[[d_]]); rownames(scRNAseqData)[GeneIndToKeep]})
#   gs2 = lapply(1:length(gs2), function(d_){ 
#     GeneIndToKeep = rownames(scRNAseqData) %in% as.character(gs2[[d_]]); rownames(scRNAseqData)[GeneIndToKeep]})
#   names(gs) = names_gs_cp; names(gs2) = names_gs_2_cp;
#   # change to all genes
#   # cell_markers_genes_score = marker_sensitivity[marker_sensitivity$gene_ %in% unique(unlist(gs)),]
#   cell_markers_genes_score = marker_sensitivity
#   
#   # z-scale if not
#   if(!scaled) Z <- t(scale(t(scRNAseqData))) else Z <- scRNAseqData
#   
#   # multiple by marker sensitivity
#   for(jj in 1:nrow(cell_markers_genes_score)){
#     Z[cell_markers_genes_score[jj,"gene_"], ] = Z[cell_markers_genes_score[jj,"gene_"], ] * cell_markers_genes_score[jj, "score_marker_sensitivity"]
#   }
#   
#   # subselect only with marker genes
#   Z = Z[unique(c(unlist(gs),unlist(gs2))), ]
#   
#   # combine scores
#   es = do.call("rbind", lapply(names(gs), function(gss_){ 
#     sapply(1:ncol(Z), function(j) {
#       gs_z = Z[gs[[gss_]], j]; gz_2 = Z[gs2[[gss_]], j] * -1
#       sum_t1 = (sum(gs_z) / sqrt(length(gs_z))); sum_t2 = sum(gz_2) / sqrt(length(gz_2));
#       if(is.na(sum_t2)){
#         sum_t2 = 0;
#       }
#       sum_t1 + sum_t2
#     })
#   })) 
#   
#   dimnames(es) = list(names(gs), colnames(Z))
#   es.max <- es[!apply(is.na(es) | es == "", 1, all),] # remove na rows
#   
#   es.max
# }

# #### run sctype database ####
# #### cell type assignment
# ### load database
# ### prepare gene set list
# ## download gene lists
# db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
# # set tissue type
# tissue = "Brain"
# # prepare gene sets
# gs_list = gene_sets_prepare(db_,
#                             tissue)
# 
# # ## test marker specificity
# # marker.specificity = marker_sensitivity_score(scRNAseqData = burtoni.snseq.combined.sct[["SCT"]]@scale.data, 
# #                                               scaled = TRUE, 
# #                                               gs = gs_list$gs_positive, 
# #                                               gs2 = gs_list$gs_negative,
# #                                               gene_names_to_uppercase = FALSE)
# 
# ### assign cell types to clusters
# # get cell-type by cell matrix
# es.max = sctype_score(scRNAseqData = burtoni.snseq.combined.sct.all[["integrated"]]@scale.data, 
#                       scaled = TRUE, 
#                       gs = gs_list$gs_positive, 
#                       gs2 = gs_list$gs_negative) 
# 
# 
# # NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
# # In case Seurat is used, it is either pbmc[["RNA"]]@scale.data (default), pbmc[["SCT"]]@scale.data, in case sctransform is used for normalization,
# # or pbmc[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.
# 
# # merge by cluster
# cL_results = do.call("rbind", lapply(unique(burtoni.snseq.combined.sct.all@meta.data$seurat_clusters), function(cl){
#   es.max.cl = sort(rowSums(es.max[ ,rownames(burtoni.snseq.combined.sct.all@meta.data[burtoni.snseq.combined.sct.all@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
#   head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(burtoni.snseq.combined.sct.all@meta.data$seurat_clusters==cl)), 10)
# }))
# sctype_scores = cL_results %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
# 
# # set low-confident (low ScType score) clusters to "unknown"
# sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
# print(sctype_scores[,1:3])
# 
# #graph umap
# burtoni.snseq.combined.sct.all@meta.data$customclassif = ""
# for(j in unique(sctype_scores$cluster)){
#   cl_type = sctype_scores[sctype_scores$cluster==j,]; 
#   burtoni.snseq.combined.sct.all@meta.data$customclassif[burtoni.snseq.combined.sct.all@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
# }
# 
# #graph umap
# DimPlot(burtoni.snseq.combined.sct.all, 
#         reduction = "umap", 
#         label = TRUE, 
#         repel = TRUE, 
#         group.by = 'customclassif') +
#   ggtitle('sctype')
# ggsave('./sctype/UMAP cell types sctype.png',
#        height = 5.5,
#        width = 8.5)
# 
# DimPlot(burtoni.snseq.combined.sct.all, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'Cell.type')   

# #### graph alluvial plot database ####
# ### compare sctype to scsorter
# burtoni.snseq.combined.sct.all@meta.data %>% 
#   select(c(customclassif,
#            Cell.type)) %>% 
#   table() %>% 
#   as.data.frame() %>% 
#   mutate(Count = sum(Freq)) %>% 
#   ungroup() %>% 
#   mutate(Freq.scale = 100*Freq/Count) %>% 
#   ggplot(aes(axis1 = reorder(Cell.type, -Freq.scale),
#              axis2 = reorder(customclassif,-Freq.scale),
#              y = Freq.scale)) +
#   geom_alluvium(aes(fill = Cell.type)) +
#   geom_stratum() +
#   geom_text(stat = "stratum",
#             aes(label = after_stat(stratum))) +
#   scale_x_discrete(limits = c("Cell.type", 
#                               "customclassif"),
#                    expand = c(0.15, 0.05)) +
#   scale_fill_viridis_d() +
#   theme_classic() +
#   theme(axis.text.x=element_blank())
# ggsave('./sctype/Alluvial cells by celltype.png',
#        height = 10,
#        width = 10)
# 
# ### compare sctype to scsorter
# burtoni.snseq.combined.sct.all@meta.data %>% 
#   select(c(customclassif,
#            seurat_clusters)) %>% 
#   table() %>% 
#   as.data.frame() %>% 
#   mutate(Count = sum(Freq)) %>% 
#   ungroup() %>% 
#   mutate(Freq.scale = 100*Freq/Count) %>% 
#   ggplot(aes(axis1 = reorder(seurat_clusters,-Freq.scale),
#              axis2 = reorder(customclassif,-Freq.scale),
#              y = Freq.scale)) +
#   geom_alluvium(aes(fill = customclassif)) +
#   geom_stratum() +
#   geom_text(size = 5,
#             stat = "stratum",
#             aes(label = after_stat(stratum))) +
#   scale_x_discrete(limits = c("customclassif", 
#                               "seurat_clusters"),
#                    expand = c(0.15, 0.05)) +
#   scale_fill_viridis_d() +
#   theme_classic() +
#   theme(axis.text.x=element_blank())
# ggsave('./sctype/Alluvial cells by celltype and cluster.png',
#        height = 5.5,
#        width = 6.2)
# 
# 
# 
# 

# #### bubble network database ####
#' # load libraries
#' lapply(c("ggraph","igraph","tidyverse", "data.tree"), library, character.only = T)
#' 
#' # prepare edges
#' cL_results=cL_results[order(cL_results$cluster),]; edges = cL_results; edges$type = paste0(edges$type,"_",edges$cluster); edges$cluster = paste0("cluster ", edges$cluster); edges = edges[,c("cluster", "type")]; colnames(edges) = c("from", "to"); rownames(edges) <- NULL
#' 
#' # prepare nodes
#' nodes_lvl1 = sctype_scores[,c("cluster", "ncells")]; nodes_lvl1$cluster = paste0("cluster ", nodes_lvl1$cluster); nodes_lvl1$Colour = "#f1f1ef"; nodes_lvl1$ord = 1; nodes_lvl1$realname = nodes_lvl1$cluster; nodes_lvl1 = as.data.frame(nodes_lvl1); nodes_lvl2 = c(); 
#' 
#' ccolss=rainbow(24)
#' # ccolss= c("#5f75ae","#92bbb8","#64a841","#e5486e","#de8e06","#eccf5a","#b5aa0f","#e4b680","#7ba39d","#b15928","#ffff99", "#6a3d9a","#cab2d6","#ff7f00","#fdbf6f","#e31a1c","#fb9a99","#33a02c","#b2df8a","#1f78b4","#a6cee3")
#' for (i in 1:length(unique(cL_results$cluster))){
#'   dt_tmp = cL_results[cL_results$cluster == unique(cL_results$cluster)[i], ]; nodes_lvl2 = rbind(nodes_lvl2, data.frame(cluster = paste0(dt_tmp$type,"_",dt_tmp$cluster), ncells = dt_tmp$scores, Colour = ccolss[i], ord = 2, realname = dt_tmp$type))
#' }
#' nodes = rbind(nodes_lvl1, nodes_lvl2); nodes$ncells[nodes$ncells<1] = 1;
#' files_db = openxlsx::read.xlsx(db_)[,c("cellName","shortName")]; files_db = unique(files_db); nodes = merge(nodes, files_db, all.x = T, all.y = F, by.x = "realname", by.y = "cellName", sort = F)
#' nodes$shortName[is.na(nodes$shortName)] = nodes$realname[is.na(nodes$shortName)]; nodes = nodes[,c("cluster", "ncells", "Colour", "ord", "shortName", "realname")]
#' 
#' # need to remove some rows with shortName:
#' #'Immune system' and 'Endothelial'
#' mygraph <- graph_from_data_frame(edges, 
#'                                  vertices=nodes %>% 
#'                                    filter(shortName != 'Immune system',
#'                                           shortName != 'Endothelial'))
#' 
#' # Make the graph
#' gggr<- ggraph(mygraph, layout = 'circlepack', weight=I(ncells)) + 
#'   geom_node_circle(aes(filter=ord==1,fill=I("#F5F5F5"), colour=I("#D3D3D3")), alpha=0.9) + geom_node_circle(aes(filter=ord==2,fill=I(Colour), colour=I("#D3D3D3")), alpha=0.9) +
#'   theme_void() + geom_node_text(aes(filter=ord==2, label=shortName, colour=I("#ffffff"), fill="white", repel = !1, parse = T, size = I(log(ncells,25)*1.5)))+ geom_node_label(aes(filter=ord==1,  label=shortName, colour=I("#000000"), size = I(3), fill="white", parse = T), repel = !0, segment.linetype="dotted")
#' 
#' labelledUMAP = DimPlot(burtoni.snseq.combined.sct.all, reduction = "umap", label = TRUE, repel = TRUE, cols = ccolss)
#' 
#' cowplot::plot_grid(labelledUMAP,
#'                    gggr)
#' 
#' cowplot::plot_grid(gggr)
#' ggsave('sctype/UMAP cell types with network.png',
#'        width = 10,
#'        height = 10)
#' 
#' # No function scater but could just use cowplot
#' # scater::multiplot(DimPlot(burtoni.snseq.combined.sct.all, reduction = "umap", label = TRUE, repel = TRUE, cols = ccolss), gggr, cols = 2)
#' 


# #### run sctype ####
# #### cell type assignment
# ### prepare gene set list
# # set tissue type
# tissue = "Brain"
# # prepare gene sets
# gs_list = gene_sets_prepare.df(annotation.data,
#                             tissue)
# 
# 
# ## test marker specificity
# marker.specificity = marker_sensitivity_score(scRNAseqData = burtoni.snseq.combined.sct[["SCT"]]@scale.data,
#                                               scaled = TRUE,
#                                               gs = gs_list$gs_positive,
#                                               gs2 = gs_list$gs_negative,
#                                               gene_names_to_uppercase = FALSE)
# 
# ### assign cell types to clusters
# # get cell-type by cell matrix
# es.max = sctype_score(scRNAseqData = burtoni.snseq.combined.sct[["SCT"]]@scale.data, 
#                       scaled = TRUE, 
#                       gs = gs_list$gs_positive, 
#                       gs2 = gs_list$gs_negative,
#                       gene_names_to_uppercase = FALSE) 
# 
# 
# # es.max = sctype_score(scRNAseqData = burtoni.snseq.combined.sct[["SCT"]]@data %>% 
# #                             as.matrix(), 
# #                       scaled = TRUE, 
# #                       gs = gs_list$gs_positive, 
# #                       gs2 = gs_list$gs_negative,
# #                       gene_names_to_uppercase = FALSE) 
# 
# # NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
# # In case Seurat is used, it is either pbmc[["RNA"]]@scale.data (default), pbmc[["SCT"]]@scale.data, in case sctransform is used for normalization,
# # or pbmc[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.
# 
# # merge by cluster
# cL_results = do.call("rbind", lapply(unique(burtoni.snseq.combined.sct@meta.data$seurat_clusters), function(cl){
#   es.max.cl = sort(rowSums(es.max[ ,rownames(burtoni.snseq.combined.sct@meta.data[burtoni.snseq.combined.sct@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
#   head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(burtoni.snseq.combined.sct@meta.data$seurat_clusters==cl)), 10)
# }))
# sctype_scores = cL_results %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
# 
# # set low-confident (low ScType score) clusters to "unknown"
# sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
# print(sctype_scores[,1:3])
# 
# #graph umap
# burtoni.snseq.combined.sct@meta.data$sctypemarkers = ""
# for(j in unique(sctype_scores$cluster)){
#   cl_type = sctype_scores[sctype_scores$cluster==j,]; 
#   burtoni.snseq.combined.sct@meta.data$sctypemarkers[burtoni.snseq.combined.sct@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
# }
# 
# #graph umap
# DimPlot(burtoni.snseq.combined.sct,
#         reduction = "umap",
#         label = TRUE,
#         repel = TRUE,
#         group.by = 'sctypemarkers') +
#   ggtitle('sctype')
# ggsave('./sctype/UMAP cell types sctype custom.png',
#        height = 5.5,
#        width = 6.5)
# 
# # set low-confident (low ScType score) clusters to "unknown"
# sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
# print(sctype_scores[,1:3])
# 
# #graph umap
# burtoni.snseq.combined.sct@meta.data$sctypemarkers = ""
# for(j in unique(sctype_scores$cluster)){
#   cl_type = sctype_scores[sctype_scores$cluster==j,]; 
#   burtoni.snseq.combined.sct@meta.data$sctypemarkers[burtoni.snseq.combined.sct@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
# }
# 
# #graph umap
# DimPlot(burtoni.snseq.combined.sct,
#         reduction = "umap",
#         label = TRUE,
#         repel = TRUE,
#         group.by = 'sctypemarkers') +
#   ggtitle('sctype')
# ggsave('./sctype/UMAP cell types sctype custom confidence.png',
#        height = 5.5,
#        width = 6.5)
# 
# # DimPlot(burtoni.snseq.combined.sct, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'Cell.type')

# #### graph alluvial plot database ####
# ### compare sctype to seurat
# burtoni.snseq.combined.sct@meta.data %>% 
#   select(c(sctypemarkers,
#            seurat_clusters)) %>% 
#   table() %>% 
#   as.data.frame() %>% 
#   mutate(Count = sum(Freq)) %>% 
#   ungroup() %>% 
#   mutate(Freq.scale = 100*Freq/Count) %>% 
#   ggplot(aes(axis1 = reorder(seurat_clusters,-Freq.scale),
#              axis2 = reorder(sctypemarkers,-Freq.scale),
#              y = Freq.scale)) +
#   geom_alluvium(aes(fill = sctypemarkers)) +
#   geom_stratum() +
#   geom_text(stat = "stratum",
#             aes(label = after_stat(stratum))) +
#   scale_x_discrete(limits = c("seurat_clusters", 
#                               "sctypemarkers"),
#                    expand = c(0.15, 0.05)) +
#   scale_fill_viridis_d() +
#   theme_classic()
# ggsave('./sctype/Alluvial cells by celltype and cluster custom.png',
#        height = 10,
#        width = 10)
# 
# ### compare sctype to scsorter
# burtoni.snseq.combined.sct@meta.data %>% 
#   select(c(sctypemarkers,
#            Cell.type)) %>% 
#   table() %>% 
#   as.data.frame() %>% 
#   mutate(Count = sum(Freq)) %>% 
#   ungroup() %>% 
#   mutate(Freq.scale = 100*Freq/Count) %>% 
#   ggplot(aes(axis1 = reorder(Cell.type, -Freq.scale),
#              axis2 = reorder(sctypemarkers,-Freq.scale),
#              y = Freq.scale)) +
#   geom_alluvium(aes(fill = Cell.type)) +
#   geom_stratum() +
#   geom_text(stat = "stratum",
#             aes(label = after_stat(stratum))) +
#   scale_x_discrete(limits = c("Cell.type", 
#                               "sctypemarkers"),
#                    expand = c(0.15, 0.05)) +
#   scale_fill_viridis_d() +
#   theme_classic() 
# ggsave('./sctype/Alluvial cells by celltype custom.png',
#        height = 10,
#        width = 10)
# 
# ## color by broad cell types
# burtoni.snseq.combined.sct@meta.data %>% 
#   select(c(sctypemarkers,
#            Cell.type)) %>% 
#   table() %>% 
#   as.data.frame() %>% 
#   mutate(Count = sum(Freq)) %>% 
#   ungroup() %>% 
#   mutate(Freq.scale = 100*Freq/Count) %>% 
#   mutate(Cell.type.broad = case_when(Cell.type == "oligodendrocytes" ~ "glia",
#                                      Cell.type =="endothelial"~ "vascular",
#                                      Cell.type =="excitatory"~ "neuron",
#                                      Cell.type =="inhibitory"~ "neuron",
#                                      Cell.type =="neurons"~ "neuron",
#                                      Cell.type =="Unknown"~ "unknown",
#                                      Cell.type =="astrocytes"~ "glia",
#                                      Cell.type =="microglia"~ "glia",
#                                      Cell.type =="mural"~ "vascular",
#                                      Cell.type =="VSM"~ "vascular",
#                                      Cell.type =="OPCs"~ "glia",
#                                      Cell.type =="ependymal"~ "vascular",
#                                      TRUE ~ "NA")) %>% 
#   mutate(sctypemarkers.broad = case_when(sctypemarkers == "oligodendrocytes" ~ "glia",
#                                      sctypemarkers =="endothelial"~ "vascular",
#                                      sctypemarkers =="excitatory"~ "neuron",
#                                      sctypemarkers =="inhibitory"~ "neuron",
#                                      sctypemarkers =="neurons"~ "neuron",
#                                      sctypemarkers =="Unknown"~ "unknown",
#                                      sctypemarkers =="astrocytes"~ "glia",
#                                      sctypemarkers =="microglia"~ "glia",
#                                      sctypemarkers =="mural"~ "vascular",
#                                      sctypemarkers =="VSM"~ "vascular",
#                                      sctypemarkers =="OPCs"~ "glia",
#                                      sctypemarkers =="ependymal"~ "vascular",
#                                      TRUE ~ "NA")) %>% 
#   ggplot(aes(axis1 = reorder(Cell.type.broad, -Freq.scale),
#              axis2 = reorder(sctypemarkers.broad,-Freq.scale),
#              y = Freq.scale)) +
#   geom_alluvium(aes(fill = Cell.type.broad)) +
#   geom_stratum(aes(fill = Cell.type.broad)) +
#   geom_stratum(aes(fill = sctypemarkers.broad)) +
#   geom_text(color = 'white',
#             size = 5,
#             stat = "stratum",
#             aes(label = after_stat(stratum))) +
#   scale_x_discrete(limits = c("Cell.type.broad", 
#                               "sctypemarkers.broad"),
#                    expand = c(0.15, 0.05)) +
#   scale_fill_viridis_d() +
#   theme_classic() 
# ggsave('./sctype/Alluvial cells by celltype custom broad.png',
#        height = 5.5,
#        width = 7)
# 
# 
# 
# 
# 

# #### bubble network database ####
#' # load libraries
#' lapply(c("ggraph","igraph","tidyverse", "data.tree"), library, character.only = T)
#' 
#' # prepare edges
#' cL_results=cL_results[order(cL_results$cluster),]; edges = cL_results; edges$type = paste0(edges$type,"_",edges$cluster); edges$cluster = paste0("cluster ", edges$cluster); edges = edges[,c("cluster", "type")]; colnames(edges) = c("from", "to"); rownames(edges) <- NULL
#' 
#' # prepare nodes
#' nodes_lvl1 = sctype_scores[,c("cluster", "ncells")]; nodes_lvl1$cluster = paste0("cluster ", nodes_lvl1$cluster); nodes_lvl1$Colour = "#f1f1ef"; nodes_lvl1$ord = 1; nodes_lvl1$realname = nodes_lvl1$cluster; nodes_lvl1 = as.data.frame(nodes_lvl1); nodes_lvl2 = c(); 
#' 
#' ccolss=rainbow(32)
#' # ccolss= c("#5f75ae","#92bbb8","#64a841","#e5486e","#de8e06","#eccf5a","#b5aa0f","#e4b680","#7ba39d","#b15928","#ffff99", "#6a3d9a","#cab2d6","#ff7f00","#fdbf6f","#e31a1c","#fb9a99","#33a02c","#b2df8a","#1f78b4","#a6cee3")
#' for (i in 1:length(unique(cL_results$cluster))){
#'   dt_tmp = cL_results[cL_results$cluster == unique(cL_results$cluster)[i], ]; nodes_lvl2 = rbind(nodes_lvl2, data.frame(cluster = paste0(dt_tmp$type,"_",dt_tmp$cluster), ncells = dt_tmp$scores, Colour = ccolss[i], ord = 2, realname = dt_tmp$type))
#' }
#' 
#' nodes = rbind(nodes_lvl1, nodes_lvl2); nodes$ncells[nodes$ncells<1] = 1;
#' files_db = annotation.data[,c("cellName")]; files_db = unique(files_db); nodes = merge(nodes, files_db, all.x = T, all.y = F, by.x = "realname", by.y = "cellName", sort = F)
#' nodes$shortName[is.na(nodes$shortName)] = nodes$realname[is.na(nodes$shortName)]; nodes = nodes[,c("cluster", "ncells", "Colour", "ord", "shortName", "realname")]
#' 
#' nodes = rbind(nodes_lvl1, nodes_lvl2); nodes$ncells[nodes$ncells<1] = 1;
#' files_db = annotation.data %>% dplyr::select(cellName); files_db = unique(files_db); nodes = merge(nodes, files_db, all.x = T, all.y = F, by.x = "realname", by.y = "cellName", sort = F)
#' nodes$shortName[is.na(nodes$shortName)] = nodes$realname[is.na(nodes$shortName)]; nodes = nodes[,c("cluster", "ncells", "Colour", "ord", "realname")]
#' 
#' #' # need to remove some rows with shortName:
#' #' #'Immune system' and 'Endothelial'
#' #' mygraph <- graph_from_data_frame(edges, 
#' #'                                  vertices=nodes %>% 
#' #'                                    filter(shortName != 'Immune system',
#' #'                                           shortName != 'Endothelial'))
#' 
#' # Make the graph
#' gggr<- ggraph(mygraph, layout = 'circlepack', weight=I(ncells)) + 
#'   geom_node_circle(aes(filter=ord==1,fill=I("#F5F5F5"), colour=I("#D3D3D3")), alpha=0.9) + geom_node_circle(aes(filter=ord==2,fill=I(Colour), colour=I("#D3D3D3")), alpha=0.9) +
#'   theme_void() + geom_node_text(aes(filter=ord==2, label=realname, colour=I("#ffffff"), fill="white", repel = !1, parse = T, size = I(log(ncells,25)*1.5)))+ geom_node_label(aes(filter=ord==1,  label=realname, colour=I("#000000"), size = I(3), fill="white", parse = T), repel = !0, segment.linetype="dotted")
#' 
#' labelledUMAP = DimPlot(burtoni.snseq.combined.sct.all, reduction = "umap", label = TRUE, repel = TRUE, cols = ccolss)
#' 
#' cowplot::plot_grid(labelledUMAP,
#'                    gggr)
#' 
#' cowplot::plot_grid(gggr)
#' ggsave('sctype/UMAP cell types with network custom.png',
#'        width = 10,
#'        height = 10)
#' 
#' # No function scater but could just use cowplot
#' # scater::multiplot(DimPlot(burtoni.snseq.combined.sct.all, reduction = "umap", label = TRUE, repel = TRUE, cols = ccolss), gggr, cols = 2)



# #### graph cell counts ####
# burtoni.snseq.combined.sct@meta.data %>%
#   mutate(orig.ident = ifelse(orig.ident == "dom_burtoni_snseq",
#                              "Dominant",
#                              "Subordinate")) %>% 
#   mutate(Count = 1) %>% 
#   group_by(Genotype.id,
#            orig.ident) %>% 
#   summarise(Total = sum(Count)) %>%
#   ggplot(aes(y = Total,
#              x = orig.ident,
#              fill = orig.ident)) + 
#   geom_point(size = 8,
#              shape = 21) +
#   scale_fill_manual(values = c('Dominant' = '#4e499e',
#                                'Subordinate' = '#60bb46')) +
#   theme_classic() +
#   xlab('') +
#   ylab('Total cells assigned')+
#   theme_classic() +
#   # theme(panel.border = element_rect(color = "black",
#   #                                   fill = NA,
#   #                                   size = 1))+
#   theme(axis.text = element_text(size = 15,
#                                  angle = 45,
#                                  hjust = 1))  +
#   theme(axis.title = element_text(size = 20))+
#   theme(legend.position = 'none',
#         legend.box.background = element_rect(colour = "black"),
#         legend.background = element_blank(),
#         legend.text=element_text(size=20),
#         legend.title = element_blank()) +
#   ylim(c(0,6000))
# ggsave('../souporcell/figures/filtered genotype id by status cell count presentation.png',
#        height = 5.5,
#        width = 5.5,
#        units = "in",
#        dpi = 300)
# 
# ### all cells
# ## status, cell type, cluster
# # color celltype
# burtoni.snseq.combined.sct@meta.data %>% 
#   select(c(orig.ident,
#            Cell.type,
#            seurat_clusters)) %>% 
#   table() %>% 
#   as.data.frame() %>% 
#   group_by(orig.ident) %>% 
#   mutate(Count = sum(Freq)) %>% 
#   ungroup() %>% 
#   mutate(Freq.scale = 100*Freq/Count) %>% 
#   ggplot(aes(axis1 = orig.ident, 
#              axis2 = reorder(Cell.type, -Freq.scale),
#              axis3 = reorder(seurat_clusters,-Freq.scale),
#              y = Freq.scale)) +
#   geom_alluvium(aes(fill = Cell.type)) +
#   geom_stratum() +
#   geom_text(stat = "stratum",
#             aes(label = after_stat(stratum))) +
#   scale_x_discrete(limits = c("orig.ident", 
#                               "seurat_clusters"),
#                    expand = c(0.15, 0.05)) +
#   scale_fill_viridis_d() +
#   theme_classic() +
#   theme(axis.text.x=element_blank())
# ggsave('./neuropeptides/avp.oxt/avp/alluvial/All cells by celltype.png',
#        height = 10,
#        width = 10)
# 
# ### graph alluvial
# # broad cell type
# # color broad cell
# burtoni.snseq.combined.sct@meta.data %>% 
#   select(c(orig.ident,
#            Cell.type,
#            seurat_clusters)) %>% 
#   mutate(Cell.type.broad = case_when(Cell.type == "oligodendrocytes" ~ "glia",
#                                      Cell.type =="endothelial"~ "vascular",
#                                      Cell.type =="excitatory"~ "neuron",
#                                      Cell.type =="inhibitory"~ "neuron",
#                                      Cell.type =="neurons"~ "neuron",
#                                      Cell.type =="Unknown"~ "unknown",
#                                      Cell.type =="astrocytes"~ "glia",
#                                      Cell.type =="microglia"~ "glia",
#                                      Cell.type =="mural"~ "vascular",
#                                      Cell.type =="VSM"~ "vascular",
#                                      Cell.type =="OPCs"~ "glia",
#                                      Cell.type =="ependymal"~ "vascular",
#                                      TRUE ~ "NA")) %>% 
#   mutate(orig.ident = ifelse(orig.ident == 'dom_burtoni_snseq',
#                              'Dom',
#                              'Sub')) %>% 
#   select(-c(Cell.type)) %>% 
#   table() %>% 
#   as.data.frame() %>% 
#   group_by(orig.ident) %>% 
#   mutate(Count = sum(Freq)) %>% 
#   ungroup() %>% 
#   mutate(Freq.scale = 100*Freq/Count) %>% 
#   ggplot(aes(axis1 = orig.ident, 
#              axis2 = reorder(Cell.type.broad, -Freq.scale),
#              axis3 = reorder(seurat_clusters,-Freq.scale),
#              y = Freq.scale)) +
#   geom_alluvium(aes(fill = Cell.type.broad)) +
#   geom_stratum(aes(fill = Cell.type.broad)) +
#   geom_text(color = 'white',
#             size = 5,
#             stat = "stratum",
#             aes(label = after_stat(stratum))) +
#   scale_x_discrete(limits = c("orig.ident", 
#                               "seurat_clusters"),
#                    expand = c(0.15, 0.05)) +
#   scale_fill_viridis_d() +
#   theme_classic() +
#   theme(axis.text.x=element_blank()) 
# ggsave('./sctype/All cells by broad cell type.png',
#        height = 5.5,
#        width = 10)
# 

#### run sctype HypoMap ####
#### cell type assignment
### prepare gene set list
# set tissue type
tissue = "Brain"
# prepare gene sets
gs_list = gene_sets_prepare.df(annotation.data.hypomap,
                               tissue)

### assign cell types to clusters
# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = burtoni.snseq.combined.sct[["SCT"]]@scale.data,
                      scaled = TRUE,
                      gs = gs_list$gs_positive,
                      gs2 = gs_list$gs_negative,
                      gene_names_to_uppercase = FALSE)


# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix.
# In case Seurat is used, it is either pbmc[["RNA"]]@scale.data (default), pbmc[["SCT"]]@scale.data, in case sctransform is used for normalization,
# or pbmc[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.


# merge by cluster
cL_results = do.call("rbind", lapply(unique(burtoni.snseq.combined.sct@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(burtoni.snseq.combined.sct@meta.data[burtoni.snseq.combined.sct@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(burtoni.snseq.combined.sct@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_results %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

#graph umap
burtoni.snseq.combined.sct@meta.data$sctypemarkers.hypo = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,];
  burtoni.snseq.combined.sct@meta.data$sctypemarkers.hypo[burtoni.snseq.combined.sct@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

#graph umap
DimPlot(burtoni.snseq.combined.sct,
        reduction = "umap",
        label = TRUE,
        repel = TRUE,
        group.by = 'sctypemarkers.hypo') +
  ggtitle('sctype')
ggsave('./sctype.hypo/UMAP cell types sctype custom.png',
       height = 5.5,
       width = 6.5)

## graph counts
# drop zeros
burtoni.snseq.combined.sct@meta.data %>%
  dplyr::select(c(Genotype.id,
                  sctypemarkers.hypo,
                  orig.ident)) %>%
  table() %>%
  as.data.frame() %>%
  filter(Freq > 0) %>%
  ggplot(aes(x = sctypemarkers.hypo,
             y = Freq,
             color = orig.ident,
             group = Genotype.id)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave('./sctype.hypo/Count cell types sctype.png',
       height = 5.5,
       width = 6.5)

### add predicted cell counts
## create predicted data frame
predicted.cell.counts = data.frame(orig.ident = 'predicted',
                                   sctypemarkers.hypo = c('C7-1: GLU',
                                                        'C7-2: GABA',
                                                        'C7-3: Astro-Ependymal',
                                                        'C7-4: Oligo+Precursor',
                                                        'C7-5: Immune',
                                                        'C7-6: ParsTuber',
                                                        'C7-7: Vascular'),
                                      percent = c(24.01,
                                                  32.89,
                                                  16.86,
                                                  18.16,
                                                  3.78,
                                                  0.19,
                                                  4.11))

## convert counts to percent
sctype.percent = burtoni.snseq.combined.sct@meta.data %>%
  dplyr::select(c(Genotype.id,
                  sctypemarkers.hypo,
                  orig.ident)) %>%
  table() %>%
  as.data.frame() %>%
  filter(Freq > 0) %>%
  group_by(Genotype.id) %>%
  mutate(Total = sum(Freq)) %>%
  ungroup() %>%
  mutate(percent = 100*Freq/Total)


## graph counts
sctype.percent %>%
  full_join(predicted.cell.counts) %>%
  ggplot(aes(x = sctypemarkers.hypo,
             y = percent,
             color = orig.ident,
             group = Genotype.id)) +
  geom_line() +
  geom_point() +
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave('./sctype.hypo/Percentage cell types sctype predicted.png',
       height = 5.5,
       width = 6.5)



#sctype score per cluster
for (i in unique(cL_results$cluster)) {
  cL_results %>%
    mutate(type.call = ncells/4,
           color.in = ifelse(scores > type.call,
                             'above',
                             'below')) %>%
    filter(cluster == i) %>%
    ggplot(aes(x = reorder(type, -100*scores/ncells),
               y = 100*scores/ncells,
               fill = color.in)) +
    geom_col() +
    theme_classic() +
    xlab('') +
    ylab('Scaled cell type score') +
    labs(fill = "Threshold: ncells/4") +
    ggtitle(paste('Cluster ',
                  i,
                  sep = ''))+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  ggsave(paste('./sctype.hypo/clusters/Cluster ',
               i,
               ' score.png'),
         height = 5,
         width = 10)
}

## graph max value for each cluster
cL_results %>%
  mutate(type.call = ncells/4,
         color.in = ifelse(scores > type.call,
                           'above',
                           'below')) %>%
  full_join(burtoni.snseq.combined.sct@meta.data %>%
              dplyr::select(seurat_clusters,
                            sctypemarkers.hypo) %>%
              distinct(),
            by = c("cluster" = "seurat_clusters")) %>%
  group_by(cluster,
           sctypemarkers.hypo,
           ncells) %>%
  summarise(max = max(scores)) %>%
  mutate(percentage = 100*max/ncells) %>%
  ggplot(aes(color = sctypemarkers.hypo,
             y = percentage,
             x= reorder(cluster,
                        -percentage))) +
  geom_point() +
  theme_classic() +
  geom_hline(yintercept = 25,
             linetype = 'dashed')
ggsave('./sctype.hypo/sctype score percentage per cluster.png',
       height = 5,
       width = 10)

## create neuron score
cL_results %>%
  mutate(type.call = ncells/4,
         color.in = ifelse(scores > type.call,
                           'above',
                           'below')) %>%
  full_join(burtoni.snseq.combined.sct@meta.data %>%
              dplyr::select(seurat_clusters,
                            sctypemarkers.hypo) %>%
              distinct(),
            by = c("cluster" = "seurat_clusters")) %>%
  view()

### call indvidual cells
## histogram
es.max %>%
  as.data.frame() %>%
  rownames_to_column('type') %>%
  pivot_longer(!type,
               names_to = 'cell',
               values_to = 'sctype.score') %>%
  filter(sctype.score > 0) %>%
  ggplot() +
  geom_histogram(aes(sctype.score,
                     fill = sctype.score > 1)) +
  theme_bw() +
  facet_wrap(~type)
ggsave('./sctype.hypo/Histogram sctype score per cell.png',
       height = 10,
       width = 10)

## Call individual cell types
# only include sctype score above 1
# select top sctype score per cell to assign type
# calculate percentage of sctype score, needs to be above 50%
es.max.cl = es.max %>%
  as.data.frame() %>%
  rownames_to_column('sctype.hypo') %>%
  pivot_longer(!sctype.hypo,
               names_to = 'cell.id',
               values_to = 'sctype.hypo.score') %>%
  mutate(sctype.hypo.score.reduce = ifelse(sctype.hypo.score > 1,
                                           sctype.hypo.score,
                                           0)) %>%
  arrange(desc(sctype.hypo.score.reduce)) %>%
  group_by(cell.id) %>%
  mutate(sctype.hypo.score.reduce.max = max(sctype.hypo.score.reduce),
         sctype.hypo.score.reduce.sum = sum(sctype.hypo.score.reduce)) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(sctype.hypo.score.reduce.percent = 100*sctype.hypo.score.reduce/sctype.hypo.score.reduce.sum,
         sctype.hypo.score.reduce.percent = ifelse(is.na(sctype.hypo.score.reduce.percent),
                                                   0,
                                                   sctype.hypo.score.reduce.percent),
         sctype.hypo.score.reduce.percent.group = case_when(sctype.hypo.score.reduce.percent >= 90 ~ '~100%',
                                                            sctype.hypo.score.reduce.percent >= 50 ~ '> 50%',
                                                            sctype.hypo.score.reduce.percent >= 0 ~ '< 50%',
                                                            TRUE ~ 'NA'),
         sctype.hypo.update = ifelse(sctype.hypo.score.reduce.percent.group == '< 50%',
                                     'Unknown',
                                     sctype.hypo))


# graph
es.max.cl %>%
  dplyr::select(sctype.hypo,
                sctype.hypo.update,
                sctype.hypo.score.reduce.percent.group) %>%
  table() %>%
  as.data.frame() %>%
  filter(Freq > 0) %>%
  ggplot() +
  geom_point(aes(x = sctype.hypo,
                 y = Freq,
                 group = sctype.hypo.update,
                 color = sctype.hypo.score.reduce.percent.group)) +
  theme_bw()
ggsave('./sctype.hypo/Sctype score per cell.png',
       height = 10,
       width = 10)


### compare results of cell to cluster sctype hypomap
burtoni.snseq.combined.sct@meta.data %>%
  rownames_to_column('cell.id') %>%
  full_join(es.max.cl %>%
              select(c(cell.id,
                       sctype.hypo,
                       sctype.hypo.score.reduce.percent.group,
                       sctype.hypo.update))) %>%
  select(c(sctypemarkers.hypo,
           sctype.hypo.update)) %>%
  table() %>%
  as.data.frame() %>%
  mutate(Count = sum(Freq)) %>%
  ungroup() %>%
  mutate(Freq.scale = 100*Freq/Count) %>%
  ggplot(aes(axis1 = reorder(sctypemarkers.hypo, -Freq.scale),
             axis2 = reorder(sctype.hypo.update,-Freq.scale),
             y = Freq.scale)) +
  geom_alluvium(aes(fill = sctype.hypo.update)) +
  geom_stratum() +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("sctypemarkers.hypo",
                              "sctype.hypo.update"),
                   expand = c(0.15, 0.05)) +
  scale_fill_viridis_d() +
  theme_classic()
ggsave('./sctype.hypo/Alluvial cells by sctype.hypo cluster vs individual.png',
       height = 10,
       width = 10)

# graph UMAP
burtoni.snseq.combined.sct@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  rownames_to_column("cell.id") %>%
  full_join(es.max.cl %>%
              select(c(cell.id,
                       sctype.hypo,
                       sctype.hypo.score.reduce.percent.group,
                       sctype.hypo.update))) %>%
  ggplot() +
  geom_point(aes(x = UMAP_1,
             y = UMAP_2,
             color = sctype.hypo)) +
  theme_classic()
ggsave('./sctype.hypo/UMAP sctype hypo individual cell.png',
       height = 10,
       width = 10)

# graph UMAP
burtoni.snseq.combined.sct@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  rownames_to_column("cell.id") %>%
  full_join(es.max.cl %>%
              select(c(cell.id,
                       sctype.hypo,
                       sctype.hypo.score.reduce.percent.group,
                       sctype.hypo.update))) %>%
  ggplot() +
  geom_point(aes(x = UMAP_1,
                 y = UMAP_2,
                 color = sctype.hypo.update)) +
  theme_classic()
ggsave('./sctype.hypo/UMAP sctype hypo individual cell unknown.png',
       height = 10,
       width = 10)

# rename to broad cell types
burtoni.snseq.combined.sct@meta.data = burtoni.snseq.combined.sct@meta.data %>%
  separate(sctypemarkers.hypo,
           c("Cell.type.level",
             "sctypemarkers.hypo.broad"),
           remove = F,
           sep = ':') %>%
  mutate(sctypemarkers.hypo.broad = ifelse(is.na(sctypemarkers.hypo.broad),
                                           sctypemarkers.hypo,
                                           sctypemarkers.hypo.broad))

#### save metadata ####
### save sctype hypo cell types
# extract metadata to save
sctypemarkers.hypo = burtoni.snseq.combined.sct@meta.data %>%
  dplyr::select(sctypemarkers.hypo,
                Cell.type.level,
                sctypemarkers.hypo.broad) %>%
  rownames_to_column('Cell.id')


#save file
write.csv(sctypemarkers.hypo,
          file = './sctype.hypo/sctypemarkers.hypo.csv')


### load sctype.hypo data
burtoni.sctypemarkers.hypo = read.csv('./sctype.hypo/sctypemarkers.hypo.csv')

## sctype.hypo
burtoni.snseq.combined.sct = AddMetaData(
  object = burtoni.snseq.combined.sct,
  metadata = burtoni.sctypemarkers.hypo %>%
    select(sctypemarkers.hypo,
           Cell.id) %>%
    column_to_rownames(var = "Cell.id"),
  col.name = 'sctypemarkers.hypo'
)

## sctype.hypo
burtoni.snseq.combined.sct = AddMetaData(
  object = burtoni.snseq.combined.sct,
  metadata = burtoni.sctypemarkers.hypo %>%
    select(sctypemarkers.hypo.broad,
           Cell.id) %>%
    column_to_rownames(var = "Cell.id"),
  col.name = 'sctypemarkers.hypo.broad'
)


# #### cell counts and umap presentation ####
# ### cell counts
# # burtoni.snseq.combined.sct@meta.data %>% 
# #   dplyr::select(sctypemarkers.hypo,
# #                 Genotype.id,
# #                 orig.ident) %>% 
# #   table() %>% 
# #   as.data.frame() %>% 
# #   filter(Freq > 0) %>% 
# #   ggplot(aes(x = sctypemarkers.hypo,
# #              y = Freq,
# #              fill = orig.ident)) +
# #   geom_point(size = 8,
# #              shape = 21) +
# #   scale_fill_manual(values = c("#60bb46", "#4e499e")) +
# #   theme_classic() + 
# #   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
# #         text = element_text(size = 20),
# #         legend.position = 'none') +
# #   ylab("Total cells assigned") +
# #   xlab('') +
# #   labs('')
# # ggsave('sctype.hypo/Count cell per status.png',
# #        height = 6,
# #        width = 6)
# # 
# # # broad
# # burtoni.snseq.combined.sct@meta.data %>% 
# #   dplyr::select(sctypemarkers.hypo.broad,
# #                 Genotype.id,
# #                 orig.ident) %>% 
# #   table() %>% 
# #   as.data.frame() %>% 
# #   filter(Freq > 0) %>% 
# #   mutate(order = Freq,
# #          order = ifelse(sctypemarkers.hypo.broad == "Unknown",
# #                         0,
# #                         order)) %>% 
# #   ggplot(aes(x = reorder(sctypemarkers.hypo.broad,
# #                          -order),
# #              y = Freq,
# #              fill = orig.ident)) +
# #   geom_point(size = 8,
# #              shape = 21) +
# #   scale_fill_manual(values = c("#60bb46", "#4e499e")) +
# #   theme_classic() + 
# #   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
# #         text = element_text(size = 20),
# #         legend.position = 'none') +
# #   ylab("Total cells assigned") +
# #   xlab('') +
# #   labs('')
# # ggsave('sctype.hypo/Count cell broad per status.png',
# #        height = 6,
# #        width = 6)
# ## UMAP for paper
# # just cell labels
# DimPlot(burtoni.snseq.combined.sct,
#         group.by = 'sctypemarkers.hypo.broad',
#         alpha = 0,
#         label = T,
#         label.box = T,
#         cols = c('#F8766D',
#                  '#B79F00',
#                  '#00BA38',
#                  '#00BFC4',
#                  '#619CFF',
#                  '#F564E3',
#                  'grey')) +
#   NoLegend() +
#   theme(
#     plot.title = element_blank()
#   )
# ggsave('sctype.hypo/UMAP cell broad paper labels.png',
#        height = 7,
#        width = 7)
# 
# # just cluster labels
# DimPlot(burtoni.snseq.combined.sct,
#         group.by = 'seurat_clusters',
#         alpha = 0,
#         label = T) +
#   NoLegend() +
#   theme(
#     plot.title = element_blank()
#   )
# ggsave('sctype.hypo/UMAP cell broad paper clusters.png',
#        height = 7,
#        width = 7)
# 
# # just color cells
# DimPlot(burtoni.snseq.combined.sct,
#         group.by = 'sctypemarkers.hypo.broad',
#         label = F,
#         cols = c('#F8766D',
#                  '#B79F00',
#                  '#00BA38',
#                  '#00BFC4',
#                  '#619CFF',
#                  '#F564E3',
#                  'grey')) +
#   NoLegend() +
#   theme(
#     plot.title = element_blank())
# ggsave('sctype.hypo/UMAP cell broad paper.png',
#        height = 7,
#        width = 7)
# # ## UMAP for presentation
# # # cell type broad
# # DimPlot(burtoni.snseq.combined.sct,
# #         group.by = 'sctypemarkers.hypo.broad',
# #         label = T) +
# #   NoLegend() +
# #   theme(
# #     plot.title = element_blank(),
# #     axis.title.x = element_blank(),
# #     axis.title.y = element_blank())
# # ggsave('sctype.hypo/UMAP cell broad presentation.png',
# #        height = 5.5,
# #        width = 5.5)
# # 
# # # grey umap
# # DimPlot(burtoni.snseq.combined.sct,
# #         group.by = 'orig.ident',
# #         cols = c('grey',
# #                  'grey')) +
# #   NoLegend()+
# #   theme(
# #     plot.title = element_blank(),
# #     axis.title.x = element_blank(),
# #     axis.title.y = element_blank())
# # ggsave('sctype.hypo/UMAP grey presentation.png',
# #        height = 5.5,
# #        width = 5.5)
# # 
# # # clusters umap
# # DimPlot(burtoni.snseq.combined.sct,
# #         group.by = 'seurat_clusters',
# #         label = T) +
# #   NoLegend() +
# #   theme(
# #     plot.title = element_blank(),
# #     axis.title.x = element_blank(),
# #     axis.title.y = element_blank())
# # ggsave('sctype.hypo/UMAP clusters presentation.png',
# #        height = 5.5,
# #        width = 5.5)
# 
# 

# #### cell type percentages ####
# ### percentages by cluster
# ## get counts of nuclei per cluster per genotype
# # create table of samples by cluster
# table_samples_by_clusters <- burtoni.snseq.combined.sct@meta.data %>% 
#   group_by(Genotype.id, 
#            seurat_clusters) %>%
#   summarize(count = n()) %>%
#   spread(seurat_clusters, 
#          count, 
#          fill = 0) %>%
#   ungroup() %>%
#   mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
#   dplyr::select(c('Genotype.id', 
#                   'total_cell_count', 
#                   everything())) %>%
#   arrange(factor(Genotype.id, levels = levels(burtoni.snseq.combined.sct@meta.data$Genotype.id)))
# 
# # create table of clusters by sample
# table_clusters_by_samples <- burtoni.snseq.combined.sct@meta.data %>%
#   group_by(Genotype.id, 
#            seurat_clusters) %>%
#   summarize(count = n()) %>%
#   spread(Genotype.id,
#          count, 
#          fill = 0) %>%
#   ungroup() %>%
#   mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
#   select(c('seurat_clusters', 
#            'total_cell_count', 
#            everything())) %>%
#   arrange(factor(seurat_clusters, levels = levels(burtoni.snseq.combined.sct@meta.data$seurat_clusters)))
# 
# ## graph percentages by cluster
# # get count per sample
# temp_labels <- burtoni.snseq.combined.sct@meta.data %>%
#   group_by(Genotype.id) %>%
#   tally()
# 
# # create graph of sample by cluster
# p1 <- table_samples_by_clusters %>%
#   dplyr::select(-c('total_cell_count')) %>%
#   reshape2::melt(id.vars = 'Genotype.id') %>%
#   # mutate(Genotype.id = factor(Genotype.id, levels = levels(burtoni.snseq.combined.sct@meta.data$Genotype.id))) %>%
#   ggplot(aes(Genotype.id, value)) +
#   geom_bar(aes(fill = variable), position = 'fill', stat = 'identity') +
#   geom_text(
#     data = temp_labels,
#     aes(x = Genotype.id, y = Inf, label = paste0(format(n, big.mark = ',', trim = TRUE)), vjust = -1),
#     color = 'black', size = 2.8
#   ) +
#   # scale_fill_manual(name = 'Cluster') +
#   scale_y_continuous(name = 'Percentage [%]', labels = scales::percent_format(), expand = c(0.01,0)) +
#   coord_cartesian(clip = 'off') +
#   theme_bw() +
#   theme(
#     legend.position = 'left',
#     plot.title = element_text(hjust = 0.5),
#     text = element_text(size = 16),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     axis.title.x = element_blank(),
#     axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
#     plot.margin = margin(t = 20, r = 0, b = 0, l = 0, unit = 'pt')
#   )
# 
# # get count per cluster
# temp_labels <- burtoni.snseq.combined.sct@meta.data %>%
#   group_by(seurat_clusters) %>%
#   tally()
# 
# # create graph of cluster by sample
# p2 <- table_clusters_by_samples %>%
#   select(-c('total_cell_count')) %>%
#   reshape2::melt(id.vars = 'seurat_clusters') %>%
#   mutate(seurat_clusters = factor(seurat_clusters, levels = levels(burtoni.snseq.combined.sct@meta.data$seurat_clusters))) %>%
#   ggplot(aes(seurat_clusters, value)) +
#   geom_bar(aes(fill = variable), position = 'fill', stat = 'identity', color = 'black') +
#   geom_text(
#     data = temp_labels, aes(x = seurat_clusters, y = Inf, label = paste0(format(n, big.mark = ',', trim = TRUE)), vjust = -1),
#     color = 'black', size = 2.8
#   ) +
#   scale_fill_manual(name = 'Sample', values = c("#60bb46", 
#                                                 "#60bb46", 
#                                                 "#60bb46", 
#                                                 "#4e499e",
#                                                 "#4e499e",
#                                                 "#4e499e"
#                                                 )) +
#   scale_y_continuous(name = 'Percentage [%]', labels = scales::percent_format(), expand = c(0.01,0)) +
#   coord_cartesian(clip = 'off') +
#   theme_bw() +
#   theme(
#     legend.position = 'right',
#     plot.title = element_text(hjust = 0.5),
#     text = element_text(size = 16),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     axis.title = element_blank(),
#     plot.margin = margin(t = 20, r = 0, b = 0, l = 10, unit = 'pt')
#   )
# 
# ggsave(
#   'sctype.hypo/Percentage genotype by clusters.png',
#   p1 + p2 +
#     patchwork::plot_layout(ncol = 2, widths = c(
#       burtoni.snseq.combined.sct@meta.data$Genotype.id %>% unique() %>% length(),
#       burtoni.snseq.combined.sct@meta.data$seurat_clusters %>% unique() %>% length()
#     )),
#   width = 18, height = 8
# )
# 
# ### percentages by celltype
# ## get counts of nuclei per celltype per genotype
# # create table of samples by celltype
# table_samples_by_celltype <- burtoni.snseq.combined.sct@meta.data %>% 
#   group_by(Genotype.id, 
#            sctypemarkers.hypo.broad) %>%
#   summarize(count = n()) %>%
#   spread(sctypemarkers.hypo.broad, 
#          count, 
#          fill = 0) %>%
#   ungroup() %>%
#   mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
#   dplyr::select(c('Genotype.id', 
#                   'total_cell_count', 
#                   everything())) %>%
#   arrange(factor(Genotype.id, levels = levels(burtoni.snseq.combined.sct@meta.data$Genotype.id)))
# 
# # create table of celltype by sample
# table_celltype_by_samples <- burtoni.snseq.combined.sct@meta.data %>%
#   group_by(Genotype.id, 
#            sctypemarkers.hypo.broad) %>%
#   summarize(count = n()) %>%
#   spread(Genotype.id,
#          count, 
#          fill = 0) %>%
#   ungroup() %>%
#   mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
#   select(c('sctypemarkers.hypo.broad', 
#            'total_cell_count', 
#            everything())) %>%
#   arrange(factor(sctypemarkers.hypo.broad, levels = levels(burtoni.snseq.combined.sct@meta.data$sctypemarkers.hypo.broad)))
# 
# ## graph percentages by celltype
# # get count per sample
# temp_labels <- burtoni.snseq.combined.sct@meta.data %>%
#   group_by(Genotype.id) %>%
#   tally()
# 
# # create graph of sample by celltype
# p1 <- table_samples_by_celltype %>%
#   dplyr::select(-c('total_cell_count')) %>%
#   reshape2::melt(id.vars = 'Genotype.id') %>%
#   # mutate(Genotype.id = factor(Genotype.id, levels = levels(burtoni.snseq.combined.sct@meta.data$Genotype.id))) %>%
#   ggplot(aes(Genotype.id, value)) +
#   geom_bar(aes(fill = variable), position = 'fill', stat = 'identity') +
#   geom_text(
#     data = temp_labels,
#     aes(x = Genotype.id, y = Inf, label = paste0(format(n, big.mark = ',', trim = TRUE)), vjust = -1),
#     color = 'black', size = 2.8
#   ) +
#   # scale_fill_manual(name = 'Cluster') +
#   scale_y_continuous(name = 'Percentage [%]', labels = scales::percent_format(), expand = c(0.01,0)) +
#   coord_cartesian(clip = 'off') +
#   theme_bw() +
#   theme(
#     legend.position = 'left',
#     plot.title = element_text(hjust = 0.5),
#     text = element_text(size = 16),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     axis.title.x = element_blank(),
#     axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
#     plot.margin = margin(t = 20, r = 0, b = 0, l = 0, unit = 'pt')
#   )
# 
# # get count per celltype
# temp_labels <- burtoni.snseq.combined.sct@meta.data %>%
#   group_by(sctypemarkers.hypo.broad) %>%
#   tally()
# 
# # create graph of celltype by sample
# p2 <- table_celltype_by_samples %>%
#   select(-c('total_cell_count')) %>%
#   reshape2::melt(id.vars = 'sctypemarkers.hypo.broad') %>% 
#   # mutate(sctypemarkers.hypo.broad = factor(sctypemarkers.hypo.broad, levels = levels(burtoni.snseq.combined.sct@meta.data$sctypemarkers.hypo.broad))) %>%
#   ggplot(aes(sctypemarkers.hypo.broad, value)) +
#   geom_bar(aes(fill = variable), position = 'fill', stat = 'identity', color = 'black') +
#   geom_text(
#     data = temp_labels, aes(x = sctypemarkers.hypo.broad, y = Inf, label = paste0(format(n, big.mark = ',', trim = TRUE)), vjust = -1),
#     color = 'black', size = 2.8
#   ) +
#   scale_fill_manual(name = 'Sample', values = c("#60bb46", 
#                                                 "#60bb46", 
#                                                 "#60bb46", 
#                                                 "#4e499e",
#                                                 "#4e499e",
#                                                 "#4e499e"
#   )) +
#   scale_y_continuous(name = 'Percentage [%]', labels = scales::percent_format(), expand = c(0.01,0)) +
#   coord_cartesian(clip = 'off') +
#   theme_bw() +
#   theme(
#     legend.position = 'right',
#     plot.title = element_text(hjust = 0.5),
#     text = element_text(size = 16),
#     panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(),
#     axis.title = element_blank(),
#     plot.margin = margin(t = 20, r = 0, b = 0, l = 10, unit = 'pt'),
#     axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
#   )
# 
# ggsave(
#   'sctype.hypo/Percentage genotype by celltype.png',
#   p1 + p2 +
#     patchwork::plot_layout(ncol = 2, widths = c(
#       burtoni.snseq.combined.sct@meta.data$Genotype.id %>% unique() %>% length(),
#       burtoni.snseq.combined.sct@meta.data$sctypemarkers.hypo.broad %>% unique() %>% length()
#     )),
#   width = 18, height = 8
# )
# 
# #### plot cell types tree ####
# # https://romanhaa.github.io/projects/scrnaseq_workflow/#clustering
# 
# ### create tree of cell types
# burtoni.snseq.combined.sct.tree <- BuildClusterTree(
#   burtoni.snseq.combined.sct,
#   dims = 1:15,
#   reorder = FALSE,
#   reorder.numeric = FALSE,
#   slot = 'data',
#   assay = "SCT"
# )
# 
# ## create tree
# tree <- burtoni.snseq.combined.sct.tree@tools$BuildClusterTree
# # tree$tip.label <- paste0("Cluster ", tree$tip.label)
# 
# # convert tree to tibble
# tree = as_tibble(tree)
# 
# ## add cell type data
# tree = full_join(tree,
#                  burtoni.snseq.combined.sct@meta.data %>% 
#                    dplyr::select(seurat_clusters,
#                                  sctypemarkers.hypo.broad) %>% 
#                    distinct(),
#                        by = c("label" = "seurat_clusters"))
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
#   ggtree::geom_tippoint(aes(color = sctypemarkers.hypo.broad),
#                         shape = 16,
#                         size = 5) +
#   coord_cartesian(clip = 'off') +
#   theme(plot.margin = unit(c(0,2.5,0,0),
#                            'cm'))
# ggsave('sctype.hypo/Cluster cell type tree.png',
#        height = 6,
#        width = 6)
# 
# 
# #### graph alluvial sctype broad vs hypomap ####
# 
# 
# ### compare sctype of hypomap markers to clusters
# burtoni.snseq.combined.sct@meta.data %>%
#   select(c(sctypemarkers.hypo,
#            seurat_clusters)) %>%
#   table() %>%
#   as.data.frame() %>%
#   mutate(Count = sum(Freq)) %>%
#   ungroup() %>%
#   mutate(Freq.scale = 100*Freq/Count) %>%
#   ggplot(aes(axis1 = reorder(sctypemarkers.hypo, -Freq.scale),
#              axis2 = reorder(seurat_clusters,-Freq.scale),
#              y = Freq.scale)) +
#   geom_alluvium(aes(fill = sctypemarkers.hypo)) +
#   geom_stratum() +
#   geom_text(stat = "stratum",
#             aes(label = after_stat(stratum))) +
#   scale_x_discrete(limits = c("sctypemarkers.hypo",
#                               "seurat_clusters"),
#                    expand = c(0.15, 0.05)) +
#   scale_fill_viridis_d() +
#   theme_classic()
# ggsave('./sctype/Alluvial cells by sctype.hypo vs cluster.png',
#        height = 10,
#        width = 10)
# 
# ### compare sctype of hypomap markers to published markers
# burtoni.snseq.combined.sct@meta.data %>%
#   select(c(sctypemarkers,
#            sctypemarkers.hypo)) %>%
#   table() %>%
#   as.data.frame() %>%
#   mutate(Count = sum(Freq)) %>%
#   ungroup() %>%
#   mutate(Freq.scale = 100*Freq/Count) %>%
#   ggplot(aes(axis1 = reorder(sctypemarkers, -Freq.scale),
#              axis2 = reorder(sctypemarkers.hypo,-Freq.scale),
#              y = Freq.scale)) +
#   geom_alluvium(aes(fill = sctypemarkers.hypo)) +
#   geom_stratum() +
#   geom_text(stat = "stratum",
#             aes(label = after_stat(stratum))) +
#   scale_x_discrete(limits = c("sctypemarkers",
#                               "sctypemarkers.hypo"),
#                    expand = c(0.15, 0.05)) +
#   scale_fill_viridis_d() +
#   theme_classic()
# ggsave('./sctype/Alluvial cells by sctype vs sctype.hypo.png',
#        height = 10,
#        width = 10)
# 
# ### compare sctype of hypomap markers to scsorter published markers
# burtoni.snseq.combined.sct@meta.data %>%
#   select(c(sctypemarkers.hypo,
#            Cell.type)) %>%
#   table() %>%
#   as.data.frame() %>%
#   mutate(Count = sum(Freq)) %>%
#   ungroup() %>%
#   mutate(Freq.scale = 100*Freq/Count) %>%
#   ggplot(aes(axis1 = reorder(sctypemarkers.hypo, -Freq.scale),
#              axis2 = reorder(Cell.type,-Freq.scale),
#              y = Freq.scale)) +
#   geom_alluvium(aes(fill = sctypemarkers.hypo)) +
#   geom_stratum() +
#   geom_text(stat = "stratum",
#             aes(label = after_stat(stratum))) +
#   scale_x_discrete(limits = c("sctypemarkers.hypo",
#                               "Cell.type"),
#                    expand = c(0.15, 0.05)) +
#   scale_fill_viridis_d() +
#   theme_classic()
# ggsave('./sctype/Alluvial cells by sctype.hypo vs scsorter.png',
#        height = 10,
#        width = 10)
# 
# 
# 
# 
# 
# 

# #### run sctype HypoMap neurons ####
# burtoni.snseq.combined.sct.neurons = burtoni.snseq.combined.sct
# 
# ## keep cluster ids names
# burtoni.snseq.combined.sct.neurons = AddMetaData(burtoni.snseq.combined.sct.neurons,
#                                                  burtoni.snseq.combined.sct.neurons@meta.data$integrated_snn_res.0.8,
#                                                         col.name = 'integrated_snn_res.0.8.broad')
# 
# #set idents
# Idents(object = burtoni.snseq.combined.sct.neurons) <- "sctypemarkers.hypo"
# 
# #subset to neurons
# burtoni.snseq.combined.sct.neurons = subset(burtoni.snseq.combined.sct.neurons,
#                                             idents = c("C7-1: GLU",
#                                                        "C7-2: GABA",
#                                                        "Unknown"))
# 
# # need to set to integrated for clustering
# DefaultAssay(burtoni.snseq.combined.sct.neurons) = 'integrated'
# 
# #check data loaded correctly
# ## run PCA, UMAP, and cluster 
# #use 0.8 resolution
# burtoni.snseq.combined.sct.neurons = burtoni.snseq.combined.sct.neurons %>% 
#   RunPCA() %>%
#   FindNeighbors(dims = 1:15) %>%
#   RunUMAP(dims = 1:15) %>%
#   FindClusters(resolution = 0.8)
# 
# ## graph 
# # idents to new clusters
# Idents(object = burtoni.snseq.combined.sct.neurons) <- "integrated_snn_res.0.8"
# 
# # graph
# DimPlot(burtoni.snseq.combined.sct.neurons,
#         group.by='sctypemarkers.hypo',
#         label=TRUE) +
#   # umap_theme() +
#   ggtitle('Neurons Hypo')
# ggsave('./sctype.hypo/neurons/UMAP cell types sctype on neurons.png',
#        height = 5.5,
#        width = 6.5)
# 
# ### cell type assignment
# ### prepare gene set list
# # set tissue type
# tissue = "Brain"
# # prepare gene sets
# gs_list = gene_sets_prepare.df(annotation.data.hypomap,
#                                tissue)
# 
# 
# ## test marker specificity
# marker.specificity = marker_sensitivity_score(scRNAseqData = burtoni.snseq.combined.sct.neurons[["SCT"]]@scale.data,
#                                               scaled = TRUE,
#                                               gs = gs_list$gs_positive,
#                                               gs2 = gs_list$gs_negative,
#                                               gene_names_to_uppercase = FALSE)
# 
# ### assign cell types to clusters
# # get cell-type by cell matrix
# es.max = sctype_score(scRNAseqData = burtoni.snseq.combined.sct.neurons[["SCT"]]@scale.data, 
#                       scaled = TRUE, 
#                       gs = gs_list$gs_positive, 
#                       gs2 = gs_list$gs_negative,
#                       gene_names_to_uppercase = FALSE) 
# 
# 
# # es.max = sctype_score(scRNAseqData = burtoni.snseq.combined.sct.neurons[["SCT"]]@data %>% 
# #                             as.matrix(), 
# #                       scaled = TRUE, 
# #                       gs = gs_list$gs_positive, 
# #                       gs2 = gs_list$gs_negative,
# #                       gene_names_to_uppercase = FALSE) 
# 
# # NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
# # In case Seurat is used, it is either pbmc[["RNA"]]@scale.data (default), pbmc[["SCT"]]@scale.data, in case sctransform is used for normalization,
# # or pbmc[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.
# 
# # merge by cluster
# cL_results = do.call("rbind", lapply(unique(burtoni.snseq.combined.sct.neurons@meta.data$seurat_clusters), function(cl){
#   es.max.cl = sort(rowSums(es.max[ ,rownames(burtoni.snseq.combined.sct.neurons@meta.data[burtoni.snseq.combined.sct.neurons@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
#   head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(burtoni.snseq.combined.sct.neurons@meta.data$seurat_clusters==cl)), 10)
# }))
# sctype_scores = cL_results %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
# 
# # set low-confident (low ScType score) clusters to "unknown"
# sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
# print(sctype_scores[,1:3])
# 
# #graph umap
# burtoni.snseq.combined.sct.neurons@meta.data$sctypemarkers.hypo.neurons = ""
# for(j in unique(sctype_scores$cluster)){
#   cl_type = sctype_scores[sctype_scores$cluster==j,]; 
#   burtoni.snseq.combined.sct.neurons@meta.data$sctypemarkers.hypo.neurons[burtoni.snseq.combined.sct.neurons@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
# }
# 
# #graph umap
# DimPlot(burtoni.snseq.combined.sct.neurons,
#         reduction = "umap",
#         label = TRUE,
#         repel = TRUE,
#         group.by = 'sctypemarkers.hypo.neurons') +
#   ggtitle('Neurons Hypo sctype')
# ggsave('./sctype.hypo/neurons/UMAP cell types sctype custom neurons.png',
#        height = 5.5,
#        width = 6.5)
# 
# ## graph counts
# # drop zeros
# burtoni.snseq.combined.sct.neurons@meta.data %>% 
#   dplyr::select(c(Genotype.id,
#                   sctypemarkers.hypo.neurons,
#                   orig.ident)) %>% 
#   table() %>% 
#   as.data.frame() %>% 
#   filter(Freq > 0) %>% 
#   ggplot(aes(x = sctypemarkers.hypo.neurons,
#              y = Freq,
#              color = orig.ident,
#              group = Genotype.id)) + 
#   geom_line() +
#   geom_point() +
#   theme_bw() + 
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# ggsave('./sctype.hypo/neurons/Count cell types sctype neurons.png',
#        height = 5.5,
#        width = 6.5)
# 
# #### graph alluvial sctype broad vs hypomap neurons ###
# ### broad sctype.hypo vs sctype.hypo.neurons
# burtoni.snseq.combined.sct.neurons@meta.data %>% 
#   select(c(sctypemarkers.hypo,
#            sctypemarkers.hypo.neurons)) %>% 
#   table() %>% 
#   as.data.frame() %>% 
#   mutate(Count = sum(Freq)) %>% 
#   ungroup() %>% 
#   mutate(Freq.scale = 100*Freq/Count) %>% 
#   ggplot(aes(axis1 = reorder(sctypemarkers.hypo, -Freq.scale),
#              axis2 = reorder(sctypemarkers.hypo.neurons,-Freq.scale),
#              y = Freq.scale)) +
#   geom_alluvium(aes(fill = sctypemarkers.hypo)) +
#   geom_stratum() +
#   geom_text(stat = "stratum",
#             aes(label = after_stat(stratum))) +
#   scale_x_discrete(limits = c("sctypemarkers.hypo", 
#                               "sctypemarkers.hypo.neurons"),
#                    expand = c(0.15, 0.05)) +
#   scale_fill_viridis_d() +
#   theme_classic() 
# ggsave('./sctype.hypo/neurons/Alluvial cells by sctype.hypo vs sctype.hypo.neurons.png',
#        height = 10,
#        width = 10)
# 
# ### broad sctype.hypo vs sctype.hypo.neurons
# burtoni.snseq.combined.sct.neurons@meta.data %>% 
#   select(c(Cell.type,
#            sctypemarkers.hypo.neurons)) %>% 
#   table() %>% 
#   as.data.frame() %>% 
#   mutate(Count = sum(Freq)) %>% 
#   ungroup() %>% 
#   mutate(Freq.scale = 100*Freq/Count) %>% 
#   ggplot(aes(axis1 = reorder(Cell.type, -Freq.scale),
#              axis2 = reorder(sctypemarkers.hypo.neurons,-Freq.scale),
#              y = Freq.scale)) +
#   geom_alluvium(aes(fill = sctypemarkers.hypo.neurons)) +
#   geom_stratum() +
#   geom_text(stat = "stratum",
#             aes(label = after_stat(stratum))) +
#   scale_x_discrete(limits = c("Cell.type", 
#                               "sctypemarkers.hypo.neurons"),
#                    expand = c(0.15, 0.05)) +
#   scale_fill_viridis_d() +
#   theme_classic() 
# ggsave('./sctype.hypo/neurons/Alluvial cells by Cell.type vs sctype.hypo.neurons.png',
#        height = 10,
#        width = 10)
# 
# 
# 
# 
# 
# 
# ### broad sctype.hypo vs sctype.hypo.neurons
# burtoni.snseq.combined.sct.neurons@meta.data %>% 
#   select(c(sctypemarkers.hypo,
#            sctypemarkers.hypo.neurons)) %>% 
#   table() %>% 
#   as.data.frame() %>% 
#   mutate(Count = sum(Freq)) %>% 
#   ungroup() %>% 
#   mutate(Freq.scale = 100*Freq/Count) %>% 
#   ggplot(aes(axis1 = reorder(sctypemarkers.hypo, -Freq.scale),
#              axis2 = reorder(sctypemarkers.hypo.neurons,-Freq.scale),
#              y = Freq.scale)) +
#   geom_alluvium(aes(fill = sctypemarkers.hypo.neurons)) +
#   geom_stratum() +
#   geom_text(stat = "stratum",
#             aes(label = after_stat(stratum))) +
#   scale_x_discrete(limits = c("sctypemarkers.hypo", 
#                               "sctypemarkers.hypo.neurons"),
#                    expand = c(0.15, 0.05)) +
#   scale_fill_viridis_d() +
#   theme_classic() 
# ggsave('./sctype.hypo/neurons/Alluvial cells by sctype.hypo vs sctype.hypo.neurons.png',
#        height = 10,
#        width = 10)
# 
# 
# 
# 

#### run sctype HypoMap neurons.reduce ####
burtoni.snseq.combined.sct.neurons.reduce = burtoni.snseq.combined.sct.neurons

## keep cluster ids names
burtoni.snseq.combined.sct.neurons.reduce = AddMetaData(burtoni.snseq.combined.sct.neurons.reduce,
                                                        burtoni.snseq.combined.sct.neurons.reduce@meta.data$integrated_snn_res.0.8,
                                                        col.name = 'integrated_snn_res.0.8.neurons')

#set idents
Idents(object = burtoni.snseq.combined.sct.neurons.reduce) <- "sctypemarkers.hypo"

#subset to neurons
burtoni.snseq.combined.sct.neurons.reduce = subset(burtoni.snseq.combined.sct.neurons.reduce,
                                            idents = c("C7-1: GLU",
                                                       "C7-2: GABA"))

# need to set to integrated for clustering
DefaultAssay(burtoni.snseq.combined.sct.neurons.reduce) = 'integrated'

#check data loaded correctly
## run PCA, UMAP, and cluster 
#use 0.8 resolution
burtoni.snseq.combined.sct.neurons.reduce = burtoni.snseq.combined.sct.neurons.reduce %>% 
  RunPCA() %>%
  FindNeighbors(dims = 1:15) %>%
  RunUMAP(dims = 1:15) %>%
  FindClusters(resolution = 0.8)

## graph 
# idents to new clusters
Idents(object = burtoni.snseq.combined.sct.neurons.reduce) <- "integrated_snn_res.0.8"

# graph
DimPlot(burtoni.snseq.combined.sct.neurons.reduce,
        group.by='sctypemarkers.hypo',
        label=TRUE) +
  # umap_theme() +
  ggtitle('Neurons Hypo')
ggsave('./sctype.hypo/neurons/UMAP cell types sctype on neurons reduce.png',
       height = 5.5,
       width = 6.5)

### cell type assignment
### prepare gene set list
# set tissue type
tissue = "Brain"
# prepare gene sets
gs_list = gene_sets_prepare.df(annotation.data.hypomap,
                               tissue)


## test marker specificity
marker.specificity = marker_sensitivity_score(scRNAseqData = burtoni.snseq.combined.sct.neurons.reduce[["SCT"]]@scale.data,
                                              scaled = TRUE,
                                              gs = gs_list$gs_positive,
                                              gs2 = gs_list$gs_negative,
                                              gene_names_to_uppercase = FALSE)

### assign cell types to clusters
# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = burtoni.snseq.combined.sct.neurons.reduce[["SCT"]]@scale.data, 
                      scaled = TRUE, 
                      gs = gs_list$gs_positive, 
                      gs2 = gs_list$gs_negative,
                      gene_names_to_uppercase = FALSE) 


# es.max = sctype_score(scRNAseqData = burtoni.snseq.combined.sct.neurons.reduce[["SCT"]]@data %>% 
#                             as.matrix(), 
#                       scaled = TRUE, 
#                       gs = gs_list$gs_positive, 
#                       gs2 = gs_list$gs_negative,
#                       gene_names_to_uppercase = FALSE) 

# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
# In case Seurat is used, it is either pbmc[["RNA"]]@scale.data (default), pbmc[["SCT"]]@scale.data, in case sctransform is used for normalization,
# or pbmc[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.

# merge by cluster
cL_results = do.call("rbind", lapply(unique(burtoni.snseq.combined.sct.neurons.reduce@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(burtoni.snseq.combined.sct.neurons.reduce@meta.data[burtoni.snseq.combined.sct.neurons.reduce@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(burtoni.snseq.combined.sct.neurons.reduce@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_results %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])

#graph umap
burtoni.snseq.combined.sct.neurons.reduce@meta.data$sctypemarkers.hypo.neurons.reduce = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  burtoni.snseq.combined.sct.neurons.reduce@meta.data$sctypemarkers.hypo.neurons.reduce[burtoni.snseq.combined.sct.neurons.reduce@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

#graph umap
DimPlot(burtoni.snseq.combined.sct.neurons.reduce,
        reduction = "umap",
        label = TRUE,
        repel = TRUE,
        group.by = 'sctypemarkers.hypo.neurons.reduce') +
  ggtitle('Neurons Hypo sctype')
ggsave('./sctype.hypo/neurons/UMAP cell types sctype custom neurons reduce.png',
       height = 5.5,
       width = 6.5)

## graph counts
# drop zeros
burtoni.snseq.combined.sct.neurons.reduce@meta.data %>% 
  dplyr::select(c(Genotype.id,
                  sctypemarkers.hypo.neurons.reduce,
                  orig.ident)) %>% 
  table() %>% 
  as.data.frame() %>% 
  filter(Freq > 0) %>% 
  ggplot(aes(x = sctypemarkers.hypo.neurons.reduce,
             y = Freq,
             color = orig.ident,
             group = Genotype.id)) + 
  geom_line() +
  geom_point() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
ggsave('./sctype.hypo/neurons/Count cell types sctype neurons reduce.png',
       height = 5.5,
       width = 6.5)

#### graph alluvial sctype broad vs hypomap neurons ###
### sctype.hypo.neurons vs sctype.hypo.neurons.reduce
burtoni.snseq.combined.sct.neurons.reduce@meta.data %>% 
  select(c(sctypemarkers.hypo.neurons,
           sctypemarkers.hypo.neurons.reduce)) %>% 
  table() %>% 
  as.data.frame() %>% 
  mutate(Count = sum(Freq)) %>% 
  ungroup() %>% 
  mutate(Freq.scale = 100*Freq/Count) %>% 
  ggplot(aes(axis1 = reorder(sctypemarkers.hypo.neurons, -Freq.scale),
             axis2 = reorder(sctypemarkers.hypo.neurons.reduce,-Freq.scale),
             y = Freq.scale)) +
  geom_alluvium(aes(fill = sctypemarkers.hypo.neurons)) +
  geom_stratum() +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("sctypemarkers.hypo.neurons", 
                              "sctypemarkers.hypo.neurons.reduce"),
                   expand = c(0.15, 0.05)) +
  scale_fill_viridis_d() +
  theme_classic() 
ggsave('./sctype.hypo/neurons/Alluvial cells by sctype.hypo.neurons vs sctype.hypo.neurons.reduce.png',
       height = 10,
       width = 10)


### broad sctype.hypo vs sctype.hypo.neurons.reduce
burtoni.snseq.combined.sct.neurons.reduce@meta.data %>% 
  select(c(sctypemarkers.hypo,
           sctypemarkers.hypo.neurons.reduce)) %>% 
  table() %>% 
  as.data.frame() %>% 
  mutate(Count = sum(Freq)) %>% 
  ungroup() %>% 
  mutate(Freq.scale = 100*Freq/Count) %>% 
  ggplot(aes(axis1 = reorder(sctypemarkers.hypo, -Freq.scale),
             axis2 = reorder(sctypemarkers.hypo.neurons.reduce,-Freq.scale),
             y = Freq.scale)) +
  geom_alluvium(aes(fill = sctypemarkers.hypo.neurons.reduce)) +
  geom_stratum() +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("sctypemarkers.hypo", 
                              "sctypemarkers.hypo.neurons.reduce"),
                   expand = c(0.15, 0.05)) +
  scale_fill_viridis_d() +
  theme_classic() 
ggsave('./sctype.hypo/neurons/Alluvial cells by sctype.hypo vs sctype.hypo.neurons.reduce.png',
       height = 10,
       width = 10)


### neuron clusters vs neuron clusters reduce
burtoni.snseq.combined.sct.neurons.reduce@meta.data %>%
  select(c(sctypemarkers.hypo.neurons,
           sctypemarkers.hypo.neurons.reduce,
           integrated_snn_res.0.8.neurons,
           integrated_snn_res.0.8)) %>%
  table() %>%
  as.data.frame() %>%
  mutate(Count = sum(Freq)) %>%
  ungroup() %>%
  mutate(Freq.scale = 100*Freq/Count) %>%
  ggplot(aes(axis1 = reorder(integrated_snn_res.0.8.neurons, -Freq.scale),
             axis2 = reorder(integrated_snn_res.0.8,-Freq.scale),
             y = Freq.scale)) +
  geom_alluvium(aes(fill = sctypemarkers.hypo.neurons)) +
  geom_stratum() +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("integrated_snn_res.0.8.neurons",
                              "integrated_snn_res.0.8"),
                   expand = c(0.15, 0.05)) +
  scale_fill_viridis_d() +
  theme_classic()
ggsave('./sctype.hypo/neurons/Alluvial cells by neurons cluster vs neurons reduce cluster.png',
       height = 10,
       width = 10)


### broad clusters vs neuron clusters reduce
burtoni.snseq.combined.sct.neurons.reduce@meta.data %>%
  select(c(sctypemarkers.hypo.neurons,
           sctypemarkers.hypo.neurons.reduce,
           integrated_snn_res.0.8.broad,
           integrated_snn_res.0.8)) %>%
  table() %>%
  as.data.frame() %>%
  mutate(Count = sum(Freq)) %>%
  ungroup() %>%
  mutate(Freq.scale = 100*Freq/Count) %>%
  ggplot(aes(axis1 = reorder(integrated_snn_res.0.8.broad, -Freq.scale),
             axis2 = reorder(integrated_snn_res.0.8,-Freq.scale),
             y = Freq.scale)) +
  geom_alluvium(aes(fill = sctypemarkers.hypo.neurons.reduce)) +
  geom_stratum() +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("integrated_snn_res.0.8.broad",
                              "integrated_snn_res.0.8"),
                   expand = c(0.15, 0.05)) +
  scale_fill_viridis_d() +
  theme_classic()
ggsave('./sctype.hypo/neurons/Alluvial cells by broad cluster vs neurons reduce cluster.png',
       height = 10,
       width = 10)

### sctype.hypo.neurons vs sctype.hypo.neurons.reduce
burtoni.snseq.combined.sct.neurons.reduce@meta.data %>%
  select(c(Cell.type,
           sctypemarkers.hypo.neurons.reduce)) %>% 
  table() %>% 
  as.data.frame() %>% 
  mutate(Count = sum(Freq)) %>% 
  ungroup() %>% 
  mutate(Freq.scale = 100*Freq/Count) %>% 
  ggplot(aes(axis1 = reorder(Cell.type, -Freq.scale),
             axis2 = reorder(sctypemarkers.hypo.neurons.reduce,-Freq.scale),
             y = Freq.scale)) +
  geom_alluvium(aes(fill = sctypemarkers.hypo.neurons.reduce)) +
  geom_stratum() +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Cell.type", 
                              "sctypemarkers.hypo.neurons.reduce"),
                   expand = c(0.15, 0.05)) +
  scale_fill_viridis_d() +
  theme_classic() 
ggsave('./sctype.hypo/neurons/Alluvial cells by Cell.type vs sctype.hypo.neurons.reduce.png',
       height = 10,
       width = 10)








#### annotate unknown cells ####
### assign cell type using highest cell type score 
## use graphs in: './sctype.hypo/clusters/'
## use graph: './sctype.hypo/Cluster cell type tree.png'

### new annotations:
#cluster,cell.type
#7,GLU
#3,GABA
#10,GABA

### pull cluster and sctype.hypo calls
data.cluster.sctype.hypo = burtoni.snseq.combined.sct@meta.data %>% 
  dplyr::select(seurat_clusters,
                sctypemarkers.hypo,
                sctypemarkers.hypo.broad) %>% 
  distinct() 
# remove rownames
row.names(data.cluster.sctype.hypo) <- NULL

### replace unknowns based on new annotations
data.cluster.sctype.hypo = data.cluster.sctype.hypo %>% 
  mutate(sctypemarkers.hypo.broad = case_when(seurat_clusters == 7 ~ " GLU",
                                              seurat_clusters == 3 ~ " GABA",
                                              seurat_clusters == 10 ~ " GABA",
                                              TRUE ~ sctypemarkers.hypo.broad)) %>% 
  mutate(sctypemarkers.hypo = case_when(seurat_clusters == 7 ~ "C7-1: GLU",
                                        seurat_clusters == 3 ~ "C7-2: GABA",
                                        seurat_clusters == 10 ~ "C7-2: GABA",
                                        TRUE ~ sctypemarkers.hypo)) 


## combine with cell.id
# extract metadata to save
sctypemarkers.hypo.unknown = burtoni.snseq.combined.sct@meta.data %>%
  dplyr::select(seurat_clusters) %>%
  rownames_to_column('Cell.id') %>% 
  full_join(data.cluster.sctype.hypo)

#### save metadata unknown ####
#save file
write.csv(sctypemarkers.hypo.unknown,
          file = './sctype.hypo/sctypemarkers.hypo.unknown.csv')

### load sctype.hypo data
sctypemarkers.hypo.unknown = read.csv('./sctype.hypo/sctypemarkers.hypo.unknown.csv')

## sctype.hypo
burtoni.snseq.combined.sct = AddMetaData(
  object = burtoni.snseq.combined.sct,
  metadata = sctypemarkers.hypo.unknown %>% 
    select(sctypemarkers.hypo,
           Cell.id) %>% 
    column_to_rownames(var = "Cell.id"),
  col.name = 'sctypemarkers.hypo'
)

## sctype.hypo
burtoni.snseq.combined.sct = AddMetaData(
  object = burtoni.snseq.combined.sct,
  metadata = sctypemarkers.hypo.unknown %>% 
    select(sctypemarkers.hypo.broad,
           Cell.id) %>% 
    column_to_rownames(var = "Cell.id"),
  col.name = 'sctypemarkers.hypo.broad'
)


#### cell counts and umap presentation unknown ####
### cell counts
burtoni.snseq.combined.sct@meta.data %>% 
  dplyr::select(sctypemarkers.hypo,
                Genotype.id,
                orig.ident) %>% 
  table() %>% 
  as.data.frame() %>% 
  filter(Freq > 0) %>% 
  ggplot(aes(x = sctypemarkers.hypo,
             y = Freq,
             fill = orig.ident)) +
  geom_point(size = 8,
             shape = 21) +
  scale_fill_manual(values = c("#60bb46", "#4e499e")) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        text = element_text(size = 20),
        legend.position = 'none') +
  ylab("Total cells assigned") +
  xlab('') +
  labs('')
ggsave('sctype.hypo/unknown/Count cell per status.png',
       height = 6,
       width = 6)

# broad
burtoni.snseq.combined.sct@meta.data %>% 
  dplyr::select(sctypemarkers.hypo.broad,
                Genotype.id,
                orig.ident) %>% 
  table() %>% 
  as.data.frame() %>% 
  filter(Freq > 0) %>% 
  mutate(order = Freq,
         order = ifelse(sctypemarkers.hypo.broad == "Unknown",
                        0,
                        order)) %>% 
  ggplot(aes(x = reorder(sctypemarkers.hypo.broad,
                         -order),
             y = Freq,
             fill = orig.ident)) +
  geom_point(size = 8,
             shape = 21) +
  scale_fill_manual(values = c("#60bb46", "#4e499e")) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        text = element_text(size = 20),
        legend.position = 'none') +
  ylab("Total cells assigned") +
  xlab('') +
  labs('')
ggsave('sctype.hypo/unknown/Count cell broad per status.png',
       height = 6,
       width = 6)

# ## UMAP for presentation
# # cell type broad
# DimPlot(burtoni.snseq.combined.sct,
#         group.by = 'sctypemarkers.hypo.broad',
#         label = T) +
#   NoLegend() + 
#   theme(
#     plot.title = element_blank(),
#     axis.title.x = element_blank(),
#     axis.title.y = element_blank())
# ggsave('sctype.hypo/unknown/UMAP cell broad presentation.png',
#        height = 5.5,
#        width = 5.5)
# 
# # grey umap
# DimPlot(burtoni.snseq.combined.sct,
#         group.by = 'orig.ident',
#         cols = c('grey',
#                  'grey')) +
#   NoLegend()+ 
#   theme(
#     plot.title = element_blank(),
#     axis.title.x = element_blank(),
#     axis.title.y = element_blank())
# ggsave('sctype.hypo/unknown/UMAP grey presentation.png',
#        height = 5.5,
#        width = 5.5)
# 
# # clusters umap
# DimPlot(burtoni.snseq.combined.sct,
#         group.by = 'seurat_clusters',
#         label = T) +
#   NoLegend() + 
#   theme(
#     plot.title = element_blank(),
#     axis.title.x = element_blank(),
#     axis.title.y = element_blank())
# ggsave('sctype.hypo/unknown/UMAP clusters presentation.png',
#        height = 5.5,
#        width = 5.5)

## UMAP for paper
# just cell labels
DimPlot(burtoni.snseq.combined.sct,
        group.by = 'sctypemarkers.hypo.broad',
        alpha = 0,
        label = T,
        label.box = T,
        cols = c('#F8766D',
                 '#B79F00',
                 '#00BA38',
                 '#00BFC4',
                 '#619CFF',
                 '#F564E3',
                 'grey')) +
  NoLegend() +
  theme(
    plot.title = element_blank()
  )
ggsave('sctype.hypo/UMAP cell broad paper labels.png',
       height = 7,
       width = 7)

# just cluster labels
DimPlot(burtoni.snseq.combined.sct,
        group.by = 'seurat_clusters',
        alpha = 0,
        label = T) +
  NoLegend() +
  theme(
    plot.title = element_blank()
  )
ggsave('sctype.hypo/UMAP cell broad paper clusters.png',
       height = 7,
       width = 7)

# just color cells
DimPlot(burtoni.snseq.combined.sct,
        group.by = 'sctypemarkers.hypo.broad',
        label = F,
        cols = c('#F8766D',
                 '#B79F00',
                 '#00BA38',
                 '#00BFC4',
                 '#619CFF',
                 '#F564E3',
                 'grey')) +
  NoLegend() +
  theme(
    plot.title = element_blank())
ggsave('sctype.hypo/UMAP cell broad paper.png',
       height = 7,
       width = 7)

#### cell type percentages unknown ####
### percentages by cluster
## get counts of nuclei per cluster per genotype
# create table of samples by cluster
table_samples_by_clusters <- burtoni.snseq.combined.sct@meta.data %>% 
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
  arrange(factor(Genotype.id, levels = levels(burtoni.snseq.combined.sct@meta.data$Genotype.id)))

# create table of clusters by sample
table_clusters_by_samples <- burtoni.snseq.combined.sct@meta.data %>%
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
  arrange(factor(seurat_clusters, levels = levels(burtoni.snseq.combined.sct@meta.data$seurat_clusters)))

## graph percentages by cluster
# get count per sample
temp_labels <- burtoni.snseq.combined.sct@meta.data %>%
  group_by(Genotype.id) %>%
  tally()

# create graph of sample by cluster
p1 <- table_samples_by_clusters %>%
  dplyr::select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'Genotype.id') %>%
  # mutate(Genotype.id = factor(Genotype.id, levels = levels(burtoni.snseq.combined.sct@meta.data$Genotype.id))) %>%
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
temp_labels <- burtoni.snseq.combined.sct@meta.data %>%
  group_by(seurat_clusters) %>%
  tally()

# create graph of cluster by sample
p2 <- table_clusters_by_samples %>%
  select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'seurat_clusters') %>%
  mutate(seurat_clusters = factor(seurat_clusters, levels = levels(burtoni.snseq.combined.sct@meta.data$seurat_clusters))) %>%
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
  'sctype.hypo/unknown/Percentage genotype by clusters.png',
  p1 + p2 +
    patchwork::plot_layout(ncol = 2, widths = c(
      burtoni.snseq.combined.sct@meta.data$Genotype.id %>% unique() %>% length(),
      burtoni.snseq.combined.sct@meta.data$seurat_clusters %>% unique() %>% length()
    )),
  width = 18, height = 8
)

## just percentage per cluster 
# add bias lines
ggsave(
  'sctype.hypo/unknown/Percentage genotype by clusters bias.png',
  p2 +
    geom_hline(yintercept = .67,
               linetype = 'dashed') +
    geom_hline(yintercept = .33,
               linetype = 'dashed'),
  width = 18, 
  height = 9
)



### percentages by celltype
## get counts of nuclei per celltype per genotype
# create table of samples by celltype
table_samples_by_celltype <- burtoni.snseq.combined.sct@meta.data %>% 
  group_by(Genotype.id, 
           sctypemarkers.hypo.broad) %>%
  summarize(count = n()) %>%
  spread(sctypemarkers.hypo.broad, 
         count, 
         fill = 0) %>%
  ungroup() %>%
  mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  dplyr::select(c('Genotype.id', 
                  'total_cell_count', 
                  everything())) %>%
  arrange(factor(Genotype.id, levels = levels(burtoni.snseq.combined.sct@meta.data$Genotype.id)))

# create table of celltype by sample
table_celltype_by_samples <- burtoni.snseq.combined.sct@meta.data %>%
  group_by(Genotype.id, 
           sctypemarkers.hypo.broad) %>%
  summarize(count = n()) %>%
  spread(Genotype.id,
         count, 
         fill = 0) %>%
  ungroup() %>%
  mutate(total_cell_count = rowSums(.[c(2:ncol(.))])) %>%
  select(c('sctypemarkers.hypo.broad', 
           'total_cell_count', 
           everything())) %>%
  arrange(factor(sctypemarkers.hypo.broad, levels = levels(burtoni.snseq.combined.sct@meta.data$sctypemarkers.hypo.broad)))

## graph percentages by celltype
# get count per sample
temp_labels <- burtoni.snseq.combined.sct@meta.data %>%
  group_by(Genotype.id) %>%
  tally()

# create graph of sample by celltype
p1 <- table_samples_by_celltype %>%
  dplyr::select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'Genotype.id') %>%
  # mutate(Genotype.id = factor(Genotype.id, levels = levels(burtoni.snseq.combined.sct@meta.data$Genotype.id))) %>%
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

# get count per celltype
temp_labels <- burtoni.snseq.combined.sct@meta.data %>%
  group_by(sctypemarkers.hypo.broad) %>%
  tally()

# create graph of celltype by sample
p2 <- table_celltype_by_samples %>%
  select(-c('total_cell_count')) %>%
  reshape2::melt(id.vars = 'sctypemarkers.hypo.broad') %>% 
  # mutate(sctypemarkers.hypo.broad = factor(sctypemarkers.hypo.broad, levels = levels(burtoni.snseq.combined.sct@meta.data$sctypemarkers.hypo.broad))) %>%
  ggplot(aes(sctypemarkers.hypo.broad, value)) +
  geom_bar(aes(fill = variable), position = 'fill', stat = 'identity', color = 'black') +
  geom_text(
    data = temp_labels, aes(x = sctypemarkers.hypo.broad, y = Inf, label = paste0(format(n, big.mark = ',', trim = TRUE)), vjust = -1),
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
    plot.margin = margin(t = 20, r = 0, b = 0, l = 10, unit = 'pt'),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )

ggsave(
  'sctype.hypo/unknown/Percentage genotype by celltype.png',
  p1 + p2 +
    patchwork::plot_layout(ncol = 2, widths = c(
      burtoni.snseq.combined.sct@meta.data$Genotype.id %>% unique() %>% length(),
      burtoni.snseq.combined.sct@meta.data$sctypemarkers.hypo.broad %>% unique() %>% length()
    )),
  width = 18, height = 8
)

#### plot cell types tree unknown ####
# https://romanhaa.github.io/projects/scrnaseq_workflow/#clustering

### create tree of cell types
burtoni.snseq.combined.sct.tree <- BuildClusterTree(
  burtoni.snseq.combined.sct,
  dims = 1:15,
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
                                 sctypemarkers.hypo.broad) %>% 
                   distinct(),
                 by = c("label" = "seurat_clusters"))

# convert back to tree
tree = treeio::as.treedata(tree)

### graph tree of celltypes
ggtree::ggtree(tree, aes(x, y)) +
  scale_y_reverse() +
  ggtree::geom_tree() +
  ggtree::theme_tree() +
  ggtree::geom_tiplab(offset = 1) +
  ggtree::geom_tippoint(aes(color = sctypemarkers.hypo.broad),
                        shape = 16,
                        size = 5) +
  coord_cartesian(clip = 'off') +
  theme(plot.margin = unit(c(0,2.5,0,0),
                           'cm'))
ggsave('sctype.hypo/unknown/Cluster cell type tree.png',
       height = 6,
       width = 6)

