#### Burtoni snseq seurat analysis
### scsorter
###Note: Seurat requires R version > 4
## use lambcomp1 to run R with command 
# > R-4.0.3
##https://cran.r-project.org/web/packages/scSorter/vignettes/scSorter.html


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
#load libraries
library(tidyverse)
library(Seurat)
library(patchwork)
library(clustree)
library(scSorter)


#### load data ####
### single cell data
# load('burtoni.snseq.combined.sct.all.RData')
load('burtoni.snseq.combined.sct.RData')

burtoni.scsorter.data = burtoni.snseq.combined.sct

## meta data
burtoni.scsorter.metadata = burtoni.snseq.combined.sct@meta.data

#extract umap data
burtoni.snseq.combined.sct.umap = burtoni.snseq.combined.sct@reductions$umap@cell.embeddings %>% 
  as.data.frame() %>% 
  rownames_to_column("Cell.id")
  
### annotation data
## subset data from 'Cell.type.gene.markers.csv'
## load annotation data
annotation.data = read.csv("../Gene.lists/Cell.type.annotation.genes.csv")
## remove NA
annotation.data = annotation.data %>% 
  na.omit() %>% 
  distinct()

# #remove neuron
# annotation.data = annotation.data %>% 
#   filter(Type != 'neurons') 
  

##rename duplicate rownames
## can't handle upper/lowercase of same genes
## ex. ‘SLC32A1’, ‘SLCO1C1’
annotation.data = annotation.data %>% 
  filter(Marker != "SLCO1C1") %>% 
  filter(Marker != "SLC32A1")

### load HypoMap annotation
annotation.data.hypomap = read.csv("../Gene.lists/Cell.type.gene.markers.HypoMap.annotation.ensemble.csv")

## rename columns
# subset data to relevant columns
# remove duplicates
annotation.data.hypomap = annotation.data.hypomap %>% 
  dplyr::rename('Marker' = 'tilapia_gene_name') %>% 
  dplyr::rename('Type' = 'cluster_name') %>% 
  dplyr::select(c(Marker,
                  Type)) %>% 
  distinct() %>% 
  na.omit()

## check for multiple genes
annotation.data.hypomap %>% 
  distinct() %>% 
  group_by(Marker) %>% 
  filter(n()>1) %>% 
  nrow()

#### preprocessing data ####
### subset highly variable genes
# identify top 2000 variable genes
topgenes <- head(VariableFeatures(burtoni.scsorter.data,
                                  assay = 'integrated'), 
                 2000)

#get assay data
#as matrix
burtoni.scsorter.data.extract = GetAssayData(burtoni.scsorter.data,
                                             assay = 'SCT',
                                             slot = 'data') %>% 
  as.matrix()

# filter genes
# filter out genes with non-zero expression in less than 10% of total cells
topgene_filter = rowSums(as.matrix(burtoni.scsorter.data.extract)[topgenes, ]!=0) > ncol(burtoni.scsorter.data.extract)*.1
#create filter gene list
topgenes = topgenes[topgene_filter]

# create list of genes to subset
picked_genes = unique(c(annotation.data$Marker, 
                        topgenes))

# subset data by gene list
burtoni.scsorter.data.extract = burtoni.scsorter.data.extract[rownames(burtoni.scsorter.data.extract) %in% picked_genes, ]

## find annotation genes not listed
non.overlap.genes = annotation.data %>% 
  pull(Marker) %>% 
  setdiff(intersect(annotation.data %>% 
                      pull(Marker), 
                    burtoni.scsorter.data.extract %>% 
                      rownames()))
#create new annotation dataframe
annotation.data.overlap = annotation.data %>% 
  filter(!Marker %in% non.overlap.genes)

### HypoMap
# create list of genes to subset
picked_genes.hypo = unique(c(annotation.data.hypomap$Marker, 
                        topgenes))

# subset data by gene list
burtoni.scsorter.data.extract.hypo = burtoni.scsorter.data.extract[rownames(burtoni.scsorter.data.extract) %in% picked_genes.hypo, ]

## find annotation genes not listed
non.overlap.genes.hypo = annotation.data.hypomap %>% 
  pull(Marker) %>% 
  setdiff(intersect(annotation.data.hypomap %>% 
                      pull(Marker), 
                    burtoni.scsorter.data.extract.hypo %>% 
                      rownames()))

#create new annotation dataframe
annotation.data.hypomap.overlap = annotation.data.hypomap %>% 
  filter(!Marker %in% non.overlap.genes.hypo)

## check for duplicated gene
duplicated.genes.hypo = burtoni.scsorter.data.extract.hypo %>% 
  as.data.frame() %>% 
  rownames_to_column('Gene') %>% 
  dplyr::select(Gene) %>% 
  mutate(Gene.UP = toupper(Gene)) %>% 
  group_by(Gene.UP) %>% 
  filter(n()>1) %>%
  mutate(keep = ifelse(Gene == Gene.UP,
                       1,
                       0)) %>% 
  filter(keep == 1) %>% 
  mutate(Gene.rename = paste(Gene,
                             '.UP',
                             sep = '')) %>% 
  as.data.frame()

## rename upper case genes in dataframe
burtoni.scsorter.data.extract.hypo.rownames = burtoni.scsorter.data.extract.hypo %>% 
  rownames() %>% 
  as.data.frame() %>% 
  dplyr::rename('Gene' = '.') %>% 
  full_join(duplicated.genes.hypo) %>% 
  mutate(Gene = ifelse(!is.na(keep),
                       Gene.rename,
                       Gene))
# rename genes
rownames(burtoni.scsorter.data.extract.hypo) = burtoni.scsorter.data.extract.hypo.rownames$Gene




#### run scsorter ####
### run scsorter command
### scsorter
burtoni.scsorter.data.scsort <- scSorter(burtoni.scsorter.data.extract, 
                                         annotation.data.overlap)

#results
print(table(burtoni.scsorter.data.scsort$Pred_Type))


#### graph results ####
# graph counts
table(burtoni.scsorter.data.scsort$Pred_Type) %>% 
  as.data.frame() %>%  
  ggplot(aes(x = reorder(Var1,
                         -Freq),
             y = Freq)) +
  geom_point() +
  theme_bw() +
  xlab("Cell.type") +
  ylab('Count') +
  ggtitle("scsorter cell type count")
ggsave('scsorter/scsorter cell type count.png',
       width = 10,
       height = 10)

## remove neurons
table(burtoni.scsorter.data.scsort$Pred_Type) %>% 
  as.data.frame() %>% 
  filter(Var1 != 'neurons') %>% 
  ggplot(aes(x = reorder(Var1,
                         -Freq),
             y = Freq)) +
  geom_point() +
  theme_bw() +
  xlab("Cell.type") +
  ylab('Count') +
  ggtitle("scsorter cell type count")
ggsave('scsorter/scsorter cell type count neuron remove.png',
       width = 10,
       height = 10)

## add cell id to output
burtoni.scsorter.data.scsort.output = data.frame(Cell.type = burtoni.scsorter.data.scsort$Pred_Type,
  Cell.id = burtoni.scsorter.data.extract %>% 
           colnames())

##combine with metadata
burtoni.scsorter.data.scsort.output = left_join(burtoni.scsorter.data.scsort.output,
                                                burtoni.scsorter.metadata %>% 
                                                  rownames_to_column("Cell.id"))

# graph counts
burtoni.scsorter.data.scsort.output %>%
  select(Cell.type,
         orig.ident) %>%
  table() %>%
  as.data.frame() %>% 
  ggplot(aes(x = reorder(Cell.type,
                         -Freq),
             y = Freq,
             color = orig.ident)) +
  geom_point(size = 5) +
  theme_bw() +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("Cell.type") +
  ylab('Count') +
  ggtitle("scsorter cell type count")
ggsave('scsorter/scsorter cell type count status.png',
       width = 10,
       height = 10)

# remove neurons
burtoni.scsorter.data.scsort.output %>%
  select(Cell.type,
         orig.ident) %>%
  table() %>%
  as.data.frame() %>% 
  filter(Cell.type != 'neurons') %>% 
  ggplot(aes(x = reorder(Cell.type,
                         -Freq),
             y = Freq,
             color = orig.ident)) +
  geom_point(size = 5) +
  theme_bw() +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("Cell.type") +
  ylab('Count') +
  ggtitle("scsorter cell type count")
ggsave('scsorter/scsorter cell type count status neuron removed.png',
       width = 10,
       height = 10)

### compare to module score
### create violin plot
## compare for each module
# create list of modules
module.scores = burtoni.scsorter.data.scsort.output %>% 
  colnames() %>% 
  str_subset(".score1")

# create graph for each module
for (i in module.scores) {
  burtoni.scsorter.data.scsort.output %>% 
    ggplot(aes_string(x = paste0("reorder(Cell.type, -",
                                 i,
                                 ")"),
                      y = i)) +
    geom_hline(yintercept = 0) +
    geom_violin(trim = FALSE) +
    theme_bw()+ 
    xlab('') +
    stat_summary(fun.data=mean_sdl,
                 geom="pointrange",
                 color="red") +
    ggtitle(paste0('Scsorter vs module ',
                   i))
  #save
ggsave(paste('scsorter/module_score/',
             i,
             ' vs scsorter.png',
             sep = ''),
       width = 10,
       height = 10)
}

##combine umap with metadata
burtoni.scsorter.data.scsort.output.umap =  full_join(burtoni.scsorter.data.scsort.output,
                                                      burtoni.snseq.combined.sct.umap)
## broad cell class
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


# for poster
#graph umap with cell markers
#broad cell class
# seperate doms and subs
burtoni.scsorter.data.scsort.output.umap %>% 
  mutate(orig.ident = ifelse(orig.ident == 'dom_burtoni_snseq',
                             'Dominant',
                             'Subordinate')) %>% 
  ggplot(aes(x=UMAP_1,
             y= UMAP_2,
             color = Cell.class.broad)) +
  geom_point() +
  theme_classic() +
  scale_color_manual(values = c('#be5927',
                                '#4e499e',
                                '#60bb46',
                                'grey'),
                     breaks = c('neurons',
                                'glia',
                                'vascular',
                                'unknown')) +
  facet_grid( ~ orig.ident) +
  labs(color='Broad cell classes') +
  theme(legend.text=element_text(size=15),
        legend.title=element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=5)))+ theme(strip.background = element_blank()) +
  theme(strip.text.x = element_text(size = 15))
ggsave('scsorter/UMAP scsorter cell cell types broad sub vs dom poster.pdf',
       width = 6,
       height = 5,
       units = "in",
       dpi = 300)

## dotplot
burtoni.scsorter.data.scsort.output.umap %>% 
  group_by(seurat_clusters,
           Cell.type) %>% 
  summarise(count = n()) %>% 
  group_by(seurat_clusters) %>% 
  mutate(total = sum(count),
         percent = 100*count/total) %>% 
  ggplot(aes(y = seurat_clusters,
             x = Cell.type,
             size = percent,
             color = percent)) +
  geom_point() +
  theme_bw()
ggsave('scsorter/dotplot scsorter cell types.png',
       width = 10,
       height = 10)

#broad cell class
burtoni.scsorter.data.scsort.output.umap %>% 
  group_by(seurat_clusters,
           Cell.class.broad) %>% 
  summarise(count = n()) %>% 
  group_by(seurat_clusters) %>% 
  mutate(total = sum(count),
         percent = 100*count/total) %>% 
  ggplot(aes(y = seurat_clusters,
             x = Cell.class.broad,
             size = percent,
             color = percent)) +
  geom_point() +
  theme_bw()
ggsave('scsorter/dotplot scsorter cell types broad.png',
       width = 10,
       height = 10)




#### compare cell types ####
#### cluster cell classes using average expression

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


### calculate average expression 
## all cell types
burtoni.snseq.combined.sct.avg.celltype = AverageExpression(burtoni.snseq.combined.sct, 
                  verbose = FALSE,
                  assays = 'SCT',
                  group.by = 'Cell.type',
                  slot = 'data')

#create datframe
#filter out all but picked genes
burtoni.snseq.combined.sct.avg.celltype.df = as_tibble(burtoni.snseq.combined.sct.avg.celltype[["SCT"]],
                                                           rownames = NA) %>% 
  rownames_to_column('gene') %>% 
  filter(gene %in% picked_genes)

### using all cell types
#create dendrogram
#remove rows with all zero
burtoni.snseq.combined.sct.avg.celltype.tree = pvclust(burtoni.snseq.combined.sct.avg.celltype.df %>% 
                                                                  select(-c(gene)) %>%                                                            mutate(sum = rowSums(.)) %>% 
                                                                  filter(sum > 0) %>% 
                                                                  select(-c(sum)) %>% 
                                                             select(-c('Unknown')))

#graph
png(file = 'scsorter/pvclust scsorter cell types.png',
    width = 12,
    height = 6,
    units = 'in',
    res = 1080)
plot(burtoni.snseq.combined.sct.avg.celltype.tree,
     cex.pv	= 1,
     cex = 2)
dev.off()



#### scsorter full ####
### create new function
scsorter.full <- function (expr, anno, default_weight = 2, n_start = 10, alpha = 0, 
                           u = 0.05, max_iter = 100, setseed = 0) 
{
  anno_processed = design_matrix_builder(anno, default_weight)
  dt = data_preprocess(expr, anno_processed)
  dat = dt[[1]]
  designmat = dt[[2]]
  weightmat = dt[[3]]
  c_cost = NULL
  c_mu = list()
  c_clus = list()
  for (i in 1:n_start) {
    set.seed(i + setseed)
    t1 = Sys.time()
    pred_ot = update_func(as.matrix(dat), designmat, weightmat, 
                          unknown_threshold1 = alpha, unknown_threshold2 = u, 
                          max_iter = max_iter)
    t2 = Sys.time()
    c_cost = c(c_cost, pred_ot[[3]])
    c_mu[[i]] = pred_ot[[1]]
    c_clus[[i]] = pred_ot[[2]]
  }
  pk = which.min(c_cost)
  pred_clus = c_clus[[pk]]
  pred_clus = c(colnames(designmat), rep("Unknown", ncol(designmat)))[pred_clus]
  pred_mu = c_mu[[pk]]
  return(list(Pred_Type = pred_clus, Pred_param = pred_mu, Full_type = c_cost))
}

#set namespace for new function
environment(scsorter.full) <- asNamespace('scSorter')

## run scsort full
#~20 min
burtoni.scsorter.data.scsort.full.cost <- scsorter.full(burtoni.scsorter.data.extract, 
                                         annotation.data.overlap)

#results
print(table(burtoni.scsorter.data.scsort.full.cost$Pred_Type))





































#### save data ####
save(burtoni.scsorter.data.scsort.output,
     file = "burtoni.scsorter.data.scsort.output.new.RData")
load("burtoni.scsorter.data.scsort.output.RData")
























#### run scsorter hypomap ####
### run scsorter command
### scsorter
burtoni.scsorter.data.scsort.hypo <- scSorter(burtoni.scsorter.data.extract.hypo, 
                                         annotation.data.hypomap.overlap)





#results
print(table(burtoni.scsorter.data.scsort.hypo$Pred_Type))


#### graph results hypomap ####
# graph counts
table(burtoni.scsorter.data.scsort$Pred_Type) %>% 
  as.data.frame() %>%  
  ggplot(aes(x = reorder(Var1,
                         -Freq),
             y = Freq)) +
  geom_point() +
  theme_bw() +
  xlab("Cell.type") +
  ylab('Count') +
  ggtitle("scsorter cell type count")
ggsave('scsorter/scsorter cell type count.png',
       width = 10,
       height = 10)

## remove neurons
table(burtoni.scsorter.data.scsort$Pred_Type) %>% 
  as.data.frame() %>% 
  filter(Var1 != 'neurons') %>% 
  ggplot(aes(x = reorder(Var1,
                         -Freq),
             y = Freq)) +
  geom_point() +
  theme_bw() +
  xlab("Cell.type") +
  ylab('Count') +
  ggtitle("scsorter cell type count")
ggsave('scsorter/scsorter cell type count neuron remove.png',
       width = 10,
       height = 10)

## add cell id to output
burtoni.scsorter.data.scsort.output = data.frame(Cell.type = burtoni.scsorter.data.scsort$Pred_Type,
                                                 Cell.id = burtoni.scsorter.data.extract %>% 
                                                   colnames())

##combine with metadata
burtoni.scsorter.data.scsort.output = left_join(burtoni.scsorter.data.scsort.output,
                                                burtoni.scsorter.metadata %>% 
                                                  rownames_to_column("Cell.id"))

# graph counts
burtoni.scsorter.data.scsort.output %>%
  select(Cell.type,
         orig.ident) %>%
  table() %>%
  as.data.frame() %>% 
  ggplot(aes(x = reorder(Cell.type,
                         -Freq),
             y = Freq,
             color = orig.ident)) +
  geom_point(size = 5) +
  theme_bw() +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("Cell.type") +
  ylab('Count') +
  ggtitle("scsorter cell type count")
ggsave('scsorter/scsorter cell type count status.png',
       width = 10,
       height = 10)

# remove neurons
burtoni.scsorter.data.scsort.output %>%
  select(Cell.type,
         orig.ident) %>%
  table() %>%
  as.data.frame() %>% 
  filter(Cell.type != 'neurons') %>% 
  ggplot(aes(x = reorder(Cell.type,
                         -Freq),
             y = Freq,
             color = orig.ident)) +
  geom_point(size = 5) +
  theme_bw() +
  theme(text = element_text(size = 20),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("Cell.type") +
  ylab('Count') +
  ggtitle("scsorter cell type count")
ggsave('scsorter/scsorter cell type count status neuron removed.png',
       width = 10,
       height = 10)

### compare to module score
### create violin plot
## compare for each module
# create list of modules
module.scores = burtoni.scsorter.data.scsort.output %>% 
  colnames() %>% 
  str_subset(".score1")

# create graph for each module
for (i in module.scores) {
  burtoni.scsorter.data.scsort.output %>% 
    ggplot(aes_string(x = paste0("reorder(Cell.type, -",
                                 i,
                                 ")"),
                      y = i)) +
    geom_hline(yintercept = 0) +
    geom_violin(trim = FALSE) +
    theme_bw()+ 
    xlab('') +
    stat_summary(fun.data=mean_sdl,
                 geom="pointrange",
                 color="red") +
    ggtitle(paste0('Scsorter vs module ',
                   i))
  #save
  ggsave(paste('scsorter/module_score/',
               i,
               ' vs scsorter.png',
               sep = ''),
         width = 10,
         height = 10)
}

##combine umap with metadata
burtoni.scsorter.data.scsort.output.umap =  full_join(burtoni.scsorter.data.scsort.output,
                                                      burtoni.snseq.combined.sct.umap)
## broad cell class
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


# for poster
#graph umap with cell markers
#broad cell class
# seperate doms and subs
burtoni.scsorter.data.scsort.output.umap %>% 
  mutate(orig.ident = ifelse(orig.ident == 'dom_burtoni_snseq',
                             'Dominant',
                             'Subordinate')) %>% 
  ggplot(aes(x=UMAP_1,
             y= UMAP_2,
             color = Cell.class.broad)) +
  geom_point() +
  theme_classic() +
  scale_color_manual(values = c('#be5927',
                                '#4e499e',
                                '#60bb46',
                                'grey'),
                     breaks = c('neurons',
                                'glia',
                                'vascular',
                                'unknown')) +
  facet_grid( ~ orig.ident) +
  labs(color='Broad cell classes') +
  theme(legend.text=element_text(size=15),
        legend.title=element_text(size=15))+
  guides(colour = guide_legend(override.aes = list(size=5)))+ theme(strip.background = element_blank()) +
  theme(strip.text.x = element_text(size = 15))
ggsave('scsorter/UMAP scsorter cell cell types broad sub vs dom poster.pdf',
       width = 6,
       height = 5,
       units = "in",
       dpi = 300)

## dotplot
burtoni.scsorter.data.scsort.output.umap %>% 
  group_by(seurat_clusters,
           Cell.type) %>% 
  summarise(count = n()) %>% 
  group_by(seurat_clusters) %>% 
  mutate(total = sum(count),
         percent = 100*count/total) %>% 
  ggplot(aes(y = seurat_clusters,
             x = Cell.type,
             size = percent,
             color = percent)) +
  geom_point() +
  theme_bw()
ggsave('scsorter/dotplot scsorter cell types.png',
       width = 10,
       height = 10)

#broad cell class
burtoni.scsorter.data.scsort.output.umap %>% 
  group_by(seurat_clusters,
           Cell.class.broad) %>% 
  summarise(count = n()) %>% 
  group_by(seurat_clusters) %>% 
  mutate(total = sum(count),
         percent = 100*count/total) %>% 
  ggplot(aes(y = seurat_clusters,
             x = Cell.class.broad,
             size = percent,
             color = percent)) +
  geom_point() +
  theme_bw()
ggsave('scsorter/dotplot scsorter cell types broad.png',
       width = 10,
       height = 10)




#### compare cell types hypomap ####
#### cluster cell classes using average expression

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


### calculate average expression 
## all cell types
burtoni.snseq.combined.sct.avg.celltype = AverageExpression(burtoni.snseq.combined.sct, 
                                                            verbose = FALSE,
                                                            assays = 'SCT',
                                                            group.by = 'Cell.type',
                                                            slot = 'data')

#create datframe
#filter out all but picked genes
burtoni.snseq.combined.sct.avg.celltype.df = as_tibble(burtoni.snseq.combined.sct.avg.celltype[["SCT"]],
                                                       rownames = NA) %>% 
  rownames_to_column('gene') %>% 
  filter(gene %in% picked_genes)

### using all cell types
#create dendrogram
#remove rows with all zero
burtoni.snseq.combined.sct.avg.celltype.tree = pvclust(burtoni.snseq.combined.sct.avg.celltype.df %>% 
                                                         select(-c(gene)) %>%                                                            mutate(sum = rowSums(.)) %>% 
                                                         filter(sum > 0) %>% 
                                                         select(-c(sum)) %>% 
                                                         select(-c('Unknown')))

#graph
png(file = 'scsorter/pvclust scsorter cell types.png',
    width = 12,
    height = 6,
    units = 'in',
    res = 1080)
plot(burtoni.snseq.combined.sct.avg.celltype.tree,
     cex.pv	= 1,
     cex = 2)
dev.off()



#### save data ####
save(burtoni.scsorter.data.scsort.output,
     file = "burtoni.scsorter.data.scsort.output.new.RData")
load("burtoni.scsorter.data.scsort.output.RData")























