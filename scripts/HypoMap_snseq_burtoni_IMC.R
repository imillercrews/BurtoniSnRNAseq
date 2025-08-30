#### Burtoni snseq seurat analysis
### cell type assignment with HypoMap
###Note: Seurat requires R version > 4
## use lambcomp1 to run R with command 
# > R-4.0.3
###hypomap
##https://github.com/lsteuernagel/hypoMap_paper
###mapscvi
# https://github.com/lsteuernagel/mapscvi


### set working directory
setwd("/stor/work/Hofmann/All_projects/A_burtoni_snseq/seurat/")

#### load libraries ####
## install seurat disk
# library(remotes)
# remotes::install_github("mojaveazure/seurat-disk")
# 
# ## install mapscvi
# library(devtools)
# devtools::install_github("lsteuernagel/mapscvi")

#load libraries
library(cowplot)
library(Seurat)
library(mapscvi)
library(tidyverse)
library(ggalluvial)
library(biomaRt)

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

## genotype
burtoni.snseq.combined.sct = AddMetaData(
  object = burtoni.snseq.combined.sct,
  metadata = burtoni.souporcell.filtered %>% 
    select(Genotype.id,
           Cell.id) %>% 
    column_to_rownames(var = "Cell.id"),
  col.name = 'Genotype.id'
)


#### Map homologous genes from human ####
### Get required gene list
reference_gene_list=mapscvi::reference_hypoMap_downsample@assays$RNA %>% 
  rownames() 
###selecting biomart database
ensembl = useMart('ensembl')
#list datasets
dataset = listDatasets(ensembl)
##select dataset
#mouse
ensembl.mouse = useDataset('mmusculus_gene_ensembl',
                           mart=ensembl)

# get list of all human attributes in tilapia 
# listAttributes(ensembl.human) %>%
#   filter(grepl('onil',
#                name)) %>%
#   pull(name)

#create human attributes
mouse.to.tilapia.attributes = c('external_gene_name',
                                'ensembl_gene_id',
                                'oniloticus_homolog_ensembl_gene',
                                'oniloticus_homolog_associated_gene_name',
                                'oniloticus_homolog_orthology_type',
                                'oniloticus_homolog_orthology_confidence')

#identify human to tilapia homology
human.to.tilapia.homology = getBM(attributes = human.to.tilapia.attributes,
                                  mart = ensembl.human,
                                  values = reference_gene_list,
                                  filter = 'external_gene_name',
                                  useCache = FALSE) # useCache has to do with version of R not being up to date?

#count number of human genes?
#1287
human.to.tilapia.homology %>% 
  na.omit() %>% 
  nrow() 

#check number of tilapia genes?
#497
human.to.tilapia.homology %>% 
  na.omit() %>% 
  pull(ensembl_gene_id) %>% 
  unique() %>% 
  length()

# add tilapia gene names
human.to.tilapia.homology = human.to.tilapia.homology %>% 
  right_join(features.tsv %>%  
               dplyr::select(-c(V3)) %>% 
               dplyr::rename('oniloticus_homolog_ensembl_gene' = 'V1')%>% 
               dplyr::rename('tilapia_gene_name' = 'V2')) %>% 
  na.omit()

# keep only 1:1 orthologs
# create list of single genes
one2one.human <- human.to.tilapia.homology %>% 
  group_by(external_gene_name) %>% 
  summarise(n()) %>% 
  filter(`n()` <= 1) %>%
  pull(external_gene_name)

# filter to only orthologs 
# remove NA 
human.to.tilapia.orthologs <- human.to.tilapia.homology %>% 
  filter(external_gene_name %in% one2one.human) %>% 
  na.omit()


# remove multiple human genes 
one2many.human <- human.to.tilapia.orthologs %>% 
  group_by(oniloticus_homolog_associated_gene_name) %>% 
  summarise(n()) %>% 
  filter(`n()` <= 1) %>%
  pull(oniloticus_homolog_associated_gene_name)

# remove duplicate human genes
human.to.tilapia.orthologs <- human.to.tilapia.orthologs %>% 
  filter(oniloticus_homolog_associated_gene_name %in% one2many.human)




#### test ####
tmp = mapscvi::map_new_seurat_hypoMap(burtoni.snseq.combined.sct,
                                       suffix="query_burtoni",
                                       max_epochs=20)

names(test@reductions)


plot_query_labels(query_seura_object=test,
                  reference_seurat=mapscvi::reference_hypoMap_downsample,
                  overlay = FALSE,
                  labelonplot = FALSE)


tmp2 = prepare_query(burtoni.snseq.combined.sct,
                     assay = 'SCT',
                                      suffix="query_burtoni",
                                      normalize=FALSE)


tmp2 = predict_query(tmp2,
                     model_path = system.file("extdata/models/hypoMap_harmonized_scVI_model/", 
                                              package = 'mapscvi'),
                     max_epochs = 20)

names(tmp2@reductions)



cluster_labels = mapscvi::reference_hypoMap_downsample@meta.data$C66_named
reference_reduction = "scvi"
tmp2 = project_query(query_seurat_object = tmp2,
                                      reference_map_reduc = mapscvi::reference_hypoMap_downsample@reductions[[reference_reduction]],
                                      reference_map_umap = mapscvi::reference_hypoMap_downsample@reductions[[paste0("umap_",reference_reduction)]],
                                      query_reduction = "scvi",
                                      label_vec =cluster_labels)
