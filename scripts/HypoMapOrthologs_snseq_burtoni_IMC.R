#### Burtoni snseq seurat analysis
### HypoMap orthologs
###Note: Seurat requires R version > 4
## use lambcomp1 to run R with command 
# > R-4.0.3
###CellCall
#https://www.nature.com/articles/s42255-022-00657-y

### set working directory
setwd("/stor/work/Hofmann/All_projects/A_burtoni_snseq/seurat/")


#### load libraries ####

#load libraries
library(biomaRt)
library(tidyverse)

#### load data ####
## gene markers from HypoMap supplementary table 5
## copied mouse marker genes for these groups:  
#C7-1: GLU,C7-2: GABA, C7-3: Astro-Ependymal, C7-4: Oligo+Precursor, C7-5: Immune, C7-6: ParsTuber, C7-7: Vascular
## load HypoMap mouse marker genes
mouse.annotation.data = read.csv("../Gene.lists/Cell.type.gene.markers.HypoMap.csv")

#get list of genes
features.tsv = read.delim("../cellranger/count/NileTilapia.reference/dom/Burtoni-dom/outs/filtered_feature_bc_matrix/features.tsv.gz", 
                          header = FALSE,
                          stringsAsFactors = FALSE)

## load HypoMap neuron mouse marker genes
mouse.annotation.data = read.csv("../Gene.lists/Cell.type.gene.markers.HypoMap.csv")

#### Map homologous genes from mouse ####
###selecting biomart database
ensembl = useMart('ensembl')
#list datasets
dataset = listDatasets(ensembl)
##select dataset
#mouse
ensembl.mouse = useDataset('mmusculus_gene_ensembl',
                           mart=ensembl)

# get list of all mouse attributes in tilapia 
listAttributes(ensembl.mouse) %>%
  filter(grepl('onil',
               name)) %>%
  pull(name)

#create mouse attributes
mouse.to.tilapia.attributes = c('external_gene_name',
                                'ensembl_gene_id',
                                'oniloticus_homolog_ensembl_gene',
                                'oniloticus_homolog_associated_gene_name',
                                'oniloticus_homolog_orthology_type',
                                'oniloticus_homolog_orthology_confidence')

#identify mouse to tilapia homology
mouse.to.tilapia.homology = getBM(attributes = mouse.to.tilapia.attributes,
                                  mart = ensembl.mouse,
                                  values = mouse.annotation.data$gene,
                                  filter = 'external_gene_name',
                                  useCache = FALSE) # useCache has to do with version of R not being up to date?

#count number of mouse genes
#306
mouse.annotation.data %>% 
  nrow() 

# number of 
#347
mouse.to.tilapia.homology %>% 
  na.omit() %>% 
  nrow() 

# summarize data
# need to annotate 77 genes! (424-347 = 77)
#NA                 ortholog_many2many ortholog_one2many ortholog_one2one 
#77                 20                 73                134 
mouse.to.tilapia.homology %>% 
  dplyr::select(c(oniloticus_homolog_orthology_type,
           external_gene_name)) %>% 
  distinct() %>% 
  dplyr::select(-c(external_gene_name)) %>% 
  table()

## add cluster information
mouse.to.tilapia.homology = mouse.to.tilapia.homology %>% 
  full_join(mouse.annotation.data %>% 
              dplyr::rename('external_gene_name' = 'gene'))

# add tilapia gene names
mouse.to.tilapia.homology = mouse.to.tilapia.homology %>% 
  left_join(features.tsv %>%  
               dplyr::select(-c(V3)) %>% 
               dplyr::rename('oniloticus_homolog_ensembl_gene' = 'V1')%>% 
               dplyr::rename('tilapia_gene_name' = 'V2')) 

### save file
write.csv(mouse.to.tilapia.homology, 
          file = "../Gene.lists/Cell.type.gene.markers.HypoMap.annotation.csv",
          row.names = F)







