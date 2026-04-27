# Project Description

Includes all relevant code for https://doi.org/10.1101/2024.11.07.622449

# Script descriptions

## Initial data processing

### command line
README_basespace_IMC.txt - download data files

README_cellranger_IMC.txt - process single cell files through 10x genomics

README_souporcell_IMC.txt - genotype-free demultiplexing

### R
sample_metadata_snseq_burtoni_IMC.R - information on samples

HypoMapOrthologs_snseq_burtoni_IMC.R - get orthologs for comparison with hypomap database

## Seurat analysis in R

### Normaliztion, integration, clustering
seurat_snseq_burtoni_IMC.R 

### cell identification
souporcell_snseq_burtoni_IMC.R - genotyping data

sctype_snseq_burtoni_IMC.R - cell type calling

HypoMap_snseq_burtoni_IMC.R - comparison with HypoMap database

### DEG analysis
AVP_DEG_snseq_burtoni_IMC.R - AVP clustering analysis and DEG

Astrocytes_DEG_snseq_burtoni_IMC.R - Astrocyte specific clustering and DEG

Oligodendrocytes_DEG_snseq_burtoni_IMC.R - Oligodendrocyte specific clustering and DEG

ParsTuber_DEG_snseq_burtoni_IMC.R - ParsTuber specific clustering and DEG

Vascular_DEG_snseq_burtoni_IMC.R - Vascular specific clustering and DEG

Cell_comparison_DEG_snseq_burtoni_IMC.R - compare results from non-neuronal cell DEG

### neuropeptides
Neuropeptide_network_snseq_burtoni_IMC.R - neuropeptide analysis and networks

# neuron
hdWGCNA_snseq_burtoni_IMC.R - neuron analysis and hdWGCNA 
