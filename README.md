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

HypoMapOrthologs_snseq_burtoni_IMC.R - get orthologs for comparison with hypomap database [Data: S5]

## Seurat analysis in R

### Normaliztion, integration, clustering
seurat_snseq_burtoni_IMC.R 

### cell identification
souporcell_snseq_burtoni_IMC.R - genotyping data

sctype_snseq_burtoni_IMC.R - cell type calling [Fig. 1B, 1C]

### AVP analysis
AVP_DEG_snseq_burtoni_IMC.R - AVP clustering analysis and DEG [Fig: 3, S7; Data: S4]

### neuropeptides
Neuropeptide_network_snseq_burtoni_IMC.R - neuropeptide analysis and networks [Fig: 4, S8]

# neuron
hdWGCNA_snseq_burtoni_IMC.R - neuron analysis and hdWGCNA [Fig: 1D, 2, S1-6 S; Data: S1-3]
