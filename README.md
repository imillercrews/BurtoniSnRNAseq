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

HypoMapOrthologs_snseq_burtoni_IMC.R - get orthologs for comparison with hypomap database [Data: S3]

## Seurat analysis in R

### Normaliztion, integration, clustering
seurat_snseq_burtoni_IMC.R 

### cell identification
souporcell_snseq_burtoni_IMC.R - genotyping data

sctype_snseq_burtoni_IMC.R - cell type calling [Fig. 1B, 1C]

### DEG analysis
AVP_DEG_snseq_burtoni_IMC.R - AVP clustering analysis and DEG [Fig: 3, S8, S9; Data: S4]

Astrocytes_DEG_snseq_burtoni_IMC.R - Astrocyte specific clustering and DEG [Fig: S2; Data: S4]

Oligodendrocytes_DEG_snseq_burtoni_IMC.R - Oligodendrocyte specific clustering and DEG [Fig: S4; Data: S4]

ParsTuber_DEG_snseq_burtoni_IMC.R - ParsTuber specific clustering and DEG [Fig: S5; Data: S4]

Vascular_DEG_snseq_burtoni_IMC.R - Vascular specific clustering and DEG [Fig: S3; Data: S4]

Cell_comparison_DEG_snseq_burtoni_IMC.R - compare results from non-neuronal cell DEG [Fig: S6]

### neuropeptides
Neuropeptide_network_snseq_burtoni_IMC.R - neuropeptide analysis and networks [Fig: 4, S10]

# neuron
hdWGCNA_snseq_burtoni_IMC.R - neuron analysis and hdWGCNA [Fig: 1D, 2, S1, S7; Data: S1, S2]
