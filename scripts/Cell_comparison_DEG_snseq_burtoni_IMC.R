#### Burtoni snseq seurat analysis
### cell type DEG comparison
###Note: Seurat requires R version > 4
## use lambcomp1 to run R with command 
# > R-4.0.3



### set working directory
setwd("/stor/work/Hofmann/All_projects/A_burtoni_snseq/seurat/")

#### load libraries ####

#load libraries
library(tidyverse)

# using the cowplot theme for ggplot
# theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

#### load data ####
### get limmatrend DEG data for each cell type
data.vascular = read.csv('vascular/limmatrend/vascular.limmatrend.results.csv')
data.parstuber= read.csv('parstuber/limmatrend/parstuber.limmatrend.results.csv')
data.astrocytes = read.csv('astrocytes/limmatrend/astrocytes.limmatrend.results.csv')
data.oligodendrocytes = read.csv('oligodendrocytes/limmatrend/oligodendrocytes.limmatrend.results.csv')

### get genes to keep
data.vascular.keep = data.vascular %>% 
  filter(Keep == 'Keep') %>% 
  mutate(vascular = 1) %>% 
  dplyr::select(Gene,
                vascular)
data.parstuber.keep = data.parstuber %>% 
  filter(Keep == 'Keep') %>% 
  mutate(parstuber = 1) %>% 
  dplyr::select(Gene,
                parstuber)
data.astrocytes.keep = data.astrocytes %>% 
  filter(Keep == 'Keep') %>% 
  mutate(astrocytes = 1) %>% 
  dplyr::select(Gene,
                astrocytes)
data.oligodendrocytes.keep = data.oligodendrocytes %>% 
  filter(Keep == 'Keep') %>% 
  mutate(oligodendrocytes = 1) %>% 
  dplyr::select(Gene,
                oligodendrocytes)
# combine
data.keep = data.vascular.keep %>% 
  full_join(data.parstuber.keep) %>% 
  full_join(data.astrocytes.keep) %>% 
  full_join(data.oligodendrocytes.keep)

# convert NA to zero
data.keep = data.keep %>% 
  replace(is.na(.), 0)

# create total
data.keep = data.keep %>% 
  mutate(total = rowSums(across(where(is.numeric))))

# get table
data.keep.table = data.keep %>% 
  dplyr::select(total) %>% 
  table()

# genes present in 3 or 4 comparisons 
data.keep.gene.list = data.keep %>% 
  filter(total >=3) %>% 
  pull(Gene)

### compare genes in heat map
data.keep.value = data.vascular %>% 
  filter(Sig_DvsS == 'Sig') %>%
  dplyr::select(Gene,
                All) %>% 
  dplyr::rename(vascular = All) %>% 
  full_join(data.parstuber %>% 
              filter(Sig_DvsS == 'Sig') %>%
              dplyr::select(Gene,
                            All) %>% 
              dplyr::rename(parstuber = All)) %>% 
  full_join(data.astrocytes %>% 
              filter(Sig_DvsS == 'Sig') %>%
              dplyr::select(Gene,
                            All) %>% 
              dplyr::rename(astrocytes = All)) %>% 
  full_join(data.oligodendrocytes %>% 
              filter(Sig_DvsS == 'Sig') %>%
              dplyr::select(Gene,
                            All) %>% 
              dplyr::rename(oligodendrocytes = All)) %>%
  replace(is.na(.), 0)

# replace infinite values
data.keep.value = data.keep.value %>% 
  mutate(vascular = ifelse(vascular == Inf,
                            max(data.keep.value$vascular[is.finite(data.keep.value$vascular)]),
                            vascular)) %>% 
  mutate(parstuber = ifelse(parstuber == Inf,
                            max(data.keep.value$parstuber[is.finite(data.keep.value$parstuber)]),
                            parstuber)) %>% 
  mutate(astrocytes = ifelse(astrocytes == Inf,
                            max(data.keep.value$astrocytes[is.finite(data.keep.value$astrocytes)]),
                            astrocytes)) %>% 
  mutate(oligodendrocytes = ifelse(oligodendrocytes == Inf,
                            max(data.keep.value$oligodendrocytes[is.finite(data.keep.value$oligodendrocytes)]),
                            oligodendrocytes))


# filter down to top genes
data.keep.value.reduce = data.keep.value %>% 
  filter(Gene %in% data.keep.gene.list) %>% 
  column_to_rownames('Gene') %>% 
  as.matrix()


#### graph cell type comparison ####
### compare gene count
## bar graph
data.keep.table %>% 
  data.frame() %>% 
  ggplot(aes(x = total,
             y = Freq)) +
  geom_bar(stat = 'identity') +
  geom_label(aes(label = Freq)) +
  theme_classic() +
  ylab('Gene count') +
  xlab('Number of cell specific DEG gene is present')
ggsave(filename = "Cell.types.deg/DEGs across celltype comparison.png",
                  width = 5,
                  height = 5)

## upset plot
data.keep

### heat map of top genes
paletteLength <- 25
myColor <- colorRampPalette(c("#60bb46", "white","#4e499e"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(data.keep.value.reduce), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(data.keep.value.reduce)/paletteLength, max(data.keep.value.reduce), length.out=floor(paletteLength/2)))

# Plot the heatmap
data.keep.value.reduce %>% 
  t() %>% 
  pheatmap(color=myColor, 
           breaks=myBreaks,
           angle_col = 45,
           filename = 'Cell.types.deg/DEGs across celltype DvsS heatmap.png',
           height = 5,
           width = 10)

# keep solid
data.keep.value.reduce %>% 
  t() %>% 
  pheatmap(color=c("#60bb46", "white","#4e499e"), 
           breaks = c(min(data.keep.value.reduce),
                      -0.01,
                      0.01,
                      max(data.keep.value.reduce)),
           angle_col = 45,
           filename = 'Cell.types.deg/DEGs across celltype DvsS heatmap solid.png',
           height = 5,
           width = 10)









