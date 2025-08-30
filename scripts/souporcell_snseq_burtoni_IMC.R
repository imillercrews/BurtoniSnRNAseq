#### Burtoni snseq genotyping analysis
### souporcell
#https://github.com/wheaton5/souporcell

### set working directory
setwd("/stor/work/Hofmann/All_projects/A_burtoni_snseq/souporcell/")

#### load libraries ####
##install libraries
# install.packages("tidyverse",
#                  repos = "https://cloud.r-project.org")
# install.packages("Seurat",
#                  repos = "https://cloud.r-project.org")
#load libraries
library(Seurat)
library(tidyverse)

###ideas
## Assignment/status vs UMI/reads 
# Do doublets have high number of reads? 
# Do unassigned have low number of UMI?
## Select filtered cells, look at assignment and status
# Do any doublets or unassigned make the cut?
# Is there a bias for a specific assignment?
## Assignment of filtered cells for cluster?
# Is there a bias for any clusters?
## compare dom to sub
# ambient RNA amount: ambient_rna.txt

#### load data ####
### load souporcell data
##dom
dom.clusters = read_tsv('../souporcell/dom_souporcell/clusters.tsv')
#rename barcode to Cell.id
dom.clusters = dom.clusters %>% 
  rename(Cell.id = barcode) %>% 
  mutate(Cell.id = paste(Cell.id, 
                         '1',
                         sep = '_'),
         orig.ident = 'dom_burtoni_snseq')

##sub
sub.clusters = read_tsv('../souporcell/sub_souporcell/clusters.tsv')
#rename barcode to Cell.id
sub.clusters = sub.clusters %>% 
  rename(Cell.id = barcode) %>% 
  mutate(Cell.id = paste(Cell.id, 
                         '2',
                         sep = '_'),
         orig.ident = 'sub_burtoni_snseq')

## merge dom and sub data
burtoni.souporcell = full_join(dom.clusters,
                               sub.clusters)

# ### select data
# ##dom
# dom.clusters.select = read_tsv('dom_souporcell_select/clusters.tsv')
# #rename barcode to Cell.id
# dom.clusters.select = dom.clusters.select %>% 
#   rename(Cell.id = barcode) %>% 
#   mutate(Cell.id = paste(Cell.id, 
#                          '1',
#                          sep = '_'),
#          orig.ident = 'dom_burtoni_snseq')
# 
# ##sub
# sub.clusters.select = read_tsv('sub_souporcell_select/clusters.tsv')
# #rename barcode to Cell.id
# sub.clusters.select = sub.clusters.select %>% 
#   rename(Cell.id = barcode) %>% 
#   mutate(Cell.id = paste(Cell.id, 
#                          '2',
#                          sep = '_'),
#          orig.ident = 'sub_burtoni_snseq')
# 
# ### merge dom and sub data
# burtoni.souporcell.select = full_join(dom.clusters.select,
#                                sub.clusters.select)

### scsorter data
load("../seurat/burtoni.scsorter.data.scsort.output.RData")

### combine data
Souporcell.cell.type = burtoni.souporcell %>% 
  left_join(burtoni.scsorter.data.scsort.output)

# #select
# Souporcell.cell.type.select = burtoni.souporcell.select %>% 
#   left_join(burtoni.scsorter.data.scsort.output)

#assign genotype names 
# burtoni.souporcell.filtered = Souporcell.cell.type %>% 
#   select(c(Cell.id, status, assignment, log_prob_singleton, log_prob_doublet, cluster0, cluster1, cluster2, orig.ident, Genotype.id, Cell.type)) %>% 
#   mutate(Filtered = ifelse(is.na(Cell.type),
#                            'Removed',
#                            'Selected')) %>% 
#   mutate(Genotype.id = case_when(Filtered == 'Selected' & assignment == '1' & orig.ident == 'dom_burtoni_snseq' ~ 'Dom.A',
#                                  Filtered == 'Selected' & assignment == '2' & orig.ident == 'dom_burtoni_snseq' ~ 'Dom.B',
#                                  Filtered == 'Selected' & assignment != '1' & assignment != '2' & orig.ident == 'dom_burtoni_snseq' ~ 'Dom.C',
#                                  Filtered == 'Selected' & assignment == '1' & orig.ident == 'sub_burtoni_snseq' ~ 'Sub.A',
#                                  Filtered == 'Selected' & assignment == '2' & orig.ident == 'sub_burtoni_snseq' ~ 'Sub.B',
#                                  Filtered == 'Selected' & assignment != '1' & assignment != '2' & orig.ident == 'sub_burtoni_snseq' ~ 'Sub.C',
#                                  TRUE ~ 'NA')) %>% 
#   filter(Filtered == 'Selected') %>% 
#   select(-c(Filtered))
# 
# # save data
# write_csv(burtoni.souporcell.filtered,
#           'burtoni.souporcell.filtered.csv')

# read in data
burtoni.souporcell.filtered = read_csv('burtoni.souporcell.filtered.csv')

#check data
burtoni.souporcell.sample.count = burtoni.souporcell.filtered %>% 
  select(c(status, 
           orig.ident,
           Genotype.id)) %>% 
  table() %>% 
  as.data.frame() %>% 
  mutate(assigned = ifelse(status == 'unassigned',
                           'unassigned',
                           'assigned')) %>% 
  select(-c(status)) %>% 
  group_by(Genotype.id) %>% 
  mutate(total = sum (Freq)) %>% 
  ungroup() %>% 
  mutate(Proportion = Freq/total)


#### Graph select ####
# ### sample status
# ##total
# # Status cell count
# Souporcell.cell.type.select %>% 
#   select(status,
#          orig.ident) %>% 
#   mutate(Count = 1) %>% 
#   group_by(status,
#            orig.ident) %>% 
#   summarise(Total = sum(Count)) %>% 
#   ggplot(aes(y = Total,
#              x = status)) + 
#   geom_point(size = 5) +
#   theme_classic() +
#   ggtitle('status souporcell select') + 
#   facet_grid(.~orig.ident)
# ggsave('figures/select/status souporcell select.png')
# 
# # status by assignment
# Souporcell.cell.type.select %>% 
#   select(c(status,
#            assignment,
#            orig.ident)) %>% 
#   mutate(Count = 1) %>% 
#   group_by(assignment,
#            status,
#            orig.ident) %>% 
#   summarise(Total = sum(Count)) %>% 
#   ggplot(aes(y = Total,
#              x = assignment,
#              color = status)) + 
#   geom_jitter(height = 0,
#               width = 0.1,
#               size = 5) +
#   theme_classic() +
#   ggtitle('status by assignment souporcell select') +
#   facet_grid(.~orig.ident)
# ggsave('figures/select/status by assignment souporcell select.png')
# 
# ## cell type
# Souporcell.cell.type.select %>% 
#   mutate(Count = 1) %>% 
#   group_by(assignment,
#            status,
#            Cell.type,
#            orig.ident) %>% 
#   summarise(Total = sum(Count)) %>% 
#   ggplot(aes(y = Cell.type,
#              x = Total,
#              shape = status,
#              color = assignment)) + 
#   geom_jitter(height = 0,
#               width = 0.05,
#               size = 5) +
#   theme_classic() +
#   ggtitle('cell type by assignment and status souporcell select') + 
#   facet_grid(.~orig.ident)
# ggsave('figures/select/cell type by assignment and status souporcell select.png')

#### Graph all ####
### sample status
##total
# Status cell count
Souporcell.cell.type %>% 
  select(status,
         orig.ident) %>% 
  mutate(Count = 1) %>% 
  group_by(status,
           orig.ident) %>% 
  summarise(Total = sum(Count)) %>% 
  ggplot(aes(y = Total,
             x = status)) + 
  geom_point(size = 5) +
  theme_classic() +
  ggtitle('status souporcell') + 
  facet_grid(.~orig.ident)
ggsave('../souporcell/figures/status souporcell.png')

# status by assignment
Souporcell.cell.type %>% 
  select(c(status,
           assignment,
           orig.ident)) %>% 
  mutate(Count = 1) %>% 
  group_by(assignment,
           status,
           orig.ident) %>% 
  summarise(Total = sum(Count)) %>% 
  ggplot(aes(y = Total,
             x = assignment,
             color = status)) + 
  geom_jitter(height = 0,
              width = 0.1,
              size = 5) +
  theme_classic() +
  ggtitle('status by assignment souporcell') +
  facet_grid(.~orig.ident)
ggsave('../souporcell/figures/status by assignment souporcell.png')


## filtering comparison
# Status cell count
Souporcell.cell.type %>% 
  mutate(Filtered = ifelse(is.na(Cell.type),
                           'Removed',
                           'Selected')) %>% 
  mutate(Count = 1) %>% 
  group_by(status,
           Filtered,
           orig.ident) %>% 
  summarise(Total = sum(Count)) %>% 
  ggplot(aes(y = Total,
             x = status,
             shape = Filtered)) + 
  geom_point(size = 5) +
  theme_classic() +
  ggtitle('status souporcell') +
  facet_grid(.~orig.ident)
ggsave('../souporcell/figures/filter across status souporcell.png')

# status by assignment
Souporcell.cell.type %>% 
  mutate(Filtered = ifelse(is.na(Cell.type),
                           'Removed',
                           'Selected')) %>% 
  mutate(Count = 1) %>% 
  group_by(assignment,
           status,
           Filtered,
           orig.ident) %>% 
  summarise(Total = sum(Count)) %>% 
  ggplot(aes(y = Total,
             x = assignment,
             color = status,
             shape = Filtered)) + 
  geom_jitter(height = 0,
              width = 0.05,
              size = 5) +
  theme_classic() +
  ggtitle('status by assignment souporcell') +
  facet_grid(.~orig.ident)
ggsave('../souporcell/figures/filter acrossstatus by assignment souporcell.png')

## filtered
# Status cell count
Souporcell.cell.type %>% 
  mutate(Filtered = ifelse(is.na(Cell.type),
                           'Removed',
                           'Selected')) %>% 
  filter(Filtered == 'Selected') %>% 
  mutate(Count = 1) %>% 
  group_by(status,
           orig.ident) %>% 
  summarise(Total = sum(Count)) %>% 
  ggplot(aes(y = Total,
             x = status)) + 
  geom_point(size = 5) +
  theme_classic() +
  ggtitle('status souporcell filtered') +
  facet_grid(.~orig.ident)
ggsave('../souporcell/figures/filtered status souporcell.png')

# status by assignment
Souporcell.cell.type %>% 
  mutate(Filtered = ifelse(is.na(Cell.type),
                           'Removed',
                           'Selected')) %>% 
  filter(Filtered == 'Selected') %>% 
  mutate(Count = 1) %>% 
  group_by(assignment,
           status,
           Filtered,
           orig.ident) %>% 
  summarise(Total = sum(Count)) %>% 
  ggplot(aes(y = Total,
             x = assignment,
             color = status)) + 
  geom_jitter(height = 0,
              width = 0.05,
              size = 5) +
  theme_classic() +
  ggtitle('status by assignment souporcell filtered') + 
  facet_grid(.~orig.ident)
ggsave('../souporcell/figures/filtered status by assignment souporcell.png')


## cell type
Souporcell.cell.type %>% 
  mutate(Filtered = ifelse(is.na(Cell.type),
                           'Removed',
                           'Selected')) %>% 
  filter(Filtered == 'Selected') %>% 
  mutate(Count = 1) %>% 
  group_by(assignment,
           status,
           Filtered,
           Cell.type,
           orig.ident) %>% 
  summarise(Total = sum(Count)) %>% 
  ggplot(aes(y = Cell.type,
             x = Total,
             shape = status,
             color = assignment)) + 
  geom_jitter(height = 0,
              width = 0.05,
              size = 5) +
  theme_classic() +
  ggtitle('cell type by assignment and status souporcell filtered') + 
  facet_grid(.~orig.ident)
ggsave('../souporcell/figures/cell type filtered status by assignment souporcell.png')

#singlets
Souporcell.cell.type %>% 
  mutate(Filtered = ifelse(is.na(Cell.type),
                           'Removed',
                           'Selected')) %>% 
  filter(Filtered == 'Selected') %>% 
  mutate(Count = 1) %>% 
  group_by(assignment,
           status,
           Filtered,
           Cell.type,
           orig.ident) %>% 
  summarise(Total = sum(Count)) %>% 
  ggplot(aes(y = Cell.type,
             x = Total,
             shape = status,
             color = assignment)) + 
  geom_jitter(height = 0,
              width = 0.05,
              size = 5) +
  theme_classic() +
  ggtitle('cell type by assignment and status souporcell filtered') +
  facet_grid(status~orig.ident)
ggsave('../souporcell/figures/cell type filtered status by assignment souporcell facet.png')

#genotype by celltype
Souporcell.cell.type %>% 
  mutate(Filtered = ifelse(is.na(Cell.type),
                           'Removed',
                           'Selected')) %>% 
  filter(Filtered == 'Selected') %>% 
  mutate(Count = 1,
         status = ifelse(status == 'unassigned',
                         'unassigned',
                         'assigned')) %>% 
  group_by(Genotype.id,
           status,
           Filtered,
           Cell.type,
           orig.ident) %>% 
  summarise(Total = sum(Count)) %>% 
  ggplot(aes(y = Cell.type,
             x = Total,
             shape = status,
             color = Genotype.id)) + 
  geom_jitter(height = 0,
              width = 0.05,
              size = 5) +
  theme_classic() +
  ggtitle('cell type by assignment and status souporcell filtered') + 
  facet_grid(.~orig.ident)
ggsave('../souporcell/figures/cell type filtered status by genotype souporcell.png')

#facet
Souporcell.cell.type %>% 
  mutate(Filtered = ifelse(is.na(Cell.type),
                           'Removed',
                           'Selected')) %>% 
  filter(Filtered == 'Selected') %>% 
  mutate(Count = 1,
         status = ifelse(status == 'unassigned',
                         'unassigned',
                         'assigned')) %>% 
  group_by(Genotype.id,
           status,
           Filtered,
           Cell.type,
           orig.ident) %>% 
  summarise(Total = sum(Count)) %>% 
  ggplot(aes(y = Cell.type,
             x = Total,
             shape = status,
             color = Genotype.id)) + 
  geom_jitter(height = 0,
              width = 0.05,
              size = 5) +
  theme_classic() +
  ggtitle('cell type by assignment and status souporcell filtered') +
  facet_grid(status~orig.ident)
ggsave('../souporcell/figures/cell type filtered status by genotype souporcell facet.png')

#facet
Souporcell.cell.type %>% 
  mutate(Filtered = ifelse(is.na(Cell.type),
                           'Removed',
                           'Selected')) %>% 
  filter(Filtered == 'Selected') %>% 
  mutate(Count = 1,
         status = ifelse(status == 'unassigned',
                         'unassigned',
                         'assigned')) %>% 
  group_by(Genotype.id,
           Filtered,
           Cell.type,
           orig.ident) %>% 
  summarise(Total = sum(Count)) %>% 
  ggplot(aes(y = Cell.type,
             x = Total,
             color = Genotype.id)) + 
  geom_jitter(height = 0,
              width = 0.05,
              size = 5) +
  theme_classic() +
  ggtitle('cell type by assignment and status souporcell filtered') +
  facet_grid(~orig.ident)
ggsave('../souporcell/figures/cell type filtered status by genotype souporcell facet simple.png')


## graph total assigned per individual
# counts
burtoni.souporcell.sample.count %>% 
  ggplot(aes(y = Genotype.id,
             x = Freq,
             fill = assigned)) +
  geom_bar(position="stack", 
           stat="identity") +
  theme_bw()
ggsave('../souporcell/figures/Genotype ID counts assignment.png')

#proportion
burtoni.souporcell.sample.count %>% 
  ggplot(aes(y = Genotype.id,
             x = Proportion,
             fill = assigned)) +
  geom_bar(position="stack", 
           stat="identity") +
  theme_bw()
ggsave('../souporcell/figures/Genotype ID proportion assignment.png')  

# for poster
burtoni.souporcell.filtered %>% 
  mutate(orig.ident = ifelse(orig.ident == "dom_burtoni_snseq",
                             "Dominant",
                             "Subordinate")) %>% 
  mutate(Count = 1) %>% 
  group_by(Genotype.id,
           orig.ident,
           Cell.type) %>% 
  summarise(Total = sum(Count)) %>%
  ggplot(aes(y = Total,
             x = reorder(Cell.type,
                         -Total),
             fill = orig.ident)) + 
  geom_point(size = 8,
             shape = 21) +
  scale_fill_manual(values = c('Dominant' = '#4e499e',
                                'Subordinate' = '#60bb46')) +
  theme_classic() +
  xlab('') +
  ylab('Total cells assigned')+
  theme_classic() +
  # theme(panel.border = element_rect(color = "black",
  #                                   fill = NA,
  #                                   size = 1))+
  theme(axis.text = element_text(size = 15,
                                 angle = 45,
                                 hjust = 1))  +
  theme(axis.title = element_text(size = 20))+
  theme(legend.position = c(0.7, 
                            0.85),
        legend.box.background = element_rect(colour = "black"),
        legend.background = element_blank(),
        legend.text=element_text(size=20),
        legend.title = element_blank())
ggsave('../souporcell/figures/filtered genotype id by cell type poster.pdf',
       height = 5,
       width = 5,
       units = "in",
       dpi = 300)

## filtered
# Status cell count
burtoni.souporcell.filtered %>% 
  mutate(Count = 1) %>% 
  group_by(Genotype.id,
           orig.ident,
           Cell.type) %>% 
  summarise(Total = sum(Count)) %>% 
  ggplot(aes(y = Total,
             x = Cell.type,
             alpha = Genotype.id,
             color = orig.ident)) + 
  geom_point(size = 5) +
  theme_classic() +
  ggtitle('status souporcell filtered') 
ggsave('../souporcell/figures/filtered status souporcell poster.png')
