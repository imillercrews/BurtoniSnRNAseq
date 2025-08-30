#### Burtoni snseq 
### neuropeptides gene lists
## use lambcomp1 to run R with command 

## load library
library(UpSetR)

#### load gene lists ####
## www.neuropeptides.nl
neuropeptides.nl.list= read.csv('../Gene.lists/Nueropeptides.nl_list.csv')
#reduce
neuropeptides.nl.list.reduce = neuropeptides.nl.list %>% 
  select(Gene) %>% 
  rename(Gene.name = Gene) %>% 
  unique()


## HGNC
## https://www.genenames.org/data/genegroup/#!/group/542
HGNC.list= read.csv('../Gene.lists/HGNC_receptor_ligands.csv')

# subset neuropeptides
# reduce
HGNC.list.reduce = HGNC.list %>% 
  filter(Group.name == 'Neuropeptides') %>% 
  select(Approved.symbol) %>% 
  rename(Gene.name = Approved.symbol) %>% 
  unique()

## NeuroPedia: Neuropeptide database and spectra library
## http://proteomics.ucsd.edu/Software/NeuroPedia/
NeuroPedia.list= read.csv('../Gene.lists/Database_NeuroPedia_063011_reduced.csv',
                          header = T)
# reduce
NeuroPedia.list.reduce = NeuroPedia.list %>% 
  select(Gene.Name) %>% 
  rename(Gene.name = Gene.Name) %>% 
  unique()



#### Combine neuropeptide lists ####

neuropeptides.list = full_join(neuropeptides.nl.list.reduce ,
                               HGNC.list.reduce) %>% 
  full_join(NeuroPedia.list.reduce)

#save
write_csv(neuropeptides.list,
          "../Gene.lists/neuropeptides.list.csv")

#### compare lists ####
### create overlap matrix 
neuropeptides.matrix = full_join(neuropeptides.nl.list.reduce %>% 
                                   mutate(neuropeptides.nl.list.reduce = 1),
                               HGNC.list.reduce %>% 
                                 mutate(HGNC.list.reduce = 1)) %>% 
  full_join(NeuroPedia.list.reduce%>% 
              mutate(NeuroPedia.list.reduce = 1)) %>% 
  select(!('Gene.name')) %>%
  mutate(across(everything(), 
                replace_na, 
                0)) 

## graph
png(file = '../Gene.lists/Neuropeptide gene list overlap.png',
    width = 10,
    height = 10,
    units = 'in',
    res = 480)
upset(neuropeptides.matrix, 
      order.by = "freq",
      text.scale = c(2, 2, 2, 2, 2, 2))
dev.off()


