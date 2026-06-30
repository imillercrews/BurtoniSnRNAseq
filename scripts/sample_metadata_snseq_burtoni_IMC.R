#### Burtoni snseq 
### sample metdata analysis
## use lambcomp1 to run R with command 
# > R-4.0.3

### set working directory
setwd("/stor/work/Hofmann/All_projects/A_burtoni_snseq/seurat/")

#### load libraries ####
#load libraries
library(tidyverse)

#### load data ####
### sample meta data
sample.data = read.csv('../Burtoni_samples/sample_metadata.csv')


#### Compare measurements across social status ####
### one sample ttest
## calculate percent change
## Setup_SL
sample.data %>% 
  dplyr::select(c('Tank',
         'Dissection_SL',
         'Status')) %>% 
  pivot_wider(id_cols = 'Tank',
              names_from = 'Status',
              values_from = 'Dissection_SL') %>% 
  mutate(Percent.change = 100*(Dom-Sub)/Sub) %>% 
  pull(Percent.change) %>% 
  t.test(mu = 0, 
         alternative = "two.sided")

## Dissection_Mass
sample.data %>% 
  dplyr::select(c('Tank',
                  'Dissection_Mass',
                  'Status')) %>% 
  pivot_wider(id_cols = 'Tank',
              names_from = 'Status',
              values_from = 'Dissection_Mass') %>% 
  mutate(Percent.change = 100*(Dom-Sub)/Sub) %>% 
  pull(Percent.change) %>% 
  t.test(mu = 0, 
         alternative = "two.sided")

## GSI
sample.data %>% 
  dplyr::select(c('Tank',
                  'GSI',
                  'Status')) %>% 
  pivot_wider(id_cols = 'Tank',
              names_from = 'Status',
              values_from = 'GSI') %>% 
  mutate(Percent.change = 100*(Dom-Sub)/Sub) %>% 
  pull(Percent.change) %>% 
  t.test(mu = 0, 
         alternative = "two.sided")

### use  t-test 
## GSI
t.test(GSI ~ Status,
       data = sample.data)

# summary stats
sample.data %>% 
  group_by(Status) %>% 
  summarise(mean = mean(GSI),
            SD = sd(GSI))


## Testosterone
t.test(Testosterone_pg.mL ~ Status,
       data = sample.data)

# summary stats
sample.data %>% 
  group_by(Status) %>% 
  summarise(mean = mean(Testosterone_pg.mL/1000),
            SD = sd(Testosterone_pg.mL/1000))











