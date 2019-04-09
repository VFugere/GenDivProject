## Vincent Fugere 2019

# This code unites population-level datasets and create a master data file for analysis
# Also explores the number of populations and times series at each scale

rm(list=ls())
library(tidyverse)
library(vegan)

load('~/Desktop/Data/DF_D.Rdata')
load('~/Desktop/Data/DF_Pi.Rdata')
load('~/Desktop/Data/sequence_metadata.RData')

DF_Pi <- select(DF_Pi, -(taxon:species))
DF_Pi$year <- as.numeric(DF_Pi$year)

alldata <- inner_join(DF_D, DF_Pi, by = c('pop','year'))

#how many pops per scale?
alldata %>% group_by(scale) %>% summarize(pops = n_distinct(pop))
#10K vs. no treshold have almost the same number of pops (32 pops diff)
#Biggest jump is from 1K to 100, loosing 5K pops

alldata <- alldata %>% group_by(pop) %>% add_tally(name = 'n.years')

#how many time series?
alldata %>% filter(n.years >= 3) %>%
  group_by(scale, taxon) %>%
  summarize('pops' = n_distinct(pop), 'species' = n_distinct(species))
#plenty! and number of species ~ number of pops even at smallest scales, indicating
#that we are not simply cutting pops into small units

#get spatial centroids of pops
sc10 <- seq %>% group_by(pop10) %>% summarize('nseqs' = n(),'lat' = mean(lat), 'long' = mean(long)) %>% rename('pop' = pop10)
sc100 <- seq %>% group_by(pop100) %>% summarize('nseqs' = n(),'lat' = mean(lat), 'long' = mean(long)) %>% rename('pop' = pop100)
sc1000 <- seq %>% group_by(pop1000) %>% summarize('nseqs' = n(),'lat' = mean(lat), 'long' = mean(long)) %>% rename('pop' = pop1000)
sc10000 <- seq %>% group_by(pop10000) %>% summarize('nseqs' = n(),'lat' = mean(lat), 'long' = mean(long)) %>% rename('pop' = pop10000)
sc100000 <- seq %>% group_by(pop100000) %>% summarize('nseqs' = n(),'lat' = mean(lat), 'long' = mean(long)) %>% rename('pop' = pop100000)
centroids <- bind_rows(sc10,sc100,sc1000,sc10000,sc100000) %>% filter(nseqs > 1) %>% select(-nseqs)
rm(sc10,sc100,sc1000,sc10000,sc100000)
alldata <- left_join(alldata, centroids, by = 'pop')

#add human impact variables
sc10 <- seq %>% group_by(pop10,year) %>% summarize('nseqs' = n(),'hd' = mean(pop_tot), 'hd.var' = var(pop_tot), 'p.lu' = mean(conv_rangeland+cropland+pasture+urban,na.rm = T), 'lu.var' = var(conv_rangeland+cropland+pasture+urban,na.rm = T)) %>% rename('pop' = pop10)
sc100 <- seq %>% group_by(pop100,year) %>% summarize('nseqs' = n(),'lat' = mean(lat), 'long' = mean(long)) %>% rename('pop' = pop100)
sc1000 <- seq %>% group_by(pop1000,year) %>% summarize('nseqs' = n(),'lat' = mean(lat), 'long' = mean(long)) %>% rename('pop' = pop1000)
sc10000 <- seq %>% group_by(pop10000,year) %>% summarize('nseqs' = n(),'lat' = mean(lat), 'long' = mean(long)) %>% rename('pop' = pop10000)
sc100000 <- seq %>% group_by(pop100000,year) %>% summarize('nseqs' = n(),'lat' = mean(lat), 'long' = mean(long)) %>% rename('pop' = pop100000)
human.impacts <- bind_rows(sc10,sc100,sc1000,sc10000,sc100000) %>% filter(nseqs > 1) %>% select(-nseqs)
rm(sc10,sc100,sc1000,sc10000,sc100000)
alldata <- left_join(alldata, human.impacts, by = 'pop')
