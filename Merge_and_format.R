## Vincent Fugere 2019

# This code unites the 3 population-level dataset: π, D, and
# land use and human density variables.
# Also explores the number of populations and times series at each scale

rm(list=ls())
library(tidyverse)

load('~/Desktop/Data/DF_D.Rdata')
load('~/Desktop/Data/DF_Pi.Rdata')

filter(DF_D, species == "Scatopsciara_atomaria")
filter(DF_Pi, species == "Scatopsciara_atomaria")

DF_Pi <- select(DF_Pi, -(taxon:species))

with(subset(DF_Pi, species == 'Scatopsciara_atomaria'), table(year,pop))
with(subset(DF_D, species == 'Scatopsciara_atomaria'), table(year,pop))

alldata <- inner_join(DF_D, DF_Pi, by = c('pop','year'))

#how many pops per scale?
alldata %>% group_by(scale) %>% summarize(pops = n_distinct(pop))
#10K vs. no treshold have almost the same number of pops (32 pops diff)
#Biggest jump is from 1K to 100, loosing 5K pops

alldata <- alldata %>% group_by(pop) %>% add_tally(name = 'nyrs')
hist(alldata$nyrs,breaks=100)

#?

lrgsc <- filter(alldata, scale == '1e+05')
  
#simulating data while waiting for Chloe's file

