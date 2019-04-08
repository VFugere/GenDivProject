## Vincent Fugere March 21 2019

# This code calculates mean genetic distance (π) from pairwise sequence comparions
# computed in Julia. Returns a population-level summary dataframe with population
# and year identifier, and π, sequence number, and comparison number for each pop

rm(list=ls())

library(tidyverse)

#output data frame
mean.pi <- data.frame('pop' = character(0), 'year' = numeric(0),
                      'taxon' = character(0), 'scale' = character(0),
                      'species' = character(0), 'div' = numeric(0),
                      'nseqs' = numeric(0), 'ncomps' = numeric(0),
                      stringsAsFactors = F)

#path containing Julia output files 
dir <- ('/Users/vincentfugere/Desktop/Data/Pairwise_sequence_comparisons/')

files <- list.files(dir)

for(file in files){
  
  filedat <- read_csv(file.path(dir, file)) %>%
    rename(pop = cell) %>%
    filter(overlap >= 0.5) %>% #removes pairwise comparisons with less than 50% sequence overlap
    filter(!is.na(num_per_bp)) %>% #removes sequences alone in their pop, which have NA for seq1:num_per_bp
    separate(col = 'species', into = c('species','year'), sep = '\\.') %>%
    as.data.frame
  
  filedat$pop_yr <- paste(filedat$pop,filedat$year,sep='_') 
  pis <- unique(filedat$pop_yr)
  fileinfo <- str_remove(file, '.csv') %>% str_split('_', simplify = T)
  
  #calculate for all pops separately
  for(pi in pis){
    
    popdat <- filedat[filedat$pop_yr == pi,]
    line <- nrow(mean.pi) + 1
    
    mean.pi[line,'pop'] <- popdat[1,3]
    mean.pi[line,'year'] <- popdat[1,2]
    mean.pi[line,'taxon'] <- fileinfo[1,2]
    mean.pi[line,'scale'] <- fileinfo[1,4]
    mean.pi[line,'species'] <- popdat[1,1]
    
    mean.pi[line,'div'] <- mean(popdat$num_per_bp)
    mean.pi[line,'nseqs'] <- n_distinct(c(popdat[,3],popdat[,4]))
    mean.pi[line,'ncomps'] <- nrow(popdat)
    
  }
  
}

DF_Pi <- mean.pi

save(DF_Pi, file = '~/Desktop/DF_Pi.RData')
