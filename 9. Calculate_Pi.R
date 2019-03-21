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
base.dir <- ('/Users/vincentfugere/Desktop/Data/Pairwise_sequence_comparisons/')

files <- list.files(base.dir)

for(file in files){
  
      filedat <- read.table(file.path(dir, file), stringsAsFactors = F)
      pops <- unique(filedat$V1)
      
      #calculate for all pops separately
      for(pop in pops){
        
        popdat <- filedat[filedat$V1 == pop,]
        line <- nrow(mean.dist) + 1
        
        mean.dist[line,'pop'] <- popdat[1,1]
        mean.dist[line,'year'] <- as.numeric(str_sub(file, start = -4))
        mean.dist[line,'species'] <- str_sub(file, end = -6)
        mean.dist[line,'taxon'] <- folder
        mean.dist[line,'scale'] <- as.numeric(str_sub(subfolder, start = +16))
        
        if(nrow(popdat)==1){
          mean.dist[line,'D'] <-  NA  
        } else {
          dists <- spDists(as.matrix(popdat[,3:2]), longlat=TRUE)
          dists <- dists[lower.tri(dists, diag=F)]
          mean.dist[line,'D'] <-  mean(dists)
        }
        
      }
      
    }
    
  }
  
}

DF_Pi <- mean.pi

save(DF_Pi, file = '~/Desktop/DF_Pi.RData')
