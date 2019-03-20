## Vincent Fugère, March 20 2019

# This code takes all coords files that belong to the same species,
# pastes them together, calculates a Great Circle distance matrix among all points,
# runs a cluster analysis to find gaps in space of minimum distance x,
# then assign groups to the coords based on clusters, 
# then re-split the lines into their original coord files now with a new
# column 'population' that is used to decide what sequences are being
# grouped together in the other analyses (pi hat, D, and LU mean and var)

library(tidyverse)
library(sp)

#master parameter: min distances among groups of sequences to be considered separate pops
#in km (Great Circle distance)

scales <- c(10,100,1000,10000,100000)

for(scale in scales){
  
  minim.dist <- scale
  
  #path containing all subfolders with .coords files
  base.dir <- ('/Users/vincentfugere/Desktop/Data/coord_files/')
  
  #name of subfolders (taxa)
  folders <- list.files(base.dir)
  
  for(folder in folders){
    
    dir <- paste0(base.dir,folder)
    output.dir <- paste0(base.dir,folder,'/grouped_mindist',minim.dist)
    dir.create(output.dir)
    dir <- paste0(base.dir,folder,'/ungrouped')
    
    #list species
    species <- unique(str_sub(list.files(dir, pattern = ".coords"), end = -13))
    
    #loop through species
    for(sp in species){
      
      #list files
      files <- list.files(dir, pattern = sp)
      
      #load files and paste them together
      outdat <- data.frame('file' = character(0), 'V1' = numeric(0), 'V2' = numeric(0), stringsAsFactors = F)
      
      for(file in files){
        dat <- read.table(file.path(dir, file))
        name.vec <- data.frame('file' = rep(file,nrow(dat)))
        dat <- cbind(name.vec,dat)
        outdat <- rbind(outdat,dat)
      }
      
      if(nrow(outdat) == 1){
        next
      }
      
      #create groups/clusters
      d <- spDists(as.matrix(outdat[,2:3]), longlat=TRUE)
      d <- as.dist(d)
      hc <- hclust(d, method = 'single')
      outdat$group <- paste(sp,scale,cutree(hc, h = minim.dist),sep='_')
      
      #output individual files
      outdat <- outdat[,c(4,2,3,1)]
      for(file in files){
        dat <- outdat[outdat$file == file,1:3]
        fname <- str_remove(file,'.coords')
        write.table(dat, file = file.path(output.dir, fname), sep = "\t", col.names = F, row.names = F)
      }
      
    }
    
  }
  
}