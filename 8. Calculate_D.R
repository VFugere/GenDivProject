## Vincent Fugere March 20 2019

# This code calculates mean spatial (Great Circle) distance among a set of coordinates
# in a .coords file, looping along all .coords files in a defined set of folders
# returns for all populations 'D' sensu Gratton et al. 2017 TREE

rm(list=ls())

library(tidyverse)
library(sp)

#output data frame
mean.dist <- data.frame('pop' = character(0), 'year' = numeric(0), 'D' = numeric(0),
                        'taxon' = character(0), 'scale' = character(0),
                        'species' = character(0), stringsAsFactors = F)

#path containing all subfolders with .coords files
base.dir <- ('/Users/vincentfugere/Desktop/Data/coord_files/')

#name of subfolders
folders <- list.files(base.dir)

for(folder in folders){
  
  #scale specific subfolders
  subfolders <- list.files(paste0(base.dir,folder))
  subfolders <- subfolders[subfolders != 'ungrouped']
  
  for(subfolder in subfolders){
    
    #list files
    dir <- paste0(base.dir,folder,'/',subfolder)
    files <- list.files(dir)
    
    #loop through files and calculate spatial distance
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

DF_D <- mean.dist

save(DF_D, file = '~/Desktop/DF_D.RData')

# pdf('~/Desktop/mean_dist.pdf',width=4,height=3.5,pointsize = 8)
# hist(mean.dist$mean.dist,breaks=100, col = 'gray', border = NA, main=NULL,
#      xlab = 'mean pairwise Great Circle distance among sequences (km)',
#      ylab = 'number of species by year combinations')
# dev.off()

DF_D %>% group_by(scale,taxon) %>% summarise(na_count = sum(is.na(D)))
