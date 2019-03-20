## Vincent Fugere, March 20 2019

## This script will go through .coords and FASTA files and 
## delete files that have bad species names

library(tidyverse)

#bad characters

unwanted <- c('\\.',0:9,'BOLD','-','#','_nsp','spnov')
unwanted <- paste(unwanted, collapse = '|')

#path containing all subfolders with .coords files
coords.dir <- ('/Users/vincentfugere/Desktop/Data/coord_files/')

#path containing all subfolders with FASTA files
fasta.dir <- ('/Users/vincentfugere/Desktop/Data/FASTA_files/')

#name of subfolders (taxa)
folders <- list.files(coords.dir)

for(folder in folders){
  
  dir <- paste0(coords.dir,folder)
  names <- str_sub(list.files(dir, pattern = ".coords"), end = -13)
  
  files <- str_sub(list.files(dir, pattern = ".coords"))
  to.rm <- files[str_count(names, '_') != 1 | str_count(names, unwanted) != 0]
  setwd(dir)
  for(badfile in to.rm){
    file.remove(badfile)
  }
  
  dir <- paste0(fasta.dir,folder)
  files <- str_sub(list.files(dir, pattern = ".fasta"))
  to.rm <- files[str_count(names, '_') != 1 | str_count(names, unwanted) != 0]
  setwd(dir)
  for(badfile in to.rm){
    file.remove(badfile)
  }
  
}
