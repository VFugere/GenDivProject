############################################
# Script to plot coordinates on map, and:
# 1) remove coord and fasta lines that land on grid cell border because 
# can't assign land use values unambiguously
# 2) remove coord and fasta lines that fall in water according to HYDE 3.2 data
#
############################################

rm(list=ls())
library(stringr)
library(dplyr)

# Chloe's function
cell_coords <- function(lat_seq, long_seq){
  ## Transform sequence coordinates
  # 5 arc minute = 1 unit
  new_lat <- lat_seq * 12
  new_long <- long_seq * 12
  
  # Bottom-left cell with bottom-left corner at c(0,0)
  new_lat <- new_lat + (90 * 12)
  new_long <- new_long  + (180 * 12)
  
  ## Does sequence fall on a cell boundary?
  if((new_lat %% 1 == 0) | (new_long %% 1 == 0)){
    border <- 1
  }else{
    border <- 0
  }
  
  ## If not, compute transformed cell coordinates
  if(border == 0){
    new_lat_cell <- (new_lat %/% 1) + (1/2)
    new_long_cell <- (new_long %/% 1) + (1/2)
  }else{
    lat_cell <- long_cell <- NA
  }
  
  ## Untransformed cell coordinates
  if(border == 0){
    # Bottom-left cell with bottom-left corner at c(-90,-180)
    new_lat_cell <- new_lat_cell - (90*12)
    new_long_cell <- new_long_cell - (180*12)
    
    # 5 arc minute = 1/12 unit
    lat_cell <- new_lat_cell/12
    long_cell <- new_long_cell/12
  }
  
  # Return result
  result <- c(lat_cell, long_cell, border)
  return(result)
}

#path containing all subfolders with .coords files
coords.dir <- ('/Users/vincentfugere/Desktop/Data/coord_files/')

#path containing all subfolders with FASTA files
fasta.dir <- ('/Users/vincentfugere/Desktop/Data/FASTA_files/')

#path to output folder for summary files of filtration process
summary.dir <- ('/Users/vincentfugere/Desktop/Data/Summary_deleted_sequences/')

#load HYDE 3.2 data
load('/Users/vincentfugere/Desktop/Data/Hyde_fileshyde32_2017_08.RData') # add directory
hyde32_2017_08 <- hyde32_2017_08[,1:3]
hyde32_2017_08$lat <- round(hyde32_2017_08$lat, 5)
hyde32_2017_08$long <- round(hyde32_2017_08$long, 5)

#name of subfolders (taxa)
folders <- list.files(coords.dir)

for(folder in folders){
  
  # working directories
  input.coords <- paste0(coords.dir,folder,'/ungrouped/')
  input.fasta <- paste0(fasta.dir,folder,'/')
  output.coords <- paste0(coords.dir,folder,'/ungrouped_filtered/')
  output.fasta <-paste0(fasta.dir,folder,'_filtered/')
  
  # grab .files (create file list of .coords)
  coordlist<-list.files(input.coords, pattern="*.coords", full.names=FALSE)
  namelist<-str_remove(coordlist, ".coords")
  
  seq_removed <- data.frame(filename=character(),
                            Nsequences_border=integer(),
                            Nsequences_water=integer(),
                            Nsequences_total=integer(),
                            entireFileRemoved = character(),
                            stringsAsFactors = F)
  
  # apply chloe's function and write new .coord and .fasta files
  for (n in namelist) {
    # Get coords and fasta files
    coordsfile<-paste0(n, ".coords")
    fastafile<-paste0(n, ".fasta")
    
    coords.df <- read.table(paste(input.coords, coordsfile, sep="/"), header = FALSE)
    fasta.df <- read.table(paste(input.fasta, fastafile, sep="/"), header = FALSE)
    
    # identify sequences that fall on cell border & get HYDE 3.2 cell coordinates
    coords.df[, c("cell_lat", "cell_long", "border")] <- t(mapply(FUN=cell_coords, lat_seq=coords.df[,1], long_seq=coords.df[,2]))
    
    # get row numbers of sequences that don't fall on cell border
    seqs_keep <- which(coords.df$border == 0)
    nseqs_removed_border <- nrow(coords.df) - length(seqs_keep)
    
    # only keep corresponding sequences in .coords and .fasta files
    coords.df <- coords.df[seqs_keep,]
    fasta.df <- fasta.df[seqs_keep,]
    
    # get row numbers of sequences that don't fall in water
    coords.df$cell_lat <- round(coords.df$cell_lat, 5)
    coords.df$cell_long <- round(coords.df$cell_long, 5)
    coords.df <- left_join(coords.df, hyde32_2017_08, by=c("cell_lat" = "lat", "cell_long" = "long"))
    seqs_keep <- which(!is.na(coords.df$conv_rangeland))
    nseqs_removed_water <- nrow(coords.df) - length(seqs_keep)
    
    # only keep corresponding sequences in .coords and .fasta files
    coords.df <- coords.df[seqs_keep,]
    fasta.df <- fasta.df[seqs_keep,]
    
    coords.df <- coords.df[,1:2]
    
    # writes new .coord and .fasta files, if non empty
    if (nrow(coords.df)>0){
      write.table(file=paste(output.coords, coordsfile, sep=""), x=coords.df, row.names = F, col.names = F)
      write.table(file=paste(output.fasta, fastafile, sep=""), x=fasta.df, row.names = F, col.names = F)
    }
    
    # track removed sequences
    seq_removed[which(namelist == n), 1] <- n
    seq_removed[which(namelist == n), 2] <- nseqs_removed_border
    seq_removed[which(namelist == n), 3] <- nseqs_removed_water
    seq_removed[which(namelist == n), 4] <- nseqs_removed_border + nseqs_removed_water
    seq_removed[which(namelist == n), 5] <- ifelse(nrow(coords.df)>0, "No", "Yes")
  }
  
  write.csv(file=paste0(summary.dir, folder,"_sequences_removed.csv"), x=seq_removed)
  
}
