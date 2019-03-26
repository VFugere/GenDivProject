##
# Script to plot coordinates on map, determine if on cell border
# remove coord and fasta lines that land on grid cell border because 
# can't assign land use values unambiguously
#
# Workflow:
#  iterate through each line of coordsdf, first column lat second column long
#  pass these arguments to chloe's func  function(coordsdf[0,0], coordsdf[0,1])
#  if (chloe(lat,long))
#  then (do this when 1) - on a cell boundary
#  else (do this when 0) - not on a cell boundary
#  write coords and fasta line to new file
############################################

rm(list=ls())
library(stringr)

# Chloe's function
cell_border <- function(lat_seq, long_seq){
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
  
  ## Return result
  return(border)
}

#path containing all subfolders with .coords files
coords.dir <- ('/Users/vincentfugere/Desktop/Data/coord_files/')

#path containing all subfolders with FASTA files
fasta.dir <- ('/Users/vincentfugere/Desktop/Data/FASTA_files/')

#path to output folder for summary files of filtration process
summary.dir <- ('/Users/vincentfugere/Desktop/Data/Summary_deleted_sequences/')

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
                            Nsequences=integer(),
                            entireFileRemoved = character(),
                            stringsAsFactors = F)
  
  # apply chloe's function and write new .coord and .fasta files
  for (n in namelist) {
    # Get coords and fasta files
    coordsfile<-paste0(n, ".coords")
    fastafile<-paste0(n, ".fasta")
    
    coords.df <- read.table(paste(input.coords, coordsfile, sep="/"), header = FALSE)
    fasta.df <- read.table(paste(input.fasta, fastafile, sep="/"), header = FALSE)
    
    # identify sequences that fall on cell border
    coords.df$border <- mapply(FUN=cell_border, lat_seq=coords.df[,1], long_seq=coords.df[,2])
    
    # gets row numbers of sequences that don't fall on cell border
    seqs_keep <- which(coords.df$border == 0)
    nseqs_removed <- nrow(coords.df) - length(seqs_keep)
    
    # only keep corresponding sequences in .coords and .fasta files
    coords.df <- coords.df[seqs_keep,]
    coords.df <- coords.df[,1:2]
    
    fasta.df <- fasta.df[seqs_keep,]
    
    # writes new .coord and .fasta files, if non empty
    if (nrow(coords.df)>0){
      write.table(file=paste(output.coords, coordsfile, sep=""), x=coords.df, row.names = F, col.names = F)
      write.table(file=paste(output.fasta, fastafile, sep=""), x=fasta.df, row.names = F, col.names = F)
    }
    
    # track removed sequences
    seq_removed[which(namelist == n), 1] <- n
    seq_removed[which(namelist == n), 2] <- nseqs_removed
    seq_removed[which(namelist == n), 3] <- ifelse(nrow(coords.df)>0, "No", "Yes")
  }
  
  write.csv(file=paste0(summary.dir, folder,"_sequences_removed.csv"), x=seq_removed)
  
}
