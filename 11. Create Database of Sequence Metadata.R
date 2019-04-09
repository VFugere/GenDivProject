#### Temporal variation in intra-specific neutral genetic diversity across anthromes
#### Gonzalez Lab project - McGill University - 2016-2019
#### Script by Chloé Debyser

#### CREATE A DATABASE OF SEQUENCE METADATA
#### This script takes all clustered .coord files and creates a unified sequence metadata table,
#### including 1) sequence lat and long, 2) HYDE 3.2 cell the sequence falls in, 3) corresponding
#### HYDE 3.2 data.

#### Load Workspace ####
# Working directories
hyde32 <- "C:/Users/Chloe/Documents/Études/2. Graduate - M.Sc/6. Lab - Meetings & Projects/Anthropocene Paper/Land Use Data/HYDE 3.2/"
dataDir <- "C:/Users/Chloe/Documents/Études/2. Graduate - M.Sc/6. Lab - Meetings & Projects/Anthropocene Paper/Data/"
coords <- paste0(dataDir, "coord_files/")

# Packages
library(plyr)
library(dplyr)

# Data
taxonomy <- get(load(paste0(dataDir, "metadata/species_taxonomy.RData")))

#### Database of Sequence Metadata ####
taxa <- list.files(file.path(coords))

for (taxon in taxa){
  # Get scales of analysis
  scales <- list.files(file.path(paste0(coords, taxon)), pattern="grouped_mindist")
  scales <- as.integer(unique(substr(scales, 16, nchar(scales))))
  
  # Create empty database of sequences
  seq <- data.frame(species = character(),
                    year = integer(),
                    FASTArow = integer(),
                    lat = numeric(),
                    long = numeric(),
                    stringsAsFactors = F)
  
  for (scale in scales){
    files <- list.files(file.path(paste0(coords, taxon, "/grouped_mindist", scale, "/")))
    
    # First scale
    if (match(scale, scales) == 1){
      seq[,paste0("pop", scale)] <- character(0)
      
      for (file in files){
        # Get data
        data <- read.table(file=paste0(coords, taxon, "/grouped_mindist", scale, "/", file))
        
        # Cleanup data
        colnames(data) <- c(paste0("pop", scale), "lat", "long") 
        data$FASTArow <- as.integer(rownames(data))
        data$species <- gsub("_", " ", substr(file, 1, nchar(file)-5))
        data$year <- as.integer(substr(file, nchar(file)-3, nchar(file)))
        data <- data[,colnames(seq)]  
        
        # Add to sequence database
        seq <- rbind(seq, data)
      }
      
      # Number of sequences for that taxon
      nseq <- nrow(seq)
    
    # Other scales
    }else{
      new_seq <- data.frame(species = character(),
                            year = integer(),
                            FASTArow = integer(),
                            lat = numeric(),
                            long = numeric(),
                            stringsAsFactors = F)
      new_seq[,paste0("pop", scale)] <- character(0)
      
      for (file in files){
        # Get data
        data <- read.table(file=paste0(coords, taxon, "/grouped_mindist", scale, "/", file))
        
        # Cleanup data
        colnames(data) <- c(paste0("pop", scale), "lat", "long") 
        data$FASTArow <- as.integer(rownames(data))
        data$species <- gsub("_", " ", substr(file, 1, nchar(file)-5))
        data$year <- substr(file, nchar(file)-3, nchar(file))
        
        # Add to sequence database
        new_seq <- rbind(new_seq, data)
      }
      
      # Join with sequence database for other scales
      seq <- join(seq, new_seq, by=c("species", "year", "FASTArow", "lat", "long"), type="full")
    }
    
    print(paste(taxon, scale, Sys.time()))
  }
  
  # Check that the number of sequences is the same across scales (i.e. that the join operation worked properly)
  if(!nrow(seq) == nseq){print("ERROR: join operation did not work as expected")}
  
  # Save sequence metadata for the focal taxon
  seq$class <- taxon
  seq <- seq[,c(ncol(seq), 1:(ncol(seq)-1))]
  assign(paste0("seq_", taxon), seq)
}
rm(seq, new_seq, data, nseq)

# Consolidate sequence metadata for all taxa into a single database
seq_list <- as.list(paste0("seq_", taxa))
seq_list <- lapply(X=seq_list, FUN=get)
seq <- bind_rows(seq_list)
rm(seq_list)

# Add unique sequence ID
seq$seqID <- 1:nrow(seq)
seq <- seq[,c(ncol(seq), 1:(ncol(seq)-1))]

#### Add Land Use and Human Density Data ####
# Function to assign cell coordinates to sequence coordinates
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

# Apply
seq[, c("cell_lat", "cell_long", "border")] <- t(mapply(FUN=cell_coords, lat_seq=seq$lat, long_seq=seq$long))
seq$border <- sapply(seq$border, function(x) ifelse(x==1, "Yes", "No"))

# Add LU and HD data to sequence metadata, for each year
years <- as.integer(substr(list.files(file.path(hyde32), pattern="_08.RData"), start=8, stop=11))
seq_list <- list()

for (year in years){
  # Load HD & LU data
  load(paste0(hyde32, "hyde32_", year, "_08.RData"))
  assign("LUHD", get(paste0("hyde32_", year, "_08")))
  rm(list=paste0("hyde32_", year, "_08"))
  
  # Select sequence collection years corresponding to the focal land use map
  if (year < 2000){
    years_included <- year:(year+9)
  }else{
    years_included <- year
  }
  
  # Only keep sequences from those years
  seq_sub <- seq[which(seq$year %in% years_included),]
  
  # Round lat and long values in same way
  seq_sub$cell_lat <- round(seq_sub$cell_lat, 5)
  seq_sub$cell_long <- round(seq_sub$cell_long, 5)
  LUHD$lat <- round(LUHD$lat, 5)
  LUHD$long <- round(LUHD$long, 5)
  
  # Join
  seq_sub <- left_join(seq_sub, LUHD, by=c("cell_lat" = "lat", "cell_long" = "long"))
  
  # Save
  seq_list[[as.name(year)]] <- seq_sub
  
  # Remove HD & LU data
  rm(LUHD)
  print(paste(year, Sys.time()))
}

# Add all years to the same database
seq <- bind_rows(seq_list)

# Remove "border" column
seq <- seq[order(seq$seqID),]
rownames(seq) <- 1:nrow(seq)
seq <- subset(seq, select=-c(border))
seq <- seq[,c(1:7, 13:14, 8:12, 15:26)]

#### Add Taxonomic Information ####
# All sequences
seq <- left_join(seq, taxonomy[,c(1,3:5)], by="species")
seq <- seq[,c(1:2,27:29,3:26)]

# Two species that are missing from the taxonomy dataset
seq[which(seq$species == "Mydaea winnemana"), "order"] <- "Diptera"
seq[which(seq$species == "Mydaea winnemana"), "family"] <- "Muscidae"
seq[which(seq$species == "Monomorium bicolor"), "order"] <- "	Hymenoptera"
seq[which(seq$species == "Monomorium bicolor"), "family"] <- "Formicidae"

#### Save Database ####
save(seq, file = paste0(dataDir, "metadata/sequence_metadata.RData"))
