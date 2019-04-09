#### Temporal variation in intra-specific neutral genetic diversity across anthromes
#### Gonzalez Lab project - McGill University - 2016-2019
#### Script by Chloé Debyser

#### HYDE 3.2 - Data Formatting
# This script takes the raw HYDE 3.2 text files and:
# 1) Aggregates land use and human density values at the 1, 2, and 4 degree spatial scales
# 2) Converts results to long format dataframes (with columns lat, long, as well as human density and land use variables)

#### Load Workspace ####
## Working directory
wd <- "C:/Users/Chloe/Documents/Études/2. Graduate - M.Sc/6. Lab - Meetings & Projects/Anthropocene Paper/Land Use Data/HYDE 3.2"
setwd(wd)

## Packages
library(raster) # For plotting and aggregating the maps
library(parallel) # To parallelize computations

## Memory Limit
memory.limit(size=10000) # Set to 10,000Mb

## Data
# Years of interest
time_points <- c(1980)

# General files
general.data <- c("garea", "im_reg", "iso", "landlake", "maxln", "sub_iso") # Names of the files
garea <- read.table(paste(wd, "/General_files/garea_cr.asc.txt", sep=""))
im_reg <- read.table(paste(wd, "/General_files/im_reg_cr.asc.txt", sep=""))
iso <- read.table(paste(wd, "/General_files/iso_cr.asc.txt", sep=""))
landlake <- read.table(paste(wd, "/General_files/landlake.asc.txt", sep=""))
maxln <- read.table(paste(wd, "/General_files/maxln_cr.asc.txt", sep=""))
sub_iso <- read.table(paste(wd, "/General_files/sub_iso_cr.asc.txt", sep=""))

# Land use
landuse.data <- c("conv_rangeland", "cropland", "grazing", "pasture", "rangeland", "tot_irri", "tot_rainfed", "tot_rice") # Names of the files selected
for (i in 1:length(time_points)){
  directory <- paste(wd, "/Land_use/", time_points[[i]], "AD_lu", sep="")
  for (j in 1:length(landuse.data)){
    assign(paste(landuse.data[[j]], time_points[[i]], sep=""), read.table(paste(directory, "/", landuse.data[[j]], time_points[[i]], "AD.asc", sep="")))
  }
}

# Population density
pop.data <- c("popc", "popd", "rurc", "uopp", "urbc") # Names of the files selected
for (i in 1:length(time_points)){
  directory <- paste(wd, "/Pop_density/", time_points[[i]], "AD_pop", sep="")
  for (j in 1:length(pop.data)){
    assign(paste(pop.data[[j]], time_points[[i]], sep=""), read.table(paste(directory, "/", pop.data[[j]], "_", time_points[[i]], "AD.asc", sep="")))
  }
}

#### Convert to Raster ####
# General files
for (i in 1:length(general.data)){
  matrix <- as.matrix(get(general.data[[i]]))
  matrix[matrix==-9999] <- NA
  raster <- raster(matrix)
  extent(raster) <- extent(-180, 180, -90, 90)
  assign(paste(general.data[[i]], "_08", sep=""), raster)
}

# Land use & Population density
hist.data <- c(landuse.data, pop.data)
for (i in 1:length(time_points)){
  for (j in 1:length(hist.data)){
    matrix <- as.matrix(get(paste(hist.data[[j]], time_points[[i]], sep="")))
    matrix[matrix==-9999] <- NA
    raster <- raster(matrix)
    extent(raster) <- extent(-180, 180, -90, 90)
    assign(paste(hist.data[[j]], time_points[[i]], "_08", sep=""), raster)
  }
}
rm(matrix)
rm(raster)

# Create list of all map datasets
all.data <- c(general.data, apply(expand.grid(landuse.data, time_points), 1, paste, collapse=""), apply(expand.grid(pop.data, time_points), 1, paste, collapse=""))

#### Plot the .08 Maps ####
for (i in 1:length(all.data)){
  plot(get(paste(all.data[[i]], "_08", sep="")), maxpixels=9331200)
  title(all.data[[i]])
}

#### Aggregate Maps at Lower Spatial Resolutions: 1, 2, and 4 ####
# List of scales
spatial.scales <- c("08", "1", "2", "4")
spatial_scales <- c(1/12, 1, 2, 4)

# Land mask
assign("maxln_1", aggregate(maxln_08, fact=12, fun=sum, na.rm=TRUE))
assign("maxln_2", aggregate(maxln_08, fact=24, fun=sum, na.rm=TRUE))
assign("maxln_4", aggregate(maxln_08, fact=48, fun=sum, na.rm=TRUE))

# Land use
for (i in 1:length(landuse.data)){
  for (j in 1:length(time_points)){
    assign(paste(landuse.data[[i]], time_points[[j]], "_1", sep=""), aggregate(get(paste(landuse.data[[i]], time_points[[j]], "_08", sep="")), fact=12, fun=sum, na.rm=TRUE))
    assign(paste(landuse.data[[i]], time_points[[j]], "_2", sep=""), aggregate(get(paste(landuse.data[[i]], time_points[[j]], "_08", sep="")), fact=24, fun=sum, na.rm=TRUE))
    assign(paste(landuse.data[[i]], time_points[[j]], "_4", sep=""), aggregate(get(paste(landuse.data[[i]], time_points[[j]], "_08", sep="")), fact=48, fun=sum, na.rm=TRUE))
  }
}

# Population density
abs.pop.data <- c("popc", "rurc", "urbc", "uopp")
for (i in 1:length(abs.pop.data)){
  for (j in 1:length(time_points)){
    assign(paste(abs.pop.data[[i]], time_points[[j]], "_1", sep=""), aggregate(get(paste(abs.pop.data[[i]], time_points[[j]], "_08", sep="")), fact=12, fun=sum, na.rm=TRUE))
    assign(paste(abs.pop.data[[i]], time_points[[j]], "_2", sep=""), aggregate(get(paste(abs.pop.data[[i]], time_points[[j]], "_08", sep="")), fact=24, fun=sum, na.rm=TRUE))
    assign(paste(abs.pop.data[[i]], time_points[[j]], "_4", sep=""), aggregate(get(paste(abs.pop.data[[i]], time_points[[j]], "_08", sep="")), fact=48, fun=sum, na.rm=TRUE))
  }
}

#### Divide Absolute Values by Land Area, in Each Cell ####
# Land use
for (k in 1:length(spatial.scales)){
  for (j in 1:length(time_points)){
    for (i in 1:length(landuse.data)){
      dataset <- get(paste(landuse.data[[i]], time_points[[j]], "_", spatial.scales[[k]], sep=""))
      landmask <- get(paste("maxln_", spatial.scales[[k]], sep=""))
      ratio <- dataset / landmask         # Portion of land area occupied by the land use type of interest (km2 landuse / km2 land)
      assign(paste(landuse.data[[i]], time_points[[j]], "_", spatial.scales[[k]], "_ratio", sep=""), ratio)
    }
  }
}

# Population density
for (k in 1:length(spatial.scales)){
  for (j in 1:length(time_points)){
    for (i in 1:length(abs.pop.data)){
      dataset <- get(paste(abs.pop.data[[i]], time_points[[j]], "_", spatial.scales[[k]], sep=""))
      landmask <- get(paste("maxln_", spatial.scales[[k]], sep=""))
      ratio <- dataset / landmask         # Population density (inh/km2) OR Portion of land area occupied by the land use type of interest (km2 landuse / km2 land)
      assign(paste(abs.pop.data[[i]], time_points[[j]], "_", spatial.scales[[k]], "_ratio", sep=""), ratio)
    }
  }
}
rm(dataset)
rm(landmask)
rm(ratio)

#### Convert Final Maps to Long Format ####
## Create a template results database for each spatial scale
for (i in 1:length(spatial.scales)){
  landmask <- get(paste("maxln_", spatial.scales[[i]], sep=""))
  coord <- coordinates(landmask)
  coord <- as.data.frame(coord)
  colnames(coord) <- c("long", "lat")
  assign(paste("hyde32_", spatial.scales[[i]], sep=""), coord)
}
rm(landmask)
rm(coord)

## Function to extract raster value for each cell coordinate
raster.value <- function(x, ratio, results){
  lat <- results[x, "lat"]
  long <- results[x, "long"]
  xy <- cbind(long, lat)
  value <- extract(ratio, xy)
  return(value)
}

## Set up a parallelized computing cluster
# Calculate the number of cores
no_cores <- detectCores() - 1

# Initiate cluster
cl <- makeCluster(no_cores)

# Export functions to cluster
clusterExport(cl, "raster.value")
clusterExport(cl, "extract")

## Land use
for (k in rev(1:length(spatial.scales))){
  for (j in 1:length(time_points)){
    assign("results", get(paste("hyde32_", spatial.scales[[k]], sep="")))
    for (i in 1:length(landuse.data)){
      ratio <- get(paste(landuse.data[[i]], time_points[[j]], "_", spatial.scales[[k]], "_ratio", sep=""))
      results[,landuse.data[[i]]] <- unlist(parSapply(cl, X=1:nrow(results), FUN=raster.value, ratio=ratio, results=results))
      print(paste("Data=", landuse.data[[i]], "Sys.time=", Sys.time()))
    }
    assign(paste("hyde32_", time_points[[j]], "_", spatial.scales[[k]], sep=""), results)
    print(paste("Scale=", spatial.scales[[k]], "Year=", time_points[[j]], "Sys.time=", Sys.time()))
  }
}

## Population density
for (k in rev(1:length(spatial.scales))){
  for (j in 1:length(time_points)){
    assign("results", get(paste("hyde32_", time_points[[j]], "_", spatial.scales[[k]], sep="")))
    for (i in 1:length(abs.pop.data)){
      ratio <- get(paste(abs.pop.data[[i]], time_points[[j]], "_", spatial.scales[[k]], "_ratio", sep=""))
      results[,abs.pop.data[[i]]] <- unlist(parSapply(cl, X=1:nrow(results), FUN=raster.value, ratio=ratio, results=results))
      print(paste("Data=", abs.pop.data[[i]], "Sys.time=", Sys.time()))
    }
    assign(paste("hyde32_", time_points[[j]], "_", spatial.scales[[k]], sep=""), results)
    print(paste("Scale=", spatial.scales[[k]], "Year=", time_points[[j]], "Sys.time=", Sys.time()))
  }
}
rm(ratio)
rm(results)

## Stop parallelized cluster
stopCluster(cl)

## Remove template results databases
rm(hyde32_08)
rm(hyde32_1)
rm(hyde32_2)
rm(hyde32_4)

#### Save Results ####
# Choose column order
col.order <- c("lat", "long", landuse.data, "uopp", "popc", "rurc", "urbc")

# Save data
for (i in rev(1:length(spatial.scales))){
  for (j in 1:length(time_points)){
    # Re-order columns
    final_dataset <- get(paste("hyde32_", time_points[[j]], "_", spatial.scales[[i]], sep=""))
    final_dataset <- final_dataset[,col.order]
    
    # Change column names
    names(final_dataset)[names(final_dataset) == "uopp"] <- "urban"
    names(final_dataset)[names(final_dataset) == "popc"] <- "pop_tot"
    names(final_dataset)[names(final_dataset) == "urbc"] <- "pop_urb"
    names(final_dataset)[names(final_dataset) == "rurc"] <- "pop_rur"
    
    # Save
    assign(paste("hyde32_", time_points[[j]], "_", spatial.scales[[i]], sep=""), final_dataset)
    save(list = paste("hyde32_", time_points[[j]], "_", spatial.scales[[i]], sep=""), file=paste("hyde32_", time_points[[j]], "_", spatial.scales[[i]], ".RData", sep=""))
  }
}
rm(final_dataset)

#### Verify Results ####
## Mutually exclusive?
# YES: "cropland", "grazing", "urban"
for (i in rev(1:length(spatial.scales))){
  for (j in 1:length(time_points)){
    final_dataset <- get(paste("hyde32_", time_points[[j]], "_", spatial.scales[[i]], sep=""))
    sum <- final_dataset$cropland + final_dataset$grazing + final_dataset$urban
    print(paste("Scale=", spatial.scales[[i]], "Year=", time_points[[j]], sep=" "))
    if(max(sum, na.rm=T)<(1+1/100000)){
      print(paste("CORRECT", max(sum, na.rm=T), sep=" "))
    } else {
      print(paste("WRONG", max(sum, na.rm=T), sep=" "))
    }
  }
}

## Subsets are absolute complements?
# YES: "tot_irri" & "tot_rainfed" IN  "cropland" -- ONE CELL HAS A MISTAKE FOR 1980! (long=-56.625, lat=5.791667)
for (i in rev(1:length(spatial.scales))){
  for (j in 1:length(time_points)){
    final_dataset <- get(paste("hyde32_", time_points[[j]], "_", spatial.scales[[i]], sep=""))
    dif <- abs(final_dataset$cropland - (final_dataset$tot_irri + final_dataset$tot_rainfed))
    print(paste("Scale=", spatial.scales[[i]], "Year=", time_points[[j]], sep=" "))
    if(max(dif, na.rm=T)<(1/100000)){
      print(paste("CORRECT", max(dif, na.rm=T), sep=" "))
    } else {
      print(paste("WRONG", max(dif, na.rm=T), sep=" "))
    }
  }
}

# YES: "pop_rur" & "pop_urb" IN "pop_tot"
for (i in rev(1:length(spatial.scales))){
  for (j in 1:length(time_points)){
    final_dataset <- get(paste("hyde32_", time_points[[j]], "_", spatial.scales[[i]], sep=""))
    dif <- abs(final_dataset$pop_tot - (final_dataset$pop_rur + final_dataset$pop_urb))
    print(paste("Scale=", spatial.scales[[i]], "Year=", time_points[[j]], sep=" "))
    if(max(dif, na.rm=T)<(1/100000)){
      print(paste("CORRECT", max(dif, na.rm=T), sep=" "))
    } else {
      print(paste("WRONG", max(dif, na.rm=T), sep=" "))
    }
  }
}

# YES: "pasture" & "rangeland" & "conv_rangeland" IN "grazing"
for (i in rev(1:length(spatial.scales))){
  for (j in 1:length(time_points)){
    final_dataset <- get(paste("hyde32_", time_points[[j]], "_", spatial.scales[[i]], sep=""))
    dif <- abs(final_dataset$grazing - (final_dataset$pasture + final_dataset$rangeland + final_dataset$conv_rangeland))
    print(paste("Scale=", spatial.scales[[i]], "Year=", time_points[[j]], sep=" "))
    if(max(dif, na.rm=T)<(1/100000)){
      print(paste("CORRECT", max(dif, na.rm=T), sep=" "))
    } else {
      print(paste("WRONG", max(dif, na.rm=T), sep=" "))
    }
  }
}

## Simple subsets?
# YES: "tot_rice" IN "cropland"
for (i in rev(1:length(spatial.scales))){
  for (j in 1:length(time_points)){
    final_dataset <- get(paste("hyde32_", time_points[[j]], "_", spatial.scales[[i]], sep=""))
    dif <- final_dataset$cropland - final_dataset$tot_rice
    print(paste("Scale=", spatial.scales[[i]], "Year=", time_points[[j]], sep=" "))
    if(min(dif, na.rm=T)>(-1/100000)){
      print(paste("CORRECT", min(dif, na.rm=T), sep=" "))
    } else {
      print(paste("WRONG", min(dif, na.rm=T), sep=" "))
    }
  }
}
rm(final_dataset)

#### Extra Code ####
# For getting raster cell numbers which meet a certain conditions
Which(conv_rangeland1980_08 > 5, cells = TRUE) 

# For extracting a raster value based on cell numbers
extract(conv_rangeland1980_08, c(1297884, 1297885))

# For plotting in raster at full resolution and with colour breaks
break_pts <- c(0, 0.0000000001, 1/80, 5/80, 10/80, 40/80, 50/80, 60/80, 70/80, 80/80)
colours <- c("papayawhip", "deepskyblue", "mediumaquamarine", "green", "yellow3", "yellow", "orange", "tan3", "red")
plot(cropland1980_08_ratio, breaks=break_pts, col=colours, maxpixels=9331200)

# For plotting without converting to raster
land_msk <- garea
dim(land_msk)

land_msk <- land_msk[rev(rownames(land_msk)),]
land_msk_matrix <- data.matrix(land_msk)

heatmap <- heatmap(land_msk_matrix, Rowv=NA, Colv=NA)

# Differences between 1980 and 2017 datasets
for (j in spatial.scales){
  list <- colnames(get(paste("hyde32_1980_", j, sep="")))[3:length(colnames(get(paste("hyde32_1980_", j, sep=""))))]
  for (i in 1:length(list)){
    data_1980 <- get(paste("hyde32_1980_", j, sep=""))
    data_2017 <- get(paste("hyde32_2017_", j, sep=""))
    if(sum(!is.na(data_1980[,list[[i]]]))==sum(data_1980[,list[[i]]]==data_2017[,list[[i]]], na.rm=T)){
      print(paste(j, " ", list[[i]], ": SAME", sep=""))
    }else{
      print(paste(j, " ", list[[i]], ": DIFFERENT", sep=""))
    }
  }
}

# Differences between old and new datasets
list <- colnames(old_2017)[3:length(colnames(old_2017))]
for (j in 1:length(time_points)){
  old <- get(paste("old_", time_points[[j]], sep=""))
  new <- get(paste("new_", time_points[[j]], sep=""))
  for (i in 1:length(list)){
    if(sum(!is.na(old[,list[[i]]]))==sum(old[,list[[i]]]==new[,list[[i]]], na.rm=T)){
      print(paste(time_points[[j]], " ", list[[i]], ": SAME", sep=""))
    }else{
      print(paste(time_points[[j]], " ", list[[i]], ": DIFFERENT", sep=""))
    }
  }
}

