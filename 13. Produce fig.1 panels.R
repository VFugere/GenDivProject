#### Temporal variation in intra-specific neutral genetic diversity across anthromes
#### Gonzalez Lab project - McGill University - 2016-2019
#### Script by Chloé Debyser

#### PRODUCE FIG.1 PANELS
#### This script produces the 5 panels of Fig.1

#### Load workspace ####
## Working directory
setwd("C:/Users/Chloe/Dropbox/GenDivProject - R Codes")
dataDir <- "Data/Metadata/"
finalDir <- "Figure1/"

## Packages
library(scales) # For % scale in histograms
library(ggthemes) # To change theme of ggplots
library(latticeExtra)
library(raster) # For creation of background maps
library(ncdf4)
library(tidyr) # For splitting a column into two & reshaping data from long to wide format
library(plyr) # For summarizing rows with same Id
library(viridis); library(readxl); library(stringr) # For bivariate color scheme
library(ggplot2) # For plotting the legend
library(colorspace) # For creating a bivariate colour scheme
library(rworldmap) # For world outline
library(countrycode) # For getting continents
library(dplyr) # For merging databases
library(forcats)
library(ggrepel)
library(rphylopic) # For symbols
windowsFonts(Helvetica=windowsFont("TT Helvetica"))

## Memory
memory.limit(size=3700000)

## Metadata
load(file = paste0(dataDir, "sequence_metadata.RData"))

## Count data
counts_dated <- get(load(file = "Data_Old/CountData_Final_Dated.RData"))
rm(fig1data)

## Settings
tax_groups <- unique(counts_dated$tax)
spatial_scales <- c(1/12, 1, 2, 4)
spatial.scales <- c("08", "1", "2", "4")
years <- 1980:2016

#### Format - Sequence count data ####
## Counts for each taxonomic group
for (i in 1:length(tax_groups)){
  for (j in 1:length(spatial.scales)){
    dataset <- counts_dated[which(counts_dated$tax == tax_groups[[i]] & counts_dated$scale == spatial.scales[[j]]),]
    dataset <- dataset[,c("cell", "year", "sequences")]
    
    dataset <- separate(dataset, col=cell, into=c("lat", "long"), sep="_")
    dataset[, c("lat", "long")] <- apply(dataset[, c("lat", "long")], 2, function(x) as.double(x))
    
    dataset <- dataset[order(dataset$lat, dataset$long, dataset$year),]
    dataset <- as.data.frame(dataset)
    rownames(dataset) <- 1:nrow(dataset)
    
    assign(paste(tax_groups[[i]], spatial.scales[[j]], sep="_"), dataset)
  }
}
rm(dataset)

## Pooled counts
for (i in 1:length(spatial.scales)){
  dataset <- counts_dated[which(counts_dated$scale == spatial.scales[[i]]),]
  dataset <- dataset[,c("cell", "year", "sequences")]
  
  dataset <- ddply(dataset, .(cell, year), summarize, sequences = sum(sequences))
  
  dataset <- separate(dataset, col=cell, into=c("lat", "long"), sep="_")
  dataset[, c("lat", "long")] <- apply(dataset[, c("lat", "long")], 2, function(x) as.double(x))
  
  dataset <- dataset[order(dataset$lat, dataset$long, dataset$year),]
  dataset <- as.data.frame(dataset)
  rownames(dataset) <- 1:nrow(dataset)
  
  assign(paste("pooled", spatial.scales[[i]], sep="_"), dataset)
  
}
rm(dataset)
tax_groups <- c(tax_groups, "pooled")

## Compile data (one matrix per scale)
for (i in 1:length(spatial.scales)){
  birds <- get(paste("birds", spatial.scales[[i]], sep="_"))
  fish <- get(paste("fish", spatial.scales[[i]], sep="_"))
  insects <- get(paste("insects", spatial.scales[[i]], sep="_"))
  mammals <- get(paste("mammals", spatial.scales[[i]], sep="_"))
  pooled <- get(paste("pooled", spatial.scales[[i]], sep="_"))
  
  final <- merge(birds, fish, by=c("lat", "long", "year"), all=T)
  colnames(final)[colnames(final)=="sequences.x"] <- "birds"
  colnames(final)[colnames(final)=="sequences.y"] <- "fish"
  
  final <- merge(final, insects, by=c("lat", "long", "year"), all=T)
  colnames(final)[colnames(final)=="sequences"] <- "insects"
  
  final <- merge(final, mammals, by=c("lat", "long", "year"), all=T)
  colnames(final)[colnames(final)=="sequences"] <- "mammals"
  
  final <- merge(final, pooled, by=c("lat", "long", "year"), all=T)
  colnames(final)[colnames(final)=="sequences"] <- "total"
  
  final[which(is.na(final$birds)),"birds"] <- 0
  final[which(is.na(final$fish)),"fish"] <- 0
  final[which(is.na(final$insects)),"insects"] <- 0
  final[which(is.na(final$mammals)),"mammals"] <- 0
  
  assign(paste("counts", spatial.scales[[i]], sep="_"), final)
}
rm(final, birds, fish, insects, mammals, pooled)
rm(birds_08, birds_1, birds_2, birds_4, fish_08, fish_1, fish_2, fish_4, insects_08, insects_1, insects_2, insects_4, mammals_08, mammals_1, mammals_2, mammals_4, pooled_08, pooled_1, pooled_2, pooled_4)

## Save
save(counts_08, file="Counts_08.RData")
save(counts_1, file="Counts_1.RData")
save(counts_2, file="Counts_2.RData")
save(counts_4, file="Counts_4.RData")

#### Add - Unsampled cells ####
## Load formatted count data
load("Counts_08.RData")
load("Counts_1.RData")
load("Counts_2.RData")
load("Counts_4.RData")

## Create a template for each spatial scale
for (i in rev(1:length(spatial.scales))){
  map <- get(load(paste(wd, "/Land Use Data/HYDE 3.2/hyde32_2016_", spatial.scales[[i]], ".RData", sep="")))
  map <- map[,c(1:3)]
  map <- map[which(!is.na(map$conv_rangeland)),]
  map <- map[,c(1:2)]
  
  template <- as.data.frame(matrix(nrow=nrow(map)*length(years), ncol=3))
  colnames(template) <- c("lat", "long", "year")
  template$lat <- rep(map$lat, times=length(years))
  template$long <- rep(map$long, times=length(years))
  template$year <- rep(years, each=nrow(map))
  assign(paste("template", spatial.scales[[i]], sep="_"), template)
  print(paste("Spatial scale =", spatial.scales[[i]], Sys.time()))
}
rm(template, map)

## Merge template and count data
for (i in 1:length(spatial.scales)){
  print(paste("Spatial scale =", spatial.scales[[i]], Sys.time()))
  
  counts <- get(paste("counts", spatial.scales[[i]], sep="_"))
  template <- get(paste("template", spatial.scales[[i]], sep="_"))
  
  final <- left_join(template, counts, by=c("lat", "long", "year"))
  final <- final[order(final$lat, final$long, final$year),]
  rownames(final) <- 1:nrow(final)
  
  final[which(is.na(final$birds)),"birds"] <- 0
  final[which(is.na(final$fish)),"fish"] <- 0
  final[which(is.na(final$insects)),"insects"] <- 0
  final[which(is.na(final$mammals)),"mammals"] <- 0
  final[which(is.na(final$total)),"total"] <- 0
  
  assign(paste("counts", spatial.scales[[i]], sep="_"), final)
}
rm(counts, template, final)

## Save
save(counts_08, file="Counts_08.RData")
save(counts_1, file="Counts_1.RData")
save(counts_2, file="Counts_2.RData")
save(counts_4, file="Counts_4.RData")

#### Add - Land use data ####
## Load formatted count data
load("Counts_08.RData")
load("Counts_1.RData")
load("Counts_2.RData")
load("Counts_4.RData")

## Add land use data
for (i in rev(1:length(spatial.scales))){
  print(paste("Spatial scale =", spatial.scales[[i]], Sys.time()))
  
  counts <- get(paste("counts", spatial.scales[[i]], sep="_"))
  counts[,"year_proxy"] <- unlist(lapply(counts$year, function(x) if(x<1990){y<-1980}else if(x<2000){y<-1990}else{y<-x}))
  
  counts_final <- as.data.frame(matrix(nrow=0, ncol=11))
  colnames(counts_final) <- c(colnames(counts), "humans", "landuse")
  
  years <- sort(unique(counts$year_proxy))
  
  for (year in years[18:19]){
    print(paste("Year = ", year, Sys.time()))
    data <- counts[which(counts$year_proxy == year),]
    
    map <- get(load(file = paste("Land Use Data/HYDE 3.2/hyde32_", paste(year, spatial.scales[[i]], sep="_"), ".RData", sep="")))
    rm(list=paste("hyde32", year, spatial.scales[[i]], sep="_"))
    
    data <- left_join(data, map, by=c("lat", "long"))
    colnames(data)[colnames(data)=="pop_tot"] <- "humans"
    data[,"landuse"] <- data[,"cropland"] + data[,"pasture"] + data[,"urban"] + data[,"conv_rangeland"]
    data <- data[,colnames(counts_final)]
    
    counts_final <- bind_rows(counts_final, data)
  }
  counts_final <- counts_final[order(counts_final$lat, counts_final$long, counts_final$year), c("lat", "long", "year", "birds", "fish", "insects", "mammals", "total", "humans", "landuse")]
  rownames(counts_final) <- 1:nrow(counts_final)
  assign(paste("counts", spatial.scales[[i]], sep="_"), counts_final)
} 
rm(counts, counts_final, i, year, years, map)

## Save
save(counts_08, file="Counts_08.RData")
save(counts_1, file="Counts_1.RData")
save(counts_2, file="Counts_2.RData")
save(counts_4, file="Counts_4.RData")

#### Update count_08 using seq_metadata ####
## Load formatted count data
load(paste0(dataDir, "Counts_08.RData"))

# Standardize lat and long columns
seq$cell_lat <- round(seq$cell_lat, digits=4)
seq$cell_long <- round(seq$cell_long, digits=4)
counts_08$lat <- round(counts_08$lat, digits=4)
counts_08$long <- round(counts_08$long, digits=4)

# Summarize seq by gridcell, class, and year
seq_summarized <- seq %>%
  group_by(cell_lat, cell_long, class, year) %>%
  summarize(count = n())

# Long to wide
seq_summarized <- spread(seq_summarized, class, count)
seq_summarized[is.na(seq_summarized)] <- 0

# Add total column
seq_summarized$total <- rowSums(seq_summarized[, c("birds", "fish", "insects", "mammals")])
seq_summarized$year <- as.numeric(seq_summarized$year)
seq_summarized <- as.data.frame(seq_summarized)

# Join
counts_08 <- left_join(counts_08[, c("lat", "long", "year", "humans", "landuse")], seq_summarized, by=c("lat" = "cell_lat", "long" = "cell_long", "year"))
counts_08[is.na(counts_08)] <- 0
counts_08 <- counts_08[,c(1:3,6:10,4:5)]

# Save
save(counts_08, file=paste0(dataDir, "Counts_08_updated.RData"))
save(seq_summarized, file=paste0(dataDir, "Seq_summarized.RData"))

#### Panel A - 3D map ####
## Load final dataset
load(file=paste0(dataDir, "Seq_summarized.RData"))

## Summarize by cell_lat and cell_long
seq_summarized_cell <- seq_summarized %>%
  group_by(cell_lat, cell_long) %>%
  summarize(birds = sum(birds), fish = sum(fish), insects = sum(insects), mammals = sum(mammals), total = sum(total))

## Format World map
map_poly <- getMap(resolution = "coarse")

      # Extract the attribute table (map_data) and the polygon coordinates (polygons)
map_data <- map_poly@data
map_data <- map_data[[4]]
polygons <- map_poly@polygons

      # Format the map as a linear table with coordinates & attribute data
map_final <- data.frame(long = double(),
                        lat = double(),
                        polygon_ID = integer(),
                        country = character(),
                        stringsAsFactors = F)

p <- 0
for (i in 1:nrow(map_poly@data)){
  ext.1 <- polygons[[i]]
  ext.2 <- ext.1@Polygons
  
  for(j in 1:length(ext.2)){
    ext.3 <- ext.2[[j]]
    ext.4 <- ext.3@coords
    
    n <- nrow(map_final)
    p <- p+1
    
    map_final[(n+1):(n+nrow(ext.4)), "long"] <- ext.4[,1]
    map_final[(n+1):(n+nrow(ext.4)), "lat"] <- ext.4[,2]
    map_final[(n+1):(n+nrow(ext.4)), "polygon_ID"] <- p
    map_final[(n+1):(n+nrow(ext.4)), "country"] <- as.character(map_data[i])
  }
}
map_final <- as.data.frame(map_final)

      # Indicate polygon separation with blank lines
map_final_sep <- map_final
map_final_sep[(nrow(map_final)+(1:max(map_final$polygon_ID))),] <- NA
map_final_sep[,"row"] <- 0

for(i in unique(map_final$polygon_ID)){
  data <- map_final[map_final$polygon_ID == i,]
  map_final_sep[which(map_final_sep$polygon_ID == i), "row"] <- max(map_final_sep$row) + 1:nrow(data)
  map_final_sep[which(is.na(map_final_sep$long) & map_final_sep$row==0), "row"][1] <- max(map_final_sep$row) + 1
}
map_final_sep <- map_final_sep[order(map_final_sep$row),]
rownames(map_final_sep) <- map_final_sep[,"row"]
map_final_sep <- map_final_sep[1:(nrow(map_final_sep)-1),]

## Plot it
      # Settings
bar_size <- 1/24
colour <- "lightsteelblue2"

      # Set up the plot panel
panel.3dmap <- function(..., rot.mat, distance, xlim, ylim, zlim, xlim.scaled, ylim.scaled, zlim.scaled) {
  scaled.val <- function(x, original, scaled) {
    scaled[1] + (x - original[1]) * diff(scaled)/diff(original)
  }
  m <- ltransform3dto3d(rbind(scaled.val(map_final_sep$long, xlim, xlim.scaled), scaled.val(map_final_sep$lat, ylim, ylim.scaled), zlim.scaled[1]), rot.mat, distance)
  panel.polygon(m[1,], m[2,], col = colour, border = colour, lwd=0.01)
}
    
      # Produce the plot - Rectangular bars
plot <- cloud(total ~ cell_long + cell_lat, seq_summarized_cell, panel.3d.cloud = function(...) {
  panel.3dmap(...)
  panel.3dbars(...)
}, xbase = bar_size, ybase = bar_size, col.facet = "black", scales = list(arrows = FALSE, col="black", z=list(cex=1), x=list(draw = FALSE), y=list(draw = FALSE)), zoom = .9, xlim = c(min(map_final$long), max(map_final$long)), ylim = c(min(map_final$lat), max(map_final$lat)), zlim = c(min(seq_summarized_cell$total), max(seq_summarized_cell$total)), xlab = NULL, ylab = NULL, zlab = NULL, aspect = c(diff(c(min(map_final$lat), max(map_final$lat)))/diff(c(min(map_final$long), max(map_final$long))), 0.3), panel.aspect = 0.75, lwd = 0.1, screen = list(z = 15, x = -60), par.settings = list(axis.line = list(col = "transparent"), box.3d = list(col = "transparent", alpha = 0)))

## Save the plot in pdf format
pdf(file=paste0(finalDir, "PanelA.pdf"), width=8, height=8, bg="white", family="Helvetica", pointsize=6)
print(plot)
dev.off()

## Number of sequences per continent
      # Assign points to a country
load(paste0(dataDir, "Counts_08_updated.RData"))
count_data <- counts_08
count_data_pts <- data.frame(count_data[rep(seq_len(dim(count_data)[1]), count_data$total), 1:2, drop = FALSE], row.names=NULL)
count_data_pts <- count_data_pts[,c(2,1)]
count_data_pts <- SpatialPoints(count_data_pts, proj4string=map_poly@proj4string)

count_country <- over(count_data_pts, map_poly)
count_country <- cbind(count_data_pts@coords, count_country)

      # Mexico as N. America
count_country[which(count_country$SOVEREIGNT == "Mexico"), "continent"] <- "North America"

      # Add level "Europe"
levels(count_country$continent) <- c(levels(count_country$continent), "Europe")

      # Handle points that fall outside countries
count_country[which(count_country$long < (-85.53) & count_country$lat > 18.6), "continent"] <- "North America"
count_country[which(is.na(count_country$continent) & count_country$lat == 18.458 & count_country$long == -88.292), "continent"] <- "South America"
count_country[which(is.na(count_country$continent) & count_country$lat == 18.208 & count_country$long == -87.875), "continent"] <- "North America"
count_country[which(is.na(count_country$continent) & count_country$lat == 18.292 & count_country$long == -87.875), "continent"] <- "North America"
count_country[which(is.na(count_country$continent) & count_country$lat == 18.292 & count_country$long == -87.792), "continent"] <- "North America"
count_country[which(is.na(count_country$continent) & count_country$lat > 25.138 & count_country$long > -82.92 & count_country$long < -37.62), "continent"] <- "North America"
count_country[which(is.na(count_country$continent) & count_country$lat == 35.875 & count_country$long == -5.375), "continent"] <- "Africa"
count_country[which(is.na(count_country$continent) & count_country$lat > 37.9 & count_country$long > -14.68 & count_country$long < 24.78), "continent"] <- "Europe"
count_country[which(is.na(count_country$continent) & count_country$lat == 45.458 & count_country$long == 36.125), "continent"] <- "Europe"
count_country[which(is.na(count_country$continent) & count_country$lat == 44.708 & count_country$long == 37.458), "continent"] <- "Europe"
count_country[which(is.na(count_country$continent) & count_country$lat == 43.375 & count_country$long == 39.958), "continent"] <- "Europe"
count_country[which(is.na(count_country$continent) & count_country$lat == 66.542 & count_country$long == 34.958), "continent"] <- "Europe"
count_country[which(is.na(count_country$continent) & count_country$lat == 39.208 & count_country$long == 25.958), "continent"] <- "Europe"
count_country[which(is.na(count_country$continent) & count_country$lat == 39.125 & count_country$long == 26.208), "continent"] <- "Europe"
count_country[which(is.na(count_country$continent) & count_country$lat == 37.708 & count_country$long == 26.792), "continent"] <- "Europe"
count_country[which(is.na(count_country$continent) & count_country$lat == 36.208 & count_country$long == 28.125), "continent"] <- "Europe"
count_country[which(is.na(count_country$continent) & count_country$lat == 35.958 & count_country$long == 14.458), "continent"] <- "Europe"
count_country[which(is.na(count_country$continent) & count_country$lat == 36.792 & count_country$long == 11.958), "continent"] <- "Europe"

no_country <- count_country[which(is.na(count_country$continent)),]
no_country <- SpatialPoints(no_country[,c(1:2)])

#plot(map_poly)
#plot(no_country, add=T, pch=16, col="red")

      # Distinguish between Europe and Asia
count_country[!is.na(count_country$SOVEREIGNT), "continent2"] <- factor(countrycode(sourcevar = count_country[which(!is.na(count_country$SOVEREIGNT)), "SOVEREIGNT"], origin = "country.name", destination = "continent"))
count_country[which(count_country$continent2 == "Europe"), "continent"] <- "Europe"

      # Percentage of sequences that are in N. America & Europe
100*nrow(count_country[which(count_country$continent == "North America" | count_country$continent == "Europe"),])/nrow(count_country)

#### Panel B - Taxonomic bias ####
# Settings
th <- 10
cols <- c("lightsteelblue4", "lightsteelblue3", "lightsteelblue2")
cols2 <- list(mammals='#636363', insects='#252525', fish='#cccccc', birds='#969696')
phylopic.ids <- c('b36a215a-adb3-445d-b364-1e63dddd6950','42fdc3cb-37fc-4340-bdf9-eed8e050137c','7a6448e5-09c4-40c8-8378-599d7f974bfe','5aeaf558-3c48-4173-83b4-dbf2846f8d75')
symbols <- sapply(X=1:length(phylopic.ids), FUN=function(x) image_data(phylopic.ids[x], size = "512")[[1]])

# Summary of taxonomic information
data <- seq %>%
  group_by(class, order) %>%
  summarize(NumSeq = n())
data <- data[order(data$class, -data$NumSeq),]

final <- data.frame(class = character(),
                    order = character(),
                    NumSeq = integer(),
                    Keep = character(),
                    Col = character(),
                    Col2 = character(),
                    stringsAsFactors = F)

classes_ordered <- c("mammals", "insects", "fish", "birds")
for (i in unique(data$class)){
  data_class <- data[which(data$class == i),]
  total <- sum(data_class$NumSeq)
  th_class <- (th/100)*total
  data_class$Keep <- sapply(data_class$NumSeq, function(x) ifelse(x >= th_class, "Y", "N"))
  data_class[nrow(data_class)+1, c("class", "order", "NumSeq", "Keep")] <- c(i, paste(nrow(data_class[which(data_class$Keep == "N"),]), "other orders", sep=" "), sum(data_class[which(data_class$Keep == "N"),"NumSeq"]), "Y")
  data_class$NumSeq <- as.integer(data_class$NumSeq)
  data_class <- data_class[which(data_class$Keep == "Y"),]
  data_class[1:(nrow(data_class)-1), "Col"] <- cols[1:nrow(data_class)-1]
  data_class[1:(nrow(data_class)-1), "Col2"] <- cols2[[i]]
  data_class[nrow(data_class), c("Col", "Col2")] <- "white"
  
  final <- bind_rows(final, data_class)
}
final <- final[, c("class", "order", "NumSeq", "Col", "Col2")]
final$order <- factor(final$order, levels=final$order)
final$order <- fct_rev(final$order)
final$class <- factor(final$class, levels=classes_ordered)

# Plot it
p <- ggplot(data = final, aes(x = class, y = NumSeq, fill=order, label=paste0("  ", order, "  "))) +
  scale_fill_manual(values=rev(final$Col)) +
  geom_bar(stat = "identity", color = "black", size = 0.1, width=0.4) + 
  coord_flip() +
  geom_text_repel(force=1, position = position_stack(vjust=0.5), vjust=3, hjust=0.5, direction="x", ylim  = c(-2100, 150000), family="Helvetica", size=6 * 1/72 * 25.4) +
  guides(fill = FALSE) +
  scale_y_continuous(expand=c(0,0), limits=c(-17000, 150000), breaks=c(0, 25000, 50000, 75000, 100000, 125000, 150000)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), axis.line.x = element_blank(), plot.margin = margin(c(0,7,5,0)), text = element_text(size=6, family="Helvetica")) +
  add_phylopic(symbols[[1]], alpha=1, x=1-0.1,  y=-11000, ysize=12000, color="black") +
  add_phylopic(symbols[[2]], alpha=1, x=2-0.1,  y=-11000, ysize=12000, color="black") +
  add_phylopic(symbols[[3]], alpha=1, x=3-0.1,  y=-11000, ysize=12000, color="black") +
  add_phylopic(symbols[[4]], alpha=1, x=4-0.1,  y=-11000, ysize=12000, color="black") +
  labs(x="", y="Number of sequences") +
  annotate(x=.3, xend=.3, y=0, yend=150000, colour="black", geom="segment")

# Save
pdf(file=paste0(finalDir, "PanelB.pdf"), width=3.231, height=1.699)
print(p)
dev.off()

#### Panel C - Temporal bias ####
cols <- c('#969696','#cccccc','#252525','#636363')
taxa <- c('birds','fish','insects','mammals')
scale <- "10"

load(paste0("Data/alldata.RData"))

pdf(paste0(finalDir, 'PanelC.pdf'), pointsize = 8, width = 4, height=3)

plodat <- alldata %>% filter(scale == scale) %>% 
  group_by(taxon,year) %>%
  summarize(npops = n(), nseqs = sum(nseqs)) %>%
  mutate(ptcx = rescale(npops,c(0.5,2.5)))

plot(nseqs~year,plodat,type='n',yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l',ylim=c(0,12),xlim=c(1980,2018))
title(ylab='number of sequences', cex.lab=1)
title(xlab='year', cex.lab=1)
axis(2,cex.axis=1,lwd=0,lwd.ticks=1,labels = c(1,10,100,1000,10000,100000),at=log(c(1,10,100,1000,10000,10000)))
axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(1985,1995,2005,2015))
for(j in 1:4){
  y <- filter(plodat, taxon == taxa[j]) %>% as.data.frame
  y$nseqs <- log(y$nseqs)
  points(nseqs~year,y,pch=16,type='l',lwd=1.2,col=alpha(cols[j],1))
}

dev.off()
   
#### Panels D and E - Histograms ####
## Load final datasets
load(paste0(dataDir, "Counts_08_updated.RData"))
#load("Counts_1.RData")
#load("Counts_2.RData")
#load("Counts_4.RData")

## Set plotting parameters
bin_width_humans <- 375 # Number of bins for "humans"
bin_width_landuse <- 2.5 # Number of bins for "landuse"
ln_humans <- "No" # Log variable "humans"?
ln_percent <- "Yes" # Log percentage (y-axis)?
scale <- "08" # Spatial scale
taxa <- "total" # Taxonomic group
humans_lim <- 15000 # X-axis limit for humans
alpha <- 0.4 # Colour transparency
#colours <- c("gray40", "hotpink")
colours <- c("gray40", "dodgerblue3")

## Format dataset
counts <- get(paste("counts", scale, sep="_"))
colnames(counts)[colnames(counts)==taxa] <- "n_sequences"
counts[,"n_cells"] <- 1
counts$landuse <- counts$landuse * 100
counts <- counts[,c("lat", "long", "year", "n_sequences", "n_cells", "humans", "landuse")]
if(!is.na(humans_lim)){counts[which(counts$humans>humans_lim),"humans"] <- humans_lim+1}

## Plot histograms
# Humans
h <- ggplot(counts, aes(x = humans)) +
  geom_histogram(aes(weight = n_cells, y=..count../sum(..count..), fill="World"), binwidth = bin_width_humans, boundary=0, alpha=alpha, closed="left") +
  geom_histogram(aes(weight = n_sequences, y=..count../sum(..count..), fill="Sample"), binwidth = bin_width_humans, boundary=0, alpha=alpha, closed="left") +
  scale_y_continuous(labels = percent_format()) +
  scale_x_continuous(breaks=c(0, 5000, 10000, 15000), label=c("0", "5000", "10000", ">15000")) +
  scale_fill_manual(name="", values=colours) +
  ylab("fraction") +
  xlab(bquote("humans (inhabitants "~km^-2* ")")) +
  ggtitle("Humans") +
  theme(text = element_text(size=6, family="Helvetica"), plot.title = element_text(hjust = 0.95, margin = margin(t=20, b=-10), face="bold"), panel.border=element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position=c(0.9,0.45), legend.key.size = unit(0.1,"inch"))

if(ln_percent=="Yes"){
  h <- h + 
    #scale_y_continuous(labels = percent_format(), trans = scales::trans_new(name="ln(1+x)", transform=function(x) log(x+1), inverse=function(x) (exp(x)-1)), breaks=c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), limits=c(NA,1))
    scale_y_continuous(labels = percent_format(), trans = scales::trans_new(name="ln(1+100x)", transform=function(x) log(100*x+1), inverse=function(x) ((exp(x)-1)/100)), breaks=c(0, 0.01, 0.1, 1), limits=c(NA,1))
}

pdf(file=paste0(finalDir, "PanelD.pdf", sep=""), width=3.25, height=2, bg="white")
print(h)
dev.off()

# Landuse
l <- ggplot(counts, aes(x = landuse)) +
  geom_histogram(aes(weight = n_cells, y=..count../sum(..count..), fill="World"), binwidth = bin_width_landuse, boundary=0, alpha=alpha, closed="left") +
  geom_histogram(aes(weight = n_sequences, y=..count../sum(..count..), fill="Sample"), binwidth = bin_width_landuse, boundary=0, alpha=alpha, closed="left") +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(name="", values=colours) +
  ylab("fraction") +
  xlab("landuse (% of gridcell land under heavy anthropogenic use)") +
  ggtitle("Landuse") +
  theme(text = element_text(size=6, family="Helvetica"), plot.title = element_text(hjust = 0.95, margin = margin(t=20, b=-10), face="bold"), panel.border=element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position=c(0.9,0.45), legend.key.size = unit(0.1,"inch"))

if(ln_percent=="Yes"){
  l <- l + 
    #scale_y_continuous(labels = percent_format(), trans = scales::trans_new(name="ln(1+x)", transform=function(x) log(x+1), inverse=function(x) (exp(x)-1)), breaks=c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7), limits=c(NA, 0.7))
    scale_y_continuous(labels = percent_format(), trans = scales::trans_new(name="ln(1+100x)", transform=function(x) log(100*x+1), inverse=function(x) ((exp(x)-1)/100)), breaks=c(0, 0.01, 0.1, 1), limits=c(NA,1))
}

pdf(file=paste0(finalDir, "PanelE.pdf", sep=""), width=3.25, height=2, bg="white")
print(l)
dev.off()

# Both plots together
grid.arrange(h, l, ncol=2)
