#### Temporal variation in intra-specific neutral genetic diversity across anthromes
#### Gonzalez Lab project - McGill University - 2016-2017
#### Script by Vincent Fugère

#### Part I: load and format data

rm(list=ls())
library(tidyverse)
library(magrittr)
#library(svMisc)
#library(readxl)

unwanted <- c('\\.',0:9,'BOLD','-','#','_nsp','spnov')
unwanted <- paste(unwanted, collapse = '|')

#### index dataframe to choose land use map closest in time to sequence

#very long  time gradient (going far back in time) to make sure program doesn't crash but in fact
#all sequences used are > 1980
map.yr.idx <- data.frame(
  'seq.yr' = 1800:2017,
  'map.yr' = c(rep(1980,190),rep(1990,10),2000:2017)
)

distinct(map.yr.idx, map.yr) %>% as.data.frame(.) -> map.yrs
map.yrs <- as.numeric(map.yrs$map.yr)

#### scale of aggregation = 0.08' grid cells ####

#load and format maps
for(i in 1:length(map.yrs)){
  fname <- paste0('hyde32_',map.yrs[i],'_08')
  path <- paste0('~/Google Drive/recherche/Intraspecific genetic diversity/data/',fname,'.RData')
  load(path)
  if(i == 1){
    get(fname) %>% filter(!is.na(pop_tot)) %>% mutate_at(vars(lat:long), funs(round(as.numeric(.),3))) %>%
      mutate('cell' = paste0(.data$lat,'_',.data$long)) %>% select(-lat, -long) %>%
      mutate('mapyear' = map.yrs[i]) %>% select(mapyear,cell,everything()) -> lu.map
  } else {
    get(fname) %>% filter(!is.na(pop_tot)) %>% mutate_at(vars(lat:long), funs(round(as.numeric(.),3))) %>%
      mutate('cell' = paste0(.data$lat,'_',.data$long)) %>% select(-lat, -long) %>%
      mutate('mapyear' = map.yrs[i]) %>% select(mapyear,cell,everything()) -> tmp.map
    lu.map <- bind_rows(lu.map,tmp.map)
  }
  rm(list = c(fname))
}
rm(tmp.map)

# load('~/Google Drive/recherche/Intraspecific genetic diversity/data/hyde32_2017_08.RData')
# load('~/Google Drive/recherche/Intraspecific genetic diversity/data/hyde32_1980_08.RData')
# old <- hyde32_1980_08 %>% mutate_at(vars(lat:long), funs(round(as.numeric(.),3))) %>%
#   mutate('cell' = paste0(.data$lat,'_',.data$long)) %>%
#   #rename('pop_density' = pop_tot, 'p_pasture' = grazing, 'p_urban' = urban, 'p_crop' = cropland) %>%
#   mutate('cell.r' = as.factor(cell)) %>%
#   select(-c(lat,long,cell))
# lu08 <- hyde32_2017_08 %>% mutate_at(vars(lat:long), funs(round(as.numeric(.),3))) %>%
#   mutate('cell' = paste0(.data$lat,'_',.data$long)) %>%
#   #rename('pop_density' = pop_tot, 'p_pasture' = grazing, 'p_urban' = urban, 'p_crop' = cropland) %>%
#   mutate('cell.r' = as.factor(cell))
# lu08 <- left_join(x = lu08, y = old, by = 'cell.r' , suffix = c('.2017','.1980')) %>%
#   select(lat,long,cell,cell.r,everything())
# rm(hyde32_1980_08, hyde32_2017_08,old)

# lu08 <- read_csv('~/Google Drive/recherche/Intraspecific genetic diversity/data/landuse_08.csv') %>%
#   select(-starts_with('X'))
# lu08$cell <- gsub(' ','',lu08$cell)
# lu08 <- lu08 %>% mutate(lat = gsub("_.*","", cell), long = gsub(".*_","", cell)) %>%
#   mutate_at(vars(lat:long), funs(round(as.numeric(.),3))) %>%
#   mutate(cell.r = paste(lat,'_',long,sep=''))
# lu08$cell.r <-as.factor(lu08$cell.r)
#lu08 <- lu08 %>% mutate_at(vars(pop_density:p_irrigated), as.numeric)
#save(lu08, file = '~/Google Drive/recherche/Intraspecific genetic diversity/data/landuse08.RData')
#load('~/Google Drive/recherche/Intraspecific genetic diversity/data/landuse08.RData')

## MAMMALS

#mam08_1: new file
mam08 <- read.csv('~/Google Drive/recherche/Intraspecific genetic diversity/data/mamm_pairwise_0.08by0.08_anth.csv')
# #mam08_2: older file
# mam08 <- read.csv('/Users/vincentfugere/Dropbox/Gonzalez lab meetings/Anthropocene_data_analysis/Vincent/raw files/mamm_pairwise_560by560_anth.csv')

colnames(mam08)[1] <- 'species.year.ID' #to be consistent with other files
colnames(mam08)[10] <- 'anth' #change new name back to old to be consistent with analysis script

# add year & species (split column)
info <- unlist(strsplit(as.character(mam08$species.year.ID), '[.]'))
mam08$year <- info[seq(2, length(info), by = 2)]
mam08$year <- as.numeric(mam08$year)
mam08$species <- info[seq(1, length(info), by = 2)]
mam08$species <- as.factor(mam08$species)
rm(info)

#mam08 <- mam08[mam08$num_per_bp != 0,] # removes all comparisons with 0 diffs = could be redundant sequences
mam08 <- na.omit(mam08) %>% filter(str_count(.$species, '_') == 1, str_count(.$species, unwanted) == 0) #removes NAs & misidentified species
mam08 <- droplevels(mam08)

#summarizing data frame
mam08.agg <- mam08 %>% group_by(species, year, cell) %>% summarize(
  div = mean(num_per_bp),
  anthrome = mean(anth),
  ncomps = length(num_per_bp),
  nseqs = length(unique(c(seq1,seq2)))) %>%
  mutate_at(vars(cell), as.character) %>%
  mutate(lat = gsub("_.*","", cell), long = gsub(".*_","", cell)) %>%
  mutate_at(vars(lat:long), funs(round(as.numeric(.),3))) %>%
  mutate(cell = paste0(lat,'_',long)) %>%
  select(-anthrome) %>% as.data.frame(.)

#checking if all cells from which we have sequences are found in land use map
if(sum(!(mam08.agg$cell %in% lu.map[lu.map$mapyear == map.yrs[1],'cell'])) > 0) {print('not all cells matched')}

#adding land use and human population density infomation
mam08.agg$mapyear <- map.yr.idx$map.yr[match(mam08.agg$year, map.yr.idx$seq.yr)]
mam08.agg <- inner_join(mam08.agg,lu.map,by=c('mapyear','cell'))

# removing groups with only one sequence comparison (poor measure of diversity)
#mam08.agg <- mam08.agg[mam08.agg$ncomps > 1,]

#creating unique identifier for population
mam08.agg$pop <- as.factor(paste(mam08.agg$species,mam08.agg$cell,sep='_'))

#creating unique identifier for species/year combination
mam08.agg$sp_yr <- as.factor(paste(mam08.agg$species,mam08.agg$year,sep='_'))

#adding length of time series
pop.dur <- aggregate(div~pop, mam08.agg, FUN=length)
mam08.agg$n.years <- pop.dur$div[match(mam08.agg$pop,pop.dur$pop)]
rm(pop.dur)

# # excluding pops in the water
# mam08.agg <- mam08.agg[mam08.agg$anthrome != 0,]
# 
# #converting anthrome to binary. natural are all 'semi-natural' and 'wild' in Ellis 2010
# mam08.agg$anth.cat <- 'anthro'
# mam08.agg$anth.cat[mam08.agg$anthrome >= 50] <- 'natural'
# mam08.agg$anth.cat <- as.factor(mam08.agg$anth.cat)
# 
# mam08.agg <- droplevels(mam08.agg)

# # mam08.short <- mam08.agg[!(duplicated(mam08.agg$pop)),]
# nlevels(mam08.agg$species) #number of species
# nlevels(mam08.agg$pop) #number of populations
# # how many pops have data from same year, different anthromes?
# mam08.tab <- with(mam08.agg,table(sp_yr, anth.cat))
# mam08.tab <- mam08.tab[mam08.tab[,1] > 0 & mam08.tab[,2] > 0,]
# mam08.tx.spatial <- row.names(mam08.tab)
# rm(mam08.tab)
# #number of data points for spatial analysis within species @ .28
# length(mam08.tx.spatial)

# #time series duration
# time.series08 <- data.frame('nb.year' = 1:20)
# time.series08$mams <- numeric(length(time.series08$nb.year))
# for(i in 1:dim(time.series08)[1]){
#   time.series08$mams[i] <- length(unique(mam08.agg[mam08.agg$n.years == time.series08$nb.year[i],]$pop))
# }

## BIRDS

aves08 <- read.csv('~/Google Drive/recherche/Intraspecific genetic diversity/data/aves_pairwise_0.08by0.08_anth.csv')
colnames(aves08)[1] <- 'species.year.ID' #to be consistent with other files
colnames(aves08)[10] <- 'anth' #change new name back to old to be consistent with analysis script

# add year & species (split column)
info <- unlist(strsplit(as.character(aves08$species.year.ID), '[.]'))
aves08$year <- info[seq(2, length(info), by = 2)]
aves08$year <- as.numeric(aves08$year)
aves08$species <- info[seq(1, length(info), by = 2)]
aves08$species <- as.factor(aves08$species)
rm(info)

#aves08 <- aves08[aves08$num_per_bp != 0,] # removes all comparisons with 0 diffs = could be redundant sequences
aves08 <- na.omit(aves08) %>% filter(str_count(.$species, '_') == 1, str_count(.$species, unwanted) == 0) #removes NAs & misidentified species
aves08 <- droplevels(aves08)

#summarizing data frame
aves08.agg <- aves08 %>% group_by(species, year, cell) %>% summarize(
  div = mean(num_per_bp),
  anthrome = mean(anth),
  ncomps = length(num_per_bp),
  nseqs = length(unique(c(seq1,seq2)))) %>%
  mutate_at(vars(cell), as.character) %>%
  mutate(lat = gsub("_.*","", cell), long = gsub(".*_","", cell)) %>%
  mutate_at(vars(lat:long), funs(round(as.numeric(.),3))) %>%
  mutate(cell = paste0(lat,'_',long)) %>%
  select(-anthrome) %>% as.data.frame(.)

#checking if all cells from which we have sequences are found in land use map
if(sum(!(aves08.agg$cell %in% lu.map[lu.map$mapyear == map.yrs[1],'cell'])) > 0) {print('not all cells matched')}

#adding land use and human population density infomation
aves08.agg$mapyear <- map.yr.idx$map.yr[match(aves08.agg$year, map.yr.idx$seq.yr)]
aves08.agg <- inner_join(aves08.agg,lu.map,by=c('mapyear','cell'))

# removing groups with only one sequence comparison (poor measure of diversity)
#aves08.agg <- aves08.agg[aves08.agg$ncomps > 1,]

#creating unique identifier for population
aves08.agg$pop <- as.factor(paste(aves08.agg$species,aves08.agg$cell,sep='_'))

#creating unique identifier for species/year combination
aves08.agg$sp_yr <- as.factor(paste(aves08.agg$species,aves08.agg$year,sep='_'))

#adding length of time series
pop.dur <- aggregate(div~pop, aves08.agg, FUN=length)
aves08.agg$n.years <- pop.dur$div[match(aves08.agg$pop,pop.dur$pop)]
rm(pop.dur)

# # excluding pops in the water
# aves08.agg <- aves08.agg[aves08.agg$anthrome != 0,]
# 
# #converting anthrome to binary
# aves08.agg$anth.cat <- 'anthro'
# aves08.agg$anth.cat[aves08.agg$anthrome >= 50] <- 'natural'
# aves08.agg$anth.cat <- as.factor(aves08.agg$anth.cat)
# 
# aves08.agg <- droplevels(aves08.agg)
# 
# # aves08.short <- aves08.agg[!(duplicated(aves08.agg$pop)),]
# nlevels(aves08.agg$species) #number of species
# nlevels(aves08.agg$pop) #number of populations
# # how many pops have data from same year, different anthromes?
# aves08.tab <- with(aves08.agg,table(sp_yr, anth.cat))
# aves08.tab <- aves08.tab[aves08.tab[,1] > 0 & aves08.tab[,2] > 0,]
# aves08.tx.spatial <- row.names(aves08.tab)
# rm(aves08.tab)
# #number of data points for spatial analysis within species @ .28
# length(aves08.tx.spatial)
# 
# #time series duration
# time.series08$aves <- numeric(length(time.series08$nb.year))
# for(i in 1:dim(time.series08)[1]){
#   time.series08$aves[i] <- length(unique(aves08.agg[aves08.agg$n.years == time.series08$nb.year[i],]$pop))
# }

## BONY FISHES

acti08 <- read.csv('~/Google Drive/recherche/Intraspecific genetic diversity/data/acti_pairwise_0.08by0.08_anth.csv')
colnames(acti08)[1] <- 'species.year.ID' #to be consistent with other files
colnames(acti08)[10] <- 'anth' #change new name back to old to be consistent with analysis script

# add year & species (split column)
info <- unlist(strsplit(as.character(acti08$species.year.ID), '[.]'))
acti08$year <- info[seq(2, length(info), by = 2)]
acti08$year <- as.numeric(acti08$year)
acti08$species <- info[seq(1, length(info), by = 2)]
acti08$species <- as.factor(acti08$species)
rm(info)

# #removing marine species
# library(rfishbase)
# sp.list <- as.character(levels(acti08$species))
# sp.list <- gsub('_', ' ', sp.list)
# ecology(sp.list, limit=1,)

#acti08 <- acti08[acti08$num_per_bp != 0,] # removes all comparisons with 0 diffs = could be redundant sequences
acti08 <- na.omit(acti08) %>% filter(str_count(.$species, '_') == 1, str_count(.$species, unwanted) == 0) #removes NAs & misidentified species
acti08 <- droplevels(acti08)

#summarizing data frame
acti08.agg <- acti08 %>% group_by(species, year, cell) %>% summarize(
  div = mean(num_per_bp),
  anthrome = mean(anth),
  ncomps = length(num_per_bp),
  nseqs = length(unique(c(seq1,seq2)))) %>%
  mutate_at(vars(cell), as.character) %>%
  mutate(lat = gsub("_.*","", cell), long = gsub(".*_","", cell)) %>%
  mutate_at(vars(lat:long), funs(round(as.numeric(.),3))) %>%
  mutate(cell = paste0(lat,'_',long)) %>%
  select(-anthrome) %>% as.data.frame(.)

#checking if all cells from which we have sequences are found in land use map
if(sum(!(acti08.agg$cell %in% lu.map[lu.map$mapyear == map.yrs[1],'cell'])) > 0) {print('not all cells matched')}

#adding land use and human population density infomation
acti08.agg$mapyear <- map.yr.idx$map.yr[match(acti08.agg$year, map.yr.idx$seq.yr)]
acti08.agg <- inner_join(acti08.agg,lu.map,by=c('mapyear','cell'))


# removing groups with only one sequence comparison (poor measure of diversity)
#acti08.agg <- acti08.agg[acti08.agg$ncomps > 1,]

#creating unique identifier for population
acti08.agg$pop <- as.factor(paste(acti08.agg$species,acti08.agg$cell,sep='_'))

#creating unique identifier for species/year combination
acti08.agg$sp_yr <- as.factor(paste(acti08.agg$species,acti08.agg$year,sep='_'))

#adding length of time series
pop.dur <- aggregate(div~pop, acti08.agg, FUN=length)
acti08.agg$n.years <- pop.dur$div[match(acti08.agg$pop,pop.dur$pop)]
rm(pop.dur)

# # excluding pops in the ocean
# acti08.agg <- acti08.agg[acti08.agg$anthrome != 0,]
# 
# #converting anthrome to binary
# acti08.agg$anth.cat <- 'anthro'
# acti08.agg$anth.cat[acti08.agg$anthrome >= 50] <- 'natural'
# acti08.agg$anth.cat <- as.factor(acti08.agg$anth.cat)
# 
# acti08.agg <- droplevels(acti08.agg)
# 
# # acti08.short <- acti08.agg[!(duplicated(acti08.agg$pop)),]
# nlevels(acti08.agg$species) #number of species
# nlevels(acti08.agg$pop) #number of populations
# # how many pops have data from same year, different anthromes?
# acti08.tab <- with(acti08.agg,table(sp_yr, anth.cat))
# acti08.tab <- acti08.tab[acti08.tab[,1] > 0 & acti08.tab[,2] > 0,]
# acti08.tx.spatial <- row.names(acti08.tab)
# rm(acti08.tab)
# #number of data points for spatial analysis within species @ .28
# length(acti08.tx.spatial)
# 
# #time series duration
# time.series08$acti <- numeric(length(time.series08$nb.year))
# for(i in 1:dim(time.series08)[1]){
#   time.series08$acti[i] <- length(unique(acti08.agg[acti08.agg$n.years == time.series08$nb.year[i],]$pop))
# }

## INSECTS

insect08 <- read.csv('~/Google Drive/recherche/Intraspecific genetic diversity/data/insect_pairwise_0.08by0.08_anth.csv')
colnames(insect08)[1] <- 'species.year.ID' #to be consistent with other files
colnames(insect08)[10] <- 'anth' #change new name back to old to be consistent with analysis script

insect08 <- insect08[insect08$overlap > 0.5,] # filter comparisons with less than 50% overlap as in Miraldo paper (for other taxa, this step was done in Excel)

row2rm <- which(insect08$species.year.ID == 'Lamprospilus_aff..2015')
if (length(row2rm) > 0){
  insect08 <- insect08[-row2rm,] #removes database entry with incorrect species name
}

#insect08 <- insect08[insect08$num_per_bp != 0,] # removes all comparisons with 0 diffs = could be redundant sequences
insect08 <- na.omit(insect08) #removes NAs

# add year & species (split column)
info <- unlist(strsplit(as.character(insect08$species.year.ID), '[.]'))
insect08$year <- info[seq(2, length(info), by = 2)]
insect08$year <- as.numeric(insect08$year)
insect08$species <- info[seq(1, length(info), by = 2)]
insect08$species <- as.factor(insect08$species)
rm(info)

insect08 <- insect08 %>% filter(str_count(.$species, '_') == 1, str_count(.$species, unwanted) == 0)

insect08 <- droplevels(insect08)

#summarizing data frame
insect08.agg <- insect08 %>% group_by(species, year, cell) %>% summarize(
  div = mean(num_per_bp),
  anthrome = mean(anth),
  ncomps = length(num_per_bp),
  nseqs = length(unique(c(seq1,seq2)))) %>%
  mutate_at(vars(cell), as.character) %>%
  mutate(lat = gsub("_.*","", cell), long = gsub(".*_","", cell)) %>%
  mutate_at(vars(lat:long), funs(round(as.numeric(.),3))) %>%
  mutate(cell = paste0(lat,'_',long)) %>%
  select(-anthrome) %>% as.data.frame(.)

#checking if all cells from which we have sequences are found in land use map
if(sum(!(insect08.agg$cell %in% lu.map[lu.map$mapyear == map.yrs[1],'cell'])) > 0) {print('not all cells matched')}

#adding land use and human population density infomation
insect08.agg$mapyear <- map.yr.idx$map.yr[match(insect08.agg$year, map.yr.idx$seq.yr)]
insect08.agg <- inner_join(insect08.agg,lu.map,by=c('mapyear','cell'))

# removing groups with only one sequence comparison (poor measure of diversity)
#insect08.agg <- insect08.agg[insect08.agg$ncomps > 1,]

#creating unique identifier for population
insect08.agg$pop <- as.factor(paste(insect08.agg$species,insect08.agg$cell,sep='_'))

#creating unique identifier for species/year combination
insect08.agg$sp_yr <- as.factor(paste(insect08.agg$species,insect08.agg$year,sep='_'))

#adding length of time series
pop.dur <- aggregate(div~pop, insect08.agg, FUN=length)
insect08.agg$n.years <- pop.dur$div[match(insect08.agg$pop,pop.dur$pop)]
rm(pop.dur)

# # excluding pops in the ocean
# insect08.agg <- insect08.agg[insect08.agg$anthrome != 0,]
# 
# #converting anthrome to binary
# insect08.agg$anth.cat <- 'anthro'
# insect08.agg$anth.cat[insect08.agg$anthrome >= 50] <- 'natural'
# insect08.agg$anth.cat <- as.factor(insect08.agg$anth.cat)
# 
# insect08.agg <- droplevels(insect08.agg)
# 
# # insect08.short <- insect08.agg[!(duplicated(insect08.agg$pop)),]
# nlevels(insect08.agg$species) #number of species
# nlevels(insect08.agg$pop) #number of populations
# # how many pops have data from same year, different anthromes?
# insect08.tab <- with(insect08.agg,table(sp_yr, anth.cat))
# insect08.tab <- insect08.tab[insect08.tab[,1] > 0 & insect08.tab[,2] > 0,]
# insect08.tx.spatial <- row.names(insect08.tab)
# rm(insect08.tab)
# #number of data points for spatial analysis within species @ .28
# length(insect08.tx.spatial)
# 
# #time series duration
# time.series08$insects <- numeric(length(time.series08$nb.year))
# for(i in 1:dim(time.series08)[1]){
#   time.series08$insects[i] <- length(unique(insect08.agg[insect08.agg$n.years == time.series08$nb.year[i],]$pop))
# }

## PLANTS (ITS2)

plant.its08 <- read.csv('~/Google Drive/recherche/Intraspecific genetic diversity/data/plants_ITS_pairwise_0.08by0.08_anth.csv')
colnames(plant.its08)[1] <- 'species.year.ID' #to be consistent with other files
colnames(plant.its08)[10] <- 'anth' #change new name back to old to be consistent with analysis script

# add year & species (split column)
info <- unlist(strsplit(as.character(plant.its08$species.year.ID), '[.]'))
plant.its08$year <- info[seq(2, length(info), by = 2)]
plant.its08$year <- as.numeric(plant.its08$year)
plant.its08$species <- info[seq(1, length(info), by = 2)]
plant.its08$species <- as.factor(plant.its08$species)
rm(info)

#plant.its08 <- plant.its08[plant.its08$num_per_bp != 0,] # removes all comparisons with 0 diffs = could be redundant sequences
plant.its08 <- na.omit(plant.its08) %>% filter(str_count(.$species, '_') == 1, str_count(.$species, unwanted) == 0) #removes NAs & misidentified species
plant.its08 <- droplevels(plant.its08)

#summarizing data frame
plant.its08.agg <- plant.its08 %>% group_by(species, year, cell) %>% summarize(
  div = mean(num_per_bp),
  anthrome = mean(anth),
  ncomps = length(num_per_bp),
  nseqs = length(unique(c(seq1,seq2)))) %>%
  mutate_at(vars(cell), as.character) %>%
  mutate(lat = gsub("_.*","", cell), long = gsub(".*_","", cell)) %>%
  mutate_at(vars(lat:long), funs(round(as.numeric(.),3))) %>%
  mutate(cell = paste0(lat,'_',long)) %>%
  select(-anthrome) %>% as.data.frame(.)

#checking if all cells from which we have sequences are found in land use map
if(sum(!(plant.its08.agg$cell %in% lu.map[lu.map$mapyear == map.yrs[1],'cell'])) > 0) {print('not all cells matched')}

#adding land use and human population density infomation
plant.its08.agg$mapyear <- map.yr.idx$map.yr[match(plant.its08.agg$year, map.yr.idx$seq.yr)]
plant.its08.agg <- inner_join(plant.its08.agg,lu.map,by=c('mapyear','cell'))

# removing groups with only one sequence comparison (poor measure of diversity)
#plant.its08.agg <- plant.its08.agg[plant.its08.agg$ncomps > 1,]

#creating unique identifier for population
plant.its08.agg$pop <- as.factor(paste(plant.its08.agg$species,plant.its08.agg$cell,sep='_'))

#creating unique identifier for species/year combination
plant.its08.agg$sp_yr <- as.factor(paste(plant.its08.agg$species,plant.its08.agg$year,sep='_'))

#adding length of time series
pop.dur <- aggregate(div~pop, plant.its08.agg, FUN=length)
plant.its08.agg$n.years <- pop.dur$div[match(plant.its08.agg$pop,pop.dur$pop)]
rm(pop.dur)

# # excluding pops in the water
# plant.its08.agg <- plant.its08.agg[plant.its08.agg$anthrome != 0,]
# 
# #converting anthrome to binary
# plant.its08.agg$anth.cat <- 'anthro'
# plant.its08.agg$anth.cat[plant.its08.agg$anthrome >= 50] <- 'natural'
# plant.its08.agg$anth.cat <- as.factor(plant.its08.agg$anth.cat)
# 
# plant.its08.agg <- droplevels(plant.its08.agg)
# 
# # plant.its08.short <- plant.its08.agg[!(duplicated(plant.its08.agg$pop)),]
# nlevels(plant.its08.agg$species) #number of species
# nlevels(plant.its08.agg$pop) #number of populations
# # how many pops have data from same year, different anthromes?
# plant.its08.tab <- with(plant.its08.agg,table(sp_yr, anth.cat))
# plant.its08.tab <- plant.its08.tab[plant.its08.tab[,1] > 0 & plant.its08.tab[,2] > 0,]
# plant.its08.tx.spatial <- row.names(plant.its08.tab)
# rm(plant.its08.tab)
# #number of data points for spatial analysis within species @ .28
# length(plant.its08.tx.spatial)
# 
# #time series duration
# time.series08$plant.its <- numeric(length(time.series08$nb.year))
# for(i in 1:dim(time.series08)[1]){
#   time.series08$plant.its[i] <- length(unique(plant.its08.agg[plant.its08.agg$n.years == time.series08$nb.year[i],]$pop))
# }
# sum(time.series08$plant.its[5:20]) #nb of time series

## PLANTS (MATK)

plant.matK08 <- read.csv('~/Google Drive/recherche/Intraspecific genetic diversity/data/plants_pairwise_matK_0.08by0.08_anth.csv')
colnames(plant.matK08)[1] <- 'species.year.ID' #to be consistent with other files
colnames(plant.matK08)[10] <- 'anth' #change new name back to old to be consistent with analysis script

#removing line with double period
#sapply(gregexpr("\\.", as.character(plant.matK08$species.year.ID)), tail, 1) == sapply(gregexpr("\\.", as.character(plant.matK08$species.year.ID)), head, 1)
row2rm <- which(plant.matK08$species.year.ID == 'Crataegus_indet..2011')
if (length(row2rm) > 0){
  plant.matK08 <- plant.matK08[-row2rm,] #removes database entry with incorrect species name
}

# add year & species (split column)
info <- unlist(strsplit(as.character(plant.matK08$species.year.ID), '[.]'))
plant.matK08$year <- info[seq(2, length(info), by = 2)]
plant.matK08$year <- as.numeric(plant.matK08$year)
plant.matK08$species <- info[seq(1, length(info), by = 2)]
plant.matK08$species <- as.factor(plant.matK08$species)
rm(info)

#plant.matK08 <- plant.matK08[plant.matK08$num_per_bp != 0,] # removes all comparisons with 0 diffs = could be redundant sequences
plant.matK08 <- na.omit(plant.matK08) %>% filter(str_count(.$species, '_') == 1, str_count(.$species, unwanted) == 0) #removes NAs & misidentified species
plant.matK08 <- droplevels(plant.matK08)

#summarizing data frame
plant.matK08.agg <- plant.matK08 %>% group_by(species, year, cell) %>% summarize(
  div = mean(num_per_bp),
  anthrome = mean(anth),
  ncomps = length(num_per_bp),
  nseqs = length(unique(c(seq1,seq2)))) %>%
  mutate_at(vars(cell), as.character) %>%
  mutate(lat = gsub("_.*","", cell), long = gsub(".*_","", cell)) %>%
  mutate_at(vars(lat:long), funs(round(as.numeric(.),3))) %>%
  mutate(cell = paste0(lat,'_',long)) %>%
  select(-anthrome) %>% as.data.frame(.)

#checking if all cells from which we have sequences are found in land use map
if(sum(!(plant.matK08.agg$cell %in% lu.map[lu.map$mapyear == map.yrs[1],'cell'])) > 0) {print('not all cells matched')}

#adding land use and human population density infomation
plant.matK08.agg$mapyear <- map.yr.idx$map.yr[match(plant.matK08.agg$year, map.yr.idx$seq.yr)]
plant.matK08.agg <- inner_join(plant.matK08.agg,lu.map,by=c('mapyear','cell'))

# removing groups with only one sequence comparison (poor measure of diversity)
#plant.matK08.agg <- plant.matK08.agg[plant.matK08.agg$ncomps > 1,]

#creating unique identifier for population
plant.matK08.agg$pop <- as.factor(paste(plant.matK08.agg$species,plant.matK08.agg$cell,sep='_'))

#creating unique identifier for species/year combination
plant.matK08.agg$sp_yr <- as.factor(paste(plant.matK08.agg$species,plant.matK08.agg$year,sep='_'))

#adding length of time series
pop.dur <- aggregate(div~pop, plant.matK08.agg, FUN=length)
plant.matK08.agg$n.years <- pop.dur$div[match(plant.matK08.agg$pop,pop.dur$pop)]
rm(pop.dur)

# # excluding pops in the water
# plant.matK08.agg <- plant.matK08.agg[plant.matK08.agg$anthrome != 0,]
# 
# #converting anthrome to binary
# plant.matK08.agg$anth.cat <- 'anthro'
# plant.matK08.agg$anth.cat[plant.matK08.agg$anthrome >= 50] <- 'natural'
# plant.matK08.agg$anth.cat <- as.factor(plant.matK08.agg$anth.cat)
# 
# plant.matK08.agg <- droplevels(plant.matK08.agg)
# 
# # plant.matK08.short <- plant.matK08.agg[!(duplicated(plant.matK08.agg$pop)),]
# nlevels(plant.matK08.agg$species) #number of species
# nlevels(plant.matK08.agg$pop) #number of populations
# # how many pops have data from same year, different anthromes?
# plant.matK08.tab <- with(plant.matK08.agg,table(sp_yr, anth.cat))
# plant.matK08.tab <- plant.matK08.tab[plant.matK08.tab[,1] > 0 & plant.matK08.tab[,2] > 0,]
# plant.matK08.tx.spatial <- row.names(plant.matK08.tab)
# rm(plant.matK08.tab)
# #number of data points for spatial analysis within species @ .28
# length(plant.matK08.tx.spatial)
# 
# #time series duration
# time.series08$plant.matK <- numeric(length(time.series08$nb.year))
# for(i in 1:dim(time.series08)[1]){
#   time.series08$plant.matK[i] <- length(unique(plant.matK08.agg[plant.matK08.agg$n.years == time.series08$nb.year[i],]$pop))
# }
# sum(time.series08$plant.matK[5:20]) #nb of time series

## PLANTS (rbcL)

plant.rbcL08 <- read.csv('~/Google Drive/recherche/Intraspecific genetic diversity/data/plants_pairwise_rbcL_0.08by0.08_anth.csv')
colnames(plant.rbcL08)[1] <- 'species.year.ID' #to be consistent with other files
colnames(plant.rbcL08)[10] <- 'anth' #change new name back to old to be consistent with analysis script

#which(sapply(gregexpr("\\.", as.character(plant.rbcL08$species.year.ID)), tail, 1) != sapply(gregexpr("\\.", as.character(plant.rbcL08$species.year.ID)), head, 1))
row2rm <- which(plant.rbcL08$species.year.ID == 'Crataegus_indet..2011')
if (length(row2rm) > 0){
  plant.rbcL08 <- plant.rbcL08[-row2rm,] #removes database entry with incorrect species name
}

# add year & species (split column)
info <- unlist(strsplit(as.character(plant.rbcL08$species.year.ID), '[.]'))
plant.rbcL08$year <- info[seq(2, length(info), by = 2)]
plant.rbcL08$year <- as.numeric(plant.rbcL08$year)
plant.rbcL08$species <- info[seq(1, length(info), by = 2)]
plant.rbcL08$species <- as.factor(plant.rbcL08$species)
rm(info)

#plant.rbcL08 <- plant.rbcL08[plant.rbcL08$num_per_bp != 0,] # removes all comparisons with 0 diffs = could be redundant sequences
plant.rbcL08 <- na.omit(plant.rbcL08) %>% filter(str_count(.$species, '_') == 1, str_count(.$species, unwanted) == 0) #removes NAs & misidentified species
plant.rbcL08 <- droplevels(plant.rbcL08)

#summarizing data frame
plant.rbcL08.agg <- plant.rbcL08 %>% group_by(species, year, cell) %>% summarize(
  div = mean(num_per_bp),
  anthrome = mean(anth),
  ncomps = length(num_per_bp),
  nseqs = length(unique(c(seq1,seq2)))) %>%
  mutate_at(vars(cell), as.character) %>%
  mutate(lat = gsub("_.*","", cell), long = gsub(".*_","", cell)) %>%
  mutate_at(vars(lat:long), funs(round(as.numeric(.),3))) %>%
  mutate(cell = paste0(lat,'_',long)) %>%
  select(-anthrome) %>% as.data.frame(.)

#checking if all cells from which we have sequences are found in land use map
if(sum(!(plant.rbcL08.agg$cell %in% lu.map[lu.map$mapyear == map.yrs[1],'cell'])) > 0) {print('not all cells matched')}

#adding land use and human population density infomation
plant.rbcL08.agg$mapyear <- map.yr.idx$map.yr[match(plant.rbcL08.agg$year, map.yr.idx$seq.yr)]
plant.rbcL08.agg <- inner_join(plant.rbcL08.agg,lu.map,by=c('mapyear','cell'))

# removing groups with only one sequence comparison (poor measure of diversity)
#plant.rbcL08.agg <- plant.rbcL08.agg[plant.rbcL08.agg$ncomps > 1,]

#creating unique identifier for population
plant.rbcL08.agg$pop <- as.factor(paste(plant.rbcL08.agg$species,plant.rbcL08.agg$cell,sep='_'))

#creating unique identifier for species/year combination
plant.rbcL08.agg$sp_yr <- as.factor(paste(plant.rbcL08.agg$species,plant.rbcL08.agg$year,sep='_'))

#adding length of time series
pop.dur <- aggregate(div~pop, plant.rbcL08.agg, FUN=length)
plant.rbcL08.agg$n.years <- pop.dur$div[match(plant.rbcL08.agg$pop,pop.dur$pop)]
rm(pop.dur)

# # excluding pops in the water
# plant.rbcL08.agg <- plant.rbcL08.agg[plant.rbcL08.agg$anthrome != 0,]
# 
# #converting anthrome to binary
# plant.rbcL08.agg$anth.cat <- 'anthro'
# plant.rbcL08.agg$anth.cat[plant.rbcL08.agg$anthrome >= 50] <- 'natural'
# plant.rbcL08.agg$anth.cat <- as.factor(plant.rbcL08.agg$anth.cat)
# 
# plant.rbcL08.agg <- droplevels(plant.rbcL08.agg)
# 
# # plant.rbcL08.short <- plant.rbcL08.agg[!(duplicated(plant.rbcL08.agg$pop)),]
# nlevels(plant.rbcL08.agg$species) #number of species
# nlevels(plant.rbcL08.agg$pop) #number of populations
# # how many pops have data from same year, different anthromes?
# plant.rbcL08.tab <- with(plant.rbcL08.agg,table(sp_yr, anth.cat))
# plant.rbcL08.tab <- plant.rbcL08.tab[plant.rbcL08.tab[,1] > 0 & plant.rbcL08.tab[,2] > 0,]
# plant.rbcL08.tx.spatial <- row.names(plant.rbcL08.tab)
# rm(plant.rbcL08.tab)
# #number of data points for spatial analysis within species @ .28
# length(plant.rbcL08.tx.spatial)
# 
# #time series duration
# time.series08$plant.rbcL <- numeric(length(time.series08$nb.year))
# for(i in 1:dim(time.series08)[1]){
#   time.series08$plant.rbcL[i] <- length(unique(plant.rbcL08.agg[plant.rbcL08.agg$n.years == time.series08$nb.year[i],]$pop))
# }
# sum(time.series08$plant.rbcL[5:20]) #nb of time series

#gather data summary statistics
summ08 <- data.frame('taxon'=c('mammals','birds','fish','insects','plants_its','plants_matK','plants_rbcL'),
                     'nb.sequences' = numeric(7),
                     'nb.species' = numeric(7),
                     'nb.pops' = numeric(7))
summ08[1,2:4] <- c(sum(mam08.agg$nseqs),length(levels(mam08.agg$species)),length(levels(mam08.agg$pop)))
summ08[2,2:4] <- c(sum(aves08.agg$nseqs),length(levels(aves08.agg$species)),length(levels(aves08.agg$pop)))
summ08[3,2:4] <- c(sum(acti08.agg$nseqs),length(levels(acti08.agg$species)),length(levels(acti08.agg$pop)))
summ08[4,2:4] <- c(sum(insect08.agg$nseqs),length(levels(insect08.agg$species)),length(levels(insect08.agg$pop)))
summ08[5,2:4] <- c(sum(plant.its08.agg$nseqs),length(levels(plant.its08.agg$species)),length(levels(plant.its08.agg$pop)))
summ08[6,2:4] <- c(sum(plant.matK08.agg$nseqs),length(levels(plant.matK08.agg$species)),length(levels(plant.matK08.agg$pop)))
summ08[7,2:4] <- c(sum(plant.rbcL08.agg$nseqs),length(levels(plant.rbcL08.agg$species)),length(levels(plant.rbcL08.agg$pop)))

#cleanup
rm(i,row2rm,acti08,aves08,insect08,mam08,plant.its08,plant.matK08,plant.rbcL08,lu.map)

mam08.agg <- as.data.frame(mam08.agg)
aves08.agg <- as.data.frame(aves08.agg)
acti08.agg <- as.data.frame(acti08.agg)
insect08.agg <- as.data.frame(insect08.agg)
plant.its08.agg <- as.data.frame(plant.its08.agg)
plant.matK08.agg <- as.data.frame(plant.matK08.agg)
plant.rbcL08.agg <- as.data.frame(plant.rbcL08.agg)

#### scale of aggregation = 1' grid cells ####

#load and format maps
for(i in 1:length(map.yrs)){
  fname <- paste0('hyde32_',map.yrs[i],'_1')
  path <- paste0('~/Google Drive/recherche/Intraspecific genetic diversity/data/',fname,'.RData')
  load(path)
  if(i == 1){
    get(fname) %>% filter(!is.na(pop_tot)) %>%
      mutate('cell' = paste0(.data$lat,'_',.data$long)) %>% select(-lat, -long) %>%
      mutate('mapyear' = map.yrs[i]) %>% select(mapyear,cell,everything()) -> lu.map
  } else {
    get(fname) %>% filter(!is.na(pop_tot)) %>%
      mutate('cell' = paste0(.data$lat,'_',.data$long)) %>% select(-lat, -long) %>%
      mutate('mapyear' = map.yrs[i]) %>% select(mapyear,cell,everything()) -> tmp.map
    lu.map <- bind_rows(lu.map,tmp.map)
  }
  rm(list = c(fname))
}
rm(tmp.map)

# load('~/Google Drive/recherche/Intraspecific genetic diversity/data/hyde32_2017_1.RData')
# load('~/Google Drive/recherche/Intraspecific genetic diversity/data/hyde32_1980_1.RData')
# old <- hyde32_1980_1 %>% mutate('cell' = paste0(.data$lat,'_',.data$long)) %>%
#   #rename('pop_density' = pop_tot, 'p_pasture' = grazing, 'p_urban' = urban, 'p_crop' = cropland) %>%
#   mutate('cell.r' = as.factor(cell)) %>% select(-c(lat,long,cell))
# lu1 <- hyde32_2017_1 %>% mutate('cell' = paste0(.data$lat,'_',.data$long)) %>%
#   #rename('pop_density' = pop_tot, 'p_pasture' = grazing, 'p_urban' = urban, 'p_crop' = cropland) %>%
#   mutate('cell.r' = as.factor(cell))
# lu1 <- left_join(x = lu1, y = old, by = 'cell.r' , suffix = c('.2017','.1980')) %>%
#   select(lat,long,cell,cell.r,everything())
# rm(hyde32_1980_1, hyde32_2017_1,old)

# lu1 <- read_csv('~/Google Drive/recherche/Intraspecific genetic diversity/data/landuse_1.csv') %>%
#   select(-starts_with('X'))
# lu1$cell <- gsub(' ','',lu1$cell)
# lu1 <- lu1 %>% mutate(lat = gsub("_.*","", cell), long = gsub(".*_","", cell)) %>%
#   mutate_at(vars(lat:long), funs(round(as.numeric(.),3))) %>%
#   mutate(cell.r = paste(lat,'_',long,sep=''))
# lu1$cell.r <-as.factor(lu1$cell.r)

## MAMMALS

mam1 <- read.csv('~/Google Drive/recherche/Intraspecific genetic diversity/data/mamm_pairwise_1by1_anth.csv')
colnames(mam1)[1] <- 'species.year.ID' #to be consistent with other files
colnames(mam1)[10] <- 'anth' #change new name back to old to be consistent with analysis script

# add year & species (split column)
info <- unlist(strsplit(as.character(mam1$species.year.ID), '[.]'))
mam1$year <- info[seq(2, length(info), by = 2)]
mam1$year <- as.numeric(mam1$year)
mam1$species <- info[seq(1, length(info), by = 2)]
mam1$species <- as.factor(mam1$species)
rm(info)

#mam1 <- mam1[mam1$num_per_bp != 0,] # removes all comparisons with 0 diffs = could be redundant sequences
mam1 <- na.omit(mam1) %>% filter(str_count(.$species, '_') == 1, str_count(.$species, unwanted) == 0) #removes NAs & misidentified species
mam1 <- droplevels(mam1)

#summarizing data frame
mam1.agg <- mam1 %>% group_by(species, year, cell) %>% summarize(
  div = mean(num_per_bp),
  anthrome = mean(anth),
  ncomps = length(num_per_bp),
  nseqs = length(unique(c(seq1,seq2)))) %>%
  mutate_at(vars(cell), as.character) %>%
  mutate(lat = gsub("_.*","", cell), long = gsub(".*_","", cell)) %>%
  mutate_at(vars(lat:long), funs(round(as.numeric(.),3))) %>%
  mutate(cell = paste0(lat,'_',long)) %>%
  select(-anthrome) %>% as.data.frame(.)

#checking if all cells from which we have sequences are found in land use map
if(sum(!(mam1.agg$cell %in% lu.map[lu.map$mapyear == map.yrs[1],'cell'])) > 0) {print('not all cells matched')}

#adding land use and human population density infomation
mam1.agg$mapyear <- map.yr.idx$map.yr[match(mam1.agg$year, map.yr.idx$seq.yr)]
mam1.agg <- inner_join(mam1.agg,lu.map,by=c('mapyear','cell'))


# removing groups with only one sequence comparison (poor measure of diversity)
#mam1.agg <- mam1.agg[mam1.agg$ncomps > 1,]

#creating unique identifier for population
mam1.agg$pop <- as.factor(paste(mam1.agg$species,mam1.agg$cell,sep='_'))

#creating unique identifier for species/year combination
mam1.agg$sp_yr <- as.factor(paste(mam1.agg$species,mam1.agg$year,sep='_'))

#adding length of time series
pop.dur <- aggregate(div~pop, mam1.agg, FUN=length)
mam1.agg$n.years <- pop.dur$div[match(mam1.agg$pop,pop.dur$pop)]
rm(pop.dur)

# # excluding pops in the water
# mam1.agg <- mam1.agg[mam1.agg$anthrome != 0,]
# 
# #converting anthrome to binary. natural are all 'semi-natural' and 'wild' in Ellis 2010
# mam1.agg$anth.cat <- 'anthro'
# mam1.agg$anth.cat[mam1.agg$anthrome >= 50] <- 'natural'
# mam1.agg$anth.cat <- as.factor(mam1.agg$anth.cat)
# 
# mam1.agg <- droplevels(mam1.agg)
# 
# # mam1.short <- mam1.agg[!(duplicated(mam1.agg$pop)),]
# nlevels(mam1.agg$species) #number of species
# nlevels(mam1.agg$pop) #number of populations
# # how many pops have data from same year, different anthromes?
# mam1.tab <- with(mam1.agg,table(sp_yr, anth.cat))
# mam1.tab <- mam1.tab[mam1.tab[,1] > 0 & mam1.tab[,2] > 0,]
# mam1.tx.spatial <- row.names(mam1.tab)
# rm(mam1.tab)
# #number of data points for spatial analysis within species @ .28
# length(mam1.tx.spatial)
# 
# #time series duration
# time.series1 <- data.frame('nb.year' = 1:20)
# time.series1$mams <- numeric(length(time.series1$nb.year))
# for(i in 1:dim(time.series1)[1]){
#   time.series1$mams[i] <- length(unique(mam1.agg[mam1.agg$n.years == time.series1$nb.year[i],]$pop))
# }
# sum(time.series1$mams[5:20]) #number of time series > 5 yrs

## BIRDS

aves1 <- read.csv('~/Google Drive/recherche/Intraspecific genetic diversity/data/aves_pairwise_1by1_anth.csv')
colnames(aves1)[1] <- 'species.year.ID' #to be consistent with other files
colnames(aves1)[10] <- 'anth' #change new name back to old to be consistent with analysis script

# add year & species (split column)
info <- unlist(strsplit(as.character(aves1$species.year.ID), '[.]'))
aves1$year <- info[seq(2, length(info), by = 2)]
aves1$year <- as.numeric(aves1$year)
aves1$species <- info[seq(1, length(info), by = 2)]
aves1$species <- as.factor(aves1$species)
rm(info)

#aves1 <- aves1[aves1$num_per_bp != 0,] # removes all comparisons with 0 diffs = could be redundant sequences
aves1 <- na.omit(aves1) %>% filter(str_count(.$species, '_') == 1, str_count(.$species, unwanted) == 0) #removes NAs & misidentified species
aves1 <- droplevels(aves1)

#summarizing data frame
aves1.agg <- aves1 %>% group_by(species, year, cell) %>% summarize(
  div = mean(num_per_bp),
  anthrome = mean(anth),
  ncomps = length(num_per_bp),
  nseqs = length(unique(c(seq1,seq2)))) %>%
  mutate_at(vars(cell), as.character) %>%
  mutate(lat = gsub("_.*","", cell), long = gsub(".*_","", cell)) %>%
  mutate_at(vars(lat:long), funs(round(as.numeric(.),3))) %>%
  mutate(cell = paste0(lat,'_',long)) %>%
  select(-anthrome) %>% as.data.frame(.)

#checking if all cells from which we have sequences are found in land use map
if(sum(!(aves1.agg$cell %in% lu.map[lu.map$mapyear == map.yrs[1],'cell'])) > 0) {print('not all cells matched')}

#adding land use and human population density infomation
aves1.agg$mapyear <- map.yr.idx$map.yr[match(aves1.agg$year, map.yr.idx$seq.yr)]
aves1.agg <- inner_join(aves1.agg,lu.map,by=c('mapyear','cell'))

# removing groups with only one sequence comparison (poor measure of diversity)
#aves1.agg <- aves1.agg[aves1.agg$ncomps > 1,]

#creating unique identifier for population
aves1.agg$pop <- as.factor(paste(aves1.agg$species,aves1.agg$cell,sep='_'))

#creating unique identifier for species/year combination
aves1.agg$sp_yr <- as.factor(paste(aves1.agg$species,aves1.agg$year,sep='_'))

#adding length of time series
pop.dur <- aggregate(div~pop, aves1.agg, FUN=length)
aves1.agg$n.years <- pop.dur$div[match(aves1.agg$pop,pop.dur$pop)]
rm(pop.dur)

# # excluding pops in the water
# aves1.agg <- aves1.agg[aves1.agg$anthrome != 0,]
# 
# #converting anthrome to binary
# aves1.agg$anth.cat <- 'anthro'
# aves1.agg$anth.cat[aves1.agg$anthrome >= 50] <- 'natural'
# aves1.agg$anth.cat <- as.factor(aves1.agg$anth.cat)
# 
# aves1.agg <- droplevels(aves1.agg)
# 
# # aves1.short <- aves1.agg[!(duplicated(aves1.agg$pop)),]
# nlevels(aves1.agg$species) #number of species
# nlevels(aves1.agg$pop) #number of populations
# # how many pops have data from same year, different anthromes?
# aves1.tab <- with(aves1.agg,table(sp_yr, anth.cat))
# aves1.tab <- aves1.tab[aves1.tab[,1] > 0 & aves1.tab[,2] > 0,]
# aves1.tx.spatial <- row.names(aves1.tab)
# rm(aves1.tab)
# #number of data points for spatial analysis within species @ .28
# length(aves1.tx.spatial)
# 
# #time series duration
# time.series1$aves <- numeric(length(time.series1$nb.year))
# for(i in 1:dim(time.series1)[1]){
#   time.series1$aves[i] <- length(unique(aves1.agg[aves1.agg$n.years == time.series1$nb.year[i],]$pop))
# }
# sum(time.series1$aves[5:20]) #nb of time series

## BONY FISHES

acti1 <- read.csv('~/Google Drive/recherche/Intraspecific genetic diversity/data/acti_pairwise_1by1_anth.csv')
colnames(acti1)[1] <- 'species.year.ID' #to be consistent with other files
colnames(acti1)[10] <- 'anth' #change new name back to old to be consistent with analysis script

# add year & species (split column)
info <- unlist(strsplit(as.character(acti1$species.year.ID), '[.]'))
acti1$year <- info[seq(2, length(info), by = 2)]
acti1$year <- as.numeric(acti1$year)
acti1$species <- info[seq(1, length(info), by = 2)]
acti1$species <- as.factor(acti1$species)
rm(info)

#acti1 <- acti1[acti1$num_per_bp != 0,] # removes all comparisons with 0 diffs = could be redundant sequences
acti1 <- na.omit(acti1) %>% filter(str_count(.$species, '_') == 1, str_count(.$species, unwanted) == 0) #removes NAs & misidentified species
acti1 <- droplevels(acti1)

#summarizing data frame
acti1.agg <- acti1 %>% group_by(species, year, cell) %>% summarize(
  div = mean(num_per_bp),
  anthrome = mean(anth),
  ncomps = length(num_per_bp),
  nseqs = length(unique(c(seq1,seq2)))) %>%
  mutate_at(vars(cell), as.character) %>%
  mutate(lat = gsub("_.*","", cell), long = gsub(".*_","", cell)) %>%
  mutate_at(vars(lat:long), funs(round(as.numeric(.),3))) %>%
  mutate(cell = paste0(lat,'_',long)) %>%
  select(-anthrome) %>% as.data.frame(.)

#checking if all cells from which we have sequences are found in land use map
if(sum(!(acti1.agg$cell %in% lu.map[lu.map$mapyear == map.yrs[1],'cell'])) > 0) {print('not all cells matched')}

#adding land use and human population density infomation
acti1.agg$mapyear <- map.yr.idx$map.yr[match(acti1.agg$year, map.yr.idx$seq.yr)]
acti1.agg <- inner_join(acti1.agg,lu.map,by=c('mapyear','cell'))

# removing groups with only one sequence comparison (poor measure of diversity)
#acti1.agg <- acti1.agg[acti1.agg$ncomps > 1,]

#creating unique identifier for population
acti1.agg$pop <- as.factor(paste(acti1.agg$species,acti1.agg$cell,sep='_'))

#creating unique identifier for species/year combination
acti1.agg$sp_yr <- as.factor(paste(acti1.agg$species,acti1.agg$year,sep='_'))

#adding length of time series
pop.dur <- aggregate(div~pop, acti1.agg, FUN=length)
acti1.agg$n.years <- pop.dur$div[match(acti1.agg$pop,pop.dur$pop)]
rm(pop.dur)

# # excluding pops in the ocean
# acti1.agg <- acti1.agg[acti1.agg$anthrome != 0,]
# 
# #converting anthrome to binary
# acti1.agg$anth.cat <- 'anthro'
# acti1.agg$anth.cat[acti1.agg$anthrome >= 50] <- 'natural'
# acti1.agg$anth.cat <- as.factor(acti1.agg$anth.cat)
# 
# acti1.agg <- droplevels(acti1.agg)
# 
# # acti1.short <- acti1.agg[!(duplicated(acti1.agg$pop)),]
# nlevels(acti1.agg$species) #number of species
# nlevels(acti1.agg$pop) #number of populations
# # how many pops have data from same year, different anthromes?
# acti1.tab <- with(acti1.agg,table(sp_yr, anth.cat))
# acti1.tab <- acti1.tab[acti1.tab[,1] > 0 & acti1.tab[,2] > 0,]
# acti1.tx.spatial <- row.names(acti1.tab)
# rm(acti1.tab)
# #number of data points for spatial analysis within species @ .28
# length(acti1.tx.spatial)
# 
# #time series duration
# time.series1$acti <- numeric(length(time.series1$nb.year))
# for(i in 1:dim(time.series1)[1]){
#   time.series1$acti[i] <- length(unique(acti1.agg[acti1.agg$n.years == time.series1$nb.year[i],]$pop))
# }
# sum(time.series1$acti[5:20]) #nb of time series

## INSECTS

insect1 <- read.csv('~/Google Drive/recherche/Intraspecific genetic diversity/data/insect_pairwise_1by1_anth.csv')
colnames(insect1)[1] <- 'species.year.ID' #to be consistent with other files
colnames(insect1)[10] <- 'anth' #change new name back to old to be consistent with analysis script

insect1 <- insect1[insect1$overlap > 0.5,] # filter comparisons with less than 50% overlap as in Miraldo paper (for other taxa, this step was done in Excel)
row2rm <- which(insect1$species.year.ID == 'Lamprospilus_aff..2015')
if (length(row2rm) > 0){
  insect1 <- insect1[-row2rm,] #removes database entry with incorrect species name
}

#insect1 <- insect1[insect1$num_per_bp != 0,] # removes all comparisons with 0 diffs = could be redundant sequences
insect1 <- na.omit(insect1) #removes NAs

# add year & species (split column)
info <- unlist(strsplit(as.character(insect1$species.year.ID), '[.]'))
insect1$year <- info[seq(2, length(info), by = 2)]
insect1$year <- as.numeric(insect1$year)
insect1$species <- info[seq(1, length(info), by = 2)]
insect1$species <- as.factor(insect1$species)
rm(info)

insect1 <- insect1 %>% filter(str_count(.$species, '_') == 1, str_count(.$species, unwanted) == 0)

insect1 <- droplevels(insect1)

#summarizing data frame
insect1.agg <- insect1 %>% group_by(species, year, cell) %>% summarize(
  div = mean(num_per_bp),
  anthrome = mean(anth),
  ncomps = length(num_per_bp),
  nseqs = length(unique(c(seq1,seq2)))) %>%
  mutate_at(vars(cell), as.character) %>%
  mutate(lat = gsub("_.*","", cell), long = gsub(".*_","", cell)) %>%
  mutate_at(vars(lat:long), funs(round(as.numeric(.),3))) %>%
  mutate(cell = paste0(lat,'_',long)) %>%
  select(-anthrome) %>% as.data.frame(.)

#checking if all cells from which we have sequences are found in land use map
if(sum(!(insect1.agg$cell %in% lu.map[lu.map$mapyear == map.yrs[1],'cell'])) > 0) {print('not all cells matched')}

#adding land use and human population density infomation
insect1.agg$mapyear <- map.yr.idx$map.yr[match(insect1.agg$year, map.yr.idx$seq.yr)]
insect1.agg <- inner_join(insect1.agg,lu.map,by=c('mapyear','cell'))

# removing groups with only one sequence comparison (poor measure of diversity)
#insect1.agg <- insect1.agg[insect1.agg$ncomps > 1,]

#creating unique identifier for population
insect1.agg$pop <- as.factor(paste(insect1.agg$species,insect1.agg$cell,sep='_'))

#creating unique identifier for species/year combination
insect1.agg$sp_yr <- as.factor(paste(insect1.agg$species,insect1.agg$year,sep='_'))

#adding length of time series
pop.dur <- aggregate(div~pop, insect1.agg, FUN=length)
insect1.agg$n.years <- pop.dur$div[match(insect1.agg$pop,pop.dur$pop)]
rm(pop.dur)

# # excluding pops in the ocean
# insect1.agg <- insect1.agg[insect1.agg$anthrome != 0,]
# 
# #converting anthrome to binary
# insect1.agg$anth.cat <- 'anthro'
# insect1.agg$anth.cat[insect1.agg$anthrome >= 50] <- 'natural'
# insect1.agg$anth.cat <- as.factor(insect1.agg$anth.cat)
# 
# insect1.agg <- droplevels(insect1.agg)
# 
# # insect1.short <- insect1.agg[!(duplicated(insect1.agg$pop)),]
# nlevels(insect1.agg$species) #number of species
# nlevels(insect1.agg$pop) #number of populations
# # how many pops have data from same year, different anthromes?
# insect1.tab <- with(insect1.agg,table(sp_yr, anth.cat))
# insect1.tab <- insect1.tab[insect1.tab[,1] > 0 & insect1.tab[,2] > 0,]
# insect1.tx.spatial <- row.names(insect1.tab)
# rm(insect1.tab)
# #number of data points for spatial analysis within species @ .28
# length(insect1.tx.spatial)
# 
# #time series duration
# time.series1$insects <- numeric(length(time.series1$nb.year))
# for(i in 1:dim(time.series1)[1]){
#   time.series1$insects[i] <- length(unique(insect1.agg[insect1.agg$n.years == time.series1$nb.year[i],]$pop))
# }
# sum(time.series1$insects[5:20]) #nb of time series

## PLANTS (ITS2)

plant.its1 <- read.csv('~/Google Drive/recherche/Intraspecific genetic diversity/data/plants_ITS_pairwise_1by1_anth.csv')
colnames(plant.its1)[1] <- 'species.year.ID' #to be consistent with other files
colnames(plant.its1)[10] <- 'anth' #change new name back to old to be consistent with analysis script

# add year & species (split column)
info <- unlist(strsplit(as.character(plant.its1$species.year.ID), '[.]'))
plant.its1$year <- info[seq(2, length(info), by = 2)]
plant.its1$year <- as.numeric(plant.its1$year)
plant.its1$species <- info[seq(1, length(info), by = 2)]
plant.its1$species <- as.factor(plant.its1$species)
rm(info)

#plant.its1 <- plant.its1[plant.its1$num_per_bp != 0,] # removes all comparisons with 0 diffs = could be redundant sequences
plant.its1 <- na.omit(plant.its1) %>% filter(str_count(.$species, '_') == 1, str_count(.$species, unwanted) == 0) #removes NAs & misidentified species
plant.its1 <- droplevels(plant.its1)

#summarizing data frame
plant.its1.agg <- plant.its1 %>% group_by(species, year, cell) %>% summarize(
  div = mean(num_per_bp),
  anthrome = mean(anth),
  ncomps = length(num_per_bp),
  nseqs = length(unique(c(seq1,seq2)))) %>%
  mutate_at(vars(cell), as.character) %>%
  mutate(lat = gsub("_.*","", cell), long = gsub(".*_","", cell)) %>%
  mutate_at(vars(lat:long), funs(round(as.numeric(.),3))) %>%
  mutate(cell = paste0(lat,'_',long)) %>%
  select(-anthrome) %>% as.data.frame(.)

#checking if all cells from which we have sequences are found in land use map
if(sum(!(plant.its1.agg$cell %in% lu.map[lu.map$mapyear == map.yrs[1],'cell'])) > 0) {print('not all cells matched')}

#adding land use and human population density infomation
plant.its1.agg$mapyear <- map.yr.idx$map.yr[match(plant.its1.agg$year, map.yr.idx$seq.yr)]
plant.its1.agg <- inner_join(plant.its1.agg,lu.map,by=c('mapyear','cell'))

# removing groups with only one sequence comparison (poor measure of diversity)
#plant.its1.agg <- plant.its1.agg[plant.its1.agg$ncomps > 1,]

#creating unique identifier for population
plant.its1.agg$pop <- as.factor(paste(plant.its1.agg$species,plant.its1.agg$cell,sep='_'))

#creating unique identifier for species/year combination
plant.its1.agg$sp_yr <- as.factor(paste(plant.its1.agg$species,plant.its1.agg$year,sep='_'))

#adding length of time series
pop.dur <- aggregate(div~pop, plant.its1.agg, FUN=length)
plant.its1.agg$n.years <- pop.dur$div[match(plant.its1.agg$pop,pop.dur$pop)]
rm(pop.dur)

# # excluding pops in the water
# plant.its1.agg <- plant.its1.agg[plant.its1.agg$anthrome != 0,]
# 
# #converting anthrome to binary
# plant.its1.agg$anth.cat <- 'anthro'
# plant.its1.agg$anth.cat[plant.its1.agg$anthrome >= 50] <- 'natural'
# plant.its1.agg$anth.cat <- as.factor(plant.its1.agg$anth.cat)
# 
# plant.its1.agg <- droplevels(plant.its1.agg)
# 
# # plant.its1.short <- plant.its1.agg[!(duplicated(plant.its1.agg$pop)),]
# nlevels(plant.its1.agg$species) #number of species
# nlevels(plant.its1.agg$pop) #number of populations
# # how many pops have data from same year, different anthromes?
# plant.its1.tab <- with(plant.its1.agg,table(sp_yr, anth.cat))
# plant.its1.tab <- plant.its1.tab[plant.its1.tab[,1] > 0 & plant.its1.tab[,2] > 0,]
# plant.its1.tx.spatial <- row.names(plant.its1.tab)
# rm(plant.its1.tab)
# #number of data points for spatial analysis within species @ .28
# length(plant.its1.tx.spatial)
# 
# #time series duration
# time.series1$plant.its <- numeric(length(time.series1$nb.year))
# for(i in 1:dim(time.series1)[1]){
#   time.series1$plant.its[i] <- length(unique(plant.its1.agg[plant.its1.agg$n.years == time.series1$nb.year[i],]$pop))
# }
# sum(time.series1$plant.its[5:20]) #nb of time series

## PLANTS (MATK)

plant.matK1 <- read.csv('~/Google Drive/recherche/Intraspecific genetic diversity/data/plants_pairwise_matK_1by1_anth.csv')
colnames(plant.matK1)[1] <- 'species.year.ID' #to be consistent with other files
colnames(plant.matK1)[10] <- 'anth' #change new name back to old to be consistent with analysis script

#removing line with double period
#sapply(gregexpr("\\.", as.character(plant.matK1$species.year.ID)), tail, 1) == sapply(gregexpr("\\.", as.character(plant.matK1$species.year.ID)), head, 1)
row2rm <- which(plant.matK1$species.year.ID == 'Crataegus_indet..2011')
if (length(row2rm) > 0){
  plant.matK1 <- plant.matK1[-row2rm,] #removes database entry with incorrect species name
}

# add year & species (split column)
info <- unlist(strsplit(as.character(plant.matK1$species.year.ID), '[.]'))
plant.matK1$year <- info[seq(2, length(info), by = 2)]
plant.matK1$year <- as.numeric(plant.matK1$year)
plant.matK1$species <- info[seq(1, length(info), by = 2)]
plant.matK1$species <- as.factor(plant.matK1$species)
rm(info)

#plant.matK1 <- plant.matK1[plant.matK1$num_per_bp != 0,] # removes all comparisons with 0 diffs = could be redundant sequences
plant.matK1 <- na.omit(plant.matK1) %>% filter(str_count(.$species, '_') == 1, str_count(.$species, unwanted) == 0) #removes NAs & misidentified species
plant.matK1 <- droplevels(plant.matK1)

#summarizing data frame
plant.matK1.agg <- plant.matK1 %>% group_by(species, year, cell) %>% summarize(
  div = mean(num_per_bp),
  anthrome = mean(anth),
  ncomps = length(num_per_bp),
  nseqs = length(unique(c(seq1,seq2)))) %>%
  mutate_at(vars(cell), as.character) %>%
  mutate(lat = gsub("_.*","", cell), long = gsub(".*_","", cell)) %>%
  mutate_at(vars(lat:long), funs(round(as.numeric(.),3))) %>%
  mutate(cell = paste0(lat,'_',long)) %>%
  select(-anthrome) %>% as.data.frame(.)

#checking if all cells from which we have sequences are found in land use map
if(sum(!(plant.matK1.agg$cell %in% lu.map[lu.map$mapyear == map.yrs[1],'cell'])) > 0) {print('not all cells matched')}

#adding land use and human population density infomation
plant.matK1.agg$mapyear <- map.yr.idx$map.yr[match(plant.matK1.agg$year, map.yr.idx$seq.yr)]
plant.matK1.agg <- inner_join(plant.matK1.agg,lu.map,by=c('mapyear','cell'))

# removing groups with only one sequence comparison (poor measure of diversity)
#plant.matK1.agg <- plant.matK1.agg[plant.matK1.agg$ncomps > 1,]

#creating unique identifier for population
plant.matK1.agg$pop <- as.factor(paste(plant.matK1.agg$species,plant.matK1.agg$cell,sep='_'))

#creating unique identifier for species/year combination
plant.matK1.agg$sp_yr <- as.factor(paste(plant.matK1.agg$species,plant.matK1.agg$year,sep='_'))

#adding length of time series
pop.dur <- aggregate(div~pop, plant.matK1.agg, FUN=length)
plant.matK1.agg$n.years <- pop.dur$div[match(plant.matK1.agg$pop,pop.dur$pop)]
rm(pop.dur)

# # excluding pops in the water
# plant.matK1.agg <- plant.matK1.agg[plant.matK1.agg$anthrome != 0,]
# 
# #converting anthrome to binary
# plant.matK1.agg$anth.cat <- 'anthro'
# plant.matK1.agg$anth.cat[plant.matK1.agg$anthrome >= 50] <- 'natural'
# plant.matK1.agg$anth.cat <- as.factor(plant.matK1.agg$anth.cat)
# 
# plant.matK1.agg <- droplevels(plant.matK1.agg)
# 
# # plant.matK1.short <- plant.matK1.agg[!(duplicated(plant.matK1.agg$pop)),]
# nlevels(plant.matK1.agg$species) #number of species
# nlevels(plant.matK1.agg$pop) #number of populations
# # how many pops have data from same year, different anthromes?
# plant.matK1.tab <- with(plant.matK1.agg,table(sp_yr, anth.cat))
# plant.matK1.tab <- plant.matK1.tab[plant.matK1.tab[,1] > 0 & plant.matK1.tab[,2] > 0,]
# plant.matK1.tx.spatial <- row.names(plant.matK1.tab)
# rm(plant.matK1.tab)
# #number of data points for spatial analysis within species @ .28
# length(plant.matK1.tx.spatial)
# 
# #time series duration
# time.series1$plant.matK <- numeric(length(time.series1$nb.year))
# for(i in 1:dim(time.series1)[1]){
#   time.series1$plant.matK[i] <- length(unique(plant.matK1.agg[plant.matK1.agg$n.years == time.series1$nb.year[i],]$pop))
# }
# sum(time.series1$plant.matK[5:20]) #nb of time series

## PLANTS (rbcL)

plant.rbcL1 <- read.csv('~/Google Drive/recherche/Intraspecific genetic diversity/data/plants_pairwise_rbcL_1by1_anth.csv')
colnames(plant.rbcL1)[1] <- 'species.year.ID' #to be consistent with other files
colnames(plant.rbcL1)[10] <- 'anth' #change new name back to old to be consistent with analysis script

#which(sapply(gregexpr("\\.", as.character(plant.rbcL1$species.year.ID)), tail, 1) != sapply(gregexpr("\\.", as.character(plant.rbcL1$species.year.ID)), head, 1))
row2rm <- which(plant.rbcL1$species.year.ID == 'Crataegus_indet..2011')
if (length(row2rm) > 0){
  plant.rbcL1 <- plant.rbcL1[-row2rm,] #removes database entry with incorrect species name
}

# add year & species (split column)
info <- unlist(strsplit(as.character(plant.rbcL1$species.year.ID), '[.]'))
plant.rbcL1$year <- info[seq(2, length(info), by = 2)]
plant.rbcL1$year <- as.numeric(plant.rbcL1$year)
plant.rbcL1$species <- info[seq(1, length(info), by = 2)]
plant.rbcL1$species <- as.factor(plant.rbcL1$species)
rm(info)

#plant.rbcL1 <- plant.rbcL1[plant.rbcL1$num_per_bp != 0,] # removes all comparisons with 0 diffs = could be redundant sequences
plant.rbcL1 <- na.omit(plant.rbcL1) %>% filter(str_count(.$species, '_') == 1, str_count(.$species, unwanted) == 0) #removes NAs & misidentified species
plant.rbcL1 <- droplevels(plant.rbcL1)

#summarizing data frame
plant.rbcL1.agg <- plant.rbcL1 %>% group_by(species, year, cell) %>% summarize(
  div = mean(num_per_bp),
  anthrome = mean(anth),
  ncomps = length(num_per_bp),
  nseqs = length(unique(c(seq1,seq2)))) %>%
  mutate_at(vars(cell), as.character) %>%
  mutate(lat = gsub("_.*","", cell), long = gsub(".*_","", cell)) %>%
  mutate_at(vars(lat:long), funs(round(as.numeric(.),3))) %>%
  mutate(cell = paste0(lat,'_',long)) %>%
  select(-anthrome) %>% as.data.frame(.)

#checking if all cells from which we have sequences are found in land use map
if(sum(!(plant.rbcL1.agg$cell %in% lu.map[lu.map$mapyear == map.yrs[1],'cell'])) > 0) {print('not all cells matched')}

#adding land use and human population density infomation
plant.rbcL1.agg$mapyear <- map.yr.idx$map.yr[match(plant.rbcL1.agg$year, map.yr.idx$seq.yr)]
plant.rbcL1.agg <- inner_join(plant.rbcL1.agg,lu.map,by=c('mapyear','cell'))

# removing groups with only one sequence comparison (poor measure of diversity)
#plant.rbcL1.agg <- plant.rbcL1.agg[plant.rbcL1.agg$ncomps > 1,]

#creating unique identifier for population
plant.rbcL1.agg$pop <- as.factor(paste(plant.rbcL1.agg$species,plant.rbcL1.agg$cell,sep='_'))

#creating unique identifier for species/year combination
plant.rbcL1.agg$sp_yr <- as.factor(paste(plant.rbcL1.agg$species,plant.rbcL1.agg$year,sep='_'))

#adding length of time series
pop.dur <- aggregate(div~pop, plant.rbcL1.agg, FUN=length)
plant.rbcL1.agg$n.years <- pop.dur$div[match(plant.rbcL1.agg$pop,pop.dur$pop)]
rm(pop.dur)

# # excluding pops in the water
# plant.rbcL1.agg <- plant.rbcL1.agg[plant.rbcL1.agg$anthrome != 0,]
# 
# #converting anthrome to binary
# plant.rbcL1.agg$anth.cat <- 'anthro'
# plant.rbcL1.agg$anth.cat[plant.rbcL1.agg$anthrome >= 50] <- 'natural'
# plant.rbcL1.agg$anth.cat <- as.factor(plant.rbcL1.agg$anth.cat)
# 
# plant.rbcL1.agg <- droplevels(plant.rbcL1.agg)
# 
# # plant.rbcL1.short <- plant.rbcL1.agg[!(duplicated(plant.rbcL1.agg$pop)),]
# nlevels(plant.rbcL1.agg$species) #number of species
# nlevels(plant.rbcL1.agg$pop) #number of populations
# # how many pops have data from same year, different anthromes?
# plant.rbcL1.tab <- with(plant.rbcL1.agg,table(sp_yr, anth.cat))
# plant.rbcL1.tab <- plant.rbcL1.tab[plant.rbcL1.tab[,1] > 0 & plant.rbcL1.tab[,2] > 0,]
# plant.rbcL1.tx.spatial <- row.names(plant.rbcL1.tab)
# rm(plant.rbcL1.tab)
# #number of data points for spatial analysis within species @ .28
# length(plant.rbcL1.tx.spatial)
# 
# #time series duration
# time.series1$plant.rbcL <- numeric(length(time.series1$nb.year))
# for(i in 1:dim(time.series1)[1]){
#   time.series1$plant.rbcL[i] <- length(unique(plant.rbcL1.agg[plant.rbcL1.agg$n.years == time.series1$nb.year[i],]$pop))
# }
# sum(time.series1$plant.rbcL[5:20]) #nb of time series

#gather data summary statistics
summ1 <- data.frame('taxon'=c('mammals','birds','fish','insects','plants_its','plants_matK','plants_rbcL'),
                     'nb.sequences' = numeric(7),
                     'nb.species' = numeric(7),
                     'nb.pops' = numeric(7))
summ1[1,2:4] <- c(sum(mam1.agg$nseqs),length(levels(mam1.agg$species)),length(levels(mam1.agg$pop)))
summ1[2,2:4] <- c(sum(aves1.agg$nseqs),length(levels(aves1.agg$species)),length(levels(aves1.agg$pop)))
summ1[3,2:4] <- c(sum(acti1.agg$nseqs),length(levels(acti1.agg$species)),length(levels(acti1.agg$pop)))
summ1[4,2:4] <- c(sum(insect1.agg$nseqs),length(levels(insect1.agg$species)),length(levels(insect1.agg$pop)))
summ1[5,2:4] <- c(sum(plant.its1.agg$nseqs),length(levels(plant.its1.agg$species)),length(levels(plant.its1.agg$pop)))
summ1[6,2:4] <- c(sum(plant.matK1.agg$nseqs),length(levels(plant.matK1.agg$species)),length(levels(plant.matK1.agg$pop)))
summ1[7,2:4] <- c(sum(plant.rbcL1.agg$nseqs),length(levels(plant.rbcL1.agg$species)),length(levels(plant.rbcL1.agg$pop)))

#cleanup
rm(i,row2rm,acti1,aves1,insect1,mam1,plant.its1,plant.matK1,plant.rbcL1,lu.map)

mam1.agg <- as.data.frame(mam1.agg)
aves1.agg <- as.data.frame(aves1.agg)
acti1.agg <- as.data.frame(acti1.agg)
insect1.agg <- as.data.frame(insect1.agg)
plant.its1.agg <- as.data.frame(plant.its1.agg)
plant.matK1.agg <- as.data.frame(plant.matK1.agg)
plant.rbcL1.agg <- as.data.frame(plant.rbcL1.agg)

#### scale of aggregation = 2' grid cells ####

#load and format maps
for(i in 1:length(map.yrs)){
  fname <- paste0('hyde32_',map.yrs[i],'_2')
  path <- paste0('~/Google Drive/recherche/Intraspecific genetic diversity/data/',fname,'.RData')
  load(path)
  if(i == 1){
    get(fname) %>% filter(!is.na(pop_tot)) %>%
      mutate('cell' = paste0(.data$lat,'_',.data$long)) %>% select(-lat, -long) %>%
      mutate('mapyear' = map.yrs[i]) %>% select(mapyear,cell,everything()) -> lu.map
  } else {
    get(fname) %>% filter(!is.na(pop_tot)) %>%
      mutate('cell' = paste0(.data$lat,'_',.data$long)) %>% select(-lat, -long) %>%
      mutate('mapyear' = map.yrs[i]) %>% select(mapyear,cell,everything()) -> tmp.map
    lu.map <- bind_rows(lu.map,tmp.map)
  }
  rm(list = c(fname))
}
rm(tmp.map)

# load('~/Google Drive/recherche/Intraspecific genetic diversity/data/hyde32_2017_2.RData')
# load('~/Google Drive/recherche/Intraspecific genetic diversity/data/hyde32_1980_2.RData')
# old <- hyde32_1980_2 %>% mutate('cell' = paste0(.data$lat,'_',.data$long)) %>%
#   #rename('pop_density' = pop_tot, 'p_pasture' = grazing, 'p_urban' = urban, 'p_crop' = cropland) %>%
#   mutate('cell.r' = as.factor(cell)) %>% select(-c(lat,long,cell))
# lu2 <- hyde32_2017_2 %>% mutate('cell' = paste0(.data$lat,'_',.data$long)) %>%
#   #rename('pop_density' = pop_tot, 'p_pasture' = grazing, 'p_urban' = urban, 'p_crop' = cropland) %>%
#   mutate('cell.r' = as.factor(cell))
# lu2 <- left_join(x = lu2, y = old, by = 'cell.r' , suffix = c('.2017','.1980')) %>%
#   select(lat,long,cell,cell.r,everything())
# rm(hyde32_1980_2, hyde32_2017_2,old)

# lu2 <- read_csv('~/Google Drive/recherche/Intraspecific genetic diversity/data/landuse_2.csv') %>%
#   select(-starts_with('X'))
# lu2$cell <- gsub(' ','',lu2$cell)
# lu2 <- lu2 %>% mutate(lat = gsub("_.*","", cell), long = gsub(".*_","", cell)) %>%
#   mutate_at(vars(lat:long), funs(round(as.numeric(.),3))) %>%
#   mutate(cell.r = paste(lat,'_',long,sep=''))
# lu2$cell.r <-as.factor(lu2$cell.r)

## MAMMALS

mam2 <- read.csv('~/Google Drive/recherche/Intraspecific genetic diversity/data/mamm_pairwise_2by2_anth.csv')
colnames(mam2)[1] <- 'species.year.ID' #to be consistent with other files
colnames(mam2)[10] <- 'anth' #change new name back to old to be consistent with analysis script

# add year & species (split column)
info <- unlist(strsplit(as.character(mam2$species.year.ID), '[.]'))
mam2$year <- info[seq(2, length(info), by = 2)]
mam2$year <- as.numeric(mam2$year)
mam2$species <- info[seq(1, length(info), by = 2)]
mam2$species <- as.factor(mam2$species)
rm(info)

#mam2 <- mam2[mam2$num_per_bp != 0,] # removes all comparisons with 0 diffs = could be redundant sequences
mam2 <- na.omit(mam2) %>% filter(str_count(.$species, '_') == 1, str_count(.$species, unwanted) == 0) #removes NAs & misidentified species
mam2 <- droplevels(mam2)

#summarizing data frame
mam2.agg <- mam2 %>% group_by(species, year, cell) %>% summarize(
  div = mean(num_per_bp),
  anthrome = mean(anth),
  ncomps = length(num_per_bp),
  nseqs = length(unique(c(seq1,seq2)))) %>%
  mutate_at(vars(cell), as.character) %>%
  mutate(lat = gsub("_.*","", cell), long = gsub(".*_","", cell)) %>%
  mutate_at(vars(lat:long), funs(round(as.numeric(.),3))) %>%
  mutate(cell = paste0(lat,'_',long)) %>%
  select(-anthrome) %>% as.data.frame(.)

#checking if all cells from which we have sequences are found in land use map
if(sum(!(mam2.agg$cell %in% lu.map[lu.map$mapyear == map.yrs[1],'cell'])) > 0) {print('not all cells matched')}

#adding land use and human population density infomation
mam2.agg$mapyear <- map.yr.idx$map.yr[match(mam2.agg$year, map.yr.idx$seq.yr)]
mam2.agg <- inner_join(mam2.agg,lu.map,by=c('mapyear','cell'))

# removing groups with only one sequence comparison (poor measure of diversity)
#mam2.agg <- mam2.agg[mam2.agg$ncomps > 1,]

#creating unique identifier for population
mam2.agg$pop <- as.factor(paste(mam2.agg$species,mam2.agg$cell,sep='_'))

#creating unique identifier for species/year combination
mam2.agg$sp_yr <- as.factor(paste(mam2.agg$species,mam2.agg$year,sep='_'))

#adding length of time series
pop.dur <- aggregate(div~pop, mam2.agg, FUN=length)
mam2.agg$n.years <- pop.dur$div[match(mam2.agg$pop,pop.dur$pop)]
rm(pop.dur)

# excluding pops in the water
# mam2.agg <- mam2.agg[mam2.agg$anthrome != 0,]
# 
# #converting anthrome to binary. natural are all 'semi-natural' and 'wild' in Ellis 2010
# mam2.agg$anth.cat <- 'anthro'
# mam2.agg$anth.cat[mam2.agg$anthrome >= 50] <- 'natural'
# mam2.agg$anth.cat <- as.factor(mam2.agg$anth.cat)
# 
# mam2.agg <- droplevels(mam2.agg)
# 
# # mam2.short <- mam2.agg[!(duplicated(mam2.agg$pop)),]
# nlevels(mam2.agg$species) #number of species
# nlevels(mam2.agg$pop) #number of populations
# # how many pops have data from same year, different anthromes?
# mam2.tab <- with(mam2.agg,table(sp_yr, anth.cat))
# mam2.tab <- mam2.tab[mam2.tab[,1] > 0 & mam2.tab[,2] > 0,]
# mam2.tx.spatial <- row.names(mam2.tab)
# rm(mam2.tab)
# #number of data points for spatial analysis within species @ .28
# length(mam2.tx.spatial)
# 
# #time series duration
# time.series2 <- data.frame('nb.year' = 1:20)
# time.series2$mams <- numeric(length(time.series2$nb.year))
# for(i in 1:dim(time.series2)[1]){
#   time.series2$mams[i] <- length(unique(mam2.agg[mam2.agg$n.years == time.series2$nb.year[i],]$pop))
# }
# sum(time.series2$mams[5:20]) #number of time series > 5 yrs

## BIRDS

aves2 <- read.csv('~/Google Drive/recherche/Intraspecific genetic diversity/data/aves_pairwise_2by2_anth.csv')
colnames(aves2)[1] <- 'species.year.ID' #to be consistent with other files
colnames(aves2)[10] <- 'anth' #change new name back to old to be consistent with analysis script

# add year & species (split column)
info <- unlist(strsplit(as.character(aves2$species.year.ID), '[.]'))
aves2$year <- info[seq(2, length(info), by = 2)]
aves2$year <- as.numeric(aves2$year)
aves2$species <- info[seq(1, length(info), by = 2)]
aves2$species <- as.factor(aves2$species)
rm(info)

#aves2 <- aves2[aves2$num_per_bp != 0,] # removes all comparisons with 0 diffs = could be redundant sequences
aves2 <- na.omit(aves2) %>% filter(str_count(.$species, '_') == 1, str_count(.$species, unwanted) == 0) #removes NAs & misidentified species
aves2 <- droplevels(aves2)

#summarizing data frame
aves2.agg <- aves2 %>% group_by(species, year, cell) %>% summarize(
  div = mean(num_per_bp),
  anthrome = mean(anth),
  ncomps = length(num_per_bp),
  nseqs = length(unique(c(seq1,seq2)))) %>%
  mutate_at(vars(cell), as.character) %>%
  mutate(lat = gsub("_.*","", cell), long = gsub(".*_","", cell)) %>%
  mutate_at(vars(lat:long), funs(round(as.numeric(.),3))) %>%
  mutate(cell = paste0(lat,'_',long)) %>%
  select(-anthrome) %>% as.data.frame(.)

#checking if all cells from which we have sequences are found in land use map
if(sum(!(aves2.agg$cell %in% lu.map[lu.map$mapyear == map.yrs[1],'cell'])) > 0) {print('not all cells matched')}

#adding land use and human population density infomation
aves2.agg$mapyear <- map.yr.idx$map.yr[match(aves2.agg$year, map.yr.idx$seq.yr)]
aves2.agg <- inner_join(aves2.agg,lu.map,by=c('mapyear','cell'))


# removing groups with only one sequence comparison (poor measure of diversity)
#aves2.agg <- aves2.agg[aves2.agg$ncomps > 1,]

#creating unique identifier for population
aves2.agg$pop <- as.factor(paste(aves2.agg$species,aves2.agg$cell,sep='_'))

#creating unique identifier for species/year combination
aves2.agg$sp_yr <- as.factor(paste(aves2.agg$species,aves2.agg$year,sep='_'))

#adding length of time series
pop.dur <- aggregate(div~pop, aves2.agg, FUN=length)
aves2.agg$n.years <- pop.dur$div[match(aves2.agg$pop,pop.dur$pop)]
rm(pop.dur)

# # excluding pops in the water
# aves2.agg <- aves2.agg[aves2.agg$anthrome != 0,]
# 
# #converting anthrome to binary
# aves2.agg$anth.cat <- 'anthro'
# aves2.agg$anth.cat[aves2.agg$anthrome >= 50] <- 'natural'
# aves2.agg$anth.cat <- as.factor(aves2.agg$anth.cat)
# 
# aves2.agg <- droplevels(aves2.agg)
# 
# # aves2.short <- aves2.agg[!(duplicated(aves2.agg$pop)),]
# nlevels(aves2.agg$species) #number of species
# nlevels(aves2.agg$pop) #number of populations
# # how many pops have data from same year, different anthromes?
# aves2.tab <- with(aves2.agg,table(sp_yr, anth.cat))
# aves2.tab <- aves2.tab[aves2.tab[,1] > 0 & aves2.tab[,2] > 0,]
# aves2.tx.spatial <- row.names(aves2.tab)
# rm(aves2.tab)
# #number of data points for spatial analysis within species @ .28
# length(aves2.tx.spatial)
# 
# #time series duration
# time.series2$aves <- numeric(length(time.series2$nb.year))
# for(i in 1:dim(time.series2)[1]){
#   time.series2$aves[i] <- length(unique(aves2.agg[aves2.agg$n.years == time.series2$nb.year[i],]$pop))
# }
# sum(time.series2$aves[5:20]) #nb of time series

## BONY FISHES

acti2 <- read.csv('~/Google Drive/recherche/Intraspecific genetic diversity/data/acti_pairwise_2by2_anth.csv')
colnames(acti2)[1] <- 'species.year.ID' #to be consistent with other files
colnames(acti2)[10] <- 'anth' #change new name back to old to be consistent with analysis script

# add year & species (split column)
info <- unlist(strsplit(as.character(acti2$species.year.ID), '[.]'))
acti2$year <- info[seq(2, length(info), by = 2)]
acti2$year <- as.numeric(acti2$year)
acti2$species <- info[seq(1, length(info), by = 2)]
acti2$species <- as.factor(acti2$species)
rm(info)

#acti2 <- acti2[acti2$num_per_bp != 0,] # removes all comparisons with 0 diffs = could be redundant sequences
acti2 <- na.omit(acti2) %>% filter(str_count(.$species, '_') == 1, str_count(.$species, unwanted) == 0) #removes NAs & misidentified species
acti2 <- droplevels(acti2)

#summarizing data frame
acti2.agg <- acti2 %>% group_by(species, year, cell) %>% summarize(
  div = mean(num_per_bp),
  anthrome = mean(anth),
  ncomps = length(num_per_bp),
  nseqs = length(unique(c(seq1,seq2)))) %>%
  mutate_at(vars(cell), as.character) %>%
  mutate(lat = gsub("_.*","", cell), long = gsub(".*_","", cell)) %>%
  mutate_at(vars(lat:long), funs(round(as.numeric(.),3))) %>%
  mutate(cell = paste0(lat,'_',long)) %>%
  select(-anthrome) %>% as.data.frame(.)

#checking if all cells from which we have sequences are found in land use map
if(sum(!(acti2.agg$cell %in% lu.map[lu.map$mapyear == map.yrs[1],'cell'])) > 0) {print('not all cells matched')}

#adding land use and human population density infomation
acti2.agg$mapyear <- map.yr.idx$map.yr[match(acti2.agg$year, map.yr.idx$seq.yr)]
acti2.agg <- inner_join(acti2.agg,lu.map,by=c('mapyear','cell'))



# removing groups with only one sequence comparison (poor measure of diversity)
#acti2.agg <- acti2.agg[acti2.agg$ncomps > 1,]

#creating unique identifier for population
acti2.agg$pop <- as.factor(paste(acti2.agg$species,acti2.agg$cell,sep='_'))

#creating unique identifier for species/year combination
acti2.agg$sp_yr <- as.factor(paste(acti2.agg$species,acti2.agg$year,sep='_'))

#adding length of time series
pop.dur <- aggregate(div~pop, acti2.agg, FUN=length)
acti2.agg$n.years <- pop.dur$div[match(acti2.agg$pop,pop.dur$pop)]
rm(pop.dur)

# # excluding pops in the ocean
# acti2.agg <- acti2.agg[acti2.agg$anthrome != 0,]
# 
# #converting anthrome to binary
# acti2.agg$anth.cat <- 'anthro'
# acti2.agg$anth.cat[acti2.agg$anthrome >= 50] <- 'natural'
# acti2.agg$anth.cat <- as.factor(acti2.agg$anth.cat)
# 
# acti2.agg <- droplevels(acti2.agg)
# 
# # acti2.short <- acti2.agg[!(duplicated(acti2.agg$pop)),]
# nlevels(acti2.agg$species) #number of species
# nlevels(acti2.agg$pop) #number of populations
# # how many pops have data from same year, different anthromes?
# acti2.tab <- with(acti2.agg,table(sp_yr, anth.cat))
# acti2.tab <- acti2.tab[acti2.tab[,1] > 0 & acti2.tab[,2] > 0,]
# acti2.tx.spatial <- row.names(acti2.tab)
# rm(acti2.tab)
# #number of data points for spatial analysis within species @ .28
# length(acti2.tx.spatial)
# 
# #time series duration
# time.series2$acti <- numeric(length(time.series2$nb.year))
# for(i in 1:dim(time.series2)[1]){
#   time.series2$acti[i] <- length(unique(acti2.agg[acti2.agg$n.years == time.series2$nb.year[i],]$pop))
# }
# sum(time.series2$acti[5:20]) #nb of time series

## INSECTS

insect2 <- read.csv('~/Google Drive/recherche/Intraspecific genetic diversity/data/insect_pairwise_2by2_anth.csv')
colnames(insect2)[1] <- 'species.year.ID' #to be consistent with other files
colnames(insect2)[10] <- 'anth' #change new name back to old to be consistent with analysis script

insect2 <- insect2[insect2$overlap > 0.5,] # filter comparisons with less than 50% overlap as in Miraldo paper (for other taxa, this step was done in Excel)
row2rm <- which(insect2$species.year.ID == 'Lamprospilus_aff..2015')
if (length(row2rm) > 0){
  insect2 <- insect2[-row2rm,] #removes database entry with incorrect species name
}

#insect2 <- insect2[insect2$num_per_bp != 0,] # removes all comparisons with 0 diffs = could be redundant sequences
insect2 <- na.omit(insect2) #removes NAs

# add year & species (split column)
info <- unlist(strsplit(as.character(insect2$species.year.ID), '[.]'))
insect2$year <- info[seq(2, length(info), by = 2)]
insect2$year <- as.numeric(insect2$year)
insect2$species <- info[seq(1, length(info), by = 2)]
insect2$species <- as.factor(insect2$species)
rm(info)

insect2 <- insect2 %>% filter(str_count(.$species, '_') == 1, str_count(.$species, unwanted) == 0)

insect2 <- droplevels(insect2)

#summarizing data frame
insect2.agg <- insect2 %>% group_by(species, year, cell) %>% summarize(
  div = mean(num_per_bp),
  anthrome = mean(anth),
  ncomps = length(num_per_bp),
  nseqs = length(unique(c(seq1,seq2)))) %>%
  mutate_at(vars(cell), as.character) %>%
  mutate(lat = gsub("_.*","", cell), long = gsub(".*_","", cell)) %>%
  mutate_at(vars(lat:long), funs(round(as.numeric(.),3))) %>%
  mutate(cell = paste0(lat,'_',long)) %>%
  select(-anthrome) %>% as.data.frame(.)

#checking if all cells from which we have sequences are found in land use map
if(sum(!(insect2.agg$cell %in% lu.map[lu.map$mapyear == map.yrs[1],'cell'])) > 0) {print('not all cells matched')}

#adding land use and human population density infomation
insect2.agg$mapyear <- map.yr.idx$map.yr[match(insect2.agg$year, map.yr.idx$seq.yr)]
insect2.agg <- inner_join(insect2.agg,lu.map,by=c('mapyear','cell'))


# removing groups with only one sequence comparison (poor measure of diversity)
#insect2.agg <- insect2.agg[insect2.agg$ncomps > 1,]

#creating unique identifier for population
insect2.agg$pop <- as.factor(paste(insect2.agg$species,insect2.agg$cell,sep='_'))

#creating unique identifier for species/year combination
insect2.agg$sp_yr <- as.factor(paste(insect2.agg$species,insect2.agg$year,sep='_'))

#adding length of time series
pop.dur <- aggregate(div~pop, insect2.agg, FUN=length)
insect2.agg$n.years <- pop.dur$div[match(insect2.agg$pop,pop.dur$pop)]
rm(pop.dur)

# # excluding pops in the ocean
# insect2.agg <- insect2.agg[insect2.agg$anthrome != 0,]
# 
# #converting anthrome to binary
# insect2.agg$anth.cat <- 'anthro'
# insect2.agg$anth.cat[insect2.agg$anthrome >= 50] <- 'natural'
# insect2.agg$anth.cat <- as.factor(insect2.agg$anth.cat)
# 
# insect2.agg <- droplevels(insect2.agg)
# 
# # insect2.short <- insect2.agg[!(duplicated(insect2.agg$pop)),]
# nlevels(insect2.agg$species) #number of species
# nlevels(insect2.agg$pop) #number of populations
# # how many pops have data from same year, different anthromes?
# insect2.tab <- with(insect2.agg,table(sp_yr, anth.cat))
# insect2.tab <- insect2.tab[insect2.tab[,1] > 0 & insect2.tab[,2] > 0,]
# insect2.tx.spatial <- row.names(insect2.tab)
# rm(insect2.tab)
# #number of data points for spatial analysis within species @ .28
# length(insect2.tx.spatial)
# 
# #time series duration
# time.series2$insects <- numeric(length(time.series2$nb.year))
# for(i in 1:dim(time.series2)[1]){
#   time.series2$insects[i] <- length(unique(insect2.agg[insect2.agg$n.years == time.series2$nb.year[i],]$pop))
# }
# sum(time.series2$insects[5:20]) #nb of time series

## PLANTS (ITS2)

plant.its2 <- read.csv('~/Google Drive/recherche/Intraspecific genetic diversity/data/plants_ITS_pairwise_2by2_anth.csv')
colnames(plant.its2)[1] <- 'species.year.ID' #to be consistent with other files
colnames(plant.its2)[10] <- 'anth' #change new name back to old to be consistent with analysis script

# add year & species (split column)
info <- unlist(strsplit(as.character(plant.its2$species.year.ID), '[.]'))
plant.its2$year <- info[seq(2, length(info), by = 2)]
plant.its2$year <- as.numeric(plant.its2$year)
plant.its2$species <- info[seq(1, length(info), by = 2)]
plant.its2$species <- as.factor(plant.its2$species)
rm(info)

#plant.its2 <- plant.its2[plant.its2$num_per_bp != 0,] # removes all comparisons with 0 diffs = could be redundant sequences
plant.its2 <- na.omit(plant.its2) %>% filter(str_count(.$species, '_') == 1, str_count(.$species, unwanted) == 0) #removes NAs & misidentified species
plant.its2 <- droplevels(plant.its2)

#summarizing data frame
plant.its2.agg <- plant.its2 %>% group_by(species, year, cell) %>% summarize(
  div = mean(num_per_bp),
  anthrome = mean(anth),
  ncomps = length(num_per_bp),
  nseqs = length(unique(c(seq1,seq2)))) %>%
  mutate_at(vars(cell), as.character) %>%
  mutate(lat = gsub("_.*","", cell), long = gsub(".*_","", cell)) %>%
  mutate_at(vars(lat:long), funs(round(as.numeric(.),3))) %>%
  mutate(cell = paste0(lat,'_',long)) %>%
  select(-anthrome) %>% as.data.frame(.)

#checking if all cells from which we have sequences are found in land use map
if(sum(!(plant.its2.agg$cell %in% lu.map[lu.map$mapyear == map.yrs[1],'cell'])) > 0) {print('not all cells matched')}

#adding land use and human population density infomation
plant.its2.agg$mapyear <- map.yr.idx$map.yr[match(plant.its2.agg$year, map.yr.idx$seq.yr)]
plant.its2.agg <- inner_join(plant.its2.agg,lu.map,by=c('mapyear','cell'))


# removing groups with only one sequence comparison (poor measure of diversity)
#plant.its2.agg <- plant.its2.agg[plant.its2.agg$ncomps > 1,]

#creating unique identifier for population
plant.its2.agg$pop <- as.factor(paste(plant.its2.agg$species,plant.its2.agg$cell,sep='_'))

#creating unique identifier for species/year combination
plant.its2.agg$sp_yr <- as.factor(paste(plant.its2.agg$species,plant.its2.agg$year,sep='_'))

#adding length of time series
pop.dur <- aggregate(div~pop, plant.its2.agg, FUN=length)
plant.its2.agg$n.years <- pop.dur$div[match(plant.its2.agg$pop,pop.dur$pop)]
rm(pop.dur)

# # excluding pops in the water
# plant.its2.agg <- plant.its2.agg[plant.its2.agg$anthrome != 0,]
# 
# #converting anthrome to binary
# plant.its2.agg$anth.cat <- 'anthro'
# plant.its2.agg$anth.cat[plant.its2.agg$anthrome >= 50] <- 'natural'
# plant.its2.agg$anth.cat <- as.factor(plant.its2.agg$anth.cat)
# 
# plant.its2.agg <- droplevels(plant.its2.agg)
# 
# # plant.its2.short <- plant.its2.agg[!(duplicated(plant.its2.agg$pop)),]
# nlevels(plant.its2.agg$species) #number of species
# nlevels(plant.its2.agg$pop) #number of populations
# # how many pops have data from same year, different anthromes?
# plant.its2.tab <- with(plant.its2.agg,table(sp_yr, anth.cat))
# plant.its2.tab <- plant.its2.tab[plant.its2.tab[,1] > 0 & plant.its2.tab[,2] > 0,]
# plant.its2.tx.spatial <- row.names(plant.its2.tab)
# rm(plant.its2.tab)
# #number of data points for spatial analysis within species @ .28
# length(plant.its2.tx.spatial)
# 
# #time series duration
# time.series2$plant.its <- numeric(length(time.series2$nb.year))
# for(i in 1:dim(time.series2)[1]){
#   time.series2$plant.its[i] <- length(unique(plant.its2.agg[plant.its2.agg$n.years == time.series2$nb.year[i],]$pop))
# }
# sum(time.series2$plant.its[5:20]) #nb of time series

## PLANTS (MATK)

plant.matK2 <- read.csv('~/Google Drive/recherche/Intraspecific genetic diversity/data/plants_pairwise_matK_2by2_anth.csv')
colnames(plant.matK2)[1] <- 'species.year.ID' #to be consistent with other files
colnames(plant.matK2)[10] <- 'anth' #change new name back to old to be consistent with analysis script

#removing line with double period
#sapply(gregexpr("\\.", as.character(plant.matK2$species.year.ID)), tail, 1) == sapply(gregexpr("\\.", as.character(plant.matK2$species.year.ID)), head, 1)
row2rm <- which(plant.matK2$species.year.ID == 'Crataegus_indet..2011')
if (length(row2rm) > 0){
  plant.matK2 <- plant.matK2[-row2rm,] #removes database entry with incorrect species name
}

# add year & species (split column)
info <- unlist(strsplit(as.character(plant.matK2$species.year.ID), '[.]'))
plant.matK2$year <- info[seq(2, length(info), by = 2)]
plant.matK2$year <- as.numeric(plant.matK2$year)
plant.matK2$species <- info[seq(1, length(info), by = 2)]
plant.matK2$species <- as.factor(plant.matK2$species)
rm(info)

#plant.matK2 <- plant.matK2[plant.matK2$num_per_bp != 0,] # removes all comparisons with 0 diffs = could be redundant sequences
plant.matK2 <- na.omit(plant.matK2) %>% filter(str_count(.$species, '_') == 1, str_count(.$species, unwanted) == 0) #removes NAs & misidentified species
plant.matK2 <- droplevels(plant.matK2)

#summarizing data frame
plant.matK2.agg <- plant.matK2 %>% group_by(species, year, cell) %>% summarize(
  div = mean(num_per_bp),
  anthrome = mean(anth),
  ncomps = length(num_per_bp),
  nseqs = length(unique(c(seq1,seq2)))) %>%
  mutate_at(vars(cell), as.character) %>%
  mutate(lat = gsub("_.*","", cell), long = gsub(".*_","", cell)) %>%
  mutate_at(vars(lat:long), funs(round(as.numeric(.),3))) %>%
  mutate(cell = paste0(lat,'_',long)) %>%
  select(-anthrome) %>% as.data.frame(.)

#checking if all cells from which we have sequences are found in land use map
if(sum(!(plant.matK2.agg$cell %in% lu.map[lu.map$mapyear == map.yrs[1],'cell'])) > 0) {print('not all cells matched')}

#adding land use and human population density infomation
plant.matK2.agg$mapyear <- map.yr.idx$map.yr[match(plant.matK2.agg$year, map.yr.idx$seq.yr)]
plant.matK2.agg <- inner_join(plant.matK2.agg,lu.map,by=c('mapyear','cell'))


# removing groups with only one sequence comparison (poor measure of diversity)
#plant.matK2.agg <- plant.matK2.agg[plant.matK2.agg$ncomps > 1,]

#creating unique identifier for population
plant.matK2.agg$pop <- as.factor(paste(plant.matK2.agg$species,plant.matK2.agg$cell,sep='_'))

#creating unique identifier for species/year combination
plant.matK2.agg$sp_yr <- as.factor(paste(plant.matK2.agg$species,plant.matK2.agg$year,sep='_'))

#adding length of time series
pop.dur <- aggregate(div~pop, plant.matK2.agg, FUN=length)
plant.matK2.agg$n.years <- pop.dur$div[match(plant.matK2.agg$pop,pop.dur$pop)]
rm(pop.dur)

# # excluding pops in the water
# plant.matK2.agg <- plant.matK2.agg[plant.matK2.agg$anthrome != 0,]
# 
# #converting anthrome to binary
# plant.matK2.agg$anth.cat <- 'anthro'
# plant.matK2.agg$anth.cat[plant.matK2.agg$anthrome >= 50] <- 'natural'
# plant.matK2.agg$anth.cat <- as.factor(plant.matK2.agg$anth.cat)
# 
# plant.matK2.agg <- droplevels(plant.matK2.agg)
# 
# # plant.matK2.short <- plant.matK2.agg[!(duplicated(plant.matK2.agg$pop)),]
# nlevels(plant.matK2.agg$species) #number of species
# nlevels(plant.matK2.agg$pop) #number of populations
# # how many pops have data from same year, different anthromes?
# plant.matK2.tab <- with(plant.matK2.agg,table(sp_yr, anth.cat))
# plant.matK2.tab <- plant.matK2.tab[plant.matK2.tab[,1] > 0 & plant.matK2.tab[,2] > 0,]
# plant.matK2.tx.spatial <- row.names(plant.matK2.tab)
# rm(plant.matK2.tab)
# #number of data points for spatial analysis within species @ .28
# length(plant.matK2.tx.spatial)
# 
# #time series duration
# time.series2$plant.matK <- numeric(length(time.series2$nb.year))
# for(i in 1:dim(time.series2)[1]){
#   time.series2$plant.matK[i] <- length(unique(plant.matK2.agg[plant.matK2.agg$n.years == time.series2$nb.year[i],]$pop))
# }
# sum(time.series2$plant.matK[5:20]) #nb of time series

## PLANTS (rbcL)

plant.rbcL2 <- read.csv('~/Google Drive/recherche/Intraspecific genetic diversity/data/plants_pairwise_rbcL_2by2_anth.csv')
colnames(plant.rbcL2)[1] <- 'species.year.ID' #to be consistent with other files
colnames(plant.rbcL2)[10] <- 'anth' #change new name back to old to be consistent with analysis script

#which(sapply(gregexpr("\\.", as.character(plant.rbcL2$species.year.ID)), tail, 1) != sapply(gregexpr("\\.", as.character(plant.rbcL2$species.year.ID)), head, 1))
row2rm <- which(plant.rbcL2$species.year.ID == 'Crataegus_indet..2011')
if (length(row2rm) > 0){
  plant.rbcL2 <- plant.rbcL2[-row2rm,] #removes database entry with incorrect species name
}

# add year & species (split column)
info <- unlist(strsplit(as.character(plant.rbcL2$species.year.ID), '[.]'))
plant.rbcL2$year <- info[seq(2, length(info), by = 2)]
plant.rbcL2$year <- as.numeric(plant.rbcL2$year)
plant.rbcL2$species <- info[seq(1, length(info), by = 2)]
plant.rbcL2$species <- as.factor(plant.rbcL2$species)
rm(info)

#plant.rbcL2 <- plant.rbcL2[plant.rbcL2$num_per_bp != 0,] # removes all comparisons with 0 diffs = could be redundant sequences
plant.rbcL2 <- na.omit(plant.rbcL2) %>% filter(str_count(.$species, '_') == 1, str_count(.$species, unwanted) == 0) #removes NAs & misidentified species
plant.rbcL2 <- droplevels(plant.rbcL2)

#summarizing data frame
plant.rbcL2.agg <- plant.rbcL2 %>% group_by(species, year, cell) %>% summarize(
  div = mean(num_per_bp),
  anthrome = mean(anth),
  ncomps = length(num_per_bp),
  nseqs = length(unique(c(seq1,seq2)))) %>%
  mutate_at(vars(cell), as.character) %>%
  mutate(lat = gsub("_.*","", cell), long = gsub(".*_","", cell)) %>%
  mutate_at(vars(lat:long), funs(round(as.numeric(.),3))) %>%
  mutate(cell = paste0(lat,'_',long)) %>%
  select(-anthrome) %>% as.data.frame(.)

#checking if all cells from which we have sequences are found in land use map
if(sum(!(plant.rbcL2.agg$cell %in% lu.map[lu.map$mapyear == map.yrs[1],'cell'])) > 0) {print('not all cells matched')}

#adding land use and human population density infomation
plant.rbcL2.agg$mapyear <- map.yr.idx$map.yr[match(plant.rbcL2.agg$year, map.yr.idx$seq.yr)]
plant.rbcL2.agg <- inner_join(plant.rbcL2.agg,lu.map,by=c('mapyear','cell'))


# removing groups with only one sequence comparison (poor measure of diversity)
#plant.rbcL2.agg <- plant.rbcL2.agg[plant.rbcL2.agg$ncomps > 1,]

#creating unique identifier for population
plant.rbcL2.agg$pop <- as.factor(paste(plant.rbcL2.agg$species,plant.rbcL2.agg$cell,sep='_'))

#creating unique identifier for species/year combination
plant.rbcL2.agg$sp_yr <- as.factor(paste(plant.rbcL2.agg$species,plant.rbcL2.agg$year,sep='_'))

#adding length of time series
pop.dur <- aggregate(div~pop, plant.rbcL2.agg, FUN=length)
plant.rbcL2.agg$n.years <- pop.dur$div[match(plant.rbcL2.agg$pop,pop.dur$pop)]
rm(pop.dur)

# # excluding pops in the water
# plant.rbcL2.agg <- plant.rbcL2.agg[plant.rbcL2.agg$anthrome != 0,]
# 
# #converting anthrome to binary
# plant.rbcL2.agg$anth.cat <- 'anthro'
# plant.rbcL2.agg$anth.cat[plant.rbcL2.agg$anthrome >= 50] <- 'natural'
# plant.rbcL2.agg$anth.cat <- as.factor(plant.rbcL2.agg$anth.cat)
# 
# plant.rbcL2.agg <- droplevels(plant.rbcL2.agg)
# 
# # plant.rbcL2.short <- plant.rbcL2.agg[!(duplicated(plant.rbcL2.agg$pop)),]
# nlevels(plant.rbcL2.agg$species) #number of species
# nlevels(plant.rbcL2.agg$pop) #number of populations
# # how many pops have data from same year, different anthromes?
# plant.rbcL2.tab <- with(plant.rbcL2.agg,table(sp_yr, anth.cat))
# plant.rbcL2.tab <- plant.rbcL2.tab[plant.rbcL2.tab[,1] > 0 & plant.rbcL2.tab[,2] > 0,]
# plant.rbcL2.tx.spatial <- row.names(plant.rbcL2.tab)
# rm(plant.rbcL2.tab)
# #number of data points for spatial analysis within species @ .28
# length(plant.rbcL2.tx.spatial)
# 
# #time series duration
# time.series2$plant.rbcL <- numeric(length(time.series2$nb.year))
# for(i in 1:dim(time.series2)[1]){
#   time.series2$plant.rbcL[i] <- length(unique(plant.rbcL2.agg[plant.rbcL2.agg$n.years == time.series2$nb.year[i],]$pop))
# }
# sum(time.series2$plant.rbcL[5:20]) #nb of time series

#gather data summary statistics
summ2 <- data.frame('taxon'=c('mammals','birds','fish','insects','plants_its','plants_matK','plants_rbcL'),
                     'nb.sequences' = numeric(7),
                     'nb.species' = numeric(7),
                     'nb.pops' = numeric(7))
summ2[1,2:4] <- c(sum(mam2.agg$nseqs),length(levels(mam2.agg$species)),length(levels(mam2.agg$pop)))
summ2[2,2:4] <- c(sum(aves2.agg$nseqs),length(levels(aves2.agg$species)),length(levels(aves2.agg$pop)))
summ2[3,2:4] <- c(sum(acti2.agg$nseqs),length(levels(acti2.agg$species)),length(levels(acti2.agg$pop)))
summ2[4,2:4] <- c(sum(insect2.agg$nseqs),length(levels(insect2.agg$species)),length(levels(insect2.agg$pop)))
summ2[5,2:4] <- c(sum(plant.its2.agg$nseqs),length(levels(plant.its2.agg$species)),length(levels(plant.its2.agg$pop)))
summ2[6,2:4] <- c(sum(plant.matK2.agg$nseqs),length(levels(plant.matK2.agg$species)),length(levels(plant.matK2.agg$pop)))
summ2[7,2:4] <- c(sum(plant.rbcL2.agg$nseqs),length(levels(plant.rbcL2.agg$species)),length(levels(plant.rbcL2.agg$pop)))

#cleanup
rm(i,row2rm,acti2,aves2,insect2,mam2,plant.its2,plant.matK2,plant.rbcL2,lu.map)

mam2.agg <- as.data.frame(mam2.agg)
aves2.agg <- as.data.frame(aves2.agg)
acti2.agg <- as.data.frame(acti2.agg)
insect2.agg <- as.data.frame(insect2.agg)
plant.its2.agg <- as.data.frame(plant.its2.agg)
plant.matK2.agg <- as.data.frame(plant.matK2.agg)
plant.rbcL2.agg <- as.data.frame(plant.rbcL2.agg)

#### scale of aggregation = 4' grid cells ####

#load and format maps
for(i in 1:length(map.yrs)){
  fname <- paste0('hyde32_',map.yrs[i],'_4')
  path <- paste0('~/Google Drive/recherche/Intraspecific genetic diversity/data/',fname,'.RData')
  load(path)
  if(i == 1){
    get(fname) %>% filter(!is.na(pop_tot)) %>%
      mutate('cell' = paste0(.data$lat,'_',.data$long)) %>% select(-lat, -long) %>%
      mutate('mapyear' = map.yrs[i]) %>% select(mapyear,cell,everything()) -> lu.map
  } else {
    get(fname) %>% filter(!is.na(pop_tot)) %>%
      mutate('cell' = paste0(.data$lat,'_',.data$long)) %>% select(-lat, -long) %>%
      mutate('mapyear' = map.yrs[i]) %>% select(mapyear,cell,everything()) -> tmp.map
    lu.map <- bind_rows(lu.map,tmp.map)
  }
  rm(list = c(fname))
}
rm(tmp.map)

# load('~/Google Drive/recherche/Intraspecific genetic diversity/data/hyde32_2017_4.RData')
# load('~/Google Drive/recherche/Intraspecific genetic diversity/data/hyde32_1980_4.RData')
# old <- hyde32_1980_4 %>% mutate('cell' = paste0(.data$lat,'_',.data$long)) %>%
#   #('pop_density' = pop_tot, 'p_pasture' = grazing, 'p_urban' = urban, 'p_crop' = cropland) %>%
#   mutate('cell.r' = as.factor(cell)) %>% select(-c(lat,long,cell))
# lu4 <- hyde32_2017_4 %>% mutate('cell' = paste0(.data$lat,'_',.data$long)) %>%
#   #rename('pop_density' = pop_tot, 'p_pasture' = grazing, 'p_urban' = urban, 'p_crop' = cropland) %>%
#   mutate('cell.r' = as.factor(cell))
# lu4 <- left_join(x = lu4, y = old, by = 'cell.r' , suffix = c('.2017','.1980')) %>%
#   select(lat,long,cell,cell.r,everything())
# rm(hyde32_1980_4, hyde32_2017_4,old)

# lu4 <- read_csv('~/Google Drive/recherche/Intraspecific genetic diversity/data/landuse_4.csv') %>%
#   select(-starts_with('X'))
# lu4$cell <- gsub(' ','',lu4$cell)
# lu4 <- lu4 %>% mutate(lat = gsub("_.*","", cell), long = gsub(".*_","", cell)) %>%
#   mutate_at(vars(lat:long), funs(round(as.numeric(.),3))) %>%
#   mutate(cell.r = paste(lat,'_',long,sep=''))
# lu4$cell.r <-as.factor(lu4$cell.r)

## MAMMALS

mam4 <- read.csv('~/Google Drive/recherche/Intraspecific genetic diversity/data/mamm_pairwise_4by4_anth.csv')
colnames(mam4)[1] <- 'species.year.ID' #to be consistent with other files
colnames(mam4)[10] <- 'anth' #change new name back to old to be consistent with analysis script

# add year & species (split column)
info <- unlist(strsplit(as.character(mam4$species.year.ID), '[.]'))
mam4$year <- info[seq(2, length(info), by = 2)]
mam4$year <- as.numeric(mam4$year)
mam4$species <- info[seq(1, length(info), by = 2)]
mam4$species <- as.factor(mam4$species)
rm(info)

#mam4 <- mam4[mam4$num_per_bp != 0,] # removes all comparisons with 0 diffs = could be redundant sequences
mam4 <- na.omit(mam4) %>% filter(str_count(.$species, '_') == 1, str_count(.$species, unwanted) == 0) #removes NAs & misidentified species
mam4 <- droplevels(mam4)

#summarizing data frame
mam4.agg <- mam4 %>% group_by(species, year, cell) %>% summarize(
  div = mean(num_per_bp),
  anthrome = mean(anth),
  ncomps = length(num_per_bp),
  nseqs = length(unique(c(seq1,seq2)))) %>%
  mutate_at(vars(cell), as.character) %>%
  mutate(lat = gsub("_.*","", cell), long = gsub(".*_","", cell)) %>%
  mutate_at(vars(lat:long), funs(round(as.numeric(.),3))) %>%
  mutate(cell = paste0(lat,'_',long)) %>%
  select(-anthrome) %>% as.data.frame(.)

#checking if all cells from which we have sequences are found in land use map
if(sum(!(mam4.agg$cell %in% lu.map[lu.map$mapyear == map.yrs[1],'cell'])) > 0) {print('not all cells matched')}

#adding land use and human population density infomation
mam4.agg$mapyear <- map.yr.idx$map.yr[match(mam4.agg$year, map.yr.idx$seq.yr)]
mam4.agg <- inner_join(mam4.agg,lu.map,by=c('mapyear','cell'))

# removing groups with only one sequence comparison (poor measure of diversity)
#mam4.agg <- mam4.agg[mam4.agg$ncomps > 1,]

#creating unique identifier for population
mam4.agg$pop <- as.factor(paste(mam4.agg$species,mam4.agg$cell,sep='_'))

#creating unique identifier for species/year combination
mam4.agg$sp_yr <- as.factor(paste(mam4.agg$species,mam4.agg$year,sep='_'))

#adding length of time series
pop.dur <- aggregate(div~pop, mam4.agg, FUN=length)
mam4.agg$n.years <- pop.dur$div[match(mam4.agg$pop,pop.dur$pop)]
rm(pop.dur)

# # excluding pops in the water
# mam4.agg <- mam4.agg[mam4.agg$anthrome != 0,]
# 
# #converting anthrome to binary. natural are all 'semi-natural' and 'wild' in Ellis 2010
# mam4.agg$anth.cat <- 'anthro'
# mam4.agg$anth.cat[mam4.agg$anthrome >= 50] <- 'natural'
# mam4.agg$anth.cat <- as.factor(mam4.agg$anth.cat)
# 
# mam4.agg <- droplevels(mam4.agg)
# 
# # mam4.short <- mam4.agg[!(duplicated(mam4.agg$pop)),]
# nlevels(mam4.agg$species) #number of species
# nlevels(mam4.agg$pop) #number of populations
# # how many pops have data from same year, different anthromes?
# mam4.tab <- with(mam4.agg,table(sp_yr, anth.cat))
# mam4.tab <- mam4.tab[mam4.tab[,1] > 0 & mam4.tab[,2] > 0,]
# mam4.tx.spatial <- row.names(mam4.tab)
# rm(mam4.tab)
# #number of data points for spatial analysis within species @ .28
# length(mam4.tx.spatial)
# 
# #time series duration
# time.series4 <- data.frame('nb.year' = 1:20)
# time.series4$mams <- numeric(length(time.series4$nb.year))
# for(i in 1:dim(time.series4)[1]){
#   time.series4$mams[i] <- length(unique(mam4.agg[mam4.agg$n.years == time.series4$nb.year[i],]$pop))
# }
# sum(time.series4$mams[5:20]) #number of time series > 5 yrs

## BIRDS

aves4 <- read.csv('~/Google Drive/recherche/Intraspecific genetic diversity/data/aves_pairwise_4by4_anth.csv')
colnames(aves4)[1] <- 'species.year.ID' #to be consistent with other files
colnames(aves4)[10] <- 'anth' #change new name back to old to be consistent with analysis script

# add year & species (split column)
info <- unlist(strsplit(as.character(aves4$species.year.ID), '[.]'))
aves4$year <- info[seq(2, length(info), by = 2)]
aves4$year <- as.numeric(aves4$year)
aves4$species <- info[seq(1, length(info), by = 2)]
aves4$species <- as.factor(aves4$species)
rm(info)

#aves4 <- aves4[aves4$num_per_bp != 0,] # removes all comparisons with 0 diffs = could be redundant sequences
aves4 <- na.omit(aves4) %>% filter(str_count(.$species, '_') == 1, str_count(.$species, unwanted) == 0) #removes NAs & misidentified species
aves4 <- droplevels(aves4)

#summarizing data frame
aves4.agg <- aves4 %>% group_by(species, year, cell) %>% summarize(
  div = mean(num_per_bp),
  anthrome = mean(anth),
  ncomps = length(num_per_bp),
  nseqs = length(unique(c(seq1,seq2)))) %>%
  mutate_at(vars(cell), as.character) %>%
  mutate(lat = gsub("_.*","", cell), long = gsub(".*_","", cell)) %>%
  mutate_at(vars(lat:long), funs(round(as.numeric(.),3))) %>%
  mutate(cell = paste0(lat,'_',long)) %>%
  select(-anthrome) %>% as.data.frame(.)

#checking if all cells from which we have sequences are found in land use map
if(sum(!(aves4.agg$cell %in% lu.map[lu.map$mapyear == map.yrs[1],'cell'])) > 0) {print('not all cells matched')}

#adding land use and human population density infomation
aves4.agg$mapyear <- map.yr.idx$map.yr[match(aves4.agg$year, map.yr.idx$seq.yr)]
aves4.agg <- inner_join(aves4.agg,lu.map,by=c('mapyear','cell'))

# removing groups with only one sequence comparison (poor measure of diversity)
#aves4.agg <- aves4.agg[aves4.agg$ncomps > 1,]

#creating unique identifier for population
aves4.agg$pop <- as.factor(paste(aves4.agg$species,aves4.agg$cell,sep='_'))

#creating unique identifier for species/year combination
aves4.agg$sp_yr <- as.factor(paste(aves4.agg$species,aves4.agg$year,sep='_'))

#adding length of time series
pop.dur <- aggregate(div~pop, aves4.agg, FUN=length)
aves4.agg$n.years <- pop.dur$div[match(aves4.agg$pop,pop.dur$pop)]
rm(pop.dur)

# # excluding pops in the water
# aves4.agg <- aves4.agg[aves4.agg$anthrome != 0,]
# 
# #converting anthrome to binary
# aves4.agg$anth.cat <- 'anthro'
# aves4.agg$anth.cat[aves4.agg$anthrome >= 50] <- 'natural'
# aves4.agg$anth.cat <- as.factor(aves4.agg$anth.cat)
# 
# aves4.agg <- droplevels(aves4.agg)
# 
# # aves4.short <- aves4.agg[!(duplicated(aves4.agg$pop)),]
# nlevels(aves4.agg$species) #number of species
# nlevels(aves4.agg$pop) #number of populations
# # how many pops have data from same year, different anthromes?
# aves4.tab <- with(aves4.agg,table(sp_yr, anth.cat))
# aves4.tab <- aves4.tab[aves4.tab[,1] > 0 & aves4.tab[,2] > 0,]
# aves4.tx.spatial <- row.names(aves4.tab)
# rm(aves4.tab)
# #number of data points for spatial analysis within species @ .28
# length(aves4.tx.spatial)
# 
# #time series duration
# time.series4$aves <- numeric(length(time.series4$nb.year))
# for(i in 1:dim(time.series4)[1]){
#   time.series4$aves[i] <- length(unique(aves4.agg[aves4.agg$n.years == time.series4$nb.year[i],]$pop))
# }
# sum(time.series4$aves[5:20]) #nb of time series

## BONY FISHES

acti4 <- read.csv('~/Google Drive/recherche/Intraspecific genetic diversity/data/acti_pairwise_4by4_anth.csv')
colnames(acti4)[1] <- 'species.year.ID' #to be consistent with other files
colnames(acti4)[10] <- 'anth' #change new name back to old to be consistent with analysis script

# add year & species (split column)
info <- unlist(strsplit(as.character(acti4$species.year.ID), '[.]'))
acti4$year <- info[seq(2, length(info), by = 2)]
acti4$year <- as.numeric(acti4$year)
acti4$species <- info[seq(1, length(info), by = 2)]
acti4$species <- as.factor(acti4$species)
rm(info)

#acti4 <- acti4[acti4$num_per_bp != 0,] # removes all comparisons with 0 diffs = could be redundant sequences
acti4 <- na.omit(acti4) %>% filter(str_count(.$species, '_') == 1, str_count(.$species, unwanted) == 0) #removes NAs & misidentified species
acti4 <- droplevels(acti4)

#summarizing data frame
acti4.agg <- acti4 %>% group_by(species, year, cell) %>% summarize(
  div = mean(num_per_bp),
  anthrome = mean(anth),
  ncomps = length(num_per_bp),
  nseqs = length(unique(c(seq1,seq2)))) %>%
  mutate_at(vars(cell), as.character) %>%
  mutate(lat = gsub("_.*","", cell), long = gsub(".*_","", cell)) %>%
  mutate_at(vars(lat:long), funs(round(as.numeric(.),3))) %>%
  mutate(cell = paste0(lat,'_',long)) %>%
  select(-anthrome) %>% as.data.frame(.)

#checking if all cells from which we have sequences are found in land use map
if(sum(!(acti4.agg$cell %in% lu.map[lu.map$mapyear == map.yrs[1],'cell'])) > 0) {print('not all cells matched')}

#adding land use and human population density infomation
acti4.agg$mapyear <- map.yr.idx$map.yr[match(acti4.agg$year, map.yr.idx$seq.yr)]
acti4.agg <- inner_join(acti4.agg,lu.map,by=c('mapyear','cell'))

# removing groups with only one sequence comparison (poor measure of diversity)
#acti4.agg <- acti4.agg[acti4.agg$ncomps > 1,]

#creating unique identifier for population
acti4.agg$pop <- as.factor(paste(acti4.agg$species,acti4.agg$cell,sep='_'))

#creating unique identifier for species/year combination
acti4.agg$sp_yr <- as.factor(paste(acti4.agg$species,acti4.agg$year,sep='_'))

#adding length of time series
pop.dur <- aggregate(div~pop, acti4.agg, FUN=length)
acti4.agg$n.years <- pop.dur$div[match(acti4.agg$pop,pop.dur$pop)]
rm(pop.dur)

# # excluding pops in the ocean
# acti4.agg <- acti4.agg[acti4.agg$anthrome != 0,]
# 
# #converting anthrome to binary
# acti4.agg$anth.cat <- 'anthro'
# acti4.agg$anth.cat[acti4.agg$anthrome >= 50] <- 'natural'
# acti4.agg$anth.cat <- as.factor(acti4.agg$anth.cat)
# 
# acti4.agg <- droplevels(acti4.agg)
# 
# # acti4.short <- acti4.agg[!(duplicated(acti4.agg$pop)),]
# nlevels(acti4.agg$species) #number of species
# nlevels(acti4.agg$pop) #number of populations
# # how many pops have data from same year, different anthromes?
# acti4.tab <- with(acti4.agg,table(sp_yr, anth.cat))
# acti4.tab <- acti4.tab[acti4.tab[,1] > 0 & acti4.tab[,2] > 0,]
# acti4.tx.spatial <- row.names(acti4.tab)
# rm(acti4.tab)
# #number of data points for spatial analysis within species @ .28
# length(acti4.tx.spatial)
# 
# #time series duration
# time.series4$acti <- numeric(length(time.series4$nb.year))
# for(i in 1:dim(time.series4)[1]){
#   time.series4$acti[i] <- length(unique(acti4.agg[acti4.agg$n.years == time.series4$nb.year[i],]$pop))
# }
# sum(time.series4$acti[5:20]) #nb of time series
# 
## INSECTS

insect4 <- read.csv('~/Google Drive/recherche/Intraspecific genetic diversity/data/insect_pairwise_4by4_anth.csv')
colnames(insect4)[1] <- 'species.year.ID' #to be consistent with other files
colnames(insect4)[10] <- 'anth' #change new name back to old to be consistent with analysis script

insect4 <- insect4[insect4$overlap > 0.5,] # filter comparisons with less than 50% overlap as in Miraldo paper (for other taxa, this step was done in Excel)
row2rm <- which(insect4$species.year.ID == 'Lamprospilus_aff..2015')
if (length(row2rm) > 0){
  insect4 <- insect4[-row2rm,] #removes database entry with incorrect species name
}

#insect4 <- insect4[insect4$num_per_bp != 0,] # removes all comparisons with 0 diffs = could be redundant sequences
insect4 <- na.omit(insect4) #removes NAs

# add year & species (split column)
info <- unlist(strsplit(as.character(insect4$species.year.ID), '[.]'))
insect4$year <- info[seq(2, length(info), by = 2)]
insect4$year <- as.numeric(insect4$year)
insect4$species <- info[seq(1, length(info), by = 2)]
insect4$species <- as.factor(insect4$species)
rm(info)

insect4 <- insect4 %>% filter(str_count(.$species, '_') == 1, str_count(.$species, unwanted) == 0)

insect4 <- droplevels(insect4)

#summarizing data frame
insect4.agg <- insect4 %>% group_by(species, year, cell) %>% summarize(
  div = mean(num_per_bp),
  anthrome = mean(anth),
  ncomps = length(num_per_bp),
  nseqs = length(unique(c(seq1,seq2)))) %>%
  mutate_at(vars(cell), as.character) %>%
  mutate(lat = gsub("_.*","", cell), long = gsub(".*_","", cell)) %>%
  mutate_at(vars(lat:long), funs(round(as.numeric(.),3))) %>%
  mutate(cell = paste0(lat,'_',long)) %>%
  select(-anthrome) %>% as.data.frame(.)

#checking if all cells from which we have sequences are found in land use map
if(sum(!(insect4.agg$cell %in% lu.map[lu.map$mapyear == map.yrs[1],'cell'])) > 0) {print('not all cells matched')}

#adding land use and human population density infomation
insect4.agg$mapyear <- map.yr.idx$map.yr[match(insect4.agg$year, map.yr.idx$seq.yr)]
insect4.agg <- inner_join(insect4.agg,lu.map,by=c('mapyear','cell'))

# removing groups with only one sequence comparison (poor measure of diversity)
#insect4.agg <- insect4.agg[insect4.agg$ncomps > 1,]

#creating unique identifier for population
insect4.agg$pop <- as.factor(paste(insect4.agg$species,insect4.agg$cell,sep='_'))

#creating unique identifier for species/year combination
insect4.agg$sp_yr <- as.factor(paste(insect4.agg$species,insect4.agg$year,sep='_'))

#adding length of time series
pop.dur <- aggregate(div~pop, insect4.agg, FUN=length)
insect4.agg$n.years <- pop.dur$div[match(insect4.agg$pop,pop.dur$pop)]
rm(pop.dur)

# # excluding pops in the ocean
# insect4.agg <- insect4.agg[insect4.agg$anthrome != 0,]
# 
# #converting anthrome to binary
# insect4.agg$anth.cat <- 'anthro'
# insect4.agg$anth.cat[insect4.agg$anthrome >= 50] <- 'natural'
# insect4.agg$anth.cat <- as.factor(insect4.agg$anth.cat)
# 
# insect4.agg <- droplevels(insect4.agg)
# 
# # insect4.short <- insect4.agg[!(duplicated(insect4.agg$pop)),]
# nlevels(insect4.agg$species) #number of species
# nlevels(insect4.agg$pop) #number of populations
# # how many pops have data from same year, different anthromes?
# insect4.tab <- with(insect4.agg,table(sp_yr, anth.cat))
# insect4.tab <- insect4.tab[insect4.tab[,1] > 0 & insect4.tab[,2] > 0,]
# insect4.tx.spatial <- row.names(insect4.tab)
# rm(insect4.tab)
# #number of data points for spatial analysis within species @ .28
# length(insect4.tx.spatial)
# 
# #time series duration
# time.series4$insects <- numeric(length(time.series4$nb.year))
# for(i in 1:dim(time.series4)[1]){
#   time.series4$insects[i] <- length(unique(insect4.agg[insect4.agg$n.years == time.series4$nb.year[i],]$pop))
# }
# sum(time.series4$insects[5:20]) #nb of time series

## PLANTS (ITS2)

plant.its4 <- read.csv('~/Google Drive/recherche/Intraspecific genetic diversity/data/plants_ITS_pairwise_4by4_anth.csv')
colnames(plant.its4)[1] <- 'species.year.ID' #to be consistent with other files
colnames(plant.its4)[10] <- 'anth' #change new name back to old to be consistent with analysis script

# add year & species (split column)
info <- unlist(strsplit(as.character(plant.its4$species.year.ID), '[.]'))
plant.its4$year <- info[seq(2, length(info), by = 2)]
plant.its4$year <- as.numeric(plant.its4$year)
plant.its4$species <- info[seq(1, length(info), by = 2)]
plant.its4$species <- as.factor(plant.its4$species)
rm(info)

#plant.its4 <- plant.its4[plant.its4$num_per_bp != 0,] # removes all comparisons with 0 diffs = could be redundant sequences
plant.its4 <- na.omit(plant.its4) %>% filter(str_count(.$species, '_') == 1, str_count(.$species, unwanted) == 0) #removes NAs & misidentified species
plant.its4 <- droplevels(plant.its4)

#summarizing data frame
plant.its4.agg <- plant.its4 %>% group_by(species, year, cell) %>% summarize(
  div = mean(num_per_bp),
  anthrome = mean(anth),
  ncomps = length(num_per_bp),
  nseqs = length(unique(c(seq1,seq2)))) %>%
  mutate_at(vars(cell), as.character) %>%
  mutate(lat = gsub("_.*","", cell), long = gsub(".*_","", cell)) %>%
  mutate_at(vars(lat:long), funs(round(as.numeric(.),3))) %>%
  mutate(cell = paste0(lat,'_',long)) %>%
  select(-anthrome) %>% as.data.frame(.)

#checking if all cells from which we have sequences are found in land use map
if(sum(!(plant.its4.agg$cell %in% lu.map[lu.map$mapyear == map.yrs[1],'cell'])) > 0) {print('not all cells matched')}

#adding land use and human population density infomation
plant.its4.agg$mapyear <- map.yr.idx$map.yr[match(plant.its4.agg$year, map.yr.idx$seq.yr)]
plant.its4.agg <- inner_join(plant.its4.agg,lu.map,by=c('mapyear','cell'))

# removing groups with only one sequence comparison (poor measure of diversity)
#plant.its4.agg <- plant.its4.agg[plant.its4.agg$ncomps > 1,]

#creating unique identifier for population
plant.its4.agg$pop <- as.factor(paste(plant.its4.agg$species,plant.its4.agg$cell,sep='_'))

#creating unique identifier for species/year combination
plant.its4.agg$sp_yr <- as.factor(paste(plant.its4.agg$species,plant.its4.agg$year,sep='_'))

#adding length of time series
pop.dur <- aggregate(div~pop, plant.its4.agg, FUN=length)
plant.its4.agg$n.years <- pop.dur$div[match(plant.its4.agg$pop,pop.dur$pop)]
rm(pop.dur)

# # excluding pops in the water
# plant.its4.agg <- plant.its4.agg[plant.its4.agg$anthrome != 0,]
# 
# #converting anthrome to binary
# plant.its4.agg$anth.cat <- 'anthro'
# plant.its4.agg$anth.cat[plant.its4.agg$anthrome >= 50] <- 'natural'
# plant.its4.agg$anth.cat <- as.factor(plant.its4.agg$anth.cat)
# 
# plant.its4.agg <- droplevels(plant.its4.agg)
# 
# # plant.its4.short <- plant.its4.agg[!(duplicated(plant.its4.agg$pop)),]
# nlevels(plant.its4.agg$species) #number of species
# nlevels(plant.its4.agg$pop) #number of populations
# # how many pops have data from same year, different anthromes?
# plant.its4.tab <- with(plant.its4.agg,table(sp_yr, anth.cat))
# plant.its4.tab <- plant.its4.tab[plant.its4.tab[,1] > 0 & plant.its4.tab[,2] > 0,]
# plant.its4.tx.spatial <- row.names(plant.its4.tab)
# rm(plant.its4.tab)
# #number of data points for spatial analysis within species @ .28
# length(plant.its4.tx.spatial)
# 
# #time series duration
# time.series4$plant.its <- numeric(length(time.series4$nb.year))
# for(i in 1:dim(time.series4)[1]){
#   time.series4$plant.its[i] <- length(unique(plant.its4.agg[plant.its4.agg$n.years == time.series4$nb.year[i],]$pop))
# }
# sum(time.series4$plant.its[5:20]) #nb of time series

## PLANTS (MATK)

plant.matK4 <- read.csv('~/Google Drive/recherche/Intraspecific genetic diversity/data/plants_pairwise_matK_4by4_anth.csv')
colnames(plant.matK4)[1] <- 'species.year.ID' #to be consistent with other files
colnames(plant.matK4)[10] <- 'anth' #change new name back to old to be consistent with analysis script

#removing line with double period
#sapply(gregexpr("\\.", as.character(plant.matK4$species.year.ID)), tail, 1) == sapply(gregexpr("\\.", as.character(plant.matK4$species.year.ID)), head, 1)
row2rm <- which(plant.matK4$species.year.ID == 'Crataegus_indet..2011')
if (length(row2rm) > 0){
  plant.matK4 <- plant.matK4[-row2rm,] #removes database entry with incorrect species name
}

# add year & species (split column)
info <- unlist(strsplit(as.character(plant.matK4$species.year.ID), '[.]'))
plant.matK4$year <- info[seq(2, length(info), by = 2)]
plant.matK4$year <- as.numeric(plant.matK4$year)
plant.matK4$species <- info[seq(1, length(info), by = 2)]
plant.matK4$species <- as.factor(plant.matK4$species)
rm(info)

#plant.matK4 <- plant.matK4[plant.matK4$num_per_bp != 0,] # removes all comparisons with 0 diffs = could be redundant sequences
plant.matK4 <- na.omit(plant.matK4) %>% filter(str_count(.$species, '_') == 1, str_count(.$species, unwanted) == 0) #removes NAs & misidentified species
plant.matK4 <- droplevels(plant.matK4)

#summarizing data frame
plant.matK4.agg <- plant.matK4 %>% group_by(species, year, cell) %>% summarize(
  div = mean(num_per_bp),
  anthrome = mean(anth),
  ncomps = length(num_per_bp),
  nseqs = length(unique(c(seq1,seq2)))) %>%
  mutate_at(vars(cell), as.character) %>%
  mutate(lat = gsub("_.*","", cell), long = gsub(".*_","", cell)) %>%
  mutate_at(vars(lat:long), funs(round(as.numeric(.),3))) %>%
  mutate(cell = paste0(lat,'_',long)) %>%
  select(-anthrome) %>% as.data.frame(.)

#checking if all cells from which we have sequences are found in land use map
if(sum(!(plant.matK4.agg$cell %in% lu.map[lu.map$mapyear == map.yrs[1],'cell'])) > 0) {print('not all cells matched')}

#adding land use and human population density infomation
plant.matK4.agg$mapyear <- map.yr.idx$map.yr[match(plant.matK4.agg$year, map.yr.idx$seq.yr)]
plant.matK4.agg <- inner_join(plant.matK4.agg,lu.map,by=c('mapyear','cell'))

# removing groups with only one sequence comparison (poor measure of diversity)
#plant.matK4.agg <- plant.matK4.agg[plant.matK4.agg$ncomps > 1,]

#creating unique identifier for population
plant.matK4.agg$pop <- as.factor(paste(plant.matK4.agg$species,plant.matK4.agg$cell,sep='_'))

#creating unique identifier for species/year combination
plant.matK4.agg$sp_yr <- as.factor(paste(plant.matK4.agg$species,plant.matK4.agg$year,sep='_'))

#adding length of time series
pop.dur <- aggregate(div~pop, plant.matK4.agg, FUN=length)
plant.matK4.agg$n.years <- pop.dur$div[match(plant.matK4.agg$pop,pop.dur$pop)]
rm(pop.dur)

# # excluding pops in the water
# plant.matK4.agg <- plant.matK4.agg[plant.matK4.agg$anthrome != 0,]
# 
# #converting anthrome to binary
# plant.matK4.agg$anth.cat <- 'anthro'
# plant.matK4.agg$anth.cat[plant.matK4.agg$anthrome >= 50] <- 'natural'
# plant.matK4.agg$anth.cat <- as.factor(plant.matK4.agg$anth.cat)
# 
# plant.matK4.agg <- droplevels(plant.matK4.agg)
# 
# # plant.matK4.short <- plant.matK4.agg[!(duplicated(plant.matK4.agg$pop)),]
# nlevels(plant.matK4.agg$species) #number of species
# nlevels(plant.matK4.agg$pop) #number of populations
# # how many pops have data from same year, different anthromes?
# plant.matK4.tab <- with(plant.matK4.agg,table(sp_yr, anth.cat))
# plant.matK4.tab <- plant.matK4.tab[plant.matK4.tab[,1] > 0 & plant.matK4.tab[,2] > 0,]
# plant.matK4.tx.spatial <- row.names(plant.matK4.tab)
# rm(plant.matK4.tab)
# #number of data points for spatial analysis within species @ .28
# length(plant.matK4.tx.spatial)
# 
# #time series duration
# time.series4$plant.matK <- numeric(length(time.series4$nb.year))
# for(i in 1:dim(time.series4)[1]){
#   time.series4$plant.matK[i] <- length(unique(plant.matK4.agg[plant.matK4.agg$n.years == time.series4$nb.year[i],]$pop))
# }
# sum(time.series4$plant.matK[5:20]) #nb of time series

## PLANTS (rbcL)

plant.rbcL4 <- read.csv('~/Google Drive/recherche/Intraspecific genetic diversity/data/plants_pairwise_rbcL_4by4_anth.csv')
colnames(plant.rbcL4)[1] <- 'species.year.ID' #to be consistent with other files
colnames(plant.rbcL4)[10] <- 'anth' #change new name back to old to be consistent with analysis script

#which(sapply(gregexpr("\\.", as.character(plant.rbcL4$species.year.ID)), tail, 1) != sapply(gregexpr("\\.", as.character(plant.rbcL4$species.year.ID)), head, 1))
row2rm <- which(plant.rbcL4$species.year.ID == 'Crataegus_indet..2011')
if (length(row2rm) > 0){
  plant.rbcL4 <- plant.rbcL4[-row2rm,] #removes database entry with incorrect species name
}

# add year & species (split column)
info <- unlist(strsplit(as.character(plant.rbcL4$species.year.ID), '[.]'))
plant.rbcL4$year <- info[seq(2, length(info), by = 2)]
plant.rbcL4$year <- as.numeric(plant.rbcL4$year)
plant.rbcL4$species <- info[seq(1, length(info), by = 2)]
plant.rbcL4$species <- as.factor(plant.rbcL4$species)
rm(info)

#plant.rbcL4 <- plant.rbcL4[plant.rbcL4$num_per_bp != 0,] # removes all comparisons with 0 diffs = could be redundant sequences
plant.rbcL4 <- na.omit(plant.rbcL4) %>% filter(str_count(.$species, '_') == 1, str_count(.$species, unwanted) == 0) #removes NAs & misidentified species
plant.rbcL4 <- droplevels(plant.rbcL4)

#summarizing data frame
plant.rbcL4.agg <- plant.rbcL4 %>% group_by(species, year, cell) %>% summarize(
  div = mean(num_per_bp),
  anthrome = mean(anth),
  ncomps = length(num_per_bp),
  nseqs = length(unique(c(seq1,seq2)))) %>%
  mutate_at(vars(cell), as.character) %>%
  mutate(lat = gsub("_.*","", cell), long = gsub(".*_","", cell)) %>%
  mutate_at(vars(lat:long), funs(round(as.numeric(.),3))) %>%
  mutate(cell = paste0(lat,'_',long)) %>%
  select(-anthrome) %>% as.data.frame(.)

#checking if all cells from which we have sequences are found in land use map
if(sum(!(plant.rbcL4.agg$cell %in% lu.map[lu.map$mapyear == map.yrs[1],'cell'])) > 0) {print('not all cells matched')}

#adding land use and human population density infomation
plant.rbcL4.agg$mapyear <- map.yr.idx$map.yr[match(plant.rbcL4.agg$year, map.yr.idx$seq.yr)]
plant.rbcL4.agg <- inner_join(plant.rbcL4.agg,lu.map,by=c('mapyear','cell'))

# removing groups with only one sequence comparison (poor measure of diversity)
#plant.rbcL4.agg <- plant.rbcL4.agg[plant.rbcL4.agg$ncomps > 1,]

#creating unique identifier for population
plant.rbcL4.agg$pop <- as.factor(paste(plant.rbcL4.agg$species,plant.rbcL4.agg$cell,sep='_'))

#creating unique identifier for species/year combination
plant.rbcL4.agg$sp_yr <- as.factor(paste(plant.rbcL4.agg$species,plant.rbcL4.agg$year,sep='_'))

#adding length of time series
pop.dur <- aggregate(div~pop, plant.rbcL4.agg, FUN=length)
plant.rbcL4.agg$n.years <- pop.dur$div[match(plant.rbcL4.agg$pop,pop.dur$pop)]
rm(pop.dur)

# # excluding pops in the water
# plant.rbcL4.agg <- plant.rbcL4.agg[plant.rbcL4.agg$anthrome != 0,]
# 
# #converting anthrome to binary
# plant.rbcL4.agg$anth.cat <- 'anthro'
# plant.rbcL4.agg$anth.cat[plant.rbcL4.agg$anthrome >= 50] <- 'natural'
# plant.rbcL4.agg$anth.cat <- as.factor(plant.rbcL4.agg$anth.cat)
# 
# plant.rbcL4.agg <- droplevels(plant.rbcL4.agg)
# 
# # plant.rbcL4.short <- plant.rbcL4.agg[!(duplicated(plant.rbcL4.agg$pop)),]
# nlevels(plant.rbcL4.agg$species) #number of species
# nlevels(plant.rbcL4.agg$pop) #number of populations
# # how many pops have data from same year, different anthromes?
# plant.rbcL4.tab <- with(plant.rbcL4.agg,table(sp_yr, anth.cat))
# plant.rbcL4.tab <- plant.rbcL4.tab[plant.rbcL4.tab[,1] > 0 & plant.rbcL4.tab[,2] > 0,]
# plant.rbcL4.tx.spatial <- row.names(plant.rbcL4.tab)
# rm(plant.rbcL4.tab)
# #number of data points for spatial analysis within species @ .28
# length(plant.rbcL4.tx.spatial)
# 
# #time series duration
# time.series4$plant.rbcL <- numeric(length(time.series4$nb.year))
# for(i in 1:dim(time.series4)[1]){
#   time.series4$plant.rbcL[i] <- length(unique(plant.rbcL4.agg[plant.rbcL4.agg$n.years == time.series4$nb.year[i],]$pop))
# }
# sum(time.series4$plant.rbcL[5:20]) #nb of time series

#gather data summary statistics
summ4 <- data.frame('taxon'=c('mammals','birds','fish','insects','plants_its','plants_matK','plants_rbcL'),
                     'nb.sequences' = numeric(7),
                     'nb.species' = numeric(7),
                     'nb.pops' = numeric(7))
summ4[1,2:4] <- c(sum(mam4.agg$nseqs),length(levels(mam4.agg$species)),length(levels(mam4.agg$pop)))
summ4[2,2:4] <- c(sum(aves4.agg$nseqs),length(levels(aves4.agg$species)),length(levels(aves4.agg$pop)))
summ4[3,2:4] <- c(sum(acti4.agg$nseqs),length(levels(acti4.agg$species)),length(levels(acti4.agg$pop)))
summ4[4,2:4] <- c(sum(insect4.agg$nseqs),length(levels(insect4.agg$species)),length(levels(insect4.agg$pop)))
summ4[5,2:4] <- c(sum(plant.its4.agg$nseqs),length(levels(plant.its4.agg$species)),length(levels(plant.its4.agg$pop)))
summ4[6,2:4] <- c(sum(plant.matK4.agg$nseqs),length(levels(plant.matK4.agg$species)),length(levels(plant.matK4.agg$pop)))
summ4[7,2:4] <- c(sum(plant.rbcL4.agg$nseqs),length(levels(plant.rbcL4.agg$species)),length(levels(plant.rbcL4.agg$pop)))

#cleanup
rm(i,row2rm,acti4,aves4,insect4,mam4,plant.its4,plant.matK4,plant.rbcL4,lu.map)

mam4.agg <- as.data.frame(mam4.agg)
aves4.agg <- as.data.frame(aves4.agg)
acti4.agg <- as.data.frame(acti4.agg)
insect4.agg <- as.data.frame(insect4.agg)
plant.its4.agg <- as.data.frame(plant.its4.agg)
plant.matK4.agg <- as.data.frame(plant.matK4.agg)
plant.rbcL4.agg <- as.data.frame(plant.rbcL4.agg)

#### saving results ####

data.summ <- rbind(summ08,summ1,summ2,summ4)
data.summ$scale <- c(rep(0.08,7),rep(1,7),rep(2,7),rep(4,7))
data.summ <- data.summ[,c(5,1:4)]
write.csv(data.summ,'~/Google Drive/Recherche/Intraspecific genetic diversity/gendiv_data_summary.csv',row.names=F)

rm(summ08,summ1,summ2,summ4)

# ts <- rbind(t(time.series08),t(time.series1)[-1,],t(time.series2)[-1,],t(time.series4)[-1,])
# ts <- ts[,1:15]
# ts2 <- data.frame('scale' = c(rep(0.08,7),rep(1,7),rep(2,7),rep(4,7)),
#                   'taxon' = row.names(ts)[-1])
# ts2 <- cbind(ts2,ts[-1,])
# write.csv(ts2,'~/Google Drive/Recherche/Intraspecific genetic diversity/gendiv_timeseries_details.csv',row.names=F)
#
 #rm(ts,ts2,time.series4,time.series2,time.series1,time.series08)

rm(fname, map.yrs, path, unwanted)

save.image('~/Google Drive/Recherche/Intraspecific genetic diversity/intrasp_gen_div.RData')
