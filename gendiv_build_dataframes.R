#### Temporal variation in intra-specific neutral genetic diversity across anthromes
#### Gonzalez Lab project - McGill University - 2016-2018
#### Script by Vincent Fugère

## This script takes pairwise sequence comparions and compute a mean genetic
## diversity value for each species x year x grid cell combination (population)
## Also adds land use and human density values to each gen div value,
## using latest land use map available

#### Part I: load and format data

rm(list=ls())
library(tidyverse)
library(magrittr)

#list of special characters which could occur in species names
#things we want to remove (ie species containing these symbols are discarded)
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

## MAMMALS

mam08 <- read.csv('~/Google Drive/recherche/Intraspecific genetic diversity/data/mamm_pairwise_0.08by0.08_anth.csv')

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

#cleanup
rm(i,row2rm,acti08,aves08,insect08,mam08,lu.map)

mam08.agg <- as.data.frame(mam08.agg)
aves08.agg <- as.data.frame(aves08.agg)
acti08.agg <- as.data.frame(acti08.agg)
insect08.agg <- as.data.frame(insect08.agg)

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

#cleanup
rm(i,row2rm,acti1,aves1,insect1,mam1,lu.map)

mam1.agg <- as.data.frame(mam1.agg)
aves1.agg <- as.data.frame(aves1.agg)
acti1.agg <- as.data.frame(acti1.agg)
insect1.agg <- as.data.frame(insect1.agg)

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

#cleanup
rm(i,row2rm,acti2,aves2,insect2,mam2,lu.map)

mam2.agg <- as.data.frame(mam2.agg)
aves2.agg <- as.data.frame(aves2.agg)
acti2.agg <- as.data.frame(acti2.agg)
insect2.agg <- as.data.frame(insect2.agg)

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

#cleanup
rm(i,row2rm,acti4,aves4,insect4,mam4,lu.map)

mam4.agg <- as.data.frame(mam4.agg)
aves4.agg <- as.data.frame(aves4.agg)
acti4.agg <- as.data.frame(acti4.agg)
insect4.agg <- as.data.frame(insect4.agg)

#### saving results ####

rm(fname, map.yrs, path, unwanted, map.yr.idx)

save.image('~/Google Drive/Recherche/Intraspecific genetic diversity/intrasp_gen_div.RData')
