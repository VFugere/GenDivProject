## Vincent Fugere 2019

# This code unites population-level datasets and create a master data file for analysis
# Also explores the number of populations and times series at each scale

rm(list=ls())
library(tidyverse)
library(vegan)

load('~/Google Drive/Recherche/Intraspecific genetic diversity/Data/DF_D.Rdata')
load('~/Google Drive/Recherche/Intraspecific genetic diversity/Data/DF_Pi.Rdata')
load('~/Google Drive/Recherche/Intraspecific genetic diversity/Data/sequence_metadata.RData')

DF_Pi <- select(DF_Pi, -(taxon:species))
DF_Pi$year <- as.numeric(DF_Pi$year)

alldata <- inner_join(DF_D, DF_Pi, by = c('pop','year'))

#how many pops per scale?
alldata %>% group_by(scale) %>% summarize(pops = n_distinct(pop))
#10K vs. no treshold have almost the same number of pops (32 pops diff)
#Biggest jump is from 1K to 100, loosing 5K pops

alldata <- alldata %>% group_by(pop) %>% add_tally(name = 'n.years')

#how many time series?
alldata %>% filter(n.years >= 4) %>%
  group_by(scale, taxon) %>%
  summarize('pops' = n_distinct(pop), 'species' = n_distinct(species))
#plenty! and number of species ~ number of pops even at smallest scales, indicating
#that we are not simply cutting pops into small units

#average great circle distance within species, for figure S1
dd <- DF_D %>% filter(scale == '1e+05') %>% group_by(species) %>% summarize('D' = mean(D,na.rm=T))
#pdf('~/Desktop/S1.pdf',width=4.5,height=3.5,pointsize = 10)
hist(dd$D, breaks=100, xlab = 'mean Great Circle distance among sequences (km)',
     ylab='species',col='grey50',border=0,main=NULL)
#dev.off()
rm(dd)

#how many sequences per pi hat?
#pdf('~/Desktop/S6.pdf',width=4.5,height=3.5,pointsize = 10)
hist(log10(DF_Pi$nseqs), breaks=100, xlab=number~of~sequences~compared,
     ylab=expression(number~of~hat(pi)~estimates),col='grey50',border=0,main=NULL,axes=F)
box(bty='l')
axis(1,lwd=0,lwd.tick=1,at=log10(c(2,5,10,50,100,1000)),labels = c(2,5,10,50,100,1000))
axis(2,lwd=0,lwd.tick=1)
#dev.off()

#get spatial centroids of pops
sc10 <- seq %>% group_by(pop10) %>% summarize('nseqs' = n(),'lat' = mean(lat), 'long' = mean(long)) %>% rename('pop' = pop10)
sc100 <- seq %>% group_by(pop100) %>% summarize('nseqs' = n(),'lat' = mean(lat), 'long' = mean(long)) %>% rename('pop' = pop100)
sc1000 <- seq %>% group_by(pop1000) %>% summarize('nseqs' = n(),'lat' = mean(lat), 'long' = mean(long)) %>% rename('pop' = pop1000)
sc10000 <- seq %>% group_by(pop10000) %>% summarize('nseqs' = n(),'lat' = mean(lat), 'long' = mean(long)) %>% rename('pop' = pop10000)
sc100000 <- seq %>% group_by(pop100000) %>% summarize('nseqs' = n(),'lat' = mean(lat), 'long' = mean(long)) %>% rename('pop' = pop100000)
centroids <- bind_rows(sc10,sc100,sc1000,sc10000,sc100000) %>% filter(nseqs > 1) %>% select(-nseqs)
rm(sc10,sc100,sc1000,sc10000,sc100000)
alldata <- left_join(alldata, centroids, by = 'pop')

#add human impact variables

seq <- mutate(seq, 'low.impact' = 1-(urban+cropland+conv_rangeland+pasture), 'managed.grazing' = conv_rangeland+pasture)

lu.div.func <- function(x1,x2,x3,x4){
  y <- cbind(x1,x2,x3,x4)
  z <- dist(y, method = "euclidean", diag = FALSE, upper = FALSE)
  q <- mean(z)
  return(q)
}

sc10 <- seq %>% group_by(pop10,year) %>% summarize('nseqs' = n(), 'hd' = mean(pop_tot), 'hd.var' = var(pop_tot), 'p.lu' = mean(1-low.impact,na.rm = T), 'lu.var' = var(1-low.impact,na.rm = T), 'lu.div' = lu.div.func(urban,low.impact,managed.grazing,cropland)) %>% rename('pop' = pop10)
sc100 <- seq %>% group_by(pop100,year) %>% summarize('nseqs' = n(), 'hd' = mean(pop_tot), 'hd.var' = var(pop_tot), 'p.lu' = mean(1-low.impact,na.rm = T), 'lu.var' = var(1-low.impact,na.rm = T), 'lu.div' = lu.div.func(urban,low.impact,managed.grazing,cropland)) %>% rename('pop' = pop100)
sc1000 <- seq %>% group_by(pop1000,year) %>% summarize('nseqs' = n(), 'hd' = mean(pop_tot), 'hd.var' = var(pop_tot), 'p.lu' = mean(1-low.impact,na.rm = T), 'lu.var' = var(1-low.impact,na.rm = T), 'lu.div' = lu.div.func(urban,low.impact,managed.grazing,cropland)) %>% rename('pop' = pop1000)
sc10000 <- seq %>% group_by(pop10000,year) %>% summarize('nseqs' = n(), 'hd' = mean(pop_tot), 'hd.var' = var(pop_tot), 'p.lu' = mean(1-low.impact,na.rm = T), 'lu.var' = var(1-low.impact,na.rm = T), 'lu.div' = lu.div.func(urban,low.impact,managed.grazing,cropland)) %>% rename('pop' = pop10000)
sc100000 <- seq %>% group_by(pop100000,year) %>% summarize('nseqs' = n(), 'hd' = mean(pop_tot), 'hd.var' = var(pop_tot), 'p.lu' = mean(1-low.impact,na.rm = T), 'lu.var' = var(1-low.impact,na.rm = T), 'lu.div' = lu.div.func(urban,low.impact,managed.grazing,cropland)) %>% rename('pop' = pop100000)
human.impacts <- bind_rows(sc10,sc100,sc1000,sc10000,sc100000) %>% filter(nseqs > 1) %>% select(-nseqs) %>% mutate('year' = as.numeric(year))
rm(sc10,sc100,sc1000,sc10000,sc100000)
alldata <- left_join(alldata, human.impacts, by = c('pop','year'))

DF <- alldata %>% ungroup %>% arrange(taxon,scale,year,species,pop) %>%
  select(taxon,scale,pop,year,species,div:long,D,hd:lu.div) %>% as.data.frame

#adding family + order
taxonomy <- seq %>% select(species,family,order) %>% distinct(species, .keep_all = T)
rm(seq)
taxonomy$species <- str_replace(taxonomy$species, ' ', '_')
DF <- left_join(DF,taxonomy, by = 'species')
DF <- select(DF, taxon,scale,pop:species,family,order,div:lu.div)

save(DF, file = '~/Google Drive/Recherche/Intraspecific genetic diversity/Data/DF_Master.RData')
