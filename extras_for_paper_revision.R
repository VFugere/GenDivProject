rm(list=ls())

library(tidyverse)
library(magrittr)
library(scales)
scale.fun <-function(x){y <- scales::rescale(log1p(x), to = c(0,1)); return(y)}
library(devtools)
source_url('https://raw.githubusercontent.com/VFugere/Rfuncs/master/utils.R')

load('~/Google Drive/Recherche/Intraspecific genetic diversity/Data/DF_Master.RData')
load('~/Google Drive/Recherche/Intraspecific genetic diversity/Data/sequence_metadata.RData')

#parameters
min.nb.seqs <- 2
taxa <- c('birds','fish','insects','mammals')
scales <- c('10','100','1000','10000')

# how many outliers are excluded?
outliers <- data.frame()
for(tax in taxa){
  for(scl in scales){
    temp <- DF %>% filter(scale == scl, taxon == tax)
    temp2 <- filter(temp, div < mean(div)+10*sd(div))
    res <- temp[1,1:2]
    res[1,3] <- nrow(temp)
    res[1,4] <- nrow(temp)-nrow(temp2)
    outliers <- rbind(outliers,res)
  }
}
rm(res,temp,temp2)
colnames(outliers)[3:4] <- c('pis','excluded.pis')
outliers$prop <- (outliers$excluded.pis/outliers$pis)*100
100*(sum(outliers$excluded.pis)/sum(outliers$pis))

#how many pi estimates don't have year-specific HYDE estimates?
yrs <- data.frame('seq.yr' = 1980:2016, 'hyde.yr' = c(rep(1980,10),rep(1990,10),2000:2016))
df.yr <- select(seq, species:year) %>% mutate(year = as.numeric(year))
df.yr <- left_join(df.yr, yrs, by = c('year' = 'seq.yr'))
df.yr %>% filter(year<2000, year %!in% c(1980,1990)) %>% nrow(.)/nrow(seq)
df.yr$lag <- df.yr$year - df.yr$hyde.yr
mean(df.yr$lag)
sd(df.yr$lag)

#how many species have more than one pop per year?
pops <- data.frame()

for(tax in taxa){
  
  for(scl in scales){
    
    temp <- DF %>% filter(scale == scl, taxon == tax)
    temp %<>% filter(nseqs >= min.nb.seqs, div < mean(div)+10*sd(div)) %>%
      mutate('year' = as.numeric(year))
    temp$sp.yr <- paste(temp$species,temp$year,sep='_')
    pops.t <- temp %>% count(sp.yr) %>% filter(n > 1)
    
    results <- data.frame('tax'=tax,'scale'=scl,'n.species'=nrow(pops.t))
    pops <- rbind(pops,results)
    
  }
}

# what is the average latitudinal variation within species with multiple populations of data?

min.nb.pops <- 5
scl <- '10'

pdf('~/Desktop/latranges.pdf',width = 6,height = 6,pointsize = 8)
par(mfrow=c(2,2),cex=1)

for(tax in taxa){
  
  temp <- DF %>% filter(scale == scl, taxon == tax)
  temp %<>% filter(nseqs >= min.nb.seqs, div < mean(div)+10*sd(div)) %>%
    mutate('year' = as.numeric(year))
  
  #retain one point per pop, but multiple data points per species
  temp %<>% group_by(pop) %>% mutate_at(vars(lat,long), median) %>%
    mutate('year' = ceiling(median(year))) %>%
    mutate_at(vars(div:ncomps,D:lu.div), mean) %>% ungroup %>%
    distinct(pop, .keep_all = T) %>%
    mutate('lat.abs' = abs(lat))
  
  temp <- add_count(temp, species)
  temp <- filter(temp, n >= min.nb.pops)
  
  lat.vars <- temp %>% group_by(species) %>% summarize('lat.var' = max(lat.abs)-min(lat.abs))
  hist(lat.vars$lat.var, breaks=50,main=tax,xlab='latitude range across populations',ylab='species')
  
}

dev.off()

# how many species show evidence of change in div within species along lat gradient?

library(cplm)

results <- data.frame()

for(tax in taxa){
  
  temp <- DF %>% filter(scale == scl, taxon == tax)
  temp %<>% filter(nseqs >= min.nb.seqs, div < mean(div)+10*sd(div)) %>%
    mutate('year' = as.numeric(year))
  
  #retain one point per pop, but multiple data points per species
  temp %<>% group_by(pop) %>% mutate_at(vars(lat,long), median) %>%
    mutate('year' = ceiling(median(year))) %>%
    mutate_at(vars(div:ncomps,D:lu.div), mean) %>% ungroup %>%
    distinct(pop, .keep_all = T) %>%
    mutate('lat.abs' = abs(lat))
  
  temp <- add_count(temp, species)
  temp <- filter(temp, n >= min.nb.pops)
  
  temp <- temp %>% group_by(species) %>%
    mutate('lu.var' = var(p.lu), 'hd.var' = var(hd), 'lat.var' = max(lat.abs)-min(lat.abs)) %>%
    #mutate('div' = scale(div, scale=T)) %>%
    ungroup %>%
    filter(lu.var > 0, hd.var > 0, lat.var >= 10)
  
  temp$species <-as.factor(temp$species)
  
  for(focalsp in levels(temp$species)){
    spdat <- filter(temp, species == focalsp)
    mod <- cpglm(div~lat.abs, data=spdat)
    rez <- as.data.frame(summary(mod)$coefficients)
    rez <- rez[2,c(1,4)]
    rez <- cbind(rez,tax,focalsp,nrow(spdat))
    results <- rbind(results,rez)
  }
  
}

