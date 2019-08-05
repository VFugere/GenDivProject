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