## Vincent Fugere 2019

# time series GAMMs testing land use impacts on gen div

rm(list=ls())
options(tibble.print_max = 100, scipen = 999)

library(tidyverse)
library(magrittr)
library(mgcv)
library(tictoc)
library(scales)
scale.fun <-function(x){y <- scales::rescale(log1p(x), to = c(0,1)); return(y)}

load('~/Google Drive/Recherche/Intraspecific genetic diversity/Data/DF_Master.RData')

#parameters
min.nb.seqs <- 2
min.nb.years <- 3
taxa <- c('birds','fish','insects','mammals')
#scales <- c('10','100','1000','10000')
scl <- '1000'

models <- list()

for(tax in taxa){
  
  temp <- DF %>% filter(scale == scl, taxon == tax)
  temp %<>% filter(nseqs >= min.nb.seqs, div < mean(div)+10*sd(div), n.years >= min.nb.years) %>%
    mutate('year' = as.numeric(year))
  
  #removing duplicate pops per species because nb species == nb pops (almost),
  popdir <- temp %>% group_by(species, pop) %>% summarize('seqs' = sum(nseqs), 'yrs' = median(n.years))
  popdir <- popdir %>% arrange(species, desc(yrs)) %>% distinct(species, .keep_all = T) %>% droplevels
  
  temp <- temp %>% filter(pop %in% popdir$pop) %>% droplevels %>%
    mutate_at(vars(D:lu.div), scale.fun) %>%
    mutate('lat.abs' = rescale(abs(lat),to=c(0,1))) %>% 
    mutate_at(vars(pop,species), as.factor) %>%
    mutate('wts' = log(nseqs)/mean(log(nseqs))) %>%
    as.data.frame
  
  tsmod <- bam(div ~ s(lat,long, bs='gp', k = 50) + s(D, k = 8, bs = 'tp') + s(year, k = 8, bs = 'tp') +
                 ti(year,hd, k = 8) + ti(year,p.lu, k = 8) + s(year,pop, bs = 'fs', k = 5, m = 1),
               data = temp, family = tw, method='fREML', discrete = T, weights = wts)
  
  # summary(tsmod)
  # gam.check(tsmod)
  
  modname <- paste('m1',tax,scl,sep='_')
  assign(modname,tsmod)
  
  mod <- list(get(modname))
  names(mod) <- modname
  models <- append(models,mod)
  
  rm(tsmod,mod,temp)
  
}

save(models, file = '~/Desktop/temporalGAMMs.Rdata')
