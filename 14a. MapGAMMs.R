## Vincent Fugere 2019

# GAMMs to model/map gen div after removing effect of spatial distance, taxonomy, and sampling year

rm(list=ls())
options(tibble.print_max = 100, scipen = 999)

library(tidyverse)
library(magrittr)
library(mgcv)
library(scales)
scale.fun <-function(x){y <- scales::rescale(log1p(x), to = c(0,1)); return(y)}

load('~/Google Drive/Recherche/Intraspecific genetic diversity/Data/DF_Master.RData')

#parameters
min.nb.seqs <- 2
taxa <- c('birds','fish','insects','mammals')
scales <- c('10','100','1000','10000')

map.models <- list()

scl <- scales[3]

for(tax in taxa){
  
    temp <- DF %>% filter(scale == scl, taxon == tax)
    temp %<>% filter(nseqs >= min.nb.seqs, div < mean(div)+10*sd(div)) %>%
      mutate('year' = as.numeric(year))
    temp %<>% group_by(pop) %>% mutate_at(vars(lat,long), median) %>%
      mutate('year' = ceiling(median(year))) %>%
      mutate_at(vars(div:ncomps,D:lu.div), mean) %>% ungroup %>%
      distinct(pop, .keep_all = T)
    temp <- temp %>% mutate_at(vars(D:lu.div), scale.fun) %>%
      mutate('lat.abs' = rescale(abs(lat),to=c(0,1))) %>% 
      mutate_at(vars(pop,year,species, family, order), as.factor)
    
    mapmod <- bam(div ~ s(D, k = 10, bs = 'tp') + s(year,bs = 're', k = 5, m=1) + s(order, bs='re',k = 5, m=1) + s(family, bs='re',k = 5, m=1), data = temp, family = tw, method='fREML', discrete = T)
    #p <- str_split(family(mapmod)[[1]], '=', simplify = T)[1,2] %>% str_remove('\\)') %>% as.numeric
    
    modname <- paste('m0',tax,scl,sep='_')
    assign(modname,mapmod)
    
    mod <- list(get(modname))
    names(mod) <- modname
    map.models <- append(map.models,mod)
    
    rm(mapmod,temp,mod)
  
}

save(map.models, file = '~/Google Drive/Recherche/Intraspecific genetic diversity/Data/mapGAMMs.Rdata')
