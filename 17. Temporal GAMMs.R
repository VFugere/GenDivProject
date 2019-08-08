## Vincent Fugere 2019

# time series GAMMs testing land use impacts on gen div

rm(list=ls())
options(tibble.print_max = 100, scipen = 999)

library(tidyverse)
library(magrittr)
library(mgcv)
library(itsadug)
library(tictoc)
library(scales)
library(rlang)
scale.fun <-function(x){y <- scales::rescale(log1p(x), to = c(0,1)); return(y)}

load('~/Google Drive/Recherche/Intraspecific genetic diversity/Data/DF_Master.RData')

#parameters
min.nb.seqs <- 2
min.nb.years <- 4
taxa <- c('birds','fish','insects','mammals')
#scales <- c('10','100','1000','10000')
scl <- '1000'

# #how many time series of each taxa at each scale?
# table <- data.frame()
# for(y in 3:10){
#   for(ns in 2:10){
#     var_name <- paste0('yrs',quo_name(y),'_seqs',quo_name(ns))
#     t1 <- DF %>% filter(nseqs >= ns, div < mean(div)+10*sd(div)) %>%
#       add_count(pop) %>%
#       filter(n >= y) %>%
#       group_by(taxon,scale) %>% summarize(!! var_name := n_distinct(pop))
#     if(y == 3 & ns == 2){
#       table <- t1}else{
#         table <- left_join(table,t1, by = c('taxon','scale'))
#       }
#   }
# }
# table[is.na(table)] <- 0
# writexl::write_xlsx(table, '~/Desktop/table.xlsx')

#Notes on data exploration using various treshold number of years and min number of sequences
#-goal is to have at least 20 time series per taxon
#-for birds, 4 yrs & 2 seqs @ 1000K scale is the most stringent criterion that has 20+ series
#-both fish and birds have < 20 4yr series at 10, 100 scales
#-tested small-scale effect (@100km scl) for mamms + insects and no significant human impacts when seqs = 5 & yrs = 5
#-so minimum scale is 1000K. most stringent criteria we can use at that scale 
#is 4yrs 2seqs. If I increase min seq or min yrs, model does not fit for birds.

models <- list()

for(tax in taxa){
  
  temp <- DF %>% filter(scale == scl, taxon == tax)
  temp %<>% filter(nseqs >= min.nb.seqs, div < mean(div)+10*sd(div)) %>%
    add_count(pop) %>%
    select(-n.years) %>%
    rename('n.years' = n) %>%
    select(taxon:ncomps,n.years,lat:lu.div) %>%
    filter(n.years >= min.nb.years) %>%
    mutate('year' = as.numeric(year))
  
  #removing duplicate pops per species because nb species == nb pops (almost),
  popdir <- temp %>% group_by(species, pop) %>% summarize('seqs' = sum(nseqs), 'yrs' = median(n.years))
  popdir <- popdir %>% arrange(species, desc(yrs)) %>% distinct(species, .keep_all = T) %>% droplevels
  
  temp <- temp %>% filter(pop %in% popdir$pop) %>% droplevels %>%
    mutate_at(vars(D:lu.div), scale.fun) %>%
    mutate('lat.abs' = rescale(abs(lat),to=c(0,1))) %>% 
    mutate_at(vars(pop,species,family,order), as.factor) %>%
    mutate('wts' = log(nseqs)/mean(log(nseqs))) %>%
    as.data.frame
  
  tsmod <- bam(div ~
                 s(lat,long, bs='gp', k = 50) +
                 s(D, k = 8, bs = 'tp') +
                 s(year, k = 8, bs = 'tp') +
                 s(hd, k = 8, bs = 'tp') +
                 s(p.lu, k = 8, bs = 'tp') +
                 ti(year,hd, k = 8) +
                 ti(year,p.lu, k = 8) +
                 s(year,pop, bs = 'fs', k = 5, m = 1) +
                 s(order, bs='re') +
                 s(family, bs='re'),
               data = temp, family = tw, method='fREML', discrete = T, weights = wts)
  
  # summary(tsmod)
  # anova(tsmod)
  # gam.check(tsmod)

  modname <- paste('m1',tax,scl,sep='_')
  assign(modname,tsmod)
  
  mod <- list(get(modname))
  names(mod) <- modname
  models <- append(models,mod)
  
  rm(tsmod,mod,temp)
  
}

save(models, file = '~/Desktop/temporalGAMMs.Rdata')
