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
  
  summary(tsmod)
  gam.check(tsmod)
  library(viridis)
  vis.gam(tsmod, view = c('year','hd'), cond = list('lat' = 0, 'long' = 0, 'p.lu' = 0), plot.type = 'contour', type = 'response', color = 'bw', main = NULL)
  fvisgam(tsmod, view = c('year','hd'), cond = list('lat' = 0, 'long' = 0, 'p.lu' = 0), ylab='human density',add.color.legend=F,hide.label=T,xlab = 'year',transform=exp,plot.type = 'contour', color = viridis(50), main = NULL)
  
  modname <- paste('m1',tax,scl,sep='_')
  assign(modname,tsmod)
  
  mods <- list(get(modnames[1]),get(modnames[2]))
  names(mods) <- modnames
  models <- append(models,mods)
  
  rm(pmod,fullmod,temp,mods)
  
}

save(models, file = '~/Desktop/spatialGAMMs.Rdata')

# summary(m1)
# gam.check(m1)
# 
# 
# #spatial autocorrelation in residuals?
# temp$R <- resid(m1)
# temp$rcol <- 4
# temp[temp$R < 0, 'rcol'] <- 2
# plot(map, xlim = c(-180,180), ylim = c(-90,90),border=NA,col='grey95',axes=F)
# points(lat~long,temp,pch=16,col=alpha(temp$rcol,0.5),cex=0.5+abs(temp$R))
# 
# spatdat <- select(temp, long,lat,R)
# coordinates(spatdat) <- c('long','lat')
# variogram(R~long+lat,spatdat,width=0.1,cutoff=50) -> var1
# scatter.smooth(x=var1$dist,y=var1$gamma,xlab='distance',ylab='semivariance',ylim=c(0,max(var1$gamma)),pch=16,col='gray')


