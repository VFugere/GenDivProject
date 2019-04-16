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
min.nb.years <- 5
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
  
  tsmod2 <- bam(div ~ s(lat,long, bs='gp', k = 50) + s(D, k = 8, bs = 'tp') + s(year, k = 8, bs = 'tp') +
                 s(hd, k = 6, bs = 'tp') + s(p.lu, k = 6, bs = 'tp') +
                 ti(year,hd, k = 6) + ti(year,p.lu, k = 6) + s(year,pop, bs = 'fs', k = 5, m = 1),
               data = temp, family = tw, method='fREML', discrete = T, weights = wts)
  
  summary(tsmod2)
  anova(tsmod2)
  gam.check(tsmod2)
  
  fvisgam(tsmod2, view = c('year','p.lu'), cond = list('hd' = 0), ylab='land use',add.color.legend=F,hide.label=T,xlab = 'year',plot.type = 'contour', color = viridis(50), main = NULL)
  pvisgam(tsmod2, view = c('year','p.lu'), ylab='land use',add.color.legend=F,hide.label=T,xlab = 'year',plot.type = 'contour', color = viridis(50), main = NULL)
  pvisgam(tsmod2, view = c('year','hd'), ylab='land use',add.color.legend=F,hide.label=T,xlab = 'year',plot.type = 'contour', color = viridis(50), main = NULL)
  
  plot_smooth(tsmod2, view="year", cond=list('hd' = 0,'p.lu'= 0.1), col=1, rm.ranef=T, se=1.96, yaxt='n',xaxt='n',ann=F, hide.label = T,lwd=2,main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA)
  plot_smooth(tsmod2, view="year", cond=list('hd' = 0,'p.lu' = 0.9), lwd=2,col='darkorchid1', add=T,rm.ranef=T, se=1.96,rug=F)
  
  inspect_random(tsmod2)
  
  pvisgam(tsmod2, view = c('year','hd'), cond = list('p.lu' = 0), ylab='human density',add.color.legend=F,hide.label=T,xlab = 'year',plot.type = 'contour', color = viridis(50), main = NULL)
  plot_smooth(tsmod2, view="year", cond=list('hd' = 0.1,'p.lu'= 0), lwd=2,col=1, rm.ranef=T, se=1.96, yaxt='n',xaxt='n',ann=F, hide.label = T,main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA)
  plot_smooth(tsmod2, view="year", cond=list('hd' = 0.9,'p.lu' = 0), lwd=2,col='darkorange1', add=T,rm.ranef=T, se=1.96,rug=F)
  
  viridis(50)[1] -> col_ln
  plot_smooth(tsmod2, view="year", plot_all=c('pop'),col=alpha(1,0.2), rm.ranef=F, se=0, yaxt='n',xaxt='n',ann=F, hide.label = T,main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA)
  plot_smooth(tsmod2, view="year", lwd=3, col=col_ln, rm.ranef=T, se=1.96, rug=F, add=T)
  
  plot_smooth(tsmod2, view="year", lwd=3, col=col_ln, rm.ranef=T, se=1.96, rug=F, add=F)
  
  plot_smooth(tsmod2, view="year", plot_all=c('pop'),add=T,col=alpha(1,0.2), rm.ranef=F, se=0, yaxt='n',xaxt='n',ann=F, hide.label = T,main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA)
  
  
  modname <- paste('m1',tax,scl,sep='_')
  assign(modname,tsmod)
  
  mod <- list(get(modname))
  names(mod) <- modname
  models <- append(models,mod)
  
  rm(tsmod,mod,temp)
  
}

save(models, file = '~/Desktop/temporalGAMMs.Rdata')
