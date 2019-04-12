## Vincent Fugere 2019

# spatial GAMM testing land use impacts on gen div

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
taxa <- c('birds','fish','insects','mammals')
scales <- c('10','100','1000','10000')

models <- list()

tax <- taxa[1]

for(scl in scales){
  
  temp <- DF %>% filter(scale == scl, taxon == tax)
  temp %<>% filter(nseqs >= min.nb.seqs, div < mean(div)+10*sd(div)) %>%
    mutate('year' = as.numeric(year))
  temp %<>% group_by(pop) %>% mutate_at(vars(lat,long), median) %>%
    mutate('year' = ceiling(median(year))) %>%
    mutate_at(vars(div:ncomps,D:lu.div), mean) %>% ungroup %>%
    distinct(pop, .keep_all = T)
  temp <- temp %>% mutate_at(vars(D:lu.div), scale.fun) %>%
    mutate('lat.abs' = rescale(abs(lat),to=c(0,1)), 'wts' = log(nseqs)/mean(log(nseqs))) %>%
    select(-taxon,-scale,-nseqs, -ncomps,-n.years, -lu.var) %>%
    mutate_at(vars(pop,year,species), as.factor) %>% as.data.frame
  
  pmod <- bam(div ~ s(lat,long, bs='gp', k = 50, m = c(1,.5)) + s(D, k = 10, bs = 'tp'), data = temp, family = tw, method='fREML', discrete = T)
  p <- str_split(family(pmod)[[1]], '=', simplify = T)[1,2] %>% str_remove('\\)') %>% as.numeric
  
  fullmod <- bam(div ~ s(lat,long, bs='tp', k = 15) + s(D, k = 10, bs = 'tp') + s(year,bs = 're', k = 5, m=1) +
                   s(lat.abs, k = 15, bs = 'tp') + s(hd, k=10, bs ='tp') +
                   s(p.lu, k=10, bs='tp') + s(species, bs='re', k = 5, m=1),
                 data = temp, family = Tweedie(p = p), method='fREML', discrete = T, weights = wts)
  
  modnames <- paste(c('m0','m1'),tax,scl,sep='_')
  
  assign(modnames[1],pmod)
  assign(modnames[2],fullmod)
  
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


