## Vincent Fugere 2019

# time series GAMMs testing land use impacts on gen div

rm(list=ls())
options(tibble.print_max = 100, scipen = 999)

library(tidyverse)
library(magrittr)
library(cplm)
library(scales)
library(RColorBrewer)
library(viridis)
scale.fun <-function(x){y <- scales::rescale(log1p(x), to = c(0,1)); return(y)}
scale.fun2 <-function(x){y <- scales::rescale(x, to = c(0,1)); return(y)}

load('~/Google Drive/Recherche/Intraspecific genetic diversity/Data/DF_Master.RData')

#parameters
min.nb.seqs <- 2
min.nb.years <- 4
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
  
  temp <- temp %>%
    filter(pop %in% popdir$pop) %>%
    droplevels %>%
    mutate_at(vars(D:lu.div), scale.fun) %>%
    mutate('lat.abs' = rescale(abs(lat),to=c(0,1))) %>%
    mutate_at(vars(year,lat,long), scale.fun2) %>% 
    mutate_at(vars(pop,species,family,order), as.factor) %>%
    mutate('wts' = log(nseqs)/mean(log(nseqs))) %>%
    as.data.frame
  
  tsmod <- cpglmm(div ~ 
                    lat*long +
                    D +
                    year +
                    hd +
                    p.lu +
                    (1+year|species) +
                    (1|order/family),
                  data = temp,
                  weights = wts,
                  control=list('max.iter' = 100000))
                    
  modname <- paste('m1',tax,scl,sep='_')
  assign(modname,tsmod)
  
  mod <- list(get(modname))
  names(mod) <- modname
  models <- append(models,mod)
  
  rm(tsmod,mod,temp)
  
}

save(models, file = '~/Desktop/temporalGLMMs.Rdata')

#### Sketching these GLMMs ####

pdf('~/Desktop/temporalGLMMs.pdf',width=5,height=5,pointsize = 8)

viridis(100)[1] -> col_ln
par(mfrow=c(2,2),cex=1,mar=c(2,2,1,1),oma=c(3,3,0,0))

xlims <- rbind(c(1985,2015),c(1994,2016),c(1988,2016),c(1984,2015))
ylims <- rbind(c(-9,-4),c(-7,-3),c(-12,-4),c(-7,-3))

for(i in 1:4){
  
  tax <- taxa[i]
  
  spmod<-get(paste('m1',tax,'1000',sep='_'))
  
  modsum <- summary(spmod)
  testres <- modsum$coefs[c('year'),c('Estimate','Std. Error')]
  low <- round(testres[1]-1.96*testres[2],2)
  high <- round(testres[1]+1.96*testres[2],2)
  
  plotfunctions::emptyPlot(xlim = range(spmod$frame$year),yaxt='n',xaxt='n',ann=F, ylim=ylims[i,],bty='l')
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=seq(-12,0,1))
  xlabs <- xlims[i,]
  xlabs <- round(seq(xlabs[1],xlabs[2],length.out = 3))
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(0,0.5,1),labels = as.character(xlabs))
  ifelse(tax == 'insects', assign('ln.alpha',0.05), assign('ln.alpha',0.2))
  for(focalsp in levels(spmod$frame$species)){
    spdat <- spmod$frame %>% filter(species == focalsp)
    xs <- seq(min(spdat$year),max(spdat$year),by=0.01)
    conditions <- expand.grid(lat=0,
                              long=0,
                              D=0,
                              hd=0,
                              p.lu=0,
                              year=xs,
                              species=focalsp,
                              order=spdat$order[1],
                              family=spdat$family[1])
    ys <- predict(spmod, newdata = conditions) 
    points(log(ys)~xs,col=alpha(1,ln.alpha),type='l')
  }
  xs <- seq(min(spmod$frame$year),max(spmod$frame$year),by=0.01)
  ys <- fixef(spmod)[1] + fixef(spmod)[5]*xs
  lines(ys~xs,lwd=3,col=col_ln)
  legend('topright', bty='n', legend = paste0('[',low,',',high,']'))
}

mtext(expression(log[e]~COI~diversity~(hat(pi))),at=.5,side=2,outer=T,cex=1.2,line=1)
mtext('year',at=.5,side=1,outer=T,cex=1.2,line=1)

dev.off()

