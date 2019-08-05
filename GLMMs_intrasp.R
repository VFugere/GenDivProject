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

load('~/Google Drive/Recherche/Intraspecific genetic diversity/Data/DF_Master.RData')

#parameters
min.nb.seqs <- 2
min.nb.pops <- 4
taxa <- c('birds','fish','insects','mammals')
#scales <- c('10','100','1000','10000')
scl <- '10'

intrasp.models <- list()

for(tax in taxa){
  
  temp <- DF %>% filter(scale == scl, taxon == tax)
  temp %<>% filter(nseqs >= min.nb.seqs, div < mean(div)+10*sd(div)) %>%
    mutate('year' = as.numeric(year))
  temp$sp.yr <- paste(temp$species,temp$year,sep='_')
  temp <- add_count(temp, sp.yr)
  temp <- filter(temp, n >= min.nb.pops)
  
  #even in insects, the vast majority of populations from species with > 1 pop available
  #have a single data point. Thus, retain one point per pop, but multiple data points per species
  temp %<>% group_by(pop) %>% mutate_at(vars(lat,long), median) %>%
    mutate('year' = ceiling(median(year))) %>%
    mutate_at(vars(div:ncomps,D:lu.div), mean) %>% ungroup %>%
    distinct(pop, .keep_all = T)
  
  temp <- temp %>% group_by(species) %>%
    mutate('lu.var' = var(p.lu), 'hd.var' = var(hd)) %>% 
    ungroup %>%
    filter(lu.var > 0, hd.var > 0)
  
  temp <- temp %>%
    mutate_at(vars(D:lu.div), scale.fun) %>%
    mutate('lat.abs' = rescale(abs(lat),to=c(0,1))) %>% 
    mutate_at(vars(sp.yr,pop,species,family,order,year), as.factor) %>%
    droplevels %>%
    mutate('wts' = log(nseqs)/mean(log(nseqs))) %>%
    as.data.frame
  
  mod <- cpglmm(div ~
               lat.abs +
               hd +
               p.lu +
                 (hd|species) +
                 (p.lu|species) +
                 (1|year) +
                 (1|family) +
                 (1|order) +
                 (1|species),
             data = temp, weights = wts)
  
  # plot(mod)
  # summary(mod)
  # 
  # re <- resid(mod)
  # plot(re~temp$lat)
  # plot(re~temp$long) # not including a lat-long smooth does not seem to be a problem
  # 
  
  modname <- paste('m1',tax,scl,sep='_')
  assign(modname,mod)
  
  mod <- list(get(modname))
  names(mod) <- modname
  intrasp.models <- append(intrasp.models,mod)
  
  rm(mod,temp)
  
}

save(intrasp.models, file = '~/Desktop/GLMMs_intrasp.Rdata')

#### Sketching these GAMMs ####

viridis(100)[1] -> col_ln

pdf('~/Desktop/GLMMs_intrasp.pdf',width=4.5,height=8.5,pointsize = 8)
par(mfrow=c(4,2),cex=1,mar=c(4,1,1,1),oma=c(1,2.8,0,0))

for(i in 1:4){
  
  tax <- taxa[i]
  
  spmod<-get(paste('m1',tax,'10',sep='_'))
  
  modsum <- summary(spmod)
  testres <- modsum$coefs[c('hd','p.lu'),c('t value')]
  testres <- round(testres,2)
  
  emptyPlot(xlim = c(0,1),yaxt='n',xaxt='n',ann=F, ylim=c(-9,-2),bty='l')
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=seq(-10,0,1))
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(0,0.25,0.5,0.75,1),labels=c('0','0.25','0.5','0.75','1'))
  if(tax == 'mammals'){title(xlab='human density (scaled)',cex=1.2)}
  #title(ylab=expression(log[e]~COI~diversity~(hat(pi))),line=2.8)
  ifelse(tax == 'insects', assign('ln.alpha',0.05), assign('ln.alpha',0.2))
  for(focalsp in levels(spmod$frame$species)){
    spdat <- spmod$frame %>% filter(species == focalsp)
    xs <- seq(min(spdat$hd),max(spdat$hd),by=0.01)
    conditions <- expand.grid(lat.abs=0,
                       hd=xs,
                       p.lu=0,
                       species=focalsp,
                       order=spdat$order[1],
                       family=spdat$family[1],
                       year=spdat$year[1])
    ys <- predict(spmod, newdata = conditions) 
    points(log(ys)~xs,col=alpha(1,ln.alpha),type='l')
  }
  xs <- seq(0,1,by=0.01)
  ys <- fixef(spmod)[1] + fixef(spmod)[3]*xs
  lines(ys~xs,lwd=3,col=col_ln)
  mtext(text=bquote(italic('t') == .(testres[1])),side=3,adj=1)
  
  ## land use
  
  emptyPlot(xlim = c(0,1),yaxt='n',xaxt='n',ann=F, ylim=c(-9,-2),bty='l')
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=seq(-10,0,1))
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(0,0.25,0.5,0.75,1),labels=c('0','0.25','0.5','0.75','1'))
  if(tax == 'mammals'){title(xlab='land use intensity (scaled)',cex=1.2)}
  #title(ylab=expression(log[e]~COI~diversity~(hat(pi))),line=2.8)
  ifelse(tax == 'insects', assign('ln.alpha',0.05), assign('ln.alpha',0.2))
  for(focalsp in levels(spmod$frame$species)){
    spdat <- spmod$frame %>% filter(species == focalsp)
    xs <- seq(min(spdat$p.lu),max(spdat$p.lu),by=0.01)
    conditions <- expand.grid(lat.abs=0,
                              hd=0,
                              p.lu=xs,
                              species=focalsp,
                              order=spdat$order[1],
                              family=spdat$family[1],
                              year=spdat$year[1])
    ys <- predict(spmod, newdata = conditions) 
    points(log(ys)~xs,col=alpha(1,ln.alpha),type='l')
  }
  xs <- seq(0,1,by=0.01)
  ys <- fixef(spmod)[1] + fixef(spmod)[4]*xs
  lines(ys~xs,lwd=3,col=col_ln)
  mtext(text=bquote(italic('t') == .(testres[2])),side=3,adj=1)
  
}

mtext(expression(log[e]~COI~diversity~(hat(pi))),at=.5,side=2,outer=T,cex=1.2,line=1)

dev.off()
