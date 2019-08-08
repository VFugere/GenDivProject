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
    filter(lu.var > 0, hd.var > 0) #this only keeps species with pops that differ in human impact variables, i.e. not all from the same level of impacts
  
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
                  (hd-1|species) +
                  (p.lu-1|species) +
                  (1|year) +
                  (1|order/family) +
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

save(intrasp.models, file = '~/Desktop/spatialGLMMs.Rdata')

#### Sketching these GLMMs ####

viridis(100)[1] -> col_ln

pdf('~/Desktop/spatialGLMMs.pdf',width=4.5,height=8.5,pointsize = 8)
par(mfrow=c(4,2),cex=1,mar=c(4,1,1,1),oma=c(1,2.8,0,0))

ylims <- rbind(c(-8,-3.5),c(-7.5,-4),c(-8.5,-1),c(-7,-2))

for(i in 1:4){
  
  tax <- taxa[i]
  
  spmod<-get(paste('m1',tax,'10',sep='_'))
  
  modsum <- summary(spmod)
  testres <- modsum$coefs[c('hd','p.lu'),c('Estimate','Std. Error')]
  low <- round(testres[,1]-1.96*testres[,2],2)
  high <- round(testres[,1]+1.96*testres[,2],2)
  
  plotfunctions::emptyPlot(xlim = c(0,1),yaxt='n',xaxt='n',ann=F, ylim=ylims[i,],bty='l')
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=seq(-10,0,1))
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(0,0.25,0.5,0.75,1),labels=c('0','0.25','0.5','0.75','1'))
  if(tax == 'mammals'){title(xlab='human density (scaled)',cex=1.2)}
  #title(ylab=expression(log[e]~COI~diversity~(hat(pi))),line=2.8)
  ifelse(tax == 'insects', assign('ln.alpha',0.05), assign('ln.alpha',0.2))
  for(focalsp in levels(spmod$frame$species)){
    spdat <- spmod$frame %>% filter(species == focalsp)
    xs <- seq(min(spdat$hd),max(spdat$hd),by=0.01)
    conditions <- expand.grid(lat.abs=median(spmod$frame$lat.abs),
                       hd=xs,
                       p.lu=median(spmod$frame$p.lu),
                       species=focalsp,
                       order=spdat$order[1],
                       family=spdat$family[1],
                       year=spdat$year[1])
    ys <- predict(spmod, newdata = conditions) 
    points(log(ys)~xs,col=alpha(1,ln.alpha),type='l')
  }
  xs <- seq(0,1,by=0.01)
  ys <- fixef(spmod)[1] + (fixef(spmod)[2]*median(spmod$frame$lat.abs)) + fixef(spmod)[3]*xs + (fixef(spmod)[4]*median(spmod$frame$p.lu))
  lines(ys~xs,lwd=3,col=col_ln)
  
  legend('topright', bty='n', legend = paste0('[',low[1],',',high[1],']'))
  
  ## land use
  
  plotfunctions::emptyPlot(xlim = c(0,1),yaxt='n',xaxt='n',ann=F, ylim=ylims[i,],bty='l')
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=seq(-10,0,1))
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(0,0.25,0.5,0.75,1),labels=c('0','0.25','0.5','0.75','1'))
  if(tax == 'mammals'){title(xlab='land use intensity (scaled)',cex=1.2)}
  #title(ylab=expression(log[e]~COI~diversity~(hat(pi))),line=2.8)
  ifelse(tax == 'insects', assign('ln.alpha',0.05), assign('ln.alpha',0.2))
  for(focalsp in levels(spmod$frame$species)){
    spdat <- spmod$frame %>% filter(species == focalsp)
    xs <- seq(min(spdat$p.lu),max(spdat$p.lu),by=0.01)
    conditions <- expand.grid(lat.abs=median(spmod$frame$lat.abs),
                              hd=median(spmod$frame$hd),
                              p.lu=xs,
                              species=focalsp,
                              order=spdat$order[1],
                              family=spdat$family[1],
                              year=spdat$year[1])
    ys <- predict(spmod, newdata = conditions) 
    points(log(ys)~xs,col=alpha(1,ln.alpha),type='l')
  }
  xs <- seq(0,1,by=0.01)
  ys <- fixef(spmod)[1] + (fixef(spmod)[2]*median(spmod$frame$lat.abs)) + fixef(spmod)[4]*xs + (fixef(spmod)[3]*median(spmod$frame$hd))
  lines(ys~xs,lwd=3,col=col_ln)
  legend('topright', bty='n', legend = paste0('[',low[2],',',high[2],']'))
  
}

mtext(expression(log[e]~COI~diversity~(hat(pi))),at=.5,side=2,outer=T,cex=1.2,line=1)

dev.off()
