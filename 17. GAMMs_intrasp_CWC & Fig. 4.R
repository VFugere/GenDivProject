## Vincent Fugere 2019

# spatial gamms with centering within cluster (CWC) to test whether there is an overall trend
# of decreasing div with lat within lat range of species, regardless of where they
# are found. I.e. removing species turnover along lat gradient
# interesting: http://philippmasur.de/blog/2018/05/23/how-to-center-in-multilevel-models/

rm(list=ls())
options(tibble.print_max = 100, scipen = 999)

library(tidyverse)
library(magrittr)
library(mgcv)
library(itsadug)
library(scales)
library(RColorBrewer)
library(viridis)
scale.fun <-function(x){y <- scales::rescale(log1p(x), to = c(0,1)); return(y)}
scl1 <-function(x){y <- scales::rescale(x, to = c(0,1)); return(y)}
scl2 <- function(x){y <- scale(x, scale=F); return(y)}

load('~/Google Drive/Recherche/Intraspecific genetic diversity/Data/DF_Master.RData')

#parameters
min.nb.seqs <- 2
min.nb.pops <- 5
taxa <- c('birds','fish','insects','mammals')
#scales <- c('10','100','1000','10000')
scl <- '10'

intrasp.models <- list()

for(tax in taxa){
  
  temp <- DF %>% filter(scale == scl, taxon == tax)
  temp %<>% filter(nseqs >= min.nb.seqs, div < mean(div)+10*sd(div)) %>%
    mutate('year' = as.numeric(year))
  
  #even in insects, the vast majority of populations from species with > 1 pop available
  #have a single data point. Thus, retain one point per pop, but multiple data points per species
  temp %<>% group_by(pop) %>% mutate_at(vars(lat,long), median) %>%
    mutate('year' = ceiling(median(year))) %>%
    mutate_at(vars(div:ncomps,D:lu.div), mean) %>% ungroup %>%
    distinct(pop, .keep_all = T) %>%
    mutate('lat.abs' = abs(lat))
  
  temp <- add_count(temp, species)
  temp <- filter(temp, n >= min.nb.pops)
  
  temp <- temp %>% group_by(species) %>%
    mutate('lu.var' = var(p.lu), 'hd.var' = var(hd), 'lat.var' = max(lat.abs)-min(lat.abs)) %>%
    #mutate('div' = scale(div, scale=T)) %>%
    ungroup %>%
    filter(lu.var > 0, hd.var > 0, lat.var >= 10)
  
  #log-transforming hd+p.lu, then centering within cluster (species), then rescaling from 0 to 1.
  temp <- temp %>% 
    mutate_at(vars(hd,p.lu), log1p) %>%
    group_by(species) %>%
    mutate_at(vars(lat.abs, hd, p.lu), scl2) %>%
    ungroup %>%
    mutate_at(vars(lat.abs, hd, p.lu), scl1) %>%
    mutate_at(vars(pop,species,family,order,year), as.factor) %>%
    droplevels %>%
    mutate('wts' = log(nseqs)/mean(log(nseqs))) %>%
    as.data.frame
  
  mod <- bam(div ~
               #s(lat,long, bs='gp', k = 50) +
               #s(D, k = 8, bs = 'tp') +
               s(lat.abs, k=8) +
               s(hd, k = 8, bs = 'tp') +
               s(p.lu, k = 8, bs = 'tp') +
               s(lat.abs, species, bs = 'fs', k = 6, m = 1) +
               s(hd, species, bs = 'fs', k = 6, m = 1) +
               s(p.lu, species, bs = 'fs', k = 6, m = 1) +
               s(year, bs = 're') +
               s(order, bs='re') +
               s(family, bs='re'),
             data = temp, method='fREML', family='tw',discrete = T, weights = wts, nthreads=2)
  
  # plot(mod)
  # summary(mod, re.test=F)

  # re <- resid(mod)
  # plot(re~temp$lat)
  # plot(re~temp$long) # not including a lat-long smooth does not seem to be a problem
  # 
  # anova(mod)
  # gam.check(mod)

  modname <- paste('m1',tax,scl,sep='_')
  assign(modname,mod)
  
  mod <- list(get(modname))
  names(mod) <- modname
  intrasp.models <- append(intrasp.models,mod)
  
  rm(mod,temp)
  
}

save(intrasp.models, file = '~/Desktop/spatialGAMMs_intrasp_CWC.Rdata')

#### Sketching these GAMMs (Fig. 4) ####

load('~/Desktop/spatialGAMMs_intrasp_CWC.RData')
list2env(intrasp.models,envir=.GlobalEnv)
rm(intrasp.models)

viridis(100)[1] -> col_ln

pdf('~/Desktop/spatialGAMMs_intrasp_CWC.pdf',width=6.75,height=8.5,pointsize = 8)
par(mfrow=c(4,3),cex=1,mar=c(4,1,1,1),oma=c(1,2.8,0,0))

for(i in 1:4){
  
  tax <- taxa[i]
  
  # GAMM fit
  
  spmod<-get(paste('m1',tax,'10',sep='_'))
  
  modsum <- summary(spmod, re.test=F)
  testres <- modsum$s.table[c('s(lat.abs)','s(hd)','s(p.lu)'),c('F','p-value')]
  testres[,1] <- round(testres[,1],2)
  testres[,2] <- round(testres[,2],4)
  testres[testres[,2] == 0,2] <- 0.0001
  rsq <- round(modsum$r.sq,2)
  
  ylims <- predict(spmod,se.fit=T)
  ylims <- range(c(ylims$fit+1.96*ylims$se.fit,ylims$fit-1.96*ylims$se.fit))
  
  #latitude
  emptyPlot(xlim = c(0,1),yaxt='n',xaxt='n',ann=F, ylim=ylims,bty='l')
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=seq(-10,0,1))
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(0,0.25,0.5,0.75,1),labels=c('0','0.25','0.5','0.75','1'))
  if(tax == 'mammals'){title(xlab='absolute latitude (scaled)',cex=1.5)}
  #title(ylab=expression(log[e]~COI~diversity~(hat(pi))),line=2.8)
  ifelse(tax == 'insects', assign('ln.alpha',0.1), assign('ln.alpha',0.2))
  for(focalsp in levels(spmod$model$species)){
    spdat <- spmod$model %>% filter(species == focalsp)
    xs <- seq(min(spdat$lat.abs),max(spdat$lat.abs),by=0.01)
    conditions <- list(lat.abs=xs,
                       hd=median(spmod$model$hd),
                       p.lu=median(spmod$model$p.lu),
                       species=focalsp,
                       order=spdat$order[1],
                       family=spdat$family[1],
                       year=spdat$year[1])
    ys <- get_predictions(spmod, cond = conditions, se=F, print.summary = F, rm.ranef=F) 
    points(fit~lat.abs,ys,col=alpha(1,ln.alpha),type='l')
  }
  plot_smooth(spmod, view="lat.abs", lwd=3, col=col_ln, rm.ranef=T, se=1.96, rug=F, add=T)
  legend('topright',bty='n',legend=bquote(atop(italic('F') == .(testres[1,1]),italic('p') == .(testres[1,2]))))
  
  #HD
  emptyPlot(xlim = c(0,1),yaxt='n',xaxt='n',ann=F, ylim=ylims,bty='l')
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=seq(-10,0,1))
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(0,0.25,0.5,0.75,1),labels=c('0','0.25','0.5','0.75','1'))
  if(tax == 'mammals'){title(xlab='human density (scaled)',cex=1.5)}
  #title(ylab=expression(log[e]~COI~diversity~(hat(pi))),line=2.8)
  ifelse(tax == 'insects', assign('ln.alpha',0.1), assign('ln.alpha',0.2))
  for(focalsp in levels(spmod$model$species)){
    spdat <- spmod$model %>% filter(species == focalsp)
    xs <- seq(min(spdat$hd),max(spdat$hd),by=0.01)
    conditions <- list(lat.abs=median(spmod$model$lat.abs),
                       hd=xs,
                       p.lu=median(spmod$model$p.lu),
                       species=focalsp,
                       order=spdat$order[1],
                       family=spdat$family[1],
                       year=spdat$year[1])
    ys <- get_predictions(spmod, cond = conditions, se=F, print.summary = F, rm.ranef=F) 
    points(fit~hd,ys,col=alpha(1,ln.alpha),type='l')
  }
  plot_smooth(spmod, view="hd", lwd=3, col=col_ln, rm.ranef=T, se=1.96, rug=F, add=T)
  legend('topright',bty='n',legend=bquote(atop(italic('F') == .(testres[2,1]),italic('p') == .(testres[2,2]))))
  #legend('bottomright',bty='n',legend=bquote(model~italic(R)^2 == .(rsq)))
  
  ## land use
  
  emptyPlot(xlim = c(0,1),yaxt='n',xaxt='n',ann=F, ylim=ylims,bty='l')
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=seq(-10,0,1))
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(0,0.25,0.5,0.75,1),labels=c('0','0.25','0.5','0.75','1'))
  if(tax == 'mammals'){title(xlab='land use intensity (scaled)',cex=1.5)}
  #title(ylab=expression(log[e]~COI~diversity~(hat(pi))),line=2.8)
  ifelse(tax == 'insects', assign('ln.alpha',0.1), assign('ln.alpha',0.2))
  for(focalsp in levels(spmod$model$species)){
    spdat <- spmod$model %>% filter(species == focalsp)
    xs <- seq(min(spdat$p.lu),max(spdat$p.lu),by=0.01)
    conditions <- list(lat.abs=median(spmod$model$lat.abs),
                       hd=median(spmod$model$hd),
                       p.lu=xs,
                       species=focalsp,
                       order=spdat$order[1],
                       family=spdat$family[1],
                       year=spdat$year[1])
    ys <- get_predictions(spmod, cond = conditions, se=F, print.summary = F, rm.ranef=F) 
    points(fit~p.lu,ys,col=alpha(1,ln.alpha),type='l')
  }
  plot_smooth(spmod, view="p.lu", lwd=3, col=col_ln, rm.ranef=T, se=1.96, rug=F, add=T)
  legend('topright',bty='n',legend=bquote(atop(italic('F') == .(testres[3,1]),italic('p') == .(testres[3,2]))))
  legend('bottomright',bty='n',legend=bquote(model~italic(R)^2 == .(rsq)))
  
}

mtext(expression(log[e]~COI~diversity~(hat(pi))),at=.5,side=2,outer=T,cex=1.2,line=1)

dev.off()
