## Vincent Fugere 2019

# spatial GAMM testing whether north vs. south hemisphere affect land use impacts on gen div

rm(list=ls())
options(tibble.print_max = 100, scipen = 999)

library(tidyverse)
library(magrittr)
library(mgcv)
library(RColorBrewer)
library(scales)
library(viridis)
library(itsadug)
scale.fun <-function(x){y <- scales::rescale(log1p(x), to = c(0,1)); return(y)}

load('~/Google Drive/Recherche/Intraspecific genetic diversity/Data/DF_Master.RData')

#parameters
min.nb.seqs <- 2
taxa <- c('birds','fish','insects','mammals')
scales <- c('10','100','1000','10000')

#adding hemisphere
DF$zone <- 'Temperate'
DF$zone[DF$lat > -30 & DF$lat < 30] <- 'Tropical'
DF$zone <- as.ordered(DF$zone)
#lat > 60 degrees excluded below

## GAMMS

tax<-taxa[4]
scl<-scales[2]

for(tax in taxa){
  
  for(scl in scales){
    
    temp <- DF %>% filter(scale == scl, taxon == tax)
    temp %<>% filter(nseqs >= min.nb.seqs, div < mean(div)+10*sd(div)) %>%
      mutate('year' = as.numeric(year))
    temp %<>% group_by(pop) %>% mutate_at(vars(lat,long), median) %>%
      mutate('year' = ceiling(median(year))) %>%
      mutate_at(vars(div:ncomps,D:lu.div), mean) %>% ungroup %>%
      distinct(pop, .keep_all = T)
    temp <- temp %>% mutate_at(vars(D:lu.div), scale.fun) %>%
      mutate('lat.abs' = abs(lat)) %>%
      filter(lat.abs < 60) %>% mutate('lat.abs' = rescale(lat.abs,c(0,1))) %>%
      mutate_at(vars(pop,year,species, family, order), as.factor)
    
    #removing duplicate π estimates per species because max number is 1.9 π per species (on average)
    #and because at largest spatial scale, nb species == nb π values, so including a species
    #RE is equivalent to have an observation-level RE (i.e. overfitting)
    temp <- temp %>% arrange(species,desc(nseqs)) %>% distinct(species, .keep_all = T) %>%
      droplevels %>% mutate('wts' = log(nseqs)/mean(log(nseqs))) %>%
      select(-taxon,-scale,-nseqs, -ncomps,-n.years, -lu.var) %>% as.data.frame
    
    fullmod <- bam(div ~ zone +
                     s(lat,long, bs='gp') +
                     s(D, k = 8, bs = 'tp') +
                     s(lat.abs, k = 8, bs = 'tp') +
                     s(hd, k=8, bs ='tp') +
                     s(hd, k=8, by=zone) +
                     s(p.lu, k=8, bs='tp') +
                     s(p.lu, k=8, by=zone) +
                     s(year,bs = 're', k = 5, m=1) +
                     s(order, bs='re',k = 5, m=1) +
                     s(family, bs='re',k = 5, m=1),
                   data = temp, family = tw, method='fREML', discrete = T, weights = wts, nthreads = 2)
    
    modname <- paste('m1',tax,scl,sep='_')
    assign(modname,fullmod)
    
    rm(fullmod,temp)
    
  }
  
}

#### sketching these GAMMs ####

pdf('~/Desktop/SuppFig_zoneGAMMs_HD.pdf',width=7,height=7,pointsize = 8,onefile = T)
layout(rbind(c(1,5,9,13),c(2,6,10,14),c(3,7,11,15),c(4,8,12,16)))
par(cex=1,mar=c(2,2,1,1),oma=c(2.8,2.8,0,0))

for(i in 1:4){
  
  tax <- taxa[i]
  
  m1<-get(paste('m1',tax,'10',sep='_'))
  m2<-get(paste('m1',tax,'100',sep='_'))
  m3<-get(paste('m1',tax,'1000',sep='_'))
  m4<-get(paste('m1',tax,'10000',sep='_'))
  
  # gettting pvals for difference
  
  m1sum <- summary(m1)
  m2sum <- summary(m2)
  m3sum <- summary(m3)
  m4sum <- summary(m4)
  
  m1res <- m1sum$s.table[c('s(hd):zoneTropical'),c('p-value')]
  m1res <- round(m1res,4)
  m1res[m1res == 0] <- 0.0001
  
  m2res <- m2sum$s.table[c('s(hd):zoneTropical'),c('p-value')]
  m2res <- round(m2res,4)
  m2res[m2res == 0] <- 0.0001
  
  m3res <- m3sum$s.table[c('s(hd):zoneTropical'),c('p-value')]
  m3res <- round(m3res,4)
  m3res[m3res == 0] <- 0.0001
  
  m4res <- m4sum$s.table[c('s(hd):zoneTropical'),c('p-value')]
  m4res <- round(m4res,4)
  m4res[m4res == 0] <- 0.0001
  
  # hd
  
  plot_smooth(m1, view="hd", plot_all="zone", cond=list(D=0,long=0), col=c('black','royalblue2'), xlim=c(0,1), rm.ranef=T, se=1.96, yaxt='n',xaxt='n',ann=F, hide.label = T,main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA)
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=-10:0,labels=-10:0)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(0,0.25,0.5,0.75,1))
  legend('topleft',legend='10 km',bty='n')
  legend('bottomleft',bty='n',legend=bquote(italic('p') == .(m1res[1])))
  
  plot_smooth(m2, view="hd", plot_all="zone", cond=list(D=0,long=0), col=c('black','royalblue2'), xlim=c(0,1), rm.ranef=T, se=1.96, yaxt='n',xaxt='n',ann=F, hide.label = T,main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA)
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=-10:0,labels=-10:0)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(0,0.25,0.5,0.75,1))
  legend('topleft',legend='100 km',bty='n')
  legend('bottomleft',bty='n',legend=bquote(italic('p') == .(m2res[1])))
  
  plot_smooth(m3, view="hd", plot_all="zone", cond=list(D=0,long=0), col=c('black','royalblue2'), xlim=c(0,1), rm.ranef=T, se=1.96, yaxt='n',xaxt='n',ann=F, hide.label = T,main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA)
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=-10:0,labels=-10:0)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(0,0.25,0.5,0.75,1))
  legend('topleft',legend='1000 km',bty='n')
  legend('bottomleft',bty='n',legend=bquote(italic('p') == .(m3res[1])))
  
  plot_smooth(m4, view="hd", plot_all="zone", cond=list(D=0,long=0), col=c('black','royalblue2'), xlim=c(0,1), rm.ranef=T, se=1.96, yaxt='n',xaxt='n',ann=F, hide.label = T,main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA)
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=-10:0,labels=-10:0)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(0,0.25,0.5,0.75,1))
  legend('topleft',legend='10 000 km',bty='n')
  legend('bottomleft',bty='n',legend=bquote(italic('p') == .(m4res[1])))
  
}

mtext('human density (scaled)',side=1,outer=T,cex=1.2,at=0.5,line=1)
mtext(expression(log[e]~COI~diversity~(hat(pi))),at=.5,side=2,outer=T,cex=1.2,line=1)

dev.off()

##LU

pdf('~/Desktop/SuppFig_zoneGAMMs_LU.pdf',width=7,height=7,pointsize = 8,onefile = T)
layout(rbind(c(1,5,9,13),c(2,6,10,14),c(3,7,11,15),c(4,8,12,16)))
par(cex=1,mar=c(2,2,1,1),oma=c(2.8,2.8,0,0))

for(i in 1:4){
  
  tax <- taxa[i]
  
  m1<-get(paste('m1',tax,'10',sep='_'))
  m2<-get(paste('m1',tax,'100',sep='_'))
  m3<-get(paste('m1',tax,'1000',sep='_'))
  m4<-get(paste('m1',tax,'10000',sep='_'))
  
  # gettting pvals for difference
  
  m1sum <- summary(m1)
  m2sum <- summary(m2)
  m3sum <- summary(m3)
  m4sum <- summary(m4)
  
  m1res <- m1sum$s.table[c('s(p.lu):zoneTropical'),c('p-value')]
  m1res <- round(m1res,4)
  m1res[m1res == 0] <- 0.0001
  
  m2res <- m2sum$s.table[c('s(p.lu):zoneTropical'),c('p-value')]
  m2res <- round(m2res,4)
  m2res[m2res == 0] <- 0.0001
  
  m3res <- m3sum$s.table[c('s(p.lu):zoneTropical'),c('p-value')]
  m3res <- round(m3res,4)
  m3res[m3res == 0] <- 0.0001
  
  m4res <- m4sum$s.table[c('s(p.lu):zoneTropical'),c('p-value')]
  m4res <- round(m4res,4)
  m4res[m4res == 0] <- 0.0001
  
  # p.lu
  
  plot_smooth(m1, view="p.lu", plot_all="zone", cond=list(D=0,long=0), col=c('black','forestgreen'), xlim=c(0,1), rm.ranef=T, se=1.96, yaxt='n',xaxt='n',ann=F, hide.label = T,main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA)
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=-10:0,labels=-10:0)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(0,0.25,0.5,0.75,1))
  legend('topleft',legend='10 km',bty='n')
  legend('bottomleft',bty='n',legend=bquote(italic('p') == .(m1res[1])))
  
  plot_smooth(m2, view="p.lu", plot_all="zone", cond=list(D=0,long=0), col=c('black','forestgreen'), xlim=c(0,1), rm.ranef=T, se=1.96, yaxt='n',xaxt='n',ann=F, hide.label = T,main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA)
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=-10:0,labels=-10:0)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(0,0.25,0.5,0.75,1))
  legend('topleft',legend='100 km',bty='n')
  legend('bottomleft',bty='n',legend=bquote(italic('p') == .(m2res[1])))
  
  plot_smooth(m3, view="p.lu", plot_all="zone", cond=list(D=0,long=0), col=c('black','forestgreen'), xlim=c(0,1), rm.ranef=T, se=1.96, yaxt='n',xaxt='n',ann=F, hide.label = T,main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA)
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=-10:0,labels=-10:0)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(0,0.25,0.5,0.75,1))
  legend('topleft',legend='1000 km',bty='n')
  legend('bottomleft',bty='n',legend=bquote(italic('p') == .(m3res[1])))
  
  plot_smooth(m4, view="p.lu", plot_all="zone", cond=list(D=0,long=0), col=c('black','forestgreen'), xlim=c(0,1), rm.ranef=T, se=1.96, yaxt='n',xaxt='n',ann=F, hide.label = T,main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA)
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=-10:0,labels=-10:0)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(0,0.25,0.5,0.75,1))
  legend('topleft',legend='10 000 km',bty='n')
  legend('bottomleft',bty='n',legend=bquote(italic('p') == .(m4res[1])))
  
}

mtext('land use intensity (scaled)',side=1,outer=T,cex=1.2,at=0.5,line=1)
mtext(expression(log[e]~COI~diversity~(hat(pi))),at=.5,side=2,outer=T,cex=1.2,line=1)

dev.off()
