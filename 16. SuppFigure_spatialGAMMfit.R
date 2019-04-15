## Vincent Fugere 2019

## Supp Figure showing GAMM fit for all taxa

rm(list=ls())
options(tibble.print_max = 100, scipen = 999)

library(tidyverse)
library(magrittr)
library(RColorBrewer)
library(mgcv)
library(itsadug)
library(scales)
library(viridis)

taxa <- c('birds','fish','insects','mammals')
scales <- c('10','100','1000','10000')
colz <- c(1,'#E69F00','#56B4E9','#009E73')

load('~/Google Drive/Recherche/Intraspecific genetic diversity/Data/spatialGAMMs.RData')
list2env(models,envir=.GlobalEnv)

##

pdf('~/Desktop/SuppFig_spatialGAMMs.pdf',width=7,height=7,pointsize = 8,onefile = T)
layout(rbind(c(1,5,9,13),c(2,6,10,14),c(3,7,11,15),c(4,8,12,16)))
par(cex=1,mar=c(2,2,1,1),oma=c(2.5,2.8,0,0))

ymins <- c(-9,-9,-10,-10)
ymaxs <-c(-4,-2,-2,-2)

for(i in 1:4){
  
  tax <- taxa[i]
  ylims <- c(ymins[i],ymaxs[i])

  m1<-get(paste('m1',tax,'10',sep='_'))
  m2<-get(paste('m1',tax,'100',sep='_'))
  m3<-get(paste('m1',tax,'1000',sep='_'))
  m4<-get(paste('m1',tax,'10000',sep='_'))
  
  # D
  
  plot_smooth(m1, view="D", col=colz[1], xlim=c(0,1),cond=list('lat' = 0, 'long' = 0), rm.ranef=F, se=1.96, yaxt='n',xaxt='n',ann=F, hide.label = T,main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA, ylim = ylims)
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=-10:0,labels=-10:0)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(0,0.25,0.5,0.75,1),labels = c('0','0.25','0.5','0.75','1'))
  legend('topleft',bty='n',legend = 'scale = 10 km', col = colz[1])
  
  plot_smooth(m2, view="D", col=colz[2], xlim=c(0,1),cond=list('lat' = 0, 'long' = 0), rm.ranef=F, se=1.96, yaxt='n',xaxt='n',ann=F, hide.label = T,main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA, ylim = ylims)
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=-10:0,labels=-10:0)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(0,0.25,0.5,0.75,1),labels = c('0','0.25','0.5','0.75','1'))
  legend('topleft',bty='n',legend = 'scale = 100 km', col = colz[2])
  
  plot_smooth(m3, view="D", col=colz[3], xlim=c(0,1),cond=list('lat' = 0, 'long' = 0), rm.ranef=F, se=1.96, yaxt='n',xaxt='n',ann=F, hide.label = T,main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA, ylim = ylims)
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=-10:0,labels=-10:0)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(0,0.25,0.5,0.75,1),labels = c('0','0.25','0.5','0.75','1'))
  legend('topleft',bty='n',legend = 'scale = 1 000 km', col = colz[3])
  
  plot_smooth(m4, view="D", col=colz[4], xlim=c(0,1),cond=list('lat' = 0, 'long' = 0), rm.ranef=F, se=1.96, yaxt='n',xaxt='n',ann=F, hide.label = T,main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA, ylim = ylims)
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=-10:0,labels=-10:0)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(0,0.25,0.5,0.75,1),labels = c('0','0.25','0.5','0.75','1'))
  legend('topleft',bty='n',legend = 'scale = 10 000 km', col = colz[4])
  
  mtext('mean spatial distance',side=1,outer=F,cex=1.2,adj=0.5,line=3)
  
  # lat

  plot_smooth(m1, view="lat.abs", col=colz[1], xlim=c(0,1),cond=list('lat' = 0, 'long' = 0), rm.ranef=F, se=1.96, yaxt='n',xaxt='n',ann=F, hide.label = T,main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA, ylim = ylims)
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=-10:0,labels=-10:0)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(0,0.25,0.5,0.75,1),labels = c('0','0.25','0.5','0.75','1'))
  
  plot_smooth(m2, view="lat.abs", col=colz[2], xlim=c(0,1),cond=list('lat' = 0, 'long' = 0), rm.ranef=F, se=1.96, yaxt='n',xaxt='n',ann=F, hide.label = T,main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA, ylim = ylims)
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=-10:0,labels=-10:0)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(0,0.25,0.5,0.75,1),labels = c('0','0.25','0.5','0.75','1'))
  
  plot_smooth(m3, view="lat.abs", col=colz[3], xlim=c(0,1),cond=list('lat' = 0, 'long' = 0), rm.ranef=F, se=1.96, yaxt='n',xaxt='n',ann=F, hide.label = T,main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA, ylim = ylims)
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=-10:0,labels=-10:0)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(0,0.25,0.5,0.75,1),labels = c('0','0.25','0.5','0.75','1'))
  
  plot_smooth(m4, view="lat.abs", col=colz[4], xlim=c(0,1),cond=list('lat' = 0, 'long' = 0), rm.ranef=F, se=1.96, yaxt='n',xaxt='n',ann=F, hide.label = T,main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA, ylim = ylims)
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=-10:0,labels=-10:0)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(0,0.25,0.5,0.75,1),labels = c('0','0.25','0.5','0.75','1'))
  
  mtext('absolute latitude',side=1,outer=F,cex=1.2,adj=0.5,line=3)
  
  # HD
  
  plot_smooth(m1, view="hd", col=colz[1], xlim=c(0,1),cond=list('lat' = 0, 'long' = 0), rm.ranef=F, se=1.96, yaxt='n',xaxt='n',ann=F, hide.label = T,main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA, ylim = ylims)
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=-10:0,labels=-10:0)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(0,0.25,0.5,0.75,1),labels = c('0','0.25','0.5','0.75','1'))
  
  plot_smooth(m2, view="hd", col=colz[2], xlim=c(0,1),cond=list('lat' = 0, 'long' = 0), rm.ranef=F, se=1.96, yaxt='n',xaxt='n',ann=F, hide.label = T,main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA, ylim = ylims)
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=-10:0,labels=-10:0)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(0,0.25,0.5,0.75,1),labels = c('0','0.25','0.5','0.75','1'))
  
  plot_smooth(m3, view="hd", col=colz[3], xlim=c(0,1),cond=list('lat' = 0, 'long' = 0), rm.ranef=F, se=1.96, yaxt='n',xaxt='n',ann=F, hide.label = T,main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA, ylim = ylims)
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=-10:0,labels=-10:0)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(0,0.25,0.5,0.75,1),labels = c('0','0.25','0.5','0.75','1'))
  
  plot_smooth(m4, view="hd", col=colz[4], xlim=c(0,1),cond=list('lat' = 0, 'long' = 0), rm.ranef=F, se=1.96, yaxt='n',xaxt='n',ann=F, hide.label = T,main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA, ylim = ylims)
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=-10:0,labels=-10:0)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(0,0.25,0.5,0.75,1),labels = c('0','0.25','0.5','0.75','1'))
  
  mtext('human density',side=1,outer=F,cex=1.2,adj=0.5,line=3)
  
  # LU
  
  plot_smooth(m1, view="p.lu", col=colz[1], xlim=c(0,1),cond=list('lat' = 0, 'long' = 0), rm.ranef=F, se=1.96, yaxt='n',xaxt='n',ann=F, hide.label = T,main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA, ylim = ylims)
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=-10:0,labels=-10:0)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(0,0.25,0.5,0.75,1),labels = c('0','0.25','0.5','0.75','1'))
  
  plot_smooth(m2, view="p.lu", col=colz[2], xlim=c(0,1),cond=list('lat' = 0, 'long' = 0), rm.ranef=F, se=1.96, yaxt='n',xaxt='n',ann=F, hide.label = T,main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA, ylim = ylims)
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=-10:0,labels=-10:0)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(0,0.25,0.5,0.75,1),labels = c('0','0.25','0.5','0.75','1'))
  
  plot_smooth(m3, view="p.lu", col=colz[3], xlim=c(0,1),cond=list('lat' = 0, 'long' = 0), rm.ranef=F, se=1.96, yaxt='n',xaxt='n',ann=F, hide.label = T,main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA, ylim = ylims)
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=-10:0,labels=-10:0)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(0,0.25,0.5,0.75,1),labels = c('0','0.25','0.5','0.75','1'))
  
  plot_smooth(m4, view="p.lu", col=colz[4], xlim=c(0,1),cond=list('lat' = 0, 'long' = 0), rm.ranef=F, se=1.96, yaxt='n',xaxt='n',ann=F, hide.label = T,main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA, ylim = ylims)
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=-10:0,labels=-10:0)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(0,0.25,0.5,0.75,1),labels = c('0','0.25','0.5','0.75','1'))
  
  mtext('land use intensity',side=1,outer=F,cex=1.2,adj=0.5,line=3)

  mtext(expression(log[e]~COI~diversity~(hat(pi))),at=.5,side=2,outer=T,cex=1.2,line=1)
  
}

dev.off()

