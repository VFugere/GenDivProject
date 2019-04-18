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

load('~/Google Drive/Recherche/Intraspecific genetic diversity/Data/spatialGAMMs_5seqs.RData')
list2env(models,envir=.GlobalEnv)

##

pdf('~/Desktop/SuppFig_spatialGAMMs.pdf',width=7,height=7,pointsize = 8,onefile = T)
layout(rbind(c(1,5,9,13),c(2,6,10,14),c(3,7,11,15),c(4,8,12,16)))
par(cex=1,mar=c(2,2,1,1),oma=c(2.5,2.8,0,0))

# #2seqs
# ymins <- c(-9,-9,-10,-10)
# ymaxs <-c(-4,-2,-2,-2)

#5seqs
ymins <- c(-8,-10,-10,-10)
ymaxs <- c(-2,-2,-2,-2)

for(i in 1:4){
  
  tax <- taxa[i]
  ylims <- c(ymins[i],ymaxs[i])

  m1<-get(paste('m1',tax,'10',sep='_'))
  m2<-get(paste('m1',tax,'100',sep='_'))
  m3<-get(paste('m1',tax,'1000',sep='_'))
  m4<-get(paste('m1',tax,'10000',sep='_'))
  
  m1sum <- summary(m1)
  m2sum <- summary(m2)
  m3sum <- summary(m3)
  m4sum <- summary(m4)
  
  m1res <- m1sum$s.table[c('s(D)','s(lat.abs)','s(hd)','s(p.lu)'),c('F','p-value')]
  m1res[,1] <- round(m1res[,1],2)
  m1res[,2] <- round(m1res[,2],4)
  m1res[m1res[,2] == 0,2] <- 0.0001
  m1rsq <- round(m1sum$r.sq,2)
  
  m2res <- m2sum$s.table[c('s(D)','s(lat.abs)','s(hd)','s(p.lu)'),c('F','p-value')]
  m2res[,1] <- round(m2res[,1],2)
  m2res[,2] <- round(m2res[,2],4)
  m2res[m2res[,2] == 0,2] <- 0.0001
  m2rsq <- round(m2sum$r.sq,2)
  
  m3res <- m3sum$s.table[c('s(D)','s(lat.abs)','s(hd)','s(p.lu)'),c('F','p-value')]
  m3res[,1] <- round(m3res[,1],2)
  m3res[,2] <- round(m3res[,2],4)
  m3res[m3res[,2] == 0,2] <- 0.0001
  m3rsq <- round(m3sum$r.sq,2)
  
  m4res <- m4sum$s.table[c('s(D)','s(lat.abs)','s(hd)','s(p.lu)'),c('F','p-value')]
  m4res[,1] <- round(m4res[,1],2)
  m4res[,2] <- round(m4res[,2],4)
  m4res[m4res[,2] == 0,2] <- 0.0001
  m4rsq <- round(m4sum$r.sq,2)
  
  
  # D
  
  plot_smooth(m1, view="D", col=colz[1], xlim=c(0,1),cond=list('lat' = 0, 'long' = 0), rm.ranef=F, se=1.96, yaxt='n',xaxt='n',ann=F, hide.label = T,main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA, ylim = ylims)
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=-10:0,labels=-10:0)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(0,0.25,0.5,0.75,1),labels = c('0','0.25','0.5','0.75','1'))
  legend('bottomleft',bty='n',legend = c(as.expression('scale = 10 km'),as.expression(bquote('model'~italic(R)^2 == .(m1rsq)))))
  legend('topright',bty='n',legend=bquote(italic('F') == .(m1res[1,1])~','~italic('p') == .(m1res[1,2])))
  
  plot_smooth(m2, view="D", col=colz[2], xlim=c(0,1),cond=list('lat' = 0, 'long' = 0), rm.ranef=F, se=1.96, yaxt='n',xaxt='n',ann=F, hide.label = T,main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA, ylim = ylims)
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=-10:0,labels=-10:0)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(0,0.25,0.5,0.75,1),labels = c('0','0.25','0.5','0.75','1'))
  legend('bottomleft',bty='n',legend = c(as.expression('scale = 100 km'),as.expression(bquote('model'~italic(R)^2 == .(m2rsq)))))
  legend('topright',bty='n',legend=bquote(italic('F') == .(m2res[1,1])~','~italic('p') == .(m2res[1,2])))
  
  plot_smooth(m3, view="D", col=colz[3], xlim=c(0,1),cond=list('lat' = 0, 'long' = 0), rm.ranef=F, se=1.96, yaxt='n',xaxt='n',ann=F, hide.label = T,main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA, ylim = ylims)
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=-10:0,labels=-10:0)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(0,0.25,0.5,0.75,1),labels = c('0','0.25','0.5','0.75','1'))
  legend('bottomleft',bty='n',legend = c(as.expression('scale = 1 000 km'),as.expression(bquote('model'~italic(R)^2 == .(m3rsq)))))
  legend('topright',bty='n',legend=bquote(italic('F') == .(m3res[1,1])~','~italic('p') == .(m3res[1,2])))
  
  plot_smooth(m4, view="D", col=colz[4], xlim=c(0,1),cond=list('lat' = 0, 'long' = 0), rm.ranef=F, se=1.96, yaxt='n',xaxt='n',ann=F, hide.label = T,main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA, ylim = ylims)
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=-10:0,labels=-10:0)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(0,0.25,0.5,0.75,1),labels = c('0','0.25','0.5','0.75','1'))
  legend('bottomleft',bty='n',legend = c(as.expression('scale = 10 000 km'),as.expression(bquote('model'~italic(R)^2 == .(m4rsq)))))
  legend('topright',bty='n',legend=bquote(italic('F') == .(m4res[1,1])~','~italic('p') == .(m4res[1,2])))
  
  mtext('mean spatial distance',side=1,outer=F,cex=1.2,adj=0.5,line=3)
  
  # lat

  plot_smooth(m1, view="lat.abs", col=colz[1], xlim=c(0,1),cond=list('lat' = 0, 'long' = 0), rm.ranef=F, se=1.96, yaxt='n',xaxt='n',ann=F, hide.label = T,main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA, ylim = ylims)
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=-10:0,labels=-10:0)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(0,0.25,0.5,0.75,1),labels = c('0','0.25','0.5','0.75','1'))
  legend('topright',bty='n',legend=bquote(italic('F') == .(m1res[2,1])~','~italic('p') == .(m1res[2,2])))
  
  plot_smooth(m2, view="lat.abs", col=colz[2], xlim=c(0,1),cond=list('lat' = 0, 'long' = 0), rm.ranef=F, se=1.96, yaxt='n',xaxt='n',ann=F, hide.label = T,main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA, ylim = ylims)
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=-10:0,labels=-10:0)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(0,0.25,0.5,0.75,1),labels = c('0','0.25','0.5','0.75','1'))
  legend('topright',bty='n',legend=bquote(italic('F') == .(m2res[2,1])~','~italic('p') == .(m2res[2,2])))
  
  plot_smooth(m3, view="lat.abs", col=colz[3], xlim=c(0,1),cond=list('lat' = 0, 'long' = 0), rm.ranef=F, se=1.96, yaxt='n',xaxt='n',ann=F, hide.label = T,main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA, ylim = ylims)
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=-10:0,labels=-10:0)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(0,0.25,0.5,0.75,1),labels = c('0','0.25','0.5','0.75','1'))
  legend('topright',bty='n',legend=bquote(italic('F') == .(m3res[2,1])~','~italic('p') == .(m3res[2,2])))
  
  plot_smooth(m4, view="lat.abs", col=colz[4], xlim=c(0,1),cond=list('lat' = 0, 'long' = 0), rm.ranef=F, se=1.96, yaxt='n',xaxt='n',ann=F, hide.label = T,main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA, ylim = ylims)
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=-10:0,labels=-10:0)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(0,0.25,0.5,0.75,1),labels = c('0','0.25','0.5','0.75','1'))
  legend('topright',bty='n',legend=bquote(italic('F') == .(m4res[2,1])~','~italic('p') == .(m4res[2,2])))
  
  mtext('absolute latitude',side=1,outer=F,cex=1.2,adj=0.5,line=3)
  
  # HD
  
  plot_smooth(m1, view="hd", col=colz[1], xlim=c(0,1),cond=list('lat' = 0, 'long' = 0), rm.ranef=F, se=1.96, yaxt='n',xaxt='n',ann=F, hide.label = T,main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA, ylim = ylims)
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=-10:0,labels=-10:0)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(0,0.25,0.5,0.75,1),labels = c('0','0.25','0.5','0.75','1'))
  legend('topright',bty='n',legend=bquote(italic('F') == .(m1res[3,1])~','~italic('p') == .(m1res[3,2])))
  
  plot_smooth(m2, view="hd", col=colz[2], xlim=c(0,1),cond=list('lat' = 0, 'long' = 0), rm.ranef=F, se=1.96, yaxt='n',xaxt='n',ann=F, hide.label = T,main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA, ylim = ylims)
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=-10:0,labels=-10:0)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(0,0.25,0.5,0.75,1),labels = c('0','0.25','0.5','0.75','1'))
  legend('topright',bty='n',legend=bquote(italic('F') == .(m2res[3,1])~','~italic('p') == .(m2res[3,2])))
  
  plot_smooth(m3, view="hd", col=colz[3], xlim=c(0,1),cond=list('lat' = 0, 'long' = 0), rm.ranef=F, se=1.96, yaxt='n',xaxt='n',ann=F, hide.label = T,main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA, ylim = ylims)
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=-10:0,labels=-10:0)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(0,0.25,0.5,0.75,1),labels = c('0','0.25','0.5','0.75','1'))
  legend('topright',bty='n',legend=bquote(italic('F') == .(m3res[3,1])~','~italic('p') == .(m3res[3,2])))
  
  plot_smooth(m4, view="hd", col=colz[4], xlim=c(0,1),cond=list('lat' = 0, 'long' = 0), rm.ranef=F, se=1.96, yaxt='n',xaxt='n',ann=F, hide.label = T,main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA, ylim = ylims)
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=-10:0,labels=-10:0)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(0,0.25,0.5,0.75,1),labels = c('0','0.25','0.5','0.75','1'))
  legend('topright',bty='n',legend=bquote(italic('F') == .(m4res[3,1])~','~italic('p') == .(m4res[3,2])))
  
  mtext('human density',side=1,outer=F,cex=1.2,adj=0.5,line=3)
  
  # LU
  
  plot_smooth(m1, view="p.lu", col=colz[1], xlim=c(0,1),cond=list('lat' = 0, 'long' = 0), rm.ranef=F, se=1.96, yaxt='n',xaxt='n',ann=F, hide.label = T,main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA, ylim = ylims)
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=-10:0,labels=-10:0)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(0,0.25,0.5,0.75,1),labels = c('0','0.25','0.5','0.75','1'))
  legend('topright',bty='n',legend=bquote(italic('F') == .(m1res[4,1])~','~italic('p') == .(m1res[4,2])))
  
  plot_smooth(m2, view="p.lu", col=colz[2], xlim=c(0,1),cond=list('lat' = 0, 'long' = 0), rm.ranef=F, se=1.96, yaxt='n',xaxt='n',ann=F, hide.label = T,main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA, ylim = ylims)
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=-10:0,labels=-10:0)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(0,0.25,0.5,0.75,1),labels = c('0','0.25','0.5','0.75','1'))
  legend('topright',bty='n',legend=bquote(italic('F') == .(m2res[4,1])~','~italic('p') == .(m2res[4,2])))
  
  plot_smooth(m3, view="p.lu", col=colz[3], xlim=c(0,1),cond=list('lat' = 0, 'long' = 0), rm.ranef=F, se=1.96, yaxt='n',xaxt='n',ann=F, hide.label = T,main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA, ylim = ylims)
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=-10:0,labels=-10:0)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(0,0.25,0.5,0.75,1),labels = c('0','0.25','0.5','0.75','1'))
  legend('topright',bty='n',legend=bquote(italic('F') == .(m3res[4,1])~','~italic('p') == .(m3res[4,2])))
  
  plot_smooth(m4, view="p.lu", col=colz[4], xlim=c(0,1),cond=list('lat' = 0, 'long' = 0), rm.ranef=F, se=1.96, yaxt='n',xaxt='n',ann=F, hide.label = T,main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA, ylim = ylims)
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=-10:0,labels=-10:0)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(0,0.25,0.5,0.75,1),labels = c('0','0.25','0.5','0.75','1'))
  legend('topright',bty='n',legend=bquote(italic('F') == .(m4res[4,1])~','~italic('p') == .(m4res[4,2])))
  
  mtext('land use intensity',side=1,outer=F,cex=1.2,adj=0.5,line=3)

  mtext(expression(log[e]~COI~diversity~(hat(pi))),at=.5,side=2,outer=T,cex=1.2,line=1)
  
}

dev.off()

