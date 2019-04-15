## Vincent Fugere 2019

# Code for Figures 2 & 3: maps and plots of spatial GAMMs

rm(list=ls())
options(tibble.print_max = 100, scipen = 999)

library(tidyverse)
library(magrittr)
library(RColorBrewer)
library(mgcv)
library(itsadug)
library(scales)
library(rworldmap)
library(gstat)
library(sp)
library(viridis)

load('~/Google Drive/Recherche/Intraspecific genetic diversity/Data/DF_Master.RData')
map <- getMap(resolution = "coarse")

#parameters
min.nb.seqs <- 2
#taxa <- c('birds','fish','insects','mammals')
taxa <- rep('birds',4)
scales <- c('10','100','1000','10000')
colz <- c(1,'#E69F00','#56B4E9','#009E73')

load('~/Desktop/spatialGAMMs.RData')
names(models) <- paste(c('m0','m1'),'birds',rep(scales, each = 2),sep='_')
list2env(models,envir=.GlobalEnv)

#### version of figure 3 with 10 seqs min ####

pdf('~/Desktop/Fig3_10seqs.pdf',width=8.5,height=9,pointsize = 8)
par(mfrow=c(4,4),cex=1)

tax <- taxa[1]

m1<-get(paste('m1',tax,'10',sep='_'))
m2<-get(paste('m1',tax,'100',sep='_'))
m3<-get(paste('m1',tax,'1000',sep='_'))
m4<-get(paste('m1',tax,'10000',sep='_'))

m1sum <- summary(m1)
m2sum <- summary(m2)
m3sum <- summary(m3)
m4sum <- summary(m4)

plot.tbl <- as.data.frame(rbind(round(m1sum$s.table[c(2,4,5,6),3:4],2),
                                round(m2sum$s.table[c(2,4,5,6),3:4],2),
                                round(m3sum$s.table[c(2,4,5,6),3:4],2),
                                round(m4sum$s.table[c(2,4,5,6),3:4],2)))

# cols <- viridis(20)[5:20]
# colfunc <- colorRampPalette(cols)
# colfunc(1000) -> cols.plot

colnames(plot.tbl)[1:2] <- c('Fval','pval')
plot.tbl$scale <- rep(scales, each=4)
plot.tbl$ln_wdth <-rescale(plot.tbl$Fval, to = c(1,4))
#plot.tbl$ln_wdth <- rep(c(0.6,1.2,1.8,2.4), each=4)
#plot.tbl$ln_col <- cols.plot[rescale(plot.tbl$Fval, to = c(1,1000))]
plot.tbl$ln_col <- rep(colz,each=4)
plot.tbl$ln_type <- 1
plot.tbl$ln_type[plot.tbl$pval > 0.05] <- 3
#plot.tbl$ln_col[plot.tbl$pval > 0.05] <- 'grey90'
plot.tbl$ctype <- rep(c('1-D','2-lat','3-hd','4-lu'),4)
plot.tbl <- arrange(plot.tbl, ctype, scale)

#ylims <-range(c(fitted(m1),fitted(m2),fitted(m3),fitted(m4)))

plot_smooth(m1, view="D", cond=list('lat' = 0, 'long' = 0), col=alpha(plot.tbl$ln_col[1],0.8), rm.ranef=T, se=0, lty=plot.tbl$ln_type[1], lwd=plot.tbl$ln_wdth[1], hide.label = T,xlab='mean spatial distance (D)',ylab=expression(COI~hat(pi)),main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA, transform = exp,ylim=c(0.0005,0.005))
plot_smooth(m2, view="D", cond=list('lat' = 0, 'long' = 0), col=alpha(plot.tbl$ln_col[2],0.8), rm.ranef=T, se=0, lty=plot.tbl$ln_type[2], lwd=plot.tbl$ln_wdth[2], add=T, rug=F, transform = exp)
plot_smooth(m3, view="D", cond=list('lat' = 0, 'long' = 0), col=alpha(plot.tbl$ln_col[3],0.8), rm.ranef=T, se=0, lty=plot.tbl$ln_type[3], lwd=plot.tbl$ln_wdth[3], add=T, rug=F, transform = exp)
plot_smooth(m4, view="D", cond=list('lat' = 0, 'long' = 0), col=alpha(plot.tbl$ln_col[4],0.8), rm.ranef=T, se=0, lty=plot.tbl$ln_type[4], lwd=plot.tbl$ln_wdth[4], add=T, rug=F, transform = exp)
legend('topleft',bty='n',legend=c('scale = 10 km','100 km','1 000 km','10 000 km'),pch=16,col=colz)

plot_smooth(m1, view="lat.abs", cond=list('lat' = 0, 'long' = 0), col=alpha(plot.tbl$ln_col[5],0.8), rm.ranef=T, se=0, lty=plot.tbl$ln_type[5], lwd=plot.tbl$ln_wdth[5], hide.label = T,xlab=expression(latitude~(absolute)),ylab=expression(COI~hat(pi)),main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA, transform = exp,ylim=c(0.0005,0.005))
plot_smooth(m2, view="lat.abs", cond=list('lat' = 0, 'long' = 0), col=alpha(plot.tbl$ln_col[6],0.8), rm.ranef=T, se=0, lty=plot.tbl$ln_type[6], lwd=plot.tbl$ln_wdth[6], add=T, rug=F, transform = exp)
plot_smooth(m3, view="lat.abs", cond=list('lat' = 0, 'long' = 0), col=alpha(plot.tbl$ln_col[7],0.8), rm.ranef=T, se=0, lty=plot.tbl$ln_type[7], lwd=plot.tbl$ln_wdth[7], add=T, rug=F, transform = exp)
plot_smooth(m4, view="lat.abs", cond=list('lat' = 0, 'long' = 0), col=alpha(plot.tbl$ln_col[8],0.8), rm.ranef=T, se=0, lty=plot.tbl$ln_type[8], lwd=plot.tbl$ln_wdth[8], add=T, rug=F, transform = exp)

plot_smooth(m1, view="hd", cond=list('lat' = 0, 'long' = 0), col=alpha(plot.tbl$ln_col[9],0.8), rm.ranef=T, se=0, lty=plot.tbl$ln_type[9], lwd=plot.tbl$ln_wdth[9], hide.label = T,xlab=expression(human~density~(people~km^-2)),ylab=expression(COI~hat(pi)),main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA, transform = exp,ylim=c(0.0005,0.005))
plot_smooth(m2, view="hd", cond=list('lat' = 0, 'long' = 0), col=alpha(plot.tbl$ln_col[10],0.8), rm.ranef=T, se=0, lty=plot.tbl$ln_type[10], lwd=plot.tbl$ln_wdth[10], add=T, rug=F, transform = exp)
plot_smooth(m3, view="hd", cond=list('lat' = 0, 'long' = 0), col=alpha(plot.tbl$ln_col[11],0.8), rm.ranef=T, se=0, lty=plot.tbl$ln_type[11], lwd=plot.tbl$ln_wdth[11], add=T, rug=F, transform = exp)
plot_smooth(m4, view="hd", cond=list('lat' = 0, 'long' = 0), col=alpha(plot.tbl$ln_col[12],0.8), rm.ranef=T, se=0, lty=plot.tbl$ln_type[12], lwd=plot.tbl$ln_wdth[12], add=T, rug=F, transform = exp)

plot_smooth(m1, view="p.lu", cond=list('lat' = 0, 'long' = 0), col=alpha(plot.tbl$ln_col[13],0.8), rm.ranef=T, se=0, lty=plot.tbl$ln_type[13], lwd=plot.tbl$ln_wdth[13], hide.label = T,xlab=expression(land~use~intensity),ylab=expression(COI~hat(pi)),main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA, transform = exp,ylim=c(0.0005,0.005))
plot_smooth(m2, view="p.lu", cond=list('lat' = 0, 'long' = 0), col=alpha(plot.tbl$ln_col[14],0.8), rm.ranef=T, se=0, lty=plot.tbl$ln_type[14], lwd=plot.tbl$ln_wdth[14], add=T, rug=F, transform = exp)
plot_smooth(m3, view="p.lu", cond=list('lat' = 0, 'long' = 0), col=alpha(plot.tbl$ln_col[15],0.8), rm.ranef=T, se=0, lty=plot.tbl$ln_type[15], lwd=plot.tbl$ln_wdth[15], add=T, rug=F, transform = exp)
plot_smooth(m4, view="p.lu", cond=list('lat' = 0, 'long' = 0), col=alpha(plot.tbl$ln_col[16],0.8), rm.ranef=T, se=0, lty=plot.tbl$ln_type[16], lwd=plot.tbl$ln_wdth[16], add=T, rug=F, transform = exp)

dev.off()

