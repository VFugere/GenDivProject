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
taxa <- c('birds','fish','insects','mammals')
scales <- c('10','100','1000','10000')
#colz <- brewer.pal(4, 'Dark2')

colz <- c(1,'#E69F00','#56B4E9','#009E73')

load('~/Desktop/spatialGAMMs.RData')
names(models) <- paste(c('m0','m1'),'birds',rep(scales, each = 2),sep='_')
list2env(models,envir=.GlobalEnv)

tax <- taxa[1]
scl <- scales[3]

temp <- DF %>% filter(scale == scl, taxon == tax)
temp %<>% filter(nseqs >= min.nb.seqs, div < mean(div)+10*sd(div)) %>%
  mutate('year' = as.numeric(year))
temp %<>% group_by(pop) %>% mutate_at(vars(lat,long), median) %>%
  mutate('year' = ceiling(median(year))) %>%
  mutate_at(vars(div:ncomps,D:lu.div), mean) %>% ungroup %>%
  distinct(pop, .keep_all = T)

m1<-m1_birds_10
m2<-m1_birds_100
m3<-m1_birds_1000
m4<-m1_birds_10000

pdf('~/Desktop/Fig3.pdf',width=8,height=10,pointsize = 8)
par(mfrow=c(4,4),cex=1)

m1sum <- summary(m1)
m2sum <- summary(m2)
m3sum <- summary(m3)
m4sum <- summary(m4)

plot.tbl <- as.data.frame(rbind(round(m1sum$s.table[c(2,4,5,6),3:4],2),
                                round(m2sum$s.table[c(2,4,5,6),3:4],2),
                                round(m3sum$s.table[c(2,4,5,6),3:4],2),
                                round(m4sum$s.table[c(2,4,5,6),3:4],2)))

cols <- viridis(20)[5:20]
colfunc <- colorRampPalette(cols)
colfunc(1000) -> cols.plot

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

pdf('~/Desktop/Fig3.pdf',width=8,height=10,pointsize = 8)
par(mfrow=c(4,4),cex=1)

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



cols <- rev(c('#d73027','#f46d43','#fdae61','#fee090','#e0f3f8','#abd9e9','#74add1','#4575b4'))
colfunc <- colorRampPalette(cols)
colfunc(1000) -> cols.plot

max.div <- quantile(fitted(pmod),0.95)
temp$fit <- fitted(pmod)
temp$fit[temp$fit > max.div] <- max.div
temp$fit.sc <- rescale(temp$fit,to=c(0,1000))

temp <- temp %>% mutate('long' = (ceiling(long)-ceiling(long)%%4), 'lat' = (ceiling(lat)-ceiling(lat)%%4)) %>% 
  group_by(long,lat) %>% summarize('fit.sc' = mean(fit.sc,na.rm=T))
pdf('~/Desktop/birds_1000.pdf',width = 8, height = 4, pointsize = 8)
plot(map, xlim = c(-180,180), ylim = c(-90,90), border=NA,col='grey95',axes=F,asp=1,cex.lab=0.5)
polygon(x=c(-180,180,180,-180),y=c(-60,-60,-90,-90),col='white',border=NA)
points(lat~long,temp,pch=21,bg='white',col=1,cex=1.4)
points(lat~long,temp,pch=16,col=alpha(cols.plot[temp$fit.sc],1),cex=1.2)
legs <- as.character(round(seq(0,max.div*7/8,length.out = 8),3))
legs[8] <- paste('>',round(max.div*7/8,3))
xseqs <- seq(-15,60,length.out = 8)
rect(xleft=xseqs,xright=xseqs+(xseqs[2]-xseqs[1]),ybottom = rep(-55,8),ytop = rep(-50,8),col=cols,border=NULL,lwd=0.2)
segments(x0=xseqs[2:8],x1=xseqs[2:8],y0=rep(-50,7),y1=c(-59,-57,-57,-59,-57,-57,-59),lwd=c(0.5,0.3,0.3,0.5,0.3,0.3,0.5))
text(x=xseqs[c(2,5,8)],y=rep(-59,3),labels = legs[c(2,5,8)],pos=1)
text(x=27.85714,y=-49,cex=1,label=smoothed~COI~hat(pi),pos=3)
dev.off()
