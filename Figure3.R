## Vincent Fugere 2019

# Code for Figure 3: spatial GAMM testing land use impacts on gen div

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
library(tictoc)

plot.s <- function(model,xvar,cond,col,xlab,ylab){
  plot_smooth(model, view=xvar, cond=cond, col=col, rm.ranef=T, hide.label = T,xlab=xlab,ylab=ylab,main=NULL,rug=F,bty='l',legend_plot_all = F, h0=NA)
}

map <- getMap(resolution = "coarse")

load('~/Google Drive/Recherche/Intraspecific genetic diversity/Data/DF_Master.RData')

#start with a single scale
scale.fun <-function(x){y <- scales::rescale(x, to = c(0,1)); return(y)}
temp <- DF %>% filter(scale == '1000', taxon == 'birds')
temp <- temp %>% mutate_at(vars(D:lu.div), scale.fun)
temp %<>% mutate('lat.sq' = scale.fun(lat^2), 'wts' = log(nseqs)/mean(log(nseqs)))
temp %<>% select(-taxon,-scale,-nseqs, -ncomps,-n.years, -lu.var) %>%
  mutate_at(vars(pop,year,species), as.factor) %>% as.data.frame
temp %<>% filter(div < mean(div)+10*sd(div))

#smoothed genetic diversity, accounting for D and Yr

m1 <- bam(div ~ s(lat,long, k = 500, m = 1) + s(D, k = 20) + s(year,bs = 're', k = 40), data = temp, family = tw, method='fREML')
gam.check(m1)

plot.s(m1,'D',list('lat' = 0, 'long' = 0),1,'D','div')

cols <- rev(c('#d73027','#f46d43','#fdae61','#fee090','#e0f3f8','#abd9e9','#74add1','#4575b4'))
colfunc <- colorRampPalette(cols)
colfunc(1000) -> cols.plot

max.div <- quantile(fitted(m1),0.95)
temp$fit <- fitted(m1)
temp$fit[temp$fit > max.div] <- max.div
temp$fit.sc <- rescale(temp$fit,to=c(0,1000))

#temp <- temp %>% arrange(fit)
plot(map, xlim = c(-180,180), ylim = c(-90,90), border=NA,col='grey95',axes=F,asp=1,cex.lab=0.5)
polygon(x=c(-180,180,180,-180),y=c(-60,-60,-90,-90),col='white',border=NA)
points(lat~long,temp,pch=21,bg='white',col=1,cex=0.5)
points(lat~long,temp,pch=16,col=alpha(cols.plot[temp$fit.sc],1),cex=0.5)
legs <- as.character(seq(0,max.div*7/8,length.out = 8))
legs[8] <- paste('>',max.div*7/8)
xseqs <- seq(-15,60,length.out = 8)
rect(xleft=xseqs,xright=xseqs+(xseqs[2]-xseqs[1]),ybottom = rep(-55,8),ytop = rep(-50,8),col=cols,border=NULL,lwd=0.2)
segments(x0=xseqs[2:8],x1=xseqs[2:8],y0=rep(-50,7),y1=c(-59,-57,-57,-59,-57,-57,-59),lwd=c(0.5,0.3,0.3,0.5,0.3,0.3,0.5))
text(x=xseqs[c(2,5,8)],y=rep(-59,3),labels = legs[c(2,5,8)],pos=1)
text(x=27.85714,y=-49,cex=1,label=smoothed~COI~pi,pos=3)

tic()
m1 <- bam(div ~ s(lat,long, bs='gp', k = 50, m = 1) + s(D, bs='cr', k = 20) + s(year,bs = 're', k = 40) + s(species,bs = 're', k = 2), data = temp, family = tw, method='fREML')
toc()

m1 <- bam(div ~ s(lat,long) + s(lat.sq) + s(D) + s(hd) + s(hd.var) + s(p.lu) + s(lu.div) + s(year, bs='re') + s(species, bs = 're'),
          data = temp, family = tw, method='fREML',cluster=cl)

summary(m1)
gam.check(m1)


#spatial autocorrelation in residuals?
temp$R <- resid(m1)
temp$rcol <- 4
temp[temp$R < 0, 'rcol'] <- 2
plot(map, xlim = c(-180,180), ylim = c(-90,90),border=1,col=NA,axes=F)
points(lat~long,temp,pch=16,col=alpha(temp$rcol,0.5),cex=0.5+abs(temp$R))

spatdat <- select(temp, long,lat,R)
coordinates(spatdat) <- c('long','lat')
variogram(R~long+lat,spatdat,width=0.1,cutoff=50) -> var1
scatter.smooth(x=var1$dist,y=var1$gamma,xlab='distance',ylab='semivariance',ylim=c(0,max(var1$gamma)),pch=16,col='gray')


