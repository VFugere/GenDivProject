## Vincent Fugere 2019

# Code for Figure 3: spatial GAMM testing land use impacts on gen div

rm(list=ls())
options(tibble.print_max = 100, scipen = 999)

library(tidyverse)
library(magrittr)
library(RColorBrewer)
library(mgcv)
library(scales)
library(rworldmap)
library(gstat)
library(sp)

map <- getMap(resolution = "coarse")

load('~/Google Drive/Recherche/Intraspecific genetic diversity/Data/DF_Master.RData')

#start with a single scale
scale.fun <-function(x){y <- scales::rescale(x, to = c(0,1)); return(y)}
temp <- DF %>% filter(scale == '10', taxon == 'birds')
temp <- temp %>% mutate_at(vars(D:lu.div), scale.fun)
temp %<>% mutate('lat.sq' = scale.fun(lat^2), 'wts' = log(nseqs)/mean(log(nseqs)))
temp %<>% select(-taxon,-scale,-nseqs, -ncomps,-n.years, -lu.var) %>%
  mutate_at(vars(pop,year,species), as.factor) %>% as.data.frame

library(tictoc)
tic()
m1 <- bam(div ~ s(lat,long, bs='gp') + s(D, k = 12, bs='cr'), data = temp, family = tw, method='fREML')
toc()

library(parallel)
cl <- makeCluster(2)
m1 <- bam(div ~ s(lat,long) + s(lat.sq) + s(D) + s(hd) + s(hd.var) + s(p.lu) + s(lu.div) + s(year, bs='re') + s(species, bs = 're'),
          data = temp, family = tw, method='fREML',cluster=cl)

tic()
m1 <- bam(div ~ s(lat,long, k = 40) + s(D, k = 12), data = temp, family = tw, method='fREML',cluster=cl)
toc()
summary(m1)
gam.check(m1)
temp$R <- resid(m1)

#spatial autocorrelation in residuals?
temp$rcol <- 4
temp[temp$R < 0, 'rcol'] <- 2
plot(map, xlim = c(-180,180), ylim = c(-90,90),border=NA,col='light gray',axes=F)
points(lat~long,temp,pch=16,col=alpha(temp$rcol,0.5),cex=0.5+abs(temp$R))
spatdat <- select(temp, long,lat,R)
coordinates(spatdat) <- c('long','lat')
variogram(R~long+lat,spatdat,width=0.1,cutoff=50) -> var1
scatter.smooth(x=var1$dist,y=var1$gamma,xlab='distance',ylab='semivariance',ylim=c(0,max(var1$gamma)),pch=16,col='gray')


