## Vincent Fugere 2019

# Code for Figures 4: time series analysis

rm(list=ls())
options(tibble.print_max = 100, scipen = 999)

library(tidyverse)
library(magrittr)
library(RColorBrewer)
library(mgcv)
library(itsadug)
library(scales)
library(viridis)
library(Kendall)

load('~/Google Drive/Recherche/Intraspecific genetic diversity/Data/DF_Master.RData')

#parameters
min.nb.seqs <- 2
min.nb.years <- 4
taxa <- c('birds','fish','insects','mammals')
scl <- '1000'

load('~/Google Drive/Recherche/Intraspecific genetic diversity/Data/temporalGAMMs.RData')
list2env(models,envir=.GlobalEnv)
rm(models)

#### Getting Mann-Kendall coefficient for each time series

MKtau <- function(x){MannKendall(x)$tau}
MKp <- function(x){MannKendall(x)$sl}

mk.coefs <- data.frame()

for(tax in taxa){
  
  temp <- DF %>% filter(scale == scl, taxon == tax) %>%
    filter(nseqs >= min.nb.seqs, div < mean(div)+10*sd(div)) %>%
    add_count(pop) %>%
    select(-n.years) %>%
    rename('n.years' = n) %>%
    select(taxon:ncomps,n.years,lat:lu.div) %>%
    filter(n.years >= min.nb.years) %>%
    mutate('year' = as.numeric(year)) %>%
    group_by(pop) %>%
    summarize('yrs' = median(n.years), 'tau' = MKtau(div), 'p' = MKp(div))
  
  temp$taxon <- tax
  mk.coefs <- bind_rows(mk.coefs,temp)
  rm(temp)
  
}

#### Figure 4 ####

ylims_all <- rbind(c(-9,-3),c(-9,-3),c(-9,-2),c(-7,-2.5))
viridis(100)[1] -> col_ln

pdf('~/Desktop/Fig4.pdf',width=8.5,height=8.5,pointsize = 8)
par(mfrow=c(4,4),cex=1,mar=c(4,4,2,2),oma=c(1,1,1,1))

for(i in 1:4){
  
  tax <- taxa[i]
  
  # MK plot
  
  dat <- filter(mk.coefs, taxon == tax)
  dat$pt.sig <- cut(dat$p, breaks=c(0,0.05,1),labels = 2:1)
  dat$pt.sig <- as.numeric(dat$pt.sig)
  dat$pt.col <- c(col_ln,'black')[dat$pt.sig]
  if(tax == 'insects'){
    dat$pt.alph <- c(0.8,0.05)[dat$pt.sig]
  }else{
  dat$pt.alph <- c(0.8,0.1)[dat$pt.sig]
  }
  
  emptyPlot(xlim = c(-1,1),yaxt='n',xaxt='n',ann=F, ylim=range(dat$yrs)+c(-0.5,0.5),bty='l')
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=seq(4,17,1))
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=seq(-1,1,0.5),labels = c('-1','','0','','1'))
  title(xlab='Mann-Kendall coefficient')
  title(ylab='years in time series')
  abline(v=0,lty=2)
  points(jitter(yrs)~tau,dat,pch=16,col=alpha(pt.col,pt.alph))
  if(tax == 'birds'){legend('topright',inset=c(-0.1,-0.18),bg='white',bty='o',box.col='white',ncol = 1,xpd=NA,pch=21,col=1,pt.bg=c('grey95',col_ln),legend = c('p > 0.05','p < 0.05'))}
  
  # GAMM fit
  
  ylims <- ylims_all[i,]
  tsmod<-get(paste('m1',tax,'1000',sep='_'))
  
  modsum <- summary(tsmod)
  testres <- modsum$s.table[c('s(year)','ti(year,hd)','ti(year,p.lu)'),c('F','p-value')]
  testres[,1] <- round(testres[,1],2)
  testres[,2] <- round(testres[,2],4)
  testres[testres[,2] == 0,2] <- 0.0001
  rsq <- round(modsum$r.sq,2)
  
  emptyPlot(xlim = range(tsmod$model$year),yaxt='n',xaxt='n',ann=F, ylim=ylims,bty='l')
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1,at=seq(-9,-2,1))
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
  title(xlab='year')
  title(ylab=expression(log[e]~COI~diversity~(hat(pi))),line=2.8)
  ifelse(tax == 'insects', assign('ln.alpha',0.05), assign('ln.alpha',0.2))
  for(focalpop in levels(tsmod$model$pop)){
    popdat <- tsmod$model %>% filter(pop == focalpop)
    xs <- seq(min(popdat$year),max(popdat$year),by=1)
    conditions <- list(lat=median(tsmod$model$lat),
                       long=median(tsmod$model$long),
                       D=median(tsmod$model$D),
                       hd=median(tsmod$model$hd),
                       p.lu=median(tsmod$model$p.lu),
                       year=xs,
                       pop=focalpop,
                       order=popdat$order[1],
                       family=popdat$family[1])
    ys <- get_predictions(tsmod, cond = conditions, se=F, print.summary = F) 
    points(fit~year,ys,col=alpha(1,ln.alpha),type='l')
  }
  plot_smooth(tsmod, view="year", lwd=3, col=col_ln, rm.ranef=T, se=1.96, rug=F, add=T)
  mtext(text=bquote(atop(italic('F') == .(testres[1,1]),italic('p') == .(testres[1,2]))),side=3,adj=1)
  legend('bottomright',bty='n',legend=bquote(model~italic(R)^2 == .(rsq)))
  
  # contour plots
  
  newD <- expand.grid(lat=median(tsmod$model$lat),long=median(tsmod$model$long),D=median(tsmod$model$D),
                      year=seq(min(tsmod$model$year),max(tsmod$model$year),length.out = 30),
                      hd=seq(0,1,length.out = 30), p.lu=seq(0,1,length.out = 30), pop=tsmod$model$pop[1],family=tsmod$model$family[1],order=tsmod$model$order[1])
  newD$fit <- predict(tsmod,newD,exclude = c('s(year,pop)','s(order)','s(family)'), discrete=F)
  zlims <- range(newD$fit)
  rm(newD)
  
  fvisgam(tsmod, view = c('year','hd'), cond = list('p.lu' = 0), zlim=zlims,ylab='human density', add.color.legend=T,hide.label=T,xlab = 'year',plot.type = 'contour', lwd=1.5,color = viridis(100), main = NULL,rm.ranef = T,dec=1)
  mtext(text=bquote(italic('F') == .(testres[2,1])~','~italic('p') == .(testres[2,2])),side=3,adj=1)
  
  fvisgam(tsmod, view = c('year','p.lu'), cond = list('hd' = 0), zlim=zlims,ylab='land use',add.color.legend=T,hide.label=T,xlab = 'year',plot.type = 'contour', lwd=1.5,nCol=100,color = viridis(100), main = NULL,rm.ranef = T,dec=1)
  mtext(text=bquote(italic('F') == .(testres[3,1])~','~italic('p') == .(testres[3,2])),side=3,adj=1)
  
}

dev.off()
