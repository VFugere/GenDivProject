#### Spatio-temporal variation in intra-specific neutral genetic diversity across anthromes
#### Gonzalez Lab project - McGill University - 2016-2018
#### Script by Vincent Fugère

rm(list=ls())
options(tibble.print_max = 100, scipen = 999)

#### Part II: figures

#load functions
library(lme4)
library(tidyverse)
library(magrittr)
library(scales)
library(RColorBrewer)
library(cplm)
library(rphylopic)
library(gstat)
library(sp)
library(rworldmap)
library(writexl)
library(broom)
#source('/Users/vincentfugere/Google Drive/Recherche/PhD/R/functions/Zuur.R')
#source('/Users/vincentfugere/Google Drive/Recherche/PhD/R/functions/utils.R')

#from https://www.r-bloggers.com/finding-the-midpoint-when-creating-intervals/
midpoints <- function(x, dp=2){
  lower <- as.numeric(gsub(',.*','',gsub('\\(|\\[|\\)|\\]','', x)))
  upper <- as.numeric(gsub('.*,','',gsub('\\(|\\[|\\)|\\]','', x)))
  return(round(lower+(upper-lower)/2, dp))
}


#load data
load('/Users/vincentfugere/Google Drive/Recherche/Intraspecific genetic diversity/intrasp_gen_div.RData')

#parameters to choose
treshold.yr <- 1980
types.lu <- c('cropland','pasture','urban','conv_rangeland') #to calculate total land use (p.lu) and land use diversity
land.use.var <- 'p.lu' #which land use variable should be included in model (div, or total, or chg - they are all collinear)
human.dens.var <- 'pop_tot' #name of human population density variable to use
file.tag <- 'newmaps' #character string to append to files

#useful for indexing and plotting
#cols <- brewer.pal(4,'Dark2')[c(1,3,4,2)]
cols <- c('#969696','#cccccc','#252525','#636363')
taxa <- c('birds','fish','insects','mammals')
scales <- c('08','1','2','4')
shortax <- c('aves','acti','insect','mam')
map <- getMap(resolution = "coarse")
phylopic.ids <- c('42fdc3cb-37fc-4340-bdf9-eed8e050137c','7a6448e5-09c4-40c8-8378-599d7f974bfe','5aeaf558-3c48-4173-83b4-dbf2846f8d75','b36a215a-adb3-445d-b364-1e63dddd6950')
# #which years to use?
# hist(mam4.agg$year, main='mammals')
# hist(aves4.agg$year, main='birds')
# hist(acti4.agg$year, main='fish')
# hist(insect4.agg$year, main='insects')
# #above 1980 for all, only excluding old fish. LU map goes back to 1980 anyways.

# #should we exclude div estimates with only a few data points?
# par(mfrow=c(4,4))
# for(i in 1:length(scales)){
#   for(j in 1:4){
#     dat <- get(paste(shortax[j],scales[i],'.agg',sep=''))
#     scatter.smooth(y=dat$div,x=dat$nseqs,col=alpha('gray',0.4),log='x',xlab='nb seqs',ylab='div')
#     legend('topright',legend=paste0(shortax[j],'_',scales[i]),bty='n')
#   }
# }
# #no strong relationship between diversity and nb sequences

#gathering all data for figs 1 & 2 time series panels
alldata <- mam08.agg[0,]
for(i in 1:length(scales)){
  for(j in 1:4){
    dat <- get(paste0(shortax[j],scales[i],'.agg')) %>%
      filter(.,!is.na(select(.,human.dens.var))) %>% 
      filter(year >= treshold.yr)
    #remove extreme outliers (10 sd greater than mean)
    dat %<>% filter(div < mean(div)+10*sd(div))
    #remove years with less than 10 populations
    # dat %>% aggregate(div~year, data = ., FUN = 'length') -> y
    # y %<>% filter(div > 9)
    # dat %<>% filter(year %in% y$year)
    dat$scale <- scales[i]
    dat$tax <- taxa[j]
    rbind(alldata,dat) -> alldata
  }
}

# #how many pixels have LU change between 1980 and 1990 and between 2000 and 1990?
# #is it a problem that we are using sometimes 9 years old land use map, when there's no hyde map for a given year?
# 
# #load 1990 map
# load('~/Google Drive/recherche/Intraspecific genetic diversity/data/hyde32_1990_08.RData')
# map1990 <- hyde32_1990_08 %>% filter(!is.na(pop_tot)) %>% mutate_at(vars(lat:long), funs(round(as.numeric(.),3))) %>%
#   mutate('cell' = paste0(.data$lat,'_',.data$long)) %>% select(-lat, -long) %>% select(cell,everything())
# rm(hyde32_1990_08)
# 
# #load 2000 map
# load('~/Google Drive/recherche/Intraspecific genetic diversity/data/hyde32_2000_08.RData')
# map2000 <- hyde32_2000_08 %>% filter(!is.na(pop_tot)) %>% mutate_at(vars(lat:long), funs(round(as.numeric(.),3))) %>%
#   mutate('cell' = paste0(.data$lat,'_',.data$long)) %>% select(-lat, -long) %>% select(cell,everything())
# rm(hyde32_2000_08)
# 
# par(mfrow=c(2,2))
# 
# pixdat <- alldata %>% filter(scale == '08', year %in% (1981:1989)) %>% distinct(cell, .keep_all = T)
# pixdat <- inner_join(pixdat,map1990,by='cell') %>% mutate('lut1' = rowSums(select(.,paste0(types.lu,'.x'))), 'lut2' = rowSums(select(.,paste0(types.lu,'.y')))) 
# plot(lut2~lut1,pixdat, main = '1981-1989',xlab='land use 1980', ylab='land use 1990'); abline(a=0,b=1,lty=2)
# plot(log1p(pop_tot.y)~log1p(pop_tot.x),pixdat, main = '1981-1989',xlab='log human density 1980', ylab='log human density 1990'); abline(a=0,b=1,lty=2)
# 
# pixdat <- alldata %>% filter(scale == '08', year %in% (1991:1999)) %>% distinct(cell, .keep_all = T)
# pixdat <- inner_join(pixdat,map2000,by='cell') %>% mutate('lut1' = rowSums(select(.,paste0(types.lu,'.x'))), 'lut2' = rowSums(select(.,paste0(types.lu,'.y')))) 
# plot(lut2~lut1,pixdat, main = '1991-1999',xlab='land use 1990', ylab='land use 2000'); abline(a=0,b=1,lty=2)
# plot(log1p(pop_tot.y)~log1p(pop_tot.x),pixdat, main = '1991-1999',xlab='log human density 1990', ylab='log human density 2000'); abline(a=0,b=1,lty=2)

alldata %<>% mutate(p.lu = rowSums(select(.,types.lu)))
alldata[alldata$p.lu > 1,'p.lu'] <- 1

alldata %<>% select(scale,tax,pop,species:long,p.lu,pop_tot)

# #for Chloe's map (Fig. 1a)
# fig1data <- alldata %>% group_by(scale,tax,cell,year) %>% summarize(sequences = sum(nseqs))
# save(fig1data, file = '~/Desktop/ChloeDat.Rdata')

#for Table S1 : note that for paper, we 
table1 <- alldata %>% group_by(scale,tax) %>%
  summarize(sequences = sum(nseqs), populations = n_distinct(pop), species = n_distinct(species)) %>%
  rename('taxon' = tax)
#write_xlsx(table1, '~/Desktop/table1.xlsx')

#### figure 1 time series panel ####

#pdf('~/Desktop/1b.pdf', pointsize = 8, width = 4, height=3)

plodat <- alldata %>% filter(scale == '08') %>% 
  group_by(tax,year) %>%
  summarize(npops = n(), nseqs = sum(nseqs)) %>%
  mutate(ptcx = rescale(npops,c(0.5,2.5)))
plot(nseqs~year,plodat,type='n',yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l',ylim=c(0,10.2),xlim=c(1980,2018))
title(ylab='number of sequences', cex.lab=1)
title(xlab='year', cex.lab=1)
axis(2,cex.axis=1,lwd=0,lwd.ticks=1,labels = c(1,10,100,1000,10000),at=log(c(1,10,100,1000,10000)))
axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(1985,1995,2005,2015))
#legs <- c(paste(as.character(min(plodat$npops)),'population'),paste(as.character(max(plodat$npops)),'populations'))
#legend('topleft',x.intersp=3,inset=0.03,bty='n',cex=0.7,pch=16,col=1,pt.cex=c(2.5,0.5),legend=rev(legs),y.intersp = 2.5)
for(j in 1:4){
  y <- filter(plodat, tax == taxa[j]) %>% as.data.frame
  y$nseqs <- log(y$nseqs)
  #points(nseqs~year,y,pch=1,type='p',cex=ptcx,col=alpha(cols[j],1),lwd=0.5)
  points(nseqs~year,y,pch=16,type='l',lwd=1.2,col=alpha(cols[j],1))
}
# xlocs<-seq(from=0.4,to=0.8,length.out = 4)
# ylocs<-log(c(2,25,42,475))/10.2
# ylocs[2] <- ylocs[2] - 0.06
# ylocs[3] <- ylocs[3] + 0.02
# for(j in 1:4){
#   label <- image_data(phylopic.ids[j], size = 128)[[1]]
#   add_phylopic_base(label, x = 0.94, y = ylocs[j], ysize = 0.15, alpha=1,color=cols[j])
# }
#legs <- c('1 population','1500 populations','3500 populations')
#legend('topleft',x.intersp=1.5,inset=c(0.05,-0.02),bty='n',cex=1,pch=16,col=1,pt.cex=rev(c(0.500000,1.345460,2.473491)),legend=rev(legs),y.intersp = 1.2)
#dev.off()

#### figure S1 ####

#pdf('~/Desktop/FigS1.pdf',  pointsize = 8, width = 3.25, height=5.5)
par(mfrow=c(3,1),cex=1,mar=c(4,4,1,1))

# (a). Correlation between sequence number and population/species number
plodat <- alldata %>% filter(scale == '08') %>% distinct(pop, .keep_all = T) %>% 
  group_by(cell) %>%
  summarize(nspecies = n(), nseqs = sum(nseqs))
plot(nspecies~nseqs,plodat,type='n',yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l',log='xy')
title(ylab='species', cex.lab=1)
title(xlab='sequences', cex.lab=1)
axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
points(nspecies~nseqs,plodat,pch=16,cex=0.3,col=alpha(1,0.4))
legend('topleft',legend = bquote(italic(r[S]) == .(round(cor(plodat$nspecies,plodat$nseqs,method='spearman'),2))),bty='n',cex=1)

# (b,c) accumulation of populations and species in dataset over time

plodat <- alldata %>% filter(scale == '08') %>% 
  group_by(tax,year) %>%
  summarize(nspecies = n_distinct(species), npops = n(), nseqs = sum(nseqs))
  
plot(nspecies~year,plodat,type='n',yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l',ylim=c(0,7.5),xlim=c(1980,2018))
title(ylab='number of species', cex.lab=1)
title(xlab='year', cex.lab=1)
axis(2,cex.axis=1,lwd=0,lwd.ticks=1,labels = c(1,10,100,1000),at=log(c(1,10,100,1000)))
axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(1985,1995,2005,2015))
for(j in 1:4){
  y <- filter(plodat, tax == taxa[j]) %>% as.data.frame
  y$nspecies <- log(y$nspecies)
  points(nspecies~year,y,pch=16,type='l',lwd=1.2,col=alpha(cols[j],1))
}

plot(npops~year,plodat,type='n',yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l',ylim=c(0,8.5),xlim=c(1980,2018))
title(ylab='number of sequences', cex.lab=1)
title(xlab='year', cex.lab=1)
axis(2,cex.axis=1,lwd=0,lwd.ticks=1,labels = c(1,10,100,1000),at=log(c(1,10,100,1000)))
axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(1985,1995,2005,2015))
for(j in 1:4){
  y <- filter(plodat, tax == taxa[j]) %>% as.data.frame
  y$npops <- log(y$npops)
  points(npops~year,y,pch=16,type='l',lwd=1.2,col=alpha(cols[j],1))
}

#dev.off()

#### figure 2 ####

# maps of genetic diversity
#pdf('~/Desktop/2a.pdf',width=5,height=10,pointsize = 6)

par(mfrow=c(4,1),oma=c(0,0,0,0),mar=c(0,0,0,0),cex=1)
#cols2<-rev(c("#7A1E07", "#DE0D0D", "#F57A08", "#F5ED03", "#9FE685", "#40B316", "#4388F0", "#0000C4"))
#cols3<-rev(c("#A60000", "#FF1500", "#FF9100", "#FAFA32", "#B9FF40", "#2EA610", "#03769C", "#003D61"))
cols3<-rev(c('#d73027','#f46d43','#fdae61','#fee090','#e0f3f8','#abd9e9','#74add1','#4575b4'))

max.div <- c(0.008,0.016,0.024,0.016)

for(j in 1:4){
  dat <- get(paste(shortax[j],scales[4],'.agg',sep=''))
  dat %<>% group_by(cell) %>% summarize_at(vars(lat,long,div), funs(mean)) %>% as.data.frame(.)
  # #8 quantiles
  #probfunc <- ecdf(quantile(dat$div[dat$div !=0], probs = seq(0,1, length.out = 7)))
  # dat$dq <- 0
  # dat$dq[dat$div != 0] <- probfunc(dat$div[dat$div!=0])
  # dat$dq <- rescale(dat$dq, to = c(1,8))
  # #8 equal intervals of 0.01
  dat$div[dat$div > max.div[j]] <- max.div[j]
  dat$dq <- as.numeric(cut_interval(dat$div, 8))
  dat$col <- cols3[dat$dq]
  plot(map, xlim = c(-180,180), ylim = c(-90,90),col=alpha('gray',0.3),border=NA,asp=1.2,axes=F,cex.lab=0.5)
  symbols(x=dat$long,y=dat$lat,squares=rep(4,length(dat$lat)),inches=F,fg=0,bg=dat$col,add=T,lwd=0.2)
  polygon(x=c(-180,180,180,-180),y=c(-60,-60,-90,-90),col='white',border=NA)
  # #legend on the lower left side
  # quants <- as.character(round(quantile(dat$div[dat$div !=0], probs = seq(0,1, length.out = 7)),4))
  # legs <- c(paste('<',quants[1]),paste(quants[1],'-',quants[2]),paste(quants[2],'-',quants[3]),paste(quants[3],'-',quants[4]),paste(quants[4],'-',quants[5]),paste(quants[5],'-',quants[6]),paste(quants[6],'-',quants[7]),paste('>',quants[7]))
  # legend(x=-175,y=-10,box.col = 'white',bg='white',cex=0.5,pt.cex=1,pch=15,col=rev(cols3),legend=rev(legs),y.intersp = 1)
  # legend(x=-180,y=-10,bty='n',cex=0.5,pt.cex=1,pch=15,col=rev(cols3),legend=rep('',8),y.intersp = 1)
  # legend(x=-194,y=2,bty='n',cex=0.7,legend='Genetic diversity')
  # #legend('bottom',horiz=T,bty='n',cex=0.7,pt.cex=seq(0.5,5.5,length.out=8),legend=c(round(min(dat$div),3),rep('',6),round(max(dat$div),3)),pch=16,inset=0.03,col=alpha(4,0.3),x.intersp=1.2)
  # label <- image_data(phylopic.ids[j], size = 128)[[1]]
  # add_phylopic_base(label, x = 0.9, y = 0.8, ysize = 0.15, alpha=1,color=1)
  # horizontal legend below africa
  legs <- as.character(seq(0,max.div[j]*7/8,length.out = 8))
  legs[8] <- paste('>',max.div[j]*7/8)
  xseqs <- seq(-15,60,length.out = 8)
  rect(xleft=xseqs,xright=xseqs+(xseqs[2]-xseqs[1]),ybottom = rep(-55,8),ytop = rep(-50,8),col=cols3,border=NULL,lwd=0.2)
  segments(x0=xseqs[2:8],x1=xseqs[2:8],y0=rep(-50,7),y1=c(-59,-57,-57,-59,-57,-57,-59),lwd=c(0.5,0.3,0.3,0.5,0.3,0.3,0.5))
  text(x=xseqs[c(2,5,8)],y=rep(-59,3),labels = legs[c(2,5,8)],pos=1)
  text(x=27.85714,y=-49,cex=1,label='Genetic diversity',pos=3)
}
#mtext(text=paste('scale =',scales[i],'x',scales[i]), side=3,outer=T,line=1)
#dev.off()

# for(j in 1:4){
#   dat <- get(paste(shortax[j],scales[2],'.agg',sep=''))
#   dat %<>% group_by(cell) %>% summarize_at(vars(lat,long,div), funs(mean)) %>% as.data.frame(.)
#   probfunc <- ecdf(quantile(dat$div[dat$div !=0], probs = seq(0,1, length.out = 7)))
#   dat$dq <- 0
#   dat$dq[dat$div != 0] <- probfunc(dat$div[dat$div!=0])
#   dat$dq <- rescale(dat$dq, to = c(1,8))
#   dat$col <- cols3[dat$dq]
#   plot(map, xlim = c(-180,180), ylim = c(-90,90),col=alpha('gray',0.3),border=NA,asp=1.2,axes=F,cex.lab=0.5)
#   symbols(x=dat$long,y=dat$lat,squares=rep(1,length(dat$lat)),inches=F,fg=0,bg=dat$col,add=T,lwd=0.1)
#   quants <- as.character(round(quantile(dat$div[dat$div !=0], probs = seq(0,1, length.out = 7)),4))
#   legs <- c(paste('<',quants[1]),paste(quants[1],'-',quants[2]),paste(quants[2],'-',quants[3]),paste(quants[3],'-',quants[4]),paste(quants[4],'-',quants[5]),paste(quants[5],'-',quants[6]),paste(quants[6],'-',quants[7]),paste('>',quants[7]))
#   legend(x=-180,y=0,box.col = 'white',bg='white',cex=0.5,pt.cex=1,pch=22,col=0,pt.bg=rev(cols3),legend=rev(legs),y.intersp = 1)
#   label <- image_data(phylopic.ids[j], size = 128)[[1]]
#   add_phylopic_base(label, x = 0.9, y = 0.8, ysize = 0.15, alpha=1,color=1)
# }

# #panel b) global time series
# #pdf('~/Desktop/Fig2b.pdf',width=2.5,height = 7, pointsize=6)
# #par(mfrow=c(4,1),oma=c(3,3,0,0),mar=c(2,2,1,1),cex=1)
# par(mfrow=c(4,1),oma=c(0,0,0,0),mar=c(4,4,1,1),cex=1)
# for(j in 1:4){
#   y <- filter(alldata, scale == '4', tax == taxa[j]) %>% group_by(year) %>%
#   summarise(mean = mean(div),se = sd(div)/sqrt(n()), n = n()) %>% filter(!is.na(se)) %>% as.data.frame()
#   plot(mean~year,y,type='n',ylim=c(0,0.025),yaxt='n',xaxt='n',ann=F,bty='n')
#   box(bty='l')
#   axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
#   axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
#   polygon(x=c(y$year,rev(y$year)),y=c((y$mean-1.96*y$se),rev(y$mean+1.96*y$se)),col=alpha(1,0.3),border=NA)
#   #points(mean~year,y,pch=1,type='o',cex=rescale(y$n,c(0.5,2.5)))
#   points(mean~year,y,pch=1,type='l')
#   #text(x=y$year,y=y$mean*1.4,labels=as.character(y$n),pos=3,cex=0.6,srt=90)
#   #legs <- paste(as.character(c(min(y$n),max(y$n))),'populations')
#   #legend('topleft',inset=0.01,box.col = alpha('white',0.5),box.lwd=0,bg='white',pch=16,col=1,pt.cex=c(0.5,2.5),legend=legs,y.intersp = 2.5)
#   #label <- image_data(phylopic.ids[j], size = 128)[[1]]
#   #add_phylopic_base(label, x = 0.5, y = 0.9, ysize = 0.25, alpha=1,color='black')
#   title(ylab='mean genetic diversity', cex.lab=1)
#   title(xlab='year', cex.lab=1)
# }
# #dev.off()

#pdf('~/Desktop/Fig2b.pdf',width=2.5,height = 7, pointsize=8)
par(mfrow=c(4,1),oma=c(0,0,0,0),mar=c(4,4,1,1),cex=1)
for(j in 1:4){
  y <- filter(alldata, scale == '4', tax == taxa[j]) %>% group_by(year) %>%
    summarise(mean = mean(div),ci = 1.96*(sd(div)/sqrt(n())), n = n()) %>% filter(!is.na(ci)) %>% as.data.frame()
  plot(mean~year,y,type='n',yaxt='n',xaxt='n',ann=F,bty='l',ylim=range(c(y$mean-y$ci,y$mean+y$ci)))
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
  y$lwr <- y$mean - y$ci
  y$upr <- y$mean + y$ci
  arrows(x0=y$year,y0=y$lwr,y1=y$upr,length=0,lwd=1.5,col=alpha(1,0.4))
  points(mean~year,y,pch=21,cex=1,col=1,bg=alpha(cols[j],0.8))
  lines(smooth.spline(y=y$mean,x=y$year,spar=0.6),lwd=3,col=alpha(1,0.4))
  label <- image_data(phylopic.ids[j], size = 128)[[1]]
  add_phylopic_base(label, x = 0.85, y = 0.9, ysize = 0.2, alpha=1,color=1)
  title(ylab='mean genetic diversity', cex.lab=1)
  title(xlab='year', cex.lab=1)
}
#dev.off()

#### Fig. S2: other drivers of gen div ####

#pdf('~/Desktop/FigS2.pdf',pointsize = 8, width=8,height=8)
par(mfrow=c(4,4),oma=c(0,2.5,0,0),cex=1,mar=c(4,2,1,1))

#latitude
for(j in 1:4){
  y <- filter(alldata, scale == '4', tax == taxa[j]) %>%
    mutate(lat = abs(lat)) %>% group_by(lat) %>%
    summarise(mean = mean(div),ci = 1.96*(sd(div)/sqrt(n())), n = n()) %>% filter(!is.na(ci)) %>% as.data.frame()
  plot(mean~lat,y,type='n',yaxt='n',xaxt='n',ann=F,bty='l',ylim=c(0,range(c(y$mean-y$ci,y$mean+y$ci))[2]))
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
  y$lwr <- y$mean - y$ci
  y$upr <- y$mean + y$ci
  arrows(x0=y$lat,y0=y$lwr,y1=y$upr,length=0,lwd=1.5,col=alpha(1,0.4))
  points(mean~lat,y,pch=21,cex=1,col=1,bg=alpha(cols[j],0.8))
  lines(smooth.spline(y=y$mean,x=y$lat,spar=0.6),lwd=3,col=alpha(1,0.4))
  label <- image_data(phylopic.ids[j], size = 128)[[1]]
  add_phylopic_base(label, x = 0.85, y = 0.9, ysize = 0.2, alpha=1,color=1)
  title(xlab='absolute latitude (degrees)', cex.lab=1.3, line=2.7)
}

#longitude
for(j in 1:4){
  y <- filter(alldata, scale == '4', tax == taxa[j]) %>% 
    mutate(long = midpoints(cut(long,30))) %>% group_by(long) %>%
    summarise(mean = mean(div),ci = 1.96*(sd(div)/sqrt(n())), n = n()) %>%
    filter(!is.na(ci)) %>% as.data.frame()
  plot(mean~long,y,type='n',yaxt='n',xaxt='n',ann=F,bty='l',ylim=c(0,range(c(y$mean-y$ci,y$mean+y$ci))[2]))
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
  y$lwr <- y$mean - y$ci
  y$upr <- y$mean + y$ci
  arrows(x0=y$long,y0=y$lwr,y1=y$upr,length=0,lwd=1.5,col=alpha(1,0.4))
  points(mean~long,y,pch=21,cex=1,col=1,bg=alpha(cols[j],0.8))
  lines(smooth.spline(y=y$mean,x=y$long,spar=0.6),lwd=3,col=alpha(1,0.4))
  label <- image_data(phylopic.ids[j], size = 128)[[1]]
  add_phylopic_base(label, x = 0.85, y = 0.9, ysize = 0.2, alpha=1,color=1)
  title(xlab='longitude (degrees)', cex.lab=1.3, line=2.7)
}

#land use
for(j in 1:4){
  y <- filter(alldata, scale == '4', tax == taxa[j]) %>% 
    mutate(p.lu = midpoints(cut(p.lu,30))) %>% group_by(p.lu) %>%
    summarise(mean = mean(div),ci = 1.96*(sd(div)/sqrt(n())), n = n()) %>%
    filter(!is.na(ci)) %>% as.data.frame()
  plot(mean~p.lu,y,type='n',yaxt='n',xaxt='n',ann=F,bty='l',ylim=c(0,range(c(y$mean-y$ci,y$mean+y$ci))[2]))
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
  y$lwr <- y$mean - y$ci
  y$upr <- y$mean + y$ci
  arrows(x0=y$p.lu,y0=y$lwr,y1=y$upr,length=0,lwd=1.5,col=alpha(1,0.4))
  points(mean~p.lu,y,pch=21,cex=1,col=1,bg=alpha(cols[j],0.8))
  lines(smooth.spline(y=y$mean,x=y$p.lu,spar=0.6),lwd=3,col=alpha(1,0.4))
  label <- image_data(phylopic.ids[j], size = 128)[[1]]
  add_phylopic_base(label, x = 0.85, y = 0.9, ysize = 0.2, alpha=1,color=1)
  title(xlab='land use (proportion)', cex.lab=1.3, line=2.7)
}

#human density
for(j in 1:4){
  y <- filter(alldata, scale == '4', tax == taxa[j]) %>% 
    mutate(pop_tot = midpoints(cut(log1p(pop_tot),30))) %>% group_by(pop_tot) %>%
    summarise(mean = mean(div),ci = 1.96*(sd(div)/sqrt(n())), n = n()) %>%
    filter(!is.na(ci)) %>% as.data.frame()
  plot(mean~pop_tot,y,type='n',yaxt='n',xaxt='n',ann=F,bty='l',ylim=c(0,range(c(y$mean-y$ci,y$mean+y$ci))[2]))
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
  y$lwr <- y$mean - y$ci
  y$upr <- y$mean + y$ci
  arrows(x0=y$pop_tot,y0=y$lwr,y1=y$upr,length=0,lwd=1.5,col=alpha(1,0.4))
  points(mean~pop_tot,y,pch=21,cex=1,col=1,bg=alpha(cols[j],0.8))
  lines(smooth.spline(y=y$mean,x=y$pop_tot,spar=0.6),lwd=3,col=alpha(1,0.4))
  label <- image_data(phylopic.ids[j], size = 128)[[1]]
  add_phylopic_base(label, x = 0.85, y = 0.9, ysize = 0.2, alpha=1,color=1)
  title(xlab='human density (log1+x)', cex.lab=1.3, line=2.7)
}

mtext(text='genetic diversity', cex=1.3,side=2,line=1,outer=T)
#dev.off()

#### Fig. SX: illustrating scale dependence ####
# 
# pdf('~/Desktop/FigS3.pdf',pointsize = 8, width=8,height=8)
# par(mfrow=c(4,4),oma=c(0,2.5,0,0),cex=1,mar=c(4,2,1,1))
# 
# #latitude
# for(i in 1:4){
#   for(j in 1:4){
#     y <- filter(alldata, scale == scales[i], tax == taxa[j]) %>%
#       mutate(lat = abs(lat)) %>% group_by(lat) %>%
#       summarise(mean = mean(div),ci = 1.96*(sd(div)/sqrt(n())), n = n()) %>% filter(!is.na(ci)) %>% as.data.frame()
#     plot(mean~lat,y,type='n',yaxt='n',xaxt='n',ann=F,bty='l',ylim=c(0,range(c(y$mean-y$ci,y$mean+y$ci))[2]))
#     axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
#     axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
#     y$lwr <- y$mean - y$ci
#     y$upr <- y$mean + y$ci
#     arrows(x0=y$lat,y0=y$lwr,y1=y$upr,length=0,lwd=1.5,col=alpha(1,0.4))
#     points(mean~lat,y,pch=21,cex=1,col=1,bg=alpha(cols[j],0.8))
#     lines(smooth.spline(y=y$mean,x=y$lat,spar=0.6),lwd=3,col=alpha(1,0.4))
#     label <- image_data(phylopic.ids[j], size = 128)[[1]]
#     add_phylopic_base(label, x = 0.85, y = 0.9, ysize = 0.2, alpha=1,color=1)
#     title(xlab='absolute latitude (degrees)', cex.lab=1.3, line=2.7)
#   }
# }
# mtext(text='genetic diversity', cex=1.3,side=2,line=1,outer=T)
# dev.off()
# 
# #land use
# for(i in 1:4){
#   for(j in 1:4){
#     y <- filter(alldata, scale == scales[i], tax == taxa[j]) %>% 
#       mutate(p.lu = midpoints(cut(p.lu,30))) %>% group_by(p.lu) %>%
#       summarise(mean = mean(div),ci = 1.96*(sd(div)/sqrt(n())), n = n()) %>%
#       filter(!is.na(ci)) %>% as.data.frame()
#     plot(mean~p.lu,y,type='n',yaxt='n',xaxt='n',ann=F,bty='l',ylim=c(0,range(c(y$mean-y$ci,y$mean+y$ci))[2]))
#     axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
#     axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
#     y$lwr <- y$mean - y$ci
#     y$upr <- y$mean + y$ci
#     arrows(x0=y$p.lu,y0=y$lwr,y1=y$upr,length=0,lwd=1.5,col=alpha(1,0.4))
#     points(mean~p.lu,y,pch=21,cex=1,col=1,bg=alpha(cols[j],0.8))
#     lines(smooth.spline(y=y$mean,x=y$p.lu,spar=0.6),lwd=3,col=alpha(1,0.4))
#     label <- image_data(phylopic.ids[j], size = 128)[[1]]
#     add_phylopic_base(label, x = 0.85, y = 0.9, ysize = 0.2, alpha=1,color=1)
#     title(xlab='land use (proportion)', cex.lab=1.3, line=2.7)
#   }
# }
# mtext(text='genetic diversity', cex=1.3,side=2,line=1,outer=T)
# dev.off()
# 
# #human density
# for(i in 1:4){
#   for(j in 1:4){
#     y <- filter(alldata, scale == scales[i], tax == taxa[j]) %>% 
#       mutate(pop_tot = midpoints(cut(log1p(pop_tot),30))) %>% group_by(pop_tot) %>%
#       summarise(mean = mean(div),ci = 1.96*(sd(div)/sqrt(n())), n = n()) %>%
#       filter(!is.na(ci)) %>% as.data.frame()
#     plot(mean~pop_tot,y,type='n',yaxt='n',xaxt='n',ann=F,bty='l',ylim=c(0,range(c(y$mean-y$ci,y$mean+y$ci))[2]))
#     axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
#     axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
#     y$lwr <- y$mean - y$ci
#     y$upr <- y$mean + y$ci
#     arrows(x0=y$pop_tot,y0=y$lwr,y1=y$upr,length=0,lwd=1.5,col=alpha(1,0.4))
#     points(mean~pop_tot,y,pch=21,cex=1,col=1,bg=alpha(cols[j],0.8))
#     lines(smooth.spline(y=y$mean,x=y$pop_tot,spar=0.6),lwd=3,col=alpha(1,0.4))
#     label <- image_data(phylopic.ids[j], size = 128)[[1]]
#     add_phylopic_base(label, x = 0.85, y = 0.9, ysize = 0.2, alpha=1,color=1)
#     title(xlab='human density (log1+x)', cex.lab=1.3, line=2.7)
#   }
# }
# mtext(text='genetic diversity', cex=1.3,side=2,line=1,outer=T)
# dev.off()

#### creating a clean dataset for model selection script for global analysis ####
 
#init dataframe to receive all data used in global model
modeldata <- mam08.agg[0,]

for(i in 1:length(scales)){
  
  for(j in 1:4){

    dat <- get(paste0(shortax[j],scales[i],'.agg')) %>%
      filter(.,!is.na(select(.,human.dens.var))) %>%
      filter(year >= treshold.yr) %>%
      mutate('p.lu' = rowSums(select(.,types.lu))) %>%
      mutate('p.wild' = 1-p.lu)

    # Adding some quality control

    # i) no less than x sequence comparisons to get measure of div
    # dat %<>% filter(nseqs > 2)

    # ii) remove extreme outliers (x sd greater than mean)
    dat %<>% filter(div < mean(div)+10*sd(div))

    # iii) remove duplicated pops because cannot fit a random effect for pop for insects (R crashes..)
    dat %<>% distinct(pop, .keep_all = T)

    # iv.a) remove yrs with < pops than 0.5 X µ nb of pops across all years (Dornelas criterion)
    #dat %>% aggregate(div~year, data = ., FUN = 'length') -> y
    #y %<>% filter(div > mean(y$div)*0.5)
    #dat %<>% filter(year %in% y$year)
    # iv.b) the dornelas criterion removes yr with 35 pops for data-poorest scale/taxon combination (mam @ 4x4). This is a little strict, so instead remove yrs < 30 pops
    #dat %>% aggregate(div~year, data = ., FUN = 'length') -> y
    #y %<>% filter(div >= 30)
    #dat %<>% filter(year %in% y$year)
    # iv.c) at a minium, remove years with less than 10 populations
    dat %>% aggregate(div~year, data = ., FUN = 'length') -> y
    y %<>% filter(div > 9)
    dat %<>% filter(year %in% y$year)

    #correcting land use with proportions exceeding 1
    dat[dat$p.lu > 1,'p.wild'] <- 0
    dat[dat$p.lu > 1,'p.lu'] <- 1

    dat$scale <- scales[i]
    dat$tax <- taxa[j]

    # #adding land use diversity  (Shannon exponent)
    # ludiv <- numeric(0)
    # for(k in 1:nrow(dat)){
    #   dtemp <- dat[k,c(types.lu,'p.wild')]
    #   dtemp <- dtemp[which(dtemp != 0)]
    #   ludiv[k] <- exp(-(sum(dtemp*log(dtemp))))
    # }
    # dat$lu.div <- ludiv
    
    rbind(modeldata,dat) -> modeldata
  }
}

# #save data
# save(modeldata, file = '~/Google Drive/Recherche/Intraspecific genetic diversity/modeldata.Rdata')

#### Model selection with sequential deletion of non significant terms ####

# see other script named 'gendiv_modelselection'

## Table S1: model selection results

#reading models/script output
load('~/Google Drive/Recherche/Intraspecific genetic diversity/CPGLMMs/models08.Rdata')
load('~/Google Drive/Recherche/Intraspecific genetic diversity/CPGLMMs/models1.Rdata')
load('~/Google Drive/Recherche/Intraspecific genetic diversity/CPGLMMs/models2.Rdata')
load('~/Google Drive/Recherche/Intraspecific genetic diversity/CPGLMMs/models4.Rdata')

#one large object with all models
models <- append(models08, c(models1, models2, models4))

#pulling out the coefficients and SEs and making a long df
modseltbl <- models %>% map_df(~ rownames_to_column(as.data.frame(summary(.)$coefs))) %>%
  rename(coef = rowname, value = Estimate, se = `Std. Error`) %>% select(coef:se) %>%
  filter(coef != '(Intercept)')

#adding a column with model name, for spread function below
nb.coef.per.mod <- unlist(models %>% map(~ nrow(summary(.)$coefs))) - 1
mod.name.vec <- data.frame('model' = rep(seq(from=1,to=length(models)), nb.coef.per.mod))
modseltbl <- bind_cols(mod.name.vec, modseltbl)

#reformat values and table to large
modseltbl <- modseltbl %>% mutate_at(vars(value:se), funs(round(.,2))) %>%
  mutate('values' = paste0(value,' (',se,')')) %>% select(-value, -se)

tableS1 <- modseltbl %>% group_by(model) %>% spread(key = coef, value = values, fill = NA) %>%
  select(1,4,9,14,2,13,6,8,5,7,12,10,11,15,16,3)
#write_xlsx(tableS1, '~/Desktop/TableS1.xlsx')

#### Caterpillar plot (Figure 3) ####

# #panel b
#
# #caterpillar plot for hypothesis-driven model (4 effects only)
#
# coefs$scale <- c(rep(0.08,nrow(results)),rep(1,nrow(results)),rep(2,nrow(results)),rep(4,nrow(results)))
#
# #coefs <- read_csv('~/Google Drive/Recherche/Intraspecific genetic diversity/models_p.lu.csv')
#
# coefs$y.lwr <- coefs$year - coefs$year.ci
# coefs$y.upr <- coefs$year + coefs$year.ci
# coefs$lat.lwr <- coefs$lat - coefs$lat.ci
# coefs$lat.upr <- coefs$lat + coefs$lat.ci
# coefs$l.lwr <- coefs$lu - coefs$lu.ci
# coefs$l.upr <- coefs$lu + coefs$lu.ci
# coefs$h.lwr <- coefs$humans - coefs$humans.ci
# coefs$h.upr <- coefs$humans + coefs$humans.ci
# # coefs$l.y.lwr <- coefs$lu.yr - coefs$lu.yr.ci
# # coefs$l.y.upr <- coefs$lu.yr + coefs$lu.yr.ci
# # coefs$h.y.lwr <- coefs$hu.yr - coefs$hu.yr.ci
# # coefs$h.y.upr <- coefs$hu.yr + coefs$hu.yr.ci
#
# #write_csv(coefs, path = paste0('~/Google Drive/Recherche/Intraspecific genetic diversity/models_',land.use.var,'_',file.tag,'.csv'))
#
# #global plotting parameters
# xlims <- c(-max(abs(coefs[,11:18])),max(abs(coefs[,11:18])))
# #labLs <- c(expression(italic(Mammalia)),expression(italic(Aves)),expression(italic(Actinopterygii)),expression(italic(Insecta)))
#
# #pdf for figure
# pdf(file = paste0('~/Google Drive/Recherche/Intraspecific genetic diversity/catplot_',land.use.var,'_',file.tag,'.pdf'),width = 6, height = 5)
# par(mfrow=c(2,2),oma=c(2,2,0,0),mar=c(2,2,1,1))
#
# for(i in 1:4){
#
#   subdat <- filter(coefs, taxon == taxa[i])
#
#   uprs <- subdat %>% select(scale, ends_with('.upr')) %>% gather(key = effect, value = value, -scale) %>% arrange(scale)
#   lwrs <- subdat %>% select(scale, ends_with('.lwr')) %>% gather(key = effect, value = value, -scale) %>% arrange(scale)
#   means <- subdat %>% select(scale, year, lat, lu, humans) %>% gather(key = effect, value = value, -scale) %>% arrange(scale)
#
#   df <- cbind(means,lwrs$value,uprs$value)
#   colnames(df)[4:5] <- c('lwr','upr')
#   df$ptcol <- cols[i]
#   df$ptcol[df$upr*df$lwr < 0] <- 'white'
#   df$lncol <- 'dark gray'
#   df$lncol[df$upr*df$lwr < 0] <- 'light gray'
#
#   plot(0,type='n',yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l',ylim=c(0.5,16.5),xlim=xlims)
#   axis(2,cex.axis=1,lwd=0,lwd.ticks=0,at=seq(2.5,17,4),labels = c('0.08','1','2','4'))
#   axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
#   abline(h=c(4.5,8.5,12.5),col=alpha(1,0.4),lty=3,lwd=0.5)
#   abline(v=0,lty=2,col='black')
#
#   arrows(x0=df$lwr,x1=df$upr,y0=1:16,length=0,lwd=1.5,col=df$lncol)
#   points(y=1:16,x=means$value,col=1,bg=df$ptcol,pch=rep(21:24,4))
#
#   label <- image_data(phylopic.ids[i], size = 128)[[1]]
#   add_phylopic_base(label, x = 0.15, y = 0.9, ysize = 0.25, alpha=0.8,color=cols[i])
#
#   if(i==1){
#     #legend('topright',box.col = 'white',cex=0.7,pt.cex=0.8,legend=rev(c('year','latitude','land use','humans','year:land use','year:humans')),pch=c(8,25:21),pt.bg=1,y.intersp=1,bg='white')
#     legend('topright',box.col = 'white',cex=0.7,pt.cex=0.7,legend=rev(c('year','latitude','land use','humans')),pch=c(24:21),pt.bg=1,y.intersp=1,bg='white')
#   }
# }
#
# mtext(text='coefficient (effect on genetic diversity)', side=1,outer=T,line=1)
# mtext(text="spatial scale (grid cell size in degrees)", side=2,outer=T,line=1)
#
# dev.off()
#
# rm(dat,lwrs,means,model,resid.ts,results,spatdat,subdat,uprs,var1,y,xlims,r2,line,label,j,i,df)
# #save.image('~/Google Drive/Recherche/Intraspecific genetic diversity/alldata_image.RData')
# #save(alldata,file = '~/Google Drive/Recherche/Intraspecific genetic diversity/alldata.RData')

# #panel b


#caterpillar plot for model with all two way interactions, using models fitted already

coefs <- data.frame('scale' = numeric(0),
                    'tax' = character(0),
                    'year' = numeric(0),
                    'lat' = numeric(0),
                    'long' = numeric(0),
                    'LU' = numeric(0),
                    'HD' = numeric(0),
                    'year:lat' = numeric(0),
                    'year:long' = numeric(0),
                    'year:LU' = numeric(0),
                    'year:HD' = numeric(0),
                    'lat:long' = numeric(0),
                    'lat:LU' = numeric(0),
                    'lat:HD' = numeric(0),
                    'long:LU' = numeric(0),
                    'long:HD' = numeric(0),
                    'scale' = numeric(0))

models <- get(paste0('models',scales[i]))
fixef(models[[1]])

coefs$scale <- c(rep(0.08,nrow(results)),rep(1,nrow(results)),rep(2,nrow(results)),rep(4,nrow(results)))

coefs$y.lwr <- coefs$year - coefs$year.ci
coefs$y.upr <- coefs$year + coefs$year.ci
coefs$lat.lwr <- coefs$lat - coefs$lat.ci
coefs$lat.upr <- coefs$lat + coefs$lat.ci
coefs$l.lwr <- coefs$lu - coefs$lu.ci
coefs$l.upr <- coefs$lu + coefs$lu.ci
coefs$h.lwr <- coefs$humans - coefs$humans.ci
coefs$h.upr <- coefs$humans + coefs$humans.ci
# coefs$l.y.lwr <- coefs$lu.yr - coefs$lu.yr.ci
# coefs$l.y.upr <- coefs$lu.yr + coefs$lu.yr.ci
# coefs$h.y.lwr <- coefs$hu.yr - coefs$hu.yr.ci
# coefs$h.y.upr <- coefs$hu.yr + coefs$hu.yr.ci

#write_csv(coefs, path = paste0('~/Google Drive/Recherche/Intraspecific genetic diversity/models_',land.use.var,'_',file.tag,'.csv'))

#global plotting parameters
xlims <- c(-max(abs(coefs[,11:18])),max(abs(coefs[,11:18])))
#labLs <- c(expression(italic(Mammalia)),expression(italic(Aves)),expression(italic(Actinopterygii)),expression(italic(Insecta)))

#pdf for figure
pdf(file = paste0('~/Google Drive/Recherche/Intraspecific genetic diversity/catplot_',land.use.var,'_',file.tag,'.pdf'),width = 6, height = 5)
par(mfrow=c(2,2),oma=c(2,2,0,0),mar=c(2,2,1,1))

for(i in 1:4){

  subdat <- filter(coefs, taxon == taxa[i])

  uprs <- subdat %>% select(scale, ends_with('.upr')) %>% gather(key = effect, value = value, -scale) %>% arrange(scale)
  lwrs <- subdat %>% select(scale, ends_with('.lwr')) %>% gather(key = effect, value = value, -scale) %>% arrange(scale)
  means <- subdat %>% select(scale, year, lat, lu, humans) %>% gather(key = effect, value = value, -scale) %>% arrange(scale)

  df <- cbind(means,lwrs$value,uprs$value)
  colnames(df)[4:5] <- c('lwr','upr')
  df$ptcol <- cols[i]
  df$ptcol[df$upr*df$lwr < 0] <- 'white'
  df$lncol <- 'dark gray'
  df$lncol[df$upr*df$lwr < 0] <- 'light gray'

  plot(0,type='n',yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l',ylim=c(0.5,16.5),xlim=xlims)
  axis(2,cex.axis=1,lwd=0,lwd.ticks=0,at=seq(2.5,17,4),labels = c('0.08','1','2','4'))
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
  abline(h=c(4.5,8.5,12.5),col=alpha(1,0.4),lty=3,lwd=0.5)
  abline(v=0,lty=2,col='black')

  arrows(x0=df$lwr,x1=df$upr,y0=1:16,length=0,lwd=1.5,col=df$lncol)
  points(y=1:16,x=means$value,col=1,bg=df$ptcol,pch=rep(21:24,4))

  label <- image_data(phylopic.ids[i], size = 128)[[1]]
  add_phylopic_base(label, x = 0.15, y = 0.9, ysize = 0.25, alpha=0.8,color=cols[i])

  if(i==1){
    #legend('topright',box.col = 'white',cex=0.7,pt.cex=0.8,legend=rev(c('year','latitude','land use','humans','year:land use','year:humans')),pch=c(8,25:21),pt.bg=1,y.intersp=1,bg='white')
    legend('topright',box.col = 'white',cex=0.7,pt.cex=0.7,legend=rev(c('year','latitude','land use','humans')),pch=c(24:21),pt.bg=1,y.intersp=1,bg='white')
  }
}

mtext(text='coefficient (effect on genetic diversity)', side=1,outer=T,line=1)
mtext(text="spatial scale (grid cell size in degrees)", side=2,outer=T,line=1)

dev.off()

#### Time series analysis ####

#load('~/Google Drive/Recherche/Intraspecific genetic diversity/alldata_image.RData')

pdf('~/Desktop/Fig4.pdf',width=3.25,pointsize = 6,height=3)

par(mfrow=c(2,2), mar=c(4,4,1,1),oma=c(0,0,0,0),cex=1)

i <- 4

ymaxes <- c(0.032,0.032,0.032,0.06)
ymins <- c(0.0005,0.0001,0.0001,0.0001)

#init dataframe to receive all data
plodat <- data.frame('p.lu' = numeric(0), 'pop_tot' = numeric(0), 'slopes' = numeric(0), 'tax' = character(0), stringsAsFactors = F)

TSmods <- list()

for(j in c(1,4)){
  
  dat <- get(paste0(shortax[j],scales[i],'.agg')) %>%
    filter(.,!is.na(select(.,human.dens.var))) %>% 
    filter(year >= treshold.yr) %>%
    mutate('p.lu' = rowSums(select(.,types.lu))) %>%
    mutate('p.wild' = 1-p.lu)
  
  # remove extreme outliers (x sd greater than mean)
  dat %<>% filter(div < mean(div)+10*sd(div))
  
  #only keeping time series with 3+ time points
  dat <- dat[dat$n.years >= 3,]
  dat <- droplevels(dat)
  
  #correcting land use with proportions exceeding 1
  dat[dat$p.lu > 1,'p.wild'] <- 0
  dat[dat$p.lu > 1,'p.lu'] <- 1
  
  #choosing relevant land use column for model
  dat$lu <- subset(dat, select = land.use.var)[,1]
  
  dat$s.yr <- as.numeric(scale(dat$year))
  dat$s.nyr <- as.numeric(scale(dat$n.years))
  dat$lu <- as.numeric(scale(dat$lu))
  dat$sc.pop <- as.numeric(scale(log1p(subset(dat, select =human.dens.var))[,1]))
  dat$salat <- as.numeric(scale(abs(dat$lat))) #absolute value of latitude
  dat$snseqs <- as.numeric(scale(log(dat$nseqs)))
  dat$slong <- as.numeric(scale(dat$long))
  
  #checking colinearity
  corvif(dat[,c('s.yr','snseqs','lu','salat','slong','sc.pop')])
  
  #checking trend per ts duration
  #model <- cpglmm(div ~ 0 + s.yr * as.factor(n.years) + (1+s.yr|pop) + (1|cell), data=dat, weights = log(nseqs))
  #write.csv(as.data.frame(summary(model)@coefs), file = '~/Desktop/insects.csv')
  
  #model: weighing or not, adding correlated intercept/slope or not, does not change anything
  #went with weights & uncoorrelated slope/intercept, which is perhaps cleanest
  model <- cpglmm(div ~ 1 + s.yr + (1+s.yr|pop) + (1|cell), data=dat, weights = log(nseqs))
  #model <- cpglmm(div ~ 1 + s.yr + (0+s.yr|pop) + (1|pop), data=dat)
  #model <- cpglmm(div ~ 1 + s.yr + (1+s.yr|pop), data=dat, weights = log(nseqs))
  #model <- cpglmm(div ~ 1 + s.yr + (1+s.yr|pop), data=dat)
  #fullmodel <- cpglmm(div ~ 1 + s.yr*(lu + sc.pop) + (0+s.yr|pop) + (1|pop) + (1|cell), data=dat, weights=log(nseqs))
  #summary(fullmodel)
  #TSmods <- c(TSmods, fullmodel)
  #model.ME <- cpglmm(div ~ 1 + s.yr + lu + sc.pop + salat + slong + (0+s.yr|pop) + (1|pop), data=dat, weights=log(nseqs))
  
  intercept <- fixef(model)[1]
  ef <-  fixef(model)[2]
  lwr <- fixef(model)[2] - summary(model)$coefs[2,2]*1.96
  upr <- fixef(model)[2] + summary(model)$coefs[2,2]*1.96
  
  yr1 <- dat$year[which(dat$s.yr == min(dat$s.yr))[1]]
  yr2 <- dat$year[which(dat$s.yr == max(dat$s.yr))[1]]
  
  #sketching model
  plot(x=0,y=0.01,dat,type='n',yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l',xlim=range(dat$s.yr),ylim=c(ymins[j],ymaxes[j]),log='y')
  title(ylab='genetic diversity', cex.lab=1)
  title(xlab='year', cex.lab=1)
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=seq(min(dat$s.yr),max(dat$s.yr),length.out = 5),labels=round(seq(yr1,yr2,length.out = 5),0))
  n <- nlevels(dat$pop)
  for(k in 1:n){
    pop.data <- dat[dat$pop == levels(dat$pop)[k],]
    #points(div~s.yr,pop.data,col=alpha(1,0.2),type='o',pch=16)
    ys <- predict(model, pop.data)
    points(ys~pop.data$s.yr,type='l',col=alpha(1,0.05))
  }
  y1 <- as.numeric(exp(intercept + ef*min(dat$s.yr)))
  y2 <- as.numeric(exp(intercept + ef*max(dat$s.yr)))
  segments(x0=min(dat$s.yr),x1=max(dat$s.yr),y0=y1,y1=y2,lwd=2)
  
  #plotting random slopes
  df <- dat %>% group_by(pop) %>% summarize_at(vars(p.lu,pop_tot,n.years), funs(mean, max, min)) %>% as.data.frame(.)
  ranef(model)$pop[,'s.yr'] + ef -> df$slopes
  
  h <- hist(df$slopes, breaks = 30, plot = F)
  histmax <- max(h$counts)
  
  plot(h, col = 1, main=NULL, xlab='temporal trend (random slope)',ylab='populations',las=1,border=0,ylim=c(-(histmax/12),histmax))
  ypos <- -(histmax/12)*0.7
  segments(x0=lwr,x1=upr,y0=ypos,y1=ypos,lty=1,col=1)
  points(x=ef,y=ypos,pch=16,col=1,cex=1.2)
  
  label <- image_data(phylopic.ids[j], size = 128)[[1]]
  add_phylopic_base(label, x = 0.85, y = 0.9, ysize = 0.25, alpha=1,color=1)
  
  df$tax <- taxa[j]
  df %<>% select(-pop)
  
  plodat <- bind_rows(plodat,df)
  
  # #validation
  # dat$R <- resid(model, type='standardized')
  # dat$Yfit <- fitted(model)
  # 
  # #how good is the model fit?
  # plot(Yfit~div,dat,ylab='fitted',xlab='observed',pch=16,col=alpha(1,0.2),cex=log(nseqs)/2)
  # abline(a=0,b=1,lty=2)
  # abline(lm(Yfit~div,dat))
  # round(summary(lm(Yfit~div,dat))$adj.r.squared,4) -> r2
  # legend('topleft',legend=paste0('(',r2,')'),bty='n')
  # 
  # #heteroscedasticity?
  # plot(R~Yfit,dat,ylab='standardized residuals',xlab='fitted',pch=16,col=alpha(1,0.2))
  
}

plodat$col.idx <- 1
plodat$col.idx[plodat$tax == 'insects'] <- 4
plodat$logHD <- log1p(plodat$pop_tot_mean)
plodat$lu.chg <- plodat$p.lu_max - plodat$p.lu_min
plodat$HD.chg <-  plodat$pop_tot_max - plodat$pop_tot_min

dev.off()

pdf('~/Desktop/TS_supp.pdf',width=3.25,pointsize = 6,height=6)
par(mfrow=c(4,2), mar=c(4,4,1,1),oma=c(0,0,0,0),cex=1)

for(z in 4:7){

  for(j in c(1,4)){
    
    dat <- get(paste0(shortax[j],scales[i],'.agg')) %>%
      filter(.,!is.na(select(.,human.dens.var))) %>% 
      filter(year >= treshold.yr) %>%
      mutate('p.lu' = rowSums(select(.,types.lu))) %>%
      mutate('p.wild' = 1-p.lu)
    
    # remove extreme outliers (x sd greater than mean)
    dat %<>% filter(div < mean(div)+10*sd(div))
    
    #only keeping time series with z time points (4-7)
    dat <- dat[dat$n.years >= z,]
    dat <- droplevels(dat)
    
    #correcting land use with proportions exceeding 1
    dat[dat$p.lu > 1,'p.wild'] <- 0
    dat[dat$p.lu > 1,'p.lu'] <- 1
    
    #choosing relevant land use column for model
    dat$lu <- subset(dat, select = land.use.var)[,1]
    
    dat$s.yr <- as.numeric(scale(dat$year))
    dat$s.nyr <- as.numeric(scale(dat$n.years))
    dat$lu <- as.numeric(scale(dat$lu))
    dat$sc.pop <- as.numeric(scale(log1p(subset(dat, select =human.dens.var))[,1]))
    dat$salat <- as.numeric(scale(abs(dat$lat))) #absolute value of latitude
    dat$snseqs <- as.numeric(scale(log(dat$nseqs)))
    dat$slong <- as.numeric(scale(dat$long))
    
    model <- cpglmm(div ~ 1 + s.yr + (1+s.yr|pop) + (1|cell), data=dat, weights = log(nseqs))
    
    ef <-  fixef(model)[2]
    lwr <- fixef(model)[2] - summary(model)$coefs[2,2]*1.96
    upr <- fixef(model)[2] + summary(model)$coefs[2,2]*1.96
    
    #plotting random slopes
    df <- dat %>% group_by(pop) %>% summarize_at(vars(p.lu,pop_tot,n.years), funs(mean, max, min)) %>% as.data.frame(.)
    ranef(model)$pop[,'s.yr'] -> df$slopes
    
    h <- hist(df$slopes, breaks = 30, plot = F)
    histmax <- max(h$counts)
    
    plot(h, col = cols[j], main=NULL, xlab='temporal trend (random slope)',ylab='populations',las=1,border=0,ylim=c(-(histmax/12),histmax))
    ypos <- -(histmax/12)*0.7
    segments(x0=lwr,x1=upr,y0=ypos,y1=ypos,lty=1,col=cols[j])
    points(x=ef,y=ypos,pch=16,col=cols[j],cex=1.2)
    
    label <- image_data(phylopic.ids[j], size = 128)[[1]]
    add_phylopic_base(label, x = 0.85, y = 0.9, ysize = 0.25, alpha=1,color=1)
    
     }
}

dev.off()

# plot(slopes~p.lu_mean,plodat,type='n',yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l')
# title(xlab='proportion of land used', cex.lab=1)
# title(ylab='population trend', cex.lab=1)
# axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
# axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
# abline(h=0,lty=2)
# points(slopes~p.lu_mean,plodat,pch=c(0,1)[col.idx^0.5],cex=0.5,col=alpha(cols[col.idx],0.2))
# legend('topleft',inset=c(0,-0.05),legend = bquote(italic(r[S]) == .(round(cor(plodat$p.lu_mean,plodat$slopes,method='spearman'),2))),bty='n',cex=1)
# 
# plot(slopes~logHD,plodat,type='n',yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l')
# title(xlab=expression(log[10](1+human~density)), cex.lab=1)
# title(ylab='population trend', cex.lab=1)
# axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
# axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
# abline(h=0,lty=2)
# points(slopes~logHD,plodat,pch=c(0,1)[col.idx^0.5],cex=0.5,col=alpha(cols[col.idx],0.2))
# legend('topleft',inset=c(0,-0.05),legend = bquote(italic(r[S]) == .(round(cor(plodat$logHD,plodat$slopes,method='spearman'),2))),bty='n',cex=1)
# 
dev.off()
# 
# pdf('~/Desktop/Fig4_extra.pdf',width=3.75,pointsize = 6,height=1.5)
# 
# par(mfrow=c(1,2), mar=c(4,4,1,1),oma=c(0,0,0,0),cex=1)
# 
# plot(slopes~lu.chg,plodat,type='n',yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l')
# title(xlab='land use change (proportion)', cex.lab=1)
# title(ylab='population trend', cex.lab=1)
# axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
# axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
# abline(h=0,lty=2)
# points(slopes~lu.chg,plodat,pch=c(0,1)[col.idx^0.5],cex=0.5,col=alpha(cols[col.idx],0.2))
# legend('topleft',inset=c(0,-0.05),legend = bquote(italic(r[S]) == .(round(cor(plodat$lu.chg,plodat$slopes,method='spearman'),2))),bty='n',cex=1)
# 
# plot(slopes~HD.chg,plodat,type='n',yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l')
# title(xlab=human~density~increase, cex.lab=1)
# title(ylab='population trend', cex.lab=1)
# axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
# axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
# abline(h=0,lty=2)
# points(slopes~HD.chg,plodat,pch=c(0,1)[col.idx^0.5],cex=0.5,col=alpha(cols[col.idx],0.2))
# legend('topleft',inset=c(0,-0.05),legend = bquote(italic(r[S]) == .(round(cor(plodat$HD.chg,plodat$slopes,method='spearman'),2))),bty='n',cex=1)
# 
# dev.off()

#### w ####

pdf('~/Desktop/TSfunnel.pdf',width=7,pointsize = 12,height=4)

par(mfrow=c(1,2), mar=c(4,4,1,1),oma=c(0,0,0,0),cex=1)

plot(slopes~n.years_mean,subset(plodat, tax == 'mammals'),type='n',yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l')
title(xlab='duration', cex.lab=1)
title(ylab='population trend', cex.lab=1)
axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
abline(h=0,lty=2)
points(slopes~jitter(n.years_mean),subset(plodat, tax == 'mammals'),pch=1,cex=1,col=alpha(1,0.2))
legend('topright',legend = 'mammals',bty='n',cex=1)

plot(slopes~n.years_mean,subset(plodat, tax == 'insects'),type='n',yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l')
title(xlab='duration', cex.lab=1)
title(ylab='population trend', cex.lab=1)
axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
abline(h=0,lty=2)
points(slopes~jitter(n.years_mean),subset(plodat, tax == 'insects'),pch=1,cex=1,col=alpha(1,0.2))
legend('topright',legend = 'insects',bty='n',cex=1)

dev.off()
