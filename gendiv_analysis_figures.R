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
library(randomcoloR)
#source('/Users/vincentfugere/Google Drive/Recherche/PhD/R/functions/Zuur.R')
#source('/Users/vincentfugere/Google Drive/Recherche/PhD/R/functions/utils.R')

#from https://www.r-bloggers.com/finding-the-midpoint-when-creating-intervals/
midpoints <- function(x, dp=2){
  lower <- as.numeric(gsub(',.*','',gsub('\\(|\\[|\\)|\\]','', x)))
  upper <- as.numeric(gsub('.*,','',gsub('\\(|\\[|\\)|\\]','', x)))
  return(round(lower+(upper-lower)/2, dp))
}

#VIF functions from Zuur et al. 2009
corvif <- function(dataz) {
  dataz <- as.data.frame(dataz)
  #correlation part
  cat("Correlations of the variables\n\n")
  tmp_cor <- cor(dataz,use="complete.obs")
  print(tmp_cor)
  #vif part
  form    <- formula(paste("fooy ~ ",paste(strsplit(names(dataz)," "),collapse=" + ")))
  dataz   <- data.frame(fooy=1,dataz)
  lm_mod  <- lm(form,dataz)
  cat("\n\nVariance inflation factors\n\n")
  print(myvif(lm_mod))
}

#Support function for corvif. Will not be called by the user
myvif <- function(mod) {
  v <- vcov(mod)
  assign <- attributes(model.matrix(mod))$assign
  if (names(coefficients(mod)[1]) == "(Intercept)") {
    v <- v[-1, -1]
    assign <- assign[-1]
  } else warning("No intercept: vifs may not be sensible.")
  terms <- labels(terms(mod))
  n.terms <- length(terms)
  if (n.terms < 2) stop("The model contains fewer than 2 terms")
  if (length(assign) > dim(v)[1] ) {
    diag(tmp_cor)<-0
    if (any(tmp_cor==1.0)){
      return("Sample size is too small, 100% collinearity is present")
    } else {
      return("Sample size is too small")
    }
  }
  R <- cov2cor(v)
  detR <- det(R)
  result <- matrix(0, n.terms, 3)
  rownames(result) <- terms
  colnames(result) <- c("GVIF", "Df", "GVIF^(1/2Df)")
  for (term in 1:n.terms) {
    subs <- which(assign == term)
    result[term, 1] <- det(as.matrix(R[subs, subs])) * det(as.matrix(R[-subs, -subs])) / detR
    result[term, 2] <- length(subs)
  }
  if (all(result[, 2] == 1)) {
    result <- data.frame(GVIF=result[, 1])
  } else {
    result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
  }
  invisible(result)
}
#END VIF FUNCTIONS


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
#pdf('~/Desktop/Fig2.pdf',width=5,height=10,pointsize = 6)

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

#### Fig 3: div as a function of drivers ####

#pdf('~/Desktop/Fig3.pdf',pointsize = 8, width=8,height=10)
par(mfrow=c(5,4),oma=c(0,2.5,0,0),cex=1,mar=c(4,2,1,1))

scl <- '4'
wigl <- 0.7
lncol <- c("#00329E")

#time
for(j in 1:4){
  y <- filter(alldata, scale == scl, tax == taxa[j]) %>% group_by(year) %>%
    summarise(mean = mean(div),ci = 1.96*(sd(div)/sqrt(n())), n = n()) %>% filter(!is.na(ci)) %>% as.data.frame()
  ymax <- max(y$mean)
  yticks <- seq(0,ymax,length.out = 5)
  ytickslab <- round(seq(0,ymax,length.out = 3),3)
  ytickslab <- c(ytickslab[1],'',ytickslab[2],'',ytickslab[3])
  plot(mean~year,y,type='n',yaxt='n',xaxt='n',ann=F,bty='l',xlim=c(1980,2016),ylim=c(0,ymax))
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1, at=yticks, labels = ytickslab)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=seq(1980,2016,length.out = 5))
  y$lwr <- y$mean - y$ci
  y$upr <- y$mean + y$ci
  arrows(x0=y$year,y0=y$lwr,y1=y$upr,length=0,lwd=1.5,col=alpha(1,0.4))
  points(mean~year,y,pch=21,cex=1,col=1,bg=alpha(cols[j],0.8))
  lines(smooth.spline(y=y$mean,x=y$year,spar=wigl,w=log(y$n)),lwd=3,col=alpha(lncol,0.6))
  #label <- image_data(phylopic.ids[j], size = 128)[[1]]
  #add_phylopic_base(label, x = 0.85, y = 0.9, ysize = 0.2, alpha=1,color=1)
  title(xlab='year', cex.lab=1.3, line = 2.7)
}

#latitude
for(j in 1:4){
  y <- filter(alldata, scale == scl, tax == taxa[j]) %>%
    mutate(lat = abs(lat)) %>% group_by(lat) %>%
    summarise(mean = mean(div),ci = 1.96*(sd(div)/sqrt(n())), n = n()) %>%
    filter(!is.na(ci)) %>% as.data.frame()
  ymax <- max(y$mean)
  yticks <- seq(0,ymax,length.out = 5)
  ytickslab <- round(seq(0,ymax,length.out = 3),3)
  ytickslab <- c(ytickslab[1],'',ytickslab[2],'',ytickslab[3])
  plot(mean~lat,y,type='n',yaxt='n',xaxt='n',ann=F,bty='l',ylim=c(0,max(y$mean)))
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1, at=yticks, labels = ytickslab)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
  y$lwr <- y$mean - y$ci
  y$upr <- y$mean + y$ci
  arrows(x0=y$lat,y0=y$lwr,y1=y$upr,length=0,lwd=1.5,col=alpha(1,0.4))
  points(mean~lat,y,pch=21,cex=1,col=1,bg=alpha(cols[j],0.8))
  lines(smooth.spline(y=y$mean,x=y$lat,spar=wigl,w=log(y$n)),lwd=3,col=alpha(lncol,0.6))
  # label <- image_data(phylopic.ids[j], size = 128)[[1]]
  # add_phylopic_base(label, x = 0.85, y = 0.9, ysize = 0.2, alpha=1,color=1)
  title(xlab='absolute latitude (degrees)', cex.lab=1.3, line=2.7)
}

#longitude
for(j in 1:4){
  y <- filter(alldata, scale == scl, tax == taxa[j]) %>% 
    mutate(long = midpoints(cut(long,30))) %>% group_by(long) %>%
    summarise(mean = mean(div),ci = 1.96*(sd(div)/sqrt(n())), n = n()) %>%
    filter(!is.na(ci)) %>% as.data.frame()
  ymax <- max(y$mean)
  yticks <- seq(0,ymax,length.out = 5)
  ytickslab <- round(seq(0,ymax,length.out = 3),3)
  ytickslab <- c(ytickslab[1],'',ytickslab[2],'',ytickslab[3])
  plot(mean~long,y,type='n',yaxt='n',xaxt='n',ann=F,bty='l',ylim=c(0,max(y$mean)))
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1, at=yticks, labels = ytickslab)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
  y$lwr <- y$mean - y$ci
  y$upr <- y$mean + y$ci
  arrows(x0=y$long,y0=y$lwr,y1=y$upr,length=0,lwd=1.5,col=alpha(1,0.4))
  points(mean~long,y,pch=21,cex=1,col=1,bg=alpha(cols[j],0.8))
  lines(smooth.spline(y=y$mean,x=y$long,spar=wigl,w=log(y$n)),lwd=3,col=alpha(lncol,0.6))
  # label <- image_data(phylopic.ids[j], size = 128)[[1]]
  # add_phylopic_base(label, x = 0.85, y = 0.9, ysize = 0.2, alpha=1,color=1)
  title(xlab='longitude (degrees)', cex.lab=1.3, line=2.7)
}

#land use
for(j in 1:4){
  y <- filter(alldata, scale == scl, tax == taxa[j]) %>% 
    mutate(p.lu = midpoints(cut(p.lu,30))) %>% group_by(p.lu) %>%
    summarise(mean = mean(div),ci = 1.96*(sd(div)/sqrt(n())), n = n()) %>%
    filter(!is.na(ci)) %>% as.data.frame()
  ymax <- max(y$mean)
  yticks <- seq(0,ymax,length.out = 5)
  ytickslab <- round(seq(0,ymax,length.out = 3),3)
  ytickslab <- c(ytickslab[1],'',ytickslab[2],'',ytickslab[3])
  plot(mean~p.lu,y,type='n',yaxt='n',xaxt='n',ann=F,bty='l',ylim=c(0,max(y$mean)))
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1, at=yticks, labels = ytickslab)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
  y$lwr <- y$mean - y$ci
  y$upr <- y$mean + y$ci
  arrows(x0=y$p.lu,y0=y$lwr,y1=y$upr,length=0,lwd=1.5,col=alpha(1,0.4))
  points(mean~p.lu,y,pch=21,cex=1,col=1,bg=alpha(cols[j],0.8))
  lines(smooth.spline(y=y$mean,x=y$p.lu,spar=wigl,w=log(y$n)),lwd=3,col=alpha(lncol,0.6))
  # label <- image_data(phylopic.ids[j], size = 128)[[1]]
  # add_phylopic_base(label, x = 0.85, y = 0.9, ysize = 0.2, alpha=1,color=1)
  title(xlab='land use (proportion)', cex.lab=1.3, line=2.7)
}

#human density
for(j in 1:4){
  y <- filter(alldata, scale == scl, tax == taxa[j]) %>% 
    mutate(pop_tot = midpoints(cut(log1p(pop_tot),30))) %>% group_by(pop_tot) %>%
    summarise(mean = mean(div),ci = 1.96*(sd(div)/sqrt(n())), n = n()) %>%
    filter(!is.na(ci)) %>% as.data.frame()
  ymax <- max(y$mean)
  yticks <- seq(0,ymax,length.out = 5)
  ytickslab <- round(seq(0,ymax,length.out = 3),3)
  ytickslab <- c(ytickslab[1],'',ytickslab[2],'',ytickslab[3])
  plot(mean~pop_tot,y,type='n',yaxt='n',xaxt='n',ann=F,bty='l',ylim=c(0,max(y$mean)))
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1, at=yticks, labels = ytickslab)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
  y$lwr <- y$mean - y$ci
  y$upr <- y$mean + y$ci
  arrows(x0=y$pop_tot,y0=y$lwr,y1=y$upr,length=0,lwd=1.5,col=alpha(1,0.4))
  points(mean~pop_tot,y,pch=21,cex=1,col=1,bg=alpha(cols[j],0.8))
  lines(smooth.spline(y=y$mean,x=y$pop_tot,spar=wigl,w=log(y$n)),lwd=3,col=alpha(lncol,0.6))
  # label <- image_data(phylopic.ids[j], size = 128)[[1]]
  # add_phylopic_base(label, x = 0.85, y = 0.9, ysize = 0.2, alpha=1,color=1)
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

rm(dat)

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

# #### Testing for some visible non-linear effects (Fig 3) with GAMMs
# 
# library(mgcv)
# d <- modeldata %>% filter(scale == '4', tax == 'mammals')
# d$lhd <- log1p(d$pop_tot)
# m <- gam(div ~ s(lhd), data=d)
# summary(m); plot(m)
# gam.check(m)

#### Caterpillar plot (Figure 4) ####

##making the coef matrices: one each for coef values, lwr bound, and upr bound

#useful vectors
mod.kp <- c(3,6,9,13,16,19,22,26,29,32,35,38,41,44,47,50)
tax.scale <- expand.grid('taxon' = c('mam','bird','fish','insect'),
                        'scale' = c(0.08,1,2,4))
colvec <- rep(rev(cols), each = 4)

#getting coefs
modseltbl <- models %>% map_df(~ rownames_to_column(as.data.frame(summary(.)$coefs))) %>%
  rename(coef = rowname, value = Estimate, se = `Std. Error`) %>% select(coef:se) %>%
  filter(coef != '(Intercept)') %>% bind_cols(mod.name.vec)

coefs <- modseltbl %>% select(-se) %>% group_by(model) %>% spread(key = coef, value = value, fill = NA) %>%
  select(1,4,9,14,2,13,6,8,5,7,12,10,11,15,16,3) %>% filter(model %in% mod.kp) %>%
  ungroup %>% select(-model)

se <- modseltbl %>% select(-value) %>% mutate(se = se*1.96) %>% group_by(model) %>%
  spread(key = coef, value = se, fill = NA) %>%
  select(1,4,9,14,2,13,6,8,5,7,12,10,11,15,16,3) %>% filter(model %in% mod.kp) %>%
  ungroup %>% select(-model)

#re-arranging in right order from birds to mammals top to bottom
coefs <- coefs[c(1,5,9,13,4,8,12,16,3,7,11,15,2,6,10,14),]
se <- se[c(1,5,9,13,4,8,12,16,3,7,11,15,2,6,10,14),]
coefs <- as.data.frame(coefs)
lwr <- coefs - se
upr <- coefs + se

##figure 4

pdf(file = '~/Desktop/Fig4.pdf',width=8.5,height=3,pointsize=7)
par(mfrow=c(1,1),mar=c(4,3,1,1),cex=1,oma=c(0,0,0,0))
plot(0,type='n',yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='n',xlim=c(0,15.1),ylim=c(0.5,16.5),yaxs='i', xaxs='i')
for(i in seq(0.5,15.5,by=2)){
  polygon(x=c(0,15,15,0),y=c(i,i,i+1,i+1),col='light gray',border=NA)
}
axis(2,cex.axis=1,lwd=0,lwd.ticks=0,at=1:16,labels = rep(c("5'",'1°','2°','4°'),4),las=1)
axis(2,cex.axis=1,lwd=1,lwd.ticks=1,at=c(0.5,4.5,8.5,12.5,16.5), labels = rep('',5))
abline(v=1:15,lty=1,lwd=0.5)
abline(v=seq(from=0.5,to=14.5,by=1),lty=2,lwd=0.5)
abline(h=c(4.5,8.5,12.5))
axis(1,cex.axis=1,lwd=0,lwd.ticks=0,at=seq(from=0.5,to=14.5,by=1),labels=rep(0,15),line=-0.75)
axis(1,cex.axis=0.7,lwd=1,lwd.ticks=0.5,at=seq(from=0,to=15,by=1),labels=rep('',16))
axis(1,cex.axis=1,lwd=0,lwd.ticks=0,at=seq(from=0.5,to=14.5,by=1),labels=c('YR','LAT','LONG','LU','HD','YR:LAT','YR:LONG','YR:LU','YR:HD','LAT:LONG','LAT:LU','LAT:HD','LONG:LU','LONG:HD','LU:HD'),line=1)
for(i in 1:15){
  for(j in 1:16){
    pt <- coefs[j,i] + i-0.5
    xg <- lwr[j,i] + i-0.5
    xd <- upr[j,i] + i-0.5
    if(is.na(pt)){
      ptfill <- NA
      lncol <- 'gray'
      lnwd <- 1
    } else if(0 < upr[j,i] & 0 > lwr[j,i]){
      ptfill <-  'white'
      lncol <- 'black'
      lnwd <- 1.25
    } else {
      ptfill <- 'black'
      lncol <- 'black'
      lnwd <- 2
    }
    arrows(x0=xg,x1=xd,y0=j,length=0,col=lncol,lwd=lnwd)
    points(x=pt,y=j,pch=21,col=1,bg=ptfill,cex=1.2)
  }
}
dev.off()

#####

rm(models, models08, models1, models2, models4, modeldata)

#### Time series analysis ####

#how many time series per scale/tax combinations?

table1$ts <- numeric(16)

for(i in 1:4){
  for(j in 1:4){
    dat <- get(paste0(shortax[j],scales[i],'.agg')) %>%
      distinct(pop, .keep_all = T) %>% as.data.frame 
    table1$ts[((i-1)*4)+j] <- sum(dat$n.years > 2)
  }
}

#### Fig. 4 ####

#init objects to receive metadata (slopes/time series properties) and models
tsdat <- data.frame()
tsmodels <- list()

#plotting and indexing helpers
i <- 4
ymaxes <- c(0.032,0.032,0.06,0.032)
ymins <- c(0.0001,0.0001,0.0001,0.0005)

#pdf('~/Desktop/Fig4.pdf',width=4,pointsize = 8,height=7)
par(mfrow=c(4,2), mar=c(4,4,1,1),oma=c(0,0,0,0),cex=1)

for(j in 1:4){
  
  dat <- get(paste0(shortax[j],scales[i],'.agg')) %>%
    filter(.,!is.na(select(.,human.dens.var))) %>% 
    filter(year >= treshold.yr) %>%
    mutate('p.lu' = rowSums(select(.,types.lu))) %>%
    mutate('p.wild' = 1-p.lu)
  
  # remove extreme outliers (x sd greater than mean)
  dat %<>% filter(div < mean(div)+10*sd(div))
  
  #recomputing number of time points after excluding outliers and old data points
  dat <- dat %>% select(-n.years) %>% add_count(pop) %>% rename(n.years = n)
  
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
  
  # #checking colinearity
  # corvif(dat[,c('s.yr','s.nyr','lu','salat','sc.pop')])
  
  #model: weighing or not, adding correlated intercept/slope or not, does not change anything
  model <- cpglmm(div ~ 1 + s.yr + (s.yr|pop) + (1|cell) + (1|species), data=dat, weights = log(nseqs))
  #model <- cpglmm(div ~ 1 + s.yr + (s.yr-1|pop) + (1|pop), data=dat)
  #model <- cpglmm(div ~ 1 + s.yr + (1+s.yr|pop), data=dat, weights = log(nseqs))
  #model <- cpglmm(div ~ 1 + s.yr + (1+s.yr|pop), data=dat)
  
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
  randcols <- distinctColorPalette(k = n)
  for(k in 1:n){
    pop.data <- dat[dat$pop == levels(dat$pop)[k],]
    #points(div~s.yr,pop.data,col=alpha(1,0.2),type='o',pch=16)
    ys <- predict(model, pop.data)
    points(ys~pop.data$s.yr,type='l',col=alpha(randcols[k],0.3))
  }
  y1 <- as.numeric(exp(intercept + ef*min(dat$s.yr)))
  y2 <- as.numeric(exp(intercept + ef*max(dat$s.yr)))
  segments(x0=min(dat$s.yr),x1=max(dat$s.yr),y0=y1,y1=y2,lwd=2)
  
  #plotting random slopes
  df <- dat %>% group_by(pop) %>% summarize_at(vars(p.lu,pop_tot,n.years), funs(mean, max, min)) %>% as.data.frame(.)
  ranef(model)$pop[,'s.yr'] + ef -> df$slopes
  
  h <- hist(df$slopes, breaks = 30, plot = F)
  histmax <- max(h$counts)
  histrange <- range(h$breaks)
  
  plot(h, col = 1, main=NULL, xlab='temporal trend (slope)',ylab='populations',las=1,border=0,ylim=c(-(histmax/12),histmax),xaxt="n")
  axis(1, at=seq(from=-2,to=2,by=0.5))
  ypos <- -(histmax/12)*0.7
  segments(x0=lwr,x1=upr,y0=ypos,y1=ypos,lty=1,col=1)
  points(x=ef,y=ypos,pch=16,col=1,cex=1.2)
  
  # label <- image_data(phylopic.ids[j], size = 128)[[1]]
  # add_phylopic_base(label, x = 0.85, y = 0.9, ysize = 0.25, alpha=1,color=1)
  # 
  df$tax <- taxa[j]
  df %<>% select(-pop)
  
  tsdat <- bind_rows(tsdat,df)
  
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
  
  #more complex model with other predictors
  model2 <- cpglmm(div ~ 1 + s.yr*salat + (s.yr|pop) + (1|cell) + (1|species), data=dat, weights = log(nseqs))
  model3 <- cpglmm(div ~ 1 + s.yr*lu + (s.yr|pop) + (1|cell) + (1|species), data=dat, weights = log(nseqs))
  model4 <- cpglmm(div ~ 1 + s.yr*sc.pop + (s.yr|pop) + (1|cell) + (1|species), data=dat, weights = log(nseqs))
  
  tsmodels <- append(tsmodels,c(model,model2,model3,model4))

  # HAD TO DITCH THIS PART BECAUSE HIGHLY COLINEAR PREDICTORS THAT CANNOT BE INCLUDED IN THE SAME MODEL
  # #getting coefs
  # coefs <- rownames_to_column(as.data.frame(summary(model2)$coefs)) %>%
  #   rename(coef = rowname, value = Estimate, se = `Std. Error`) %>% select(coef:se) %>%
  #   filter(coef != '(Intercept)')
  # coefs %<>% mutate('lwr' = value-1.96*se, 'upr' = value+1.96*se)
  # coefs <- coefs[c(4,3,5,2,1),]
  # 
  # xmin <- min(coefs$lwr)
  # xmax <- max(coefs$upr)
  # 
  # #caterpillar plot of time series model
  # plot(0,type='n',yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l',xlim=c(xmin,xmax),ylim=c(0.5,5.5))
  # # for(i in seq(0.5,5.5,by=2)){
  # #   polygon(x=c(xmin*1.1,xmax*1.1,xmax*1.1,xmin*1.1),y=c(i,i,i+1,i+1),col='#BACBED',border=NA)
  # # }
  # axis(2,cex.axis=1,lwd=0,lwd.ticks=0,at=1:5,labels = rev(c('YR','DUR','LAT','LU','HD')),las=1)
  # abline(v=0,lty=2,lwd=1)
  # abline(h=1:5,lty=3,lwd=0.25)
  # title(xlab='parameter estimate')
  # axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
  # for(w in 1:5){
  #   pt <- coefs[w,'value']
  #   lwr <- coefs[w,'lwr']
  #   upr <- coefs[w,'upr']
  #   if(0 < upr & 0 > lwr){
  #     ptfill <-  'white'
  #     lncol <- '#CACACA'
  #     lnwd <- 3
  #   } else {
  #     ptfill <- 'black'
  #     lncol <- 'black'
  #     lnwd <- 3
  #   }
  #   segments(x0=lwr,x1=upr,y0=w,y1=w,col=lncol,lwd=lnwd)
  #   points(x=pt,y=w,pch=21,col=1,bg=ptfill,cex=1.5)
  # }
}

# tsdat$logHD <- log1p(tsdat$pop_tot_mean)
# tsdat$lu.chg <- tsdat$p.lu_max - tsdat$p.lu_min
# tsdat$HD.chg <-  tsdat$pop_tot_max - tsdat$pop_tot_min

#dev.off()

#### Table S2 ####

#pulling out the coefficients and SEs and making a long df
tblS2 <- tsmodels %>% map_df(~ rownames_to_column(as.data.frame(summary(.)$coefs))) %>%
  rename(coef = rowname, value = Estimate, se = `Std. Error`) %>% select(coef:se) %>%
  filter(coef != '(Intercept)')

#adding a column with model name, for spread function below
nb.coef.per.mod <- unlist(tsmodels %>% map(~ nrow(summary(.)$coefs))) - 1
mod.name.vec <- data.frame('model' = rep(seq(from=1,to=length(tsmodels)), nb.coef.per.mod))
tblS2 <- bind_cols(mod.name.vec, tblS2)

#reformat values and table to large
tblS2 <- tblS2 %>% mutate_at(vars(value:se), funs(round(.,2))) %>%
  mutate('values' = paste0(value,' (',se,')')) %>% select(-value, -se)

tblS2 <- tblS2 %>% group_by(model) %>% spread(key = coef, value = values, fill = NA) %>%
  select(1,3,7,5,2,4,8,6)
#write_xlsx(tblS2, '~/Desktop/TableS2.xlsx')

#### re-doing analysis, but excluding time series less than x years ####

#getting coefficients

results <- data.frame('tax' = character(0), 
                      'min.yr' = numeric(0),
                      'es' = numeric(0),
                      'lwr' = numeric(0),
                      'upr' = numeric(0), stringsAsFactors = F)

for(z in 3:7){
  
  for(j in 1:4){
    
    dat <- get(paste0(shortax[j],scales[i],'.agg')) %>%
      filter(.,!is.na(select(.,human.dens.var))) %>% 
      filter(year >= treshold.yr)
    
    # remove extreme outliers (x sd greater than mean)
    dat %<>% filter(div < mean(div)+10*sd(div))
    
    #recomputing number of time points after excluding outliers and old data points
    dat <- dat %>% select(-n.years) %>% add_count(pop) %>% rename(n.years = n)
    
    #only keeping time series with z time points (4-7)
    dat <- dat[dat$n.years >= z,]
    dat <- droplevels(dat)
    
    dat$s.yr <- as.numeric(scale(dat$year))
    dat$s.nyr <- as.numeric(scale(dat$n.years))
    dat$snseqs <- as.numeric(scale(log(dat$nseqs)))
    
    model <- cpglmm(div ~ 1 + s.yr + (s.yr|pop) + (1|cell) + (1|species), data=dat, weights = log(nseqs))
    
    ef <-  fixef(model)[2]
    lwr <- fixef(model)[2] - summary(model)$coefs[2,2]*1.96
    upr <- fixef(model)[2] + summary(model)$coefs[2,2]*1.96
    
    line2add <- results[0,]
    line2add[1,1] <- shortax[j]
    line2add[1,2] <- z
    line2add[1,3:5] <- c(ef,lwr,upr)

    results <- rbind(results,line2add)
    
     }
}

#### Figure S3 ####

results$tax <- str_replace(results$tax, 'acti', 'fish')
results <- results %>% arrange(tax)

pdf('~/Desktop/FigS3.pdf',width=5,pointsize = 8,height=6)

layout(rbind(c(1,1),c(2,3),c(4,5)))
par(mar=c(4,4,1,1),oma=c(0,0,0,0),cex=1)

plot(0,type='n',yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l',xlim=c(1,16),ylim=range(c(results$lwr,results$upr)))
title(xlab='minimum number of years in time series', cex.lab=1)
title(ylab='fixed effect of time', cex.lab=1)
axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
abline(h=0,lty=3)
abline(v=c(4.5,8.5,12.5))
axis(1,cex.axis=1,lwd=0,lwd.ticks=0,at=1:16,labels=rep(3:6,4))
axis(1,cex.axis=1,lwd=0,lwd.ticks=1,at=c(4.5,8.5,12.5),labels=rep('',3))
arrows(x0=1:16,y0=results$lwr,y1=results$upr,length=0,lwd=1.5,col=1)
points(x=1:16,y=results$es,pch=20,cex=2,col=1)

for(i in 1:4){
  plot(slopes~n.years_mean,subset(tsdat, tax == taxa[i]),type='n',yaxt='n',xaxt='n',cex.axis=1,ann=F,bty='l')
  axis(2,cex.axis=1,lwd=0,lwd.ticks=1)
  axis(1,cex.axis=1,lwd=0,lwd.ticks=1)
  title(xlab='number of years', cex.lab=1)
  title(ylab='population trend', cex.lab=1)
  abline(h=0,lty=2)
  points(slopes~jitter(n.years_mean),subset(tsdat, tax == taxa[i]),pch=1,cex=1,col=alpha(1,0.2))
  abline(lm(slopes~n.years_mean,subset(tsdat, tax == taxa[i])),lwd=2)
}

dev.off()

### how many time series, species etc.?

alldata_TS <- data.frame()

for(j in 1:4){
  dat <- get(paste0(shortax[j],scales[i],'.agg')) %>%
    filter(.,!is.na(select(.,human.dens.var))) %>% 
    filter(year >= treshold.yr)
  dat %<>% filter(div < mean(div)+10*sd(div))
  dat <- dat %>% select(-n.years) %>% add_count(pop) %>% rename(n.years = n)
  dat <- dat[dat$n.years >= 3,]
  dat <- droplevels(dat)
  alldata_TS <- bind_rows(alldata_TS,dat)
}

n_distinct(alldata_TS$pop)
n_distinct(alldata_TS$species)
n_distinct(alldata_TS$cell)
alldata_TS %>% group_by(pop) %>% summarize(nb.yr = n())
alldata_TS %>% group_by(pop) %>% summarize(nb.yr = n()) -> temp
mean(temp$nb.yr)
sd(temp$nb.yr)/sqrt(nrow(temp))
range(temp$nb.yr)

