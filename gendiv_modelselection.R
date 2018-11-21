#### Spatio-temporal variation in intra-specific neutral genetic diversity across anthromes
#### Gonzalez Lab project - McGill University - 2016-2018
#### Script by Vincent Fugère

rm(list=ls())
options(tibble.print_max = 100, scipen = 999)

#### Part III: model fitting for global model using all data (not only time series)

library(lme4)
library(tidyverse)
library(magrittr)
library(scales)
library(RColorBrewer)
library(cplm)

## notes 

# Cannot add random slopes as most species have a single data point.
# Adding (1|species/pop) makes it impossible to fit a model for insects. 
# Therefore, removed duplicated rows (sample 2+ for pops with repeated sampling) and kept a single data point per pop
# Bayesian models never fit as of Dec 2017

#load data
load('~/Google Drive/Recherche/Intraspecific genetic diversity/modeldata.Rdata')

#parameters
treshold.yr <- 1980
types.lu <- c('cropland','pasture','urban','conv_rangeland') #to calculate total land use (p.lu) and land use diversity
land.use.var <- 'p.lu' #which land use variable should be included in model (div, or total, or chg - they are all collinear)
human.dens.var <- 'pop_tot' #name of human population density variable to use
cols <- brewer.pal(4,'Dark2') 
taxa <- c('mammals','birds','fish','insects')
scales <- c('08','1','2','4')
shortax <- c('mam','aves','acti','insect')

#### model selection ####

#adding/transforming variables - scaling is further below to scale within sub-dataset
modeldata$lu <- subset(modeldata, select = land.use.var)[,1]
modeldata$s.yr <- modeldata$year
modeldata$sc.pop <- log1p(subset(modeldata, select =human.dens.var))[,1]
modeldata$salat <- abs(modeldata$lat)
modeldata$snseqs <- log(modeldata$nseqs)
modeldata$slong <- modeldata$long

# model formula
f1 <- formula(div ~ s.yr+salat+slong+lu+sc.pop)
f2 <- update(f1, .~.*.)
f3 <- update(f2, . ~ . + (1|species) + (1|cell))

#mam08
dat <- filter(modeldata, scale == '08' & tax == 'mammals') %>% select(div,lu,s.yr,sc.pop,salat,slong,species,cell,nseqs) %>% mutate_at(vars(lu,s.yr,sc.pop,salat,slong),funs(as.numeric(scale(.))))
M1 <- cpglmm(f3, data=dat, weights=log(nseqs))
summary(M1)
M2 <- cpglmm(div ~ s.yr+salat+slong+lu+sc.pop+s.yr:lu+salat:sc.pop+(1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M2)
M3 <- cpglmm(div ~ s.yr+salat+lu+sc.pop+s.yr:lu+salat:sc.pop+(1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M3)
models <- list(M1,M2,M3)

#aves08
dat <- filter(modeldata, scale == '08' & tax == 'birds') %>% select(div,lu,s.yr,sc.pop,salat,slong,species,cell,nseqs) %>% mutate_at(vars(lu,s.yr,sc.pop,salat,slong),funs(as.numeric(scale(.))))
M1 <- cpglmm(f3, data=dat, weights=log(nseqs))
summary(M1)
M2 <- cpglmm(div ~ s.yr+salat+slong+lu+sc.pop +(1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M2)
M3 <- cpglmm(div ~ s.yr+lu+(1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M3)
models <- append(models, c(M1,M2,M3))

#fish08
dat <- filter(modeldata, scale == '08' & tax == 'fish') %>% select(div,lu,s.yr,sc.pop,salat,slong,species,cell,nseqs) %>% mutate_at(vars(lu,s.yr,sc.pop,salat,slong),funs(as.numeric(scale(.))))
M1 <- cpglmm(f3, data=dat, weights=log(nseqs))
summary(M1)
M2 <- cpglmm(div ~ s.yr+salat+slong+lu+sc.pop + (1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M2)
M3 <- cpglmm(div ~ s.yr + salat + (1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M3)
models <- append(models, c(M1,M2,M3))

#insect08
dat <- filter(modeldata, scale == '08' & tax == 'insects') %>% select(div,lu,s.yr,sc.pop,salat,slong,species,cell,nseqs) %>% mutate_at(vars(lu,s.yr,sc.pop,salat,slong),funs(as.numeric(scale(.))))
M1 <- cpglmm(f3, data=dat, weights=log(nseqs))
summary(M1)
M2 <- cpglmm(div ~ s.yr + salat + slong + lu + sc.pop +s.yr:salat+salat:slong+salat:lu+slong:lu +(1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M2)
M3 <- cpglmm(div ~ s.yr + salat + slong + lu + sc.pop + s.yr:salat+salat:slong + (1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M3)
M4 <- cpglmm(div ~ s.yr + salat + slong + s.yr:salat+salat:slong + (1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M4)
models <- append(models, c(M1,M2,M3,M4))

models08 <- models
save(models08, file='~/Desktop/models08.Rdata')
rm(models, models08, M1, M2, M3, M4)

#mam1
dat <- filter(modeldata, scale == '1' & tax == 'mammals') %>% select(div,lu,s.yr,sc.pop,salat,slong,species,cell,nseqs) %>% mutate_at(vars(lu,s.yr,sc.pop,salat,slong),funs(as.numeric(scale(.))))
M1 <- cpglmm(f3, data=dat, weights=log(nseqs))
summary(M1)
M2 <- cpglmm(div ~ s.yr + salat + slong + lu + sc.pop +(1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M2)
M3 <- cpglmm(div ~ salat + (1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M3)
models <- list(M1,M2,M3)

#aves1
dat <- filter(modeldata, scale == '1' & tax == 'birds') %>% select(div,lu,s.yr,sc.pop,salat,slong,species,cell,nseqs) %>% mutate_at(vars(lu,s.yr,sc.pop,salat,slong),funs(as.numeric(scale(.))))
M1 <- cpglmm(f3, data=dat, weights=log(nseqs))
summary(M1)
M2 <- cpglmm(div ~ s.yr + salat + slong + lu + sc.pop + lu:sc.pop +(1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M2)
M3 <- cpglmm(div ~ s.yr + lu + sc.pop + lu:sc.pop +(1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M3)
models <- append(models, c(M1,M2,M3))

#fish1
dat <- filter(modeldata, scale == '1' & tax == 'fish') %>% select(div,lu,s.yr,sc.pop,salat,slong,species,cell,nseqs) %>% mutate_at(vars(lu,s.yr,sc.pop,salat,slong),funs(as.numeric(scale(.))))
M1 <- cpglmm(f3, data=dat, weights=log(nseqs))
summary(M1)
M2 <- cpglmm(div ~ s.yr+salat+slong+lu+sc.pop + s.yr:sc.pop+slong:sc.pop +(1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M2)
M3 <- cpglmm(div ~ s.yr+salat+slong+sc.pop + s.yr:sc.pop+slong:sc.pop + (1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M3)
models <- append(models, c(M1,M2,M3))

#insect1
dat <- filter(modeldata, scale == '1' & tax == 'insects') %>% select(div,lu,s.yr,sc.pop,salat,slong,species,cell,nseqs) %>% mutate_at(vars(lu,s.yr,sc.pop,salat,slong),funs(as.numeric(scale(.))))
M1 <- cpglmm(f3, data=dat, weights=log(nseqs))
summary(M1)
M2 <- cpglmm(div ~ s.yr + salat + slong + lu + sc.pop + s.yr:salat +(1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M2)
M3 <- cpglmm(div ~ s.yr + salat + slong + lu + sc.pop + (1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M3)
M4 <- cpglmm(div ~ salat + (1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M4)
models <- append(models, c(M1,M2,M3,M4))

models1 <- models
save(models1, file='~/Desktop/models1.Rdata')
rm(models, models1, M1, M2, M3, M4)

#mam2
dat <- filter(modeldata, scale == '2' & tax == 'mammals') %>% select(div,lu,s.yr,sc.pop,salat,slong,species,cell,nseqs) %>% mutate_at(vars(lu,s.yr,sc.pop,salat,slong),funs(as.numeric(scale(.))))
M1 <- cpglmm(f3, data=dat, weights=log(nseqs))
summary(M1)
M2 <- cpglmm(div ~ s.yr+salat+slong+lu+sc.pop + salat:sc.pop +(1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M2)
M3 <- cpglmm(div ~ salat + sc.pop + salat:sc.pop +(1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M3)
models <- list(M1,M2,M3)

#aves2
dat <- filter(modeldata, scale == '2' & tax == 'birds') %>% select(div,lu,s.yr,sc.pop,salat,slong,species,cell,nseqs) %>% mutate_at(vars(lu,s.yr,sc.pop,salat,slong),funs(as.numeric(scale(.))))
M1 <- cpglmm(f3, data=dat, weights=log(nseqs))
summary(M1)
M2 <- cpglmm(div ~ s.yr+salat+slong+lu+sc.pop +(1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M2)
M3 <- cpglmm(div ~ s.yr+salat +(1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M3)
models <- append(models, c(M1,M2,M3))

#fish2
dat <- filter(modeldata, scale == '2' & tax == 'fish') %>% select(div,lu,s.yr,sc.pop,salat,slong,species,cell,nseqs) %>% mutate_at(vars(lu,s.yr,sc.pop,salat,slong),funs(as.numeric(scale(.))))
M1 <- cpglmm(f3, data=dat, weights=log(nseqs))
summary(M1)
M2 <- cpglmm(div ~ s.yr+salat+slong+lu+sc.pop + s.yr:lu+s.yr:sc.pop+slong:sc.pop + (1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M2)
M3 <- cpglmm(div ~ s.yr+salat+lu+sc.pop + s.yr:lu+s.yr:sc.pop+ (1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M3)
models <- append(models, c(M1,M2,M3))

#insect2
dat <- filter(modeldata, scale == '2' & tax == 'insects') %>% select(div,lu,s.yr,sc.pop,salat,slong,species,cell,nseqs) %>% mutate_at(vars(lu,s.yr,sc.pop,salat,slong),funs(as.numeric(scale(.))))
M1 <- cpglmm(f3, data=dat, weights=log(nseqs))
summary(M1)
M2 <- cpglmm(div ~ s.yr + salat + slong + lu + sc.pop +s.yr:salat+s.yr:sc.pop+salat:lu+lu:sc.pop +(1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M2)
M3 <- cpglmm(div ~ s.yr + salat + lu + sc.pop +s.yr:salat+s.yr:sc.pop+salat:lu+lu:sc.pop + (1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M3)
models <- append(models, c(M1,M2,M3))

models2 <- models
save(models2, file='~/Desktop/models2.Rdata')
rm(models, models2, M1, M2, M3)

#mam4
dat <- filter(modeldata, scale == '4' & tax == 'mammals') %>% select(div,lu,s.yr,sc.pop,salat,slong,species,cell,nseqs) %>% mutate_at(vars(lu,s.yr,sc.pop,salat,slong),funs(as.numeric(scale(.))))
M1 <- cpglmm(f3, data=dat, weights=log(nseqs))
summary(M1)
M2 <- cpglmm(div ~ s.yr+salat+slong+lu+sc.pop+(1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M2)
M3 <- cpglmm(div ~ salat+(1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M3)
models <- list(M1,M2,M3)

#aves4
dat <- filter(modeldata, scale == '4' & tax == 'birds') %>% select(div,lu,s.yr,sc.pop,salat,slong,species,cell,nseqs) %>% mutate_at(vars(lu,s.yr,sc.pop,salat,slong),funs(as.numeric(scale(.))))
M1 <- cpglmm(f3, data=dat, weights=log(nseqs))
summary(M1)
M2 <- cpglmm(div ~ s.yr+salat+slong+lu+sc.pop + s.yr:sc.pop+slong:lu +(1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M2)
M3 <- cpglmm(div ~ s.yr+slong+lu+sc.pop + s.yr:sc.pop+slong:lu + (1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M3)
models <- append(models, c(M1,M2,M3))

#fish4
dat <- filter(modeldata, scale == '4' & tax == 'fish') %>% select(div,lu,s.yr,sc.pop,salat,slong,species,cell,nseqs) %>% mutate_at(vars(lu,s.yr,sc.pop,salat,slong),funs(as.numeric(scale(.))))
M1 <- cpglmm(f3, data=dat, weights=log(nseqs))
summary(M1)
M2 <- cpglmm(div ~ s.yr+salat+slong+lu+sc.pop +s.yr:lu+slong:sc.pop+ (1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M2)
M3 <- cpglmm(div ~ s.yr + salat + (1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M3)
models <- append(models, c(M1,M2,M3))

#insect4
dat <- filter(modeldata, scale == '4' & tax == 'insects') %>% select(div,lu,s.yr,sc.pop,salat,slong,species,cell,nseqs) %>% mutate_at(vars(lu,s.yr,sc.pop,salat,slong),funs(as.numeric(scale(.))))
M1 <- cpglmm(f3, data=dat, weights=log(nseqs))
summary(M1)
M2 <- cpglmm(div ~ s.yr + salat + slong + lu + sc.pop +s.yr:salat+lu:sc.pop +(1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M2)
M3 <- cpglmm(div ~ s.yr + salat + lu + sc.pop +s.yr:salat+lu:sc.pop + (1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M3)
models <- append(models, c(M1,M2,M3))

models4 <- models
save(models4, file='~/Desktop/models4.Rdata')
rm(models, models4, M1, M2, M3)

# ## Model diagnostic code. 
# ## Need a sub dataset with one scale+taxon, cleaned up already (as when grabbed from modeldata)
# ## and a cpglmm object called 'model'
# 
# #pdf to store plots
# pdf(file = paste0('~/Google Drive/Recherche/Intraspecific genetic diversity/diagplots_',land.use.var,'_',file.tag,'.pdf'),width = 8, height = 8, onefile=T)
# par(mfrow=c(3,3),oma=c(0,0,3,0),mar=c(4,4,1,1))
#
# #checking colinearity
# corvif(dat[,c('s.yr','snseqs','lu','salat','slong','sc.pop')])
#
# #validation
# dat$R <- resid(model, type='standardized')
# dat$Yfit <- fitted(model)
# 
# # #how good is the model fit?
# # plot(Yfit~div,dat,ylab='fitted',xlab='observed',pch=16,col=alpha(1,0.2))
# # abline(a=0,b=1,lty=2)
# # abline(lm(Yfit~div,dat))
# # round(summary(lm(Yfit~div,dat))$adj.r.squared,4) -> r2
# # legend('topleft',legend=paste0('(',r2,')'),bty='n')
# #
# #heteroscedasticity?
# plot(R~Yfit,dat,ylab='standardized residuals',xlab='fitted',pch=16,col=alpha(1,0.2))
# 
# #nuisance pattern in residuals as function of predictors?
# scatter.smooth(x=dat$s.yr,y=dat$R,xlab='year',ylab='standardized residuals',pch=16,col=alpha(1,0.2))
# scatter.smooth(x=dat$salat,y=dat$R,xlab='latitude',ylab='standardized residuals',pch=16,col=alpha(1,0.2))
# scatter.smooth(x=as.numeric(scale(dat$long)),y=dat$R,xlab='longitude',ylab='standardized residuals',pch=16,col=alpha(1,0.2))
# scatter.smooth(x=dat$lu,y=dat$R,xlab='land use',ylab='standardized residuals',pch=16,col=alpha(1,0.2))
# scatter.smooth(x=dat$sc.pop,y=dat$R,xlab='human density',ylab='standardized residuals',pch=16,col=alpha(1,0.2))
# # scatter.smooth(x=as.numeric(scale(log(dat$nseqs))),y=dat$R,xlab='number of sequences',ylab='standardized residuals',pch=16,col=alpha(1,0.2))
# # boxplot(R~cell,dat,names=NULL,xlab='cell')
# # abline(h=0,lty=3)
# # boxplot(R~species,dat,names=NULL,xlab='species')
# # abline(h=0,lty=3)
# # boxplot(R~pop,dat,names=NULL,xlab='population')
# # abline(h=0,lty=3)
# # scatter.smooth(x=scale(dat$lu.chg),y=dat$R,xlab='land use chg',ylab='standardized residuals')
# # scatter.smooth(x=scale(dat$pop.chg),y=dat$R,xlab='human dens chg',ylab='standardized residuals')
# 
# #temporal autocorrelation in residuals?
# dat %>% group_by(year) %>% summarize('mean' = mean(R)) -> resid.ts
# acf(resid.ts$mean,main='')
# 
# #spatial autocorrelation in residuals?
# dat$rcol <- 4
# dat[dat$R < 0, 'rcol'] <- 2
# plot(map, xlim = c(-180,180), ylim = c(-90,90),border=NA,col='light gray',axes=F)
# points(lat~long,dat,pch=16,col=alpha(dat$rcol,0.5),cex=0.5+abs(dat$R))
# spatdat <- select(dat, long,lat,R)
# coordinates(spatdat) <- c('long','lat')
# variogram(R~long+lat,spatdat,width=0.1,cutoff=50) -> var1
# scatter.smooth(x=var1$dist,y=var1$gamma,xlab='distance',ylab='semivariance',ylim=c(0,max(var1$gamma)),pch=16,col='gray')
