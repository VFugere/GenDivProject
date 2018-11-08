#### Spatio-temporal variation in intra-specific neutral genetic diversity across anthromes
#### Gonzalez Lab project - McGill University - 2016-2018
#### Script by Vincent Fugère

rm(list=ls())
options(tibble.print_max = 100, scipen = 999)

#### Part III: model fitting on AWS EC2

library(lme4)
library(tidyverse)
library(magrittr)
library(scales)
library(RColorBrewer)
library(cplm)

#load data
load('~/Google Drive/Recherche/Intraspecific genetic diversity/alldata.Rdata')

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
alldata$lu <- subset(alldata, select = land.use.var)[,1]
alldata$s.yr <- alldata$year
alldata$sc.pop <- log1p(subset(alldata, select =human.dens.var))[,1]
alldata$salat <- abs(alldata$lat)
alldata$snseqs <- log(alldata$nseqs)
alldata$slong <- alldata$long

# model formula
f1 <- formula(div ~ s.yr+salat+slong+lu+sc.pop)
f2 <- update(f1, .~.*.)
f3 <- update(f2, . ~ . + (1|species) + (1|cell))

#mam08
dat <- filter(alldata, scale == '08' & tax == 'mammals') %>% select(div,lu,s.yr,sc.pop,salat,slong,species,cell,nseqs) %>% mutate_at(vars(lu,s.yr,sc.pop,salat,slong),funs(as.numeric(scale(.))))
M1 <- cpglmm(f3, data=dat, weights=log(nseqs))
summary(M1)
M2 <- cpglmm(div ~ lu+s.yr+sc.pop+salat+slong+lu:s.yr+sc.pop:salat+(1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M2)
M3 <- cpglmm(div ~ lu+s.yr+sc.pop+salat+lu:s.yr+sc.pop:salat+(1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M3)
models <- list(M1,M2,M3)

#aves08
dat <- filter(alldata, scale == '08' & tax == 'birds') %>% select(div,lu,s.yr,sc.pop,salat,slong,species,cell,nseqs) %>% mutate_at(vars(lu,s.yr,sc.pop,salat,slong),funs(as.numeric(scale(.))))
M1 <- cpglmm(f3, data=dat, weights=log(nseqs))
summary(M1)
M2 <- cpglmm(div ~ s.yr+salat+slong+lu+sc.pop +(1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M2)
M3 <- cpglmm(div ~ s.yr+lu+(1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M3)
models <- append(models, c(M1,M2,M3))

#fish08
dat <- filter(alldata, scale == '08' & tax == 'fish') %>% select(div,lu,s.yr,sc.pop,salat,slong,species,cell,nseqs) %>% mutate_at(vars(lu,s.yr,sc.pop,salat,slong),funs(as.numeric(scale(.))))
M1 <- cpglmm(f3, data=dat, weights=log(nseqs))
summary(M1)
M2 <- cpglmm(div ~ s.yr+salat+slong+lu+sc.pop + (1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M2)
M3 <- cpglmm(div ~ s.yr + salat + (1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M3)
models <- append(models, c(M1,M2,M3))

#insect08
dat <- filter(alldata, scale == '08' & tax == 'insects') %>% select(div,lu,s.yr,sc.pop,salat,slong,species,cell,nseqs) %>% mutate_at(vars(lu,s.yr,sc.pop,salat,slong),funs(as.numeric(scale(.))))
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
dat <- filter(alldata, scale == '1' & tax == 'mammals') %>% select(div,lu,s.yr,sc.pop,salat,slong,species,cell,nseqs) %>% mutate_at(vars(lu,s.yr,sc.pop,salat,slong),funs(as.numeric(scale(.))))
M1 <- cpglmm(f3, data=dat, weights=log(nseqs))
summary(M1)
M2 <- cpglmm(div ~ s.yr + salat + slong + lu + sc.pop +(1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M2)
M3 <- cpglmm(div ~ salat + (1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M3)
models <- list(M1,M2,M3)

#aves1
dat <- filter(alldata, scale == '1' & tax == 'birds') %>% select(div,lu,s.yr,sc.pop,salat,slong,species,cell,nseqs) %>% mutate_at(vars(lu,s.yr,sc.pop,salat,slong),funs(as.numeric(scale(.))))
M1 <- cpglmm(f3, data=dat, weights=log(nseqs))
summary(M1)
M2 <- cpglmm(div ~ s.yr + salat + slong + lu + sc.pop + lu:sc.pop +(1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M2)
M3 <- cpglmm(div ~ s.yr + lu + sc.pop + lu:sc.pop +(1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M3)
models <- append(models, c(M1,M2,M3))

#fish1
dat <- filter(alldata, scale == '1' & tax == 'fish') %>% select(div,lu,s.yr,sc.pop,salat,slong,species,cell,nseqs) %>% mutate_at(vars(lu,s.yr,sc.pop,salat,slong),funs(as.numeric(scale(.))))
M1 <- cpglmm(f3, data=dat, weights=log(nseqs))
summary(M1)
M2 <- cpglmm(div ~ s.yr+salat+slong+lu+sc.pop + s.yr:sc.pop+slong:sc.pop +(1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M2)
M3 <- cpglmm(div ~ s.yr+salat+slong+sc.pop + s.yr:sc.pop+slong:sc.pop + (1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M3)
models <- append(models, c(M1,M2,M3))

#insect1
dat <- filter(alldata, scale == '1' & tax == 'insects') %>% select(div,lu,s.yr,sc.pop,salat,slong,species,cell,nseqs) %>% mutate_at(vars(lu,s.yr,sc.pop,salat,slong),funs(as.numeric(scale(.))))
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
dat <- filter(alldata, scale == '2' & tax == 'mammals') %>% select(div,lu,s.yr,sc.pop,salat,slong,species,cell,nseqs) %>% mutate_at(vars(lu,s.yr,sc.pop,salat,slong),funs(as.numeric(scale(.))))
M1 <- cpglmm(f3, data=dat, weights=log(nseqs))
summary(M1)
M2 <- cpglmm(div ~ s.yr+salat+slong+lu+sc.pop + salat:sc.pop +(1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M2)
M3 <- cpglmm(div ~ salat + sc.pop + salat:sc.pop +(1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M3)
models <- list(M1,M2,M3)

#aves2
dat <- filter(alldata, scale == '2' & tax == 'birds') %>% select(div,lu,s.yr,sc.pop,salat,slong,species,cell,nseqs) %>% mutate_at(vars(lu,s.yr,sc.pop,salat,slong),funs(as.numeric(scale(.))))
M1 <- cpglmm(f3, data=dat, weights=log(nseqs))
summary(M1)
M2 <- cpglmm(div ~ s.yr+salat+slong+lu+sc.pop +(1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M2)
M3 <- cpglmm(div ~ s.yr+salat +(1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M3)
models <- append(models, c(M1,M2,M3))

#fish2
dat <- filter(alldata, scale == '2' & tax == 'fish') %>% select(div,lu,s.yr,sc.pop,salat,slong,species,cell,nseqs) %>% mutate_at(vars(lu,s.yr,sc.pop,salat,slong),funs(as.numeric(scale(.))))
M1 <- cpglmm(f3, data=dat, weights=log(nseqs))
summary(M1)
M2 <- cpglmm(div ~ s.yr+salat+slong+lu+sc.pop + s.yr:lu+s.yr:sc.pop+slong:sc.pop + (1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M2)
M3 <- cpglmm(div ~ s.yr+salat+lu+sc.pop + s.yr:lu+s.yr:sc.pop+ (1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M3)
models <- append(models, c(M1,M2,M3))

#insect2
dat <- filter(alldata, scale == '2' & tax == 'insects') %>% select(div,lu,s.yr,sc.pop,salat,slong,species,cell,nseqs) %>% mutate_at(vars(lu,s.yr,sc.pop,salat,slong),funs(as.numeric(scale(.))))
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
dat <- filter(alldata, scale == '4' & tax == 'mammals') %>% select(div,lu,s.yr,sc.pop,salat,slong,species,cell,nseqs) %>% mutate_at(vars(lu,s.yr,sc.pop,salat,slong),funs(as.numeric(scale(.))))
M1 <- cpglmm(f3, data=dat, weights=log(nseqs))
summary(M1)
M2 <- cpglmm(div ~ s.yr+salat+slong+lu+sc.pop+(1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M2)
M3 <- cpglmm(div ~ salat+(1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M3)
models <- list(M1,M2,M3)

#aves4
dat <- filter(alldata, scale == '4' & tax == 'birds') %>% select(div,lu,s.yr,sc.pop,salat,slong,species,cell,nseqs) %>% mutate_at(vars(lu,s.yr,sc.pop,salat,slong),funs(as.numeric(scale(.))))
M1 <- cpglmm(f3, data=dat, weights=log(nseqs))
summary(M1)
M2 <- cpglmm(div ~ s.yr+salat+slong+lu+sc.pop + s.yr:sc.pop+slong:lu +(1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M2)
M3 <- cpglmm(div ~ s.yr+slong+lu+sc.pop + s.yr:sc.pop+slong:lu + (1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M3)
models <- append(models, c(M1,M2,M3))

#fish4
dat <- filter(alldata, scale == '4' & tax == 'fish') %>% select(div,lu,s.yr,sc.pop,salat,slong,species,cell,nseqs) %>% mutate_at(vars(lu,s.yr,sc.pop,salat,slong),funs(as.numeric(scale(.))))
M1 <- cpglmm(f3, data=dat, weights=log(nseqs))
summary(M1)
M2 <- cpglmm(div ~ s.yr+salat+slong+lu+sc.pop +s.yr:lu+slong:sc.pop+ (1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M2)
M3 <- cpglmm(div ~ s.yr + salat + (1|species)+(1|cell), data=dat, weights=log(nseqs))
summary(M3)
models <- append(models, c(M1,M2,M3))

#insect4
dat <- filter(alldata, scale == '4' & tax == 'insects') %>% select(div,lu,s.yr,sc.pop,salat,slong,species,cell,nseqs) %>% mutate_at(vars(lu,s.yr,sc.pop,salat,slong),funs(as.numeric(scale(.))))
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

##### analysis across all scales ####

load('~/Google Drive/Recherche/Intraspecific genetic diversity/alldata.Rdata')

alldata$tax %<>% factor(.,levels = c('mammals','birds','fish','insects'))
alldata$scale %<>% as.factor

alldata$lu <- subset(alldata, select = land.use.var)[,1]

#scaling variables. log transform for human dens and nseqs because very skewed distribution
alldata$s.yr <- as.numeric(scale(alldata$year))
alldata$lu <- as.numeric(scale(alldata$lu))
alldata$sc.pop <- as.numeric(scale(log1p(subset(alldata, select = human.dens.var)[,1])))
alldata$salat <- as.numeric(scale(abs(alldata$lat))) #absolute value of latitude
alldata$snseqs <- as.numeric(scale(log(alldata$nseqs)))
alldata$slong <- as.numeric(scale(alldata$long))

dat <- alldata %>% select(div,lu,s.yr,sc.pop,salat,slong,species,cell,nseqs,scale,tax)

# model formula
f1 <- formula(div ~ lu+s.yr+sc.pop+salat+slong)
f2 <- update(f1, .~.*.)
f3 <- update(f2, . ~ .*scale)
f4 <- update(f3, . ~ . + (1|species) + (1|cell))

M1 <- cpglmm(f4, data=alldata, weights=log(nseqs))
AIC(M1); summary(M1)

M1 <- lmer(log1p(div) ~ lu+s.yr+sc.pop+salat+slong+scale + (1|species) + (1|cell), data=alldata, weights=log(nseqs))

models <- append(models, c(M1,M2,M3))

#reading models
load('~/Desktop/CPGLMMs/models08.Rdata')
load('~/Desktop/CPGLMMs/models1.Rdata')
load('~/Desktop/CPGLMMs/models2.Rdata')
load('~/Desktop/CPGLMMs/models4.Rdata')

