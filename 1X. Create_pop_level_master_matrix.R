## Vincent Fugere March 22 2019

# This code combines and format the 3 independent pop-level 
# matrices (π, D, and LUHD). The output of this script is 
# the dataframe used in all statistical models, in which
# populations are the unit of replication

rm(list=ls())
library(tidyverse)

#load files
dir <- ('/Users/vincentfugere/Desktop/Data/')
load(file.path(dir,'DF_Pi.RData'))
load(file.path(dir,'DF_D.RData'))
load(file.path(dir,'DF_LUHD.RData'))

alldata <- inner_join(DF_Pi,DF_D) %>% inner_join(DF_LUHD)