## setwd('~/Dropbox/sunflower/')
rm(list=ls())
setwd('analysis/parasiteCommunity')
source('src/initialize.R')
source("src/misc.R")
source("src/predictIntervals.R")
source("src/CIplotting.R")
source("src/diagnostics.R")

load('saved/parMods.RData')

## ************************************************************
## network metrics
## ************************************************************

dd.totalAbund <- expand.grid(TotalAbundance=seq(
                                 from= min(spec.wild$TotalAbundance),
                                 to= max(spec.wild$TotalAbundance),
                                 length=10),
                             Doy=  mean(spec.wild$Doy),
                             Nat350=mean(spec.wild$Nat350),
                             TransectType="SF",
                             HR1000=mean(spec.wild$HR1000),
                             SunflowerCurrent1000=mean(spec.wild$SunflowerCurrent1000)

                             ParasitePresence= 0)
