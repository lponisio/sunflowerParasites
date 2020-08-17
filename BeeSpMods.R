## setwd("~/Dropbox/sunflower")
setwd("analysis/parasiteCommunity")
rm(list=ls())
library(lme4)
library(lmerTest)
library(car)
library(MuMIn)

## focal bee set by deafult to all if not specified
source("src/initialize.R")
print(focal.bee)

## *************************************************************
## model selection: bee abundunace
## *************************************************************

## full model
bee.abund.mod <- glmer(TotalAbundance~
                              scale(Doy)+
                              scale(log(Nat350))*SFspecialist +
                              scale(log(Nat1000))*SFspecialist +
                              scale(log(HR350))*SFspecialist +
                              scale(log(HR1000))*SFspecialist +
                              scale(log(SunflowerCurrent1000))*SFspecialist +
                              scale(log(SunflowerLastYr1000))*SFspecialist +
                              TransectType*scale(SFBloom) +
                              scale(FloralAbundance) +
                              scale(FloralRichness) +
                              scale(FloralDiv) +
                              (1|Site),
                          na.action = "na.fail",
                       data=sp.by.site.wild,
                       family="poisson")
## poisson has much lower aic than negative binominal

## exclude the different gaussian decays from being included in the
## same model
ms.bee.abund <- dredge(bee.abund.mod,
                       subset =
                !("scale(log(Nat1000))" && "scale(log(Nat350))") &&
                !("scale(log(HR350))" && "scale(log(HR1000))") &&
                !("scale(FloralAbundance)" && "scale(FloralRichness)")&&
                !("scale(FloralRichness)" && "scale(FloralDiv)"))

ma.bee.abund <- model.avg(ms.bee.abund, subset= delta < 2,
                          revised.var = TRUE)

summary(ma.bee.abund)

## print(ms.bee.abund, abbrev.names=FALSE)


save(ma.bee.abund,
     ms.bee.abund,
     file=sprintf("saved/%s_SpbeeMods.RData",
                  gsub(" ", "", focal.bee)))
