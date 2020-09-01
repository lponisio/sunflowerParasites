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

## are bee/floral richness and abundance correlated?
## cor.test(by.site$Richness, by.site$TotalAbundance)
## cor.test(by.site$FloralRichness, by.site$FloralAbundance)
## cor.test(by.site$FloralDiv, by.site$FloralAbundance)
## cor.test(by.site$FloralDiv, by.site$FloralRichness)
## SHUCKS!

## *************************************************************
## model selection: bee abundunace
## *************************************************************

## full model
bee.abund.mod <- glmer.nb(TotalAbundance~
                              scale(Doy)*TransectType+
                              scale(I(Doy^2))*TransectType+
                              scale(log(Nat350)) +
                              scale(log(Nat1000)) +
                              scale(log(HR350)) +
                              scale(log(HR1000)) +
                              scale(log(SunflowerCurrent1000)) +
                              scale(log(SunflowerLastYr1000)) +
                              scale(log(SunflowerCurrent350)) +
                              scale(log(SunflowerLastYr350)) +
                              TransectType*scale(SFBloom) +
                              scale(FloralAbundance) +
                              scale(FloralRichness) +
                              scale(FloralDiv) +
                              (1|Site),
                          na.action = "na.fail",
                          data=by.site)

## exclude the different gaussian decays from being included in the
## same model
ms.bee.abund <- dredge(bee.abund.mod,
                       subset =
                !("scale(log(Nat1000))" && "scale(log(Nat350))") &&
                !("scale(log(HR350))" && "scale(log(HR1000))") &&
                !("scale(log(SunflowerCurrent1000))" &&
                  "scale(log(SunflowerCurrent350))") &&
                  !("scale(log(SunflowerLastYr1000))" &&
                  "scale(log(SunflowerLastYr350))") &&
                !("scale(FloralAbundance)" && "scale(FloralRichness)")&&
                !("scale(FloralRichness)" && "scale(FloralDiv)"))

## model average within 2 AICc of the min
ma.bee.abund <- model.avg(ms.bee.abund, subset= delta < 2,
                          revised.var = TRUE)


##richness
bee.rich.mod <- lmer(Richness~
                              scale(Doy)*TransectType+
                              scale(I(Doy^2))*TransectType+
                              scale(log(Nat350)) +
                              scale(log(Nat1000)) +
                              scale(log(HR350)) +
                              scale(log(HR1000)) +
                              scale(log(SunflowerCurrent1000)) +
                              scale(log(SunflowerLastYr1000)) +
                              TransectType*scale(SFBloom) +
                              scale(FloralAbundance) +
                              scale(FloralRichness) +
                              scale(FloralDiv) +
                              (1|Site),
                          na.action = "na.fail",
                          data=by.site)

## exclude the different gaussian decays from being included in the
## same model
ms.bee.rich <- dredge(bee.rich.mod,
                       subset =
                !("scale(log(Nat1000))" && "scale(log(Nat350))") &&
                !("scale(log(HR350))" && "scale(log(HR1000))") &&
                !("scale(FloralAbundance)" && "scale(FloralRichness)")&&
                !("scale(FloralRichness)" && "scale(FloralDiv)"))

ma.bee.rich <- model.avg(ms.bee.rich, subset= delta < 2,
                          revised.var = TRUE)


save(ma.bee.abund,
     ms.bee.abund,
     ma.bee.rich,
     ms.bee.rich,
     file=sprintf("saved/%s_beeMods.RData",
                  gsub(" ", "", focal.bee)))
