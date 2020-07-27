## setwd("~/Dropbox/sunflower")
setwd("analysis/parasiteCommunity")
rm(list=ls())
library(lme4)
library(lmerTest)
library(car)
library(MuMIn)

source("src/initialize.R")


## are bee richness and abundance correlated?
cor.test(by.site$Richness, by.site$TotalAbundance)
## SHUCKS!

## *************************************************************
## model selection: bee abundunace
## *************************************************************

## full model
bee.abund.mod <- glmer.nb(TotalAbundance~
                              scale(Doy)+
                              scale(Nat350) +
                              scale(Nat1000) +
                              scale(HR350) +
                              scale(HR1000) +
                              scale(SunflowerCurrent1000) +
                              scale(SunflowerLastYr1000) +
                              TransectType*scale(SFBloom) +
                              (1|Site),
                          na.action = "na.fail",
                          data=by.site)
## exclude the different gaussian decays from being included in the
## same model
ms.bee.abund <- dredge(bee.abund.mod,
                       subset =  !("scale(Nat1000)" && "scale(Nat350)") &&
                           !("scale(HR350)" && "scale(HR1000)"))

subset(ms.bee.abund, delta <2)
ma.bee.abund <- model.avg(ms.bee.abund, subset= delta < 2,
                          revised.var = TRUE)

summary(ma.bee.abund)
## print(ms.bee.abund, abbrev.names=FALSE)

## *************************************************************
## bee richness
## *************************************************************

bee.rich.mod <- lmer(Richness~ scale(Doy) +
                         scale(Nat1000) +
                           scale(Nat350) +
                           scale(HR350) +
                           scale(HR1000) +
                           scale(SunflowerCurrent1000) +
                           scale(SunflowerLastYr1000) +
                           TransectType*scale(SFBloom) +
                           (1|Site),
                        na.action = "na.fail",
                       data=by.site)
## exclude the different gaussian decays from being included in the
## same model
ms.bee.rich <- dredge(bee.rich.mod,
   subset =  !("scale(Nat1000)" && "scale(Nat350)") &&
            !("scale(HR350)" && "scale(HR1000)"))

subset(ms.bee.rich, delta <2)
ma.bee.rich <- model.avg(ms.bee.rich, subset= delta < 2.2,
                          revised.var = TRUE)

summary(ma.bee.rich)

## *************************************************************
## parasite presence (any parasite)
## *************************************************************
parasite.pres.mod <- glmer(ParasitePresence~
                           TransectType*scale(SFBloom) +
                           scale(TotalAbundance) +
                            scale(Richness) +
                           scale(r.degree) +
                           scale(MeanITD)+
                           Sociality +
                           (1|Site),
                           family="binomial",
                           glmerControl(optimizer="bobyqa"),
                           data=spec.wild.sub,
                           na.action = "na.fail")

## include richness and abundaunce from the same model as they are
## very colinear
ms.parasite.pres <- dredge(parasite.pres.mod,
           subset =  !("scale(Richness)" && "scale(TotalAbundance)"))

subset(ms.parasite.pres, delta <2)
ma.parasite.pres <- model.avg(ms.parasite.pres, subset= delta < 2,
                          revised.var = TRUE)

summary(ma.parasite.pres)

## simplified model sociality, degree and body size are highly
## colinear. Small bees tend to be social and generalized


## the rate of parasitism in generalist bees < specialist bees
## large bees > small bees
## social bees > solitary
## higher abundaunce lowers parasitism rates, suggesting dilution

## *************************************************************
## parasite richness within a bee
## *************************************************************

parasite.rich.mod <- glmer(cbind(ParasiteRichness, PossibleParasite)~
                           TransectType*scale(SFBloom) +
                               scale(TotalAbundance) +
                               scale(Richness) +
                           scale(r.degree) +
                           scale(MeanITD)+
                           Sociality +
                           (1|Site),
                           family="binomial",
                           glmerControl(optimizer="bobyqa"),
                           data=spec.wild.sub,
                           na.action = "na.fail")

## include richness and abundaunce from the same model as they are
## very colinear
ms.parasite.rich <- dredge(parasite.rich.mod,
           subset =  !("scale(Richness)" && "scale(TotalAbundance)"))

subset(ms.parasite.rich, delta <2)
ma.parasite.rich <- model.avg(ms.parasite.rich, subset= delta < 2,
                          revised.var = TRUE)

summary(ma.parasite.rich)
## same issue of colinearily between degree and body size

## the number of  parasites  in generalist bees < specialist bees
## large bees > small bees
## higher abundaunce lowers parasite richness, suggesting dilution



## *************************************************************
## parasite presence honey bees
## *************************************************************

hb <- spec[spec$GenusSpecies == "Apis mellifera",]
hb$SFBloom <- as.numeric(hb$SFBloom)

parasite.pres.mod.hb <- glmer(ParasitePresence~
                           TransectType*scale(SFBloom) +
                           scale(TotalAbundance) +
                            scale(Richness) +
                           (1|Site),
                           family="binomial",
                           glmerControl(optimizer="bobyqa"),
                           data=hb,
                           na.action = "na.fail")

ms.parasite.pres.hb <- dredge(parasite.pres.mod.hb,
           subset =  !("scale(Richness)" && "scale(TotalAbundance)"))

subset(ms.parasite.pres.hb, delta <2)
ma.parasite.pres.hb <- model.avg(ms.parasite.pres.hb, subset= delta < 2,
                          revised.var = TRUE)

summary(ma.parasite.pres.hb)

## native bee abundaunce has a positive effect on honey bee parasitism
## rates

parasite.rich.mod.hb <- glmer(cbind(ParasiteRichness, PossibleParasite)~
                           TransectType*scale(SFBloom) +
                           scale(TotalAbundance) +
                            scale(Richness) +
                           (1|Site),
                           family="binomial",
                           glmerControl(optimizer="bobyqa"),
                           data=hb,
                           na.action = "na.fail")

ms.parasite.rich.hb <- dredge(parasite.rich.mod.hb,
           subset =  !("scale(Richness)" && "scale(TotalAbundance)"))

subset(ms.parasite.rich.hb, delta <2)
ma.parasite.rich.hb <- model.avg(ms.parasite.rich.hb, subset= delta < 2,
                          revised.var = TRUE)

summary(ma.parasite.rich.hb)

## parasite richness is - related to native bee richness

save(ma.bee.abund, ma.bee.rich, ma.parasite.pres,
ma.parasite.rich.hb, ma.parasite.pres.hb,
    file="saved/parMods.RData")



coeffs <- summary(ma.parasite.pres)$coefmat.subset
ci <- confint(ma.parasite.pres)

plot(by.site$Parasitism ~ by.site$TotalAbundance)
