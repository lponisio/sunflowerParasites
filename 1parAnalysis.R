## setwd("~/Dropbox/sunflower")
setwd("analysis/parasiteCommunity")
rm(list=ls())
library(lme4)
library(lmerTest)
library(car)
library(MuMIn)

source("src/initialize.R")


## are bee/floral richness and abundance correlated?
cor.test(by.site$Richness, by.site$TotalAbundance)
cor.test(by.site$FloralRichness, by.site$FloralAbundance)
cor.test(by.site$FloralDiv, by.site$FloralAbundance)
cor.test(by.site$FloralDiv, by.site$FloralRichness)
## SHUCKS!

## *************************************************************
## model selection: bee abundunace
## *************************************************************

## full model
bee.abund.mod <- glmer.nb(TotalAbundance~
                              scale(Doy)+
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

## *************************************************************
## bee richness
## *************************************************************

bee.rich.mod <- lmer(Richness~ scale(Doy) +
                         scale(log(Nat1000)) +
                         scale(log(Nat350)) +
                         scale(log(HR350)) +
                         scale(log(HR1000)) +
                         scale(log(SunflowerCurrent1000)) +
                         TransectType*scale(SFBloom) +
                         scale(FloralAbundance) +
                         scale(FloralDiv) +
                         scale(FloralRichness) +
                         (1|Site),
                     na.action = "na.fail",
                     data=by.site)

## exclude the different gaussian decays from being included in the
## same model
ms.bee.rich <- dredge(bee.rich.mod,
       subset =  !("scale(log(Nat1000))" && "scale(log(Nat350))") &&
       !("scale(log(HR350))" && "scale(log(HR1000))") &&
       !("scale(FloralAbundance)" && "scale(FloralRichness)")&&
       !("scale(FloralRichness)" && "scale(FloralDiv)"))

ma.bee.rich <- model.avg(ms.bee.rich, subset= delta < 2,
                         revised.var = TRUE)

## just use the top model since no other models are withing 2 AIC
ma.bee.rich <- get.models(ms.bee.rich, 1)[[1]]

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
                               scale(FloralAbundance) +
                               scale(FloralDiv) +
                               scale(FloralRichness) +
                               (1|Site),
                           family="binomial",
                           glmerControl(optimizer="bobyqa"),
                           data=spec.wild.sub,
                           na.action = "na.fail")

## include richness and abundaunce from the same model as they are
## very colinear
ms.parasite.pres <- dredge(parasite.pres.mod,
            subset =  !("scale(Richness)" && "scale(TotalAbundance)") &&
            !("scale(FloralAbundance)" && "scale(FloralRichness)")&&
            !("scale(FloralRichness)" && "scale(FloralDiv)"))

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
                               scale(FloralAbundance) +
                               scale(FloralDiv) +
                               scale(FloralRichness) +
                               (1|Site),
                           family="binomial",
                           glmerControl(optimizer="bobyqa"),
                           data=spec.wild.sub,
                           na.action = "na.fail")

## include richness and abundaunce from the same model as they are
## very colinear
ms.parasite.rich <- dredge(parasite.rich.mod,
         subset =  !("scale(Richness)" && "scale(TotalAbundance)") &&
            !("scale(FloralAbundance)" && "scale(FloralRichness)")&&
            !("scale(FloralRichness)" && "scale(FloralDiv)"))


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

parasite.pres.mod.hb <- glmer(ParasitePresence~
                                  TransectType*scale(SFBloom) +
                                  scale(TotalAbundance) +
                                  scale(Richness) +
                                  scale(FloralAbundance) +
                                  scale(FloralRichness) +
                                  scale(FloralDiv) +
                                  (1|Site),
                              family="binomial",
                              glmerControl(optimizer="bobyqa"),
                              data=hb,
                              na.action = "na.fail")

ms.parasite.pres.hb <- dredge(parasite.pres.mod.hb,
       subset =  !("scale(Richness)" && "scale(TotalAbundance)")&&
            !("scale(FloralRichness)" && "scale(FloralDiv)"))

ma.parasite.pres.hb <- model.avg(ms.parasite.pres.hb,
                                 subset= delta < 2,
                                 revised.var = TRUE)

summary(ma.parasite.pres.hb)

parasite.rich.mod.hb <- glmer(cbind(ParasiteRichness,
                                    PossibleParasite)~
                                  TransectType*scale(SFBloom) +
                                  scale(TotalAbundance) +
                                  scale(Richness) +
                                  scale(FloralAbundance) +
                                  scale(FloralDiv) +
                                  scale(FloralRichness) +
                                  (1|Site),
                              family="binomial",
                              glmerControl(optimizer="bobyqa"),
                              data=hb,
                              na.action = "na.fail")

ms.parasite.rich.hb <- dredge(parasite.rich.mod.hb,
         subset =  !("scale(Richness)" && "scale(TotalAbundance)")&&
            !("scale(FloralRichness)" && "scale(FloralDiv)"))

ma.parasite.rich.hb <- model.avg(ms.parasite.rich.hb,
                                 subset= delta < 2,
                                 revised.var = TRUE)

summary(ma.parasite.rich.hb)

save(ma.bee.abund,
     ms.bee.abund,
     ma.bee.rich,
     ms.bee.rich,
     ma.parasite.pres,
     ms.parasite.pres,
     ma.parasite.rich,
     ms.parasite.rich,
     ma.parasite.rich.hb,
     ms.parasite.rich.hb,
     ma.parasite.pres.hb,
     ms.parasite.pres.hb,
     file="saved/parMods.RData")


