## setwd("~/Dropbox/sunflower")
setwd("analysis/parasiteCommunity")
rm(list=ls())
library(lme4)
library(lmerTest)
library(car)
library(MuMIn)

focal.bee <- "all"

source("src/initialize.R")

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

## *************************************************************
## parasite richness honey bees
## *************************************************************

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

save(ma.parasite.rich.hb,
     ms.parasite.rich.hb,
     ma.parasite.pres.hb,
     ms.parasite.pres.hb,
     file="saved/HBMods.RData")

## plotting
source("src/misc.R")
source("src/predictIntervals.R")
source("src/plotPanels.R")
source("src/CIplotting.R")
source("src/diagnostics.R")
library(viridis)
library(boot)
library(pals)

## Total abundance
summary(ma.parasite.pres.hb)
top.mod.hb <- get.models(ms.parasite.pres.hb, 1)[[1]]
summary(top.mod.hb)
r.squaredGLMM(top.mod.hb)


## dd.hb.rich <- expand.grid(FloralRichness=seq(
##                             from= min(by.site$FloralRichness, na.rm=TRUE),
##                             to= max(by.site$FloralRichness, na.rm=TRUE),
##                             length=20),
##                         ParasitePresence=0)

## hb.rich.pi <- predict.int(top.mod.hb,
##                         dd.hb.rich,
##                         "ParasitePresence",
##                         "binomial")

## by.site2 <- by.site
## colnames(by.site2)[colnames(by.site2) == "ParasitePresence"] <-
##     "ParasitePresence1"

## colnames(by.site2)[colnames(by.site2) == "HBParasitism"] <-
##     "ParasitePresence"

## cols.var <- add.alpha("goldenrod3", 0.5)
## names(cols.var) <- "all"

## ## sig only modelsx
## pdf.f(plotSigModelsHBPres,
##       file="figures/HBParasitism.pdf",
##       height=5, width=6)
