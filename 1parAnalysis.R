## setwd("~/Dropbox/sunflower")
setwd("analysis/parasiteCommunity")
rm(list=ls())
## library(devtools)
## install_github("jslefche/piecewiseSEM")
library(piecewiseSEM)
library(lme4)
library(lmerTest)
library(car)
library(MuMIn)

load("../../data/spec.Rdata")
by.site <- read.csv("../../data/bySite.csv")
load('../../../sunflower_saved/data/nat_HR_buffers.Rdata')

source("src/initialize.R")

## are bee richness and abundance correlated?
cor.test(by.site$Richness, by.site$TotalAbundance)
## SHUCKS!

## *************************************************************
## make formulas for path analyses
## *************************************************************

spec.wild.sub <- spec.wild[!is.na(spec.wild$MeanITD),]

## formula for site effects on the bee community
formula.bee <- formula(s.TotalAbundance ~ s.Nat1000 +
                           s.HR1000 +
                           AdjHR +
                           s.SFBloom   +
                           (1|Site))

## formulas for the site effects on parasites

formula.par <- formula(ParasitePresence ~
                           AdjHR +
                           s.SFBloom   +
                           s.TotalAbundance +
                           s.r.degree +
                           (1|Site))


## don't know why this doesn't work!!!!!!!!!!!!!!
mods  <- psem(
        BeeAbund = lmer(formula.bee,
                         data = by.site),
        ParPath = glmer(formula.par,
                       data = spec.wild.sub,
                       family="binomial"))

summary(mods)
## *************************************************************
## trying withing sem in pieces
## *************************************************************

library(car)
bee.abund.mod <- glmer(TotalAbundance~ scale(Nat1000) +
                           scale(HR1000) +
                           TransectType*scale(SFBloom) +
                           AdjHR +
                           (1|Site),
                       data=by.site, family="poisson")

vif(bee.abund.mod)
summary(bee.abund.mod)
## plot(density(residuals(bee.abund.mod)))

bee.rich.mod <- lmer(Richness~ scale(Nat1000) +
                           scale(HR350) +
                           TransectType*scale(SFBloom) +
                           AdjHR +
                           (1|Site),
                       data=by.site)

vif(bee.rich.mod)
summary(bee.rich.mod)
## plot(density(residuals(bee.rich.mod)))

## parasite presence (any parasite)
bee.parasite.pres.mod <- glmer(ParasitePresence~
                           AdjHR +
                           TransectType +scale(SFBloom) +
                           scale(TotalAbundance) +
                           scale(r.degree) +
                           scale(MeanITD)+
                           Sociality +
                           NestLoc +
                           (1|Site),
                           family="binomial",
                           glmerControl(optimizer="bobyqa"),
                           data=spec.wild.sub,
                           na.action = "na.fail")

vif(bee.parasite.pres.mod)
summary(bee.parasite.pres.mod)

ms.parasite.pres <- dredge(bee.parasite.pres.mod)
print(ms.parasite.pres, abbrev.names=FALSE)

## simplified model sociality, degree and body size are highly
## colinear. Small bees tend to be social and generalized

bee.parasite.pres.mod <- glmer(ParasitePresence~
                           scale(TotalAbundance) +
                           scale(r.degree) +
                           ## scale(MeanITD)+
                           ## Sociality +
                           (1|Site),
                           family="binomial",
                           glmerControl(optimizer="bobyqa"),
                           data=spec.wild.sub,
                           na.action = "na.fail")

vif(bee.parasite.pres.mod)
summary(bee.parasite.pres.mod)

## the rate of parasitism in generalist bees < specialist bees
## large bees > small bees
## social bees > solitary
## higher abundaunce lowers parasitism rates, suggesting dilution

## parasite richness within a bee
bee.parasite.rich.mod <- glmer(cbind(ParasiteRichness,
                                     PossibleParasite)~
                           AdjHR +
                           TransectType +scale(SFBloom) +
                           scale(TotalAbundance) +
                           scale(r.degree) +
                           scale(MeanITD)+
                           Sociality +
                           NestLoc +
                           (1|Site),
                           family="binomial",
                           glmerControl(optimizer="bobyqa"),
                           data=spec.wild.sub,
                           na.action = "na.fail")

vif(bee.parasite.rich.mod)
summary(bee.parasite.rich.mod)

ms.parasite.rich <- dredge(bee.parasite.rich.mod)
print(ms.parasite.rich, abbrev.names=FALSE)

## simplified model
bee.parasite.rich.mod <- glmer(cbind(ParasiteRichness,
                                     PossibleParasite)~
                                   scale(TotalAbundance) +
                                   ## scale(r.degree) +
                                   scale(MeanITD)+
                                   (1|Site),
                               family="binomial",
                               glmerControl(optimizer="bobyqa"),
                               data=spec.wild.sub,
                               na.action = "na.fail")

vif(bee.parasite.rich.mod)
summary(bee.parasite.rich.mod)

## same issue of colinearily between degree and body size

## the number of  parasites  in generalist bees < specialist bees
## large bees > small bees
## higher abundaunce lowers parasite richness, suggesting dilution
