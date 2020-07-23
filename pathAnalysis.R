setwd("~/Dropbox/sunflower")
setwd("analysis/parasiteCommunity")
rm(list=ls())
library(devtools)
install_github("jslefche/piecewiseSEM")
library(piecewiseSEM)
library(lme4)
library(lmerTest)

source("src/initialize.R")

## are bee richness and abundance correlated?
cor.test(by.site$Richness, by.site$TotalAbundance)
## SHUCKS!

## *************************************************************
## make formulas for path analyses
## *************************************************************

spec.wild.sub <- spec.wild[!is.na(spec.wild$MeanITD),]

## formula for site effects on the bee community
formula.bee <- formula(log.TotalAbundance ~ s.Nat1000 +
                           s.HR1000 +
                           s.SunflowerCurrent1000+
                            s.SunflowerLastYr1000+
                           AdjHR +
                           s.SFBloom   +
                           (1|Site))

## formulas for the site effects on parasites

formula.par <- formula(ParasitePresence ~
                           ## AdjHR +
                           ## s.SFBloom   +
                           s.SunflowerCurrent1000+
                           log.TotalAbundance +
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
