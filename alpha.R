## setwd("~/Dropbox/sunflower")
rm(list=ls())
setwd("analysis/parasiteCommunity")
load('../../data/spec.Rdata')
source("src/plotCoeff.R")
library(lme4)
library(RColorBrewer)

parasites <- c("Apicystis",         "Ascosphaera",
               "Aspergillus",       "CrithidiaSpp",
               "CrithidiaBombi",    "CrithidiaExpoeki",
               "NosemaCeranae",     "NosemaBombi" )

## not.enough.data <- c("Melissodes tepida timberlakei",
##                      "Melissodes stearnsi",
##                      "Melissodes communis alopex",
##                      "Melissodes lupina")

not.path.screen <- apply(spec[, parasites], 1,
                         function(x) all(is.na(x)))

spec <- spec[!not.path.screen,]

spec[, parasites][spec[, parasites] == ""] <- 0
spec[, parasites] <- apply(spec[, parasites], 2,  as.numeric)

## spec <- spec[!spec$GenusSpecies %in% not.enough.data,]

spec$ParasiteRichness <- rowSums(spec[, parasites],
                                 na.rm=TRUE)
spec$PossibleParasite <- apply(spec[, parasites], 1,
                               function(x) sum(!is.na(x)))
spec$ParasitePresence <- (spec$ParasiteRichness >= 1)*1

## The richness of pathogens, how many pathogens out of total possible
## successful screenings

spec$SiteType <- factor(spec$SiteType,
                        levels= c("HR", "HR + SF", "WM", "WM + SF"))

spec.wild <- spec[spec$GenusSpecies != "Apis mellifera",]

mod.rich  <- glmer(cbind(spec.wild$ParasiteRichness,
                         spec.wild$PossibleParasite)~
                       scale(spec.wild$r.degree)*spec.wild$SiteType +
                       (1|Site),
                   family="binomial",
                   data=spec.wild)

## year did not have a sig effect

summary(mod.rich)

## plotCoeffs(mod.rich, spec, "Richness",  "Parasite Richness", binom=TRUE)

## are you sick at all?
mod.pres  <- glmer(cbind(spec.wild$ParasitePresence,
                         spec.wild$PossibleParasite)~
                       scale(spec.wild$r.degree)*spec.wild$SiteType +
                       (1|Site),
                   family="binomial",
                   data=spec.wild)


summary(mod.pres)

## plotCoeffs(mod.pres, spec, "Presence",  "Parasite prevalence",
##            adj1=0.02, binom=TRUE)
