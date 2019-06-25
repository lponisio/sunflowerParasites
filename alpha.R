## setwd("~/Dropbox/sunflower")
rm(list=ls())
setwd("analysis/parasiteCommunity")
load('../../data/spec.Rdata')
source("src/plotCoeff.R")
library(lme4)

parasites <- c("Aspergillus",
               "Ascosphaera", "Apicystis", "Crithidiaspp", "Cexpoeki",
               "Cbombi",  "Nbombi", "Nceranae")

not.path.screen <- apply(spec[, parasites], 1,
                         function(x) all(is.na(x)))
spec <- spec[!not.path.screen,]

spec$ParasiteRichness <- rowSums(spec[, parasites],
                                 na.rm=TRUE)
spec$PossibleParasite <- apply(spec[, parasites], 1,
                               function(x) sum(!is.na(x)))
spec$ParasitePresence <- (spec$ParasiteRichness >= 1)*1


## The richness of pathogens, how many pathogens out of total possible
## successful screenings
mod.rich  <- glm(cbind(spec$ParasiteRichness,
          spec$PossibleParasite)~spec$GenusSpecies,
    family="binomial")

summary(mod.rich)

plotCoeffs(mod.rich, spec, "Richness",  "Parasite Richness", binom=TRUE)

## are you sick at all?
mod.pres  <- glm(spec$ParasitePresence~spec$GenusSpecies,
    family="binomial")

summary(mod.pres)

plotCoeffs(mod.pres, spec, "Presence",  "P(Parasite)",
           adj1=0.02, binom=TRUE)
