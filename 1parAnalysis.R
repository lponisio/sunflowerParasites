## setwd("~/Dropbox/sunflower")
setwd("analysis")
rm(list=ls())
library(piecewiseSEM)
library(lme4)

load("../data/spec.Rdata")
by.site <- read.csv("../data/bySite.csv")
load('../../sunflower_saved/data/nat_HR_buffers.Rdata')

colnames(hr.area.sum) <- paste0("HR", colnames(hr.area.sum))
hr.area.sum <- as.data.frame(hr.area.sum)
hr.area.sum$Site <- rownames(hr.area.sum)
colnames(nat.area.sum) <- paste0("Nat", colnames(nat.area.sum))
nat.area.sum <- as.data.frame(nat.area.sum)
nat.area.sum$Site <- rownames(nat.area.sum)

by.site <- merge(by.site, hr.area.sum, by="Site")
by.site <- merge(by.site, nat.area.sum, by="Site")

spec <- merge(spec, by.site)

parasites <- c("Apicystis", "Ascosphaera", "CrithidiaSpp",
               "CrithidiaBombi", "CrithidiaExpoeki",
               "NosemaCeranae", "NosemaBombi" )

not.path.screen <- apply(spec[, parasites], 1,
                         function(x) all(is.na(x)))

spec <- spec[!not.path.screen,]
spec.wild <- spec[spec$GenusSpecies != "Apis mellifera",]


## are bee richness and abundance correlated?
cor.test(by.site$Richness, by.site$TotalAbundance)
## SHUCKS!

## *************************************************************
## make formulas for path analyses
## *************************************************************

## formula for site effects on the bee community
formula.bee <- formula(TotalAbundance ~ scale(Nat2500) +
                           scale(HR1000) +
                           TransectType +
                           AdjHR +
                           scale(SFBloom) +
                           (1|Site))

## formulas for the site effects on parasites

xvar.par.path <- c("AdjHR",
                   "TransectType",
                   "scale(SFBloom)",
                   "scale(TotalAbundance)",
                   "scale(r.degree)")

## *************************************************************

calcMods <- function(this.formula,
                     formula.bee,
                     col.trials,
                     dats,
                     site.char){
    colnames(dats)[colnames(dats) == col.trials] <- "trials"
    mod = psem(
        BeeAbund = glmer(formula.bee,
                         data = site.char,
                         family="poisson"),
        ParPath = glmer(this.formula,
                       data = dats,

library(car)
bee.abund.mod <- glmer(TotalAbundance~ scale(Nat2500) +
                           scale(HR1000) +
                           TransectType +
                           AdjHR +
                           scale(SFBloom) +
                           (1|Site),
                       data=by.site, family="poisson")

vif(bee.abund.mod)
summary(bee.abund.mod)
plot(density(residuals(bee.abund.mod)))

library(lmerTest)
bee.abund.mod <- lmer(Richness~ scale(Nat2500) +
                           scale(HR350) +
                           TransectType +
                           AdjHR +
                           scale(SFBloom) +
                           (1|Site),
                       data=by.site)

vif(bee.abund.mod)
summary(bee.abund.mod)
plot(density(residuals(bee.abund.mod)))



bee.abund.mod <- glmer(cbind(ParasiteRichness, PossibleParasite)~
                           AdjHR +
                           TransectType +
                           scale(SFBloom) +
                           scale(TotalAbundance) +
                           scale(r.degree) +
                           (1|Site),
                       family="binomial",
                       data=spec.wild)

vif(bee.abund.mod)
summary(bee.abund.mod)
plot(density(residuals(bee.abund.mod)))
