## setwd('~/Dropbox/sunflower/')
rm(list=ls())
setwd('analysis/parasiteCommunity')
source("src/misc.R")
source("src/predictIntervals.R")
source("src/plotPanels.R")
source("src/CIplotting.R")
source("src/diagnostics.R")
library(viridis)
library(MuMIn)
library(lme4)
library(boot)
library(pals)
library(car)

source('src/initialize.R')
print(focal.bee)


sp.by.site <- read.csv("~/Dropbox/sunflower/data/SpbySite.csv")
sp.by.site <- sp.by.site[!is.na(sp.by.site$MeanITD),]
sp.by.site$Sociality <- factor(sp.by.site$Sociality,
                               levels=c("social", "solitary"))

by.site$TransectType <- factor(by.site$TransectType,
                               levels=c("HR", "SF", "WM"))

cols.sp <- polychrome(length(unique(sp.by.site$GenusSpecies)))
names(cols.sp) <- sort(unique(sp.by.site$GenusSpecies))

## ************************************************************
## bee abund models
## ************************************************************
load(sprintf('saved/%s_BeeMods.RData', gsub(" ", "", focal.bee)))

summary(ma.bee.abund)

top.mod.abund <- get.models(ms.bee.abund, 1)[[1]]
summary(top.mod.abund)
r.squaredGLMM(top.mod.abund)
vif(top.mod.abund)

cols.var <- add.alpha(viridis(3), 0.5)
names(cols.var) <- levels(by.site$TransectType)

## day of the year
dd.doy <- expand.grid(Doy=seq(
                          from= min(by.site$Doy),
                          to= max(by.site$Doy),
                          length=20),
                      Nat1000=mean(by.site$Nat1000),
                      SFBloom=mean(by.site$SFBloom),
                      SunflowerCurrent1000=
                          mean(by.site$SunflowerCurrent1000),
                      TransectType=unique(by.site$TransectType),
                      TotalAbundance=0)

doy.pi <- predict.int(top.mod.abund,
                      dd.doy,
                      "TotalAbundance",
                      "poisson")

## ## natural habitat
## dd.nat <- expand.grid(Nat1000=seq(
##                           from= min(by.site$Nat1000),
##                           to= max(by.site$Nat1000),
##                           length=20),
##                       Doy=mean(by.site$Doy),
##                       SFBloom=mean(by.site$SFBloom),
##                       SunflowerCurrent1000=
##                           mean(by.site$SunflowerCurrent1000),
##                       TransectType=unique(by.site$TransectType),
##                       TotalAbundance=0)

## nat.pi <- predict.int(top.mod.abund,
##                       dd.nat,
##                       "TotalAbundance",
##                       "poisson")


## current sunflower surrounding
dd.sfprox <- expand.grid(SunflowerCurrent1000=seq(
                          from= min(by.site$SunflowerCurrent1000),
                          to= max(by.site$SunflowerCurrent1000),
                          length=20),
                      Doy=mean(by.site$Doy),
                      SFBloom=mean(by.site$SFBloom),
                      Nat1000=mean(by.site$Nat1000),
                      TransectType=unique(by.site$TransectType),
                      TotalAbundance=0)

sfprox.pi <- predict.int(top.mod.abund,
                      dd.sfprox,
                      "TotalAbundance",
                      "poisson")

##  sunflower bloom
## dd.sf <- expand.grid(SFBloom=seq(
##                           from= min(by.site$SFBloom),
##                           to= max(by.site$SFBloom),
##                           length=20),
##                       Doy=mean(by.site$Doy),
##                      SunflowerCurrent1000=
##                          mean(by.site$SunflowerCurrent1000),
##                       Nat1000=mean(by.site$Nat1000),
##                       TransectType=unique(by.site$TransectType),
##                       TotalAbundance=0)

## sf.pi <- predict.int(top.mod.abund,
##                       dd.sf,
##                       "TotalAbundance",
##                       "poisson")

## sig only models
pdf.f(plotSigModelsBeeAbund,
      file=sprintf("figures/mods/%s_beeAbund.pdf", focal.bee),
      height=5, width=10)

## ************************************************************
## bee richness models
## ************************************************************

if(focal.bee == "all" |focal.bee == "NotLasioMel"){
    top.mod.rich <- get.models(ms.bee.rich,1)[[1]]
    summary(top.mod.rich)
    r.squaredGLMM(top.mod.rich)
    vif(top.mod.rich)

    dd.rich.doy <- expand.grid(Doy=seq(
                                   from= min(by.site$Doy),
                                   to= max(by.site$Doy),
                                   length=20),
                               TransectType=unique(by.site$TransectType),
                               SFBloom=mean(by.site$SFBloom),
                               Richness=0)

    doy.rich.pi <- predict.int(top.mod.rich,
                               dd.rich.doy,
                               "Richness")

    pdf.f(plotSigModelsBeeRich,
          file=sprintf("figures/mods/%s_beeRich.pdf", focal.bee),
          height=5, width=6)
}

## ************************************************************
## Parasite presence models wild bees
## ************************************************************
load(sprintf('saved/%s_parMods.RData', gsub(" ", "", focal.bee)))

## Total abundance
summary(ma.parasite.pres)
top.mod <- get.models(ms.parasite.pres, 1)[[1]]
summary(top.mod)
r.squaredGLMM(top.mod)
vif(top.mod)

## cols.var <- add.alpha(viridis(2), 0.5)
## names(cols.var) <- levels(sp.by.site$Sociality)

cols.var <- add.alpha(viridis(1), 0.5)
names(cols.var) <- "all"

by.site1 <- sp.by.site
colnames(by.site1)[colnames(by.site1) == "SiteParasitismRate"] <-
    "ParasitePresence"


dd.abund <- expand.grid(TotalAbundance=seq(
                            from= min(sp.by.site$TotalAbundance),
                            to= max(sp.by.site$TotalAbundance),
                            length=20),
                        MeanITD=mean(sp.by.site$MeanITD),
                        FloralDiv=mean(sp.by.site$FloralDiv, na.rm=TRUE),
                        FloralAbundance=mean(sp.by.site$FloralAbundance, na.rm=TRUE),
                        cats="all",
                        ParasitePresence=0)

abund.pi <- predict.int(top.mod,
                        dd.abund,
                        "ParasitePresence",
                        "binomial")

## floral diversity
dd.bloom <- expand.grid(FloralDiv=seq(
                            from= min(sp.by.site$FloralDiv, na.rm=TRUE),
                            to= max(sp.by.site$FloralDiv, na.rm=TRUE),
                            length=20),
                        TotalAbundance=mean(sp.by.site$TotalAbundance),
                        MeanITD=mean(sp.by.site$MeanITD),
                        FloralAbundance=mean(sp.by.site$FloralAbundance, na.rm=TRUE),
                          cats="all",
                        ParasitePresence=0)

bloom.pi <- predict.int(top.mod,
                        dd.bloom,
                        "ParasitePresence",
                        "binomial")


if(focal.bee == "all" |focal.bee == "NotLasioMel"){
## ITD
dd.itd <- expand.grid(MeanITD=seq(
                          from= min(sp.by.site$MeanITD),
                          to= max(sp.by.site$MeanITD),
                          length=20),
                      TotalAbundance=mean(sp.by.site$TotalAbundance),
                                              FloralAbundance=mean(sp.by.site$FloralAbundance, na.rm=TRUE),
                      cats="all",
                      FloralDiv=mean(sp.by.site$FloralDiv, na.rm=TRUE),
                      ParasitePresence=0)

itd.pi <- predict.int(top.mod,
                      dd.itd,
                      "ParasitePresence",
                      "binomial")

## sig only models
pdf.f(plotSigModelsParPresSite,
      file=sprintf("figures/mods/%s_wildParasitismSite.pdf", focal.bee),
      height=5, width=10)

pdf.f(plotSigModelsParPresSp,
      file=sprintf("figures/mods/%s_wildParasitismSp.pdf", focal.bee),
      height=7, width=6)

}
