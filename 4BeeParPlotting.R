rm(list=ls())
source("src/misc.R")
source("src/predictIntervals.R")
source("src/plotPanels.R")
source("src/CIplotting.R")
source("src/diagnostics.R")
source("src/plotInteractions.R")
library(viridis)
library(boot)
library(car)
library(ggplot2)
library(performance)

focal.bee <- "all"
source('src/initialize.R')

sp.by.site <- read.csv("../sunflower/data/SpbySite.csv")
traits <- read.csv("../sunflower/data/traits.csv")

sp.by.site <- merge(sp.by.site, traits)

sp.by.site <- sp.by.site[!is.na(sp.by.site$MeanITD)
                         & !is.na(sp.by.site$Lecty)
                         & !is.na(sp.by.site$Sociality),]

sp.by.site$Sociality[sp.by.site$Sociality == 1] <- "solitary"
sp.by.site$Sociality[sp.by.site$Sociality == 0] <- "social"


sp.by.site$Lecty[sp.by.site$Lecty == 1] <- "polylectic"
sp.by.site$Lecty[sp.by.site$Lecty == 0] <- "oligolectic"
sp.by.site$Lecty <- factor(sp.by.site$Lecty,
                           levels=c("oligolectic", "polylectic"))

sp.by.site <- sp.by.site[sp.by.site$Year == "2019",]

by.site$TransectType <- factor(by.site$TransectType,
                               levels=c("SF", "HR", "WM"))

cols.sp <- magma(length(unique(sp.by.site$GenusSpecies)))
names(cols.sp) <- sort(unique(sp.by.site$GenusSpecies))

## ************************************************************
## bee abund models
## ************************************************************
load(sprintf('saved/%s_BeeMods.RData', gsub(" ", "", focal.bee)))

top.mod.abund <- bee.abund.mod2
summary(top.mod.abund)

pdf("figures/diagnostics/beeAbund.pdf",
    width = 12, height = 8)
check_model(top.mod.abund)
dev.off()

cols.var <- add.alpha(viridis(3), 0.5)[c(3,1, 2)]
names(cols.var) <- levels(by.site$TransectType)

## day of the year
dd.doy <- expand.grid(Doy=seq(
                          from= min(by.site$Doy),
                          to= max(by.site$Doy),
                          length=20),
                      TransectType=unique(by.site$TransectType),
                      FloralAbundance=mean(by.site$FloralAbundance),
                      SunflowerCurrent350=mean(by.site$SunflowerCurrent350),
                      SunflowerLastYr350=mean(by.site$SunflowerLastYr350),
                      TotalAbundance=0)

doy.pi <- predict.int(top.mod.abund,
                      dd.doy,
                      "TotalAbundance",
                      "poisson")

## sf last year
dd.sflyrprox <- expand.grid(SunflowerLastYr350=seq(
                                from= min(by.site$SunflowerLastYr350),
                                to= max(by.site$SunflowerLastYr350),
                                length=20),
                            Doy=mean(by.site$Doy),
                            TransectType=unique(by.site$TransectType),
                            FloralAbundance=mean(by.site$FloralAbundance),
                            SunflowerCurrent350=mean(by.site$SunflowerCurrent350),
                            TotalAbundance=0)

sflyrprox.pi <- predict.int(top.mod.abund,
                            dd.sflyrprox,
                            "TotalAbundance",
                            "poisson")

## sig only models
pdf.f(plotSigModelsBeeAbund,
      file=sprintf("figures/mods/%s_beeAbund.pdf", focal.bee),
      height=4, width=9)

pdf.f(plotDiagAbund,
      file=file.path('figures/diagnostics/beeAbund.pdf'),
      height=7, width=3)

## ************************************************************
## bee richness models
## ************************************************************
top.mod.rich <- bee.rich.mod2
summary(top.mod.rich)

dd.rich.sflyrprox <- expand.grid(SunflowerLastYr350=seq(
                                     from= min(by.site$SunflowerLastYr350),
                                     to= max(by.site$SunflowerLastYr350),
                                     length=20),
                                 Doy=mean(by.site$Doy),
                                 TransectType=unique(by.site$TransectType),
                                 FloralAbundance=mean(by.site$FloralAbundance),
                                 SunflowerCurrent350=mean(by.site$SunflowerCurrent350),
                                 Richness=0)

sflyrprox.rich.pi <- predict.int(top.mod.rich,
                                 dd.rich.sflyrprox,
                                 "Richness")

pdf("figures/diagnostics/beeRich.pdf",
    width = 12, height = 8)
check_model(top.mod.rich)
dev.off()

## ************************************************************
## Parasite presence models wild bees
## ************************************************************
load(sprintf('saved/%s_parMods.RData', gsub(" ", "", focal.bee)))

## Total abundance
top.mod <- parasite.pres.mod
summary(top.mod)

cols.var <- add.alpha(c("darkolivegreen", "dodgerblue"), 0.5)
names(cols.var) <- unique(sp.by.site$Lecty)

by.site1 <- sp.by.site
colnames(by.site1)[colnames(by.site1) == "SiteParasitismRate"] <-
    "ParasitePresence"

pdf.f(plotInteractions, file="figures/interactions.pdf",
      height=8, width=6)

## ITD
sp.pchs  <- as.numeric(by.site1$Lecty) + 15
names(sp.pchs) <- by.site1$GenusSpecies

dd.itd <- expand.grid(MeanITD=seq(
                          from= min(sp.by.site$MeanITD),
                          to= max(sp.by.site$MeanITD),
                          length=20),
                      TotalAbundance=mean(sp.by.site$TotalAbundance),
                      FloralAbundance=mean(sp.by.site$FloralAbundance,
                                          na.rm=TRUE),
                      Lecty=unique(sp.by.site$Lecty),
                      Sociality= unique(sp.by.site$Sociality),
                      cats="all",
                      ParasitePresence=0)

itd.pi <- predict.int(top.mod,
                      dd.itd,
                      "ParasitePresence",
                      "binomial")

itd.pi <- itd.pi[itd.pi$Sociality == "solitary",]


pdf.f(plotSigModelsParPresSp,
      file=sprintf("figures/mods/%s_wildParasitismSp.pdf", focal.bee),
      height=6, width=6.7)


## richness diagnostics
pdf("figures/diagnostics/parRichLn.pdf",
    width = 12, height = 8)
check_model(parasite.rich.mod.ln)
dev.off()

## richness diagnostics
pdf("figures/diagnostics/parRichPoi.pdf",
    width = 12, height = 8)
check_model(parasite.rich.mod.poi)
dev.off()

