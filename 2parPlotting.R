## setwd('~/Dropbox/sunflower/')
rm(list=ls())
setwd('analysis/parasiteCommunity')
source('src/initialize.R')
source("src/misc.R")
source("src/predictIntervals.R")
source("src/plotPanels.R")
source("src/CIplotting.R")
source("src/diagnostics.R")
library(viridis)

load('saved/parMods.RData')
sp.by.site <- read.csv("~/Dropbox/sunflower/data/SpbySite.csv")
sp.by.site <- sp.by.site[!is.na(sp.by.site$MeanITD),]
sp.by.site$Sociality <- factor(sp.by.site$Sociality,
                               levels=c("social", "solitary"))

cols.sp <- magma(length(unique(sp.by.site$GenusSpecies)))
names(cols.sp) <- unique(sp.by.site$GenusSpecies)

cols.var <- add.alpha(c("grey47", "grey87"), 0.8)
names(cols.var) <- unique(sp.by.site$Sociality)

## ************************************************************
## Parasite presence models
## ************************************************************

## Total abundance
top.mod <- get.models(ms.parasite.pres, 6)[[1]]

dd.abund <- expand.grid(TotalAbundance=seq(
                            from= min(sp.by.site$TotalAbundance),
                            to= max(sp.by.site$TotalAbundance),
                            length=20),
                        r.degree=mean(sp.by.site$r.degree, na.rm=TRUE),
                        MeanITD=mean(sp.by.site$MeanITD, na.rm=TRUE),
                        SFBloom=mean(sp.by.site$SFBloom, na.rm=TRUE),
                        Sociality=unique(sp.by.site$Sociality),
                        ParasitePresence=0)

abund.pi <- predict.int(top.mod,
                        dd.abund,
                        "ParasitePresence",
                        "binomial")
by.site1 <- sp.by.site
colnames(by.site1)[colnames(by.site1) == "ParasitePresence"] <-
    "ParasitePresence1"

colnames(by.site1)[colnames(by.site1) == "Parasitism"] <-
    "ParasitePresence"

## SF bloom
dd.bloom <- expand.grid(SFBloom=seq(
                            from= min(sp.by.site$SFBloom),
                            to= max(sp.by.site$SFBloom),
                            length=20),
                        TotalAbundance=mean(sp.by.site$TotalAbundance),
                        MeanITD=mean(sp.by.site$MeanITD),
                        Sociality=unique(sp.by.site$Sociality),
                        r.degree=mean(sp.by.site$r.degree, na.rm=TRUE),
                        ParasitePresence=0)

bloom.pi <- predict.int(top.mod,
                        dd.bloom,
                        "ParasitePresence",
                        "binomial")
## ITD

dd.itd <- expand.grid(MeanITD=seq(
                          from= min(sp.by.site$MeanITD),
                          to= max(sp.by.site$MeanITD),
                          length=20),
                      r.degree=mean(sp.by.site$r.degree, na.rm=TRUE),
                      TotalAbundance=mean(sp.by.site$TotalAbundance),
                      Sociality=unique(sp.by.site$Sociality),
                      SFBloom=mean(sp.by.site$SFBloom, na.rm=TRUE),
                      ParasitePresence=0)

itd.pi <- predict.int(top.mod,
                      dd.itd,
                      "ParasitePresence",
                      "binomial")

## r.degree

dd.degree <- expand.grid(r.degree=seq(
                             from= min(sp.by.site$r.degree),
                             to= max(sp.by.site$r.degree),
                             length=20),
                    TotalAbundance=mean(sp.by.site$TotalAbundance),
                         MeanITD=mean(sp.by.site$MeanITD),
                         Sociality=unique(sp.by.site$Sociality),
                         SFBloom=mean(sp.by.site$SFBloom, na.rm=TRUE),
                         ParasitePresence=0)

degree.pi <- predict.int(top.mod,
                         dd.degree,
                         "ParasitePresence",
                         "binomial")

## sig only models
pdf.f(plotSigModels, file="figures/wildParasitism.pdf",
      height=5, width=10)

## "important" but not sig
pdf.f(plotNonSigModels, file="figures/wildParasitism_nosig.pdf",
      height=5, width=10)

## ************************************************************
## Parasite richness models
## ************************************************************

top.mod.rich <- get.models(ms.parasite.rich, 1)[[1]]
