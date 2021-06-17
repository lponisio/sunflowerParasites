library(lme4)
library(car)
library(MuMIn)

load("data/spec.Rdata")
by.site <- read.csv("data/bySite.csv")
sp.by.site <- read.csv("data/SpbySite.csv")

source("src/misc.R")
source("src/calcCoeffTable.R")

args <- commandArgs(trailingOnly=TRUE)

spec.raw <- spec

parasites <- c("Apicystis", "Ascosphaera", "CrithidiaSpp",
               "CrithidiaBombi", "CrithidiaExpoeki",
               "NosemaCeranae", "NosemaBombi" )

spec <- spec[spec$Apidae == 1,]

## drop managed bees
managed.bees  <- c("Apis mellifera", "Osmia californica")

spec.wild <- spec[!spec$GenusSpecies %in% managed.bees,]

spec.wild.sub <- spec.wild[spec.wild$Year == "2019",]


## sadly drops all 2018 data
by.site <- by.site[!is.na(by.site$FloralRichness),]

by.site$TransectType <- factor(by.site$TransectType,
                               levels=c("HR", "SF", "WM"))


## site by species data for plotting
sp.by.site <- sp.by.site[!is.na(sp.by.site$MeanITD),]
sp.by.site$Sociality <- factor(sp.by.site$Sociality,
                               levels=c("social", "solitary"))

sp.by.site <- sp.by.site[sp.by.site$Year == "2019", ]

sp.by.site.wild <-
    sp.by.site[!sp.by.site$GenusSpecies %in% managed.bees,]

spec.wild$SFBloom <- as.numeric(spec.wild$SFBloom)
spec.wild.sub$SFBloom <- as.numeric(spec.wild.sub$SFBloom)

spec.wild.sub$MelissodesAbundance[is.na(spec.wild.sub$MelissodesAbundance)] <- 0

## honey bee data
hb <- spec[spec$GenusSpecies == "Apis mellifera",]

hb$SFBloom <- as.numeric(hb$SFBloom)
hb <- hb[hb$Year == "2019", ]

