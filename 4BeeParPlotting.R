## setwd('~/Dropbox/sunflowerParasites/')
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

source('src/initialize.R')
print(focal.bee)


sp.by.site <- read.csv("~/Dropbox/sunflower/data/SpbySite.csv")
traits <- read.csv("~/Dropbox/sunflower/data/traits.csv")

sp.by.site <- merge(sp.by.site, traits)

sp.by.site <- sp.by.site[!is.na(sp.by.site$MeanITD),]

sp.by.site$Sociality <- factor(sp.by.site$Sociality,
                               levels=c("social", "solitary"))

sp.by.site <- sp.by.site[sp.by.site$Year == "2019",]

by.site$TransectType <- factor(by.site$TransectType,
                               levels=c("SF", "HR", "WM"))

cols.sp <- viridis(length(unique(sp.by.site$GenusSpecies)))
names(cols.sp) <- sort(unique(sp.by.site$GenusSpecies))

## ************************************************************
## bee abund models
## ************************************************************
load(sprintf('saved/%s_BeeMods.RData', gsub(" ", "", focal.bee)))

top.mod.abund <- bee.abund.mod2
summary(top.mod.abund)

plotDiagAbund <- function(){
    plotDiagnostics(top.mod.abund, by.site)
}
pdf.f(plotDiagAbund,
      file=file.path('figures/diagnostics/beeAbund.pdf'),
height=7, width=3)


cols.var <- add.alpha(viridis(3), 0.5)[c(3,1, 2)]
names(cols.var) <- levels(by.site$TransectType)

## day of the year
dd.doy <- expand.grid(Doy=seq(
                          from= min(by.site$Doy),
                          to= max(by.site$Doy),
                          length=20),
                      TransectType=unique(by.site$TransectType),
                      FloralAbundance=mean(by.site$FloralAbundance),
                      FloralDiv=mean(by.site$FloralDiv),
                      SunflowerCurrent350=mean(by.site$SunflowerCurrent350),
                      SunflowerLastYr350=mean(by.site$SunflowerLastYr350),
                      TotalAbundance=0)

doy.pi <- predict.int(top.mod.abund,
                      dd.doy,
                      "TotalAbundance",
                      "poisson")

## current sunflower surrounding
dd.sfprox <- expand.grid(SunflowerCurrent350=seq(
                          from= min(by.site$SunflowerCurrent350),
                          to= max(by.site$SunflowerCurrent350),
                          length=20),
                         Doy=mean(by.site$Doy),
                 TransectType=unique(by.site$TransectType),
                      FloralAbundance=mean(by.site$FloralAbundance),
                      FloralDiv=mean(by.site$FloralDiv),
                      SunflowerLastYr350=mean(by.site$SunflowerLastYr350),
                      TotalAbundance=0)

sfprox.pi <- predict.int(top.mod.abund,
                      dd.sfprox,
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
                      FloralDiv=mean(by.site$FloralDiv),
                      SunflowerCurrent350=mean(by.site$SunflowerCurrent350),
                      TotalAbundance=0)

sflyrprox.pi <- predict.int(top.mod.abund,
                      dd.sflyrprox,
                      "TotalAbundance",
                      "poisson")

## sig only models
pdf.f(plotSigModelsBeeAbund,
      file=sprintf("figures/mods/%s_beeAbund.pdf", focal.bee),
      height=4, width=11)

pdf.f(plotDiagAbund,
      file=file.path('figures/diagnostics/beeAbund.pdf'),
      height=7, width=3)

## ************************************************************
## bee richness models
## ************************************************************

if(focal.bee == "all" |focal.bee == "NotLasioMel"){
  top.mod.rich <- bee.rich.mod2
summary(top.mod.rich)

    dd.rich.sflyrprox <- expand.grid(SunflowerLastYr350=seq(
                          from= min(by.site$SunflowerLastYr350),
                          to= max(by.site$SunflowerLastYr350),
                          length=20),
                         Doy=mean(by.site$Doy),
                 TransectType=unique(by.site$TransectType),
                      FloralAbundance=mean(by.site$FloralAbundance),
                      FloralDiv=mean(by.site$FloralDiv),
                      SunflowerCurrent350=mean(by.site$SunflowerCurrent350),
                      Richness=0)

    sflyrprox.rich.pi <- predict.int(top.mod.rich,
                               dd.rich.sflyrprox,
                               "Richness")

    pdf.f(plotSigModelsBeeRich,
          file=sprintf("figures/mods/%s_beeRich.pdf", focal.bee),
          height=5, width=6)

    plotDiagRich <- function(){
        plotDiagnostics(top.mod.rich, by.site)
    }
    pdf.f(plotDiagRich,
          file=file.path('figures/diagnostics/beeRich.pdf'),
          height=7, width=3)

}

## ************************************************************
## Parasite presence models wild bees
## ************************************************************
load(sprintf('saved/%s_parMods.RData', gsub(" ", "", focal.bee)))

## Total abundance
top.mod <- parasite.pres.mod
summary(top.mod)

cols.var <- add.alpha("grey", 0.5)
names(cols.var) <- "all"

by.site1 <- sp.by.site
colnames(by.site1)[colnames(by.site1) == "SiteParasitismRate"] <-
    "ParasitePresence"

pdf.f(plotInteractions, file="figures/interactions.pdf",
      height=8, width=6)


## dd.abund <- expand.grid(TotalAbundance=seq(
##                             from= min(sp.by.site$TotalAbundance),
##                             to= max(sp.by.site$TotalAbundance),
##                             length=20),
##                         MeanITD=mean(sp.by.site$MeanITD),
##                         Sociality= levels(sp.by.site$Sociality),
##                         FloralAbundance=mean(sp.by.site$FloralAbundance,
##                                              na.rm=TRUE),
##                         FloralDiv=mean(sp.by.site$FloralDiv,
##                                        na.rm=TRUE),
##                           Lecty=levels(sp.by.site$Lecty)),
##                         cats="all",
##                         ParasitePresence=0)

## abund.pi <- predict.int(top.mod,
##                         dd.abund,
##                         "ParasitePresence",
##                         "binomial")
## abund.pi <- abund.pi[abund.pi$Sociality == "solitary",]

## ## floral diversity
## dd.bloom <- expand.grid(FloralAbundance=seq(
##                             from= min(sp.by.site$FloralAbundance, na.rm=TRUE),
##                             to= max(sp.by.site$FloralAbundance, na.rm=TRUE),
##                             length=20),
##                        MeanITD=mean(sp.by.site$MeanITD),
##                         Sociality= levels(sp.by.site$Sociality),
##                         TotalAbundance=mean(sp.by.site$TotalAbundance,
##                                              na.rm=TRUE),
##                         FloralDiv=mean(sp.by.site$FloralDiv,
##                                        na.rm=TRUE),
##                           Lecty=levels(sp.by.site$Lecty),
##                         cats="all",
##                         ParasitePresence=0)

## bloom.pi <- predict.int(top.mod,
##                         dd.bloom,
##                         "ParasitePresence",
##                         "binomial")

## bloom.pi <- bloom.pi[bloom.pi$Sociality == "solitary",]


## if(focal.bee == "all" |focal.bee == "NotLasioMel"){
## ## ITD
## dd.itd <- expand.grid(MeanITD=seq(
##                           from= min(sp.by.site$MeanITD),
##                           to= max(sp.by.site$MeanITD),
##                           length=20),
##                       TotalAbundance=mean(sp.by.site$TotalAbundance),
##                       FloralAbundance=mean(sp.by.site$FloralAbundance,
##                                            na.rm=TRUE),
##                       FloralDiv=mean(sp.by.site$FloralDiv,
##                                        na.rm=TRUE),
##                           Lecty=levels(sp.by.site$Lecty),
##                         Sociality= levels(sp.by.site$Sociality),
##                       cats="all",
##                       ParasitePresence=0)

## itd.pi <- predict.int(top.mod,
##                       dd.itd,
##                       "ParasitePresence",
##                       "binomial")

## itd.pi <- itd.pi[itd.pi$Sociality == "solitary",]

## ## sig only models
## pdf.f(plotSigModelsParPresSite,
##       file=sprintf("figures/mods/%s_wildParasitismSite.pdf", focal.bee),
##       height=5, width=12)

## pdf.f(plotSigModelsParPresSp,
##       file=sprintf("figures/mods/%s_wildParasitismSp.pdf", focal.bee),
##       height=7, width=6)

## }



## plotDiagPres <- function(){
##     plotDiagnostics(top.mod, spec.wild.sub)
## }
## pdf.f(plotDiagPres,
##       file=file.path('figures/diagnostics/parasitePres.pdf'),
##       height=7, width=3)



## ************************************************************
## Floral div historgrams
## ************************************************************

## by.site$SiteType <- as.factor(by.site$SiteType)
## by.site$AdjHR <- as.factor(by.site$AdjHR)

## pltfloralAbund <- ggplot(by.site, aes(FloralAbundance,
##                                       fill = AdjHR)) +
##     geom_histogram(alpha = 0.5,
##                    position = 'identity') +
##     labs(fill = "", x="Floral abundance") +
##     ylab("Count") +
##     theme(legend.position="top")+
##     scale_colour_viridis_d() +
##     scale_fill_viridis_d() +
##     geom_vline(aes(xintercept =
##                        median(FloralAbundance[AdjHR == "HR"]),
##                    col = AdjHR=="HR"), show.legend = FALSE) +
##     geom_vline(aes(xintercept =
##                        median(FloralAbundance[AdjHR =="no HR"]),
##                    col = AdjHR =="no HR"), show.legend = FALSE) +
##     theme(axis.title.x = element_text(vjust=-1),
##           axis.title.y = element_text(vjust=3),
##           plot.margin = margin(t = 10, r = 10, b = 10, l = 20))

## ggsave("figures/floralAbundHist.pdf", width = 4, height = 4)



## pltfloralAbund <- ggplot(by.site, aes(FloralDiv,
##                                       fill = AdjHR)) +
##     geom_histogram(alpha = 0.5,
##                    position = 'identity') +
##     labs(fill = "", x="Floral diversity (Shannon's)") +
##     ylab("Count") +
##     theme(legend.position="top")+
##     scale_colour_viridis_d() +
##     scale_fill_viridis_d() +
##     geom_vline(aes(xintercept =
##                        median(FloralDiv[AdjHR == "HR"]),
##                    col = AdjHR=="HR"), show.legend = FALSE) +
##     geom_vline(aes(xintercept =
##                        median(FloralDiv[AdjHR =="no HR"]),
##                    col = AdjHR =="no HR"), show.legend = FALSE) +
##     theme(axis.title.x = element_text(vjust=-1),
##           axis.title.y = element_text(vjust=3),
##           plot.margin = margin(t = 10, r = 10, b = 10, l = 20))

## ggsave("figures/floralDivHist.pdf", width = 4, height = 4)

