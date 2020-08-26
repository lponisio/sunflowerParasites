## setwd("~/Dropbox/sunflower")
setwd("analysis/parasiteCommunity")
rm(list=ls())
library(lme4)
library(lmerTest)
library(car)
library(MuMIn)

focal.bee <- "all"

source("src/initialize.R")

ys <- c("ParasitePresence",
        "cbind(ParasiteRichness, PossibleParasite)",
        parasites)

xvars <-   c("TransectType*scale(SFBloom)",
             "scale(TotalAbundance)",
             "scale(Richness)",
             "scale(FloralAbundance)",
             "scale(FloralDiv)",
             "scale(FloralRichness)")

formulas <-lapply(ys, function(y) {
    as.formula(paste(y, "~",
                     paste(paste(xvars,
                                 collapse="+"),
                           "(1|Site)", sep="+")))
})

names(formulas) <- c("Presence", "Richness", parasites)

## *************************************************************
## parasite presence honey bees
## *************************************************************

parasite.pres.mod.hb <- glmer(formulas[[1]],
                              family="binomial",
                              glmerControl(optimizer="bobyqa"),
                              data=hb,
                              na.action = "na.fail")

ms.parasite.pres.hb <- dredge(parasite.pres.mod.hb,
       subset =  !("scale(Richness)" && "scale(TotalAbundance)")&&
            !("scale(FloralRichness)" && "scale(FloralDiv)"))

ma.parasite.pres.hb <- model.avg(ms.parasite.pres.hb,
                                 subset= delta < 2,
                                 revised.var = TRUE)

summary(ma.parasite.pres.hb)

## *************************************************************
## parasite richness honey bees
## *************************************************************

parasite.rich.mod.hb <- glmer(formulas[[2]],
                              family="binomial",
                              glmerControl(optimizer="bobyqa"),
                              data=hb,
                              na.action = "na.fail")

ms.parasite.rich.hb <- dredge(parasite.rich.mod.hb,
         subset =  !("scale(Richness)" && "scale(TotalAbundance)")&&
            !("scale(FloralRichness)" && "scale(FloralDiv)"))

ma.parasite.rich.hb <- model.avg(ms.parasite.rich.hb,
                                 subset= delta < 2,
                                 revised.var = TRUE)

summary(ma.parasite.rich.hb)

save(ma.parasite.rich.hb,
     ms.parasite.rich.hb,
     ma.parasite.pres.hb,
     ms.parasite.pres.hb,
     file="saved/HBMods.RData")

## plotting
source("src/misc.R")
source("src/predictIntervals.R")
source("src/plotPanels.R")
source("src/CIplotting.R")
source("src/diagnostics.R")
library(viridis)
library(boot)
library(pals)

## Total abundance
summary(ma.parasite.pres.hb)
top.mod.hb <- get.models(ms.parasite.pres.hb, 1)[[1]]
summary(top.mod.hb)
r.squaredGLMM(top.mod.hb)

## *************************************************************
## parasite specific models
## *************************************************************
runParModel <- function(parasite){
    print(parasite)
    parasite.mod <- glmer(formulas[[parasite]],
                               family="binomial",
                               glmerControl(optimizer="bobyqa"),
                               data=hb,
                               na.action = "na.fail")

    ## include richness and abundaunce from the same model as they are
    ## very colinear
    ms.parasite <- dredge(parasite.mod,
         subset =  !("scale(Richness)" && "scale(TotalAbundance)") &&
         !("scale(FloralAbundance)" && "scale(FloralRichness)")&&
         !("scale(FloralRichness)" && "scale(FloralDiv)"))

    ma.parasite <- model.avg(ms.parasite, subset= delta < 2,
                                  revised.var = TRUE)

    return(list(  ma.parasite=  ma.parasite,
                ms.parasite=ms.parasite))
}

parasite.mods <- lapply(parasites, runParModel)
names(parasite.mods) <- parasites

par.sums <- lapply(parasite.mods, function(x)
    summary(x$ma.parasite)$coefmat.subset)

mapply(function(a,b)
    write.csv(a, file=sprintf("saved/tables/hb_%s.csv", b)),
                 a=par.sums,
                 b=names(par.sums))

save(parasite.mods,
     file=sprintf("saved/HB_%s_parasiteSpecific_parMods.RData",
                  focal.bee))

## dd.hb.rich <- expand.grid(FloralRichness=seq(
##                             from= min(by.site$FloralRichness, na.rm=TRUE),
##                             to= max(by.site$FloralRichness, na.rm=TRUE),
##                             length=20),
##                         ParasitePresence=0)

## hb.rich.pi <- predict.int(top.mod.hb,
##                         dd.hb.rich,
##                         "ParasitePresence",
##                         "binomial")

## by.site2 <- by.site
## colnames(by.site2)[colnames(by.site2) == "ParasitePresence"] <-
##     "ParasitePresence1"

## colnames(by.site2)[colnames(by.site2) == "HBParasitism"] <-
##     "ParasitePresence"

## cols.var <- add.alpha("goldenrod3", 0.5)
## names(cols.var) <- "all"

## ## sig only modelsx
## pdf.f(plotSigModelsHBPres,
##       file="figures/HBParasitism.pdf",
##       height=5, width=6)
