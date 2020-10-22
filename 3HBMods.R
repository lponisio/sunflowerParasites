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
             "scale(FloralDiv)")

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
       subset =  !("scale(Richness)" && "scale(TotalAbundance)"))

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
         subset =  !("scale(Richness)" && "scale(TotalAbundance)"))

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
         subset =  !("scale(Richness)" && "scale(TotalAbundance)"))

    ma.parasite <- model.avg(ms.parasite, subset= delta < 2,
                                  revised.var = TRUE)

    return(list(  ma.parasite=  ma.parasite,
                ms.parasite=ms.parasite))
}

parasite.mods <- lapply(parasites, runParModel)
names(parasite.mods) <- parasites
par.sums <- lapply(parasite.mods,
                   function(x) sumMSdredge(x$ma.parasite))


save(parasite.mods,
     file=sprintf("saved/HB_%s_parasiteSpecific_parMods.RData",
                  focal.bee))


mapply(function(a,b){
    write.csv(a, file=sprintf("saved/tables/HB_%s.csv", b))
    write.table(a, sep="&",  file=sprintf("saved/tables/HB_%s.txt", b))
   },
                 a=par.sums,
    b=names(par.sums))

