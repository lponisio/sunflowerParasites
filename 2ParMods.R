## setwd("~/Dropbox/sunflower")
setwd("analysis/parasiteCommunity")
rm(list=ls())
library(lme4)
library(lmerTest)
library(car)
library(MuMIn)
source("src/initialize.R")
print(focal.bee)

ys <- c("ParasitePresence",
        "cbind(ParasiteRichness, PossibleParasite)",
        parasites)

if(focal.bee == "all" |focal.bee == "NotLasioMel"){
    xvars <-   c("TransectType",
                 "SFBloom",
                 "scale(TotalAbundance)",
                 "Sociality",
                 "scale(Richness)",
                 "scale(r.degree)",
                 "scale(MeanITD)",
                 "scale(FloralAbundance)",
                 "scale(FloralDiv)",
                 "(1|GenusSpecies)",
                 "(1|Genus)")

} else{
    ## don't include traits in species specific models, add sp
    ## abundance
    xvars <-   c("TransectType", "scale(SFBloom)",
                 "scale(TotalAbundance)",
                 "scale(Richness)",
                 "scale(FloralAbundance)",
                 "scale(FloralDiv)")
}

formulas <-lapply(ys, function(y) {
        as.formula(paste(y, "~",
                         paste(paste(xvars,
                                     collapse="+"),
                               "(1|Site)", sep="+")))
})

names(formulas) <- c("Presence", "Richness", parasites)


## *************************************************************
## parasite presence (any parasite)
## *************************************************************
parasite.pres.mod <- glmer(formulas[["Presence"]],
                           family="binomial",
                           glmerControl(optimizer="bobyqa"),
                           data=spec.wild.sub,
                           na.action = "na.fail")

vif(parasite.pres.mod)

## include richness and abundaunce from the same model as they are
## very colinear
ms.parasite.pres <- dredge(parasite.pres.mod,
            subset =  !("scale(Richness)" && "scale(TotalAbundance)"))

ma.parasite.pres <- model.avg(ms.parasite.pres, subset= delta < 2,
                              revised.var = TRUE)

summary(ma.parasite.pres)

## simplified model sociality, degree and body size are highly
## colinear. Small bees tend to be social and generalized


## the rate of parasitism in generalist bees < specialist bees
## large bees > small bees
## social bees > solitary
## higher abundaunce lowers parasitism rates, suggesting dilution

## *************************************************************
## parasite richness within a bee
## *************************************************************

parasite.rich.mod <- glmer(formulas[["Richness"]],
                           family="binomial",
                           glmerControl(optimizer="bobyqa"),
                           data=spec.wild.sub,
                           na.action = "na.fail")

## include richness and abundaunce from the same model as they are
## very colinear
ms.parasite.rich <- dredge(parasite.rich.mod,
         subset =  !("scale(Richness)" && "scale(TotalAbundance)"))


ma.parasite.rich <- model.avg(ms.parasite.rich, subset= delta < 2,
                              revised.var = TRUE)

summary(ma.parasite.rich)
## same issue of colinearily between degree and body size

## the number of  parasites  in generalist bees < specialist bees
## large bees > small bees
## higher abundaunce lowers parasite richness, suggesting dilution



mapply(function(x, y)
    write.csv(x,
              file=sprintf("saved/tables/parasiteMod_%s.csv",
                           y)),
    x=list(summary(ma.parasite.pres)$coefmat.subset,
           summary(ma.parasite.rich)$coefmat.subset),
    y=ys[1:2]
    )

mapply(function(x, y)
    write.table(x,
              file=sprintf("saved/tables/parasiteMod_%s.txt",
                           y), sep="&"),
    x=list(round(summary(ma.parasite.pres)$coefmat.subset, 3),
           round(summary(ma.parasite.rich)$coefmat.subset, 3)),
    y=ys[1:2]
    )



save(ma.parasite.pres,
     ms.parasite.pres,
     ma.parasite.rich,
     ms.parasite.rich,
     file=sprintf("saved/%s_parMods.RData",
                  focal.bee))



## *************************************************************
## parasite specific models
## *************************************************************
runParModel <- function(parasite){
    print(parasite)
    parasite.mod <- glmer(formulas[[parasite]],
                               family="binomial",
                               glmerControl(optimizer="bobyqa"),
                               data=spec.wild.sub,
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

par.sums <- lapply(parasite.mods, function(x)
    round(summary(x$ma.parasite)$coefmat.subset ,2))

mapply(function(a,b)
    write.csv(a, file=sprintf("saved/tables/%s.csv", b)),
                 a=par.sums,
    b=names(par.sums))

mapply(function(a,b)
    write.table(a, sep="&",  file=sprintf("saved/tables/%s.txt", b)),
                 a=par.sums,
                 b=names(par.sums))


save(parasite.mods,
     file=sprintf("saved/%s_parasiteSpecific_parMods.RData",
                  focal.bee))


## across the entire community
## Apicystis:
## Ascosphaera: total abundance negative relationship
## CrithidiaSpp:  degree and total abund negative relationship, solitary>social
## CrithidiaBombi:  total abundance negative relationship, + SF bloom
## CrithidiaExpoeki: degree (marginal) and total abund negative relationship
## NosemaCeranae: foral diversity negative, body size +, abundance -
## NosemaBombi: floral richness/floral diversity +, body size +, total abundance -, solitary < social
