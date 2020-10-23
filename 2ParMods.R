## In this script we run the parasite presence models.

rm(list=ls())
source("src/initialize.R")
print(focal.bee)

## parasite response variables
ys <- c("ParasitePresence",
        "cbind(ParasiteRichness, PossibleParasite)",
        parasites)

## parasite explanatory variables
    xvars <-   c("SFBloom",
                 "scale(TotalAbundance)",
                 "Sociality",
                 "scale(Richness)",
                 "scale(r.degree)",
                 "scale(MeanITD)",
                 "scale(FloralAbundance)",
                 "scale(FloralDiv)",
                 "(1|GenusSpecies)",
                 "(1|Genus)")

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


## exclude richness and abundaunce from the same model as they are
## very colinear
ms.parasite.pres <- dredge(parasite.pres.mod,
            subset =  !("scale(Richness)" && "scale(TotalAbundance)"))

ma.parasite.pres <- model.avg(ms.parasite.pres, subset= delta < 2,
                              revised.var = TRUE)

summary(ma.parasite.pres)

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

mods <- list(ma.parasite.pres,
             ma.parasite.rich)
coeffs <- lapply(mods, sumMSdredge)

## make some nice tables for ms
mapply(function(x, y){
    write.csv(x,
              file=sprintf("saved/tables/parasiteMod_%s.csv",
                           y))
      write.table(x,
              file=sprintf("saved/tables/parasiteMod_%s.txt",
                           y), sep="&")
    },
    x=coeffs,
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

par.sums <- lapply(parasite.mods,
                   function(x) sumMSdredge(x$ma.parasite))

save(parasite.mods,
     file=sprintf("saved/%s_parasiteSpecific_parMods.RData",
                  focal.bee))

mapply(function(a,b){
    write.csv(a, file=sprintf("saved/tables/%s.csv", b))
    write.table(a, sep="&",  file=sprintf("saved/tables/%s.txt", b))
   },
                 a=par.sums,
    b=names(par.sums))
