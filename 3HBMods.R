## This script runs the parasite models for the honey bess

rm(list=ls())
focal.bee <- "all"
source("src/initialize.R")

## response variables
ys <- c("ParasitePresence",
        "cbind(ParasiteRichness, PossibleParasite)",
        parasites)


## explanatory variables
xvars <-   c("scale(TotalAbundance)*scale(FloralRichness)")

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
                              data=hb)

vif(parasite.pres.mod.hb)
summary(parasite.pres.mod.hb)
r.squaredGLMM(parasite.pres.mod.hb)

## *************************************************************
## parasite richness honey bees
## *************************************************************

parasite.rich.mod.hb <- glmer(formulas[[2]],
                              family="binomial",
                              glmerControl(optimizer="bobyqa"),
                              data=hb)

vif(parasite.rich.mod.hb)
summary(parasite.rich.mod.hb)
r.squaredGLMM(parasite.rich.mod.hb)


save(parasite.rich.mod.hb,
     parasite.rich.mod.hb,
     file="saved/HBMods.RData")

## *************************************************************
## parasite specific models
## *************************************************************
runParModel <- function(parasite){
    print(parasite)
    parasite.mod <- glmer(formulas[[parasite]],
                               family="binomial",
                               glmerControl(optimizer="bobyqa"),
                               data=hb)

    return(parasite.mod)
}

parasite.mods <- lapply(parasites, runParModel)
names(parasite.mods) <- parasites

par.sums <- lapply(parasite.mods,
                   function(x) round(coefficients(summary(x)),3))

par.mod.select <- lapply(parasite.mods,
                         drop1, test="Chisq")


save(parasite.mods,
     file=sprintf("saved/HB_%s_parasiteSpecific_parMods.RData",
                  focal.bee))


mapply(function(a,b){
    write.csv(a, file=sprintf("saved/tables/HB_%s.csv", b))
    write.table(a, sep="&",  file=sprintf("saved/tables/HB_%s.txt", b))
   },
                 a=par.sums,
    b=names(par.sums))

