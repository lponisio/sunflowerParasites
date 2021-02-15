## In this script we run the parasite presence models.

rm(list=ls())
source("src/initialize.R")
print(focal.bee)

## parasite response variables
ys <- c("ParasitePresence",
        "cbind(ParasiteRichness, PossibleParasite)",
        parasites)

## parasite explanatory variables
xvars <-   c(
    "scale(TotalAbundance)*scale(FloralAbundance)",
    "Sociality",
    "Lecty",
    "scale(MeanITD)",
    "(1|GenusSpecies)")

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
                           data=spec.wild.sub)

summary(parasite.pres.mod)
vif(parasite.pres.mod)
r.squaredGLMM(parasite.pres.mod)

drop1(parasite.pres.mod, test="Chisq")

## *************************************************************
## parasite richness within a bee
## *************************************************************

parasite.rich.mod <- glmer(formulas[["Richness"]],
                           family="binomial",
                           glmerControl(optimizer="bobyqa"),
                           data=spec.wild.sub)

vif(parasite.rich.mod)
summary(parasite.rich.mod)
r.squaredGLMM(parasite.rich.mod)


drop1(parasite.rich.mod, test="Chisq")

mods <- list(parasite.rich.mod,
             parasite.pres.mod)


coeffs <- lapply(mods, function(x) round(coefficients(summary(x)),3))


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

save(parasite.rich.mod,
     parasite.pres.mod,
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
                          data=spec.wild.sub)
    return(parasite.mod)

}

parasite.mods <- lapply(parasites, runParModel)
names(parasite.mods) <- parasites

par.select <- lapply(parasite.mods, drop1, test="Chisq")

par.sums <- lapply(parasite.mods,
                   function(x) round(coefficients(summary(x)),3))

save(parasite.mods,
     file=sprintf("saved/%s_parasiteSpecific_parMods.RData",
                  focal.bee))

mapply(function(a,b){
    write.csv(a, file=sprintf("saved/tables/%s.csv", b))
    write.table(a, sep="&",  file=sprintf("saved/tables/%s.txt", b))
},
a=par.sums,
b=names(par.sums))
