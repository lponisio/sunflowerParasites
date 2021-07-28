## In this script we run the parasite presence models.

rm(list=ls())
source("src/initialize.R")

## parasite response variables
ys <- c("ParasitePresence",
        "cbind(ParasiteRichness, PossibleParasite)",
        "log(ParasiteRichness +1)", ## reviwer 3 suggestion
        "ParasiteRichness", ## reviwer 3 suggestion
        parasites)

## parasite explanatory variables
xvars <-   c(
    "scale(TotalAbundance)*scale(FloralAbundance)",
    "Sociality",
    "Lecty",
    "scale(MeanITD)",
    "(1|GenusSpecies)")

## *************************************************************
## reviewer suggestions

## ## reviewer 1 suggestion to match Graystock et al. 2020
## ## https://doi.org/10.1038/s41559-020-1247-x
## ## no evidence for a strong effect of Doy on parasitism

## xvars <-   c(
##     "scale(Doy)",
##     "scale(I(Doy^2))",
##     "scale(FloralAbundance)",
##     "Sociality",
##     "Lecty",
##     "scale(MeanITD)",
##     "(1|GenusSpecies)")

## ## reviewer 2 suggestion 2

## xvars <-   c(
##     "scale(HBParasitismRate)",
##     "scale(MelissodesAbundance)",
##     "scale(TotalAbundance)*scale(FloralAbundance)",
##     "Sociality",
##     "Lecty",
##     "scale(MeanITD)",
##     "(1|GenusSpecies)")

## ## reviewer 2 suggestion 1
## xvars <-   c(
##     "scale(HBNosemaCeranae)",
##     "scale(TotalAbundance)*scale(FloralAbundance)",
##     "Sociality",
##     "Lecty",
##     "scale(MeanITD)",
##     "(1|GenusSpecies)")


## ## reviewer 2 suggestion 2
## xvars <-   c(
##     "scale(BeeDivSimp)*scale(FloralAbundance)",
##     "Sociality",
##     "Lecty",
##     "scale(MeanITD)",
##     "(1|GenusSpecies)")

## *************************************************************

formulas <-lapply(ys, function(y) {
    as.formula(paste(y, "~",
                     paste(paste(xvars,
                                 collapse="+"),
                           "(1|Site)", sep="+")))
})

names(formulas) <- c("Presence", "Richness", "ln_Richness",
                     "Richness_poi", parasites)

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

## for double checking pvalues with a more conservative approach
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
AIC(parasite.rich.mod)


## natural log transformation
parasite.rich.mod.ln <- lmer(formulas[["ln_Richness"]],
                             data=spec.wild.sub)

summary(parasite.rich.mod.ln)
AIC(parasite.rich.mod.ln)

## not transformation, using poission family
parasite.rich.mod.poi <- glmer(formulas[["Richness_poi"]],
                               family="poisson",
                               glmerControl(optimizer="bobyqa"),
                               data=spec.wild.sub)

summary(parasite.rich.mod.poi)
AIC(parasite.rich.mod.poi)

####

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
     parasite.rich.mod.ln,
     parasite.rich.mod.poi,
     parasite.pres.mod,
     file="saved/all_parMods.RData")

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

## par.select <- lapply(parasite.mods, drop1, test="Chisq")

par.sums <- lapply(parasite.mods,
                   function(x) round(coefficients(summary(x)),3))

save(parasite.mods,
     file=("saved/all_parasiteSpecific_parMods.RData"))

mapply(function(a,b){
    write.csv(a, file=sprintf("saved/tables/%s.csv", b))
    write.table(a, sep="&",  file=sprintf("saved/tables/%s.txt", b))
},
a=par.sums,
b=names(par.sums))
