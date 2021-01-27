## This script runs the models for bee abundance and richness
rm(list=ls())

## focal bee set by deafult to all if not specified
source("src/initialize.R")
print(focal.bee)

## *************************************************************
## model selection: bee abundunace
## *************************************************************
ys <- c("TotalAbundance",
        "Richness")

## HR proximity and nat hab proximity deleted

by.site$TransectType <- factor(by.site$TransectType,
                               levels=c("SF", "HR", "WM"))


all.mod.vars <-  c("TransectType",
                "scale(Doy)",
                "scale(I(Doy^2))",
                "scale(FloralAbundance)",
                "scale(FloralDiv)",
                "(1|Site)")

## to choose between, 350 or 1000 buffer
 xvars1 <-      c("scale(log(SunflowerCurrent1000))",
                  "scale(log(SunflowerLastYr1000))")
xvars2 <-      c("scale(log(SunflowerCurrent350))",
                "scale(log(SunflowerLastYr350))")

formula1 <-lapply(ys, function(y) {
        as.formula(paste(y, "~",
                         paste(paste(c(all.mod.vars, xvars1),
                                     collapse="+"))))
})

formula2 <-lapply(ys, function(y) {
        as.formula(paste(y, "~",
                         paste(paste(c(all.mod.vars, xvars2),
                                     collapse="+"))))
})

names(formula1) <- names(formula2) <- ys

## *********************************************************************
## full model of wild bee abundance
bee.abund.mod2 <- glmer.nb(formula2[["TotalAbundance"]],
                           glmerControl(optimizer="bobyqa",
                                        optCtrl=list(maxfun=1e6)),
                          data=by.site)

vif(bee.abund.mod2)
summary(bee.abund.mod2)
AIC(bee.abund.mod2)


## having some convergence issues, though sunflower1000 are not
## colinear according to VIF
bee.abund.mod1 <- glmer.nb(formula1[["TotalAbundance"]],
                           glmerControl(optimizer="bobyqa",
                                        optCtrl=list(maxfun=1e6)),
                          data=by.site)

vif(bee.abund.mod1)
summary(bee.abund.mod1)
AIC(bee.abund.mod1)

## 350 buffers have the lower AIC

## *********************************************************************
## full model of wild bee richness
bee.rich.mod2 <- glmer.nb(formula2[["Richness"]],
                           glmerControl(optimizer="bobyqa",
                                        optCtrl=list(maxfun=1e6)),
                          data=by.site)

vif(bee.rich.mod2)
summary(bee.rich.mod2)
AIC(bee.rich.mod2)

## having some convergence issues
bee.rich.mod1 <- glmer.nb(formula1[["Richness"]],
                           glmerControl(optimizer="bobyqa",
                                        optCtrl=list(maxfun=1e6)),
                          data=by.site)

vif(bee.rich.mod1)
summary(bee.rich.mod1)
AIC(bee.rich.mod1)

## 350 buffers have the lower AIC

## *********************************************************************
save(bee.abund.mod2,
    bee.rich.mod2,
     file=sprintf("saved/%s_beeMods.RData",
                  gsub(" ", "", focal.bee)))



mods <- list(bee.abund.mod2,
             bee.rich.mod2)

coeffs <- lapply(mods, function(x) round(coefficients(summary(x)),3))


mapply(function(x, y){
    write.csv(x,
              file=sprintf("saved/tables/beeMods_%s.csv",
                           y))
      write.table(x,
              file=sprintf("saved/tables/beeMods_%s.txt",
                           y), sep="&")
    },
    x=coeffs,
    y=ys[1:2]
)
