## setwd("~/Dropbox/sunflower")
setwd("analysis/parasiteCommunity")
rm(list=ls())
library(lme4)
library(lmerTest)
library(car)
library(MuMIn)

## focal bee set by deafult to all if not specified
source("src/initialize.R")
print(focal.bee)

## *************************************************************
## model selection: bee abundunace
## *************************************************************

ys <- c("TotalAbundance",
        "Richness")

xvars <-      c("TransectType",
                "SFBloom",
                "scale(Doy)",
                "scale(I(Doy^2))",
                "scale(log(Nat350))",
                "scale(log(Nat1000))",
                "scale(log(HR350))",
                "scale(log(HR1000))",
                "scale(log(SunflowerCurrent1000))",
                "scale(log(SunflowerLastYr1000))",
                "scale(log(SunflowerCurrent350))",
                "scale(log(SunflowerLastYr350))",
                "scale(FloralAbundance)",
                "scale(FloralDiv)",
                "(1|Site)")


formulas <-lapply(ys, function(y) {
        as.formula(paste(y, "~",
                         paste(paste(xvars,
                                     collapse="+"))))
})

names(formulas) <- ys

## *********************************************************************
## full model of wild bee abundance
bee.abund.mod <- glmer.nb(formulas[["TotalAbundance"]],
                          na.action = "na.fail",
                          data=by.site)

vif(bee.abund.mod)
## exclude the different gaussian decays from being included in the
## same model
ms.bee.abund <- dredge(bee.abund.mod,
                       subset =
                !("scale(log(Nat1000))" && "scale(log(Nat350))") &&
                !("scale(log(HR350))" && "scale(log(HR1000))") &&
                !("scale(log(SunflowerCurrent1000))" &&
                  "scale(log(SunflowerCurrent350))") &&
                  !("scale(log(SunflowerLastYr1000))" &&
                    "scale(log(SunflowerLastYr350))"))
## model average within 2 AICc of the min
ma.bee.abund <- model.avg(ms.bee.abund, subset= delta < 2,
                          revised.var = TRUE)
## *********************************************************************
## full model of wild bee richness
bee.rich.mod <- lmer(formulas[["Richness"]],
                          na.action = "na.fail",
                          data=by.site)
## exclude the different gaussian decays from being included in the
## same model
ms.bee.rich <- dredge(bee.rich.mod,
                       subset =
                 !("scale(log(Nat1000))" && "scale(log(Nat350))") &&
                !("scale(log(HR350))" && "scale(log(HR1000))") &&
                !("scale(log(SunflowerCurrent1000))" &&
                  "scale(log(SunflowerCurrent350))") &&
                  !("scale(log(SunflowerLastYr1000))" &&
                  "scale(log(SunflowerLastYr350))"))
ma.bee.rich <- model.avg(ms.bee.rich, subset= delta < 2,
                          revised.var = TRUE)



mapply(function(x, y)
    write.csv(x,
              file=sprintf("saved/tables/beeMod_%s.csv",
                           y)),
    x=list(summary(ma.bee.abund)$coefmat.subset,
           summary(ma.bee.rich)$coefmat.subset),
    y=ys
    )

mapply(function(x, y)
    write.table(x,
              file=sprintf("saved/tables/beeMod_%s.txt",
                           y), sep="&"),
    x=list(round(summary(ma.bee.abund)$coefmat.subset, 3),
           round(summary(ma.bee.rich)$coefmat.subset, 3)),
    y=ys
    )

save(ma.bee.abund,
     ms.bee.abund,
     ma.bee.rich,
     ms.bee.rich,
     file=sprintf("saved/%s_beeMods.RData",
                  gsub(" ", "", focal.bee)))
