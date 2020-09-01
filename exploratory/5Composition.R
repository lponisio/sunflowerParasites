setwd("~/Dropbox/sunflower")
setwd("analysis/parasiteCommunity")
rm(list=ls())
library(lme4)
library(lmerTest)
library(car)
library(MuMIn)
library(ggplot2)
library(dplyr)
library(gplots)
library(pals)

source("src/initialize.R")

## abundance through time
sp.by.site.wild <- sp.by.site.wild[sp.by.site.wild$Year == 2019,]


plotAbund <- function(){
    spp <- sort(unique(sp.by.site.wild$GenusSpecies))
    cols <- polychrome(length(unique(spp)))
    names(cols) <- sort(unique(spp))

    layout(matrix(1:2, ncol=1),
           heights=c(3, 8))
        par(oma=c(4,4,1,2), mar=c(1,1,1,1),
            mgp=c(1.5,0.5,0))

    plot(NA, ylim=c(0,1), xlim=c(0,1),
         xaxt="n",
         yaxt="n",
         bty="n",
         ylab="",
         xlab="")

    legend("center", legend=spp, col=cols, pch=16, ncol=3, cex=0.6,
           bty="n")


    plot(NA, xlim=range(sp.by.site.wild$Doy, na.rm=TRUE),
         ylim=range(log(sp.by.site.wild$Abundance),
                    na.rm=TRUE),
         xlab="",
         ylab="",
         las=1)

        mtext("Day of the year", 1, line=2)
         mtext("Abundance (log)", 2, line=2)

    for(sp in unique(sp.by.site.wild$GenusSpecies)){
        this.sp <- sp.by.site.wild[sp.by.site.wild$GenusSpecies == sp,]
        points(log(this.sp$Abundance) ~ this.sp$Doy,
               col=cols[sp], pch=16)

    }

}

pdf.f(plotAbund, file="figures/AbundDoy.pdf",
      height=8, width=7)


sp.prop <- sp.by.site.wild  %>%
  group_by(Site, Doy, GenusSpecies) %>%
  summarise(n = sum(Abundance)) %>%
  mutate(Proportion = n / sum(n))

plotPorpAbundSite <- function(){
    sites <- unique(sp.prop$Site)
    cols <- polychrome(length(unique(sp.prop$GenusSpecies)))
    names(cols) <- unique(sp.prop$GenusSpecies)

    panels <- vector("list", length(sites))

    for(i in 1:length(sites)){
        ## these.cols <-
        ##     cols[unique(sp.prop$GenusSpecies[sp.prop$Site== sites[i]])]
        panels[[i]] <-  ggplot(sp.prop[sp.prop$Site== sites[i],],
                               aes(x=Doy, y=Proportion,
                                   fill=GenusSpecies)) +
                               labs(main=sites[i]) +
            theme(legend.position = "none") +
            scale_fill_manual(values = cols) +
        geom_area(alpha=0.6 , size=1, colour="black")
    }

    do.call(grid.arrange, panels)

}

pdf.f(plotPorpAbundSite, file="figures/PropAbund.pdf",
      height=11, width=8.5)


sp.prop.all.sites <- sp.by.site.wild  %>%
    group_by(Doy, GenusSpecies) %>%
    summarise(n = sum(Abundance)) %>%
    mutate(Proportion = n / sum(n))

cols <- polychrome(length(unique(sp.prop.all.sites$GenusSpecies)))
names(cols) <- unique(sp.prop.all.sites$GenusSpecies)
ggplot(sp.prop.all.sites,
       aes(x=Doy, y=Proportion,
           fill=GenusSpecies)) +
    ## scale_fill_manual(values = cols) +
    geom_area(alpha=0.6 , size=1, colour="black")
