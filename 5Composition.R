## setwd('~/Dropbox/sunflower')
rm(list=ls())
setwd('analysis/parasiteCommunity')
library(vegan)
library(bipartite)
library(tidyr)
library(dplyr)
library(ecodist)
library(fields)
library(viridis)
library(ggplot2)
library(gplots)

source("src/commPrep.R")
source("src/plotPcoa.R")
source("src/initialize.R")
source("src/runMRM.R")

## **********************************************************
##  distance matrices for parasites, bees, plants
## **********************************************************
spec.raw <- spec.raw[spec.raw$Year == "2019",]
spec.raw.wild <- spec[spec$GenusSpecies != "Apis mellifera",]
spec <- spec[spec$Year == "2019",]

## geo
geo <- unique(spec[,c("Site", "Lat", "Long")])

## bees
bee.comm <- makeCommStruct(spec.raw.wild, "GenusSpecies")
dist.bee <- as.matrix(vegdist(bee.comm$comm,
                              "gower"))

## plants
plant.comm <- makeCommStruct(spec.raw, "PlantGenusSpecies")
dist.plant <- as.matrix(vegdist(plant.comm$comm,
                                "gower"))

## parasites
parasite.comm <- makeStructParasite(spec, parasites)
dist.parasite <- as.matrix(vegdist(parasite.comm$comm,
                                   "gower"))

## geographic distance
dist.geo <- rdist.earth(geo[,c("Long", "Lat")])
rownames(dist.geo) <- colnames(dist.geo)  <- geo$Site

## **********************************************************
##  Perm anova
## **********************************************************
## Permutational Multivariate Analysis of Variance Using Distance
## Matrices

plotCommDistbyGroup(dist.parasite, parasite.comm,
                    "adjsf", "parasite")

plotCommDistbyGroup(dist.bee, bee.comm,
                    "adjsf", "bees")

## not sure if this is necessary given it is rather obvious
plotCommDistbyGroup(dist.plant, plant.comm,
                    "adjsf", "plants")

## make a long table to plot parasitism by site type,
## parasite species
parasite.pre.long <- as.data.frame(parasite.comm$comm)
parasite.pre.long$Site <- rownames(parasite.pre.long)
parasite.long <- pivot_longer(parasite.pre.long, cols=parasites,
                              names_to = "Parasite",
                              values_to = "Parasitism")
parasite.long$adjSF <- spec$AdjSF[match(parasite.long$Site,
                                              spec$Site)]

p <- ggplot(parasite.long, aes(x=Parasite, y=Parasitism,
                               fill=adjSF)) +
  geom_boxplot()
p <- p + scale_fill_viridis_d() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("figures/pcas/parasitism_bySF.pdf", height=5, width=7)

## **********************************************************
##  MRM
## **********************************************************
## multiple regression on distance matrices
dist.bee <- dist.bee[rownames(dist.geo), rownames(dist.geo)]
dist.plant <- dist.plant[rownames(dist.geo), rownames(dist.geo)]
dist.parasite <- dist.parasite[rownames(dist.geo), rownames(dist.geo)]

MRM(as.dist(dist.parasite) ~ as.dist(dist.bee) +
        as.dist(dist.plant) + as.dist(dist.geo),  nperm=10^4)

## within each site type
site.types <- unique(spec[,c("Site", "SiteType")])
hr <- site.types$Site[site.types$SiteType == "HR" |
                      site.types$SiteType == "HR + SF"]

wm <- site.types$Site[site.types$SiteType == "WM" |
                      site.types$SiteType == "WM + SF"]

sf <- site.types$Site[site.types$SiteType == "HR + SF" |
                      site.types$SiteType == "WM + SF"]

runMantelSiteType(sf,
                  dis.bee, dist.plant,
                  dist.geo, dist.parasite)


runMantelSiteType(wm,
                  dis.bee, dist.plant,
                  dist.geo, dist.parasite)


runMantelSiteType(hr,
                  dis.bee, dist.plant,
                  dist.geo, dist.parasite)



## **********************************************************
##  parasitism by sp, site
## **********************************************************

sp.n <- c(table(spec$GenusSpecies))
sp.n <- sp.n[sp.n >=4]

spec <- spec[spec$GenusSpecies %in% names(sp.n),]

## heat maps of # of infected individuals

parasite.comms <- lapply(parasites, getParComm)
names(parasite.comms) <- parasites

for(parasite in parasites){
    pdf.f(plotParasiteMap,
          file=sprintf("figures/heatmaps/%s.pdf", parasite),
          width=8, height=7)
}
