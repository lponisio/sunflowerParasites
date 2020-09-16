load("../../data/spec.Rdata")
by.site <- read.csv("../../data/bySite.csv")
sp.by.site <- read.csv("../../data/SpbySite.csv")

source("src/misc.R")

args <- commandArgs(trailingOnly=TRUE)

if(length(args) != 0){
    focal.bee <- args[1]
} else{
    focal.bee <- "all"
}
spec.raw <- spec

parasites <- c("Apicystis", "Ascosphaera", "CrithidiaSpp",
               "CrithidiaBombi", "CrithidiaExpoeki",
               "NosemaCeranae", "NosemaBombi" )

## not.path.screen <- apply(spec[, parasites], 1,
##                          function(x) all(is.na(x)))

spec <- spec[spec$Apidae == 1,]

## drop managed bees
managed.bees  <- c("Apis mellifera", "Osmia californica")

spec.wild <- spec[!spec$GenusSpecies %in% managed.bees,]

spec.wild.sub <- spec.wild[!is.na(spec.wild$Sociality) &
                           !is.na(spec.wild$MeanITD) &
                           !is.na(spec.wild$r.degree),]

spec.wild.sub <- spec.wild.sub[spec.wild.sub$Year == "2019",]

## sadly drops all 2018 data
by.site <- by.site[!is.na(by.site$FloralRichness),]

by.site$TransectType <- factor(by.site$TransectType,
                               levels=c("HR", "SF", "WM"))


## site by species data for plotting
sp.by.site <- sp.by.site[!is.na(sp.by.site$MeanITD),]
sp.by.site$Sociality <- factor(sp.by.site$Sociality,
                               levels=c("social", "solitary"))

sp.by.site <- sp.by.site[sp.by.site$Year == "2019", ]

sp.by.site.wild <-
    sp.by.site[!sp.by.site$GenusSpecies %in% managed.bees,]



if(focal.bee != "all" & focal.bee != "NotLasioMel"){
    by.site  <- sp.by.site.wild[sp.by.site.wild$GenusSpecies ==
                                focal.bee,]
    colnames(by.site)[colnames(by.site) == "Abundance"] <-
        "TotalAbundance"

    spec.wild.sub <-  spec.wild.sub[spec.wild.sub$GenusSpecies ==
                                focal.bee,]
}

if(focal.bee == "NotLasioMel"){
    most.abund <- c("Melissodes agilis", "Lasioglossum incompletum")
    sub.by.sp  <- sp.by.site.wild[!sp.by.site.wild$GenusSpecies %in%
                                  most.abund,]

    total.abund <- aggregate(list(TotalAbundance=sub.by.sp$Abundance),
                             list(Site=sub.by.sp$Site,
                                  Doy=sub.by.sp$Doy,
                                  Year=sub.by.sp$Year),
                             sum)
    by.site$TotalAbundance <- NULL

    by.site$TotalAbundance <-   total.abund$TotalAbundance[
                                                match(paste(by.site$Doy,
                                                            by.site$Site),
                                                      paste(total.abund$Doy,
                                                            total.abund$Site))]
    by.site$TotalAbundance[is.na(by.site$TotalAbundance)] <- 0
}


spec.wild$SFBloom <- as.numeric(spec.wild$SFBloom)
spec.wild.sub$SFBloom <- as.numeric(spec.wild.sub$SFBloom)


## honey bee data
hb <- spec[spec$GenusSpecies == "Apis mellifera",]

hb$SFBloom <- as.numeric(hb$SFBloom)
hb <- hb[hb$Year == "2019", ]

