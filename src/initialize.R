load("../../data/spec.Rdata")
by.site <- read.csv("../../data/bySite.csv")
load('../../../sunflower_saved/data/nat_HR_buffers.Rdata')
load('../../../sunflower_saved/data/sunflower_buffers.Rdata')

colnames(hr.area.sum) <- paste0("HR", colnames(hr.area.sum))
hr.area.sum <- as.data.frame(hr.area.sum)
hr.area.sum$Site <- rownames(hr.area.sum)

colnames(nat.area.sum) <- paste0("Nat", colnames(nat.area.sum))
nat.area.sum <- as.data.frame(nat.area.sum)
nat.area.sum$Site <- rownames(nat.area.sum)

colnames(sunflower.2017.area.sum) <- paste0("SF2017_",
                               colnames(sunflower.2017.area.sum))
sunflower.2017.area.sum <- as.data.frame(sunflower.2017.area.sum)
sunflower.2017.area.sum$Site <- rownames(sunflower.2017.area.sum)

colnames(sunflower.2018.area.sum) <- paste0("SF2018_",
                               colnames(sunflower.2018.area.sum))
sunflower.2018.area.sum <- as.data.frame(sunflower.2018.area.sum)
sunflower.2018.area.sum$Site <- rownames(sunflower.2018.area.sum)

colnames(sunflower.2019.area.sum) <- paste0("SF2019_",
                               colnames(sunflower.2019.area.sum))
sunflower.2019.area.sum <- as.data.frame(sunflower.2019.area.sum)
sunflower.2019.area.sum$Site <- rownames(sunflower.2019.area.sum)


by.site <- merge(by.site, hr.area.sum, by="Site")
by.site <- merge(by.site, nat.area.sum, by="Site")
by.site <- merge(by.site, sunflower.2019.area.sum, by="Site")
by.site <- merge(by.site, sunflower.2018.area.sum, by="Site")
by.site <- merge(by.site, sunflower.2017.area.sum, by="Site")

by.site$SunflowerCurrent350 <- by.site$SF2019_350
by.site$SunflowerCurrent1000 <- by.site$SF2019_1000
by.site$SunflowerCurrent2500 <- by.site$SF2019_2500

by.site$SunflowerCurrent350[by.site$Year == "2018"] <-
    by.site$SF2018_350[by.site$Year == "2018"]
by.site$SunflowerCurrent1000[by.site$Year == "2018"] <-
    by.site$SF2018_1000[by.site$Year == "2018"]
by.site$SunflowerCurrent2500[by.site$Year == "2018"] <-
    by.site$SF2018_2500[by.site$Year == "2018"]


by.site$SunflowerLastYr350 <- by.site$SF2018_350
by.site$SunflowerLastYr1000 <- by.site$SF2018_1000
by.site$SunflowerLastYr2500 <- by.site$SF2018_2500

by.site$SunflowerLastYr350[by.site$Year == "2018"] <-
    by.site$SF2017_350[by.site$Year == "2018"]
by.site$SunflowerLastYr1000[by.site$Year == "2018"] <-
    by.site$SF2017_1000[by.site$Year == "2018"]
by.site$SunflowerLastYr2500[by.site$Year == "2018"] <-
    by.site$SF2017_2500[by.site$Year == "2018"]

by.site$SFBloom <- as.numeric(by.site$SFBloom)

by.site$LogTotalAbundance <- log(by.site$TotalAbundance)
## by.site$Doy <- as.numeric(strftime(by.site$Date, format = "%j"))


spec <- merge(spec, by.site)

parasites <- c("Apicystis", "Ascosphaera", "CrithidiaSpp",
               "CrithidiaBombi", "CrithidiaExpoeki",
               "NosemaCeranae", "NosemaBombi" )

not.path.screen <- apply(spec[, parasites], 1,
                         function(x) all(is.na(x)))

spec <- spec[!not.path.screen,]

spec$r.degree <- as.numeric(spec$r.degree)

## spec$s.TotalAbundance <- scale(spec$TotalAbundance)

## drop Apis
spec.wild <- spec[spec$GenusSpecies != "Apis mellifera",]
spec.wild <- spec.wild[spec.wild$GenusSpecies != "Osmia californica",]
spec.wild$SFBloom <- as.numeric(spec.wild$SFBloom)

spec.wild.sub <- spec.wild[!is.na(spec.wild$Sociality) &
                           !is.na(spec.wild$MeanITD) &
                           !is.na(spec.wild$r.degree),]



## sadly drops all 2018 data
by.site <- by.site[!is.na(by.site$FloralRichness),]


hb <- spec[spec$GenusSpecies == "Apis mellifera",]
hb$SFBloom <- as.numeric(hb$SFBloom)
