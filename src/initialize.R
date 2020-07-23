
colnames(hr.area.sum) <- paste0("HR", colnames(hr.area.sum))
hr.area.sum <- as.data.frame(hr.area.sum)
hr.area.sum$Site <- rownames(hr.area.sum)
colnames(nat.area.sum) <- paste0("Nat", colnames(nat.area.sum))
nat.area.sum <- as.data.frame(nat.area.sum)
nat.area.sum$Site <- rownames(nat.area.sum)

by.site <- merge(by.site, hr.area.sum, by="Site")
by.site <- merge(by.site, nat.area.sum, by="Site")


by.site$s.Nat350 <- scale(by.site$Nat350)
by.site$s.Nat1000 <- scale(by.site$Nat1000)
by.site$s.Nat2500 <- scale(by.site$Nat2500)

by.site$s.HR350 <- scale(by.site$HR350)
by.site$s.HR1000 <- scale(by.site$HR1000)
by.site$s.HR2500 <- scale(by.site$HR2500)

by.site$SFBloom <- as.numeric(by.site$SFBloom)
by.site$s.SFBloom <- scale(by.site$SFBloom)

by.site$s.TotalAbundance <- log(by.site$TotalAbundance)

spec <- merge(spec, by.site)

parasites <- c("Apicystis", "Ascosphaera", "CrithidiaSpp",
               "CrithidiaBombi", "CrithidiaExpoeki",
               "NosemaCeranae", "NosemaBombi" )

not.path.screen <- apply(spec[, parasites], 1,
                         function(x) all(is.na(x)))

spec <- spec[!not.path.screen,]

spec$r.degree <- as.numeric(spec$r.degree)
spec$s.r.degree <- scale(spec$r.degree)

## spec$s.TotalAbundance <- scale(spec$TotalAbundance)

## drop Apis
spec.wild <- spec[spec$GenusSpecies != "Apis mellifera",]
spec.wild$SFBloom <- as.numeric(spec.wild$SFBloom)
