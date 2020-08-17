
samp2site.spp <- function(site,spp,abund) {
    x <- tapply(abund, list(site=site,spp=spp), sum)
    x[is.na(x)] <- 0
    return(x)
}

makeCommStruct <- function(spec.dat, type){
    ## prep site by species matrix
    prep.comm <- aggregate(spec.dat[, type],
                           list(site= spec.dat$Site,
                                status= spec.dat$SiteType,
                                sp= spec.dat[, type]), length)

    comm <-  samp2site.spp(site= prep.comm$site,
                           spp= prep.comm$sp, abund=
                                                  prep.comm$x)
    sites <- rownames(comm)
    site.type <- spec.dat$SiteType[match(rownames(comm),
                                         spec.dat$Site)]
    adjsf <- spec.dat$AdjSF[match(rownames(comm),
                                  spec.dat$Site)]

    comm <- bipartite::empty(comm)

    return(list(comm=comm,
                ## sites=sites,
                site.type = site.type,
                adjsf = adjsf))
}


makeStructParasite <- function(spec, parasites){
    parasite.pre.comm <- spec[, c("Site", parasites)]

    parasite.pre.comm <- parasite.pre.comm  %>%
        group_by(Site) %>%
        summarise_each(list(mean))

    site.type <- spec$SiteType[match(parasite.pre.comm$Site,
                                     spec$Site)]

    adjsf <- spec$AdjSF[match(parasite.pre.comm$Site,
                             spec$Site)]

    sites <- parasite.pre.comm$Site

    comm <- parasite.pre.comm
    comm$Site <- NULL
    comm <- as.matrix(comm)
    rownames(comm) <- parasite.pre.comm$Site
    list(comm=comm,
         ## sites=sites,
         site.type=site.type,
         adjsf = adjsf)
}


getParComm <- function(parasite){
    parasite <- aggregate(list(Parasite=spec[, parasite]),
                          list(GenusSpecies=spec$GenusSpecies,
                               Site=spec$AltSiteName,
                               SiteType=spec$SiteType),
                          function(x) sum(x) / length(x))

    parasite.comm <- samp2site.spp(parasite$Site,
                                   parasite$GenusSpecies,
                                   parasite$Parasite)
    return(parasite.comm)
}
