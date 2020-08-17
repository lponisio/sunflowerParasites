runMantelSiteType <- function(site.in.type,
                              dis.bee, dist.plant,
                              dist.geo, dist.parasite){

    dist.bee.type <- dist.bee[rownames(dist.bee) %in% site.in.type,
                              colnames(dist.bee) %in% site.in.type]

    dist.plant.type <- dist.plant[rownames(dist.plant) %in% site.in.type,
                                  colnames(dist.plant) %in% site.in.type]

    dist.geo.type <- dist.geo[rownames(dist.geo) %in% site.in.type,
                              colnames(dist.geo) %in% site.in.type]

    dist.parasite.type <- dist.parasite[rownames(dist.parasite) %in% site.in.type,
                                        colnames(dist.parasite) %in% site.in.type]

    out <- MRM(as.dist(dist.parasite.type) ~ as.dist(dist.bee.type) +
                   as.dist(dist.plant.type) + as.dist(dist.geo.type),  nperm=10^4)
    return(out)

}
