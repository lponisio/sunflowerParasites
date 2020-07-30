
plotSigModelsParPres <- function(){
    layout(matrix(1:3, nrow=1))
    par(oma=c(4,4,1,2), mar=c(1,2,1,1),
        mgp=c(1.5,0.5,0))
    plot.panel(dats=by.site1,
               new.dd=abund.pi,
               y1="ParasitePresence",
               xs="TotalAbundance",
               col.lines="black",
               col.fill=cols.var,
               ylabel="Parasitism rate",
               plot.x=TRUE,
               factor.var=unique(sp.by.site$Sociality),
               factor.var.col="Sociality",
               cols.points ="black")
    mtext("Wild bee abundance", 1, line=3)
    legend("topright", legend=names(cols.var),
           col=cols.var, pch=c(15),
           bty="n", cex=2)


     plot.panel(dats=by.site1,
               new.dd=bloom.pi,
               y1="ParasitePresence",
               xs="FloralDiv",
               col.lines="black",
               col.fill=cols.var,
               ylabel="",
               plot.x=TRUE,
               factor.var=unique(sp.by.site$Sociality),
               factor.var.col="Sociality",
               cols.points ="black")
    mtext("Floral diversity", 1, line=3)
       legend("topright", legend=names(cols.var),
           col=cols.var, pch=c(15),
           bty="n", cex=2)


    plot.panel(dats=sp.by.site,
               new.dd=itd.pi,
               y1="ParasitePresence",
               xs="MeanITD",
               col.lines="black",
               col.fill=cols.var,
               ylabel="",
               plot.x=TRUE,
               factor.var=unique(sp.by.site$Sociality),
               factor.var.col="Sociality",
               cols.points =cols.sp)

    mtext("Body size (ITD mm)", 1, line=3)

      legend("bottomright", legend=names(cols.sp),
           col=cols.sp, pch=c(15),
           bty="n", cex=.65, ncol=1)



}



plotSigModelsBeeAbund <- function(){
    layout(matrix(1:4, nrow=2, byrow=TRUE))
    par(oma=c(4,4,1,2), mar=c(4,2,1,1),
        mgp=c(1.5,0.5,0))
    plot.panel(dats=by.site,
               new.dd=doy.pi,
               y1="TotalAbundance",
               xs="Doy",
               col.lines="black",
               col.fill=cols.var,
               ylabel="Wild bee abundance",
               plot.x=TRUE,
               factor.var=unique(by.site$TransectType),
               factor.var.col="TransectType",
               cols.points ="black")
    mtext("Day of the year", 1, line=3)

    plot.panel(dats=by.site,
               new.dd=nat.pi,
               y1="TotalAbundance",
               xs="Nat1000",
               col.lines="black",
               col.fill=cols.var,
               ylabel="",
               plot.x=TRUE,
               factor.var=unique(by.site$TransectType),
               factor.var.col="TransectType",
               cols.points ="black")
    legend("topright",
           legend=c("Hedgerow", "Sunflower", "Weedy margin"),
           col=cols.var, pch=c(15),
           bty="n", cex=2)

    mtext("Natural habitat proximity", 1, line=3)


       plot.panel(dats=by.site,
               new.dd=sfprox.pi,
               y1="TotalAbundance",
               xs="SunflowerCurrent1000",
               col.lines="black",
               col.fill=cols.var,
               ylabel="Wild bee abundance",
               plot.x=TRUE,
               factor.var=unique(by.site$TransectType),
               factor.var.col="TransectType",
               cols.points ="black")
    mtext("Sunflower proximity", 1, line=3)


       plot.panel(dats=by.site,
               new.dd=sf.pi,
               y1="TotalAbundance",
               xs="SFBloom",
               col.lines="black",
               col.fill=cols.var,
               ylabel="",
               plot.x=TRUE,
               factor.var=unique(by.site$TransectType),
               factor.var.col="TransectType",
               cols.points ="black")
    mtext("Proportion sunflower field in bloom", 1, line=3)
}






plotSigModelsBeeRich <- function(){
    par(oma=c(4,4,1,2), mar=c(1,2,1,1),
        mgp=c(1.5,0.5,0))
    plot.panel(dats=by.site,
               new.dd=doy.rich.pi,
               y1="Richness",
               xs="Doy",
               col.lines="black",
               col.fill=cols.var,
               ylabel="Wild bee richness",
               plot.x=TRUE,
               factor.var=unique(by.site$TransectType),
               factor.var.col="TransectType",
               cols.points ="black")
    mtext("Day of the year", 1, line=3)

    legend("topleft",
           legend=c("Hedgerow", "Sunflower", "Weedy margin"),
           col=cols.var, pch=c(15),
           bty="n")

}




plotSigModelsHBPres <- function(){
    par(oma=c(4,4,1,2), mar=c(1,2,1,1),
        mgp=c(1.5,0.5,0))
    by.site2$cats <- "all"
    hb.rich.pi$cats <- "all"
    plot.panel(dats=by.site2,
               new.dd=hb.rich.pi,
               y1="ParasitePresence",
               xs="FloralRichness",
               col.lines="black",
               col.fill=cols.var,
               ylabel="Honey bee parasitism",
               plot.x=TRUE,
               factor.var="all",
               factor.var.col="cats",
               cols.points ="black")
    mtext("FloralRichness", 1, line=3)

}
