
plotSigModelsParPresSite <- function(){
    layout(matrix(1:2, nrow=1))
    par(oma=c(4,4,1,2), mar=c(1,2,1,1),
        mgp=c(1.5,0.5,0))

    by.site1$cats <- "all"

    plot.panel(dats=by.site1,
               new.dd=abund.pi,
               y1="ParasitePresence",
               xs="TotalAbundance",
               col.lines="black",
               col.fill=cols.var,
               ylabel="Parasitism rate",
               plot.x=TRUE,
               factor.var="all",
               factor.var.col="cats",
               cols.points ="black")
    mtext("Wild bee abundance", 1, line=3)


     plot.panel(dats=by.site1,
               new.dd=bloom.pi,
               y1="ParasitePresence",
               xs="FloralDiv",
               col.lines="black",
               col.fill=cols.var,
               ylabel="",
               plot.x=TRUE,
               factor.var="all",
               factor.var.col="cats",
               cols.points ="black")
    mtext("Floral diversity", 1, line=3)

       ## legend("topright", legend=names(cols.var),
       ##     col=cols.var, pch=c(15),
       ##     bty="n", cex=2)


    ## plot.panel(dats=by.site1,
    ##            new.dd=itd.pi,
    ##            y1="ParasitePresence",
    ##            xs="MeanITD",
    ##            col.lines="black",
    ##            col.fill=cols.var,
    ##            ylabel="",
    ##            plot.x=TRUE,
    ##            factor.var="all",
    ##            factor.var.col="cats",
    ##            cols.points =cols.sp)

    ## mtext("Body size (ITD mm)", 1, line=3)

    ##   legend("bottomright", legend=names(cols.sp),
    ##        col=cols.sp, pch=c(15),
    ##        bty="n", cex=.65, ncol=1)



}


plotSigModelsParPresSp <- function(){
      by.site1$cats <- "all"
     layout(matrix(1:2, ncol=1),
           heights=c(3, 8))
        par(oma=c(4,4,1,2), mar=c(0,0,0,0),
            mgp=c(1.5,0.5,0))

    plot(NA, ylim=c(0,1), xlim=c(0,1),
         xaxt="n",
         yaxt="n",
         bty="n",
         ylab="",
         xlab="")

     legend("center", legend=names(cols.sp),
            col=cols.sp,
            pch=16, ncol=3, cex=0.6,
           bty="n")

    plot.panel(dats=by.site1,
               new.dd=itd.pi,
               y1="ParasitePresence",
               xs="MeanITD",
               col.lines="black",
               col.fill=cols.var,
               ylabel="",
               plot.x=TRUE,
               factor.var="all",
               factor.var.col="cats",
               cols.points =cols.sp)

      mtext("Body size (ITD mm)", 1, line=3)
      mtext("Parasitism rate", 2, line=3)

      ## legend("bottomright", legend=names(cols.sp),
      ##      col=cols.sp, pch=c(15),
      ##      bty="n", cex=.65, ncol=1)

}


plotSigModelsBeeAbund <- function(){
    layout(matrix(1:2, nrow=1, byrow=TRUE))
    ## par(oma=c(4,4,1,2), mar=c(4,2,1,1),
    ##     mgp=c(1.5,0.5,0))
      par(oma=c(2,4,1,2), mar=c(4,2,1,1),
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

    ## plot.panel(dats=by.site,
    ##            new.dd=nat.pi,
    ##            y1="TotalAbundance",
    ##            xs="Nat1000",
    ##            col.lines="black",
    ##            col.fill=cols.var,
    ##            ylabel="",
    ##            plot.x=TRUE,
    ##            factor.var=unique(by.site$TransectType),
    ##            factor.var.col="TransectType",
    ##            cols.points ="black")

    ## mtext("Natural habitat proximity", 1, line=3)


       plot.panel(dats=by.site,
               new.dd=sfprox.pi,
               y1="TotalAbundance",
               xs="SunflowerCurrent1000",
               col.lines="black",
               col.fill=cols.var,
               ylabel="",
               plot.x=TRUE,
               factor.var=unique(by.site$TransectType),
               factor.var.col="TransectType",
               cols.points ="black")
    mtext("Sunflower weighted area", 1, line=3)

        legend("topright",
           legend=c("Hedgerow", "Sunflower", "Weedy margin"),
           col=cols.var, pch=c(15),
           bty="n", cex=1)


    ##    plot.panel(dats=by.site,
    ##            new.dd=sf.pi,
    ##            y1="TotalAbundance",
    ##            xs="SFBloom",
    ##            col.lines="black",
    ##            col.fill=cols.var,
    ##            ylabel="",
    ##            plot.x=TRUE,
    ##            factor.var=unique(by.site$TransectType),
    ##            factor.var.col="TransectType",
    ##            cols.points ="black")
    ## mtext("Proportion sunflower field in bloom", 1, line=3)
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
