
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
    mtext("Wild bee abundance", 1, line=3, cex=1.2)


    plot.panel(dats=by.site1,
               new.dd=bloom.pi,
               y1="ParasitePresence",
               xs="FloralAbundance",
               col.lines="black",
               col.fill=cols.var,
               ylabel="",
               plot.x=TRUE,
               plot.y=FALSE,
               factor.var="all",
               factor.var.col="cats",
               cols.points ="black")
    mtext("Floral abundance", 1, line=3, cex=1.2)

}


plotSigModelsParPresSp <- function(){
    ## by.site1$cats <- "all"

    layout(matrix(1:2, ncol=1),
           heights=c(2, 7))
    par(oma=c(4,5,1,2), mar=c(1,0,0,0),
        mgp=c(1.5,0.5,0), xpd = NA)

    plot(NA, ylim=c(0,1), xlim=c(0,1),
         xaxt="n",
         yaxt="n",
         bty="n",
         ylab="",
         xlab="")

    savefont <- par(font=3)
    legend("center", legend=names(cols.sp),
           col=cols.sp,
           pch=sp.pchs[names(cols.sp)],
           ncol=4, cex=0.545,
           bty="n")
    par(savefont)
    par( xpd = FALSE)
    plot.panel(dats=by.site1,
               new.dd=itd.pi,
               y1="ParasitePresence",
               xs="MeanITD",
               col.lines="black",
               col.fill=cols.var,
               ylabel="",
               plot.x=TRUE,
               factor.var=unique(by.site1$Lecty),
               factor.var.col="Lecty",
               cols.points =cols.sp,
               pchs.points=sp.pchs)

    mtext("Body size (ITD mm)", 1, line=3, cex=1.2)
    mtext("Parasitism rate", 2, line=3, cex=1.2)

    legend("bottomright", legend=names(cols.var),
         col=cols.var, pch=c(2, 1),
         bty="n", cex=1.2, ncol=1)

}


plotSigModelsBeeAbund <- function(){
    options(scipen=-1)
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
               plot.y=TRUE,
               factor.var=unique(by.site$TransectType),
               factor.var.col="TransectType",
               cols.points ="black",
               other.lim=400)
    mtext("Day of the year", 1, line=3, cex=1.2)

       legend("topright",
           legend=c("Sunflower", "Hedgerow", "Weedy margin"),
           col=cols.var, pch=c(15),
           bty="n", cex=1.2)

    plot.panel(dats=by.site,
               new.dd=sflyrprox.pi,
               y1="TotalAbundance",
               xs="SunflowerLastYr350",
               col.lines="black",
               col.fill=cols.var,
               ylabel="",
               plot.x=TRUE,
               plot.y=FALSE,
               factor.var=unique(by.site$TransectType),
               factor.var.col="TransectType",
               cols.points ="black",
               other.lim=400)

    mtext("Previous year's \n  sunflower weighted area", 1, line=3,  cex=1.2)

}






plotSigModelsBeeRich <- function(){
    par(oma=c(4,4,1,2), mar=c(1,2,1,1),
        mgp=c(1.5,0.5,0))
    plot.panel(dats=by.site,
               new.dd=sflyrprox.rich.pi,
               y1="Richness",
               xs="SunflowerLastYr350",
               col.lines="black",
               col.fill=cols.var,
               ylabel="Wild bee richness",
               plot.x=TRUE,
               factor.var=unique(by.site$TransectType),
               factor.var.col="TransectType",
               cols.points ="black")
    mtext("Previous year's sunflower weighted area", 1, line=3)

    legend("topleft",
           legend=c("Sunflower", "Hedgerow", "Weedy margin"),
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
