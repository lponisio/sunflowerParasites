
plotSigModels <- function(){
    layout(matrix(1:2, nrow=1))
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
               pchs=c(16),
               factor.var=unique(sp.by.site$Sociality),
               factor.var.col="Sociality",
               cols.points ="black")
    mtext("Wild bee abundance", 1, line=3)
    legend("topright", legend=names(cols.var),
           col=cols.var, pch=16,
           bty="n")

    plot.panel(dats=sp.by.site,
               new.dd=itd.pi,
               y1="ParasitePresence",
               xs="MeanITD",
               col.lines="black",
               col.fill=cols.var,
               ylabel="",
               plot.x=TRUE,
               pchs=c(16),
               factor.var=unique(sp.by.site$Sociality),
               factor.var.col="Sociality",
               cols.points =cols.sp)

    mtext("Body size (ITD mm)", 1, line=3)

}




plotNonSigModels <- function(){
    layout(matrix(1:2, nrow=1))
    par(oma=c(4,4,1,2), mar=c(1,2,1,1),
        mgp=c(1.5,0.5,0))
    plot.panel(dats=sp.by.site,
               new.dd=degree.pi,
               y1="ParasitePresence",
               xs="r.degree",
               col.lines="black",
               col.fill=cols.var,
               ylabel="Parasitism rate",
               plot.x=TRUE,
               pchs=c(16),
               factor.var=unique(sp.by.site$Sociality),
               factor.var.col="Sociality",
               cols.points =cols.sp)
    mtext("Degree", 1, line=3)

 plot.panel(dats=by.site1,
           new.dd=bloom.pi,
           y1="ParasitePresence",
           xs="SFBloom",
           col.lines="black",
           col.fill=cols.var,
           ylabel="",
           plot.x=TRUE,
           pchs=c(16),
           factor.var=unique(sp.by.site$Sociality),
           factor.var.col="Sociality",
           cols.points ="black")
  legend("top", legend=names(cols.var),
           col=cols.var, pch=16,
           bty="n")

    mtext("Sunflower bloom", 1, line=3)

}

