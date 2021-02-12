

plotInteractions <- function(){

    layout(matrix(c(1,1, 2,2), nrow=2, byrow=TRUE),
           heights=c(1,1))

    par(oma=c(2, 5, 2, 1),
        mar=c(5, 2, 1, 1.5), cex.axis=1.5)
    nbreaks <- 20

    ## top historgram
    h <- hist(sp.by.site$FloralAbundance, breaks=nbreaks, plot=FALSE)
    cols <- rev(viridis(length(h$density)))
    plot(h, col=cols,
         xlab="", main="", ylab="", las=1)

    ## shade hedgerow data
    site.hr <- unique(spec$Site[spec$AdjHR == "HR"])
    h.hedgerow <-
        hist(sp.by.site$FloralAbundance[sp.by.site$Site %in% site.hr],
             breaks=h$breaks, plot=FALSE)
    plot(h.hedgerow,
         xlab="", main="", ylab="", yaxy="n", xaxt="n", add=TRUE, density=15)
    abline(v=mean(sp.by.site$FloralAbundance), lty=2, col="red", lwd=3)

    legend("topright", pch=c(0,7), legend=c("Weedy margin", "Hedgerow"),
           bty="n", cex=1.3)

    real.mids <- c(h$mids[h$mids < 1675], max(h$mids))
      cols <- cols[h$mids %in% real.mids]

    ## real.mids <- h$mids
    abund.hist.mids <- (real.mids - mean(sp.by.site$FloralAbundance))/
        sd(sp.by.site$FloralAbundance)

    mtext("Frequency", 2, line=4, cex=1.3)
    mtext("Non-crop floral abundance", 1, line=3, cex=1.3)

    ## interactions of floral resources and bee abundance
    plot(NA, ylim=c(0, 1), xlim=range(sp.by.site$TotalAbundance), las=1,
         ylab="", xlab="")
    mtext("Parasitism rate", 2, line=4, cex=1.3)
    mtext("Wild bee abundance", 1, line=3, cex=1.3)
    ## legend("topleft", legend="(a)", bty="n", cex=1.2)

    mod.coeffs <-  coefficients(summary(top.mod))[,1]

    for(i in 1:length(abund.hist.mids)){
        curve(inv.logit(mod.coeffs["(Intercept)"] +
                        mod.coeffs["scale(TotalAbundance)"] * x +
                        mod.coeffs["scale(FloralAbundance)"] * abund.hist.mids[i] +
                        mod.coeffs["scale(TotalAbundance):scale(FloralAbundance)"] * x * abund.hist.mids[i]),
              from=range(sp.by.site$TotalAbundance)[1],
              to=range(sp.by.site$TotalAbundance)[2],
              col=cols[i],
              lwd=2,
              add=TRUE)
    }


       curve(inv.logit(mod.coeffs["(Intercept)"] +
                        mod.coeffs["scale(TotalAbundance)"] * x +
                        mod.coeffs["scale(FloralAbundance)"] * 0 +
                        mod.coeffs["scale(TotalAbundance):scale(FloralAbundance)"] * x * 0),
              from=range(sp.by.site$TotalAbundance)[1],
              to=range(sp.by.site$TotalAbundance)[2],
             col="red",
             lty="dashed",
              lwd=2,
              add=TRUE)

}



pdf.f(plotInteractions, file="figures/interactions.pdf",
      height=8, width=6)
