plot.panel <- function(dats,
                       new.dd,
                       y1,
                       other.lim=0,
                       xs,
                       col.lines,
                       col.fill,
                       ylabel,
                       plot.x=TRUE,
                       plot.y=TRUE,
                       factor.var,
                       factor.var.col,
                       cols.points, ...){
    plotting.loop <- function(){
        for(var in factor.var){
            these.dats <- dats[dats[, factor.var.col]== var,]
            pred.dats <- new.dd[new.dd[, factor.var.col]== var,]
            ys <- data.frame(y=these.dats[,y1], x=these.dats[,xs])
            ## plots means
            if(nrow(pred.dats) > 1){
                ## add fill from ci.up to ci.lb
                polygon(c(pred.dats[,xs],
                          rev(pred.dats[,xs])),
                        c(pred.dats$phi,
                          rev(pred.dats$plo)),
                        col=col.fill[var], border=NA)
                ## plots CI
                lines(x=pred.dats[,xs],
                      y=pred.dats[,y1],
                      col=col.lines,
                      lwd=2)
                lines(x=pred.dats[,xs],
                      y=pred.dats$plo,
                      col=col.lines,
                      lty="dashed")
                lines(x=pred.dats[,xs],
                      y=pred.dats$phi,
                      col=col.lines,
                      lty="dashed")
            }
        }

        ys <- data.frame(y=dats[,y1], x=dats[,xs])
        these.cols <- cols.points[dats$GenusSpecies]
        if(all(is.na(these.cols))) these.cols <- cols.var[dats[, factor.var.col]]
        ## pchs <- 1:length(unique(dats[, factor.var.col])) +14
        ## names(pchs) <- unique(dats[, factor.var.col])
        ## these.pchs <- pchs[dats[, factor.var.col]]

        points(x=ys$x,
               y=ys$y,
               pch=16,
               col=these.cols,
               cex=1.2)
        points(x=ys$x,
               y=ys$y,
               col="black",
               cex=1.2)
    }
    plot(NA, xlim=range(dats[, xs], na.rm=TRUE),
         ylim=range(c(new.dd$plo, new.dd$phi,
                      dats[, y1], other.lim) ,
                    na.rm=TRUE),
         xlab="",
         ylab="",
         xaxt="n",
         yaxt="n",
         las=1, ...)

    if(plot.y){
        axis(2, pretty(range(c(other.lim, new.dd$plo, new.dd$phi,
                               dats[, y1]),
                             na.rm=TRUE),
                       n = 5),
             las=1, cex.axis=1.4)
        mtext(ylabel, 2, line=4, cex=1.4)
    }
    if(plot.x){
        axis(1, pretty(dats[,xs], 5), cex.axis=1.3)
    }
    plotting.loop()
}

