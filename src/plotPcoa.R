library(RColorBrewer)


pdf.f <- function(f, file, ...) {
    cat(sprintf('Writing %s\n', file))
    pdf(file, ...)
    on.exit(dev.off())
    f()
}

## add transparency to named colors
add.alpha <- function(col, alpha=0.2){
    apply(sapply(col, col2rgb)/255, 2,
          function(x)
              rgb(x[1], x[2], x[3],
                  alpha=alpha))
}


plotCommDist  <- function(dist.mat, spec, cols){

    f.pcoa <- function(){
        dist.mat <- as.matrix(dist.mat)
        dist.mat.na <- apply(dist.mat, 1, function(x) any(is.na(x)))
        dist.mat <- dist.mat[!dist.mat.na, !dist.mat.na]
        pcoa.comm <- cmdscale(dist.mat)
        GenusSpecies <-  as.vector(spec$GenusSpecies[match(
                                            rownames(pcoa.comm),
                                            spec$UniqueID)])
        pcoa.mod <- adonis(dist.mat~GenusSpecies)
        plot(NA, asp=1,  cex=1.5,
             ylim=range(pcoa.comm[,2]),
             xlim=range(pcoa.comm[,1]),
             xlab='',
             ylab='',
             xaxt='n',
             yaxt='n',
             cex.lab=1.5)
        pcoa.comm.jit <- jitter(jitter(pcoa.comm))
        for(species in uniq.gensp){
            ## all points sitting on top of eachother so triple jitter
            points(pcoa.comm.jit[GenusSpecies == species,],
                   col=cols[species], pch=16, cex=1.5)
            points(pcoa.comm.jit[GenusSpecies == species,],
                   col="black", pch=1, cex=1.5)
        }
        ordihull(pcoa.comm, GenusSpecies)

        legend("topleft", legend=uniq.gensp,
               bty="n", cex=0.8, col=cols[uniq.gensp],
               pch=16)
        legend("topleft", legend=uniq.gensp,
               bty="n", cex=0.8, col="black", pch=1)

        mtext('PCoA1', 1, line=2, cex=1.5)
        mtext('PCoA2', 2, line=2, cex=1.5)
        return(pcoa.mod)
    }
    ## function for plotting PcoA axes
    path <- 'figures'
    pdf.f(f.pcoa, file= file.path(path,
                                  "pcoa.pdf"),
          width=7, height=7)
}




plotBetaDiv  <- function(pcoa.res, cols){
    f.plotBetaDiv <- function(){
        par(mar=c(12,5,0,1))
        mp1 <- boxplot(pcoa.res,
                       las=2, col=cols, main="")
        mtext("Beta-diversity", 3, line=0.5, cex=1.2)
    }
    ## function for plotting Pcoa dispersion, i.e., beta diversity
    path <- 'figures'
    pdf.f(f.plotBetaDiv, file= file.path(path,
                                         "betadiv.pdf"),
          width=7, height=7)
}
