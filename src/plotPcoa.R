source("src/misc.R")

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

        legend("topright", legend=uniq.gensp,
               bty="n", cex=0.6, col=cols[uniq.gensp],
               pch=16)
        legend("topright", legend=uniq.gensp,
               bty="n", cex=0.6, col="black", pch=1)

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
                       las=2, col=cols, main="", xlab="")
        mtext("Beta-diversity", 3, line=0.5, cex=1.2)
    }
    ## function for plotting Pcoa dispersion, i.e., beta diversity
    path <- 'figures'
    pdf.f(f.plotBetaDiv, file= file.path(path,
                                         "betadiv.pdf"),
          width=7, height=7)
}



plotCommDistbyGroup  <- function(dist.mat, comm,
                                 comm.groups,
                                 species.type){

    f.pcoa <- function(){
        ## dist.mat.na <- apply(dist.mat, 1, function(x) any(is.na(x)))
        ## dist.mat <- dist.mat[!dist.mat.na, !dist.mat.na]
        pcoa.comm <- cmdscale(dist.mat)
        groups <- comm[[comm.groups]]
        pcoa.mod <- adonis(dist.mat~groups)
        plot(NA, asp=1,  cex=1.5,
             ylim=range(pcoa.comm[,2]),
             xlim=range(pcoa.comm[,1]),
             xlab='',
             ylab='',
             xaxt='n',
             yaxt='n',
             cex.lab=1.5)

        cols <- viridis(length(unique(groups)))
        names(cols) <- unique(groups)
        for(group in unique(groups)){
            ## all points sitting on top of eachother so triple jitter
            points(pcoa.comm[groups == group,],
                   col=cols[group], pch=16, cex=1.5)
            points(pcoa.comm[groups == group,],
                   col="black", pch=1, cex=1.5)
        }
        ordihull(pcoa.comm, groups)

        legend("topleft", legend= unique(groups),
               bty="n", col=cols[unique(groups)],
               pch=16, cex=1.2)
        legend("topleft", legend= unique(groups),
               bty="n", col="black", pch=1,
               cex=1.2)

        mtext('PCoA1', 1, line=2, cex=1.5)
        mtext('PCoA2', 2, line=2, cex=1.5)
        return(pcoa.mod)
    }
    ## function for plotting PcoA axes
    path <- 'figures/pcas'
    pdf.f(f.pcoa,
          file= file.path(path,
                          sprintf("%s_%s_pcoa.pdf",
                                  species.type, comm.groups)),
          width=7, height=7)
}



plotParasiteMap <- function(){
    colfunc <- colorRampPalette(c("white", "red"))
    par(oma=c(6,4,3,2), mar=c(1,2,2,1),
        mgp=c(1.5,0.5,0))
    heatmap.2(parasite.comms[[parasite]],
              trace="none",
              col=colfunc)
    mtext(parasite, 3, line=2.5)
}
