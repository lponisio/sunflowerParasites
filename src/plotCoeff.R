source("src/misc.R")

## inverse logit
inv.logit <- function(a){
    exp(a)/(exp(a) + 1)
}

plotCoeffs <- function(mod, spec, name, ylabel,
                       adj1=0.01, binom){
    f.coeff <- function(){
        uniq.gensp <- unique(spec$GenusSpecies)
        uniq.gensp <- sort(uniq.gensp)
        cols <- brewer.pal(length(uniq.gensp), "Set3")
        names(cols) <- uniq.gensp
        coeffs.mod <- summary(mod)$coefficients
        means <- c(coeffs.mod[1,1],
                   coeffs.mod[1,1] + coeffs.mod[2:nrow(coeffs.mod), 1])

        ci.ub <- means +
            coeffs.mod[1:nrow(coeffs.mod), 2]

        ci.lb <- means -
            coeffs.mod[1:nrow(coeffs.mod), 2]

        if(binom){
            means <- inv.logit(means)
            ci.ub <- inv.logit(ci.ub)
            ci.lb <- inv.logit(ci.lb)
            }

        par(oma=c(10, 6, 0.5, 1),
            mar=c(0.5, 0, 2.5, 1), cex.axis=1.5)
        plot(x=1:nrow(coeffs.mod),
             y=means,
             col=cols,
             pch=16,
             ylim=range(0, ci.ub, means, na.rm=TRUE),
             ## xlim=c(0.5,3 .5),
             xlab='',
             xaxt='n',
             cex=1.5, las=1)
        arrows(x0=1:nrow(coeffs.mod),
               y0= ci.lb,
               y1=ci.ub,
               angle=90,
               length=0, code=3,
               col=cols,
               lwd=1.5)
        text(1:nrow(coeffs.mod),
             par('usr')[3] - adj1,
             srt = 45, adj = 1,
             labels = uniq.gensp,
             xpd = NA,
             cex=1)
        points(x=1:nrow(coeffs.mod),
               y=means,
               pch=1,
               cex=1.5)
        mtext(ylabel, 2,
              line=4.5 , cex=1.5)

    }
    f.path <- 'figures/'
    pdf.f(f.coeff,
          file= file.path(f.path, sprintf("%s.pdf", name)),
          width=6, height=6)
}

