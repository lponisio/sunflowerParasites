
inv.logit <- function(x){
    exp(x)/(1+exp(x))
}


sumMSdredge <- function(res){
    res <- summary(res)
    mmi <- as.data.frame(res$coefmat.subset)

    intercept <- mmi[ "(Intercept)", "Estimate"]
    next.levs <- c(0, mmi[, "Estimate"][-1])
    SEs <- mmi[, "Std. Error"]

    mmi$P1 <- inv.logit(intercept + next.levs)
    mmi$P1.ci.ub <-     mmi$P1 +  qnorm(0.975)*SEs
    mmi$P1.ci.lb <-     mmi$P1 -  qnorm(0.975)*SEs

    mmi$OR  <- exp(mmi[, "Estimate"])
    mmi$OR.ci.ub <-   exp(mmi[, "Estimate"] +  qnorm(0.975)*SEs)
    mmi$OR.ci.lb <-   exp(mmi[, "Estimate"] -  qnorm(0.975)*SEs)

    mmi$ORdelta  <- (mmi$OR-1)*100
    mmi$ORdelta.ci.ub <-    (mmi$OR.ci.ub -1)*100
    mmi$ORdelta.ci.lb <-   (mmi$OR.ci.lb -1)*100
    mmi <- round(mmi, 3)
    return(mmi)
}
