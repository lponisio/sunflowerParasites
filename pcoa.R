## setwd("~/Dropbox/sunflower")
rm(list=ls())
setwd("analysis/parasiteCommunity")
library(vegan)
source("src/plotPcoa.R")

load('../../data/spec.Rdata')

parasites <- c("Aspergillus",
               "Ascosphaera", "Apicystis", "Crithidiaspp", "Cexpoeki",
               "Cbombi",  "Nbombi", "Nceranae")

pre.comm.mat <- spec[,c("UniqueID", parasites, "GenusSpecies")]

## find rows where we did not screen for pathogens
not.path.screen <- apply(pre.comm.mat[, parasites], 1,
                         function(x) any(is.na(x)))

pre.comm.mat <- pre.comm.mat[!not.path.screen,]

## converting dataframe into a community matrix with only "sites" and
## "species"
comm.mat <- pre.comm.mat
GenSp <- pre.comm.mat$GenusSpecies
rownames(comm.mat)  <-  comm.mat$UniqueID
comm.mat$UniqueID <- NULL
comm.mat$GenusSpecies <- NULL

## make distance matrix, need a method that is binary i.e. incidence
## based, because our parasite data is 0/1

dist.mat <- vegdist(comm.mat, method= "jaccard",
                    na.rm=TRUE, diag=TRUE)

beta.disper.result <- betadisper(dist.mat, GenSp, type="centroid")


## Permutation test for F and simulate missing values to compare the
## differences in the variances of the community composition of
## parasites between bee species

## perm.test <- permutest(beta.disper.result,
##           control = permControl(nperm = 100),
##           pairwise = TRUE)

## ## Draw a boxplot of the distances to centroid for each group
## boxplot(beta.disper.result)

uniq.gensp <- unique(GenSp)
uniq.gensp <- sort(uniq.gensp)
cols <- add.alpha(brewer.pal(length(uniq.gensp), "Set3"),
                  alpha=0.8)
names(cols) <- uniq.gensp

plotCommDist(dist.mat, spec, cols)
plotBetaDiv(pcoa.res=beta.disper.result, cols)
