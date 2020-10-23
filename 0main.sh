#!/usr/bin/env bash

## prep data
## only needed for original analysis, the saved .Rdata files should
## all in in github

Rscript dataPrep/dataPrep.R

##***************************************************************
## Wild Bee abundance/richness models
##***************************************************************

## the entire wild bee community
Rscript analysis/parasiteCommunity/1BeeMods.R 
Rscript analysis/parasiteCommunity/2ParMods.R 
Rscript analysis/parasiteCommunity/3HBMods.R 
Rscript analysis/parasiteCommunity/4BeeParPlotting.R 

