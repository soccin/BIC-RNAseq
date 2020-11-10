### R code from vignette source 'marray.Rnw'

###################################################
### code chunk number 1: marray.Rnw:86-88
###################################################
library(marray)
dir(system.file("swirldata", package="marray")) 


###################################################
### code chunk number 2: marray.Rnw:95-97
###################################################
datadir <- system.file("swirldata", package="marray")
swirlTargets <- read.marrayInfo(file.path(datadir, "SwirlSample.txt"))


###################################################
### code chunk number 3: marray.Rnw:104-105
###################################################
mraw <- read.Spot(targets = swirlTargets, path=datadir)


###################################################
### code chunk number 4: marray.Rnw:118-121
###################################################
galinfo <- read.Galfile("fish.gal", path=datadir)
mraw@maLayout <- galinfo$layout
mraw@maGnames <- galinfo$gnames


