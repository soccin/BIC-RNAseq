### R code from vignette source 'marrayInput.Rnw'

###################################################
### code chunk number 1: marrayInput.Rnw:120-122
###################################################
library(marray)
data(swirl)


###################################################
### code chunk number 2: marrayInput.Rnw:153-155
###################################################
datadir <- system.file("swirldata", package="marray")
dir(datadir)


###################################################
### code chunk number 3: marrayInput.Rnw:183-185
###################################################
swirlTargets <- read.marrayInfo(file.path(datadir, "SwirlSample.txt"))
summary(swirlTargets)


###################################################
### code chunk number 4: marrayInput.Rnw:203-205
###################################################
galinfo <- read.Galfile("fish.gal", path=datadir) 
names(galinfo)


###################################################
### code chunk number 5: marrayInput.Rnw:245-248
###################################################
swirl.gnames <- read.marrayInfo(file.path(datadir, "fish.gal"),
info.id=4:5, labels=5, skip=21) 
summary(swirl.gnames) 


###################################################
### code chunk number 6: marrayInput.Rnw:255-262
###################################################
swirl.layout <- read.marrayLayout(fname=file.path(datadir, "fish.gal"),
                                  ngr=4, ngc=4, nsr=22, nsc=24,
                                  skip=21,ctl.col=4)
ctl<-rep("Control",maNspots(swirl.layout))
ctl[maControls(swirl.layout)!="control"]  <- "probes"
maControls(swirl.layout)<-factor(ctl)
summary(swirl.layout)


###################################################
### code chunk number 7: marrayInput.Rnw:278-283
###################################################
mraw <- read.Spot(path=datadir, 
                  layout=galinfo$layout, 
                  gnames=galinfo$gnames, 
                  target=swirlTargets)
summary(mraw)


