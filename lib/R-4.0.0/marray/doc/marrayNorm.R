### R code from vignette source 'marrayNorm.Rnw'

###################################################
### code chunk number 1: marrayNorm.Rnw:293-296
###################################################
library("marray", verbose=FALSE)
data(swirl)
maPlate(swirl)<-maCompPlate(swirl,n=384)


###################################################
### code chunk number 2: marrayNorm.Rnw:314-316
###################################################
swirl.norm <- maNorm(swirl, norm="p")
summary(swirl.norm)


###################################################
### code chunk number 3: marrayNorm.Rnw:340-342 (eval = FALSE)
###################################################
## swirl.norm1 <- maNorm(swirl, norm="p")
## swirl.norm2 <- maNormScale(swirl.norm1, norm="p")


###################################################
### code chunk number 4: maBoxplot1pre
###################################################
boxplot(swirl[,3], xvar="maPrintTip", yvar="maM", main="Swirl array 93: pre--normalization")


###################################################
### code chunk number 5: maBoxplot2pre
###################################################
boxplot(swirl, yvar="maM", main="Swirl arrays: pre--normalization")


###################################################
### code chunk number 6: maBoxplot1post
###################################################
boxplot(swirl.norm[,3], xvar="maPrintTip", yvar="maM", main="Swirl array 93: post--normalization")


###################################################
### code chunk number 7: maBoxplot2post
###################################################
boxplot(swirl.norm, yvar="maM", main="Swirl arrays: post--normalization")


###################################################
### code chunk number 8: maPlot1pre
###################################################
plot(swirl[,3], main="Swirl array 93: pre--normalization MA--plot")


###################################################
### code chunk number 9: maPlot1post
###################################################
plot(swirl.norm[,3], main="Swirl array 93: post--normalization MA--plot")


