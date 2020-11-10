### R code from vignette source 'marrayPlots.Rnw'

###################################################
### code chunk number 1: marrayPlots.Rnw:105-108
###################################################
library(marray)
data(swirl) 
maPlate(swirl)<-maCompPlate(swirl,n=384) 


###################################################
### code chunk number 2: marrayPlots.Rnw:128-131
###################################################
Gcol<- maPalette(low="white", high="green",k=50) 
Rcol<- maPalette(low="white", high="red", k=50) 
RGcol<-maPalette(low="green", high="red", k=50) 


###################################################
### code chunk number 3: maImageGb
###################################################
tmp<-image(swirl[,3], xvar="maGb", subset=TRUE, col=Gcol,contours=FALSE, bar=FALSE) 


###################################################
### code chunk number 4: maImageRb
###################################################
tmp<-image(swirl[,3], xvar="maRb", subset=TRUE, col=Rcol, contours=FALSE, bar=FALSE) 


###################################################
### code chunk number 5: maImageMraw1
###################################################
tmp<-image(swirl[,3], xvar="maM", bar=FALSE, main="Swirl array 93: image of pre--normalization M") 


###################################################
### code chunk number 6: maImageMraw2
###################################################
tmp<-image(swirl[,3], xvar="maM", subset=maTop(maM(swirl[,3]), h=0.10,
l=0.10), col=RGcol, contours=FALSE, bar=FALSE,main="Swirl array 93:
image of pre--normalization M for  10 %  tails")  


###################################################
### code chunk number 7: maImageSpotCol
###################################################
tmp<- image(swirl[,3], xvar="maSpotCol", bar=FALSE) 


###################################################
### code chunk number 8: maImagePrintTip
###################################################
tmp<- image(swirl[,3], xvar="maPrintTip", bar=FALSE) 


###################################################
### code chunk number 9: maImageControls
###################################################
tmp<- image(swirl[,3], xvar="maControls",col=heat.colors(10),bar=FALSE) 


###################################################
### code chunk number 10: maImagePlate
###################################################
tmp<- image(swirl[,3], xvar="maPlate",bar=FALSE) 


###################################################
### code chunk number 11: maBoxplot1pre
###################################################
boxplot(swirl[,3], xvar="maPrintTip", yvar="maM", main="Swirl array 93: pre--normalization") 


###################################################
### code chunk number 12: maBoxplot2pre
###################################################
boxplot(swirl, yvar="maM", main="Swirl arrays: pre--normalization") 


###################################################
### code chunk number 13: marrayPlots.Rnw:343-344
###################################################
swirl.norm <- maNorm(swirl, norm="p")


###################################################
### code chunk number 14: maBoxplot1post
###################################################
boxplot(swirl.norm[,3], xvar="maPrintTip", yvar="maM",
	main="Swirl array 93: post--normalization") 


###################################################
### code chunk number 15: maBoxplot2post
###################################################
boxplot(swirl.norm, yvar="maM", col="green", main="Swirl arrays: post--normalization") 


###################################################
### code chunk number 16: maPlot1pre
###################################################
defs<-maDefaultPar(swirl[,3],x="maA",y="maM",z="maPrintTip")

# Function for plotting the legend
legend.func<-do.call("maLegendLines",defs$def.legend)

# Function for performing and plotting lowess fits
lines.func<-do.call("maLowessLines",c(list(TRUE,f=0.3),defs$def.lines))

plot(swirl[,3], xvar="maA", yvar="maM", zvar="maPrintTip",
		      lines.func,
		      text.func=maText(),
		      legend.func,
		      main="Swirl array 93: pre--normalization MA--plot") 


###################################################
### code chunk number 17: maPlot1post
###################################################
plot(swirl.norm[,3], xvar="maA", yvar="maM", zvar="maPrintTip",
		      lines.func,
		      text.func=maText(),
		      legend.func,
		      main="Swirl array 93: post--normalization MA--plot") 


