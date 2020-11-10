### R code from vignette source 'marrayClasses.Rnw'

###################################################
### code chunk number 1: marrayClasses.Rnw:65-67
###################################################
library("marray")
data(swirl)


###################################################
### code chunk number 2: R
###################################################
getClassDef("marrayLayout")


###################################################
### code chunk number 3: R
###################################################
getClassDef("marrayInfo")


###################################################
### code chunk number 4: R
###################################################
getClassDef("marrayRaw")


###################################################
### code chunk number 5: R
###################################################
getClassDef("marrayNorm")


###################################################
### code chunk number 6: marrayClasses.Rnw:258-266
###################################################
zebra.RG<-as.data.frame(cbind(c("swirl","WT","swirl","WT"),
c("WT","swirl","WT","swirl")))
dimnames(zebra.RG)[[2]]<-c("Cy3","Cy5")
zebra.samples<-new("marrayInfo",
		    maLabels=paste("Swirl array ",1:4,sep=""), 
		    maInfo=zebra.RG,
		    maNotes="Description of targets for Swirl experiment")
zebra.samples


###################################################
### code chunk number 7: marrayClasses.Rnw:283-285
###################################################
L<-slot(swirl, "maLayout")
L@maNgr


###################################################
### code chunk number 8: marrayClasses.Rnw:290-292
###################################################
slotNames("marrayLayout")
slotNames(swirl)


###################################################
### code chunk number 9: marrayClasses.Rnw:303-304
###################################################
validObject(maLayout(swirl), test=TRUE)


###################################################
### code chunk number 10: marrayClasses.Rnw:314-316 (eval = FALSE)
###################################################
## showMethods(classes="marrayLayout")
## showMethods("show",classes="marrayLayout") 


###################################################
### code chunk number 11: marrayClasses.Rnw:324-325
###################################################
showMethods("summary")


###################################################
### code chunk number 12: marrayClasses.Rnw:331-332
###################################################
summary(swirl)


###################################################
### code chunk number 13: marrayClasses.Rnw:346-347
###################################################
swirl[1:100,2:3]


###################################################
### code chunk number 14: marrayClasses.Rnw:376-382
###################################################
swirl.layout<-maLayout(swirl)
maNspots(swirl)
maNspots(swirl.layout)
maNgr(swirl)
maNgc(swirl.layout)
maPrintTip(swirl[1:10,3])


###################################################
### code chunk number 15: marrayClasses.Rnw:395-398
###################################################
maNotes(swirl.layout)
maNotes(swirl.layout)<- "New value"
maNotes(swirl.layout)


###################################################
### code chunk number 16: marrayClasses.Rnw:403-406
###################################################
L<-new("marrayLayout")
L
maNgr(L)<-4


###################################################
### code chunk number 17: marrayClasses.Rnw:422-423
###################################################
swirl.norm<-as(swirl, "marrayNorm")    


###################################################
### code chunk number 18: marrayClasses.Rnw:438-439
###################################################
maPlate(swirl)<-maCompPlate(swirl,n=384)


