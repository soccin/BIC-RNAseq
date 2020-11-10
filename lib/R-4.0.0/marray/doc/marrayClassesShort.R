### R code from vignette source 'marrayClassesShort.Rnw'

###################################################
### code chunk number 1: marrayClassesShort.Rnw:53-55
###################################################
library(marray)
data(swirl)


###################################################
### code chunk number 2: marrayClassesShort.Rnw:98-100
###################################################
slotNames("marrayLayout")
slotNames(swirl)


###################################################
### code chunk number 3: marrayClassesShort.Rnw:111-118
###################################################
zebra.RG<-as.data.frame(cbind(c("swirl","WT","swirl","WT"), c("WT","swirl","WT","swirl"))) 
dimnames(zebra.RG)[[2]]<-c("Cy3","Cy5")
zebra.samples<-new("marrayInfo",
		    maLabels=paste("Swirl array ",1:4,sep=""), 
		    maInfo=zebra.RG,
		    maNotes="Description of targets for Swirl experiment")
zebra.samples


###################################################
### code chunk number 4: marrayClassesShort.Rnw:134-136
###################################################
L<-slot(swirl, "maLayout")
L@maNgr


###################################################
### code chunk number 5: marrayClassesShort.Rnw:146-148 (eval = FALSE)
###################################################
## showMethods(classes="marrayLayout")
## showMethods("summary",classes="marrayLayout") 


###################################################
### code chunk number 6: marrayClassesShort.Rnw:159-160
###################################################
showMethods("print")


###################################################
### code chunk number 7: marrayClassesShort.Rnw:166-167
###################################################
summary(swirl)


###################################################
### code chunk number 8: marrayClassesShort.Rnw:181-182
###################################################
swirl[1:100,2:3]


###################################################
### code chunk number 9: marrayClassesShort.Rnw:210-216
###################################################
swirl.layout<-maLayout(swirl)
maNspots(swirl)
maNspots(swirl.layout)
maNgr(swirl)
maNgc(swirl.layout)
maPrintTip(swirl[1:10,3])


###################################################
### code chunk number 10: marrayClassesShort.Rnw:228-231
###################################################
maNotes(swirl.layout)
maNotes(swirl.layout)<- "New value"
maNotes(swirl.layout)


###################################################
### code chunk number 11: marrayClassesShort.Rnw:236-239
###################################################
L<-new("marrayLayout")
L
maNgr(L)<-4


###################################################
### code chunk number 12: marrayClassesShort.Rnw:253-254
###################################################
swirl.norm<-as(swirl, "marrayNorm")    


###################################################
### code chunk number 13: marrayClassesShort.Rnw:269-270
###################################################
maPlate(swirl)<-maCompPlate(swirl,n=384)


