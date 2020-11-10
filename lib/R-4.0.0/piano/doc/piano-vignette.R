### R code from vignette source 'piano-vignette.Rnw'

###################################################
### code chunk number 1: piano-vignette.Rnw:55-56
###################################################
library(piano)


###################################################
### code chunk number 2: piano-vignette.Rnw:129-130
###################################################
library(piano)


###################################################
### code chunk number 3: piano-vignette.Rnw:179-192
###################################################
# An example setup:
oxygen <- c("aerobic","anaerobic")
limitation <- c("Clim","Nlim")
mySetup <- cbind(oxygen[c(1,1,1,2,2,2,1,1,1,2,2,2)],
               limitation[c(1,1,1,1,1,1,2,2,2,2,2,2)])

# The rownames correspond to the CEL-file names (CAE1.CEL etc):
rownames(mySetup) <- c("CAE1","CAE2","CAE3","CAN1","CAN2","CAN3",
                     "NAE1","NAE2","NAE3","NAN1","NAN2","NAN3")
colnames(mySetup) <- c("oxygen","limitation")

# The final setup object can look like this:
mySetup


###################################################
### code chunk number 4: piano-vignette.Rnw:219-226
###################################################
# Get path to example data and setup files:
dataPath <- system.file("extdata", package="piano")

# Load pre-normalized data:
myArrayData <- loadMAdata(datadir=dataPath, 
                          dataNorm="norm_data.txt.gz", 
                          platform="yeast2")


###################################################
### code chunk number 5: piano-vignette.Rnw:231-232
###################################################
myArrayData


###################################################
### code chunk number 6: piano-vignette.Rnw:237-241
###################################################
# Check the setup:
myArrayData$setup
# Check the annotation (top 10 rows):
myArrayData$annotation[1:10,]


###################################################
### code chunk number 7: piano-vignette.Rnw:250-251
###################################################
runQC(myArrayData)


###################################################
### code chunk number 8: piano-vignette.Rnw:256-262
###################################################
# To only run the PCA:
runQC(myArrayData, rnaDeg=FALSE, nuseRle=FALSE, hist=FALSE,
      boxplot=FALSE, pca=TRUE)
# Additionally, for the PCA you can specify other colors:
runQC(myArrayData, rnaDeg=FALSE, nuseRle=FALSE, hist=FALSE,
      boxplot=FALSE, pca=TRUE, colors=c("cyan","orange"))


###################################################
### code chunk number 9: piano-vignette.Rnw:279-280
###################################################
extractFactors(myArrayData)


###################################################
### code chunk number 10: piano-vignette.Rnw:285-287
###################################################
pfc <- diffExp(myArrayData, contrasts=c("aerobic_Clim - anaerobic_Clim", 
                                        "aerobic_Nlim - anaerobic_Nlim"))


###################################################
### code chunk number 11: piano-vignette.Rnw:297-303
###################################################
# Sort genes in "aerobic_Clim - anaerobic_Clim" according to adjusted p-value:
ii <- sort(pfc$resTable[[1]]$adj.P.Val, index.return=TRUE)$ix
pfc$resTable[[1]][ii[1:5],]
# Sort genes in "aerobic_Nlim - anaerobic_Nlim" according to adjusted p-value:
ii <- sort(pfc$resTable[[2]]$adj.P.Val, index.return=TRUE)$ix
pfc$resTable[[2]][ii[1:5],]


###################################################
### code chunk number 12: piano-vignette.Rnw:332-336
###################################################
# Get p-values from the aerobic_Clim vs anaerobic_Clim comparison:
myPval <- pfc$pValues["aerobic_Clim - anaerobic_Clim"]
# Display the first values and gene IDs:
head(myPval)


###################################################
### code chunk number 13: piano-vignette.Rnw:343-353
###################################################
# Custom gene to gene set mapping:
genes2genesets <- cbind(paste("gene",c("A","A","A","B","B","C","C","C","D"),sep=""),
                        paste("set",c(1,2,3,1,3,2,3,4,4),sep=""))
genes2genesets
# Load into correct format:
myGsc <- loadGSC(genes2genesets)
# View summary:
myGsc
# View all gene sets:
myGsc$gsc


###################################################
### code chunk number 14: piano-vignette.Rnw:358-361
###################################################
myStats <- c(-1.5,-0.5,1,2)
names(myStats) <- paste("gene",c("A","B","C","D"),sep="")
myStats


###################################################
### code chunk number 15: piano-vignette.Rnw:395-396
###################################################
gsaRes <- runGSA(myStats, gsc=myGsc)


###################################################
### code chunk number 16: piano-vignette.Rnw:402-403
###################################################
gsaRes


###################################################
### code chunk number 17: piano-vignette.Rnw:406-407
###################################################
names(gsaRes)


###################################################
### code chunk number 18: piano-vignette.Rnw:414-415
###################################################
c(-4, -3, 2, 3.5)


###################################################
### code chunk number 19: piano-vignette.Rnw:418-419
###################################################
mean(c(-4, -3, 2.5, 4.5))


###################################################
### code chunk number 20: piano-vignette.Rnw:422-423
###################################################
mean(abs(c(-4, -3, 2.5, 4.5)))


###################################################
### code chunk number 21: piano-vignette.Rnw:426-428
###################################################
mean(abs(c(-4, -3)))
mean(abs(c(2.5, 4.5)))


###################################################
### code chunk number 22: piano-vignette.Rnw:433-436
###################################################
tmp<-c(0.0001,0.0001,0.9999,0.0001,0.005)
names(tmp) <- c("p (non.dir)","p (dist.dir.up)","p (dist.dir.dn)","p (mix.dir.up)","p (mix.dir.dn)")
tmp


###################################################
### code chunk number 23: piano-vignette.Rnw:439-442
###################################################
tmp<-c(100,95,5)
names(tmp)<-c("Genes (tot)","Genes (up)","Genes (down)")
tmp


###################################################
### code chunk number 24: piano-vignette.Rnw:473-479
###################################################
# Get the p-values for the contrast aerobic_Clim - anaerobic_Clim
myPval <- pfc$pValues["aerobic_Clim - anaerobic_Clim"]
head(myPval)
# Get the fold-changes for the contrast aerobic_Clim - anaerobic_Clim
myFC <- pfc$foldChanges["aerobic_Clim - anaerobic_Clim"]
head(myFC)


###################################################
### code chunk number 25: piano-vignette.Rnw:482-494 (eval = FALSE)
###################################################
## library("biomaRt")
## # Select ensembl database and S. cerevisiae dataset:
## ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset="scerevisiae_gene_ensembl" ,host="www.ensembl.org")
## # Map Yeast 2.0 microarray probeset IDs to GO:
## mapGO <- getBM(attributes=c('affy_yeast_2', 'name_1006'), 
##                #filters = 'affy_yeast_2', 
##                #values = rownames(myPval), 
##                mart = ensembl)
## # Remove blanks ("")
## mapGO <- mapGO[mapGO[,2]!="",]
## # Check the 10 first rows to see what we got:
## mapGO[1:10,]


###################################################
### code chunk number 26: piano-vignette.Rnw:497-498 (eval = FALSE)
###################################################
## myGsc <- loadGSC(mapGO)


###################################################
### code chunk number 27: piano-vignette.Rnw:535-536 (eval = FALSE)
###################################################
## GSAsummaryTable(gsaRes, save=TRUE, file="gsaResTab.xls")


###################################################
### code chunk number 28: piano-vignette.Rnw:569-570 (eval = FALSE)
###################################################
## nw <- networkPlot(gsaRes,class="non")


###################################################
### code chunk number 29: piano-vignette.Rnw:579-580 (eval = FALSE)
###################################################
## nw$geneSets


###################################################
### code chunk number 30: piano-vignette.Rnw:586-588 (eval = FALSE)
###################################################
## # An example usage:
## myGsc <- loadGSC("myModel.xml")


###################################################
### code chunk number 31: piano-vignette.Rnw:591-600
###################################################
# To load iTO977:
metMap <- system.file("extdata", "probe2metabolites_iTO977.txt.gz", 
                      package="piano")
# To load iIN800:
metMap <- system.file("extdata", "probe2metabolites_iIN800.txt.gz", 
                      package="piano")
# Convert into piano format:
metMap <- read.delim(metMap)
myGsc <- loadGSC(metMap)


###################################################
### code chunk number 32: piano-vignette.Rnw:603-607
###################################################
gsaRes <- runGSA(myPval,myFC,gsc=myGsc,
                 geneSetStat="reporter",
                 signifMethod="nullDist", nPerm=1000,
                 gsSizeLim=c(5,100))


###################################################
### code chunk number 33: piano-vignette.Rnw:609-610
###################################################
gsaRes


###################################################
### code chunk number 34: piano-vignette.Rnw:615-617
###################################################
nw <- networkPlot(gsaRes, class="distinct", direction="both",
                  significance=0.005, label="numbers")


###################################################
### code chunk number 35: piano-vignette.Rnw:626-627
###################################################
nw$geneSets


###################################################
### code chunk number 36: piano-vignette.Rnw:631-635
###################################################
gsaResTab <- GSAsummaryTable(gsaRes)
# Which columns contain p-values:
grep("p \\(",colnames(gsaResTab),value=TRUE)
grep("p \\(",colnames(gsaResTab))


###################################################
### code chunk number 37: piano-vignette.Rnw:638-640
###################################################
ii <- which(gsaResTab[,10]<0.0001)
gsaResTab$Name[ii]


###################################################
### code chunk number 38: piano-vignette.Rnw:643-650
###################################################
# Get minimum p-value for each gene set:
minPval <- apply(gsaResTab[,c(4,7,10,14,18)],1,min,na.rm=TRUE)
# Select significant gene sets:
ii <- which(minPval<0.0001)
gsaResTabSign <- gsaResTab[ii,c(1,4,7,10,14,18)]
# Look at the first 10 gene sets:
gsaResTabSign[1:10,]


###################################################
### code chunk number 39: piano-vignette.Rnw:666-668 (eval = FALSE)
###################################################
## myTval <- pfc$resTable[["aerobic_Clim - anaerobic_Clim"]]$t
## names(myTval) <- pfc$resTable[["aerobic_Clim - anaerobic_Clim"]]$ProbesetID


###################################################
### code chunk number 40: piano-vignette.Rnw:671-686 (eval = FALSE)
###################################################
## myGsc <- loadGSC(mapGO)
## gsaRes1 <- runGSA(myTval,geneSetStat="mean",gsc=myGsc,
##                   nPerm=1000,gsSizeLim=c(10,800))
## gsaRes2 <- runGSA(myTval,geneSetStat="median",gsc=myGsc,
##                   nPerm=1000,gsSizeLim=c(10,800))
## gsaRes3 <- runGSA(myTval,geneSetStat="sum",gsc=myGsc,
##                   nPerm=1000,gsSizeLim=c(10,800))
## gsaRes4 <- runGSA(myTval,geneSetStat="maxmean",gsc=myGsc,
##                   nPerm=1000,gsSizeLim=c(10,800))
## gsaRes5 <- runGSA(myPval,myFC,geneSetStat="fisher",gsc=myGsc,
##                   nPerm=1000,gsSizeLim=c(10,800))
## gsaRes6 <- runGSA(myPval,myFC,geneSetStat="stouffer",gsc=myGsc,
##                   nPerm=1000,gsSizeLim=c(10,800))
## gsaRes7 <- runGSA(myPval,myFC,geneSetStat="tailStrength",gsc=myGsc,
##                   nPerm=1000,gsSizeLim=c(10,800))


###################################################
### code chunk number 41: piano-vignette.Rnw:689-693 (eval = FALSE)
###################################################
## resList <- list(gsaRes1,gsaRes2,gsaRes3,gsaRes4,gsaRes5,gsaRes6,gsaRes7)
## names(resList) <- c("mean","median","sum","maxmean","fisher",
##                     "stouffer","tailStrength")
## ch <- consensusHeatmap(resList,cutoff=30,method="mean")


###################################################
### code chunk number 42: piano-vignette.Rnw:696-697 (eval = FALSE)
###################################################
## ch$pMat


