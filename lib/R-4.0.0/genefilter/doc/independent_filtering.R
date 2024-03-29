## ----knitr, echo=FALSE, results="hide"----------------------------------------
library("knitr")
opts_chunk$set(tidy=FALSE,dev="png",fig.show="hide",
               fig.width=4,fig.height=4.5,dpi=240,
               message=FALSE,error=FALSE,warning=FALSE)

## ----style, eval=TRUE, echo=FALSE, results="asis"--------------------------
BiocStyle:::latex()

## ----options,results='hide',echo=FALSE------------------------------------------------------------
options(digits=3, width=100)
library("pasilla") # make sure this is installed, since we need it in the next section

## ----libraries,results='hide'---------------------------------------------------------------------
library("pasilla")
data("pasillaGenes")

## ----DESeq1,results='hide'------------------------------------------------------------------------
library("DESeq")

## ----DESeq2,cache=TRUE,results='hide'-------------------------------------------------------------
cds  = estimateSizeFactors( pasillaGenes )
cds  = estimateDispersions( cds )
fit1 = fitNbinomGLMs( cds, count ~ type + condition )
fit0 = fitNbinomGLMs( cds, count ~ type  )

## ----DESeq3,cache=TRUE----------------------------------------------------------------------------
res = data.frame(
filterstat = rowMeans(counts(cds)),
pvalue    = nbinomGLMTest( fit1, fit0 ),
row.names = featureNames(cds) )

## ----headres--------------------------------------------------------------------------------------
dim(res)
head(res)

## ----pass,echo=FALSE,cache=TRUE-------------------------------------------------------------------
theta = 0.4
pass = with(res, filterstat > quantile(filterstat, theta))

## ----figscatterindepfilt--------------------------------------------------------------------------
with(res,
  plot(rank(filterstat)/length(filterstat), -log10(pvalue), pch=16, cex=0.45))

## ----figecdffilt----------------------------------------------------------------------------------
trsf = function(n) log10(n+1)
plot(ecdf(trsf(res$filterstat)), xlab=body(trsf), main="")

## ----badfilter1,cache=TRUE------------------------------------------------------------------------
badfilter = as.numeric(gsub("[+]*FBgn", "", rownames(res)))

## ----badfilter2,echo=FALSE------------------------------------------------------------------------
stopifnot(!any(is.na(badfilter)))

## ----figbadfilter---------------------------------------------------------------------------------
plot(rank(badfilter)/length(badfilter), -log10(res$pvalue), pch=16, cex=0.45)

## ----genefilter,results='hide'--------------------------------------------------------------------
library("genefilter")

## ----pBH1,cache=TRUE------------------------------------------------------------------------------
theta = seq(from=0, to=0.5, by=0.1)
pBH = filtered_p(filter=res$filterstat, test=res$pvalue, theta=theta, method="BH")

## ----pBH2-----------------------------------------------------------------------------------------
head(pBH)

## ----figrejection,fig.width=5.5,fig.height=5.5----------------------------------------------------
rejection_plot(pBH, at="sample",
               xlim=c(0, 0.5), ylim=c(0, 2000),
               xlab="FDR cutoff (Benjamini & Hochberg adjusted p-value)", main="")

## ----filtered_R1,cache=TRUE-----------------------------------------------------------------------
theta = seq(from=0, to=0.8, by=0.02)
rejBH = filtered_R(alpha=0.1, filter=res$filterstat, test=res$pvalue, theta=theta, method="BH")

## ----fignumreject,fig.width=5.5,fig.height=5.5----------------------------------------------------
plot(theta, rejBH, type="l",
     xlab=expression(theta), ylab="number of rejections")

## ----differentstats,cache=TRUE--------------------------------------------------------------------
filterChoices = data.frame(
  `mean`   = res$filterstat,
  `geneID` = badfilter,
  `min`    = rowMin(counts(cds)),
  `max`    = rowMax(counts(cds)),
  `sd`     = rowSds(counts(cds))
)
rejChoices = sapply(filterChoices, function(f)
  filtered_R(alpha=0.1, filter=f, test=res$pvalue, theta=theta, method="BH"))

## ----colours,results='hide'-----------------------------------------------------------------------
library("RColorBrewer")
myColours = brewer.pal(ncol(filterChoices), "Set1")

## ----figdifferentstats,fig.width=5.5,fig.height=5.5-----------------------------------------------
matplot(theta, rejChoices, type="l", lty=1, col=myColours, lwd=2,
        xlab=expression(theta), ylab="number of rejections")
legend("bottomleft", legend=colnames(filterChoices), fill=myColours)

## ----histindepfilt, fig.width=7, fig.height=5-----------------------------------------------------
h1 = hist(res$pvalue[!pass], breaks=50, plot=FALSE)
h2 = hist(res$pvalue[pass], breaks=50, plot=FALSE)
colori <- c(`do not pass`="khaki", `pass`="powderblue")

## ----fighistindepfilt, dev="pdf"------------------------------------------------------------------
barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
        col = colori, space = 0, main = "", ylab="frequency")
text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
     adj = c(0.5,1.7), xpd=NA)
legend("topright", fill=rev(colori), legend=rev(names(colori)))

## ----sortP, cache=TRUE----------------------------------------------------------------------------
resFilt = res[pass,]
orderInPlot = order(resFilt$pvalue)
showInPlot = (resFilt$pvalue[orderInPlot] <= 0.06)
alpha = 0.1

## ----sortedP, fig.width=4.5, fig.height=4.5-------------------------------------------------------
plot(seq(along=which(showInPlot)), resFilt$pvalue[orderInPlot][showInPlot],
     pch=".", xlab = expression(rank(p[i])), ylab=expression(p[i]))
abline(a=0, b=alpha/length(resFilt$pvalue), col="red3", lwd=2)

## ----doBH, echo=FALSE, results='hide'-------------------------------------------------------------
whichBH = which(resFilt$pvalue[orderInPlot] <= alpha*seq(along=resFilt$pvalue)/length(resFilt$pvalue))
## Test some assertions:
## - whichBH is a contiguous set of integers from 1 to length(whichBH)
## - the genes selected by this graphical method coincide with those
##   from p.adjust (i.e. padjFilt)
stopifnot(length(whichBH)>0,
          identical(whichBH, seq(along=whichBH)),
          resFilt$FDR[orderInPlot][ whichBH] <= alpha,
          resFilt$FDR[orderInPlot][-whichBH]  > alpha)

## ----SchwSpjot, echo=FALSE, results='hide'--------------------------------------------------------
j  = round(length(resFilt$pvalue)*c(1, .66))
px = (1-resFilt$pvalue[orderInPlot[j]])
py = ((length(resFilt$pvalue)-1):0)[j]
slope = diff(py)/diff(px)

## ----SchwederSpjotvoll, fig.width=4.5, fig.height=4.5---------------------------------------------
plot(1-resFilt$pvalue[orderInPlot],
     (length(resFilt$pvalue)-1):0, pch=".", xaxs="i", yaxs="i",
     xlab=expression(1-p[i]), ylab=expression(N(p[i])))
abline(a=0, slope, col="red3", lwd=2)
abline(h=slope)
text(x=0, y=slope, labels=paste(round(slope)), adj=c(-0.1, 1.3))

## ----sessionInfo, results='asis', echo=FALSE------------------------------------------------------
si = as.character( toLatex( sessionInfo() ) )
cat( si[ -grep( "Locale", si ) ], sep = "\n" )

