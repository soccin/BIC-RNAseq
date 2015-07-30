scale.factor.quantile<-function(counts.v,percentile="75%")
{
	q=quantile(counts.v,probs=seq(0,1,0.05))
	q=q[grep(percentile,names(q))]
	return(sum(counts.v[which(counts.v<=q)]))
} 
 
#######
make.cds<-function(counts.dat,conds,count.cut=10,libsizeQ=F,percentile="100%",fitType="parametric",method="per-condition", sharingMode="maximum")
{
### libsizeQ=F by default scale factors are estimated using DESeq procedure ###

### remove genes/bacteria/species with total counts.dat less than count.cut in all samples... ###
sum.bact=apply(counts.dat,1,sum)
counts.tmp.dat=counts.dat[which(sum.bact>=count.cut),]
cds = newCountDataSet(counts.tmp.dat, conds)
#cds = estimateSizeFactors(cds)
if(!libsizeQ){
    cds = estimateSizeFactors(cds)
}
if(libsizeQ | length(which(is.na(sizeFactors(cds))))>0){
	       print("estimating sizeFactors using library size")
	       if(percentile=="100%")
	     		libsizes=apply(counts.tmp.dat,2,sum)
				
				if(percentile!="100%")
						libsizes=apply(counts.tmp.dat,2,scale.factor.quantile,percentile=percentile)

 
     		libsizes=libsizes/median(libsizes)
     		sizeFactors(cds)=libsizes
}
cds=estimateDispersions(cds,fit=fitType,method=method,sharingMode=sharingMode)

return(cds)
}

######
run.nbinomTest<-function(cds, conds, condA, condB, q.cut=0.05, fc.cut=2, count.cut=10,zeroaddQ=F,libsizeQ=F,percentile="100%",fitType="parametric",gns=c())
{
ng=NULL
ans=NULL
res=NULL

#cat(c("length(gns)",length(gns),"\n"))

###start###
res = nbinomTest(cds, condA, condB)


save(res,file="resTest.Rdata")


#cat(c("dim(res)",dim(res),"\n"))

rownames(res)=res[,1]
if (length(gns)>0){
    gns = gns[rownames(res)]
    res = cbind(res,gns)
    colnames(res)[length(colnames(res))] = "GeneSymbol"
}
#cat(c("colnames(res)",colnames(res),"\n"))
#cat(c("dim(res)",dim(res),"\n"))
res=res[which(!is.nan(res[,6])),]

ng = rownames(res)[which(res$padj <= q.cut & abs(res$log2FoldChange)>=log2(fc.cut))]
norm.factors=sizeFactors(cds)
counts.scaled=counts(cds,norm=T)

### just adding genes that have zero counts in one condition and greater than count.cut/mean(sizeFactors) in another condition BUT not significant p-value as DESeq get it wrong here... ###
if(zeroaddQ)
{
#print("adding zero counts with non-significant p-values? ")
jj1=rownames(res)[which(res[,"baseMeanA"]==0 & res[,"baseMeanB"]>=count.cut/mean(norm.factors))]
jj2=rownames(res)[which(res[,"baseMeanB"]==0 & res[,"baseMeanA"]>=count.cut/mean(norm.factors))]
ng=unique(c(ng,jj1,jj2))
}

if(length(ng)>0)
{
	#print(paste("there are", length(ng), "differentially expressed genes"))
	resSig=res[ng,]
	 
	jj=unique(rownames(resSig)[which(resSig[,"baseMeanA"]>= count.cut/mean(norm.factors) | resSig[,"baseMeanB"]>= count.cut/mean(norm.factors))])
	if(length(jj)>=1)
		ng=ng[ng%in%jj]
	
	if(length(jj)==0)
	   ng=NULL    

}	

#print("moving on...")
if(length(ng)==0){
		cat("\n===================================================\n")
		cat("Nothing passes the significant cutoff of",q.cut,"\n")
		cat("and has sufficient mean number of reads. At least", count.cut, "in one condition \n")
		m=max(res[which(res$padj==min(res$padj)),"baseMeanA"],res[which(res$padj==min(res$padj)),"baseMeanB"])
		cat("Best corrected p.value =",min(res$padj),"\n")
        cat("with max number of counts", m, "\n")
		cat("\n\n")
	} else {
		#print("getting ans ready")
                if ("GeneSymbol" %in% colnames(resSig)){
                    ans=resSig[ng,c("id","GeneSymbol","padj", "log2FoldChange","baseMeanA","baseMeanB")]
                    colnames(ans)=c("GeneID","GeneSymbol", "P.adj", paste("log2[",condB,"/",condA,"]",sep=""), paste("Mean_at_cond_",condA, sep=""), paste("Mean_at_cond_",condB, sep=""))
                } else {
                    ans=resSig[ng,c("id","padj", "log2FoldChange","baseMeanA","baseMeanB")]
                    colnames(ans)=c("GeneID", "P.adj", paste("log2[",condB,"/",condA,"]",sep=""), paste("Mean_at_cond_",condA, sep=""), paste("Mean_at_cond_",condB, sep=""))
                }
	        rownames(ans)=ans[,1]

		
#	pdf(paste("diffDESeq",DATE(),condA,condB,"FDR",q.cut,"FC",fc.cut,".pdf",sep="_"))
#			plot(res$baseMean, res$log2FoldChange)
#			points(resSig$baseMean, resSig$log2FoldChange, col="red")
#	dev.off()
}

### getting FC and p-values for all genes ###
#print("getting FC and p-values for all genes")
if ("GeneSymbol" %in% colnames(res)){
    res=res[,c("id","GeneSymbol","pval","padj", "log2FoldChange","baseMeanA","baseMeanB")]
    colnames(res)=c("GeneID","GeneSymbol", "pval","P.adj", paste("log2[",condB,"/",condA,"]",sep=""), paste("Mean_at_cond_",condA, sep=""), paste("Mean_at_cond_",condB, sep=""))
} else {
    res=res[,c("id","pval","padj", "log2FoldChange","baseMeanA","baseMeanB")]
    colnames(res)=c("GeneID", "pval","P.adj", paste("log2[",condB,"/",condA,"]",sep=""), paste("Mean_at_cond_",condA, sep=""), paste("Mean_at_cond_",condB, sep=""))
}
rownames(res)=res[,1]

Res=list()
Res$ans=ans
Res$counts.scaled=counts.scaled
Res$DEgenes=ng
Res$all.res=res
Res$q.cut=q.cut
Res$fc.cut=fc.cut
Res$count.cut=count.cut
Res$norm.factors=norm.factors

   
return(Res)
}


run.DESeq<-function(counts.dat, conds, condA, condB, q.cut=0.05, fc.cut=2, count.cut=10,zeroaddQ=F, libsizeQ=F,percentile="100%",fitType="parametric",orderPvalQ=T,method="per-condition",sharingMode="maximum",gns=c())
{

cds=make.cds(counts.dat=counts.dat,conds=conds,count.cut=count.cut,libsizeQ=libsizeQ,percentile=percentile,fitType=fitType,method=method,sharingMode=sharingMode)
#vsd = getVarianceStabilizedData(cds)
#dists=dist(t(vsd))

Res=run.nbinomTest(cds=cds, conds=conds, condA=condA, condB=condB, q.cut=q.cut, fc.cut=fc.cut, count.cut=count.cut,zeroaddQ=zeroaddQ,libsizeQ=libsizeQ,percentile=percentile,fitType=fitType,gns=gns)

#print("run.nbinomTest done...")
### done with run.nbinomTest() ###
###################################
### Deal with results for ALL genes ! ###
all.res=Res$all.res
### make sure there are no Inf FC ###
jj=which(abs(all.res[,paste("log2[",condB,"/",condA,"]",sep="")])==Inf)
if(length(jj)>=1)
	all.res[jj,paste("log2[",condB,"/",condA,"]",sep="")]=log2(as.numeric(all.res[jj,paste("Mean_at_cond_",condB, sep="")])+1)-log2(as.numeric(all.res[jj,paste("Mean_at_cond_",condA, sep="")])+1)

### Remove those that have less than count.cut normalized mean counts in either group... ###
jj=unique(which(all.res[,paste("Mean_at_cond_",condA, sep="")]>=count.cut | all.res[,paste("Mean_at_cond_",condB, sep="")]>=count.cut))
all.res=all.res[jj,]

### Order results for ALL genes ###
if(orderPvalQ)
	all.res=all.res[order(all.res[,"P.adj"],decreasing=F),]
else
	all.res=all.res[order(abs(all.res[,grep("log",colnames(all.res))]),decreasing=T),]

write.dat(all.res,file.name=paste("ALLResDESeq_",condB,"_vs_",condA,".xls",sep=""))
Res$all.res=all.res

#################################
### deal with results for DE genes ... ###
#print("dealing with FC for DE genes... ")
rr=Res$ans
if(!is.null(rr))
{
	jj=which(abs(rr[,paste("log2[",condB,"/",condA,"]",sep="")])==Inf)
if(length(jj)>=1)
	rr[jj,paste("log2[",condB,"/",condA,"]",sep="")]=log2(as.numeric(rr[jj,paste("Mean_at_cond_",condB, sep="")])+1)-log2(as.numeric(rr[jj,paste("Mean_at_cond_",condA, sep="")])+1)

if(orderPvalQ)
	rr=rr[order(rr[,"P.adj"],decreasing=F),]
else
	rr=rr[order(abs(rr[,grep("log",colnames(rr))]),decreasing=T),]
	
Res$ans=rr
write.dat(Res$ans,file=paste("ResDESeq_", condB, "_vs_", condA, ".xls",sep=""))
}
return(Res)
}
