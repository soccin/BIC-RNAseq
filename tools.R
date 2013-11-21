## Nick's functions ###
`cc` <-
function(...) {
    paste(...,sep='_')
}

`cq` <-
function(...) {
    paste(...,sep='')
}

`DATE` <-
function() {
gsub("-","",Sys.Date())
}

`ldensity` <-
function(x){
dx=density(log(x))
dx$x=exp(dx$x)
dx$y=dx$y/sum(dx$y)
return(dx)
}

`len` <-
function(x) {
length(x)
}

plus.paste<-function(str.v,plus.str="+")
{
	if(length(str.v)>1)	
	{
		tmp=str.v[1]
		for(nn in 2:length(str.v))
        	 tmp=paste(tmp,str.v[nn],sep=plus.str)
	
	return(tmp)
	}
	return(str.v)
}

split.plus<-function(x.str)
{
		if(length(grep("\\+",x.str))>=1)
	  		return(unlist(strsplit(x.str,"\\+")))  

return(x.str)
}

remove.dupl.entries<-function(InputDat,str="")
{
	tmp=apply(InputDat,1,plus.paste,plus.str=str)
	return(InputDat[-which(duplicated(tmp)),])	
}




drawlogaxis <- function(side,range)
{
par(tck=0.02)
#	d <- log(range,10)
d <- range
mlog <- floor(min(d))
Mlog <- ceiling(max(d))
SeqLog <- c(mlog:Mlog)
Nlog <- (Mlog-mlog)+1
axis(side,at=SeqLog,labels=10^SeqLog)
ats <- log(seq(from=2,to=9,by=1),10)
mod <- NULL
for(i in SeqLog)
{
mod <- c(mod,rep(i,length(ats)))
}
ats <- rep(ats,Nlog)
ats <- ats+mod
par(tck=0.02/3)
axis(side,at=ats,labels=NA)
}

logplot <- function(x,y,log='xy',...,forceylim=c(0,0),forcexlim=c(0,0))
{
par(tck=0.02)
xlg <- FALSE
ylg <- FALSE
if('x'%in%strsplit(log,'')[[1]]){x <- log(x,10);xlg=TRUE}
if('y'%in%strsplit(log,'')[[1]]){y <- log(y,10);ylg=TRUE}
yl <- ifelse(forceylim==c(0,0),range(y),forceylim)
xl <- ifelse(forcexlim==c(0,0),range(x),forcexlim)
plot(x,y,...,axes=FALSE,ylim=yl,xlim=xl)
if(xlg){drawlogaxis(1,xl)}else{axis(1,at=pretty(xl),labels=pretty(xl))}
if(ylg){drawlogaxis(2,yl)}else{axis(2,at=pretty(yl),labels=pretty(yl))}
box()
}

addlog <- function(x,y,log='xy',...)
{
xlg <- FALSE
ylg <- FALSE
if('x'%in%strsplit(log,'')[[1]]){x <- log(x,10);xlg=TRUE}
if('y'%in%strsplit(log,'')[[1]]){y <- log(y,10);ylg=TRUE}
points(x,y,...)

}

### my functions ###
rr<-function(dat,header = FALSE, sep = "", skip = 0)
       return(as.matrix(read.delim(dat,header=header,sep=sep,skip=skip)))

read.counts<-function(dat,header = FALSE, sep = "", skip = 0)
{
	tmp=rr(dat,header = header, sep = sep, skip = skip)
	tmp=make.rownames(tmp)
	return(matrix2numeric(tmp))
}


pdf.hclust<-function(dat,file.name="tmp.pdf",title="",width=26,height=16,lwd=3,cex.main=3,cex.lab=3,cex=3,xlab="",ylab="")
{
	pdf(file.name, width=width,height=height)
	plot(hclust(dist(t(dat))), lwd=lwd, main=paste(title), cex.main=cex.main, xlab=xlab,ylab=ylab, cex.lab=cex.lab, cex=cex)
dev.off()
}

pdf.hclust.cor<-function(dat,file.name="tmp.pdf",title="")
{
	pdf(file.name, width=26,height=16)
	plot(hclust(dist(t(dist(1-cor(dat))))), lwd=3, main=paste(title), cex.main=3, xlab="",ylab="", cex.lab=3, cex=3)
dev.off()
}

pdf.sam<-function(dat,file.name="tmp.pdf",title="")
{
	pdf(file.name)
	
	sam1=sammon(dist(t(dat)),trace=FALSE)
	plot(sam1$points,lwd=2,cex.main=1.5, xlab="",ylab="",cex.lab=1.5,cex=1.5, main=paste(title))
	text(sam1$points, labels = colnames(dat), cex=1.5)

    dev.off()
}
    	
    	make.matrix<-function(mat,row.names="",col.names="")
{
	if(length(mat)==0)
	   return(NULL)
	   
	if(is.matrix(mat))
	  {
	  	mat1=mat
	  	#if(row.names!="")
	  	#   rownames(mat)=row.names
	  	#if(!col.names!="")
	  	#  colnames(mat)=col.names   
	  }
	  if(is.vector(mat))#!is.matrix(mat))
	   	mat1=matrix(mat,nrow=1)
	   		
	   if(length(col.names)==ncol(mat1))
		   			colnames(mat1)=col.names
		   			
	   
	   if(length(row.names)==nrow(mat1))
			   		rownames(mat1)=row.names
		   	
	   return(mat1)
		

}

rbind2mat<-function(mat1,mat2)
{
		if(is.null(mat1) & !is.null(mat2))
		  return(mat2)
		  
	    if(!is.null(mat1) & is.null(mat2))
		  return(mat1)   		
    	
    	if(!is.null(mat1) & !is.null(mat2))
  	     { 
  	     	if(ncol(mat1)!=ncol(mat2))
  	     		return("Error")
  	     	tmp=as.matrix(rbind(mat1,mat2))
  	     	rownames(tmp)=c(rownames(mat1),rownames(mat2))
  	     	return(tmp) 		
  	     } 		
}





## take a matrix, where columns correspond to samples and rows to value of each species
## generate all products [and ratios ?!] from the rows of this matrix

generate.prod.mat<-function(dat)
{
	inds=combn(1:nrow(dat),2)
  labls=NULL
  vals=colnames(dat)
for(nn in 1:ncol(inds))
{
	pair=inds[,nn]
	vals=as.matrix(rbind(vals,dat[pair[1],]*dat[pair[2],]))
	labls=c(labls,paste(rownames(dat)[pair[1]],"_",rownames(dat)[pair[2]],sep="")) 
}
vals=matrix2numeric(make.colnames(vals))
rownames(vals)=labls

return(vals)	
}


len<-function(x)
  length(x)
  
  
remove.NAs<-function(v)
{
	return(v[!is(v)])
}


countNA<-function(v)
{
	length(which(!is.na(v)))
}



### Nick's function for CGH data ###
S<-function(x){
    px=table(x)
    px=px/sum(px)
    return(-sum(px*log2(px)))
}



len.which<-function(v, cut.off,alt="l")
{
	if(alt=="g")
   	{
   		return(length(which(v>=cut.off)))
   	}	
   	if(alt=="l")
   	{
   		return(length(which(v<=cut.off)))
   		}
}



### remove hsa ### 
#remove.sub<-function(x, substr="hsa"){
#	 tmp=NULL
#	 for(n in 1:length(x)){
#	 	   if(length(grep(substr,x[n]))>0){
#	 	   	 tmp = c(tmp, unlist(strsplit(x[n],substr))[[2]]) 
#	 	   	}
#	 	   	if(length(grep(substr,x[n]))==0){
#           tmp = c(tmp, x[n])
#	 	}
#	}
#return(tmp)
#}

remove.sub<-function(x, substr="hsa", part2take=1){
tmp=NULL
for(n in 1:length(x)){
#print(paste("n=",n))
if(length(grep(substr,x[n]))>0)
{ 
  #print(paste("new name=",unlist(strsplit(x[n],substr))))
  tmp = c(tmp, unlist(strsplit(x[n],substr))[[part2take]]) 
}
if(length(grep(substr,x[n]))==0)
{
  tmp = c(tmp, x[n])
}
}
return(tmp)
}



get.this.samples<-function(the.sample,runs)
{
	return(runs[grep(the.sample,runs)])
}
	
	
get.samples <- function(samples, runs){
	       return(apply(as.matrix(samples), 1, get.this.samples, runs=runs))
	}


get.this.sample.nu<-function(the.sample,runs)
{
	n1 = grep(the.sample,runs) 
	if(length(n1)==1)
	{return(n1)}
	if(length(n1)!=1)
	{return(NA)}
}
	
	
get.samples.nu <- function(samples, runs){
	       return(apply(as.matrix(samples), 1, get.this.sample.nu, runs=runs))
	}




freq.compute.gn<-function(v, val=1)
{
	 if(val==1)
	 {
	  return(length(which(v==1 | v== 2))/length(v))
  }
  if(val==-1)
	 {
	  return(length(which(v==-1 | v== -2))/length(v))
  }
}

freq.compute<-function(v, val=1)
{
	length(which(v==val))/length(v)
}



mytable<-function(v)
{
	c(length(which(v==-1)), length(which(v==0)), length(which(v==1)))
}
	
	
	how.many.gr<-function(v, cut.off=0.6)
{
	length(which(v>=cut.off))
}

how.many.ls<-function(v, cut.off=-0.6)
{
	length(which(v<= cut.off))
}


rescale<-function(v)
{
	m1=min(v,na.rm=TRUE) 
	m2=max(v,na.rm=TRUE) 
	
	return((v-m1)/(m2-m1))
}

rescale2<-function(v,m1=0,m2=1)
{
 return((v-m1)/(m2-m1))
}




`DATE` <-
function() {
    gsub("-","",Sys.Date())
}

add.str<-function(x, str)
{
	tmp=NULL
for(i in 1:length(x))
{
	  tmp=c(tmp, paste(x[i], str, sep=""))
}
return(tmp)
	
	}
	
	
	
### from a vector X return numbers from the grep(X, Y[n]) ###

grepvector<-function(X, Y,no.doubles=TRUE){
	  res=NULL
	  for(n in 1:length(Y)){
	  	    rr = grep(Y[n], X)
         if(no.doubles)
	  	    {
	  	    	if(length(rr)>=2){
	  	    	   print(paste("n=", n, " problem with", Y[n], "double entries..."))
	  	    	   rr=rr[1]
	  	    	}
	  	    	}
	  	    	res=c(res, rr)
	  	}
	  	return(res)
	}


make.subset.colnames<-function(dat, vars)
{
     col.num = grepvector(colnames(dat),vars)
     
     	   return(dat[,col.num])

	}   


matrix2round<-function(dat){
	res=matrix(round(dat),nrow=nrow(dat))
	rownames(res)=rownames(dat)
	colnames(res)=colnames(dat)
	return(res)
}



matrix2integer<-function(dat){
	res=matrix(round(dat),nrow=nrow(dat))
	rownames(res)=rownames(dat)
	colnames(res)=colnames(dat)
	return(res)
	
}

### look for missing values ###
vals2NA<-function(dat){
	 cols.with.NAs=NULL
	 rows.with.NAs=NULL
	 for(i in 1:ncol(dat)){
	 	   ll=which(is.na(dat[,i]))
	 	   if(length(ll)>=1){
	 	   	  cols.with.NAs=c(cols.with.NAs, i)
	 	   	  rows.with.NAs=c(rows.with.NAs,ll)
	 	   	}
	 	}
	return(list(cols.with.NAs=cols.with.NAs, rows.with.NAs=unique(rows.with.NAs)))
	}




write.dat<-function(dat, file.name, add.rownames=F,col.names=T,quote=F,sep="\t")
{
	 	write.table(dat, file=file.name, col.names=col.names, row.names=add.rownames, quote=quote, sep=sep)
}

make.conds<-function(conds,condA.str,condB.str,condA,condB)
{
	h=conds
	h[grep(condA.str,h)]=condA
	h[grep(condB.str,h)]=condB
	
	return(h)
}


make.rownames<-function(dat, nn=1)
{
	rownames(dat)=dat[,nn]
	return(dat[,-nn])
}

add.rownames2dat<-function(dat,h="")
{
	tmp=as.matrix(cbind(rownames(dat),dat))
	colnames(tmp)[1]=h
    return(tmp)
}

make.colnames<-function(dat,nn=1)
{
	colnames(dat)=dat[nn,]
	return(dat[-nn,])
}



abs.mean<-function(x)
{
	return(mean(abs(x)))
	}

abs.sum<-function(x)
{
   sum(abs(x), na.rm=TRUE)
}



how.many.present<-function(v){
	length(which(v!="A"))
	}

transform4heatmap<-function(x)
{
	return(x - apply(x,1,mean))
}

average.by.name<-function(dat){
	print(paste("dim",dim(dat)))
	ugn <- unique(dat[,1])
	dat.mean <- rep(0, ncol(dat)-1)
	for(nn in 1:length(ugn))
	{
		kk <- which(dat[,1]==ugn[nn])
		if(length(kk)>=2){
		    #val <- mean(as.numeric(dat[kk,2]))
		    val<-apply(matrix(as.numeric(dat[kk,-1]),nrow=length(kk)), 2, mean)
		}
		if(length(kk)==1){
			val=as.numeric(dat[kk,-1])
	   	}
	dat.mean<-as.matrix(rbind(dat.mean, val))  	}
	dat.mean<-dat.mean[-1,]
	return(list(dat.mean=dat.mean, ugn=ugn))
}

### this is averaging for a two-column matrix:
### 1st column has names and the other values to be averaged=Fold-changes 
average.by.name2<-function(dat,func2use="mean")
{
	print(paste("dim",dim(dat)))
	ugn <- unique(dat[,1])
	dat.mean <- rep(0, ncol(dat)-1)
	for(nn in 1:length(ugn))
	{
		kk <- which(dat[,1]==ugn[nn])
		if(length(kk)>=2){
			 if(func2use=="mean")
           val <- mean(as.numeric(dat[kk,2]))
		   if(func2use=="min")
		      val <- min(as.numeric(dat[kk,2])) #for p-values! 
		     if(func2use=="max")
		      val <- max(as.numeric(dat[kk,2])) #for p-values!   

		}
		if(length(kk)==1){
			val=as.numeric(dat[kk,-1])
	   	}
	dat.mean<-as.matrix(rbind(dat.mean, val))  	}
	dat.mean<-dat.mean[-1,]
	return(list(dat.mean=dat.mean, ugn=ugn))
}




 log2_2fc<-function(fc){
 # this function trasnfers log2 FC
 # to the form that Nick uses for investigators
 
    if(fc<0){return(-1/2^fc)}
    if(fc>0){return(1/2^fc)} 
   }


remove.last.space<-function(str){
res=str

if(length(grep(" ", str))>0)	
		res=unlist(strsplit(str," "))[[1]]
		
return(res)		
}


get.array.nu<-function(s){
	tt<-unlist(strsplit(unlist(strsplit(s,".CEL")),"_"))
	return(tt[length(tt)]) 
}




remove.first.space<-function(str){
	s<-strsplit(str," ")[[1]]
	if(length(s)==2 & s[1]=="")
	{return(s[[2]])}
	if(length(s)){return(s)}
}

matrix2numeric<-function(dat){
	res=matrix(as.numeric(dat),nrow=nrow(dat))
	rownames(res)=rownames(dat)
	colnames(res)=colnames(dat)
	return(res)
}
matrix2floor<-function(dat)
{
	res=matrix(floor(dat),nrow=nrow(dat))
	rownames(res)=rownames(dat)
	colnames(res)=colnames(dat)
	return(res)
}

maploc<-function(probe_id)
{
    x=toTable(hgu133a2CHRLOC[probe_id])    
    if(nrow(x)>1){
    	   x = x[1,]
    	}
    	x$start_location=abs(x$start_location)   
    	return(as.matrix(x))
 }


### limma analysis ###
limma.analysis<-function(dat, groups, comp.groups,q.cut=0.05, lfc=log2(2), lib2use="lumiHumanAll.db", write.tables=T, write.gnlists=T, all.genes=TRUE)
{
### all.genes=FALSE 
### is a flag to write down only differentially expressed genes
### if all.genes=TRUE
### all genes, and their fold-changes are written 

design <- model.matrix(~0+groups)
colnames(design) = levels(groups)
fit=lmFit(dat,design)
 contrast.matrix <- makeContrasts(contrasts=comp.groups, levels=design)
 fit2 <- contrasts.fit(fit, contrast.matrix) 
 fit2 <- eBayes(fit2)
 
 res=decideTests(fit2,p.value=q.cut,lfc=lfc)
 coefs=colnames(res)
 
 ### go to geneSymbols in res ###
 TABs = vector(length=length(coefs),mode="list")
 Gene.Lists.UP = vector(length=length(coefs),mode="list")
 Gene.Lists.DOWN = vector(length=length(coefs),mode="list")
 
 for(nn in 1:length(coefs)){
 print(paste("analyzing", coefs[nn], sep=" "))
 if(all.genes){ #to output all genes 
 	   tab<-topTable(fit2, coef = coefs[nn], adjust = "BH", number=nrow(fit2), sort.by='logFC')
 	   }
 	   if(!all.genes) #to output only DE genes
 	   tab<-topTable(fit2, coef = coefs[nn], adjust = "BH", number=nrow(fit2), sort.by='logFC',p.value=q.cut, lfc=lfc)

 h=colnames(tab)
 rownames(tab)=tab[,1]
 gns = as.vector(getSYMBOL(rownames(tab), lib2use))
 jj=which(!is.na(gns))
 tab2 = as.matrix(cbind(gns,tab[,c(2,5,6)]))
 rownames(tab2)=tab2[,1]
 tab2=tab2[jj,]
 TABs[[nn]]<- tab2
 colnames(TABs[[nn]])=c("GeneSymbol",h[c(2,5,6)])

 ### add here the cut-offs ##3 
 tmp=coefs[nn]
 tmp=unlist(strsplit(tmp, " - "))
 file.name=paste(tmp[1], "_vs_", tmp[2], ".txt",sep="")
 if(write.tables)
 {  
 	  write.table(as.matrix(TABs[[nn]]), file=file.name,col.names=T,row.names=F,quote=F,sep="\t")
     
  }
  
 ###GeneLists###
 if(write.gnlists)
 {
 	gg1=rownames(res)[which(res[,nn]==1)]
  if(length(gg1)>=1)
  {
  gns.up= as.vector(getSYMBOL(gg1, lib2use))
  gns.up=unique(gns.up[!is.na(gns.up)])
  Gene.Lists.UP[[nn]]=gns.up
  tmp=coefs[nn]
 tmp=unlist(strsplit(tmp, " - "))
 file.name=paste("UP_gns", tmp[1], "_vs_", tmp[2], ".txt",sep="")
 write.table(as.vector(Gene.Lists.UP[[nn]]), file=file.name,col.names=F,row.names=F,quote=F,sep="\t")
 
  }else
  {
  	  print(paste("no UP DE genes"))
  	}
  
 gg2=rownames(res)[which(res[,nn]==-1) ]
 
 if(length(gg2)>=1){
 gns.down= as.vector(getSYMBOL(gg2, lib2use))
 gns.down=unique(gns.down[!is.na(gns.down)])
 Gene.Lists.DOWN[[nn]]=gns.down
 tmp=coefs[nn]
 tmp=unlist(strsplit(tmp, " - "))
 file.name=paste("DOWN_gns", tmp[1], "_vs_", tmp[2], ".txt",sep="")
 write.table(as.vector(Gene.Lists.DOWN[[nn]]), file=file.name,col.names=F,row.names=F,quote=F,sep="\t")

 }else
  {
  	  print(paste("no DOWN DE genes"))
  	}
 
 
 
  
  }
}
return(res)
}


limma.analysis4Sam<-function(dat, groups, comp.groups,q.cut=0.05, lfc=log2(2), lib2use="lumiHumanAll.db", write.tables=T, write.gnlists=T, all.genes=TRUE, affyFC=TRUE, more.infoQ=TRUE, add.expr.conds=TRUE)
{
### all.genes=FALSE 
### is a flag to write down only differentially expressed genes
### if all.genes=TRUE
### all genes, and their fold-changes are written 

design <- model.matrix(~0+groups)
colnames(design) = levels(groups)
fit=lmFit(dat,design)
 contrast.matrix <- makeContrasts(contrasts=comp.groups, levels=design)
 fit2 <- contrasts.fit(fit, contrast.matrix) 
 fit2 <- eBayes(fit2)
 
 res=decideTests(fit2,p.value=q.cut,lfc=lfc)
 
 coefs=colnames(res)
 
 ### go to geneSymbols in res ###

 TABs = vector(length=length(coefs),mode="list")
 Gene.Lists.UP = vector(length=length(coefs),mode="list")
 Gene.Lists.DOWN = vector(length=length(coefs),mode="list")
 
 for(nn in 1:length(coefs)){
 print(paste("analyzing", coefs[nn], sep=" "))
 condA = unlist(strsplit(coefs[nn]," - "))[[1]]
 condB = unlist(strsplit(coefs[nn]," - "))[[2]]
 if(all.genes){ #to output all genes 
 	   tab<-topTable(fit2, coef = coefs[nn], adjust = "BH", number=nrow(fit2), sort.by='logFC')
 	   }
 	   if(!all.genes) #to output only DE genes
 	   {
 	   	tab<-topTable(fit2, coef = coefs[nn], adjust = "BH", number=nrow(fit2), sort.by='logFC',p.value=q.cut, lfc=lfc)
     }
 h=colnames(tab)
 rownames(tab)=tab[,1]
 gns = as.vector(getSYMBOL(rownames(tab), lib2use))
 jj=which(!is.na(gns))
 ### remove probes with no geneSymbols ###
 tab=tab[jj,]
 gns=gns[jj]
 #############
 tab2 = as.matrix(cbind(gns,tab[,c(2,6)]))
 #rownames(tab2)=tab2[,1]
 if(affyFC)# ### if log2FC format use Nick's function ###
 {
 	   tab2[,2]=rn(log2FC(as.numeric(tab2[,2])))
 }
  ### adding gene description and chromosome location ###
      	maploc=mget(rownames(tab2),hgu133a2CHRLOC,ifnotfound=NA)
      	pr=NULL
      	chrID=NULL
      	start.loc=NULL
  	     for(ii in 1:length(maploc))
  	       {
  	       	
  	     	   pr=c(pr,names(maploc)[ii])
  	     	   strt.tmp=maploc[[ii]][1]
  	     	   chr.tmp=NA
  	     	   if(!is.na(strt.tmp)){
  	     	   	  chr.tmp=names(maploc[[ii]])[1]
  	     	     }
  	     	     
            chrID = c(chrID, chr.tmp)
            start.loc=c(start.loc, abs(strt.tmp))
  	     	} 
  	     gene.info=toTable(hgu133a2GENENAME[rownames(tab2)])
  	     rownames(gene.info)=gene.info[,1]
  	     gene.info=gene.info[rownames(tab2),2]
  	     gene.info=cbind(gene.info, chrID, start.loc)
  	
 	#add mean expression in each condition
   meanA = apply(dat[,which(colnames(dat)==condA)],1,mean)[rownames(tab2)]
 	 meanB = apply(dat[,which(colnames(dat)==condB)],1,mean)[rownames(tab2)]
   tab2 = as.matrix(cbind(rownames(tab2),tab2, meanA, meanB,gene.info))
   colnames(tab2)=c("ProbeID", "GeneSymbol", "FC", "FDR", paste("mean_expr_",condA,sep=""),paste("mean_expr_",condB,sep=""),"GeneDescr", "ChrID", "StartLoc")
    
 TABs[[nn]]<- tab2
 
 ### add here the cut-offs ##3 
 tmp=coefs[nn]
 tmp=unlist(strsplit(tmp, " - "))
 file.name=paste(tmp[1], "_vs_", tmp[2], ".txt",sep="")
 if(write.tables)
 {  
 	  write.table(as.matrix(TABs[[nn]]), file=file.name,col.names=T,row.names=F,quote=F,sep="\t")
     
  }
  
 ###GeneLists###
 if(write.gnlists)
 {
 	gg1=rownames(res)[which(res[,nn]==1)]
  if(length(gg1)>=1)
  {
  gns.up= as.vector(getSYMBOL(gg1, lib2use))
  gns.up=unique(gns.up[!is.na(gns.up)])
  Gene.Lists.UP[[nn]]=gns.up
  tmp=coefs[nn]
 tmp=unlist(strsplit(tmp, " - "))
 file.name=paste("UP_gns", tmp[1], "_vs_", tmp[2], ".txt",sep="")
 write.table(as.vector(Gene.Lists.UP[[nn]]), file=file.name,col.names=F,row.names=F,quote=F,sep="\t")
 
  }else
  {
  	  print(paste("no UP DE genes"))
  	}
  
 gg2=rownames(res)[which(res[,nn]==-1) ]
 
 if(length(gg2)>=1){
 gns.down= as.vector(getSYMBOL(gg2, lib2use))
 gns.down=unique(gns.down[!is.na(gns.down)])
 Gene.Lists.DOWN[[nn]]=gns.down
 tmp=coefs[nn]
 tmp=unlist(strsplit(tmp, " - "))
 file.name=paste("DOWN_gns", tmp[1], "_vs_", tmp[2], ".txt",sep="")
 write.table(as.vector(Gene.Lists.DOWN[[nn]]), file=file.name,col.names=F,row.names=F,quote=F,sep="\t")

 }else
  {
  	  print(paste("no DOWN DE genes"))
  	}
  }
}
return(res)
}

log2FC <- function(x) {
  return(ifelse(x<0,-2^(-x),2^x))
}

rn <- function(x, digits = 3){
  if (is.null(x)) 
    NULL
  else {
    if (is.matrix(x) && ncol(x) == 1) 
      x <- x[, 1]
    round(x, digits = digits)
  }
}



##############################################################
### write the names of genes in different parts of venn diagram to .xls file ###
### res is an otput from limma in gene symbols ###
 make.gnlists4venn<-function(res, gns, k1, k2, fl.name="genes_in_venn_diagram")
{
s1 = colnames(res)[k1]
s2 = colnames(res)[k2]
print(paste("comparing", s1, " and ", s2, sep=""))
### maxumim number of genes in one of the parts of the venn diagram ###
jj1 = which(res[,k1]==1 & res[,k2]==1)
jj2 = which(res[,k1]==-1 & res[,k2]==-1)
jj3=which(res[,k1]==1 & res[,k2]!=1)
jj4=which(res[,k1]==-1 & res[,k2]!=-1)
jj5 = which(res[,k1]!=1 & res[,k2]==1)
jj6= which(res[,k1]!=-1 & res[,k2]==-1)

ll = max(c(length(jj1), length(jj2), length(jj3), length(jj4), length(jj5), length(jj6)))

nn = rep(" ", ll) #max number of genes in this venn diagram... 
tmp = nn
tmp[1:length(jj1)]=gns[jj1]
out=tmp

print(gns[jj1])

tmp = nn
tmp[1:length(jj2)]=gns[jj2]
out=as.matrix(cbind(out,tmp))

tmp = nn
tmp[1:length(jj3)]=gns[jj3]
out=as.matrix(cbind(out,tmp))

tmp = nn
tmp[1:length(jj4)]=gns[jj4]
out=as.matrix(cbind(out,tmp))

tmp = nn
tmp[1:length(jj5)]=gns[jj5]
out=as.matrix(cbind(out,tmp))


tmp = nn
tmp[1:length(jj6)]=gns[jj6]
out=as.matrix(cbind(out,tmp))


colnames(out)=c(paste("UP", s1,"AND", s2,sep="_"), paste("DOWN", s1, "AND", s2, sep="_"), paste("UP", s1,"NOT", s2,sep="_"), paste("DOWN", s1,"NOT", s2,sep="_"), 
paste("NOT", s1,"UP", s2,sep="_"), paste("NOT", s1,"DOWN", s2,sep="_"))

write.table(out, file=paste(fl.name,".xls",sep=""), col.names=T,row.names=F,quote=F,sep="\t")

}







run.limma<-function(dat, groups, comp.groups,q.cut=q.cut, lfc=lfc)
{
	#comp.groups=paste(levels(groups)[[1]]," - ", levels(groups)[[2]],sep="")
	#print(comp.groups)
	#print(colnames(dat))
	design <- model.matrix(~0+groups)
	colnames(design) = levels(groups)
	fit=lmFit(dat,design)
	contrast.matrix <- makeContrasts(comp.groups,levels=design)
	fit2 <- contrasts.fit(fit, contrast.matrix) 
	fit2 <- eBayes(fit2)
	res=decideTests(fit2,p.value=q.cut,lfc=lfc)
	tab<-topTable(fit2, coef = comp.groups, adjust = "fdr", number=nrow(fit2))
	
	return(list(res=res,fit=fit2,tab=tab))
}


### auxiliary function to write results into file ###
write.limmaresults2file<-function(tab,groups,lib.file,file2write){
	nu=tab[,1]
  gene.symbs=as.vector(getSYMBOL(nu, lib.file))

  tabb = as.matrix(cbind(gene.symbs,tab$logFC,tab$adj.P.Val))
  colnames(tabb)=c("GeneSymb", paste("log[expr_", levels(groups)[[1]],"/expr_",levels(groups)[[2]],"]",sep=""), "adj.Pvalue")

  write.table(tabb, file=file2write,row.names=FALSE, col.names=TRUE,sep="\t")

	}
	
	
	 make.table <- function(tab, gene.symbs,col2choose=c(2,5)){
	 g1 <- tab[,1]
	 for(i in 1:length(g1)){
	 the.gene <- g1[i]
	 the.symb <- gene.symbs[nu.IDs==the.gene]
	 entry <- as.vector(c(the.gene, the.symb, as.vector(tab[i,col2choose])))
	 if(i==1){tab.dat <-  entry}
	 if(i>1){tab.dat<-as.matrix(rbind(tab.dat, entry))}
	 }
	 colnames(tab.dat) <- c("ID", "GeneSymb", colnames(tab)[col2choose])
	 rownames(tab.dat)=rep("", nrow(tab.dat))
	 
	 return(tab.dat)
	 }
	 
	 make.table.from.two <- function(tab1, tab2, gene.symbs, title1, title2){
	 g1 <- tab1[,1]
	 for(i in 1:length(g1)){
	 the.gene <- g1[i]
	 #see if this gene exists in table2
	 k <- which(tab2[,1]==the.gene)
	 if(length(k)==1){
	 the.symb <- gene.symbs[nu.IDs==the.gene]
	 entry <- as.vector(c(the.gene, the.symb, tab1[i,2], tab2[k,2], tab1[i,6], tab2[k,6]))
	 if(i==1){tab.dat <-  entry}
	 if(i>1){tab.dat<-as.matrix(rbind(tab.dat, entry))}
	 }
	 }
	 colnames(tab.dat) <- c("ID", "GeneSymb", paste("logFC",title1), paste("logFC",title2), paste("adjPval",title1), paste("adjPval",title2))
	 rownames(tab.dat)=rep("", nrow(tab.dat))
	 
	 return(tab.dat)
	 }


### auxiliary functions ###
### this is a little auxiliary function that takes things/names in the vector
### and puts them into a string 
vector2str<-function(vv)
{
	res=""
	for(nn in 1:length(vv))
		res=paste(res,vv[nn],sep=",")		
		
	return(res)
}

## this function is for Pathway Enrichment Analysis ###
getgns4pathway<-function(gnsetofinterest,pathway.name,geneSetCollection,universeGns,universeGnsIds)
{
	kk=which(names(geneSetCollection)==pathway.name)
	gns.indxs=which(universeGnsIds%in%as.numeric(geneIds(geneSetCollection[[kk]])))
	gns.names=universeGns[gns.indxs]
	return(gnsetofinterest[which(gnsetofinterest%in%gns.names)])
}

