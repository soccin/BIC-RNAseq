make.binary.mat<-function(comp,pval.cutoff)
{
	fls=grep("up|down",grep(comp,dir(),value=T),value=T)
	res1=rr(fls[1],header=T)
	res2=rr(fls[2],header=T)
	res1=res1[which(as.numeric(res1[,2])<pval.cutoff),]
	res2=res2[which(as.numeric(res2[,2])<pval.cutoff),]
	res1[,2]=-1
	res2[,2]=1
	res=rbind(res1,res2)
	colnames(res)=c("GeneSet","Dir")
	
	return(res)
}	


getGeneSets<-function(gsaRes,gsa.tab,fc2keep=log2(1.5),frac2keep=2,dir="Up",fcQ=T)
{
	tmp=cbind(gsa.tab,rep("",nrow(gsa.tab)),rep("",nrow(gsa.tab)))
    colnames(tmp)=c(colnames(gsa.tab),"GenesInGeneSet","GenesAndFC")
	geneSets2remove=NULL
	for(nn in 1:nrow(gsa.tab))
	{
		fc=geneSetSummary(gsaRes, gsa.tab[nn,1])$geneLevelStats
		o=order(abs(fc),decreasing=T)
		fc=fc[o]
		gn.names=rownames(as.matrix(fc))
		gn.names.str=add.str(gn.names,",")
		gn.names.fc.str=add.str(paste(gn.names,"=",fc),";")
		tmp[nn,ncol(gsa.tab)+1]=gn.names.str
		tmp[nn,ncol(gsa.tab)+2]=gn.names.fc.str

    	## Add a condition that at least half of the genes have to have 
    	## FC >= log2(1.5) ! 
    	if(dir=="Dn")
	    	jj=length(which(fc<= -fc2keep))
	    if(dir=="Up")
	    	jj=length(which(fc>= fc2keep))

    	if(jj <= length(fc)/frac2keep)
	       	geneSets2remove=c(geneSets2remove,nn)    	
	}
	if(length(geneSets2remove)>=1 & fcQ)
			tmp=tmp[-geneSets2remove,]

return(tmp)	
}

processgsaRes<-function(gsaRes,pval.cutoff=0.01,fc2keep=log2(1.5),frac2keep=2,fcQ=T)
{
	gsa.tab.dn=NULL
	gsa.tab.up=NULL
	tmp=as.matrix(GSAsummaryTable(gsaRes))
	if(!is.null(tmp))
	{
		print(tmp[1:3,])
		
		gsa.tab.dn=NULL
		gsa.tab.up=NULL
	
		indxDn=which(tmp[,"p adj (dist.dir.dn)"]<=pval.cutoff)
		indxUp=which(tmp[,"p adj (dist.dir.up)"]<=pval.cutoff)
		colsDn=c("Name","Genes (tot)", "Stat (dist.dir)","p adj (dist.dir.dn)","Genes (down)")
		colsUp=c("Name","Genes (tot)", "Stat (dist.dir)","p adj (dist.dir.up)","Genes (up)")
		if(length(indxDn)>=1)
		{
                        if(length(indxDn) == 1){
                            tmpDn=t(as.matrix(tmp[indxDn,colsDn]))
                        } else {
			    tmpDn=as.matrix(tmp[indxDn,colsDn])
                        }
			if(length(colsDn)>=2)
		    	tmpDn=as.matrix(tmpDn[order(abs(as.numeric(tmpDn[,"Stat (dist.dir)"])),decreasing=T),])
	                        if(length(indxDn) == 1){
                                    tmpDn = t(tmpDn)
                                }
	   			gsa.tab.dn=getGeneSets(gsaRes,gsa.tab=tmpDn,fc2keep=fc2keep,frac2keep=frac2keep,dir="Dn",fcQ=fcQ)
		}
		if(length(indxUp)>=1)
		{
                        if(length(indxUp) == 1){
                            tmpUp=t(as.matrix(tmp[indxUp,colsUp]))
                        } else {
			    tmpUp=as.matrix(tmp[indxUp,colsUp])
                        }
	    	if(length(colsUp)>=2)
	    		tmpUp=as.matrix(tmpUp[order(abs(as.numeric(tmpUp[,"Stat (dist.dir)"])),decreasing=T),])
   			if(length(indxUp) == 1){
                            tmpUp = t(tmpUp)
                        }
   			gsa.tab.up=getGeneSets(gsaRes,gsa.tab=tmpUp,fc2keep=fc2keep,frac2keep=frac2keep,dir="Up",fcQ=fcQ)
		}
	}	      
  tab.res = list()
  tab.res$gsa.tab.dn = gsa.tab.dn
  tab.res$gsa.tab.up = gsa.tab.up 
  #return(list(gsa.tab.dn=gsa.tab.dn,gsa.tab.up=gsa.tab.up))
  return(tab.res)
}


run.pianoGSA<-function(fls,dir4res,gs.names=c("c1","c2","c3","c5-1", "c6","c7"),dir.gs="/home/raya/Tools/Data/",min.gns.nu=10,max.gns.nu=500,nPerm=100,pval.cutoff=0.01,ncol.gns=NULL,fc2keep=log2(1.5),frac2keep=2,fcQ=T,microarrayQ=F)
{
	gs.names2load=pastec(gs.names,".all.v4.0.symbols.gmt")
	print(gs.names2load)
	
	#c1.all.v4.0.symbols.gmt
	#c2.all.v4.0.symbols.gmt
	#c3.all.v4.0.symbols.gmt
	#c5-1.all.v4.0.symbols.gmt": pathways... Gene sets are named by GO term and contain genes annotated by that term. GSEA users: Gene set enrichment analysis identifies gene sets consisting of co-regulated genes; GO gene sets are based on ontologies and do not necessarily comprise co-regulated genes
	#c6-1.all.v4.0.symbols.gmt") Gene sets represent signatures of cellular pathways which are often dis-regulated in cancer. The majority of signatures were generated directly from microarray data from NCBI GEO or from internal unpublished profiling experiments which involved perturbation of known cancer genes. In addition, a small number of oncogenic signatures were curated from scientific publications. msigdb.v4.0.symbols.gmt")

	tmp<-capture.output(suppressMessages(library("piano")))

   ## dir.gs="/home/raya/Tools/Data/" : directory where gene sets from the Broad are 
	## dir4res = directory where results will be stored 
	### ncol.gns: this is the column with gene symbols ###

	directions=c("up","down")
	res.gs=c("GeneSets","Pval","Dir")
	for(the.fl in fls)
	{
	   res=as.matrix(read.csv(the.fl ,sep="\t"))
   
	if(!microarrayQ)
	{
		 name2compare=remove.sub(remove.sub(remove.sub(the.fl,".xls",part2take=1),"ALLResDESeq_",part2take=2),"_Genes",part2take=1)
		 print(paste("working on ", name2compare))
 
 		if(ncol.gns==2)
			res=res[,-1]

		res=make.rownames(res)
		res=matrix2numeric(res)
	
		fc=as.matrix(res[,3])
		rownames(fc)=rownames(res)
		pvals=as.matrix(res[,1])
		rownames(pvals)=rownames(res)
	}
	if(microarrayQ)
	{
		name2compare=remove.sub(the.fl,".txt",part2take=1)
		 print(paste("working on ", name2compare))
 
		res=make.rownames(res)
		res=matrix2numeric(res)
		fc=as.matrix(res[,1])
		rownames(fc)=rownames(res)
		pvals=as.matrix(res[,2]) #adjusted p-value
		rownames(pvals)=rownames(res)
	}

	res.gs=c("GeneSets","Pval","Direction")
	
	for(the.gs in gs.names2load)
	{	
		gsc=loadGSC(pastec(dir.gs,the.gs))
		gs.name=remove.sub(remove.sub(the.gs,".all.v4.0.symbols.gmt",part2take=1),"-1",part2take=1)

		gsaRes <- runGSA(fc,geneSetStat="mean",gsc=gsc, nPerm=nPerm,gsSizeLim=c(min.gns.nu,max.gns.nu))
		tab.res=processgsaRes(gsaRes,pval.cutoff=pval.cutoff,fc2keep=fc2keep,frac2keep=frac2keep,fcQ=fcQ)
		if(!is.null(tab.res$gsa.tab.dn))
		{	
			file.name=paste("GeneSet_Dn_",remove.sub(remove.sub(the.fl, "ALLResDESeq_",part2take=2),"_codingGenes",part2take=1),"_",gs.name,".xls",sep="")
			write.dat(tab.res$gsa.tab.dn,file.name=file.name)
		}

		if(!is.null(tab.res$gsa.tab.up))
		{
			file.name=paste("GeneSet_Up_",remove.sub(remove.sub(the.fl, "ALLResDESeq_",part2take=2),"_codingGenes",part2take=1),"_",gs.name,".xls",sep="")
			write.dat(tab.res$gsa.tab.up,file.name=file.name)
		}


		#gsaRes2 <- runGSA(fc,geneSetStat="median",gsc=gsc,nPerm=nPerm,gsSizeLim=c(min.gns.nu,max.gns.nu))
		#gsaRes3 <- runGSA(fc,geneSetStat="sum",gsc=gsc,nPerm=nPerm,gsSizeLim=c(min.gns.nu,max.gns.nu))
		#resList <- list(gsaRes1,gsaRes2)#,gsaRes3)	
		#names(resList) <- c("mean","median")
	
		# #file4path=comps[nn]#remove.sub(the.fl,".txt",part2take=1)
		# chd=consensusScores(resList,class="distinct",direction="down")
		# jj=which(chd$pMat[,1]<=pval.cutoff)
		# tmpD=cbind(rownames(chd$pMat),chd$pMat[,1])[jj,]
		# #write.dat(tmp,file.name=paste("pathways",gs.name,"_",name2compare,"_down.xls",sep=""))
		# if(is.matrix(tmpD))
			# res.gs=rbind(res.gs,cbind(tmpD,rep(-1,nrow(tmpD))))
		# if(is.vector(tmpD))
				# res.gs=rbind(res.gs,c(tmpD,-1))

		
		
 		
}
}
}

organize.pathways<-function(wd)
{
	#just to organize all pathways ###
	setwd(wd)   
	fls=dir()
	fls=remove.sub(fls,"Up_",part2take=2)
	fls=remove.sub(fls,"Dn_",part2take=2)
	fls=remove.sub(fls,"_c1",part2take=1)
	fls=remove.sub(fls,"_c2",part2take=1)
	fls=remove.sub(fls,"_c3",part2take=1)
	fls=remove.sub(fls,"_c4",part2take=1)
	fls=remove.sub(fls,"_c5",part2take=1)
	fls=remove.sub(fls,"_c6",part2take=1)
	fls=remove.sub(fls,"_c7",part2take=1)
	fls=remove.sub(fls,".txt",part2take=1)
	fls=remove.sub(fls,".xls",part2take=1)
	fls=unique(fls)

	for(the.fl in fls)
	{
	print(the.fl)
	
	nn=1
	res=NULL
	Res=NULL
	fls2use=gr("Dn",gr(the.fl,dir()))
	comp=remove.sub(remove.sub(fls2use[1],".xls",part2take=1),".txt",part2take=1)
	
	for(the.fl2use in fls2use)
	{
		print(the.fl2use)
		res=as.matrix(read.csv(the.fl2use,sep="\t",header=T))
	print(res)
		if(nrow(res)>0)
		{		if(ncol(res)==1)
						res=t(res)
				#print(res[1:2,1:3])
			if(nn>1)
				{
					Res=rbind(Res,res)	
					nn=nn+1
				}	
			if(nn==1)
				{
					Res=res
					nn=nn+1
				}			
				}	
	}
	if(!is.null(Res))
	{
		Res=as.matrix(Res)
		if(nrow(Res)>1)
		{
			Res=Res[order(Res[,3],decreasing=T),]
			Res=Res[which(!duplicated(Res[,1])),]
		}	
		write.dat(Res,pastec(comp,"_pathways.xls"))
	}
	
    
    nn=1
    res=NULL
    Res=NULL
	fls2use=gr("Up",gr(the.fl,dir()))
	comp=remove.sub(remove.sub(fls2use[1],".xls",part2take=1),".txt",part2take=1)
	
	for(the.fl2use in fls2use)
	{
		#print(the.fl2use)
		res=as.matrix(read.csv(the.fl2use,sep="\t",header=T))
	
		if(nrow(res)>0)
		{		if(ncol(res)==1)
						res=t(res)
			#print(res[1:2,1:3])
			if(nn>1)
				{
											
					Res=rbind(Res,res)	
					nn=nn+1
				}	
			if(nn==1)
				{
					Res=res		
					nn=nn+1
				}	
		#print(Res[1:2,1:3])
		
		}	
	}
	if(!is.null(Res))
	{
		Res=as.matrix(Res)
	if(nrow(Res)>1)
	{
		Res=Res[order(Res[,3],decreasing=T),]
		Res=Res[which(!duplicated(Res[,1])),]
	}
	write.dat(Res,pastec(comp,"_pathways.xls"))
	}
	
}
}
