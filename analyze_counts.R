normalize.counts <- function(counts.file,output.dir=output.dir,conds=conds,count.cut=count.cut,libsizeQ=libsizeQ,percentile=percentile,method=method,fitType=fitType,bin=bin,key=NULL,norm.quantile=FALSE){

    source(paste(bin,"tools.R",sep="/"))
    source(paste(bin,"run_DESeq.R",sep="/"))

    sink("/dev/null")
    library("DESeq")
    library("limma")
    sink()

    #########################################
    ## read and reformat raw counts file
    ## use unique identifiers as rownames (e.g., ensembl IDs) 
    ## but keep gene symbols (if given) for later use
    #########################################
    cat("    Reformatting raw counts...\n")
    HTSeq.dat=rr(counts.file,header=T)
    HTSeq.dat=make.rownames(HTSeq.dat)
    gns = NULL
cat(c(colnames(HTSeq.dat),"\n"))
    if ("GeneSymbol" %in% colnames(HTSeq.dat)){
        gns = HTSeq.dat[,1]
        idsAndGns = as.matrix(gns)
        rownames(idsAndGns) = rownames(HTSeq.dat)
        HTSeq.dat = HTSeq.dat[,-1]
    }

    if(!is.null(key)){
        HTSeq.dat = HTSeq.dat[,key[,1]]
    }

cat(c(colnames(HTSeq.dat),"\n"))

    samps = colnames(HTSeq.dat)
    counts.dat=matrix2numeric(HTSeq.dat)

    ## if no conditions given, create a vector of one mock 
    ## condition needed for make.cds  
    if(is.null(conds)){
        conds=rep('s',length(samps))
    }
cat(conds)
cat("\n")
    ########################################
    ## write scaled data using DESeq method
    ########################################
    cat("    Making CDS...\n")

    cds <- tryCatch({
        cds=make.cds(counts.dat=counts.dat,conds=conds,count.cut=count.cut,libsizeQ=libsizeQ,percentile=percentile,method=method)
        }, error = function(err){
            warning(paste("method='",method,"' did not work.. trying method='pooled'",sep=""))
            cds=make.cds(counts.dat=counts.dat,conds=conds,count.cut=count.cut,libsizeQ=libsizeQ,percentile=percentile,method='pooled')
            return(cds)
        }
    )

    cat("    Getting DESeq scaled counts...\n")
    counts.scaled=counts(cds,norm=T)

    #cat(c("        Dim counts.scaled:",dim(counts.scaled),"\n"))
    counts.log.dat=log2(counts.scaled+1)
    dat=2^counts.log.dat
    if (!is.null(key)){
        #cat(c("        Key is not null\n"))
        key[,1] = make.names(key[,1]) ## remove any invalid characters
        dat = dat[,key[,1]]
        colnames(dat)=paste(key[,1],key[,2],sep="__")
    }
    if (exists("idsAndGns")){
        idsAndGns = as.matrix(idsAndGns[rownames(counts.scaled),])
        #cat(c("        'gns' exists and is not null\n"))
        #cat(c("        Dim 'idsAndGns'",dim(idsAndGns),"\n"))
        #cat(c("        Dim 'dat'",dim(dat),"\n"))
        dat=as.matrix(cbind(rownames(dat),idsAndGns,dat))
        colnames(dat)[1]="GeneID"
        colnames(dat)[2]="GeneSymbol"
    } else {
        #cat(c("        'idsAndgns' does not exist or is null\n"))
        gns = c()
        dat=as.matrix(cbind(rownames(dat),dat))
        colnames(dat)[1]="GeneSymbol"
    }

    cat("    Writing DESeq scaled counts file...\n")
    setwd(output.dir)
    write.dat(dat,file.name="counts_scaled_DESeq.xls")

    ######################################### 
    ## if specified, run quantile method of normalization
    ######################################### 
    if(norm.quantile){
        cat("    Getting quantile-normalized counts...\n")
        counts.log.dat=log2(counts.dat+1)
        counts.log.norm.dat=normalizeBetweenArrays(counts.log.dat,method='quantile')
        dat=2^counts.log.norm.dat
        if (exists("gns")){
            dat=as.matrix(cbind(rownames(dat),gns,dat))
            colnames(dat)[1]="GeneID"
            colnames(dat)[2]="GeneSymbol"
        } else {
            dat=as.matrix(cbind(rownames(dat),dat))
            colnames(dat)[1]="GeneSymbol"
        }
        cat("    Writing quantile-normalized counts...\n")
        write.dat(dat,file.name="counts_scaled_quantile.xls")
    } 

    ## counts.scaled will be used for clustering and 
    ## counts.dat will be used differential expression
    counts = list()
    counts$scaled = counts.scaled
    counts$raw = counts.dat
    counts$gns = gns
    return(counts)

}


cluster.samples <- function(counts.scaled,output.dir=output.dir,conds=conds){

    setwd(output.dir)

    if(is.null(conds)){
        conds=rep('s',length(colnames(counts.scaled)))
    }

    counts2hclust=log2(counts.scaled+1)

    ########################################
    ## Hierarchical clustering
    ########################################
    sink("/dev/null")
    pdf.hclust(counts2hclust,file.name="counts_DESeqscaled_clust.pdf",title="All counts scaled using DESeq method")
    sink()

    ########################################
    ###MDS clustering
    ########################################
    md=cmdscale(dist(t(counts2hclust)),2)
    pdf("MDS_plot.pdf",width=18,height=12)
    plot(md, col=as.factor(conds),main="",lwd=2.5,cex=1.5)
    legend(-200, 70, levels(as.factor(conds)),col=as.factor(levels(as.factor(conds))),pch=1,cex=1.2)
    text(md,colnames(counts2hclust),cex=0.9)
    dev.off()

    return

}


run.diff.exp <- function(counts.raw,comps,key,output.dir,gns=NULL,q.cut=0.5,lfc=1,fc.cut=2,count.cut=15,zeroaddQ=T,libsizeQ=F,fitType='parametric',method='per-condition',sharingMode='maximum'){

    setwd(output.dir)
    AllRes=vector(length=length(comps),mode="list")

    for (nn in 1:length(comps)){

        comp = comps[nn]
        cat(c("    ",comp,"\n"))
        condA = unlist(strsplit(comp," - "))[1]
        condB = unlist(strsplit(comp," - "))[2]

        conds = key[grep(paste(condA,condB,sep="|"),key[,2]),2]

        dat = counts.raw[,grep(paste(condA,condB,sep="|"),key[,2])]

        ResDESeq=run.DESeq(counts.dat=dat, conds=conds, condA=condA, condB=condB, q.cut=q.cut, fc.cut=fc.cut, count.cut=count.cut,zeroaddQ=zeroaddQ,libsizeQ=F,fitType=fitType,sharingMode=sharingMode,method=method,gns=gns)
        AllRes[[nn]]=ResDESeq

    }

    save(AllRes,comps,file="All_Results.RData")

    return
}


run.gene.set.analysis <- function(species,bin,gsa.dir,deseq.res.dir,min.gns.nu=5,max.gns.nu=1000,pval.cutoff=0.1,nPerm=1e4,fcQ=T,fc2keep=log2(1.5),frac2keep=4){

    if (species %in% c("hg19","hg18")){
        species = "human"
    }
    if (species %in% c("mm9","mm10")){
        species = "mouse"
    }
    if (!species %in% c("human","mouse")){
        cat("    Error: Only human and mouse data supported. Can not run Gene Set Analysis\n")
        q()
    }

    tmp<-capture.output(suppressMessages(library("piano")))
    tmp<-capture.output(suppressMessages(library("MASS")))
    tmp<-capture.output(suppressMessages(source(paste(bin,"/tools.R",sep=""))))
    tmp<-capture.output(suppressMessages(source(paste(bin,"/run_DESeq.R",sep=""))))
    tmp<-capture.output(suppressMessages(source(paste(bin,"/script2runGSApiano.R",sep=""))))
    gmt.dir=paste(bin,"/data",sep="")


    if(species == "human"){
        gs.names.list=c("c1.all.v4.0.symbols.gmt","c2.all.v4.0.symbols.gmt","c3.all.v4.0.symbols.gmt","c5-1.all.v4.0.symbols.gmt", "c6-1.all.v4.0.symbols.gmt","c7.all.v4.0.symbols.gmt")[-1]
    } else {
        gs.names.list=c("mouse_c1.gmt","mouse_c2.gmt","mouse_c3.gmt","mouse_c4.gmt","mouse_c5.gmt")
    }


    pd=getwd()
    setwd(deseq.res.dir)
    fls=grep("ALL",dir(),value=T)
    AllGSARes = list()
    for (gs in gs.names.list){
        gs.name = remove.sub(gs,"\\.gmt",part2take=1)
        gsc = loadGSC(paste(gmt.dir,gs,sep="/"))

        for(fl in fls){
            comp.name=remove.sub(fl,"ALLResDESeq_",part2take=2)
            res=as.matrix(read.csv(fl,sep="\t"))
            res=res[,-1] #remove ensemble IDs ###
            res=make.rownames(res)
            res=matrix2numeric(res)
            fc=as.matrix(res[,3])
            rownames(fc)=rownames(res)
            pvals=as.matrix(res[,2])
            rownames(pvals)=rownames(res)

            setwd(gsa.dir)
            res.gs=c("GeneSets","Pval","Direction")
            gsaRes <- runGSA(fc,geneSetStat="mean",gsc=gsc,gsSizeLim=c(min.gns.nu,max.gns.nu),nPerm=nPerm)
            AllGSARes[[comp.name]][[gs.name]] = gsaRes
            tab.res=processgsaRes(gsaRes,pval.cutoff=pval.cutoff,fc2keep=fc2keep,frac2keep=frac2keep,fcQ=fcQ)
            file.name=paste("GeneSet_Dn_",remove.sub(remove.sub(fl, "ALLResDESeq_",part2take=2),"_Genes",part2take=1),"_",gs.name,".xls",sep="")
            write.dat(tab.res$gsa.tab.dn,file.name=file.name)

            file.name=paste("GeneSet_Up_",remove.sub(remove.sub(fl, "ALLResDESeq_",part2take=2),"_Genes",part2take=1),"_",gs.name,".xls",sep="")
            write.dat(tab.res$gsa.tab.up,file.name=file.name)
            setwd(deseq.res.dir)
        }
    }
    setwd(gsa.dir)
    save(AllGSARes,file="AllGSARes.Rdata",compress=T)

    setwd(pd)
    return
}
