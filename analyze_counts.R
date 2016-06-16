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
    if ("GeneSymbol" %in% colnames(HTSeq.dat)){
        gns = HTSeq.dat[,1]
        idsAndGns = as.matrix(gns)
        rownames(idsAndGns) = rownames(HTSeq.dat)
        HTSeq.dat = HTSeq.dat[,-1]
    }

    ## TEMP FOR KALLISTO TESTING
    ## filter out any rows with zero counts 
    #HTSeq.dat=HTSeq.dat[apply(HTSeq.dat, 1, function(row) all(row >10 )),]

    if(!is.null(key)){
        HTSeq.dat = HTSeq.dat[,key[,1]]
    }

    samps = colnames(HTSeq.dat)

    ## for kallisto counts, round floats to integers
    counts.dat=round(matrix2numeric(HTSeq.dat))

    ## if no conditions given, create a vector of one mock 
    ## condition needed for make.cds  
    if(is.null(conds)){
        conds=rep('s',length(samps))
    }
    ########################################
    ## write scaled data using DESeq method
    ########################################
    cat("    Making CDS...\n")

    cds <- tryCatch({
        cds <- make.cds(counts.dat=counts.dat,conds=conds,count.cut=count.cut,libsizeQ=libsizeQ,percentile=percentile,method=method)
        }, error = function(err){
            warning(paste("method='",method,"' did not work.. trying method='pooled'",sep=""))
            cds <- tryCatch({
                      cds=make.cds(counts.dat=counts.dat,conds=conds,count.cut=count.cut,libsizeQ=libsizeQ,percentile=percentile,method='pooled')
                     }, error = function(x){
                      stop(paste("Both methods '",method,"' and 'pooled' failed.",sep=""))
                   }
            )
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

    if(length(colnames(counts.scaled))<3){
        cat("    Less than three samples; can not run cluster analysis\n")
        return()
    }

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

make.generic.heatmaps <- function(diff.exp.dir,norm.counts.file){

    ## load normalized counts
    dat = read.delim(norm.counts.file,sep="\t",header=T)

    ## assign IDs as rownames
    rownames(dat) = dat[,1]
    
    if ("GeneSymbol" %in% colnames(dat)){
        idHeader = "GeneSymbol"
    } else {
        idHeader = "GeneID"
    }

    comp_results = dir(diff.exp.dir)[grep("^ResDESeq.*\\.xls$",dir(diff.exp.dir))]

    for (compres in comp_results){
        file = sub("ResDESeq_","",compres)
        file = sub(".xls","",file)
        conds = unlist(strsplit(file,"_vs_"))
        condA = conds[1]
        condB = conds[2]

        genes = as.matrix(read.delim(paste(diff.exp.dir,compres,sep="/"),sep="\t",header=T))[,1]

        ## extract data for those genes
        htmp.dat = dat[genes,grep(paste(idHeader,condA,condB,sep="|"),colnames(dat))]

        ## remove duplicate genes
        if(length(htmp.dat[which(duplicated(htmp.dat[,idHeader])),idHeader])>0){
            htmp.dat = htmp.dat[-which(duplicated(as.vector(htmp.dat[,idHeader]))),]
        }

        ## remove any rows that have NAs
        htmp.dat = htmp.dat[complete.cases(htmp.dat),]
        ## as long as gene symbols are unique, assign them 
        ## as rownames; if not, we'll have to average values
        ## for genes that occur multiple times
        rownames(htmp.dat)=htmp.dat[,idHeader]
        htmp.dat = htmp.dat[,-1]

        ## from matrix
        htmp.dat = as.matrix(log2(htmp.dat))

        ## take only the top 100 genes in order to make heatmap legible
        if(dim(htmp.dat)[1]>100){
            htmp.dat = htmp.dat[1:100,]
        }

        if(dim(htmp.dat)[1]>1 && dim(htmp.dat)[2]>1){
            ## make heatmap pdf
            pdf(paste(diff.exp.dir,"/ResDESeq_",condA,"_vs_",condB,"_heatmap.pdf",sep=""),width=25,height=16)
            heatmap.2(htmp.dat - apply(htmp.dat, 1, mean), trace='none', col=colorpanel(16,"green","black","red"),cexRow=0.9,cexCol=1.2, dendrogram="both",main=paste("Top Differentially Expressed Genes ",condA," vs ",condB,sep=""), symbreaks=TRUE, keysize=0.5, margin=c(20,10))
            dev.off()
        }
    }

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
