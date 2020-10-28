#' Read in raw HTSeq count file and format for analysis
#'
#' Reformatting includes removing samples that are not in 
#' sample key (if key is provided), rounding counts to
#' nearest integer, assigning GeneIDs to rownames of matrix,
#' and removing GeneID column from matrix if it also contains
#' a GeneSymbol column.
#'
#' @param htseq.files tab-delimited file containing raw HTseq
#'                    counts; columns are samples and rows are
#'                    genes; must contain GeneID column but may 
#'                    also contain a GeneSymbol column
#' @param key         two-column matrix where column one contains
#'                    all samples to be included in analysis and
#'                    column two contains condition groups to which
#'                    the samples belong
#' @return a list of two items: 1) a matrix of raw counts formatted 
#'         for BIC differential expression analysis and 2) a vector
#'         of gene symbols named by the GeneID in htseq file (this
#'         vector will be NULL if there is no GeneSymbol column in
#'         the htseq file.
#' @export
bic.format.htseq.counts <- function(htseq.file,key=NULL){

  cat("Reformatting raw counts...")

  idsAndGns <- matrix()

  HTSeq.dat=bic.file2matrix(htseq.file,header=T, sep="\t")
  HTSeq.dat=bic.make.rownames(HTSeq.dat)

  if ("GeneSymbol" %in% colnames(HTSeq.dat)){
    idsAndGns <- as.matrix(HTSeq.dat[,1])
    rownames(idsAndGns) <- rownames(HTSeq.dat)
    HTSeq.dat <- HTSeq.dat[,-1]
  }

  ## subset/sort data according to key
  if(!is.null(key)){
    smps <- key[,1]
    if(any(!smps %in% colnames(HTSeq.dat))){
        print(paste0("WARNING: These samples are NOT in count data... ", 
                     paste(smps[!smps %in% colnames(HTSeq.dat)], collapse=", ")))
        smps <- smps[smps %in% colnames(HTSeq.dat)] 
    }
    HTSeq.dat = HTSeq.dat[,smps]
  }

  HTSeq.dat = HTSeq.dat[!rownames(HTSeq.dat) %in% c("alignment_not_unique","ambiguous","no_feature","not_aligned","too_low_aQual"),]

  ## for kallisto counts, round floats to integers
  raw.counts=round(bic.matrix2numeric(HTSeq.dat))
  cat("Done.\n")

  formatted.counts <- list()
  formatted.counts$raw <- raw.counts
  formatted.counts$ids <- idsAndGns

  return(formatted.counts)
}

#' Non-deseq way to normalize counts(?)
#'
#' TO DO: INSERT DESCRIPTION HERE
#' 
#' @param raw.counts  Raw count matrix, where a column is one sample and a row is one gene
#' @param percentile  A string representing percentile; e.g., "100\%". Default: "75\%"
#' @return TO DO
#' @export
bic.scale.factor.quantile <- function(raw.counts,percentile="75%"){
  q=quantile(raw.counts,probs=seq(0,1,0.05))
  q=q[grep(percentile,names(q))]
  return(sum(raw.counts[which(raw.counts<=q)]))
}

#' Run a differential expression comparison between two
#' conditions.
#'
#' Compare expression in condition B vs expression in condition B. Results include
#' log2FC, adjusted p-val, and mean expression in each condition. Two result sets
#' are included in returned list: one including ALL genes, and one including only 
#' genes differentially expressed as determined by \code{max.p} and \code{min.abs.fc}. 
#' User has the option to include in the results those genes with zero counts in one
#' condition and count greater than \code{min.count/mean(sizeFactors)} in the other
#' condition, even if the p-value is greater than \code{max.p}, as those genes may still
#' be of interest.
#'
#' @param countDataSet  A countDataSet object created by running DESeq's newCountDataSet()
#' @param conds         Vector of sample conditons, in the same order as samples 
#'                      in raw count matrix
#' @param condA         The first condition in comparison (must be in conditions vector)
#' @param condB         The second condition in comparison (must be in conditions vector)
#' @param max.p         A number, the max p-value cutoff to determine which genes are 
#'                      considered differentially expressed; Default: 0.05
#' @param min.abs.fc    A number, the minimum absolute fold change cutoff to determine 
#'                      which genes are considered differentially expressed; Default: 2
#' @param min.count     A number, minimum average number of reads in one condition; 
#'                      Default: 10
#' @param zeroaddQ      Logical, indicates whether to include genes that have zero 
#'                      counts in one condition and insignificant p-value, but at 
#'                      least min.count/mean(sizeFactors) in the other condition; 
#'                      Default: FALSE
#' @param genes         A named vector of gene symbols, where the names are equivalent 
#'                      to the IDs in the original raw count matrix (e.g., Ensembl IDs); 
#'                      Default: NULL
#' @return List containing a matrix of ALL results, a matrix containing results
#'         for differentially expressed genes (as determined by cutoffs), a vector
#'         of DE genes, the cutoffs used in analysis, and norm factors (TO DO: reword this) 
#' @export
bic.run.deseq.comparison <- function(countDataSet, conds, condA, condB,
                                     max.p=0.05, min.abs.fc=2, min.count=10,
                                     zeroaddQ=F, genes=c()){

  cat(paste0("Comparing condB [", condB, "] vs [", condA, "]"))

  cds <- countDataSet
  ng <- NULL ## num genes
  ans <- NULL ## filtered results

  
  res <- tryCatch({ 
             nbinomTest(cds, condA, condB) 
           }, error = function(e){
             print(e)
             return(NULL)
         })
  if(is.null(res)){
      return(NULL)
  }
  DESeqRes <- res

  ## add gene symbol to result in case 'id' is ensembl id (or other id)
  rownames(res) <- res[,1]
  if(length(genes) > 0){
    genes <- genes[rownames(res),1]
    res <- cbind(res,genes)
    colnames(res)[length(colnames(res))] = "GeneSymbol"
  }

  ng <- rownames(res)[which(res$padj <= max.p & abs(res$log2FoldChange)>=log2(min.abs.fc))]
  norm.factors <- sizeFactors(cds)

  ## add genes that have zero counts in one condition but 
  ## count greater than min.count/mean(sizeFactors) in
  ## in another condition but 
  
  if(zeroaddQ){
    jj1 <- rownames(res)[which(res[,"baseMeanA"]==0 & res[,"baseMeanB"]>=min.count/mean(norm.factors))]
    jj2 <- rownames(res)[which(res[,"baseMeanB"]==0 & res[,"baseMeanA"]>=min.count/mean(norm.factors))]
    ng <- unique(c(ng,jj1,jj2))
  }

  if(length(ng) > 0){
    resSig <- res[ng,]
    jj <- unique(rownames(resSig)[which(resSig[,"baseMeanA"]>= min.count/mean(norm.factors) | 
                                        resSig[,"baseMeanB"]>= min.count/mean(norm.factors))])
    if(length(jj) == 0){
      ng <- NULL
    } else {
      ng <- ng[ng %in% jj]
    }
  }

  ## prepare filtered output
  if(length(ng) > 0){
    if ("GeneSymbol" %in% colnames(resSig)){
      ans=resSig[ng,c("id","GeneSymbol","padj", "log2FoldChange",
                      "baseMeanA","baseMeanB")]
      colnames(ans)=c("GeneID","GeneSymbol", "P.adj",
                      paste("log2[",condB,"/",condA,"]",sep=""),
                      paste("Mean_at_cond_",condA, sep=""),
                      paste("Mean_at_cond_",condB, sep=""))
    } else {
      ans=resSig[ng,c("id","padj", "log2FoldChange","baseMeanA","baseMeanB")]
      colnames(ans)=c("GeneID", "P.adj", paste("log2[",condB,"/",condA,"]",sep=""), 
                      paste("Mean_at_cond_",condA, sep=""), 
                      paste("Mean_at_cond_",condB, sep=""))
    }
    ## sort by log2FC
    ans = ans[order(abs(ans[,paste0("log2[",condB,"/",condA,"]")]), decreasing = TRUE),]
    rownames(ans)=ans[,1]
  } else {
    cat(paste0("\n===================================================\nNo genes pass the significant cutoff of ", max.p, "\n  AND have sufficient mean number of reads across samples\n AND at least ", min.count, " reads in one condition\n===================================================\n"))
    #m=max(res[which(res$padj==min(res$padj)),"baseMeanA"],
    #      res[which(res$padj==min(res$padj)),"baseMeanB"])
    #cat("Best corrected p.value =",min(res$padj),"\n")
    #cat("with max number of counts", m, "\n")
    cat("\n\n")
  }

  #############
  ## prepare unfiltered results (results for all genes)
  if ("GeneSymbol" %in% colnames(res)){
    res=res[,c("id","GeneSymbol","pval","padj", "log2FoldChange","baseMeanA","baseMeanB")]
    colnames(res)=c("GeneID","GeneSymbol", "pval","P.adj", 
                    paste("log2[",condB,"/",condA,"]",sep=""), 
                    paste("Mean_at_cond_",condA, sep=""), 
                    paste("Mean_at_cond_",condB, sep="")
                   )
  } else {
    res=res[,c("id","pval","padj", "log2FoldChange","baseMeanA","baseMeanB")]
    colnames(res)=c("GeneID", "pval","P.adj", 
                     paste("log2[",condB,"/",condA,"]",sep=""), 
                     paste("Mean_at_cond_",condA, sep=""), 
                     paste("Mean_at_cond_",condB, sep="")
                   )
  }
  rownames(res)=res[,1]

  ## make sure there are no Inf FC
  jj <- which(abs(res[,paste("log2[",condB,"/",condA,"]",sep="")]) == Inf)
  if(length(jj) >= 1){
    res[jj, paste("log2[",condB,"/",condA,"]",sep="")] <-
          log2(as.numeric(res[jj,paste("Mean_at_cond_",condB, sep="")]) + 1) -
          log2(as.numeric(res[jj,paste("Mean_at_cond_",condA, sep="")]) + 1)
  }

  ## Keep only those that have at least min.count normalized mean counts in at least one group
  jj <- unique(which(res[,paste("Mean_at_cond_",condA, sep="")] >= min.count |
                     res[,paste("Mean_at_cond_",condB, sep="")] >= min.count
                    )
              )
  res <- res[jj,]


  ######################
  ##### clean up results filtered by FC and pval

  ## make sure there are no Inf FC
  if(!is.null(ans)){
    jj=which(abs(ans[,paste("log2[",condB,"/",condA,"]",sep="")])==Inf)  
    if(length(jj)>=1){
        ans[jj, paste("log2[",condB,"/",condA,"]",sep="")] <- 
              log2(as.numeric(ans[jj,paste("Mean_at_cond_",condB, sep="")]) + 1) - 
              log2(as.numeric(ans[jj,paste("Mean_at_cond_",condA, sep="")]) + 1)
    }
  }

  ## DO NOT REMOVE THOSE WITH LESS THAN MIN COUNT AS WE SAVED ZERO-COUNT
  ## GENES BEFORE, IN CASE THEY ARE OF INTEREST TO THE INVESTIGATOR

  Res=list()
  Res$DESeq         <- DESeqRes 
  Res$filtered      <- ans
  Res$DEgenes       <- rownames(ans)
  Res$all.res       <- res
  Res$max.p         <- max.p
  Res$min.abs.fc    <- min.abs.fc
  Res$min.count     <- min.count
  Res$norm.factors  <- norm.factors

  return(Res)
  
} 



#' Normalize raw count matrix using DESeq scaling method
#' 
#' @param raw.counts  Raw count matrix, where each column is a sample and each
#'                    row is a gene 
#' @param conds       vector of sample conditions, in the same order as
#'                    column names in \code{raw.counts}
#' @param min.count   minimum total number of reads for a gene across samples; Default: 10
#' @param libsizeQ    logical to indicate we should use library size to normalize counts 
#'                    rather than DESeq's method [TO DO: RECHECK THIS]
#' @param percentile  string indicating percentile to use when normalizing by library size
#' @param fitType     See DESeq documentation. Default: "parameteric"
#' @param method      See DESeq documentation. Default: "per-condition"
#' @param sharingMode See DESeq documentation. Default: "maximum"
#' @return DESeq's countDataSet object
#' @export
bic.get.deseq.cds <- function(raw.counts, conds,
                           min.count=10, libsizeQ=F, percentile="100%",
                           fitType="parametric", method="per-condition",
                           sharingMode="maximum"){

  cat("Getting DESeq countDataSet...")
  ## filter raw counts: remove genes with total 
  ## counts less than min.count
  x <- apply(raw.counts,1,sum)
  counts.tmp <- raw.counts[which(x >= min.count),]

  ## create CountDataSet, DESeq's main data structure
  cds <- newCountDataSet(counts.tmp,conds)

  ## by default, normalize counts using DESeq's method
  ## if otherwise specified, normalize by library size
  if(!libsizeQ){
    cds = estimateSizeFactors(cds)
  }

  if(libsizeQ | length(which(is.na(sizeFactors(cds))))>0){
    cat("\nEstimating sizeFactors using library size\n")
    if(percentile == "100%"){
      libsizes <- apply(counts.tmp,2,sum)
    } else {
      libsizes <- apply(counts.tmp,2,util.scale.factor.quantile,percentile=percentile)
    }
    libsizes <- libsizes / median(libsizes)
    sizeFactors(cds) <- libsizes
  }
  
  cds <- tryCatch({
           estimateDispersions(cds,fit=fitType,method=method,sharingMode=sharingMode)
         }, error = function(err){
           if(!method == "pooled"){
             warning(paste("method='",method,"' did not work. Trying method='pooled'",sep=""))
             cds <- tryCatch({
                  estimateDispersions(cds,fit=fitType,method='pooled',sharingMode=sharingMode)
               }, error = function(x){
                  if(!method == "blind"){
                    warning(paste("method='",method,"' did not work. Trying method='blind'",sep=""))
                    cds <- tryCatch({
                         estimateDispersions(cds,fit=fitType,method='blind',sharingMode='fit-only')
                    }, error = function(x){
                         warning(paste("All methods '",method,"', 'pooled' and 'blind' failed.",sep=""))
                         #quit(save="no",status=15,runLast=TRUE)
                         return(NULL)
                    })
                  }
               })
           } else {
             warning("method 'pooled' did not work. Please try another method.")
             return(NULL)
             #quit(save="no",status=15,runLast=TRUE)
           }   
         })
  cat("Done.\n")
  return(cds) 
}


#' Write normalized counts matrix to file
#'
#' @param  norm.counts.mat     matrix containing normalized counts
#' @param  file.name           File to which results should be written
#' @export
bic.write.normalized.counts.file <- function(norm.counts.mat,file.name){
    write.table(norm.counts.mat,file=file.name,quote=F,col.names=T,row.names=F,sep="\t")
    #openxlsx::write.xlsx(as.data.frame(norm.counts.mat), file=file.name)
}

#' Write results for ALL genes and for DE genes
#'
#' Two files are written: One containing log2FC, adjusted p-val, mean(condA),
#' mean(condB) for ALL genes, and one containing all of that same information
#' for ONLY DE genes. Results are sorted by abs(log2FC) by default but may be
#' sorted by adjusted p-val. 
#'
#' @param res         Matrix to write to file, either all.res or filtered, 
#'                    from \code{bic.get.norm.deseq.counts}
#' @param file.name   File to which results should be written
#' @param orderPvalQ  Logical indicating results should be sorted by pValue instead of FC 
#' @return Nothing.   Files are written.
#' @export
bic.write.deseq.results <- function(res, file.name, orderPvalQ=FALSE){

  head(res)

  ## Sort results  
  if(orderPvalQ){
    res <- res[order(res[,"P.adj"],decreasing=F),]
  } else {
    res <- res[order(abs(res[,grep("log",colnames(res))]),decreasing=T),]
  }
  bic.write.dat(res, file=file.name) 
}

#' Use DESeq to normalize raw HTSeq counts 
#' 
#' Use DESeq scaling method to normalize raw HTSeq counts and return a list containing both 
#' raw and normalized counts. Takes in either formatted.counts matrix from bic.format.htseq.counts()
#' OR an HTSeq counts file. Optionally takes in previously generated countDataSet. If not given,
#' also optionally takes in parameters for generating one. Sample key will be used in generation of
#' countDataSet but if not given, it will be assumed that all samples belong to one condition group.
#'
#' @param formatted.counts  matrix containing raw counts from HTSeq. One column per sample,
#'                          one row per gene; may be used INSTEAD of htseq.file, but NOT in
#'                          conjunction with
#' @param htseq.file        file containing raw HTSeq counts; One column per sample, one row
#'                          per gene; may be used INSTEAD of formatted counts, but NOT in 
#'                          conjunction with
#' @param cds               DESeq countDataSet object; Default: NULL, if null will compute
#'                          on the fly
#' @param key               sample key where column one contains sample ID that matches
#'                          that in counts file, and column two contains the condition 
#'                          group to which it belongs. Only those samples that are in the 
#'                          key matrix will be included in the returned matrices
#' @param min.count         minimum total number of reads for a gene across samples; Default: 10
#' @param libsizeQ          logical to indicate we should use library size to normalize counts 
#'                          rather than DESeq's method [TO DO: RECHECK THIS]
#' @param percentile        string indicating percentile to use when normalizing by library size
#' @param fitType           See DESeq documentation. Default: "parameteric"
#' @param method            See DESeq documentation. Default: "per-condition"
#' @param sharingMode       See DESeq documentation. Default: "maximum"#' 
#' @return list containing two matrices: one of raw counts and one of
#'         DESeq-scaled, log2 transformed counts. Column names will have condition appended
#' @export
bic.deseq.normalize.htseq.counts <- function(formatted.counts=NULL,htseq.file=NULL,cds=NULL,key=NULL,
                                             min.count=10, libsizeQ=F, percentile="100%",
                                             fitType="parametric", method="per-condition",
                                             sharingMode="maximum"){

  if(!is.null(formatted.counts) & !is.null(htseq.file)){ 
    warning("Both formatted.counts and htseq.file are not NULL. Using htseq file for raw counts")
    formatted.counts=NULL
  }  

  if(is.null(formatted.counts)){
    formatted.counts <- bic.format.htseq.counts(htseq.file,key)
  }
  counts.dat <- formatted.counts$raw
  idsAndGns <- formatted.counts$ids
  
  if(is.null(cds)){
    if(is.null(key)){
      warning("No sample key given. Creating countDataSet using one condition for all samples.")
      nsamp <- ncol(formatted.counts$raw)
      if(length(grep("Gene",colnames(formatted.counts$raw)))>0){
        nsamp <- ncol(formatted.counts$raw) - length(grep("Gene",colnames(formatted.counts$raw)))
      }
      conds <- rep("s",nsamp)
    }
    cds <- bic.get.deseq.cds(formatted.counts$raw, conds, min.count=min.count, libsizeQ=libsizeQ, 
                             percentile=percentile, fitType=fitType, method=method,
                             sharingMode=sharingMode)
  }

  cat("Getting DESeq scaled counts...")
  counts.scaled <- counts(cds,norm=T)
  cat("Done.\n")

  cat("Formatting scaled counts...")
  counts.log.dat <- log2(counts.scaled+1)
  dat <- 2^counts.log.dat-1
  if (!is.null(key)){
    key[,1] <- make.names(key[,1]) ## remove any invalid characters
    dat <- dat[,key[,1]]
    colnames(dat) <- paste(key[,1],key[,2],sep="__")
  }
  if (!is.null("idsAndGns")){
    idsAndGns <- as.matrix(idsAndGns[rownames(counts.scaled),])
    dat <- as.matrix(cbind(rownames(dat),idsAndGns,dat))
    colnames(dat)[1:2] <- c("GeneID","GeneSymbol")
  } else {
    dat <- as.matrix(cbind(rownames(dat),dat))
    colnames(dat)[1] <- "GeneSymbol"
  }
  cat("Done.\n")
 
  counts <- list()
  counts$scaled <- dat
  counts$raw <- counts.dat
  counts$ids <- idsAndGns
  return(counts) 
}


#' Use quantile normalization to normalize raw HTSeq counts 
#' 
#' Run \code{normalizeBetweenArrays()} function to normalize
#' raw HTSeq counts and return a list containing both raw
#' and normalized counts
#'
#' @param htseq.counts  matrix containing raw counts from HTSeq. One column per 
#'                      sample, one row per gene
#' @param key           matrix where column one contains sample IDs that matches 
#'                      those in counts file, and column two contains the condition 
#'                      groups to which they belong. Only those samples that are in 
#'                      the key matrix will be included in the normalization
#' @return list containing two matrices (one of raw counts and one of quantile 
#'         normalized counts) and a vector of gene symbols if one was included in 
#'         the raw counts matrix
#' @export 
bic.quantile.normalize.htseq.counts <- function(htseq.counts,key=NULL){
  ## read and reformat raw counts file
  ## use unique identifiers as rownames (e.g., ensembl IDs)
  ## but keep gene symbols (if given) for later use)
  cat("Reformatting raw counts...")
  HTSeq.dat <- bic.file2matrix(htseq.counts,header=T, sep="\t")
  HTSeq.dat <- bic.make.rownames(HTSeq.dat)

  gns <- NULL
  idsAndGns <- NULL

  if ("GeneSymbol" %in% colnames(HTSeq.dat)){
    gns <- HTSeq.dat[,1]
    idsAndGns <- as.matrix(gns)
    rownames(idsAndGns) <- rownames(HTSeq.dat)
    HTSeq.dat <- HTSeq.dat[,-1]
  }

  ## only keep columns for samples in key
  if(!is.null(key)){
    HTSeq.dat <- HTSeq.dat[,key[,1]]
  }

  samps <- colnames(HTSeq.dat)

  ## for kallisto counts, round floats to integers
  counts.dat <- round(bic.matrix2numeric(HTSeq.dat))

  cat("Done.\n")

  cat("Getting quantile-normalized counts...")
  counts.log.dat=log2(counts.dat+1)
  counts.log.norm.dat=normalizeBetweenArrays(counts.log.dat,method='quantile')
  dat=2^counts.log.norm.dat
  ## append gene symbols to matrix if there are any
  if (exists("gns")){
    dat=as.matrix(cbind(rownames(dat),gns,dat))
    colnames(dat)[1:2]=c("GeneID","GeneSymbol")
  } else {
    dat=as.matrix(cbind(rownames(dat),dat))
    colnames(dat)[1]="GeneSymbol"
  }
  cat("Done.\n")

  counts <- list()
  counts$normalized <- dat
  counts$raw <- counts.dat
  counts$gns <- gns
  return(counts)   
}


#' Get the gene sets that meet criteria set by cut offs
#'
#' @param gsaRes      results from \code{piano::runGSA()}
#' @param gsa.tab     BIC-ified GSA results
#' @param fc2keep     minimum absolute value fold change; Default: log2(1.5)
#' @param frac2keep   at least 1/frac2keep genes in set must have FC >= fc2keep
#' @param dir         direction ["Up"|"Dn"] of the results in gsaRes
#'                     
#' @return filtered gsa.tab
bic.get.gene.sets <- function(gsaRes,gsa.tab,fc2keep=log2(1.5),frac2keep=2,dir="Up"){
  tmp <- cbind(gsa.tab,rep("",nrow(gsa.tab)),rep("",nrow(gsa.tab)))
  colnames(tmp) <- c(colnames(gsa.tab),"GenesInGeneSet","GenesAndFC")
  geneSets2remove=NULL
  for(nn in 1:nrow(gsa.tab)){
    fc <- geneSetSummary(gsaRes, gsa.tab[nn,1])$geneLevelStats
    o <- order(abs(fc),decreasing=T)
    fc <- fc[o]
    gn.names <- rownames(as.matrix(fc))
    gn.names.str <- bic.join.strings(gn.names,",")
    gn.names.fc.str <- bic.join.strings(paste(gn.names,"=",fc),";")
    tmp[nn,ncol(gsa.tab)+1] <- gn.names.str
    tmp[nn,ncol(gsa.tab)+2] <- gn.names.fc.str

    ## Add a condition that at least half of the genes have to have 
    ## FC >= log2(1.5) ! 
    if(dir == "Dn"){ jj <- length(which(fc <= -fc2keep)) }
    if(dir == "Up"){ jj <- length(which(fc >= fc2keep)) }
    if(jj <= length(fc)/frac2keep){ geneSets2remove=c(geneSets2remove,nn) }
  }
  if(length(geneSets2remove) >= 1){ tmp <- tmp[-geneSets2remove,] }

  return(tmp)
}

#' Format and filter results from piano::runGSA()
#' 
#' @param gsaRes    result object from piano::runGSA()
#' @param max.p     maximum p-value to use as cutoff; Default: 0.1
#' @param fc2keep   minimum fold change
#' @param frac2keep at least 1/frac2keep genes in set must have FC >= fc2keep
#' @param fcQ       remove gene sets that do not meet cutoffs set 
#' 
#' @return list of two matrices: one with up-regulated gene sets and one with down-regulated
#'         gene sets
bic.process.gsa.res <- function(gsaRes,max.p=0.1,fc2keep=log2(1.5),frac2keep=2,fcQ=T)
{
  gsa.tab.dn <- NULL
  gsa.tab.up <- NULL
  tmp=as.matrix(GSAsummaryTable(gsaRes))
  if(!is.null(tmp)){
    gsa.tab.dn <- NULL
    gsa.tab.up <- NULL

    indxDn <- which(tmp[,"p adj (dist.dir.dn)"] <= max.p)
    indxUp <- which(tmp[,"p adj (dist.dir.up)"] <= max.p)
    colsDn <- c("Name","Genes (tot)", "Stat (dist.dir)",
             "p adj (dist.dir.dn)","Genes (down)")
    colsUp <- c("Name","Genes (tot)", "Stat (dist.dir)",
             "p adj (dist.dir.up)","Genes (up)")
    if(length(indxDn) >= 1){
      if(length(indxDn) == 1){
        tmpDn <- t(as.matrix(tmp[indxDn,colsDn]))
      } else {
        tmpDn <- as.matrix(tmp[indxDn,colsDn])
      }
      tmpDn <- as.matrix(tmpDn[order(abs(as.numeric(tmpDn[,"Stat (dist.dir)"])),decreasing=T),])
      if(length(indxDn) == 1){
        tmpDn <- t(tmpDn)
      }
      if(fcQ){
        gsa.tab.dn <- bic.get.gene.sets(gsaRes,gsa.tab=tmpDn,fc2keep=fc2keep,frac2keep=frac2keep,dir="Dn")
      } else {
        gsa.tab.dn <- tmpDn
      }
    }
    if(length(indxUp)>=1){
      if(length(indxUp) == 1){
        tmpUp <- t(as.matrix(tmp[indxUp,colsUp]))
      } else {
        tmpUp <- as.matrix(tmp[indxUp,colsUp])
      }
      tmpUp <- as.matrix(tmpUp[order(abs(as.numeric(tmpUp[,"Stat (dist.dir)"])),decreasing=T),])
      if(length(indxUp) == 1){
        tmpUp <- t(tmpUp)
      }
      if(fcQ){
        gsa.tab.up <- bic.get.gene.sets(gsaRes,gsa.tab=tmpUp,fc2keep=fc2keep,frac2keep=frac2keep,dir="Up")
      } else {
        gsa.tab.up <- tmpUp
      }
    }
  }
  tab.res <- list()
  tab.res$gsa.tab.dn <- gsa.tab.dn
  tab.res$gsa.tab.up <- gsa.tab.up
  return(tab.res)
}

#' Run gene set analysis using the package 'piano'
#'
#' Run gene set analysis for one pairwise comparison of sample groups. 
#' GSA is run using the R package 'piano'. For each gene set, means 
#' are calculated for the genes with positive fold changes and for
#' genes with negative fold changes. P-values and adjusted p-values 
#' are then calculated and gene sets with a significant number of
#' of up- or down-regulated genes are reported in the output. Gene Sets
#' included are from the MSigDB Collections:
#' 
#' C1: positional gene sets
#' C2: curated gene sets
#' C3: motif gene sets
#' C5: GO gene sets
#' C6: oncogenic signatures
#' C7: immunologic signatures
#' 
#' Note: Mouse gene sets were downloaded from 
#'       http://bioinf.wehi.edu.au/software/MSigDB/
#'
#' @param species   currently only human and mouse are supported
#' @param deseq.res the matrix of ALL DESeq results (unfiltered)
#'                  in BIC format (GeneID,GeneSymbol,pval,P.adj,
#'                  log2[condB/condA],Mean_at_cond_condA,
#'                  Mean_at_cond_condB
#' @param min.gns   minimum number of genes a set must have in 
#'                  order to be included in analysis; Default: 5
#' @param max.gns   maximum number of genes a set must have in
#'                  order to be included in analysis; Default: 1000
#' @param max.p     maximum p-value to use as cut off; Default: 0.1
#' @param nPerm     the number of times the genes are randomized;
#'                  Default: 1e4
#' @param fc2keep   fold change cutoff (see frac2keep); Default: 
#'                  log2(1.5) 
#' @param frac2keep at least \code{1/frac2keep} of genes in a set
#'                  must have \code{fold change >= abs(fc2keep)};
#'                  Default: 4
#' @param fcQ       filter gene sets based on fc2keep and frac2keep       
#'                  Default: TRUE
#' @return a list of two matrices: one with gene sets found to be
#'         significantly UP-regulated and one with gene sets found
#'         to be significantly DOWN-regulated
#' @export
bic.run.gsa <- function(species,deseq.res,gmt.dir,min.gns=5,max.gns=1000,
                        max.p=0.1,nPerm=1e4,fcQ=T,
                        fc2keep=log2(1.5),frac2keep=4){

  if (species %in% c("hg19","hg18","hybrid")){
    species = "human"
  }
  if (species %in% c("mm9","mm10")){
    species = "mouse"
  }
  if (!species %in% c("human","mouse")){
    cat("    Error: Only human and mouse data supported. 
                    Can not run Gene Set Analysis\n")
    q()
  }

  if(species == "human"){
    gs.names.list <- c("h.all.v7.1.symbols.gmt", "c2.all.v7.1.symbols.gmt",
                       "c3.all.v7.1.symbols.gmt","c5.all.v7.1.symbols.gmt",
                       "c6.all.v7.1.symbols.gmt","c7.all.v7.1.symbols.gmt")
  } else {
    gs.names.list <- c("Mm.h.all.v7.1.gmt", "Mm.c2.all.v7.1.gmt",
                       "Mm.c3.all.v7.1.gmt", "Mm.c5.all.v7.1.gmt",
                       "Mm.c6.all.v7.1.gmt", "Mm.c7.all.v7.1.gmt")
  }                 

  all.res <- list()
  all.up.res <- NULL
  all.dn.res <- NULL
  for (gs in gs.names.list){
    print(gs)
    gs.name <- bic.remove.sub(gs,"\\.gmt",part2take=1)
    gs.cat <- toupper(bic.remove.sub(bic.remove.sub(gs,".all.v4.0.symbols.gmt",part2take=1),"-1",part2take=1))
    gsc <- loadGSC(file.path(gmt.dir, gs)) 

    res <- deseq.res[,-1]  ## remove ensembl IDs; since this function only supports
                           ## human and mouse, we can be confident that the first
                           ## column always has ensembl IDs
    fc  <- as.matrix(res[,c(1,4)])
    fc <- bic.average.by.name(fc)

    gsa.res <- runGSA(fc,geneSetStat="mean",gsc=gsc,
                     gsSizeLim=c(min.gns,max.gns),nPerm=nPerm)

    tab.res <- bic.process.gsa.res(gsa.res,max.p=max.p,
                                   fc2keep=fc2keep,frac2keep=frac2keep,fcQ=fcQ)
    if(!is.null(tab.res$gsa.tab.up)){
      if(is.null(dim(tab.res$gsa.tab.up))){
          tab.res$gsa.tab.up = t(as.matrix(tab.res$gsa.tab.up))
      }
      up.res <- cbind(gs.cat,tab.res$gsa.tab.up)
      colnames(up.res)[1] <- "Gene Set Category" 
      if(is.null(all.up.res)){
        all.up.res <- up.res
      } else {
        all.up.res <- rbind(all.up.res,up.res)
      }
    }
    if(!is.null(tab.res$gsa.tab.dn)){
      if(is.null(dim(tab.res$gsa.tab.dn))){
        tab.res$gsa.tab.dn = t(as.matrix(tab.res$gsa.tab.dn))
      }
      dn.res <- cbind(gs.cat,tab.res$gsa.tab.dn)
      colnames(dn.res)[1] <- "Gene Set Category"
      if(is.null(all.dn.res)){
        all.dn.res <- dn.res
      } else {
        all.dn.res <- rbind(all.dn.res,dn.res)
      }
    }
  }
  all.res$up <- all.up.res
  all.res$dn <- all.dn.res
  return(all.res)
}
