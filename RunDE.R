#! /opt/common/CentOS_6/R/R-3.2.0/bin/R
.libPaths("/opt/common/CentOS_6/R/R-3.2.0/lib64/R/library")

usage <- function(){

    usage.str = "\nUsage: Rscript RunDE.R
    \"counts.file    = '[Required: absolute path to htseq counts file]'\"

    \"key.file       = '[Required for differential expression analysis: 
                        absolute path to key file]'\"
    \"comps          = [Required for differential expression analysis: vector 
                       containing all comparisons to be made based on conditions 
                       in key file. Must be in the following format: c('CondA - CondB',
                       'CondA - CondC','CondB - CondC')]\"
    \"species        = '[Required for gene set analysis: (hg19|human|mm9|mm10|mouse)] 
                        Note: only human or mouse currently supported; specific
                        build does not matter, as long as it is clearly human or mouse'\"

    \"Rlibs          = [Optional (default=NULL): path to R libraries\"
    \"pre            = [Optional (default='TEMP'): prefix for all output files]\"
    \"diff.exp       = [Optional (default=TRUE): run differential expression analysis]\"
    \"GSA            = [Optional (default=TRUE): run gene set analysis; if running GSA 
                       but not DESeq, BE SURE TO SET diff.exp.dir (see below) to point
                       to existing DESeq results directory.]\"
    \"heatmaps       = [Optional (default=TRUE): generate a heatmap for each comparison
                       showing relative expression of the top ~100 DE genes\"

    \"counts.dir     = [Optional (default='$PWD/counts_gene')]\"
    \"clustering.dir = [Optional (default='$PWD/clustering')]\"
    \"diff.exp.dir   = [Optional (default='$PWD/differentialExpression_gene')]: if 
                        running DESeq, this is where output will go; if NOT running 
                        DESeq, this is wherever DESeq results already exist. Must 
                        be correct in order for GSA to run properly.\"
    \"gsa.dir        = [Optional (default='$PWD/GSA')]\"

    \"no.replicates  = [Optional (default=F): T|F, automatically sets fitType, 
                        method, and sharingMode to accommodate comparison
                        of single samples]\"
    \"q.cut          = [Optional (default=0.05): insert description here]\"
    \"lfc            = [Optional (default=1): insert description here]\"
    \"fc.cut         = [Optional (default=2): Fold change cutoff]\"
    \"count.cut      = [Optional (default=15): minimum total reads required for a gene
                        to be included in analysis]\"
    \"fitType        = [Optional (default='parametric'): 'parametric'|'local', See DESeq docs]\"
    \"orderPvalQ     = [Optional (default=T): T|F, sort final results by adjusted Pvalue]\"

    \"method         = [Optional (default='per-condition'): 'pooled'|'per-condition'|'blind', 
                        See DESeq docs.]\"
    \"sharingMode    = [Optional (default='maximum'): 'maximum'|'fit-only'|'gene-est-only', 
                        See DESeq docs.]\"
    \"zeroaddQ       = [Optional (default=T): T|F, include genes with zero counts in one
                        condition even if their pValues are insignificant]\"
    \"libsizeQ       = [Optional (default=F): T|F, normalize data using library size]\"
    \"test           = [Optional (default=F): T|F, run test data provided with bicrnaseq package]\"
    \n\n"
    cat(usage.str)
}

cat(c("\n++++++++++++++++ BIC RNA-Seq Counts Analysis ++++++++++++++++\n\n"))

pd = getwd()

################################################################################
####   Set defaults
################################################################################
test              <- FALSE
counts.dir        <- file.path(pd,"counts_gene")
clustering.dir    <- file.path(pd,"clustering")
diff.exp.dir      <- file.path(pd,"differentialExpression_gene")
gsa.dir           <- file.path(pd,"GSA")
max.p             <- 0.05
lfc               <- 1   #0.57#log2(fc.cut)
min.abs.fc        <- 2   #2^0.57
min.count         <- 15
percentile        <- "100%"
fitType           <- "parametric"
orderPvalQ        <- T
method            <- "per-condition"
sharingMode       <- "maximum"
zeroaddQ          <- T
libsizeQ          <- F
no.replicates     <- F
key               <- NULL
conds             <- NULL
GSA               <- T
diff.exp          <- T
heatmaps          <- T
pre               <- "TEMP"

key.file          <- NULL
counts.file       <- NULL

norm.counts       <- NULL
cds               <- NULL
mds.labels        <- FALSE
Rlibs             <- NULL
 
################################################################################
####   get user input
################################################################################
args=(commandArgs(TRUE))

if(length(args)==0){
    usage()
    q()
}
for(i in 1:length(args)){
    eval(parse(text=args[i]))
}

if(!is.null(Rlibs)){
  .libPaths(c(Rlibs,.libPaths()))
}

diff.exp.fig.dir  <- file.path(diff.exp.dir,"figures")
diff.exp.rdat.dir <- file.path(diff.exp.dir,"Rdata")
gsa.rdat.dir      <- file.path(gsa.dir,"Rdata")

if(test){
  ## run everything to test
  counts.file <- system.file("extdata","htseq_counts.txt",package="bicrnaseq")
  key.file <- system.file("extdata","sample_key.txt",package="bicrnaseq")
  comps=c('GroupA - GroupB', 'GroupB - GroupC')
  GSA=TRUE
  species='human'
}

################################################################################
####   validate input
################################################################################
if (!exists("counts.file")){
    cat("Error: Please specify a counts file. See usage for details\n")
    q()
}
if (!file.exists(counts.file)){
    cat(c("Error: counts file",counts.file,"doesn't exist.\n"))
    q()
}
if (exists("key.file") && !is.null(key.file) && !file.exists(key.file)){
    cat(c("Error: key file",key.file,"doesn't exist.\n"))
    q()
}
if (GSA && exists("key.file") && exists("comps") && !exists("species")){
    cat ("Error: Gene set analysis is turned ON. Please either turn it OFF or specify which species this data is from\n")
    q()
}

################################################################################
####   only load lib if input is valid
################################################################################
cat("Loading bicrnaseq library...")
#tmp <- capture.output(suppressMessages(library(bicrnaseq,lib.loc=Rlibs)))
#bin <- "/ifs/work/byrne/pipelines/rnaseq_pipeline"
source(file.path(bin, "bicrnaseqR/source_bicrnaseq.R"))
#tmp <- capture.output(suppressMessages(library(tidyverse))

cat("Done.\n")


################################################################################
## if key file is given, create key 
## and check for replicates
################################################################################
#key = as.matrix(read.delim(key.file,header=F,strip.white=T,sep="\t"))

if (exists("key.file") && !is.null(key.file)){
    dir.create(diff.exp.dir,showWarnings=FALSE,mode="0755")
    dir.create(diff.exp.fig.dir,showWarnings=FALSE,mode="0755")
    dir.create(diff.exp.rdat.dir,showWarnings=FALSE,mode="0755")
    key = as.matrix(read.delim(key.file,header=F,strip.white=T,sep="\t"))
    ##remove samples to be excluded
    ex = grep("_EXCLUDE_",key[,2])
    if(length(ex)>0){
        key = key[-ex,]
    }
    
    key[,1] = make.names(key[,1])
    conds=key[,2]
    
    t=table(conds)
    for(cond in conds){
        if(t[names(t)==cond]<2){
            cat(c("No replicates found for condition \"",
                  cond,"\". Setting no.replicates=TRUE\n"))
            no.replicates=TRUE
        }
    }
} else {
    key = NULL
}


## format counts
formatted.counts <- bic.format.htseq.counts(counts.file,key=key)

## if no conditions given, create a vector of one mock condition, 
## needed for normalization
if(is.null(conds)){
  conds <- rep("s",length(colnames(formatted.counts$raw)))
  mds.labels <- TRUE
} else if(length(colnames(formatted.counts$raw)) < 20){
  mds.labels <- TRUE
}

## get DESeq countDataSet
cds <- bic.get.deseq.cds(formatted.counts$raw,
                         conds,
                         min.count = min.count,
                         libsizeQ = libsizeQ,
                         percentile = percentile,
                         fitType = fitType,
                         method = method,
                         sharingMode = sharingMode)
if(is.null(cds)){
  for(c in comps){
      condA <- unlist(strsplit(c, " - "))[1]
      condB <- unlist(strsplit(c, " - "))[2]
      file.name <- paste(diff.exp.dir,
                         paste("ResDESeq_",condB,"_vs_",condA,".xls",sep=""),
                         sep="/")
      file.name2 <- gsub("Res","ALLRes",file.name)
      cat("\nZERO RESULTS\n")
      res <- as.data.frame(matrix(NA, nrow=1, ncol=6))
      names(res) = c("GeneID","GeneSymbol", "P.adj",
                  paste("log2[",condB,"/",condA,"]",sep=""),
                  paste("Mean_at_cond_",condA, sep=""),
                  paste("Mean_at_cond_",condB, sep=""))
      write_csv(res, file.name)
      write_csv(res, file.name2)
  }
  quit(save="no", status=0)
}
save(cds,file=file.path(diff.exp.rdat.dir,"cds.Rdata"),compress=TRUE)

################################################################################
####                           clustering & QC                              ####
################################################################################
cat("Running clustering and QC...")
setwd(pd)
tmp <- capture.output(
         suppressMessages(
           bic.deseq.heatmap(cds,
                             file=file.path(clustering.dir,
                                            paste0(pre,"_heatmap_50_most_highly_expressed_genes.pdf")),
                             transform=TRUE,
                             num.gns=50)
         )
       )
tmp <- capture.output(
         suppressMessages(
           bic.sample.to.sample.distances(cds,
                                         conds,
                                         file=file.path(sub("/$","",clustering.dir),paste0(pre,"_heatmap_sample_to_sample_distances.pdf"))
           )
         )
       )
tmp <- capture.output(
         suppressMessages(
           bic.deseq.plot.pca(cds,
                              file=file.path(clustering.dir,paste0(pre,"_PCA.pdf"))
           )
         )
       )
tmp <- capture.output(
         suppressMessages(
           bic.plot.dispersion.estimates(cds,
                                         out.dir=clustering.dir,file.prefix=pre)))
cat("Done.\n")


################################################################################
####                  Normalize raw HTSeq counts (always)                   ####
################################################################################
setwd(pd)
dir.create(counts.dir,showWarnings=FALSE,mode="0755")
## normalize and format counts
norm.counts <- bic.deseq.normalize.htseq.counts(formatted.counts=formatted.counts, cds=cds, key=key) 
norm.counts.mat <- norm.counts$scaled
counts.raw <- norm.counts$raw
counts.ids <- norm.counts$ids

cat("Writing normalized counts to file...")
file.name <- paste(counts.dir,"counts_scaled_DESeq.xls",sep="/")
dat <- cbind(counts.ids,norm.counts.mat)
bic.write.normalized.counts.file(dat,file.name)
cat("Done.\n")

################################################################################
####                         Cluster samples (always)                       ####
################################################################################
setwd(pd)
if(exists("norm.counts.mat") && !is.null("norm.counts.mat")){

  cat("Clustering all samples...")
  file.name <- file.path(clustering.dir,paste0(pre,"_counts_scaled_hclust.pdf"))
  tryCatch({
      bic.hclust.samples(norm.counts.mat,file.name=file.name,conds=conds,title="All counts scaled using DESeq method")
    }, error = function(err){
      conditionMessage(err)
    })
  file.name <- file.path(clustering.dir,paste0(pre,"_counts_scaled_MDS.pdf"))
  bic.mds.clust.samples(norm.counts.mat,file=file.name,conds=conds,labels=mds.labels)
  cat("Done.\n")
} else {
  cat("ERROR: Can not find normalized counts matrix. Can not cluster samples.\n")
}

################################################################################
####        run differential expression analysis if comps and key given     ####
################################################################################
all_results <- list()
all_results$norm.counts <- norm.counts
all_results$cds <- cds
all_results$de <- list()

if(diff.exp == T){
  if(exists("comps") && !is.null(comps) && !is.null(key)){

    setwd(pd)
    if(GSA){
      dir.create(gsa.dir,showWarnings=FALSE,mode="0755")
      dir.create(gsa.rdat.dir,showWarnings=FALSE,mode="0755")
    }

    cat("Running differential expression analysis...\n")

    #load("/ifs/work/byrne/pipelines/rnaseq_pipeline/all_DE_results.Rdata")

    for (comp in comps){
      condA = unlist(strsplit(comp," - "))[1]
      condB = unlist(strsplit(comp," - "))[2]
      conds = key[grep(paste(condA,condB,sep="|"),key[,2]),2]

      ## run DESeq comparison
      cat(paste("  Comparing: ",condA," vs ",condB,"...", sep="") )
      de.res <- bic.run.deseq.comparison(cds, conds, condA, condB,
                                         max.p = max.p, 
                                         min.abs.fc = min.abs.fc, 
                                         min.count = min.count, 
                                         zeroaddQ = zeroaddQ, 
                                         genes = counts.ids)
      cat("Done.\n")

      all_results$de[[paste(condA,"_vs_",condB,sep="")]] <- de.res

      #de.res <- all_results$de[[paste(condA,"_vs_",condB,sep="")]]

      ##
      ## DESeq visualization
      ##
      cat("    Drawing DESeq plots...")
      bic.plot.ma(de.res$DESeq,
                  file=file.path(diff.exp.fig.dir,
                                 paste0(pre,"_MAplot_",condA,"_vs_",condB,".pdf")
                       )
                 )
      bic.pval.histogram(de.res$DESeq,
                         file=file.path(diff.exp.fig.dir,
                                        paste0(pre,"_pval_histogram_",condA,"_vs_",condB,".pdf",sep="")
                              )      
                        )
      cat("Done.\n")
      ##
      ## format and write DE results to file
      ##
      cat("    Writing results for DE genes to file...")
      file.name <- paste(diff.exp.dir,
                         paste("ResDESeq_",condB,"_vs_",condA,".xls",sep=""),
                         sep="/")
      if(is.null(de.res$filtered)){
          cat(paste0("\n    ZERO RESULTS AFTER FILTERING\n  "))
          write(paste(c("GeneID","GeneSymbol", "P.adj",
                      paste("log2[",condB,"/",condA,"]",sep=""),
                      paste("Mean_at_cond_",condA, sep=""),
                      paste("Mean_at_cond_",condB, sep="")),collapse="\t"),file=file.name,append=FALSE)
      } else {
          bic.write.deseq.results(de.res$filtered,file.name=file.name,orderPvalQ=orderPvalQ)
          cat("Done.\n")
      }
      ##
      ## write All DE results to file
      ##
      cat("    Writing results for ALL genes to file...")
      file.name <- paste(diff.exp.dir,
                         paste("ALLResDESeq_",condB,"_vs_",condA,".xls",sep=""),
                         sep="/")
      if(is.null(de.res$all.res)){
          cat(paste0("\n    ZERO RESULTS GENERATED BY DESEQ\n  "))
          write(paste(c("GeneID","GeneSymbol", "P.adj",
                      paste("log2[",condB,"/",condA,"]",sep=""),
                      paste("Mean_at_cond_",condA, sep=""),
                      paste("Mean_at_cond_",condB, sep="")),collapse="\t"),file=file.name,append=FALSE)
      } else {
          bic.write.deseq.results(de.res$all.res,file.name=file.name,orderPvalQ=orderPvalQ)
          cat("Done.\n")
      }
      ##
      ## heatmap of top DE genes
      ##
      if(heatmaps & !is.null(de.res$filtered) & length(rownames(de.res$filtered)) > 0){
        genes <- de.res$DEgenes
        out.file <- file.path(diff.exp.fig.dir,
                              paste(pre,"_heatmap_",condA,"_vs_",condB,".pdf",sep="")
                    )
        tryCatch({
            bic.standard.heatmap(norm.counts.mat,condA,condB,genes=genes,file=out.file)
          }, error = function(e){
              warning(paste0("Could not generate heatmap for ",condA," vs ",condB))
          })
      }

      ##
      ## gene set analysis
      ##
      if(GSA){
        if(exists("species") & !is.null(species)){
          gsa.res <- bic.run.gsa(species,de.res$all.res) 
          if(!is.null(gsa.res$dn)){
            out.file <- file.path(gsa.dir,
                              paste("GeneSet_Dn_",condB,"_vs_",condA,".xls",sep="")
                        )
            bic.write.dat(gsa.res$dn,file=out.file)
          }
          if(!is.null(gsa.res$up)){
            out.file <- file.path(gsa.dir,
                              paste("GeneSet_Up_",condB,"_vs_",condA,".xls",sep="")
                        )
            bic.write.dat(gsa.res$up,file=out.file)
          }
          all_results$gsa[[paste(condA,"_vs_",condB,sep="")]] <- gsa.res
        }
      }
    }
    #save(all_results,file=file.path(diff.exp.rdat.dir,"all_DE_results.Rdata"),compress=T)
  } else{
     cat("No sample key or comparisons found. Can not run differential 
          expression analysis.\n")
  }
} else {
    cat("Differential expression analysis is turned OFF.\n")
}

