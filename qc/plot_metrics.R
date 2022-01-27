#! /juno/opt/common/bic/R/R-4.0.2/bin/R 

usage <- function(){
    usage.str = "\nUsage: Rscript plotMetrics.R
    \"metrics.dir    = '[Required: directory containing metrics files]'\"
    \"qc.dir         = '[Required: directory to which PDF files should be saved]'\"
    \"pre            = '[Required: file prefix for all output files]'\"
    \"Rlibs          = '[Required: path to local R libraries]'
    \"bin            = '[Required: qc pipeline root]'
    \n\n"
    cat(usage.str)
}

cat(c("\n++++++++++++++++ Plot BIC RNA-Seq QC Metrics ++++++++++++++++\n\n"))

pd <- getwd()
pre <- "" ## by default, no file prefix
col.pal <- "Spectral"
err.code <- 0
Rlibs <- NULL

## testing
#if(interactive()){
#    args <- c("metrics",
#             "metrics/images",
#             "TEST",
#             "/ifs/work/byrne/pipelines/rnaseq_pipeline/lib/R")
#} else {
    ###################################
    ## get user input
    ###################################
    args=(commandArgs(TRUE))

    if(length(args)==0){
       usage()
        q()
    }
#    for(i in 1:length(args)){
#      eval(parse(text=args[i]))
#    }
#}
metrics.dir <- args[1]
qc.dir <- args[2]
if(length(args)>=3){
  pre <- args[3]
  if(length(args)>=4){
    Rlibs <- args[4]
  }
  if(length(args)>=5){
    bin <- args[5]
  }
} 

print(paste("metrics.dir",metrics.dir))
print(paste("qc.dir",qc.dir))
print(paste("pre",pre))
print(paste("Rlibs",Rlibs))
print(paste("bin",bin))
 
if(!is.null(Rlibs)){
  .libPaths(c(Rlibs,.libPaths()))
}

#####################################
## validate input
#####################################
if (!exists("metrics.dir")){
  cat("Error: Please specify full path to metrics directory. See usage for details\n")
  usage()
  q()
}
if (!exists("qc.dir")){
  cat("Error: Please specify full path to QC directory. See usage for details\n")
  usage()
  q()
}

#####################################
### Load BIC's rnaseq library
#####################################
cat("Loading bicrnaseq library...\n")
cat("Rlibs: ",Rlibs,"\n")
source(file.path(bin,"bicrnaseqR/source_bicrnaseq.R"))
cat("Done.\n")

dir.create(qc.dir,showWarnings=FALSE)
setwd(qc.dir)

crm <- NULL
asm <- NULL
cov <- NULL
ip  <- NULL

#####################################
#### Plot data from PICARD CollectRNASeqMetrics
#####################################
pattern = '_CollectRnaSeqMetrics.txt'
cat(paste("Searching for *",pattern,"* file...",sep=""))
file <- bic.non.empty.file.found(metrics.dir,pattern)
if(!is.null(file)){
  crm <- read.delim(file)
  cat("Done.\n")

  cat("Plotting alignment distribution...")
  file.name <- paste(pre,"picard_alignment_distribution.pdf",sep="_")
  tryCatch({
      bic.plot.alignment.distribution(crm,col.pal=col.pal,file=file.name)
    }, error = function(e) {
        message(paste0("\n\tERROR: ",e,"\n"))
        err.code <- 1
    }
  )
  file.name <- paste(pre,"picard_alignment_distribution_percentage.pdf",sep="_")
  tryCatch({
      bic.plot.alignment.distribution(crm,pct=TRUE,col.pal=col.pal,file=file.name)
    }, error = function(e){
      message(paste0("\n\tERROR: ",e,"\n"))
      err.code <- 1
    }
  )  

  cat("Done.\n")

  cat("Plotting coverage bias...")
  file.name <- paste(pre,"picard_5prime3prime_bias.pdf",sep="_")
  tryCatch({
      bic.plot.5prime3prime.bias(crm,col.pal=col.pal,file=file.name)
    }, error = function(e){
        message(paste0("\n\tERROR: ",e,"\n"))
        err.code <- 1
    }
  )
  cat("Done.\n")
} else {
  cat(paste("WARNING: No file matching pattern '",pattern,"' found.\n",sep=""))
}

#####################################
#### Plot data from PICARD AlignmentSummaryMetrics
#####################################
pattern <- '_AlignmentSummaryMetrics.txt'
cat(paste("Searching for *",pattern,"* file...",sep=""))
file <- bic.non.empty.file.found(metrics.dir,pattern)
if(!is.null(file)){
  asm <- read.delim(file)
  cat("Done.\n")

  cat("Plotting alignment summary...")
  file.name <- paste(pre,"picard_alignment_summary.pdf",sep="_")
  tryCatch({
      bic.plot.alignment.summary(asm,col.pal=col.pal,file=file.name)
    }, error = function(e){
        message(paste0("\n\tERROR: ",e,"\n"))
        err.code <- 1
    }
  )
  file.name <- paste(pre,"picard_alignment_summary_percentage.pdf",sep="_")
  tryCatch({
      bic.plot.alignment.summary(asm,pct=TRUE,col.pal=col.pal,file=file.name)
    }, error = function(e){
        message(paste0("\n\tERROR: ",e,"\n"))
        err.code <- 1
    }
  )
  cat("Done.\n")
} else {
  cat(paste("WARNING: No file matching pattern '",pattern,"' found.\n",sep=""))
}

#####################################
#### Plot data from merged CollectRNASeqMetrics Histogram files
#####################################
pattern <- '_CollectRnaSeqHistograms.txt'
cat(paste("Searching for *",pattern,"* file...",sep=""))
file <- bic.non.empty.file.found(metrics.dir,pattern)
if(!is.null(file)){
  cov <- read.delim(file)
  cat("Done.\n")

  cat("Plotting normalized coverage...")
  file.name <- paste(pre,"picard_coverage.pdf",sep="_")
  #bic.plot.coverage(cov,file=file.name)
  tryCatch({
      bic.plot.rseqc.line.chart(cov,"Normalized Coverage",file=file.name)
    }, error = function(e){
      message(paste0("\n\tERROR: ",e,"\n"))
      err.code <- 1
    }
  )
  cat("Done.\n")
} else {
  cat(paste("WARNING: No file matching pattern '",pattern,"' found.\n",sep=""))
}

#####################################
#### Plot insertion profiles for R1 and R2 
#####################################
pattern = "_insertion_profiles_Read-1.txt"
cat(paste("Searching for *",pattern,"* file...",sep=""))
file <- bic.non.empty.file.found(metrics.dir,pattern)
if(!is.null(file)){
  ip <- read.delim(file,header=T,sep="\t")
  cat("Done.\n")
  pattern = "_insertion_profiles_Read-2.txt"
  cat(paste("Searching for *",pattern,"* file...",sep=""))
file2 <- bic.non.empty.file.found(metrics.dir,pattern)
  if(!is.null(file2)){
    ip2 <- read.delim(file2,header=T,sep="\t")
    cat("Done.\n")

    cat("Plotting RSeQC Insertion profiles...")
    file.name <- file.path(qc.dir,paste(pre,"rseqc_insertion_profiles.pdf",sep="_"))
    tryCatch({
        bic.plot.rseqc.line.chart(ip,"Read 1",ip2,"Read 2",main="Insertion Profiles",file=file.name)
      }, error = function(e){
          message(paste0("\n\tERROR: ",e,"\n"))
          err.code <- 1
      }
    )
    cat("Done.\n")
  } else {
    cat(paste("WARNING: No file matching pattern '",pattern,"' found.\n",sep=""))
  }
} else {
  cat(paste("WARNING: No file matching pattern '",pattern,"' found.\n",sep=""))
}

#####################################
#### Plot deletion profile
#####################################
pattern = "_deletion_profiles.txt"
cat(paste("Searching for *",pattern,"* file...",sep=""))
file <- bic.non.empty.file.found(metrics.dir,pattern)
if(!is.null(file)){
  dat <- read.delim(file,header=T,sep="\t")
  cat("Done.\n")

  cat("Plotting RSeQC Deletion profiles...")
  file.name <- file.path(qc.dir,paste(pre,"rseqc_deletion_profiles.pdf",sep="_"))
  tryCatch({
      bic.plot.rseqc.line.chart(dat,"Deletion Profiles",file=file.name)
    }, error = function(e){
        message(paste0("\n\tERROR: ",e,"\n"))
        err.code <- 1
    }
  )  
  cat("Done.\n")
} else {
  cat(paste("WARNING: No file matching pattern '",pattern,"' found.\n",sep=""))
}


#####################################
#### Plot clipping profiles for R1 and R2 
#####################################
pattern = "_clipping_profiles_Read-1.txt"
cat(paste("Searching for *",pattern,"* file...",sep=""))
file <- bic.non.empty.file.found(metrics.dir,pattern)
if(!is.null(file)){
  dat <- read.delim(file,header=T,sep="\t")
  cat("Done.\n")
  pattern = "_clipping_profiles_Read-2.txt"
  cat(paste("Searching for ",pattern,"...",sep=""))
  file2 <- bic.non.empty.file.found(metrics.dir,pattern)
  if(!is.null(file2)){
    dat2 <- read.delim(file2,header=T,sep="\t")
    cat("Done.\n")

    cat("Plotting RSeQC Insertion profiles...")
    file.name <- file.path(qc.dir,paste(pre,"rseqc_clipping_profiles.pdf",sep="_"))
    tryCatch({
        bic.plot.rseqc.line.chart(dat,"Read 1",dat2,"Read 2",main="Clipping Profiles",file=file.name)
      }, error = function(e){
          message(paste0("\n\tERROR: ",e,"\n"))
          err.code <- 1
      }
    )  
    cat("Done.\n")
  } else {
    cat(paste("WARNING: No file matching pattern '",pattern,"'' found.\n",sep=""))
  }
} else {
    cat(paste("WARNING: No file matching pattern '",pattern,"' found.\n",sep=""))
}

#####################################
#### Plot GC content
#####################################
pattern = ".gc_content.txt"
cat(paste("Searching for *",pattern,"* file...",sep=""))
file <- bic.non.empty.file.found(metrics.dir,pattern)
if(!is.null(file)){
  dat <- read.delim(file,header=T,sep="\t")
  dat <- dat[order(dat[,1]),]
  cat("Done.\n")
    
  cat("Plotting RSeQC GC content...")
  file.name <- file.path(qc.dir,paste(pre,"rseqc_gc_content.pdf",sep="_"))
  tryCatch({
      bic.plot.rseqc.line.chart(dat,"GC content",file=file.name)
    }, error = function(e){
        message(paste0("\n\tERROR: ",e,"\n"))
        err.code <- 1
    }
  )  
  cat("Done.\n")
} else {
    cat(paste("WARNING: No file matching pattern '",pattern,"' found.\n",sep=""))
}

#####################################
#### Plot RSeQC read distribution metrics
#####################################
pattern = "_read_distribution.txt"
cat(paste("Searching for *",pattern,"* file...",sep=""))
file <- bic.non.empty.file.found(metrics.dir,pattern)
if(!is.null(file)){
  rd <- as.matrix(read.delim(file,header=T,sep="\t"))
  cat("Done.\n")

  if(nrow(rd) == 1){ 
      smps <- rd[,1]
      cols <- colnames(rd)[-1]
      rd <- as.numeric(rd[,-1])
      names(rd) <- cols 
      rd <- as.data.frame(t(as.data.frame(rd)))
      rd$Samples <- smps
  } else {
      rownames(rd) <- rd[,1]
      rd <- rd[,-1]
      rd <- apply(rd,2,as.numeric) 
      rd <- as.data.frame(rd)
      rd$Samples <- rownames(rd)
  }
  cat("Plotting RSeQC read distribution metrics...")
  file.name <- file.path(qc.dir,paste(pre,"rseqc_read_distribution.pdf",sep="_"))
  tryCatch({
      bic.plot.read.distribution(rd,file=file.name)
    }, error = function(e){
        message(paste0("\n\tERROR: ",e,"\n"))
        err.code <- 1
    }
  )  
  cat("Done.\n")

  cat("Plotting RSeQC read distribution by percentage...")
  file.name <- paste(pre,"rseqc_read_distribution_percentage.pdf",sep="_")
  tryCatch({
      bic.plot.read.distribution(rd,file=file.name,pct=TRUE)
    }, error = function(e){
        message(paste0("\n\tERROR: ",e,"\n"))
        err.code <- 1
    }
  )  
  cat("Done.\n")

} else {
    cat(paste("WARNING: No file matching pattern '",pattern,"' found.\n",sep=""))
}

cat("\n")
