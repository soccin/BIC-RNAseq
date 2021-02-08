usage <- function(){

  cat("\nUsage: Rscript qcDEinput.R

    [REQUIRED]
      --pre       project prefix
      --bin       pipeline directory

    [OPTIONAL]
      --pdir      project directory (default = $PWD)
      --Rlibs     absolute path to htseq counts file
    \n"
  )

}

tmpargs <- commandArgs(trailingOnly = TRUE)

if(any(grepl("Rlibs", tmpargs))){
    Rlibs = tmpargs[grep("Rlibs$", tmpargs) + 1]
    #print(paste0("Loading libraries from ", Rlibs))
    .libPaths(Rlibs)
    suppressPackageStartupMessages(library(R.utils))
    suppressPackageStartupMessages(library(logger))
    suppressPackageStartupMessages(library(tidyverse))
    log_threshold(DEBUG)
}

#
## after loading R.utils, get command args again in list format
#
args <- commandArgs(trailingOnly = TRUE, asValue = TRUE, adhoc = TRUE)
#
## check required args
#
if(length(args) == 0 || !all(c("pre", "bin") %in% names(args))){ 
    usage() 
    q("no", 1, FALSE)
}
for(i in 1:length(args)){
    #print(paste0(names(args)[i], " = ", args[[i]]))
    assign(names(args)[i], args[[i]])
}
if(!exists("pdir") || is.null(pdir)){ pdir <- "." }

#
## load BIC rnaseq library
#
source(file.path(bin, "bicrnaseqR/bic_util.R"))

#
## parse & qc key and comparisons files
#
keysAndComps <- bic.get.keys.and.comparisons(pre, path = pdir) 

# make sure all samples in key files are in mapping file
mapping <- file.path(pdir, dir(pdir)[grep("_sample_mapping.txt", dir(pdir))])
pSamples <- unique(read.csv(mapping, header = F, sep = "\t")$V2)

ec <<- 0
tmp <- lapply(1:length(keysAndComps$keys), function(z){ 
         x <- keysAndComps$keys[[z]]
         if(is.null(x)){ return(NULL) }
         if(any(!x$Sample %in% pSamples)){
             log_error(paste0("The following samples from key file [", z, "] are not in mapping file: ", 
                              paste(setdiff(x$Sample, pSamples), collapse = ", ")))
             ec <<- 1
         }
      }) %>%
unlist()

q(save = "no", status = ec, runLast = FALSE)

