bic.setup.directories <- function(counts.dir, clustering.dir, all.gene.dir = NULL, 
                                  diff.exp.dir = NULL, gsa.dir = NULL, 
                                  multi.comps = FALSE, comp.set.nums = NULL){
    
    outDirs <- list()
    all.gene <- all.gene.dir 
    de       <- diff.exp.dir
    gsa      <- gsa.dir
    counts   <- counts.dir
    clusters <- clustering.dir
    
    outDirs[[1]] <- list(countsDir  = counts.dir,
                         allGeneDir = all.gene.dir, 
                         DEdir      = diff.exp.dir,
                         GSAdir     = gsa.dir, 
                         clusterDir = clustering.dir)

    if(!is.null(diff.exp.dir)){
        outDirs[[1]]$DEfigDir   <- file.path(outDirs[[1]]$DEdir, "figures")
    }

    if(multi.comps){
        for(i in comp.set.nums){
            sfx <- paste0("_", i)
            setName <- paste0("comparisons_", i)
            dirs <- list(countsDir = file.path(dirname(counts), setName, paste0(basename(counts), sfx)),
                         clusterDir = file.path(dirname(clusters), setName, paste0(basename(clusters), sfx)))

            if(!is.null(diff.exp.dir)){
                dirs$allGeneDir <- file.path(dirname(all.gene), setName, paste0(basename(all.gene), sfx))
                dirs$DEdir      <- file.path(dirname(de), setName, paste0(basename(de), sfx))
                dirs$DEfigDir   <- file.path(dirs$DEdir, "figures")
            }
            if(!is.null(gsa.dir)){
                dirs$GSAdir <- file.path(dirname(gsa), setName, paste0(basename(gsa), sfx))
            }
            outDirs[[i]] <- dirs
        }
    }
    
    for(dr in unlist(outDirs)){
        if(!file.exists(dr)){
            cat(paste0("Creating directory: ", dr, "\n"))
            dir.create(dr, recursive = T, mode = "0755")
        } else {
            cat(paste0("Directory already exists: ", dr, "\n"))
        }
    }

    outDirs
} 

bic.empty.deseq.table <- function(condA, condB){
    res <- as.data.frame(matrix(NA, nrow=1, ncol=6))
    names(res) = c("GeneID","GeneSymbol", "P.adj",
                   paste("log2[",condB,"/",condA,"]",sep=""),
                   paste("Mean_at_cond_",condA, sep=""),
                   paste("Mean_at_cond_",condB, sep=""))
    res
}

bic.write.empty.results <- function(file.name, condA, condB){
    if(grepl(".txt$", file.name)){
        write_csv(bic.empty.deseq.table(condA, condB),
                  file = file.name, 
                  append = FALSE,
                  sep = "\t")
    } else if(grepl(".xlsx", file.name)){
        openxlsx::write.xlsx(bic.empty.deseq.table(condA, condB), file.name)
    } else {
        stop(paste0("Unrecognized file type: .", gsub(".*\\.", "", file.name))) 
    }
}


bic.write.all.empty.results <- function(comps, diff.exp.dir){
    for(cmp in comps){
        condA <- unlist(strsplit(cmp, " - "))[1]
        condB <- unlist(strsplit(cmp, " - "))[2]

        flNm <- file.path(diff.exp.dir, paste0("ALLResDESeq_", condB, "_vs_", condA, ".xlsx"))
        bic.write.empty.results(flNm, condA, condB)

        flNm <- file.path(diff.exp.dir,paste0("ResDESeq_", condB, "_vs_", condA, ".xlsx"))
        bic.write.empty.results(flNm, condA, condB)
    }
}

bic.get.keys.and.comparisons <- function(pre, key.file = NULL, comp.file = NULL, path = "."){

    if(!is.null(key.file) && !is.null(comp.file)){
        bic.setup.comparisons.from.multi.column.files(key.file, comp.file)
    } else {
        bic.setup.comparisons.from.individual.files(pre, path = path)
    }
}

bic.qc.comparison.set <- function(key, comps){

    ec <- 0
    conds <- unlist(strsplit(comps, " - "))
    keyNotComps <- setdiff(key$Group, conds)
    compsNotKey <- setdiff(conds, key$Group)

    if(any(duplicated(key$Sample))){
        log_error(paste0("Duplicate samples found in key: ", 
                         paste0(unique(key$Sample[duplicated(key$Sample)]), collapse = ", ")))
        ec <- 1
    }   
    if(length(compsNotKey) > 0){
        log_error(paste0("The following sample groups found in comparisons file but not sample key: ",
                         paste(compsNotKey, collapse = ", ")))    
        ec <- 1
    }
    groupCounts <- key %>% group_by(Group) %>% summarize(Count = n())
    if(any(groupCounts$Count < 2)){
        log_warn(paste0("The following groups apply to only one sample (no replicates): ",
                        paste0(groupCounts %>% filter(Count < 2) %>% pull(Group), collapse = ", ")))
    }
    if(length(keyNotComps) > 0){    
        log_warn(paste0("The following sample groups found in key but not comparisons: ", 
                         paste(keyNotComps, collapse = ", ")))
    }
    if(ec == 1){
        stop("Errors found in key/comparison files. Please revise and rerun.")
    }
}

bic.setup.comparisons.from.standard.files <- function(keyFile, compFile){

    key   <- read.csv(keyFile, header = F, sep = "\t", check.names = F) %>%
             as_tibble() %>%
             select_if(~!all(is.na(.))) %>%
             rename(Sample = V1, Group = V2) 
    if(any(grepl("EXCLUDE", key$Group))){
        log_warn(paste0("Excluding samples: ", 
                        paste0(key %>% filter(grepl("EXCLUDE", Group)) %>% pull(Sample), collapse = ", ")))
        key <- key %>%
               filter(!grepl("EXCLUDE", Group))
    }
    comps <- read.csv(compFile, header = F, sep = "\t", check.names = F)
    comps <- paste(comps[,1], comps[,2], sep = " - ")

    bic.qc.comparison.set(key, comps)

    list(keys = key, comparisons = comps)

}


bic.setup.comparisons.from.multi.column.files <- function(keyFile, compFile){
    allKeys <- list()
    allComps <- list()
    
    keys <- read.csv(keyFile, header = T, sep = "\t", check.names = F)   ## Sample Comp_1 Comp_2
    comps <- read.csv(compFile, header = T, sep = "\t", check.names = F) ## CompSet  Group1   Group2
                                                        ## Comp_1   GroupA   GroupB
                                                        ## Comp_1   GroupD   GroupB 
                                                        ## Comp_2   GroupX   GroupY
    for(i in 2:ncol(keys)){
        key <- keys[,c(1,i)]
        if(any(grepl("EXCLUDE", key[,2]))){
            log_warn(paste0("Excluding samples: ", 
                             paste(key[grep("EXCLUDE", key[,2]), 1], collapse = ", ")))
            key <- key[!grepl("EXCLUDE", key[,2]),]
        }

        key <- key[!is.na(key$Sample) & key$Sample != "",]
        allKeys[[i]] <- key

        cmps <- comps[comps$Comparison == paste0("Comp_", i), c(2,3)]
        allComps[[i]] <- paste(cmps[,1], cmps[,2], sep = " - ")

        bic.qc.comparison.set(allkeys[[i]], allComps[[i]])
    }

    return(list(keys = allKeys, comparisons = allComps))
}


bic.setup.comparisons.from.individual.files <- function(pre, path = "."){

    keys  <- file.path(path, dir(path)[grepl("sample_key", dir(path))])
    comps <- file.path(path, dir(path)[grepl("comparisons", dir(path))])

    keyNums  <- gsub(".txt", "", gsub(".*_sample_key", "", keys))
    compNums <- gsub(".txt", "", gsub(".*_sample_comparisons", "", comps))

    if(!identical(keyNums, compNums)){
        stop("Sample key and comparison files must have 1:1 relationship indicated by a single number appended to each file name") 
    } 

    if(length(keys) == 1 && keyNums == ""){
        keyNcomp <- bic.setup.comparisons.from.standard.files(keys, comps)
        return(list(keys = list(keyNcomp$keys), comparisons = list(keyNcomp$comparisons)))
    } 

    allKeys <- list()
    allComps <- list()

    for(num in sort(as.numeric(keyNums))){
        log_debug(paste0("Setting up comparison set [", num, "]"))
        keyFile  <- keys[grep(paste0(pre, "_sample_key", num, ".txt"), keys)]
        compFile <- comps[grep(paste0(pre, "_sample_comparisons", num, ".txt"), comps)]
        
        log_debug(paste0("  key file: ", keyFile))
        log_debug(paste0("  comparisons file: ", compFile))        

        keyNcomp <- bic.setup.comparisons.from.standard.files(keyFile, compFile)

        allKeys[[num]] <- keyNcomp$keys
        allComps[[num]] <- keyNcomp$comparisons 
    }

    return(list(keys = allKeys, comparisons = allComps))

}

#' Convert a character vector to a delimited string
#' 
#' Given a character vector and a delimiter, join all
#' items in the vector, separated by the given delimiter
#' 
#' @param x          character vector
#' @param delim      delimiter to join vector elements
#' 
#' @return a string
#'
#' @export
bic.join.strings<-function(x, delim){
  tmp=""
  for(i in 1:length(x)){
    tmp=paste(tmp,x[i], delim, sep="")
  }
  return(tmp)
}


#' Average all values in a matrix that have the same 'name'
#'
#' Given a two-column matrix where the first column may contain duplicate 
#' values ('names'), for every one of those duplicate values find the average 
#' of the corresponding numbers in the second column
#' 
#' @param dat  a matrix where the first column may contain duplicate values
#'             (strings or numbers) and the second through Nth columns contain 
#'             only numbers
#' @return a matrix containing averages of values whose names were the
#'         the same in the original matrix, plus the original values in the 
#'         second through Nth columns of the input matrix whose 'names' were 
#'         unique
#'
#' @export
bic.average.by.name <- function(dat){
  ugn <- unique(dat[,1])
  dat.mean <- matrix(NA,ncol=ncol(dat)-1,nrow=length(ugn))
  for(nn in 1:length(ugn)){
    kk <- which(dat[,1]==ugn[nn])
    if(length(kk)>=2){
      val<-apply(matrix(as.numeric(dat[kk,-1]),nrow=length(kk)), 2, mean)
    }
    if(length(kk)==1){
      val=as.numeric(dat[kk,-1])
    }
    dat.mean[nn] <- val       
  }
  rownames(dat.mean) <- ugn
  return(dat.mean)
}


#' Read data file and return as a matrix 
#' 
#' @param dat    file to read
#' @param header logical to indicate whether file contains a header. Default: FALSE
#' @param sep    delimiter in file (Default: "")
#' @param skip   number of lines to skip before reading file (Default: 0)
#'
#' @export
bic.file2matrix <- function(dat,header = FALSE, sep = "", skip = 0){
  return(as.matrix(read.delim(dat,header=header,sep=sep,skip=skip)))
}

#' Search for a non-empty file by pattern
#'
#' In the specified directory, search files for given pattern
#' and if found, determine if file is not empty
#' 
#' @param path     full path to directory in which to search for file
#' @param pattern  pattern to search
#'
#' @return NULL if file is not found or is empty, otherwise return 
#'         FIRST file name that matches pattern and is not empty. NOTE: if 
#'         multiple files match pattern, ONLY the FIRST file is returned,
#'         without warning
#'
#' @export
bic.non.empty.file.found <- function(path,pattern){
  fileName = dir(path)[grep(pattern,dir(path))]
  if(length(fileName)==0 || length(fileName)==0){
    return(NULL)
  } 
  fl = paste(path,fileName[1],sep="/")
  if(!file.exists(fl) || file.info(fl)$size == 0){
    return(NULL)
  }
  return(fl)
}

#' Write file containing tab delimited data 
#' 
#' Wrapper around \code{write.table()} with several defaults set.
#' 
#' @param dat           matrix containing data to be writtedn
#' @param file.name     file name
#' @param add.rownames  include row names in file; Default: FALSE
#' @param col.names     include column names in file; Default: TRUE
#' @param quote         include quotes around each matrix cell; Default: FALSE
#' @param sep           delimiter; Default: "\\t"
#'
#' @export
bic.write.dat<-function(dat, file.name, add.rownames=F,col.names=T,quote=F,sep="\t")
{
  if(grepl("\\.xlsx$",file.name)){
     openxlsx::write.xlsx(as.data.frame(dat), file=file.name)
  } else {
     write.table(dat, file=file.name, col.names=col.names, row.names=add.rownames, quote=quote, sep=sep)
  }
}

#' Convert all values in matrix to numeric
#'
#' Quick way to ensure you are working with a numeric matrix 
#' 
#' @param dat  matrix containing data to be converted
#' @return numeric matrix
#'
#' @export
bic.matrix2numeric<-function(dat){
  res=matrix(as.numeric(dat),nrow=nrow(dat))
  rownames(res)=rownames(dat)
  colnames(res)=colnames(dat)
  return(res)
}

#' Remove a column from a matrix and use its values as the matrix's
#' rownames.
#'
#' Shortcut to assigning row names to a matrix of data read in from 
#' a file. By default, the first column of data is used.
#'
#' @param dat    matrix containing extra column to be used as row names
#' @param nn     column to be removed and used as row names
#' @return matrix with new rownames and one less column
#'
#' @export
bic.make.rownames<-function(dat, nn=1)
{       
  rownames(dat)=dat[,nn]
  return(dat[,-nn])
}

#' Add a column to matrix containing its row names
#' 
#' Create a new column containing the row names of the matrix
#' and add it as the first column in the matrix. Takes in optional
#' second argument to set the name of the new column.
#' 
#' @param dat  matrix
#' @param h    optional header for new column
#' @return matrix containing new extra column made of existing matrix's row names
#'
#' @export
bic.add.rownames2dat<-function(dat,h=""){       
  tmp=as.matrix(cbind(rownames(dat),dat))
  colnames(tmp)[1]=h
  return(tmp)
}

#' Remove a row from matrix to use the values as column names 
#'
#' Quick way to assign column headers to a matrix of data read
#' in from file. By default, the first row is used.
#' 
#' @param dat  matrix containing extra row to be used as column names
#' @param nn   optional row number to use as column names
#' @return data matrix with new column names and one less row
#'
#' @export
bic.make.colnames<-function(dat,nn=1){       
  colnames(dat)=dat[nn,]
  return(dat[-nn,])
}

#' Remove a substring from a string
#' 
#' For every string in vector x, split the string
#' on substring \code{substr} and return either
#' the string preceding trailing \code{substr}
#'
#' @param x         character vector containing strings to 
#'                  be split
#' @param substr    substring to be used as delimiter for
#'                  splitting
#' @param part2take if 1, return string preceding \code{substr}
#'                  if 2, return string trailing \code{substr}
#'
#' @return either the string preceding or trailing \code{substr}
#' @export
bic.remove.sub<-function(x, substr="hsa", part2take=1){
  tmp=NULL
  for(n in 1:length(x)){
    if(length(grep(substr,x[n]))>0){
      tmp = c(tmp, unlist(strsplit(x[n],substr))[[part2take]])
    }
    if(length(grep(substr,x[n]))==0){
      tmp = c(tmp, x[n])
    }
  }
  return(tmp)
}

#' Get info for one or multiple packages installed in current R version
#' 
#' Given a vector of package names and fields to look up, return a 
#' data frame with the requested information
#' 
#' @param pkgs    a vector of package names to look up
#' @param fields  a vector of information fields to include in results; Default:
#'                all available fields
#'
#' @return A data frame of all packages and fields for which information
#'         could be found, or NULL if no installed packages match
#'         package name(s)
#'
#' @export
bic.get.package.info <- function(pkgs,fields=NULL){

  if(is.null(fields)){
    ## defaults
    fields = c("Package","LibPath","Version",
               "Priority","Depends","Imports",
               "LinkingTo","Suggests","Enhances",
               "License","License_is_FOSS","License_restricts_use",
               "OS_type","MD5sum","NeedsCompilation",
               "Built")
  }
  ip <- as.data.frame(installed.packages())
  psearch <- paste(pkgs,collapse="|")
  fsearch <- paste(fields,collapse="|")

  res <- as.data.frame(ip[grep(psearch,rownames(ip),ignore.case=TRUE),grep(fsearch,names(ip),ignore.case=TRUE)])
  names(res) <- names(ip)[grep(fsearch,names(ip),ignore.case=TRUE)]
  return(res)

}
