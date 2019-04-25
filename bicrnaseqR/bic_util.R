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
  #if(grepl("\\.xlsx$",file.name)){
  #   openxlsx::write.xlsx(as.data.frame(dat), file=file.name)
  #} else {
     write.table(dat, file=file.name, col.names=col.names, row.names=add.rownames, quote=quote, sep=sep)
  #}
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
