#!/usr/bin/R

usage <- function(){

    usage.str = "\nUsage: Rscript packageInfo.R
    \"pkgs='[vector packages for which to fetch info; REQUIRED]'\"
    \"fields='[vector of information fields to fetch; DEFAULT is all fields returned by installed.packages(): 
             c(\"Package\",\"LibPath\",\"Version\",
           \"Priority\",\"Depends\",\"Imports\",
           \"LinkingTo\",\"Suggests\",\"Enhances\",
           \"License\",\"License_is_FOSS\",\"License_restricts_use\",
           \"OS_type\",\"MD5sum\",\"NeedsCompilation\",
           \"Built\")]'\"
    \n
    Example: Rscript packageInfo.R \"pkgs=c('gplots','scales')\" \"fields=c('Version','LibPath')\"

    \n\n"
    cat(usage.str)
}

## defaults
fields = c("Package","LibPath","Version",
           "Priority","Depends","Imports",
           "LinkingTo","Suggests","Enhances",
           "License","License_is_FOSS","License_restricts_use",
           "OS_type","MD5sum","NeedsCompilation",
           "Built")

## get user input
args=(commandArgs(TRUE))

if(length(args)==0){
    ## print usage
    usage()
    q()
}

for(i in 1:length(args)){
    eval(parse(text=args[i]))
}

if(!exists("pkgs")){
    usage()
    q()
}

ip <- as.data.frame(installed.packages())
psearch <- paste(pkgs,collapse="|")
fsearch <- paste(fields,collapse="|")

cat("\n")
res <- as.data.frame(ip[grep(psearch,rownames(ip),ignore.case=TRUE),grep(fsearch,names(ip),ignore.case=TRUE)])
names(res) <- names(ip)[grep(fsearch,names(ip),ignore.case=TRUE)]
res
cat("\n")
