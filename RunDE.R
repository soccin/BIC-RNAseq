#!/opt/R-2.15.0/bin/Rscript

usage <- function(){

    usage.str = "\nUsage: Rscript RunDE.R 
    \"bin='[Required: directory containing source R code]'\"
    \"counts.file='[Required: absolute path to htseq counts file]'\"

    \"key.file='[Required for differential expression analysis: absolute path to key file]'\"
    \"comps=[Required for differential expression analysis: vector containing all comparisons to be made based on conditions in key file
             Must be in the following format: c('CondA - CondB','CondA - CondC','CondB - CondC')]\"
    \"species='[Required for gene set analysis: (hg19|human|mm9|mm10|mouse)] Note: only human or mouse currently supported; specific 
             build does not matter, as long as it is clearly human or mouse'\"
    
    \"diff.exp=[Optional (default=TRUE): run differential expression analysis]\"
    \"GSA=[Optional (default=TRUE): run gene set analysis; if running GSA but not DESeq, BE SURE TO SET diff.exp.dir (see below) to point
              to existing DESeq results directory.]\"
    \"counts.dir=[Optional (default='$PWD/htseq')]\"
    \"clustering.dir=[Optional (default='$PWD/clustering')]\"
    \"diff.exp.dir=[Optional (default='$PWD/DESeq')]: if running DESeq, this is where output will go; if NOT running DESeq, this is
              wherever DESeq results already exist. Must be correct in order for GSA to run properly.\"
    \"gsa.dir=[Optional (default='$PWD/GSA')]\"

    \"no.replicates=[Optional (default=F): T|F, automatically sets fitType, method, and sharingMode to accommodate comparison
                     of single samples]\"
    \"q.cut=[Optional (default=0.05): insert description here]\"
    \"lfc=[Optional (default=1): insert description here]\"
    \"fc.cut=[Optional (default=2): Fold change cutoff]\"
    \"count.cut=[Optional (default=15): minimum count to be included in analysis]\"
    \"fitType=[Optional (default='parametric'): 'parametric'|'local', insert description here]\"
    \"orderPvalQ=[Optional (default=T): T|F, insert description here]\"
    \"method=[Optional (default='per-condition'): 'pooled'|'per-condition'|'blind', insert description here]\"
    \"sharingMode=[Optional (default='maximum'): 'maximum'|'fit-only'|'gene-est-only', insert description here]\"
    \"zeroaddQ=[Optional (default=T): T|F, insert description here]\"
    \"libsizeQ=[Optional (default=F): T|F, insert description here]\"
    
    \n\n"
    cat(usage.str)
}

cat(c("\n++++++++++++++++ BIC RNA-Seq Differential Expression Analysis ++++++++++++++++\n\n"))

pd = getwd()

## defaults
counts.dir="htseq"
clustering.dir="clustering"
diff.exp.dir="DESeq"
gsa.dir="GSA"
q.cut = 0.05
lfc = 1   #0.57#log2(fc.cut)
fc.cut = 2   #2^0.57
count.cut = 15
fitType = "parametric"
orderPvalQ = T
method = "per-condition"
sharingMode = "maximum"
zeroaddQ = T
libsizeQ = F
no.replicates = F
key = NULL
conds = NULL
GSA = T
diff.exp = T


################
## get user input
################
args=(commandArgs(TRUE))

if(length(args)==0){
    ## print usage
    usage()
    q()
}

for(i in 1:length(args)){
    eval(parse(text=args[i]))
}

## validate input
if(!exists("bin")){
    cat("Error: Please specify a bin directory. See usage for details\n")
    q()
}
if (!exists("counts.file")){
    cat("Error: Please specify a counts file. See usage for details\n")
    q()
}
if (!file.exists(counts.file)){
    cat(c("Error: counts file",counts.file,"doesn't exist.\n"))
    q()
}
if (exists("key.file") && !file.exists(key.file)){
    cat(c("Error: key file",key.file,"doesn't exist.\n"))
    q()
}
if (GSA && exists("key.file") && exists("comps") && !exists("species")){
    cat ("Error: Gene set analysis is turned ON. Please either turn it OFF or specify which species this data is from\n")
    q()
}
if (no.replicates){
    fitType = 'local'
    method = 'blind'
    sharingMode = 'fit-only'
}


tmp<-capture.output(suppressMessages(source(paste(bin,"tools.R",sep="/"))))
tmp<-capture.output(suppressMessages(source(paste(bin,"run_DESeq.R",sep="/"))))
tmp<-capture.output(suppressMessages(source(paste(bin,"analyze_counts.R",sep="/"))))

tmp<-capture.output(suppressMessages(library("DESeq")))
tmp<-capture.output(suppressMessages(library("limma")))

setwd(pd)
## create key if key.file is given
if (exists("key.file")){
    key = as.matrix(read.delim(key.file,header=F,strip.white=T,sep="\t"))
    key[,1] = make.names(key[,1])    
    conds=key[,2]
} 
    

############################
## normalize data (always)
############################
setwd(pd)
dir.create(counts.dir,showWarnings=FALSE)
cat("Normalizing raw counts...\n")
counts=normalize.counts(counts.file=counts.file,
                        output.dir=counts.dir,
                        conds=conds,
                        count.cut=count.cut,
                        libsizeQ=libsizeQ,
                        percentile=percentile,
                        method=method,
                        bin=bin,
                        key=key)
cat("    Done!\n\n")

############################
## cluster samples (always)
############################
setwd(pd)
dir.create(clustering.dir,showWarnings=FALSE)
if(exists("counts") && !is.null(counts) && !is.null(counts$scaled)){
    cat("Clustering all samples...\n")
    cluster.samples(counts.scaled=counts$scaled,
                    output.dir=clustering.dir,
                    conds=conds)
    cat("    Done!\n\n")
} else{
    cat("ERROR: Normalized counts do not exist. Can not run analysis.\n")
    q()
}


############################
## if comps and key, run diff.exp
############################
if(diff.exp){
    if(exists("comps") && !is.null(comps) && !is.null(key)){
        setwd(pd)
        dir.create(diff.exp.dir,showWarnings=FALSE)
        cat("Running differential expression analysis...\n")
        run.diff.exp(counts.raw=counts$raw,
                     comps=comps,
                     key=key,
                     output.dir=diff.exp.dir,
                     gns=counts$gns,
                     q.cut=q.cut,
                     lfc=lfc,
                     fc.cut=fc.cut,
                     count.cut=count.cut,
                     zeroaddQ=zeroaddQ,
                     libsizeQ=libsizeQ,
                     fitType=fitType,
                     method=method,
                     sharingMode=sharingMode)
        cat("    Done!\n\n")
    } else{
        cat("No sample key or comparisons found. Can not differential expression analysis.\n")
    }
} else {
    cat("Differential expression analysis is turned OFF.\n")
}

############################
## if GSA and species, run gene set analysis
############################
if(GSA){
    if(exists("species") && !is.null(species)){
        setwd(pd)
        dir.create(gsa.dir,showWarnings=FALSE)
        deseq.res.dir=paste(pd,diff.exp.dir,sep="/")
        cat("Running gene set analysis...\n")
        run.gene.set.analysis(species=species,
                              bin=bin,
                              deseq.res.dir=deseq.res.dir)
        cat("    Done!\n\n")
    } else{
        cat("No species specified. Can not run gene set analysis.\n")
    }
} else {
    cat("Gene set analysis is turned OFF.\n")
}

