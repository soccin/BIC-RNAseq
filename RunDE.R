#!/opt/R-2.15.0/bin/Rscript

usage <- function(){
    
    usage.str = "\nUsage: Rscript RunDE.R 
    \"bin='[Required: directory containing source R code]'\"
    \"proj.id='[Required: 4-digit project ID created by GCL]'\"
    \"output.dir='[Required: absolute path to output directory]'\"
    \"counts.file='[Required: absolute path to htseq counts file]'\"
    \"key.file='[Required: absolute path to key file]'\"
    \"comps=[Required: vector containing all comparisons to be made based on conditions in key file
             Must be in the following format: c('CondA - CondB','CondA - CondC','CondB - CondC')]\"
    
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

################
## set defaults
################
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
if(!exists("proj.id")){
    cat("Error: Please specify a Project ID number. See usage for details\n")
    q()
}
if(!exists("output.dir")){
    cat("Error: Please specify an output directory. See usage for details\n")
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
if (!exists("key.file")){
    cat("Error: Please specify a key file. See usage for details\n")
    q()
}
if (!file.exists(key.file)){
    cat(c("Error: key file",key.file,"doesn't exist.\n"))
    q() 
} 

#################
## organizing Rayas script4DEanalysis.R here
#################
#source("/home/byrne/RNAseqPipe/trunk/bin/tools.R") 
#source("/home/byrne/RNAseqPipe/trunk/bin/run_DESeq.R") 
source(paste(bin,"tools.R",sep="/"))
source(paste(bin,"run_DESeq.R",sep="/"))

sink("/dev/null")
library("DESeq")
library("limma")
sink()

HTSeq.dat=rr(counts.file,header=T)
HTSeq.dat=make.rownames(HTSeq.dat)
if ("GeneSymbol" %in% colnames(HTSeq.dat)){
    samps = colnames(HTSeq.dat)[-1]
    gns = HTSeq.dat[,1]
    idsAndGns = as.matrix(gns)
    rownames(idsAndGns) = rownames(HTSeq.dat)
    HTSeq.dat = HTSeq.dat[,-1]
} else {
    samps = colnames(HTSeq.dat)
}

counts.dat=matrix2numeric(HTSeq.dat)

### remove samples that are of no interest in this analysis ###
### reorder columns based on order of sample names in key file
#keys = as.matrix(read.delim(key.file,header=T,strip.white=T,sep="\t"))
keys = as.matrix(read.delim(key.file,header=F,strip.white=T,sep="\t"))

counts.dat = counts.dat[,keys[,1]]
conds = keys[,2]
colnames(counts.dat)=paste(keys[,1],keys[,2],sep="__")

setwd(output.dir)

### clustering
counts.log.dat=log2(counts.dat+1)
counts.log.norm.dat=normalizeBetweenArrays(counts.log.dat,method='quantile')
sink("/dev/null")
pdf.hclust(counts.log.norm.dat, file.name="hclust_quantile.pdf",title="All counts, quantile norm")
sink()

### write scaled data
dat=2^counts.log.norm.dat
if (exists("gns")){
    dat=as.matrix(cbind(rownames(dat),gns,dat))
    colnames(dat)[1]="ID"
    colnames(dat)[2]="GeneSymbol"
} else {
    dat=as.matrix(cbind(rownames(dat),dat))
    colnames(dat)[1]="GeneSymbol"
}
write.dat(dat,file.name="counts_scaled_quantile.xls")

cat("trying make.cds...\n")
cds=make.cds(counts.dat=counts.dat,conds=conds,count.cut=count.cut,libsizeQ=libsizeQ,percentile=percentile,method=method)
cat("getting scaled counts...\n")

####### this is where some features are eliminated!!!#####
counts.scaled=counts(cds,norm=T)
##########################################################

###MDS clustering###
counts2hclust=log2(counts.scaled+1)
md<- cmdscale(dist(t(counts2hclust)),2)
pdf("MDS_plot.pdf",width=18,height=12)
plot(md, col=as.factor(keys[,2]),main="",lwd=2.5,cex=1.5)
legend(-200, 70, levels(as.factor(keys[,2])),col=as.factor(levels(as.factor(keys[,2]))),pch=1,cex=1.2)
text(md,colnames(counts2hclust),cex=0.7)
dev.off()
####################

if (exists("idsAndGns")){
    idsAndGns = as.matrix(idsAndGns[rownames(counts.scaled),])
}

counts.log.dat=log2(counts.scaled+1)
cat("clustering with DESeq scaling method...\n")
sink("/dev/null")
pdf.hclust(counts.log.dat,file.name="counts_DESeqscaled_clust.pdf",title="All counts scaled using DESeq method")
sink()

cat("writing DESeq scaled counts...\n")
dat=2^counts.log.dat
if (exists("gns")){
    dat=as.matrix(cbind(rownames(dat),idsAndGns,dat))
    colnames(dat)[1]="ID"
    colnames(dat)[2]="GeneSymbol"
} else {
    gns = c()
    dat=as.matrix(cbind(rownames(dat),dat))    
    colnames(dat)[1]="GeneSymbol"
}

write.dat(dat,file.name="counts_scaled_DESeq.xls")

cat("removing low variance genes...\n")
### remove low variance genes ###
#counts.scaled2=counts.scaled.dat[jj,]
ss=apply(counts.log.dat,1,sd)
mm=apply(counts.log.dat,1,mean)
q=quantile(ss/mm)
jj=which(ss/mm>=q[2]) # genes with variance in the upper quantile ###
sink("/dev/null")
pdf.hclust(counts.log.dat[jj,],file.name="counts_scaled_DESeq_clust_low_var_removed.pdf", title="Low var counts removed; scaled using DESeq default method")
sink()

jj=which(ss/mm>=q[4]) # genes with variance in the upper quantile ###
sink("/dev/null")
pdf.hclust(counts.log.dat[jj,],file.name="counts_scaled_DESeq_clust_high_var_genes.pdf", title="Genes with var in upper quantile; scaled using DESeq default method")
sink()

colnames(counts.dat)=keys[,1]
###################################
## Differential Expression Analysis 
###################################
AllRes=vector(length=length(comps),mode="list")

for (nn in 1:length(comps)){
    
    comp = comps[nn]
    cat(c(comp,"\n"))
    condA = unlist(strsplit(comp," - "))[1]
    condB = unlist(strsplit(comp," - "))[2]

    conds = keys[grep(paste(condA,condB,sep="|"),keys[,2]),2]
    #cat(c(conds,"\n\n\n"))

    dat = counts.dat[,grep(paste(condA,condB,sep="|"),keys[,2])]
    #cat(c(colnames(dat),"\n\n\n"))

    ResDESeq=run.DESeq(counts.dat=dat, conds=conds, condA=condA, condB=condB, q.cut=q.cut, fc.cut=fc.cut, count.cut=count.cut,zeroaddQ=zeroaddQ,libsizeQ=F,method=method,gns=gns)
    AllRes[[nn]]=ResDESeq
    
}

save(AllRes,comps,file="All_Results.RData")

    
        
