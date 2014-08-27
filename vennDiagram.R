load("All_Results.RData")

## results of ALL comparisons that have been run are stored
## in a list called "AllRes"

## first just print to screen the order of comparisons to find the
## indexes of the comparisons to be included in the venn diagram
## COPY AND PASTE this for loop only:
for (x in 1:length(AllRes)){
   cat(c(x,": ",colnames(AllRes[[x]]$ans)[grep("log2",colnames(AllRes[[x]]$ans))],"\n"))
}

###################################
#######  EDIT THIS SECTION  #######

## assign results of each individual comparison
## to its own variable (with a name indicative of what the comparison
## is)

## use the indexes of each comparison of interest you just found
## and assign the results of those comparisons to meaningful
## variable names (i.e., change these variable names and indexes
## according to what was printed by the loop above)
gfpVSpd1 = AllRes[[1]]$ans
pd1VSdn = AllRes[[2]]$ans
gfpVSdn = AllRes[[3]]$ans

## name your comparison of comparisons
comp = "GFP_plus_vs_PD1_plus___vs___GFP_plus_vs_double_neg___vs___PD1_plus_vs_double_neg"

## assign human-readable names to each comparison
comp.names=c("GFP_plus_vs_PD1_plus","GFP_plus_vs_double_neg","PD1_plus_vs_double_neg")

## create a list of results of each of the pairwise comparisons
## IMPORTANT: MUST BE IN THE SAME ORDER AS COMP.NAMES ABOVE
all.dat = list()
all.dat[[1]] = gfpVSpd1
all.dat[[2]] = gfpVSdn
all.dat[[3]] = pd1VSdn

###### END SECTION TO BE EDITED #####
#####################################


#######################################
### COPY AND PASTE EVERYTHING BELOW ###


## create output directory
out.dir = comp
dir.create(out.dir)

## name the pdf file
venn.file = paste("venn",comp,".pdf",sep="")

## get a unique list of genes that are UP/DOWN in 
## one or more of the three comparisons
all.ids = c()

for (x in 1:length(all.dat)){
    all.ids = c(all.ids,as.vector(all.dat[[x]][,"ID"]))
}
all.ids = unique(all.ids)

## create a matrix where columns are comparisons
## and rows are geneIDs that were found to be UP/DOWN in 
## one or more of the two comparisons. Assign a 0 if a gene is not up
## or down in a comparison, 1 if it is up, -1 if it is down
res = matrix(0,nrow=length(all.ids),ncol=length(all.dat))
rownames(res)=all.ids
colnames(res)=comp.names

for (x in 1:length(all.dat)){
    dat = all.dat[[x]]
    res[as.vector(dat[which(dat[,grep("log2",colnames(dat))]>0),"ID"]),colnames(res)[x]] = 1
    res[as.vector(dat[which(dat[,grep("log2",colnames(dat))]<0),"ID"]),colnames(res)[x]] = -1
}

## generate venn diagram based on matrix
library(limma)
pdf(paste(out.dir,"/",venn.file,sep=""),width=25,height=15)
vennDiagram(res,include=c("up","down"),counts.col=c("red","blue"),lwd=2,cex=1.5)
dev.off()

make.multi.comp.matrix <- function(ids,all.dat,compsToExclude=NULL){

    u = 1
    if(!is.null(compsToExclude) && compsToExclude==1){
        u = 2
    } 

    gns = as.vector(all.dat[[u]][ids,"GeneSymbol"])
    final.dat = matrix(NA,nrow=length(ids),ncol=2*length(all.dat)+2)
    rownames(final.dat) = ids
    colnames(final.dat)  = rep("tmp",2*length(all.dat)+2)

    for(i in 1:length(all.dat)){
        dat.to.write = as.matrix(all.dat[[i]][ids,-grep("Mean|ID|GeneSymbol",colnames(all.dat[[i]]))])
        colnames(dat.to.write)[1] = paste("P.adj",comp.names[i],sep="__")
        final.dat[,c(i*2-1,i*2)] = dat.to.write
        colnames(final.dat)[c(i*2-1,i*2)] = colnames(dat.to.write)
    }

    if(is.null(compsToExclude)){
        fc.to.avg = colnames(final.dat)[grep("log2",colnames(final.dat))]
        fdr.to.avg = colnames(final.dat)[grep("adj",colnames(final.dat))]
    } else {
        fc.to.avg = colnames(final.dat)[grep("log2",colnames(final.dat))][-compsToExclude]
        fdr.to.avg = colnames(final.dat)[grep("adj",colnames(final.dat))][-compsToExclude]
    }

    meanFC = apply(final.dat[,fc.to.avg],1,mean)
    meanFDR = apply(final.dat[,fdr.to.avg],1,mean)
    final.dat = cbind(gns,final.dat)
    final.dat[,ncol(final.dat)-1] = meanFDR
    final.dat[,ncol(final.dat)] = meanFC
    final.dat = final.dat[order(meanFDR),]
    colnames(final.dat)[c(ncol(final.dat)-1,ncol(final.dat))]=c("meanFDR_common_comps","meanFC_common_comps")
    final.dat = cbind(rownames(final.dat),final.dat)
    colnames(final.dat)[1]="ID"
    return(final.dat)
}

dirs = c("UP","DOWN")

for(dir in dirs){

    if(dir=="UP") x=1 else x=-1
    
    ### write data for genes that are UP/DOWN in BOTH comparisons

    ids_all_comparisons = rownames(res)[which(rowSums(res)==x*ncol(res))]

    final.dat = make.multi.comp.matrix(ids_all_comparisons,all.dat)
    colnames(final.dat)[c(ncol(final.dat)-1,ncol(final.dat))]=c("meanFDR","meanFC")
    file.name = paste("genes_",dir,"_in_",paste(colnames(res),collapse="_AND_"),".txt",sep="")
    write.table(final.dat,file=paste(out.dir,"/",file.name,sep=""),quote=F,row.names=F,col.names=T,sep="\t")


    ### write data for genes that are uniquely UP/DOWN in one comparison
    for (i in 1:length(colnames(res))){

        ## write lists of genes that are uniquely UP/DOWN in each condition
        uniquely_in_dir = rownames(res)[intersect(which(res[,i]==x),which(apply(as.matrix(res[,-i]),1,function(z)!x %in% z)))]
        final.dat = all.dat[[i]][uniquely_in_dir,-grep("Mean",colnames(all.dat[[i]]))]
        file.name=paste("genes_",dir,"_in_",colnames(res)[i],"_ONLY",".txt",sep="")
        write.table(final.dat,file=paste(out.dir,"/",file.name,sep=""),quote=F,row.names=F,col.names=T,sep="\t")

        if (length(colnames(res)) > 2){
            ## write lists of genes that are up/down in the other two comparisons but not the current one
            if (i==length(colnames(res))) n=1  else n=length(colnames(res)) 

            uniquely_not_in_dir = rownames(res)[intersect(which(res[,i]!=x),which(apply(as.matrix(res[,-i]),1,sum)==((length(colnames(res))-1)*x)))]
            final.dat = make.multi.comp.matrix(uniquely_not_in_dir,all.dat,compsToExclude=i)
            file.name=paste("genes_",dir,"_in_",paste(colnames(res)[-i],collapse="__AND__"),"__NOT__",colnames(res)[i],".txt",sep="")
            write.table(final.dat,file=paste(out.dir,"/",file.name,sep=""),quote=F,row.names=F,col.names=T,sep="\t")

         }
    }
}

