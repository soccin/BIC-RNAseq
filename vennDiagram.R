load("All_Results.RData")

##################################
###       MANUAL CONFIG        ###

## figure out the order of comparisons in 
## AllRes (done manually for now)
Ren713vsA1A1803 = AllRes[[1]]$ans
A1A6421vsA1A1803 = AllRes[[2]]$ans
Ren713vsA1A6421 = AllRes[[3]]$ans

## name your three-way comparison
comp = "Ren713_vs_A1A1803___vs___Ren713_vs_A1A6421"
out.dir=comp
venn.file = paste("venn",comp,".pdf",sep="")

## assign human-readable names to each comparison
## IMPORTANT: MAKE SURE THESE NAMES ARE IN THE SAME ORDER
## AS THEY ARE IN all.dat (below)
comp.names=c("Ren_713_vs_A1A_1803","Ren_713_vs_A1A_6421")

## list of results of each of the three pairwise comparisons
## IMPORTANT: MUST BE IN THE SAME ORDER AS COMP.NAMES ABOVE
all.dat = list()
all.dat[[1]] = Ren713vsA1A1803
all.dat[[2]] = Ren713vsA1A6421

#################################

### COPY AND PASTE EVERYTHING BELOW ###

## create output directory
dir.create(out.dir)

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

### write lists of genes UP/DOWN in ALL three comparisons
ids_DOWN_all_comparisons = rownames(res)[which(res[,1]==-1 & res[,2]==-1)]
gns = all.dat[[1]][which(all.dat[[1]][,"ID"] %in% ids_DOWN_all_comparisons),"GeneSymbol"]
write.table(gns,file=paste(out.dir,"/","genes_DOWN_in_",paste(colnames(res),collapse="_AND_"),".txt",sep=""),quote=F,row.names=F,col.names=F)

ids_UP_all_comparisons = rownames(res)[which(res[,1]==1 & res[,2]==1)]
gns = all.dat[[1]][which(all.dat[[1]][,"ID"] %in% ids_UP_all_comparisons),"GeneSymbol"]
write.table(gns,file=paste(out.dir,"/","genes_UP_in_",paste(colnames(res),collapse="_AND_"),".txt",sep=""),quote=F,row.names=F,col.names=F)

## write lists of genes that are uniquely UP/DOWN in each condition
for (i in 1:length(colnames(res))){

    ## write lists of genes that are uniquely UP/DOWN in each condition
    uniquely_up = rownames(res)[intersect(which(res[,i]==1),which(apply(as.matrix(res[,-i]),1,function(x)!1 %in% x)))]
    gns = all.dat[[i]][which(all.dat[[i]][,"ID"] %in% uniquely_up),"GeneSymbol"]
    write.table(gns,file=paste(out.dir,"/","genes_UP_in_",colnames(res)[i],"_ONLY",".txt",sep=""),quote=F,row.names=F,col.names=F)

    uniquely_down = rownames(res)[intersect(which(res[,i]==-1),which(apply(as.matrix(res[,-i]),1,function(x)!-1 %in% x)))]
    gns = all.dat[[i]][which(all.dat[[i]][,"ID"] %in% uniquely_down),"GeneSymbol"]
    write.table(gns,file=paste(out.dir,"/","genes_DOWN_in_",colnames(res)[i],"_ONLY",".txt",sep=""),quote=F,row.names=F,col.names=F)

    if (length(colnames(res)) > 2){
        ## write lists of genes that are up/down in the other two comparisons but not the current one
        if (i==length(colnames(res))){ n=1 }
        else { n=length(colnames(res)) }

        uniquely_not_up = rownames(res)[intersect(which(res[,i]!=1),which(apply(as.matrix(res[,-i]),1,sum)==(length(colnames(res))-1)))]
        gns = all.dat[[n]][which(all.dat[[n]][,"ID"] %in% uniquely_not_up),"GeneSymbol"]
        write.table(gns,file=paste(out.dir,"/","genes_UP_in_",paste(colnames(res)[-i],collapse="__AND__"),"__NOT__",colnames(res)[i],".txt",sep=""),quote=F,row.names=F,col.names=F)

        uniquely_not_down = rownames(res)[intersect(which(res[,i]!=-1),which(apply(as.matrix(res[,-i]),1,sum)==((length(colnames(res))-1)*-1)))]
        gns = all.dat[[n]][which(all.dat[[n]][,"ID"] %in% uniquely_not_down),"GeneSymbol"]
        write.table(gns,file=paste(out.dir,"/","genes_DOWN_in_",paste(colnames(res)[-i],collapse="__AND__"),"__NOT__",colnames(res)[i],".txt",sep=""),quote=F,row.names=F,col.names=F)


}



