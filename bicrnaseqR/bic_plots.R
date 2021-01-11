c3 <- c("#5395b4", "#f1592a", "#85c440")
c48 <- c("#1d915c","#5395b4",
                 "#964a48",
                 "#2e3b42",
                 "#b14e72",
                 "#402630","#f1592a",
                 "#81aa90","#f79a70", # lt pink
                 "#b5ddc2",
                 "#8fcc8b", # lt purple
                 "#9f1f63", # lt orange
                 "#865444", "#a7a9ac",
                 "#d0e088","#7c885c","#d22628","#343822","#231f20",
                 "#f5ee31","#a99fce","#54525e","#b0accc",
                 "#5e5b73","#efcd9f", "#68705d", "#f8f391", "#faf7b6", "#c4be5d", "#764c29", "#c7ac74", "#8fa7aa", "#c8e7dd", "#766a4d", "#e3a291", "#5d777a", "#299c39", "#4055a5", "#b96bac", "#d97646", "#cebb2d", "#bf1e2e", "#d89028", "#85c440", "#36c1ce", "#574a9e")


#' Arrange multiple plots on one page with one shared legend
#' 
#' Place multiple plots on one page, with a shared legend either at the right or at the bottom
#'
#' @param ...       ggplot objects
#' @param ncol      number of columns
#' @param nrow      number of rows
#' @param main      main title
#' @param position  position of legend ["bottom"|"right"]
bic.grid.arrange.shared.legend <- function(..., ncol = length(list(...)), nrow = 1,
                                           main = "", position = c("bottom", "right")) {
  ## copied and pasted from 
  ## https://github.com/tidyverse/ggplot2/wiki/share-a-legend-between-two-ggplot2-graphs

  plots <- list(...)
  position <- match.arg(position)
  top <- textGrob(main,gp=gpar(fontsize=20))
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            top = top,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           top = top,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  grid.newpage()
  grid.draw(combined)

  # return gtable invisibly
  invisible(combined)
}

#' Check formatting of input to bic.plot.rseqc.line.chart()
#'
#' Ensure object(s) passed to 'dat' and/or 'dat2' arguments of bic.plot.rseqc.line.chart()
#' are data frames where column one contains integers in order from least to greatest,
#' and each remaining column contains metric values for one sample
#' 
#' @param  dat   data frame passed to plot function
#' @param  dat2  second data frame passed to plot function (optional)
bic.check.line.chart.data <- function(dat,dat2=NULL){
  if(is.unsorted(dat[,1],strictly=TRUE)){
    stop("first data set must be sorted by first column")
  }
  tryCatch({apply(dat,2,as.numeric)},
    warning = function(w){
      stop("first data set contains non-numeric values")
    },
    error = function(e){
        stop("first data set contains non-numeric values")
    }
  )
  if(!is.null(dat2)){
    if(is.unsorted(dat2[,1],strictly=TRUE)){
      stop("second data set must be sorted by first column")
    }
    tryCatch({apply(dat2,2,as.numeric)},
      warning = function(w){
        stop("second data set contains non-numeric values")
      },
      error = function(e){
        stop("second data set contains non-numeric values")
      }
    )
    #if(!all(dim(dat) == dim(dat2))){
    #  stop(paste("first and second data sets have different dimensions\n\tdim(dat) = ",
    #             dim(dat),"\tdim(dat2) = ",dim(dat2),sep="")
    #      )
    #}
  }
  invisible(NULL)
} 

#' Plot RSeQC sequential metrics in a line graph 
#'
#' Create a line chart for merged RSeQC metrics where each line represents
#' a sample with the x axis showing sequential numeric value (e.g., read
#' read position or GC content value) and the y axis is a count
#' 
#' @param dat         data frame where column1 contains sequential numeric values 
#'                    and the remaining columns are metric values for each sample
#' @param title       plot title
#' @param dat2        second data frame with the same structure as dat1, to be plotted
#'                    adjacent to the plot of dat (optional; used to plot data for read 
#'                    1 and read 2 on one page)
#' @param title2      title of second plot (optional)
#' @param main        if plotting two datasets, this is the main title for the page
#' @param col.pal     color palette (Default: "Set3")
#' @param file        PDF file to which plot should be saved (optional)
#'
#' @export
bic.plot.rseqc.line.chart <- function(dat,title,dat2=NULL,title2=NULL,main=NULL,col.pal="Set3",file=NULL){
  width <- 10
  height <- 10
  p <- NULL
  p2 <- NULL

  ## validate input 
  bic.check.line.chart.data(dat,dat2=dat2)

  dat.m <- melt(dat,id.vars=colnames(dat)[1])
  p <- ggplot(dat.m,
              aes(x = dat.m[,1], y = as.numeric(value), color = variable)
       ) +
       geom_line(size=0.75) +
       labs(title=title) +
       xlab(colnames(dat)[1]) +
       ylab("") +
       theme(legend.position="none",aspect.ratio=4/3) +
       scale_colour_brewer(direction=-1,name="Sample",palette=col.pal) +
       scale_x_continuous()
  ## only add legend to first plot if it is the only one
  if(!is.null(dat2)){
      width <- 14
      dat2.m <- melt(dat2, id.vars=colnames(dat2)[1])
      p2 <- ggplot(dat2.m,
              aes(x = dat2.m[,1], y = as.numeric(value), color = variable)
         ) +
         geom_line(size=0.75) +
         theme(legend.position="none",aspect.ratio=4/3) +
         labs(title=title2) +
         xlab(colnames(dat)[1]) +
         ylab("") +
         scale_colour_brewer(direction=-1,name="Sample",palette=col.pal) +
         scale_x_continuous() 
  } else {
    p <- p + theme(legend.position="right")
  }
  if(!is.null(file)){ 
    pdf(file,onefile=FALSE,width=width,height=height)
  }
  if(!is.null(p2)){
    bic.grid.arrange.shared.legend(p=p,p2=p2,ncol=2,main=main,position="right")
  } else {
    print(p) 
  }
  if(!is.null(file)){
    dev.off()
  }
}

#' Check formatting of input to bic.plot.read.distribution()
#'
#' Ensure object(s) passed to 'dat' argument of bic.plot.read.distribution()
#' is a data frame containing a 'Samples' slot and at least one other slot
#' containing data to plot
#' 
#' @param  dat  data frame passed to plot function 
bic.check.read.distribution.data <- function(dat){
  if(!"Samples" %in% colnames(dat)){
    stop("data set must contain a 'Samples' column")
  }
  if(dim(dat)[2] < 2){
    stop("Data frame does not contain any data")
  }
  invisible(NULL)  
}

#' Plot RSeQC read distribution stats
#' 
#' Plot distribution of reads across different genomic features
#' like exons, introns, TSS, etc. for all samples in a project
#' 
#' @param dat        data frame containing merged output from multiple runs of
#'                   RSeQC's read_distribution.py script, where rows are
#'                   samples and columns are metrics; must contain "Samples" slot
#'                   and may contain any other slots 
#' @param pct        logical indicating that plot should show percentages
#' @param stack      logical indicating that bar chart should be stacked; Default: TRUE
#' @param horizontal logical indicating that bars should be horizontal; Default: TRUE
#' @param col.pal    name of color palette; must be from list of RColorBrewer palettes
#'                   Default: "Set3"
#' @param file       PDF file to which plot should be saved (optional)
#' @export
bic.plot.read.distribution <- function(dat,file=NULL,stack=TRUE,pct=FALSE,col.pal="Set3",horizontal=TRUE){

  ## validate input
  bic.check.read.distribution.data(dat)

  suppressMessages(dat.m <- melt(dat))

  position <- "stack"
  if(!stack){
    position <- position_dodge(width=0.7)
  }
  if(pct){
    position <- "fill"
  }

  p <- ggplot(dat.m, aes(x = Samples, y = value, fill = variable)) +
    geom_bar(stat = "identity", position = position, width = 0.7, color = "black") +
    theme(axis.text.x = element_text(angle=45,size=9,hjust=1,color="black"),
          legend.position="right",
          legend.title = element_blank()
       ) +
    labs(title="RSeQC Read Distribution") +
    scale_fill_brewer(direction=1,palette=col.pal) +
    xlab("") +
    ylab("") 
  if(horizontal){
    p <- p + coord_flip()
  }
  if(pct){
    p <- p + scale_y_continuous(labels=percent)
  }
  if(horizontal){
    p <- p + coord_flip()
  }


  if(!is.null(file)){
    pdf(file)
  }
  print(p)
  if(!is.null(file)){
    dev.off()
  }

}


#' Plot dispersion estimates
#'
#' Plot dispersion estimates for all conditions in experiment. Draw one
#' plot per PDF file.
#' 
#' @param cds          DESeq countsDataSet
#' @param out.dir      Directory in which plots should be saved; Default: $PWD
#' @param file.prefix  String to prepend to PDF file names (optional)
#'
#' @export
bic.plot.dispersion.estimates <- function(cds,out.dir=getwd(),file.prefix=""){
  if(file.prefix != ""){
    file.prefix <- paste(file.prefix,"_",sep="")
  }
  for(cond in ls(cds@fitInfo)){
    pdf(file.path(out.dir,paste(file.prefix,"dispersion_estimates_",cond,".pdf",sep="")))
    plotDispEsts(cds,name=cond)
    dev.off()
  }  
}

#' Draw a heatmap of counts table
#' 
#' Draw heatmap of either raw counts or variance stabilization
#' transformed (?) data for all samples. Write plot to PDF file.
#' 
#' @param cds       DESeq countsDataSet object
#' @param file      PDF file to which heatmap should be saved (optional)
#' @param num.gns   Include this number of most highly expressed genes in 
#'                  the heatmap; Default: 100
#' @param transform logical indicating whether to transform counts; 
#'                  Default: FALSE
#' @export
bic.deseq.heatmap <- function(cds,file=NULL,num.gns=100,transform=FALSE){
  hmcol <- colorRampPalette(brewer.pal(9,"GnBu"))(100)
  dat <- NULL
  if(transform){
    cdsBlind <- estimateDispersions(cds, method="blind")
    vst <- varianceStabilizingTransformation(cdsBlind)
    dat <- exprs(vst)
    dat <- dat[!grepl("alignment_not_unique|ambiguous|no_feature|not_aligned|too_low_aQual",
                           rownames(dat)),]
    select <- order(rowMeans(dat),decreasing=TRUE)[1:num.gns]
    dat <- dat[select,]
  } else {
    counts.cds <- counts(cds)
    counts.cds <- counts.cds[!grepl("alignment_not_unique|ambiguous|no_feature|not_aligned|too_low_aQual",
                           rownames(counts.cds)),]
    select <- order(rowMeans(counts.cds), decreasing=TRUE)[1:num.gns]
    dat <- counts.cds[select,]
  }
  if(!is.null(file)){
    pdf(file)
  }
  par(cex.main=0.8)
  heatmap.2(dat, col=hmcol, trace="none", margin=c(10,10),
            main=paste("Expression of the \n",num.gns," most highly expressed genes",sep=""),
            cexRow=0.6,cexCol=1.0,keysize=1.5,key.title=NA)
  if(!is.null(file)){
    dev.off()
  }
}

#' Draw a heatmap showing sample-to-sample distances using variance 
#' stabilized data
#' 
#' Find potential sample mislabeling errors by viewing distances
#' between every pair of samples. Save plot to PDF file.
#' 
#' @param cds    DESeq countDataSet
#' @param conds  vector of sample conditions that was used to generate \code{cds}
#' @param file   PDF file to which heatmap whould be saved (optional)
#' 
#' @export
bic.sample.to.sample.distances <- function(cds,conds,file=NULL){
  #if(is.null(file)){
  #  file = "heatmap_sample_to_sample_distances.pdf"
  #}
  hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
  cdsBlind <- estimateDispersions(cds,method="blind")
  vst <- varianceStabilizingTransformation(cdsBlind)
  distances <- dist(t(exprs(vst)))
  dist.mat <- as.matrix(distances)
  rownames(dist.mat) = colnames(dist.mat) = with(pData(cdsBlind), paste(rownames(dist.mat), " : " ,conds,sep=""))
  if(!is.null(file)){
    pdf(file)
  }
  par(cex.main=0.8)
  heatmap.2(dist.mat, main="Sample to Sample Distances", trace="none", col=rev(hmcol), margin=c(13,13), key.title=NA)
  if(!is.null(file)){
    dev.off()
  }
}

#' Plot PCA
#' 
#' Run DESeq \code{plotPCA()} on variance stabilised tranformed data
#' 
#' @param cds  DESeq countDataSet
#' @param file PDF file to which plot should be saved
#' @export
bic.deseq.plot.pca <- function(cds,file=NULL){
  #if(is.null(file)){
  #  file="pca.pdf"
  #}
  cdsBlind <- estimateDispersions(cds,method="blind")
  vst <- varianceStabilizingTransformation(cdsBlind)
  if(!is.null(file)){
    pdf(file)
  }
  plt <- DESeq::plotPCA(vst)
  print(plt)
  if(!is.null(file)){
    dev.off()
  }
}

#' Plot MA
#'
#' Plot log2 fold changes versus mean of normalized counts for a DESeq
#' comparison
#' 
#' @param res   object returned from DESeq function nbinomTest
#' @param file  PDF file to which plot should be saved (optional)
#' @export
bic.plot.ma <- function(res, file=NULL){
  #if(is.null(file)){
  #  file = "ma_plot.pdf"
  #}
  if(!is.null(file)){
    pdf(file)
  }
  DESeq::plotMA(res)
  if(!is.null(file)){
    dev.off()
  }
}

#' Histogram of p-values
#' 
#' @param dat  data frame with at least one column with name 'pvals'
#' @param file PDF file to which plot should be saved (optional)
#' @export
bic.pval.histogram <- function(dat,file=NULL){
  #if(is.null(file)){
  #  file="pval_histogram.pdf"
  #}
  if(!is.null(file)){
    pdf(file)
  }
  hist(dat$pval, breaks=100, col="blue", border="slateblue", main="",xlab="p-value",ylab="Frequency")
  if(!is.null(file)){
    dev.off()
  }
}

#' Validate data from PICARD CollectRnaSeqMetrics or AlignmentSummaryMetrics
#'
#' Check that data frame contains the correct column names and that
#' all data is numeric
#'
#' @param dat    data frame passed to plot function
#' @param name   name of plot (Options: ["alignment.distribution" |
#'                                       "5prime3prime.bias" |
#'                                       "alignment.summary"]
bic.check.picard.data <- function(dat,name,max.rows=50){

  required.slots <- switch(name,
                           "alignment.distribution" = c("SAMPLE","RIBOSOMAL_BASES","CODING_BASES",
                                                      "UTR_BASES","INTRONIC_BASES","INTERGENIC_BASES"),
                           "5prime3prime.bias" = c("SAMPLE","MEDIAN_CV_COVERAGE","MEDIAN_5PRIME_BIAS",
                                                 "MEDIAN_3PRIME_BIAS","MEDIAN_5PRIME_TO_3PRIME_BIAS"),
                           "alignment.summary" = c("SAMPLE","CATEGORY","PF_READS","PF_READS_ALIGNED")
                          )
  missing.slots = c()
  for(rs in required.slots){
    if (!rs %in% colnames(dat)){
      missing.slots = c(missing.slots,rs)
    }
  }
  if(length(missing.slots) > 0){
    stop(paste("data is missing the following required column(s): ",
               paste(missing.slots,collapse=", "),
               sep="")
        )
  }
  tryCatch({apply(dat[,required.slots[-grep("SAMPLE|CATEGORY",required.slots)]],1,as.numeric)},
    warning = function(w){
      stop("data set contains non-numeric values")
    },
    error = function(e){
      stop("data set contains non-numeric values")
    }
  )
  if(name=="alignment.summary"){
    max.rows <- max.rows * 3
    if(!"FIRST_OF_PAIR" %in% dat$CATEGORY | !"SECOND_OF_PAIR" %in% dat$CATEGORY){
      stop("'CATEGORY' column must contain at least one instance of 'FIRST_OF_PAIR' and one of 'SECOND_OF_PAIR'")
    }
  }

  if(nrow(dat) > max.rows){
      warning("Too many samples to plot. Skipping.")
      return(1)
  }

  invisible(NULL)
}


#' Plot distibution of alignment across different genomic locations
#'
#' Bar chart showing for each sample the distribution of read
#' alignments across ribosomal, coding, UTR, intronic and 
#' intergenic bases
#' 
#' @param dat        data frame consisting of data from CollectRNASeqMetrics
#'                   output
#' @param col.pal    name of color palette; must be from list of RColorBrewer palettes
#'                   Default: "Set3"
#' @param file       PDF file to which plot should be saved (optional)
#' @param horizontal logical indicating that bars should be horizontal; Default: TRUE
#' @param pct        plot percentages
#'
#' @export
bic.plot.alignment.distribution <- function(dat,pct=FALSE,horizontal=TRUE,col.pal="Set3",file=NULL,max.rows=50){
  ## validate data
  chk <- bic.check.picard.data(dat,"alignment.distribution")
  if(!is.null(chk)){ return(NULL) }
  y <- data.frame(Sample = dat$SAMPLE, 
                 Ribosomal = dat$RIBOSOMAL_BASES, 
                 Coding = dat$CODING_BASES,
                 UTR = dat$UTR_BASES,
                 Intronic = dat$INTRONIC_BASES,
                 Intergenic = dat$INTERGENIC_BASES)

  if(nrow(y) > max.rows){
      warning("Too many samples to plot. Skipping.")
      return(1)
  }

  y.m <- melt(y,id.var="Sample")

  position <- "stack"
  y.lbl <- "Total Bases (millions)"
  if(pct){ 
    position <- "fill" 
    y.lbl <- "Percent Bases"
  } 

  ### counts
  p <- ggplot(y.m, aes(x = Sample, y = value/1000000, fill = variable)) + 
    geom_bar(stat="identity", position=position, width=0.7, color="black")+
    theme(axis.text.x = element_text(angle=45,size=9,hjust=1,color="black"),
        legend.position="right",
        legend.title = element_blank()
     ) + 
    scale_fill_brewer(direction=1,palette=col.pal) +      
    labs(title="Alignment Distribution") + 
    xlab("") +  
    ylab(y.lbl)
    if(horizontal){
      p <- p + coord_flip()
    }
    if(pct){
      p <- p + scale_y_continuous(labels=percent)
    } 
  if(horizontal){
    p <- p + coord_flip()
  }
  if(!is.null(file)){
    pdf(file)
  }
  print(p)
  if(!is.null(file)){
    dev.off()
  }
}

#' Plot coverage uniformity (5' and 3' bias)
#' 
#' Bar chart showing median CV coverage, median 5' and 3' bias,
#' and median 5' to 3' bias for each sample.
#'
#' @param dat      data frame consisting of data from CollectRNASeqMetrics
#'                 output
#' @param col.pal  name of color palette; must be from list of RColorBrewer palettes
#'                 Default: "Set3"
#' @param horizontal logical indicating that bars should be horizontal; Default: TRUE
#' @param file     PDF file to which plot should be saved (optional)
#'
#' @export
bic.plot.5prime3prime.bias <- function(dat,col.pal="Set3",horizontal=TRUE,file=NULL){
  ## validate data
  dat[dat == "?"] <- NA
  chk <- bic.check.picard.data(dat,"5prime3prime.bias")
  if(!is.null(chk)){ return(NULL) }

  y <- data.frame(Sample = dat$SAMPLE,
                  cvCoverage = dat$MEDIAN_CV_COVERAGE,
                  fivePrimeBias = dat$MEDIAN_5PRIME_BIAS,
                  threePrimeBias = dat$MEDIAN_3PRIME_BIAS,
                  fivetoThreePrimeBias = dat$MEDIAN_5PRIME_TO_3PRIME_BIAS)

  suppressMessages(y.m <- melt(y))

  position <- position_dodge(width=0.7)

  p <- ggplot(y.m, aes(x = Sample, y = value, fill = variable)) + 
    geom_bar(stat = "identity", position = position, width = 0.7, color = "black") + 
    theme(axis.text.x = element_text(angle=45,size=9,hjust=1,color="black"),
          legend.position="right",
          legend.title = element_blank()
       ) +    
    scale_fill_brewer(direction=1,palette=col.pal, labels = c("Median CV of Coverage","Median 5\' Bias","Median 3\' Bias","Median 5\' to 3\' Bias")) +
    labs(title="Coverage Uniformity") + 
    xlab("") + 
    ylab("")
  
  if(horizontal){
    p <- p + coord_flip()
  }
  if(!is.null(file)){
    pdf(file)
  }
  print(p)
  if(!is.null(file)){
    dev.off()
  }

}

#' Plot normalized coverage
#'
#' Line chart showing normalized coverage for each sample
#'
#' @param dat  data frame containing combined histograms from collectrnaseqmetrics
#' @param col.pal  name of color palette; must be from list of RColorBrewer palettes
#'                 Default: "Set3"
#' @param file PDF file to which plot should be saved (optional)
#' @export
bic.plot.coverage <- function(dat,col.pal="Set3",file=NULL){

  if(nrow(dat) > 50){ warn("Too many samples. skipping."); return(NULL) }

  suppressMessages(x.m <- melt(dat, id.vars="position"))
  x.m$position <- as.integer(as.character(x.m$position))
  p <-  ggplot(x.m, aes(x = position, y = value, color = variable)) +  
        geom_line(size=0.7) + 
        theme(legend.position="right",
            legend.text = element_text(size=9),
            legend.title = element_blank()
            ) +
        labs(title="Normalized Coverage") + 
         xlab("Read position") + 
         ylab("") +
         scale_colour_brewer(direction=-1,name="Sample",palette=col.pal) 
         #scale_color_hue(colnames(dat)[-1]) + 
         #guides(colour = guide_legend(override.aes = list(size=5)))
  if(!is.null(file)){
    pdf(file)
  }
  print(p)
  if(!is.null(file)){
    dev.off()
  }
}


#' Plot numbers of mapped and unmapped reads for each sample
#'
#' Stacked bar chart showing either absolute values or percentages of mapped 
#' and unmapped reads for R1 and R2 of each sample. Takes in data frame
#' containing PICARD AlignmentSummaryMetrics, or at minimum a data frame with
#' columns: CATEGORY, PF_READS, PF_READS_ALIGNED, SAMPLE, where CATEGORY column
#' contains at least categories FIRST_OF_PAIR and SECOND_OF_PAIR.
#'
#' @param dat      data frame of PICARD AlignmentSummaryMetrics, or at least 
#'                 containing columns CATEGORY, PF_READS, PF_READS_ALIGNED and 
#'                 SAMPLE, where CATEGORY contains at least values FIRST_OF_PAIR
#'                 and SECOND_OF_PAIR
#' @param pct      logical indicating whether to show percentages rather than absolute
#'                 read counts
#' @param col.pal  name of color palette; must be from list of RColorBrewer palettes
#'                 Default: "Set3"
#' @param position position of grouped bars; Default: "stack", option: "dodge"
#' @param file     PDF file to which plot should be saved (optional)
#'
#' @export
bic.plot.alignment.summary <- function(dat,position="stack",pct=FALSE,col.pal="Set3",file=NULL){
  ## validate input data
  chk <- bic.check.picard.data(dat,"alignment.summary")
  if(!is.null(chk)){ return(NULL) }

  dat <- dat[-which(dat$CATEGORY=="PAIR"),]
  dat$UNMAPPED <- dat$PF_READS-dat$PF_READS_ALIGNED
  dat <- dat[,c("CATEGORY","PF_READS","UNMAPPED","SAMPLE")]
  dat$CATEGORY <- revalue(dat$CATEGORY,c("FIRST_OF_PAIR"="R1","SECOND_OF_PAIR"="R2"))  
  colnames(dat) <- c("Category","MappedReads","UnmappedReads","Sample")
  suppressMessages(dat.m <- melt(dat))
  
  cat.lbl.y <- -2
  y.lbl <- "Reads (xMillion)"
  if(pct){ 
    position <- "fill" 
    cat.lbl.y <- -0.04
    y.lbl <- ""
  }

  p <- ggplot(dat.m[which(dat.m$variable!="Total"),], aes(x=Category, y=value/1000000, fill=variable)) + 
   geom_bar(stat="identity", position=position, colour="black") +
   facet_wrap( ~ Sample, nrow=1, strip.position="bottom") + 
   geom_text(data=dat.m[which(dat.m$variable=="MappedReads"),], mapping=aes(x=Category, y=cat.lbl.y, label=Category), vjust=0) +
   theme(axis.text.x = element_blank(),
      axis.ticks=element_blank(),
      legend.position="right",
      legend.title = element_blank(),
      strip.text.x = element_text(size = 9,color="black",angle=90, hjust=0.5, vjust=1)
    ) +
  scale_fill_brewer(direction=1,palette=col.pal) + 
  labs(title="Alignment Summary") + 
  xlab("") +
  ylab(y.lbl)
  if(pct){
    p = p + scale_y_continuous(labels=percent) +
    ylab(y.lbl)
  } 
  if(!is.null(file)){
    pdf(file)
  }
  print(p)
  if(!is.null(file)){
    dev.off()
  }
}

#' Write PDF file containing heierarchical clustering tree
#'
#' Quickly plot heirarchical clustering of data with several default
#' parameters like width, height, font size, etc.
#' 
#' @param dat            matrix containing data to plot
#' @param file.name      file name. Default: "tmp.pdf"
#' @param title          Title of plot
#' @param sample.labels  Vector of sample labels. Default: column names in matrix
#' @param conds          vector of sample conditions, in same order as column names
#'                       in matrix or sample.labels; if given, nodes will be colored
#'                       according to conditions; Default: NULL
#' @param width          plot width (see plot docs)
#' @param height         plot height (see plot docs)
#' @param lwd            line width
#' @param cex.main       size of title
#' @param cex.lab        size of labels
#' @param cex            font size
#' @param xlab           label of x axis
#' @param ylab           label of y axis
#' @export
bic.pdf.hclust<-function(dat,conds=NULL,file.name="tmp.pdf",title="",width=20,height=14,lwd=3,
    cex.main=2.5,cex.lab=3.0,cex=3.0,xlab="",ylab="",sample.labels=NULL)
{
  if(!is.null(sample.labels)){
    colnames(dat) <- sample.labels
  }
  h <- hclust(dist(t(dat)))
  dend <- as.dendrogram(h)

  lab.cex = 2.0
  max.name = 1
  ## get longest column name in order to ensure margins are correct
  for(n in colnames(dat)){
    if (nchar(n) > max.name){
      max.name <- nchar(n)
    }
  }
  ## adjust label size if number of samples is 20
  if(ncol(dat) > 20){
    lab.cex <- 30/ncol(dat)
  }

  dend <- dendextend::set(dend,"labels_cex",lab.cex)

  ## if conditions are given, color points on leaves accordingly
  if(!is.null(conds)){
    names(conds) <- colnames(dat)
    conds <- conds[labels(dend)]
    pch <- c(19)
    col=as.numeric(as.factor(conds))
    dend <- dendextend::set(dend,"leaves_pch",pch)
    dend <- dendextend::set(dend,"leaves_col",col)
    dend <- dendextend::set(dend,"leaves_cex",lab.cex)
  }

  if(!is.null(file.name)){
    pdf(file.name,width=width,height=height)
  }
  rmar <- max.name*(lab.cex/2)
  par(oma=c(3,2,3,1))
  par(mar=c(3,2,2,rmar))
  plot(dend, horiz=T,lwd=lwd, main=title, cex.main=cex.main,xlab=xlab,ylab=ylab, cex.lab=lab.cex, cex=cex)
  if(!is.null(file.name)){
    dev.off()
  }
}


#' Plot heirachical clustering of all samples
#'
#' Plot clustering of samples using the heierarchical clustering method
#' 
#' @param norm.counts  data matrix containing normalized counts, where
#'                     a column represents a sample and a row represents
#'                     a gene. may or may not contain "GeneID" and "GeneSymbol"
#'                     columns.
#' @param conds        vector of sample conditions, in same order as column names 
#'                     in matrix or sample.labels; if given, nodes will be colored
#'                     according to conditions; Default: NULL
#' @param title        main title of plot
#' @param log2         logical indicating whether data is on log2 scale (Default: FALSE)
#' @param file.name    plot will be saved in this file; Default: 
#'                     $PWD/counts_scaled_hclust.pdf
#' @export
bic.hclust.samples <- function(norm.counts, conds = NULL, log2 = FALSE, 
                               file.name = NULL, title = ""){

  if("GeneID" %in% colnames(norm.counts) | "GeneSymbol" %in% colnames(norm.counts)){
    norm.counts <- norm.counts[,-grep("GeneID|GeneSymbol",colnames(norm.counts))]
  }
  norm.counts <- bic.matrix2numeric(as.matrix(norm.counts))

  if(length(colnames(norm.counts)) < 3){
    cat("Less than three samples; can not run cluster analysis\n")
    stop("Can not cluster less than three samples")
  }

  if(is.null(file.name)){
    file.name <- "counts_scaled_hclust.pdf"
  }
  
  if(log2==FALSE){
    pseudo <- min(norm.counts[norm.counts > 0])
    counts2hclust <- log2(norm.counts + pseudo)
  } else {
    counts2hclust <- norm.counts
  }
  tryCatch({
   bic.pdf.hclust(counts2hclust, conds = conds, file.name = file.name, title = title)
     }, error = function(err) {
        stop(err)
        traceback()
     }, warning= function(war){
        warn(war)
     }
  )

}

#' Plot MDS clustering of all samples
#'
#' Plot clustering of samples using the MDS method
#' 
#' @param norm.counts  data matrix containing normalized counts, where
#'                     a column represents a sample and a row represents
#'                     a gene
#' @param conds        vector of conditions to be appended to sample IDs
#'                     for labeling
#' @param file         PDF file to which plot should be saved (optional)
#' @param log2         logical indicating whether norm.counts is on log2 scale 
#'                     Default: FALSE
#' @param labels       logical indicating whether to include sample labels on plot; 
#'                     Default: TRUE
#' @export
bic.mds.clust.samples <- function(norm.counts, log2 = FALSE, file = NULL, 
                                  conds = NULL, labels = TRUE){

  if("GeneID" %in% colnames(norm.counts) |
     "GeneSymbol" %in% colnames(norm.counts)){
    norm.counts <- norm.counts[,-grep("GeneID|GeneSymbol",colnames(norm.counts))]
  }

  norm.counts <- bic.matrix2numeric(as.matrix(norm.counts))

  if(length(colnames(norm.counts)) < 3){
    cat("Less than three samples; can not run cluster analysis\n")
    stop("Can not cluster less than three samples")
  }

  if(is.null(conds)){
    conds=rep('s',length(colnames(norm.counts)))
  }

  if(log2==FALSE){
    pseudo <- min(norm.counts[norm.counts > 0])
    counts2hclust <- log2(norm.counts + pseudo)
  } else {
    counts2hclust <- norm.counts
  }

  md <- cmdscale(dist(t(counts2hclust)),2)

  if(!is.null(file)){
    pdf(file,width=18,height=12)
  }

  if(labels){
    plot(md, col=as.factor(conds),lwd=2.5,cex=1.5, main="Multidimensional Scaling (MDS)")
    legend("topleft", levels(as.factor(conds)),col=as.factor(levels(as.factor(conds))),pch=1,cex=1.2)
    text(md,colnames(counts2hclust),cex=0.9)
  } else {
    plot(md, col=as.factor(conds),pch=21,main="Multidimensional Scaling (MDS)",lwd=2.5,cex=1.5)
    legend("topleft", levels(as.factor(conds)),col=seq(along=levels(as.factor(conds))),pch=21,cex=1.2)
  }
  if(!is.null(file)){
    dev.off()
  }
}

#' Create PDF of generic red/green/black heatmap showing relative
#' expression values
#'
#' Given the results of DESeq comparison of two conditions and the
#' full normalized count matrix, generate a heatmap showing the relative 
#' expression values of "top 100" differentially expressed genes. If
#' there are fewer than 100 genes in the results file, all of them will
#' be included in the heatmap. Color scheme is red/black/green.
#'
#' @param norm.counts         tibble of normalized counts, where a row is a gene
#'                            and a column is a sample. May contain "GeneID" and/or
#'                            "GeneSymbol" column. If GeneSymbol column is present, 
#'                            it will be used for heatmap labeling. If not, GeneID 
#'                            column will be used.
#' @param condA               the first condition named in the bic.deseq.file name
#' @param condB               the second condition named in the bic.deseq.file name
#' @param genes               a vector of genes to be included in the heatmap;
#'                            if normalized counts matrix includes "GeneSymbol" column,
#'                            genes must be from that column. Otherwise, they must
#'                            be from the "GeneID" column.                      
#' @param file                name of PDF file to which heatmap should be written. (optional) 
#' 
#' @export
bic.standard.heatmap <- function(norm.counts, condA, condB, genes = NULL, file = NULL){

  dat <- norm.counts

  tryCatch({ 
    htmp.dat <- dat %>%
                filter(GeneID %in% genes, complete.cases(.)) %>%
                filter(row_number() < 101) %>%
                select(dplyr::matches("Gene"), grep(paste(c(condA, condB), collapse = "|"), names(.)))
                
    gns <- htmp.dat$GeneID
    if("GeneSymbol" %in% names(htmp.dat)){
        htmp.dat <- htmp.dat %>%
                    filter(!duplicated(GeneSymbol))
        gns <- htmp.dat$GeneSymbol
    }

    htmp.dat <- htmp.dat %>%
                select_if(is.numeric) %>%
                as.matrix()
    rownames(htmp.dat) <- gns

    ## replace any zeros with ones before taking log2
    htmp.dat[htmp.dat==0] <- 1
    htmp.dat <- as.matrix(log2(htmp.dat))

    if(dim(htmp.dat)[1] > 1 && dim(htmp.dat)[2] > 1){

      if(!is.null(file)){ pdf(file,width=16,height=16) }

      par(cex.main=1.4)
      heatmap.2(htmp.dat - apply(htmp.dat, 1, mean), 
                trace='none', 
                col=colorpanel(16,"green","black","red"),
                cexRow=0.9,
                dendrogram="both",
                main=paste("Top Differentially Expressed Genes ",condA," vs ",condB,sep=""), 
                symbreaks=TRUE, 
                keysize=0.5, 
                margin=c(20,10))

      if(!is.null(file)){ dev.off() }
    }
   }, error = function(e) {
      print(e)
      warning(paste0("Can not generate heatmap for ",condA," vs ",condB,". Need at least 2 rows and 2 columns to plot."))
  })
 
}


bic.deseq.qc <- function(cds, clustering.dir){
    tmp <- tryCatch({
             capture.output(
               bic.deseq.heatmap(cds,
                                 file = file.path(clustering.dir,
                                                  paste0(pre,"_heatmap_50_most_highly_expressed_genes.pdf")),
                                transform=TRUE,
                                num.gns=50))
             }, error = function(e){
                print("WARNING: There was an error while creating heatmap of the 50 most highly expressed genes.")
           })

    tmp <- tryCatch({
             capture.output(
               bic.sample.to.sample.distances(cds,
                                              conds,
                                              file = file.path(sub("/$","",clustering.dir),
                                                               paste0(pre,"_heatmap_sample_to_sample_distances.pdf"))
             ))
            }, error = function(e){
                print("WARNING: There was an error while creating heatmap of sample to sample distances")
           })

    tmp <- tryCatch({
             capture.output(
               bic.deseq.plot.pca(cds,
                                  file = file.path(clustering.dir,paste0(pre,"_PCA.pdf"))
             ))
            }, error = function(e){
               print("WARNING: There was an error while creating PCA plot.")
          })

    tmp <- tryCatch({
             capture.output(
               bic.plot.dispersion.estimates(cds,
                                             out.dir = clustering.dir,
                                             file.prefix = pre))
            }, error = function(e){
               print("WARNING: There was an error while creating dispersion estimate plot.")
          })

}
