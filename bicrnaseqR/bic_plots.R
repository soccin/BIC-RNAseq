font <- "Roboto"
#showtext_auto()

bicTheme <- theme_minimal() +
            theme(text = element_text(size = 16, family = font, color = "black"),
                  plot.title = element_text(size = 20, margin = unit(c(0.2, 0, 0.1, 0), "in")),
                  plot.subtitle = element_text(size = 13, face = "italic", margin = unit(c(0.1, 0, 0.2, 0), "in")),
                  axis.text.x = element_text(color = "black"),
                  axis.text.y = element_text(color = "black"),
                  axis.line.x = element_line(size = 0.5),
                  axis.line.y = element_line(size = 0.5))

clrPal <- c("#e5f5e0", "royalblue4", "#fae8fc", "#c4c2c3", "indianred3", "lightsalmon", "yellowgreen", "slategray1"> "turquoise4", "lightyellow2", "steelblue3", "orange2", "violetred3", "#a1788d", "gold", "orchid4", "burlywood1", "da> seagreen3", "royalblue3", "lightpink1", "plum", "snow2", "#7e8c76")

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
#' @param file.prefix  prefix for each ouput PDF file (optional)
#'
#' @export
bic.plot.dispersion.estimates <- function(cds, out.dir = NULL, file.prefix = NULL){

  plotList <- list()

  for(cond in ls(cds@fitInfo)){
      ## data prep taken directly from plotDispEsts()
      px <- rowMeans(counts(cds, normalized = TRUE))
      sel = (px > 0)
      px = px[sel]
      py = fitInfo(cds, name = cond)$perGeneDispEsts[sel]
      ymin = 10^floor(log10(min(py[py > 0], na.rm = TRUE)) - 0.1)
      xg = 10^seq(-0.5, 5, length.out = 100)

      dat <- tibble(Gene = names(px), Vals = px, Y = py) %>% 
             mutate(Y = pmax(Y, ymin))
      dat2 <- tibble(x = xg, y = fitInfo(cds, name = cond)$dispFun(xg))
      p <- ggplot(dat, aes(x = Vals, y = Y)) +
           geom_point(color = "black", alpha = 0.7, size = 0.5) +
           scale_y_continuous(trans = "log10") +
           scale_x_continuous(trans = "log10", breaks = c(1, 100, 10000)) + 
           geom_line(data = dat2, aes(x = x, y = y), color = "red", size = 1) +
           xlab('mean of normalized counts') +
           ylab('dispersion') +
           labs(title = cond) +
           bicTheme +
           theme(plot.margin = margin(t = 0.25, b = 0.25, r = 0.25, l = 0.25, unit = "cm"))
      plotList[[cond]] <- p 
 
      if(!is.null(file.prefix)){
          if(is.null(out.dir)){ out.dir <- getwd() }
          file <- paste0(file.prefix, "_", cond, ".pdf")
          if(out.dir != ""){ file <- file.path(out.dir, file) }
          pdf(file, height = 7, width = 7) 
          print(p)
          dev.off()
      }
   }

   plotList
}


ensemblIDtoGeneSymbol <- function(ensemblIDs, species){

    library('biomaRt')
    ds <- NULL 
    if(species == 'human'){
        ds <- 'hsapiens_gene_ensembl'
    } else if (species == 'mouse') {
        ds <- 'mmusculus_gene_ensembl'
    } else {
        stop("don't know which biomart dataset to use")
    }
    mart <- useMart("ensembl", dataset=ds)
    gnList <- getBM(filters= "ensembl_gene_id", 
                    attributes= c("ensembl_gene_id","external_gene_name"), 
                    values=ensemblIDs, 
                    mart= mart)
    if(nrow(gnList) > 0 && !any(is.na(gnList[,2]))){
        return(gnList)
    } else {
        return(NULL)
    }
}

matchPanelHeights <- function(g1, g2){

    panelIDs <- unique(g1$layout[grepl("panel", g1$layout$name), "t"])
    heights <- g1$heights[panelIDs]

    g2$heights[panelIDs] <- unit.c(do.call(unit, list(heights, 'null')))

    return(list(g1,g2))
}

matchPanelWidths <- function(g1, g2){

    panelIDs <- unique(g1$layout[grepl("panel", g1$layout$name), "t"])
    widths <- g1$widths[panelIDs]

    g2$widths[panelIDs] <- unit.c(do.call(unit, list(widths, 'null')))

    return(list(g1,g2))
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
bic.high.expression.heatmap <- function(cds, num.gns=100, transform = TRUE, idMap = NULL, 
                                        height = 8.5, width = 11, file = NULL, key = NULL,
                                        clrs = NULL){

    ##
    ## PREP DATA
    ##
    if(transform){
        dat <- cds %>% estimateDispersions(method="blind") %>%
               varianceStabilizingTransformation() %>%
               exprs()
    } else {
        dat <- log(counts(cds) + 1)
    }
    dat <- dat %>%
           as_tibble(rownames = "ID") %>%
           filter(!grepl("alignment_not_unique|ambiguous|no_feature|not_aligned|too_low_aQual", "ID")) %>%
           arrange(desc(rowMeans(.[,-1]))) %>%
           filter(row_number() <= num.gns) %>%
           gather(2:ncol(.), key = "Sample", value = "Val")

    if(!is.null(idMap)){
        dat <- dat %>% left_join(idMap, by = c("ID")) %>% select(ID = GeneSymbol, Sample, Val) 
    }

    if(is.null(clrs) && !is.null(key)){
        clrs <- bic.get.colors(length(unique(key$Group)))
        names(clrs) <- unique(key$Group)
    }

    htmp <- bic.rnaseq.heatmap(dat, cols = "Sample", rows = "ID", clusterCols = TRUE, clusterRows = FALSE,
                                colAnnot = key, annClrs = clrs,
                                title = "Normalized Gene Expression", 
                                subtitle = paste(num.gns, "most highly expressed genes"),
                                width = 11, height = 8.5, file = file)
    htmp
}

#' Draw a heatmap showing sample-to-sample distances using variance 
#' stabilized data
#' 
#' Find potential sample mislabeling errors by viewing distances
#' between every pair of samples. Save plot to PDF file.
#' 
#' @param cds    DESeq countDataSet
#' @param key    tibble containing sample names and group assigments (optional) 
#' @param file   PDF file to which heatmap whould be saved (optional)
#' 
#' @export
bic.sample.to.sample.distances <- function(cds, key = NULL, clrs = NULL, file=NULL, width = 8.5, height = 11){
  
    dat <- cds %>%
           estimateDispersions(method = "blind") %>%
           varianceStabilizingTransformation() %>%
           exprs %>% t %>% dist %>% as.matrix %>% as_tibble %>%
           mutate(Sample = names(.)) %>%
           select(Sample, everything()) %>%
           gather(2:ncol(.), key = "Sample2", value = "Val")

    key2 <- NULL 
    if(!is.null(key)){
        if(is.null(clrs)){
            clrs <- bic.get.colors(length(unique(key$Group)))
            names(clrs) <- unique(key$Group)
        }
        key2 <- key %>% select(Sample2 = Sample, Group)
    }

    htmp <- bic.rnaseq.heatmap(dat, cols = "Sample", rows = "Sample2", clusterCols = TRUE, 
                                clusterRows = TRUE, 
                                rowAnnot = key2, colAnnot = key, annClrs = clrs,
                                title = "Sample to Sample Distances",
                                width = 11, height = 8.5, file = file)

    htmp
}


bic.get.colors <- function(numClrs, sort = FALSE){

    aClrs <- c("green4","hotpink","gold","mediumorchid3", "mediumturquoise", 
               "yellowgreen", "royalblue3", "violetred3", "chocolate2", 
               "plum", "deepskyblue1")

    clrsDist <- c("#7e8c76", "hotpink", "gold3", "steelblue3",
                  "gold", "burlywood1", "#a1788d", "mediumturquoise", 
                  "#fcba03", "turquoise4", "lightpink1", "royalblue3", 
                  "plum", "firebrick2", "green4", "slategray1", 
                  "chocolate2", "mediumorchid3", "royalblue4", "yellowgreen", 
                  "lightsalmon", "#c4c2c3", "deepskyblue1", "violetred3")

    clrsGradient <- c("firebrick2", "chocolate2", "#fcba03", "lightsalmon", 
                      "burlywood1", "gold", "gold3", "#7e8c76", 
                      "turquoise4", "green4", "yellowgreen", "mediumturquoise", 
                      "royalblue4", "royalblue3", "steelblue3", "deepskyblue1",
                      "slategray1", "#a1788d", "mediumorchid3", "plum", 
                      "violetred3", "hotpink", "lightpink1", "#c4c2c3") 

    lens <- list(a = length(aClrs), dist = length(clrsDist), grad = length(clrsGradient))
    idxs <- seq(numClrs)

    if(numClrs > max(unlist(lens))){
        ## return the uglier, but larger palette
        if(sort){ log_warn("Color sorting is not supported when using >", max(unlist(lens)), " colors.") }
        return(c48[idxs])
    }
    if(numClrs <= lens$a){
        if(sort){ 
            return(clrsGradient[clrsGradient %in% aClrs][idxs])
        }
        return(aClrs[idxs])
    }
    if(sort){
        return(clrsGradient[idxs])
    }
    return(clrsDist[idxs])
}

#' Plot PCA
#' 
#' Run DESeq \code{plotPCA()} on variance stabilised tranformed data
#' 
#' @param cds    DESeq countDataSet
#' @param file   PDF file to which plot should be saved
#' @param ntop   number of top genes to include in PCA analysis; default=500
#' @param labels 
#' @export
bic.deseq.plot.pca <- function(cds, file=NULL, ntop = 500, labels = FALSE, clrs = NULL, idMap = NULL,
                               height = 8.5, width = 11){
  cdsBlind <- estimateDispersions(cds,method="blind")
  vst <- varianceStabilizingTransformation(cdsBlind)

  grps <- tibble(Sample = rownames(phenoData(vst)@data), 
                 Group = phenoData(vst)@data[,2]) %>%
          mutate(Group = factor(Group, levels = unique(Group)))

  if(is.null(clrs)){
    clrs <- bic.get.colors(length(unique(grps$Group)))
    names(clrs) <- sort(unique(grps$Group))
  }

  rv = rowVars(exprs(vst))
  select = order(rv, decreasing=TRUE)[seq_len(ntop)]
  pca = prcomp(t(exprs(vst)[select,]))

  smps <- rownames(pca$x)
  plotDat <- pca$x %>%
             as_tibble() %>%
             mutate(Sample = smps) %>%
             left_join(grps, by = "Sample")

  xtitle <- paste0("PC1: ", round(summary(pca)$importance[2,"PC1"]*100, 1), "% variance")
  ytitle <- paste0("PC2: ", round(summary(pca)$importance[2,"PC2"]*100, 1), "% variance")

  marg <- unit(c(0.75, 1.5, 0.75, 1.5), "in")
  if(length(levels(plotDat$Group)) > 5){
      marg <- unit(c(0.25, 1.5, 0.25, 1.5), "in")
  }
  fontPts <- getFontSizeInPoints(height, marg, nrow(plotDat) * 0.5)
  lgndRows <- ceiling(length(levels(plotDat$Group))/5)
  lgndFS <- min(14, getFontSizeInPoints(((unit(height, "in") - (marg[1] + marg[3]))/20) * lgndRows, unit(c(0,0,0,0),"in"), lgndRows ))
  pointSz <- min(6, ptToMM(fontPts))

  plt <- ggplot(data = plotDat, aes(x = PC1, y = PC2, fill = Group, label = Sample), color = "black") +
           geom_point(size = pointSz, shape=21, alpha = 0.7) +
           #geom_point(size = 4, shape = 1, color = "black") +
           xlab(xtitle) +
           ylab(ytitle) +
           labs(title = "PCA", subtitle = "") +
           scale_fill_manual(values = clrs) +
           theme_minimal() +
           bicTheme +
           theme(text = element_text(fontPts, color = "black"),
                 axis.text.x = element_text(color = "black"),
                 axis.text.y = element_text(color = "black"),
                 axis.title.y = element_blank(),
                 axis.title.x = element_blank(),
                 plot.margin = marg,
                 legend.key.size = unit(lgndFS, "pt"),
                 legend.title = element_blank(),
                 legend.position = "bottom",
                 legend.text = element_text(size = lgndFS, hjust = 0, margin = margin(l = -0.25, unit = "in")),
                 legend.spacing.x = unit(0.25, "in"))

  if(labels){
      plt <- plt + 
             geom_text_repel(color = "black", size = 4, point.padding = 0.5)
  }

  if(!is.null(file)){ 
      pdf(file, width = 11, height = 8.5) 
      print(plt)
      dev.off() 
  }
  return(plt)
}

bic.plot.pc.loading <- function(cds, plotPCs = c("PC1", "PC2"), ntop = 500, file = NULL,
                                 idMap = NULL, extraTheme = NULL){

  cdsBlind <- estimateDispersions(cds,method="blind")
  vst <- varianceStabilizingTransformation(cdsBlind)

  grps <- tibble(Sample = rownames(phenoData(vst)@data),
                 Group = phenoData(vst)@data[,2])

  rv = rowVars(exprs(vst))
  select = order(rv, decreasing=TRUE)[seq_len(ntop)]
  pca = prcomp(t(exprs(vst)[select,]))

  margins <- list(unit(c(1, 0.25, 1.5, 1.5), unit = "in"),
                  unit(c(1, 1.5, 1.5, 0.25), unit = "in"),
                  unit(c(1, 0.25, 1.5, 0.25), unit = "in"))

  plotList <- list()
  widths <- c(1)
  arrowX <- 0.1

  breaks <- list()

  for(x in seq(plotPCs)){
    pc <- plotPCs[x]
    loadDat <- sort(pca$rotation[,pc])
    plotDat <- c(loadDat[1:10], loadDat[(length(loadDat)-9):length(loadDat)])
  
    plotDat <- tibble(ID = factor(names(plotDat), levels = names(plotDat)), 
                      Value = plotDat) %>%
               mutate(Color = ifelse(Value < 0, "blue", "red"))

    ### since arrowX is not in plotDat, when using it to build the plot, only
    ### a reference is stored in the ggplot object until it is added to a list
    ### so even if arrowX is updated, to fit data in second plot, the new value
    ### will still apply to the first plot, helping match plot widths
    arrowX <- max(arrowX, max(plotDat$Value) * 1.01)

    if(!is.null(idMap)){
      plotDat <- plotDat %>% 
                 left_join(idMap, by = "ID")
      lvls <- plotDat$GeneSymbol
      plotDat <- plotDat %>%
                 mutate(ID = factor(GeneSymbol, levels = lvls))
    }

    widths[x] <- sapply(as.vector(plotDat$ID), function(x) { nchar(x) }) %>% max 

    clrs <- c(red = "red3", blue = "blue3")

    breaks[[x]] <- seq(-round(max(abs(plotDat$Value)), 1) - 0.05, round(max(abs(plotDat$Value)), 1) + 0.05, 0.05)

    plt <- NULL
    plt <- ggplot(plotDat, aes(x = Value, y = ID, fill = Color)) +
           geom_bar(stat = "identity", color = "black") +
           geom_segment(aes(x = arrowX, y = nrow(plotDat) * 0.15, xend = arrowX, yend = nrow(plotDat) * 0.85), 
                        color = "black", linetype = "dashed") +
           geom_segment(aes(x = arrowX, xend = arrowX, y = nrow(plotDat) * 0.84, yend = nrow(plotDat) * 0.85), 
                        color = "black", lineend = "butt", linejoin = "mitre", arrow = arrow(length = unit(0.2, "in"))) +
           annotate("text", x = arrowX * 1.2, y = nrow(plotDat)/2, 
                    label = paste("Top/bottom loadings", pc), angle = 270, size = 4.5) +
           scale_fill_manual(values = clrs) +
           scale_x_continuous(breaks = breaks[[x]], position = "top") +
           geom_segment(aes(x = breaks[[x]][1], y = nrow(plotDat) + 0.7,  
                            xend = breaks[[x]][length(breaks[[x]])], yend = nrow(plotDat) + 0.7), 
                        size = 0.5, color="black") +
           bicTheme +
           theme(text = element_text(family = font, size = 16),
                 legend.position = "none",
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 axis.line.x = element_blank(),
                 axis.line.y = element_blank(),
                 axis.title.x = element_blank(),
                 axis.title.y = element_blank(),
                 axis.text.x = element_text(color = "black", size = 10),
                 axis.text.y = element_text(color = "black", size = 12),
                 plot.margin = margins[[x]]) 

    if(x == 1){
        plt <- plt + labs(title = "PC Loading")
    }

    if(!is.null(extraTheme)){
        plt <- plt + extraTheme
    }
    plotList[[pc]] <- plt
  }

  pg <- list()
  for(x in 1:length(plotList)){
      gt <- ggplot_gtable(ggplot_build(plotList[[x]]))
      gt$layout$clip[gt$layout$name == "panel"] <- "off"
      pg[[x]] <- gt
  }

  pg <- plot_grid(plotlist = pg, align = 'h', axis = 'bt', rel_widths = widths/sum(widths) + 2)  

  if(!is.null(file)){ 
      pdf(file, height = 8.5, width = 5.5 * length(plotPCs))
      grid.draw(pg)
      dev.off() 
  }

  return(pg)
}

## check if data is sorted by a given column
data.sorted.by <- function(dat, sortCol){
  groups <- seq(length(unique(dat[[sortCol]])))
  names(groups) <- unique(dat[[sortCol]])
  tmp <- dat %>% 
         mutate(GroupNum = groups[Group], row = row_number())

  for(g in names(groups)){
      tmp2 <- tmp %>% filter(!!as.name(sortCol) == g) %>% pull(row)
      if(!identical(seq(tmp2[1], tmp2[length(tmp2)]), tmp2)){
          return(FALSE)
      }
  }
  TRUE
}

bic.plot.cibersort <- function(cbs, file = NULL, annotate = FALSE, sortByGroup = FALSE, extraTheme = NULL, width = 11, height = 8.5){

  clrs <- c('B cells naive' = 'gold', 'B cells memory' = 'lightyellow2', ## yellows
            'Plasma cells' = '#c4c2c3', # gray
            'T cells CD8' = 'snow2', 'T cells CD4 naive' = '#7e8c76', 'T cells CD4 memory resting' = 'turquoise4',
            'T cells CD4 memory activated' = 'darkseagreen3', 'T cells follicular helper' = 'yellowgreen',
            'T cells regulatory (Tregs)' = '#e5f5e0', 'T cells gamma delta' = '#4D6B50', # greens
            'NK cells resting' = 'lightpink1', 'NK cells activated' = 'violetred3', # pinks
            'Monocytes' = 'royalblue4', 'Macrophages M0' = 'royalblue3', 'Macrophages M1' = 'steelblue3', 'Macrophages M2' = 'slategray1', # blues
            'Dendritic cells resting' = 'orange2', 'Dendritic cells activated' = 'lightsalmon', # oranges 
            'Mast cells resting' = 'orchid4', 'Mast cells activated' = '#a1788d', 'Eosinophils' = 'plum', 'Neutrophils' = '#fae8fc' # purples
           )

  plotDat <- cbs %>%
             gather(names(clrs), key="CellType", value="PCT") %>%
             mutate(CellType = factor(CellType, levels = unique(CellType))) %>%
             separate(Sample, c("Sample", "Group"), sep="__") 

  sortedByGroup <- data.sorted.by(plotDat %>% select(Sample, Group) %>% unique, "Group")
  if(!sortedByGroup && sortByGroup){
    log_debug("Sorting data by sample group")
    plotDat <- plotDat %>% 
               arrange(Group) %>%
               mutate(Group = factor(Group, levels = unique(Group)))
  }

  marg   <- unit(c(0.5,0.5,0.25,0.5), "in")
  pWidth <- getPlotWidthInPoints(width, marg)
  numSmps <- length(unique(plotDat$Sample))
  xtext  <- min(14, (pWidth - inchToPt(3))/(numSmps + 2))  # subtract 3in (estimated) for legend, which will not change

  plt <- ggplot(plotDat, aes(Sample, PCT, fill = CellType)) +
         geom_bar(stat = "identity", position = "stack", color = "black", width=0.85) +
         scale_fill_manual(values = clrs) +
         scale_y_continuous(labels = scales::percent, limits = c(0,1.01), expand = c(0,0)) +
         scale_x_discrete(labels = gsub("s_", "", gsub("__.*", "", plotDat$Sample))) +
         xlab("") +
         ylab("Relative Percent") +
         labs(title = "CIBERSORT") +
         bicTheme +
         theme(axis.text.x = element_text(angle = 90, size = xtext, vjust = 0.5, hjust = 1.1, color = "black"),
               axis.line.x = element_blank(),
               axis.line.y = element_blank(),
               axis.ticks.y = element_line(size = 0.5),
               axis.text.y = element_text(color = "black"),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               legend.title = element_blank(),
               legend.text = element_text(size = ifelse(numSmps > 30, 12, 14)),
               legend.key.size = unit(0.2, 'in'),
               legend.position = "right",
               plot.margin = marg) +
         guides(fill=guide_legend(ncol=1))

  if(!is.null(extraTheme)){
      plt <- plt + extraTheme
  }

  pg <- plt

  if(annotate) {
    if(!sortedByGroup && !sortByGroup){
      log_error("Data not sorted by group. Can not annotate CIBERSORT plot.")
    } else {
  
      annDat <- plotDat %>% 
                select(Sample, Group) %>% unique %>%
                mutate(Sample = gsub("s_", "", Sample), Group = gsub("_", "\n", Group), 
                       x = row_number()) %>%
                group_by(Group) %>%
                mutate(first = ifelse(x == min(x), x, NA),
                       last = ifelse(x == min(x), max(x), NA),
	     	       mid = ifelse(x == min(x), (last + first)/2, NA),
                       label = ifelse(x == min(x), Group, NA))
 
      plt2 <- ggplot(annDat, mapping = aes(x = x, y = 1), color = "black", fill = "white") +
                geom_segment(mapping = aes(x = first - 0.4, y = 1, xend = last + 0.4, yend = 1)) +
                scale_x_continuous(labels = annDat$Sample, breaks = annDat$x, 
                                   limits = c(min(annDat$x) - 0.5, max(annDat$x) + 0.5), expand = c(0,0.1)) +
                scale_y_continuous(limits = c(0, 1.05), expand = c(0,0)) +
                xlab("") + ylab("") +
                bicTheme +
                theme(axis.text.x = element_blank(), 
                      axis.line.x = element_blank(),
                      axis.line.y = element_blank(),
                      axis.ticks.y = element_blank(),
                      axis.text.y = element_blank(),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      plot.margin = margin(t = -0.75, b = 1.5, unit = "lines"))
    
      for(x in which(!is.na(annDat$label))){
          plt2 <- plt2 +
                  annotate("text", x = annDat$mid[x], y = 0.5, label = annDat$label[x], hjust = 0.5, 
                           lineheight = 1, size = 3)
      }
      pg <- plot_grid(plt, plt2, align = 'v', axis = 'lr', rel_heights = c(1, 0.15), ncol=1)
    }
  }
  if(!is.null(file)){
    pdf(file, width = 11, height = 8.5)
    grid.draw(pg)
    dev.off()
  } 
  return(pg)
}

prep.for.volcano <- function(res, pvalCol, fcCol, pCut, fcCut, labelCol = NULL){
    vdat <- res %>%
            mutate(`-Log10P` = log10(!!as.name(pvalCol)) * -1,
            passP = ifelse(!!as.name(pvalCol) < pCut, TRUE, FALSE),
            passFC = ifelse(abs(!!as.name(fcCol)) >= fcCut, TRUE, FALSE),
            plotLabel = "")
    if(!is.null(labelCol)){
        if(!labelCol %in% names(vdat)){
            stop(paste0("'", labelCol, "' does not exist in data."))
        }
        vdat <- vdat %>% mutate(plotLabel = ifelse(passP & passFC, !!as.name(labelCol), NA))
    }
    vdat
}

bic.volcano.plot <- function(dat, fcCol, pvalCol = "pval", labelGenes = T, maxSig = 0.05, fcCut = 1, title = "",
                              yLabel = "Significance", xLabel = "Fold Change", pointLabelCol = "plotLabel", 
                              maxUpLabels = 25, maxDnLabels = 25, showCutoffs = TRUE, file = NULL){

    clrs <- c('#cc0000', '#008000', 'darkgray')
    lgndFClbl <- gsub("\\[", "\\(", gsub("\\]", "\\)", fcCol))
    names(clrs) <- c(paste0(lgndFClbl, ' and ', pvalCol), lgndFClbl, "NS") 

    #lgndLabels <- list(bquote(.(xLabel) ~ and ~ .(yLabel)), xLabel, bquote('NS') )
    lgndLabels <- names(clrs)
    names(lgndLabels) <- names(clrs) 

    plotDat <- dat %>%
               mutate(Significance = factor(ifelse(passP & passFC, names(clrs)[1], 
                                             ifelse(passFC, lgndFClbl, "NS")),
                                           levels = names(clrs)))

    ### if not showing cutoffs, no colors needed
    if(!showCutoffs){
        plotDat$Significance <- "NS"
        clrs <- c(NS = "black")
        lgndLabels <- c(NS = "NS")
    }

    cutoffLineSize <- 0.5
    fontSize       <- 12
    pointSize      <- 1.5
    alpha          <- 0.5
    xmax           <- max(abs(plotDat[[fcCol]]))
    ymax           <- max(plotDat[[pvalCol]])
    legendPos      <- ifelse(showCutoffs, "top", "none")
    legendDir      <- "horizontal"
    numPoints <- nrow(plotDat)

    ### make plot
    p <- ggplot(plotDat, aes(x = !!as.name(fcCol), y = !!as.name(pvalCol), 
                             color = Significance, label=!!as.name(pointLabelCol))) +
         geom_point(size = pointSize, alpha = alpha)

    if(showCutoffs){
        p <- p +
             geom_vline(xintercept = fcCut * -1, linetype = "dashed", size = cutoffLineSize) +
             geom_vline(xintercept = fcCut, linetype = "dashed", size = cutoffLineSize) +
             geom_hline(yintercept = maxSig, linetype = "dashed", size = cutoffLineSize)
    }

    p <- p +
         scale_y_continuous(limits = c(0, ymax)) +
         scale_x_continuous(limits = c(xmax * -1, xmax)) +
         scale_color_manual(values = clrs, labels = lgndLabels) 

    p <- p + 
         labs(title = title) +
         ylab(yLabel) +
         xlab(xLabel) +
         bicTheme +
         theme(legend.position = legendPos,
               legend.direction = legendDir,
               legend.justification = c(0, 1),
               legend.title = element_blank(),
               legend.spacing.x = unit(0.2, "cm"),
               legend.text = element_text(margin = margin(l = -7)),
               plot.margin = unit(c(1.5,0.25,1.5,0.25), "in"))

    if(labelGenes && !is.na(pointLabelCol)){
        upLabelDat <- plotDat %>% 
                      filter(!is.na(!!as.name(pointLabelCol)), !!sym(fcCol) > 0) %>%
                      arrange(desc(!!sym(fcCol))) %>%
                      mutate(fcRank = row_number()) %>%
                      arrange(desc(!!sym(pvalCol))) %>%
                      mutate(pRank = row_number()) %>%
                      mutate(min = ifelse(fcRank <= pRank, fcRank, pRank)) %>%
                      arrange(min) %>%
                      filter(row_number() <= maxUpLabels)

        dnLabelDat <- plotDat %>% 
                      filter(!is.na(!!as.name(pointLabelCol)), !!sym(fcCol) < 0) %>%
                      arrange(!!sym(fcCol)) %>%
                      mutate(fcRank = row_number()) %>%
                      arrange(desc(!!sym(pvalCol))) %>%
                      mutate(pRank = row_number()) %>%
                      mutate(min = ifelse(fcRank <= pRank, fcRank, pRank)) %>%
                      arrange(min) %>%
                      filter(row_number() <= maxDnLabels)

        labelDat <- bind_rows(upLabelDat, dnLabelDat)
        p <- p +
             geom_text_repel(data = labelDat,
                             aes(x = !!as.name(fcCol), y = !!as.name(pvalCol), label=!!as.name(pointLabelCol)),
                             color = "black",
                             min.segment.length = unit(2, "mm"))
    }

    ### annotate with total number of genes (points) included on plot
    total <- paste0("n=", formatC(numPoints, big.mark=","))
    p <- p + labs(title = title, caption = total)
    
    ### increase size of points in legend slightly
    p <- p + guides(shape = guide_legend(override.aes = list(size = 6)))

    if(!is.null(file)){ 
        pdf(file, width = 8.5, height = 8.5) 
        print(p)
        dev.off()
    }
    return(p)
}




#' Plot MA
#'
#' Plot log2 fold changes versus mean of normalized counts for a DESeq
#' comparison
#' 
#' @param res   object returned from DESeq function nbinomTest
#' @param file  PDF file to which plot should be saved (optional)
#' @export
bic.plot.ma <- function(res, title = "", file=NULL){

    ## data prep taken directly from DESeq::plotMA

    py <- res %>% filter(baseMean != 0) %>% pull(log2FoldChange)
    ylim <- c(-1, 1) * quantile(abs(py[is.finite(py)]), probs = 0.99) * 1.1
    res <- res %>% 
           mutate(col = ifelse(res$padj >= 0.1, "gray32", "red3"),
                  shape = ifelse(log2FoldChange < ylim[1], 6, ifelse(log2FoldChange > ylim[2], 2, 16)),
                  y = pmax(ylim[1], pmin(ylim[2], log2FoldChange))) %>%
           filter(baseMean != 0)

    clrs <- c("red3", "gray32")
    names(clrs) <- clrs
    pointSize <- 0.5
    p <- ggplot(res, aes(x = baseMean, y = y, color = col, fill = col)) +
         geom_point(data = res %>% filter(shape == 6), shape = 6, size = pointSize) +
         geom_point(data = res %>% filter(shape == 16), shape = 16, size = pointSize) +
         geom_point(data = res %>% filter(shape == 2), shape = 2, size = pointSize) +
         geom_hline(yintercept = 0, size = 1, color = "red3", alpha = 0.5) +
         scale_x_continuous(trans = "log", 
                            labels = scales::scientific, breaks = c(0.1, 10, 1000, 100000)) +
         ylab(expression(log[2] ~ fold ~ change)) +
         xlab("mean of normalized counts") +
         labs(title = title) +
         scale_color_manual(values = clrs) +
         bicTheme +
         theme(legend.position = "none",
               axis.text.y = element_text(color = "black"),
               axis.text.x = element_text(color = "black"),
               plot.margin = margin(t = 0.25, b = 0.25, r = 0.25, l = 0.25, unit = "cm"))

    if(!is.null(file)){
       pdf(file, height = 8.5, width = 11)
       print(p)
       dev.off()
    }
    return(p)
}

#' Histogram of p-values
#' 
#' @param dat  data frame with at least one column with name 'pvals'
#' @param file PDF file to which plot should be saved (optional)
#' @export
bic.pval.histogram <- function(dat, title = "", file=NULL){
  
  p <- ggplot(dat, aes(x = pval)) +
       geom_histogram(binwidth = 0.01, fill = "royalblue3", color = "darkblue") +
       scale_x_continuous(breaks = seq(0, 1, 0.2), expand = c(0.01, 0.01)) +
       scale_y_continuous(expand = c(0.01, 0.01)) +
       ylab("Frequency") +
       xlab("p-value") +
       labs(title = title) +
       bicTheme +
       theme(panel.grid.major = element_blank(),
             panel.grid.minor.x = element_blank(),
             plot.margin = margin(t = 0.25, b = 0.25, r = 0.25, l = 0.25, unit = "cm")) ## these margins help for thumbnails in report

  if(!is.null(file)){
    pdf(file, height = 8.5, width = 11)
    print(p)
    dev.off()
  }
  return(p)

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

#' Size font according to size of device and plot margins
#'
#' Calculate font size in pt using device height and plot margin
#' in inches, and the number of lines needed to fit on plot
#' 
#' @param h       device height in inches
#' @param marg    unit object; plot margins in inches
#' @param nLines  number of lines needed to fit into plot
getFontSizeInPoints <- function(h, marg, nLines){
    inH <- as.numeric(unit(h, 'in') - (marg[1] + marg[3]))
    inchToPt(inH)/nLines
}

getFontSizeInMM <- function(h, marg, nLines){
    inH <- as.numeric(unit(h, 'in') - (marg[1] + marg[3]))
    inchToMM(inH)/nLines
}

inchToPt <- function(inches){
    inches * 72.27
}

ptToInch <- function(pt){
    pt / 72.27
}

inchToMM <- function(inches){
    inches * 25.4
}

mmToPt <- function(mm){
    mm / 0.35
}

ptToMM <- function(pt){
    pt * 0.35
}

getPlotHeightInPoints <- function(h, marg){
    inH <- as.numeric(unit(h, 'in') - (marg[1] + marg[3]))
    inchToPt(inH)
}

getPlotWidthInPoints <- function(w, marg){
    inW <- as.numeric(unit(w, 'in') - (marg[2] + marg[4]))
    inchToPt(inW)
}
#' Write PDF file containing heierarchical clustering tree
#'
#' Quickly plot heirarchical clustering of data with several default
#' parameters like width, height, font size, etc.
#' 
#' @param norm.counts    matrix containing data to plot
#' @param file           file name. Default: NULL
#' @param width          plot width (see plot docs)
#' @param height         plot height (see plot docs)
#'
#' @export
bic.hclust <- function(norm.counts, clrs, file = NULL, height = 8.5, width = 11){

    library(ggdendro)
    dat <- norm.counts %>% dplyr::select(starts_with("s_")) %>% as.matrix %>% t
    dat <- log2(dat + 1)
    hc <- hclust(dist(dat))

    dendDat <- dendro_data(as.dendrogram(hc))

    grps <- sort(unique(gsub(".*__", "", dendDat$labels$label)))
    dendDat$labels <- dendDat$labels %>% 
                      mutate(Group = factor(gsub(".*__", "", label), levels = grps))

    if(is.null(clrs)){
        clrs <- bic.get.colors(length(grps))
        names(clrs) <- sort(grps)
    }

    marg <- unit(c(0.25,1,0.25,1), "in")
    lblFS <- min(14, getFontSizeInPoints(height, marg, nrow(dendDat$labels) * 1.2))
    lgndFS <- min(14, getFontSizeInPoints(height*.75, unit(c(0,0,0,0),"in"), length(grps) ))
    pointSz <- min(4.8, ptToMM(lblFS))
    p <- ggplot(dendDat$segments) + 
           geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
           geom_point(data = dendDat$labels, aes(x = x, y = 0, fill = Group), 
                      size = pointSz, color = "black", shape = 21) +
           labs(title = "Unsupervised clustering", subtitle = "All counts scaled using DESeq method") +
           scale_y_reverse() +
           scale_x_continuous(breaks = dendDat$labels$x, labels = dendDat$labels$label, position = "top") +
           scale_fill_manual(values = clrs) +
           theme_minimal() +
           bicTheme +
           theme(text = element_text(color = "black"), 
                 axis.line.y = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 axis.text.x = element_text(size = lblFS, angle = 90, hjust = 1, color = "black"),
                 axis.text.y = element_text(size = lblFS, color = "black"),
                 axis.title.y = element_blank(),
                 axis.title.x = element_blank(),
                 plot.margin = marg,
                 legend.key.height = unit(lgndFS, "pt"), 
                 legend.text = element_text(size = lgndFS, color = "black"),
                 legend.title = element_blank()) +
           guides(fill = guide_legend(ncol = max(1, ceiling(length(grps)/30)))) +
           coord_flip()

    if(!is.null(file)){ 
        pdf(file, height = height, width = width)
        print(p)
        dev.off()
    }
    return(p)
}


#' Plot MDS clustering of all samples
#'
#' Plot clustering of samples using the MDS method
#' 
#' @param norm.counts  data matrix containing normalized counts, where
#'                     a column represents a sample and a row represents
#'                     a gene
#' @param file         PDF file to which plot should be saved (optional)
#' @param log2         logical indicating whether norm.counts is on log2 scale 
#'                     Default: FALSE
#' @param labels       logical indicating whether to include sample labels on plot; 
#'                     Default: TRUE
#' @export
bic.mds.clust.samples <- function(norm.counts, clrs = NULL, plotLog2 = FALSE, file = NULL, 
                                  labels = TRUE, width = 11, height = 8.5){

  if("ID" %in% names(norm.counts) |
     "GeneSymbol" %in% names(norm.counts)){
    norm.counts <- norm.counts[,-grep("^ID$|^GeneSymbol$", names(norm.counts))]
  }

  norm.counts <- bic.matrix2numeric(as.matrix(norm.counts))

  if(ncol(norm.counts) < 3){
    log_warn("Less than three samples; can not run cluster analysis\n")
    return(NULL) 
  }

  if(plotLog2){
    pseudo <- min(norm.counts[norm.counts > 0])
    counts2hclust <- log2(norm.counts + pseudo)
  } else {
    counts2hclust <- norm.counts
  }

  md <- cmdscale(dist(t(counts2hclust)),2)
  plotDat <- as_tibble(md) %>% 
             mutate(Sample = gsub("__.*", "", rownames(md)), 
                    Group = gsub(".*__", "", rownames(md)))
  grps <- unique(plotDat$Group)
  if(!is.null(clrs)){
      grps <- names(clrs)
  }
  plotDat$Group <- factor(plotDat$Group, levels = grps)

  if(is.null(clrs)){
    clrs <- bic.get.colors(length(levels(plotDat$Group)))
    names(clrs) <- sort(levels(plotDat$Group))
  }

  marg <- unit(c(0.75, 1.5, 0.75, 1.5), "in")
  if(length(levels(plotDat$Group)) > 5){
      marg <- unit(c(0.25, 1.5, 0.25, 1.5), "in")
  }
  fontPts <- getFontSizeInPoints(height, marg, nrow(plotDat))
  lgndRows <- ceiling(length(levels(plotDat$Group))/5)
  lgndFS <- min(14, getFontSizeInPoints(((unit(height, "in") - (marg[1] + marg[3]))/20) * lgndRows, unit(c(0,0,0,0),"in"), lgndRows ))
  pointSz <- min(4, ptToMM(fontPts))

  p <- ggplot(plotDat, aes(x = V1, y = V2, fill = Group, label = Sample), color = "black") +
       geom_point(size = pointSz, shape = 21, alpha = 0.7) +
       scale_fill_manual(values = clrs, labels = names(clrs)) +
       xlab("mds[,1]") +
       ylab("mds[,2]") +
       labs(title = "Multidimensional Scaling (MDS)", subtitle = "All counts scaled using DESeq method") +
       theme_minimal() +
       bicTheme +
       theme(text = element_text(fontPts, color = "black"),
             axis.text.x = element_text(angle = 90, hjust = 1, color = "black"),
             axis.text.y = element_text(color = "black"),
             axis.title.y = element_blank(),
             axis.title.x = element_blank(),
             plot.margin = marg,
             legend.key.size = unit(lgndFS, "pt"), 
             legend.title = element_blank(),
             legend.position = "bottom",
             legend.text = element_text(size = lgndFS, hjust = 0, margin = margin(l = -0.25, unit = "in")),
             legend.spacing.x = unit(0.25, "in")) 
       
  if(labels){
      p <- p +
           geom_text_repel(color = "black", size = 4, point.padding = 0.5)
  }

  if(!is.null(file)){ 
     pdf(file, width = width, height = height) 
     print(p)
     dev.off()
  } 
  return(p)
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
#'                            and a column is a sample. May contain "ID" and/or
#'                            "GeneSymbol" column. If GeneSymbol column is present, 
#'                            it will be used for heatmap labeling. If not, ID 
#'                            column will be used. All data will be included in heatmap.
#' @param genes               a vector of genes to be included in the heatmap;
#'                            if normalized counts matrix includes "GeneSymbol" column,
#'                            genes must be from that column. Otherwise, they must
#'                            be from the "ID" column.                      
#' @param file                name of PDF file to which heatmap should be written. (optional) 
#' 
#' @export
bic.expression.heatmap <- function(norm.counts, genes, file = NULL,
                                  key = NULL, annClrs = NULL, height = 8.5, width = 11, fillFun = NULL,
                                  colorFun = NULL, maxLabels = 100, title = ""){

    dat <- norm.counts

    gnCol <- "GeneSymbol"
    if(!gnCol %in% names(dat) || !any(genes %in% dat[[gnCol]])){
        gnCol <- "ID"
    }

    htmp.dat <- dat %>%
                select(all_of(gnCol), dplyr::matches("^s_")) %>%
                filter_at(vars(all_of(gnCol)), all_vars(. %in% genes)) %>%
                gather(2:ncol(.), key = "Sample", value = "Count") %>%
                group_by_at(c("Sample", gnCol)) %>%
                summarize(Count = mean(Count)) %>%  ## for multiple transcripts, take average count
                separate(Sample, c("Sample", "Group"), sep="__") %>%
                mutate(logCount = log2(Count + 1)) %>%
                group_by_at(gnCol) %>% 
                mutate(MN = mean(logCount)) %>%
                ungroup() %>% 
                mutate(Val = logCount - MN) %>%
                select_at(c(all_of(gnCol), "Sample", "Val"))

    if(is.null(fillFun)){
        fillFun <- scale_fill_gradient2(low = "green", mid = "black", high = "red")
    }

    if(is.null(colorFun)){
        colorFun <- scale_color_gradient2(low = "green", mid = "black", high = "red")
    }

    colAnnot = NULL
    if(is.null(annClrs) && !is.null(key)){
        annClrs <- bic.get.colors(length(unique(key$Group)))
        names(annClrs) <- unique(key$Group)
        colAnnot <- key %>% filter(Sample %in% htmp.dat$Sample)
    }

    htmp <- bic.rnaseq.heatmap(htmp.dat, "Sample", gnCol, 
                               fillFun = fillFun, colorFun = colorFun, annClrs = annClrs,
                               colAnnot = colAnnot, 
                               title = title, 
                               width = width, height = height, file = file,
                               maxLabels = maxLabels)

    return(htmp)
}


bic.rnaseq.heatmap <- function(htmpDat, cols, rows, clusterCols = TRUE, clusterRows = TRUE, colAnnot = NULL, 
                                rowAnnot = NULL, title = "", subtitle = "", width = 11, height = 8.5,
                                fillFun = NULL, colorFun = NULL, annClrs = NULL, file = NULL,
                                maxLabels = 100, labelRows = TRUE, labelCols = TRUE){

    smryDat  <- NULL
    cDat     <- NULL
    cDendDat <- NULL
    rDendDat <- NULL
    cAnnLbls <- NULL
    cAnnP    <- NULL
    rAnnp    <- NULL

    bot_row <- NULL
    right_col <- NULL

    dat <- htmpDat

    cOrder   <- sort(unique(dat[[cols]]))
    rOrder   <- sort(unique(dat[[rows]]))
  
    if(length(rOrder) > maxLabels){
        labelRows <- FALSE
    } else if(length(rOrder) < 3){
        clusterRows <- FALSE
    }
    if(length(cOrder) > maxLabels){
        labelCols <- FALSE
    } else if(length(unique(dat[[cols]])) < 3){
        clusterCols <- FALSE
    }

    if(clusterRows){
        rDat <- dat %>% spread(cols, Val)
        rDendDat <- rDat %>% dist %>% hclust %>% as.dendrogram %>% dendro_data
        rOrder <- rDat[[rows]][as.numeric(rDendDat$labels$label)]
    }

    if(clusterCols){
        cDat     <- dat %>% spread(rows, Val) 
        cDendDat <- cDat %>% dist %>% hclust %>% as.dendrogram %>% dendro_data
        cOrder   <- cDat[[cols]][as.numeric(cDendDat$labels$label)]
    }

    dat[[cols]] <- factor(dat[[cols]], levels = cOrder)
    dat[[rows]] <- factor(dat[[rows]], levels = rOrder)

    smryDat <- dat %>%
               mutate(bin = cut(Val, seq(min(0, min(Val)), round_any(max(Val) + 1, 0.1), 1), include.lowest = T)) %>%
               group_by(bin) %>%
               summarize(Val = n()) %>%
               mutate(x = as.numeric(gsub("\\(|\\[", "", gsub(",.*", "", bin))))

    bgDat <- tibble(x = seq(min(dat$Val), round(max(dat$Val)), 0.1))


    ##########################
    #  Calculate dimensions  #
    ##########################

    ## semi-fixed dimensions first
    rAnnot   <- ifelse(is.null(rowAnnot), 0, ncol(rowAnnot) - 1)
    cAnnot   <- ifelse(is.null(colAnnot), 0, ncol(colAnnot) - 1)

    nRows    <- length(rOrder) + cAnnot
    nCols    <- length(cOrder) + rAnnot
    wPts     <- inchToPt(width)
    hPts     <- inchToPt(height)
    rDendW   <- min(wPts, hPts) * 0.15
    cDendH   <- rDendW * 1.5  ## need to leave room for plot title
    lgndW    <- ifelse(sum(rAnnot, cAnnot) > 0, wPts * 0.15, 0.01)
    smryW    <- rDendW 
    smryH    <- cDendH * 0.5 
    minMargW <- min(inchToPt(0.5), wPts * 0.1)
    minMargH <- min(inchToPt(0.5), hPts * 0.1) 

    maxHtmpW <- wPts - rDendW - lgndW - smryW - (minMargW*2) ### total space width available for heatmap and annotations
    maxHtmpH <- hPts - cDendH - (minMargH * 2)

    boxH     <- min(maxHtmpH/nRows, hPts/50)
    boxW     <- min(maxHtmpW/nCols, wPts/65) 
    maxRlbl  <- ifelse(labelRows, max(nchar(rOrder)), 0)
    maxClbl  <- ifelse(labelCols, max(nchar(cOrder)), 0)
 
    htmpWpts <- max(wPts * 0.18, boxW * nCols + maxRlbl)
    htmpHpts <- boxH * nRows + maxClbl

    ### determine final margins 
    lrMarg <- (wPts - (rDendW + htmpWpts + lgndW))/2
    tbMarg <- (hPts - (cDendH + htmpHpts))/2
    tMarg  <- tbMarg * 0.75
    bMarg  <- tbMarg * 0.25
    smryLR <- (lrMarg + rDendW - smryW)/2 
    smryT  <- tMarg + cDendH - smryH

    top <- tMarg + cDendH

    ## label font sizes
    cFS <- min(14, htmpWpts/(nCols + 2))
    rFS <- min(14, (htmpHpts/(nRows + 2)) * 0.8)
    if(cFS > rFS){
        cFS <- cFS * 0.75
        rFS <- rFS * 1.25
    }
    lgndKeys <- c()

    if(is.null(fillFun)){
        fillFun <- scale_fill_distiller(type = "seq", palette = "GnBu", trans = "reverse")
    }
    if(is.null(colorFun)){
        colorFun <- scale_color_distiller(type = "seq", palette = "GnBu", trans = "reverse")
    }

    blank <- ggplot(data.frame(1)) + theme_void()
    rDend <- cDend <- cAnnP <- rAnnP <- cAnnLbls <- rAnnLbls <- lgnd <- blank 

    ## plots
    smry <- ggplot() +
            geom_vline(data = bgDat,  aes(xintercept = x, color = x), size = 1) +
            geom_line(data = smryDat, aes(x = x, y = Val, group = 1), color = "cyan") +
            scale_x_continuous(expand = c(0,0)) +
            colorFun +
            xlab("Value") +
            ylab("Count") +
            bicTheme +
            theme(text = element_text(size = top * 0.045, color = "black", family = font),
                  axis.text.x = element_text(size = top * 0.045, color = "black", family = font),
                  axis.text.y = element_text(size = top * 0.045, color = "black", family = font),
                  legend.position = "none",
                  panel.border = element_rect(fill = NA, color = "black", size = 0.25),
                  plot.margin = unit(c(smryT, smryLR, 0, smryLR), "pt"))

    htmp <- ggplot(dat, aes_string(x = cols, y = rows)) +
            geom_tile(aes(fill = Val)) +
            scale_y_discrete(position = "right", expand = c(0,0)) +
            scale_x_discrete(expand = c(0,0)) +
            fillFun +
            bicTheme +
            theme(axis.line.x = element_blank(),
                  axis.line.y = element_blank(),
                  axis.text.x = element_text(size = cFS, color = "black", angle = 90, hjust = 1, vjust = 0.5),
                  axis.text.y = element_text(size = rFS, color = "black"),
                  axis.title.x = element_blank(),
                  axis.title.y = element_blank(),
                  legend.position = "none",
                  legend.direction = "vertical",
                  legend.title = element_blank(),
                  plot.margin = unit(c(0, 0, bMarg, 0), "pt"),
                  panel.border = element_rect(size = 1, color = "black", fill = NA))
    if(!labelRows){
        htmp <- htmp + theme(axis.text.y = element_blank())
    }
    if(!labelCols){
        htmp <- htmp + theme(axis.text.x = element_blank())
    }

    ### this is a hacky way of keeping the title on the plot even if not clustering columns
    cDend <- ggplot(cDendDat$segments) +
             labs(title = title, subtitle = subtitle) +
             theme_void() +
             theme(plot.margin = unit(c(tMarg, 0, 0, 0), "pt"),
                   plot.title = element_text(size = hPts * 0.032),
                   plot.subtitle = element_text(face = "italic",
                                                size = hPts * 0.023)) 
    if(clusterCols){
        cDend <- cDend + 
                 geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
                 scale_x_continuous(expand = expansion(add = 0.5)) 
    } 

    if(clusterRows){
        rDend <- ggplot(rDendDat$segments) +
                 geom_segment(aes(x = -y, y = x, xend = -yend, yend = xend)) +
                 scale_y_continuous(expand = expansion(add = 0.5)) +
                 theme_void() +
                 theme(plot.margin = unit(c(0, 0, bMarg, lrMarg), "pt"))
    } 

    lgndDat <- tibble()  ## legend dat will consist of row AND col annotation

    if(!is.null(colAnnot)){
        allAnn <- setdiff(names(colAnnot), cols)

        ann <- colAnnot %>% 
               gather(all_of(allAnn), key = "Annotation", value = "Val") %>%
               mutate(Val = factor(Val, levels = unique(Val)),
                      !!as.name(cols) := factor(!!as.name(cols), levels = cOrder))
        lgndDat <- lgndDat %>% bind_rows(ann)

        lvls <- levels(ann$Val)
        if(is.null(annClrs)){
            annClrs <- bic.get.colors(length(lvls))
            names(annClrs) <- lvls
        }

        legendCols <- 1 #ifelse(lrMarg > inchToPt(2) & length(annClrs) > 4, 2, 1)
        legendFS <- ifelse(legendCols == 1, rFS, rFS * 0.75)

        cAnnP <- ggplot(ann, aes(x = !!as.name(cols), y = Annotation)) +
                 geom_tile(aes(fill = Val), width = 0.95, height = 0.95) + 
                 scale_fill_manual(values = annClrs) +
                 scale_y_discrete(expand = c(0,0)) +
                 scale_x_discrete(expand = c(0,0)) +
                 theme_void() +
                 theme(panel.background = element_rect(fill = "gray"),
                       panel.border = element_rect(size = 0.75, color = "black", fill = NA),
                       legend.position = "none",
                       plot.margin = unit(c(0,lrMarg,0,0), unit = "pt"))

        if(length(allAnn) > 1){
            ## replace blank plot with labels
            cLblDat <- tibble(x = 1, Annotation = allAnn)
            cAnnLbls <- ggplot(cLblDat, aes(x = x, y = Annotation, label = Annotation)) + 
                        scale_x_continuous(expand = c(0,0), limits = c(0.5, 1.01)) +
                        geom_text(size = ptToMM(rFS), hjust = 1, color = "black") + 
                        theme_void() +
                        theme(plot.margin = unit(c(0,0,0,lrMarg), "pt"))
        } 
    }

    if(!is.null(rowAnnot)){
        allRann <- setdiff(names(rowAnnot), rows)

        rAnn <- rowAnnot %>%
                gather(all_of(allRann), key = "Annotation", value = "Val") %>%
                mutate(Val = factor(Val, levels = unique(Val)),
                       !!as.name(rows) := factor(!!as.name(rows), levels = rOrder))
        lgndDat <- lgndDat %>% bind_rows(rAnn)               
        
        lvls <- levels(rAnn$Val)
        if(is.null(annClrs)){
            annClrs <- bic.get.colors(length(lvls))  ### TODO: CHECK FOR THE NEED OF DIFFERENT COLORS THAN COL ANNOTATION
            names(annClrs) <- lvls
        }
        legendCols <- ifelse(lrMarg > inchToPt(2) & length(annClrs) > 4, 2, 1) 

        rAnnP <- ggplot(rAnn, aes(x = Annotation, y = !!as.name(rows))) +
                   geom_tile(aes(fill = Val), width = 0.95, height = 0.95) +
                   scale_fill_manual(values = annClrs) +
                   scale_y_discrete(expand = c(0,0)) +
                   scale_x_discrete(expand = c(0,0)) +
                   theme_void() +
                   theme(panel.background = element_rect(fill = "gray"),
                         panel.border = element_rect(size = 0.75, color = "black", fill = NA),
                         axis.text.x = element_blank(), 
                         legend.position = "none", # c(1.01,1.01),
                         plot.margin = unit(c(0,0,0,0), unit = "pt")) 

        if(length(allRann) > 1){
            rAnnP <- rAnnP + 
                     theme(axis.text.x = element_text(size = cFS, color = "black", angle = 90, hjust = 1))
        }
    }

    if(nrow(lgndDat) > 0){
        lgndDat <- lgndDat %>% select(Val) %>% unique %>% mutate(x = 1, y = rev(row_number()))
        numKeys <- length(levels(lgndDat$Val))
        lFS     <- min(14, (htmpHpts - bMarg)/(numKeys + 2))
        numRows <- (htmpHpts - bMarg)/lFS
        lgndDat <- lgndDat %>% mutate(y = y + (numRows - numKeys))

        lgnd <- ggplot(lgndDat, aes(x = 0, y = y, fill = Val, label = Val)) +
                geom_point(shape = 22, size = lFS * 0.3) +
                geom_text(aes(x = 0.05, label = Val), hjust = 0, size = ptToMM(lFS)) +
                scale_fill_manual(values = annClrs) +
                scale_x_continuous(limits = c(0, 0.5), expand = c(0.1, 0)) +
                scale_y_continuous(limits = c(0, numRows + 0.5), expand = c(0,0)) +
                theme_void() +
                theme(legend.position = "none",
                      plot.margin = margin(0, lrMarg * 0.9, bMarg, lrMarg * 0.1, "pt"))

        ### adjust margins on other plots and apply right margin to legend instead
        htmp <- htmp + theme(plot.margin = unit(c(0,0,bMarg,0), "pt"))
        if(!is.null(cDend)){
            cDend <- cDend + theme(plot.margin = unit(c(tMarg, 0, 0, 0), "pt"))
        }
        if(!is.null(cAnnP)){
            cAnnP <- cAnnP + theme(plot.margin = unit(c(0,0,0,0), "pt"))
        }
    }

    ###### assemble the parts!
    rowAnnW <- ifelse(rAnnot > 0, rAnnot * boxW, 0.01)
    colAnnH <- ifelse(cAnnot > 0, cAnnot * max(boxH,boxW), 0.01)

    rel_widths  <- c(lrMarg + rDendW, rowAnnW, htmpWpts, lgndW + lrMarg)
    rel_heights <- c(tMarg + cDendH, colAnnH, htmpHpts + bMarg)

    #log_debug("rel_widths = ", paste(rel_widths, collapse = ","))
    #log_debug("rel_heights = ", paste(rel_heights, collapse = ","))

    ### align bottom row and right column 
    ## important to start alignment with the column including main heatmap
    mainCol <- align_plots(plotlist = list(cDend, cAnnP, htmp), axis = 'lr', align = 'v') 
    btr <- align_plots(plotlist = list(rDend, rAnnP, mainCol[[3]], blank), axis = 'bt', align = 'h')
    tpr <- align_plots(plotlist = list(smry, blank, mainCol[[1]], blank), axis = 't', align = 'h')

    plotlist = list(tpr[[1]], tpr[[2]], tpr[[3]],     tpr[[4]],
                    blank,    cAnnLbls, mainCol[[2]], blank,
                    btr[[1]], btr[[2]], mainCol[[3]], lgnd)

    pg <- plot_grid(plotlist = plotlist, rel_widths = rel_widths, rel_heights = rel_heights, ncol = 4, nrow = 3) 
                           
    if(!is.null(file)){
        pdf(file, height = height, width = width)
        grid.draw(pg)
        dev.off()
    }
    return(pg)

}




bic.deseq.qc <- function(cds, conds, clustering.dir, idMap = NULL){
    tmp <- tryCatch({
             capture.output(
               bic.high.expression.heatmap(cds,
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
               bic.plot.pc.loading(cds, 
                                  file = file.path(clustering.dir,paste0(pre,"_PC_loadings.pdf")), 
                                  idMap = idMap)
             )
           }, error = function(e){
               log_error("There was an error while plotting PC loadings.")
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
