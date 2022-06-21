library(gtable)

REPORT_CONTENTS <- list(study = c('Hierarchical clustering', 
                                  'MDS', 
                                  'MDS (samples labeled)',
                                  'PCA', 
                                  'PCA (samples labeled)',  
                                  'PC loadings',
                                  'Expression of top 50 genes',  
                                  'Sample to sample distances',
                                  'Cell composition'),
                         comparisons = c('Volcano plots', 'DE heatmap'),
                         comparisons_qc = c('Dispersion Estimates', 'P-value histograms', 'MA plots'))

## some calculations of cell sizes
row_heights <- function(m){
    do.call(unit.c, apply(m, 1, function(l)
      max(do.call(unit.c, lapply(l, grobHeight)))))
}

col_widths <- function(m){
    do.call(unit.c, apply(m, 2, function(l)
      max(do.call(unit.c, lapply(l, grobWidth)))))
}


deseq.cover.page <- function(request, svnrev){
    req <- read.delim(request, header=F, sep=":") %>% deframe %>% trimws %>% as.list

    header <- "Memorial Sloan Kettering Cancer Center\nBioinformatics Core"
    title <- gsub("Proj_", "Project ", req$ProjectID)
    subtitle <- "Differential Expression Analysis Report"
    dt <- format(Sys.Date(), "%Y-%m-%d")

    info <- list(#assay <- paste('Assay:', req$Assay)
                 pipeline <- paste('Pipeline:', req$Pipelines),
                 version <- paste('Version:', svnrev),
                 run <- paste('Pipeline run number:', req$RunNumber),
                 pi <- paste('PI:', req$PI_Name),
                 piEmail <- paste('PI email:', req$`PI_E-mail`),
                 inv <- paste('Investigator:', req$Investigator_Name),
                 invEmail <- paste('Investigator email: ', req$`Investigator_E-mail`))
 
    H1 <- 28
    H2 <- 20
    body <- 14

    bicBlue <- "#3B6AA0"

    spacerGrob <- function(height = unit(1, 'in')){
                      rectGrob(gp = gpar(col = NA), height = height)
                  }
 
    headerGrob   <- grobTree(rectGrob(gp=gpar(fill="#3B6AA0")), 
                             textGrob(header, 
                                      gp=gpar(fontfamily = font, fontsize = H2 * 0.75, col = "white")))
    titleGrob    <- textGrob(title, gp = gpar(fontfamily = font, fontsize = H1))
    subtitleGrob <- textGrob(subtitle, gp = gpar(fontfamily = font, fontsize = H2))
    dateGrob     <- textGrob(dt, gp = gpar(fontfamily = font, fontsize = body))
    infoGrob     <- lapply(seq(length(info)), function(ii)
                             textGrob(info[ii], gp=gpar(fontfamily = font, fontisze = body)))
    footer       <- rectGrob(gp = gpar(fill = bicBlue), height = unit(0.5, 'in')) 

    titleMat <- matrix(list(headerGrob, 
                            spacerGrob(unit(1, 'in')), 
                            titleGrob, 
                            subtitleGrob, 
                            spacerGrob(unit(0.5, 'in')), 
                            dateGrob, 
                            spacerGrob(unit(1, 'in'))), 
                      ncol=1)
    infoMat <- matrix(infoGrob, nc = 1)
    footerMat <- matrix(list(spacerGrob(unit(0.5, 'in')), footer), nc = 1)

    grobs <- rbind(titleMat, infoMat, footerMat)

    linePadding <- unit(3,"mm")
    widths      <- col_widths(grobs) + linePadding
    heights     <- row_heights(grobs)
    heights[-1] <- heights[-1] + linePadding

    ## make height for title grob slightly larger
    heights[3] <- heights[3] * 3

    ## place labels in a gtable
    g <-  gtable_matrix("table", grobs = grobs, widths = widths, heights = heights)

    g
}

squareFigureMargin <- function(figWidth, figHeight, plotHW, unit = 'in'){
    tb = (figHeight - plotHW)/2
    lr = (figWidth - plotHW)/2
    margin(t = tb, b = tb, l = lr, r = lr, unit = unit)
}

sharedPageTitle <- function(title, size = size){
    ggplot() +
    annotate("text", x = 0, y = 0.1, label = title, size = size, hjust = 0, vjust = 1) +
    scale_y_continuous(expand = c(0.1, 1), limits = c(0, 1)) +
    scale_x_continuous(expand = c(0.02, 0), limits = c(0, 1)) +
    theme_void() 
}

sharedXtitle <- function(title, size){
    ggplot() +
    annotate("text", x = 0, y = 0, label = title, size = size, hjust = 0.5, vjust = 0) +
    scale_y_continuous(expand = c(0, 0.2), limits = c(-0.5, 0.1)) +
    theme_void() 
}

sharedYtitle <- function(title, size){
    ggplot() +
    annotate("text", x = 0, y = 0, label = title, size = size, angle = 90, vjust = 1) +
    scale_x_continuous(expand = c(0.1, 1), limits = c(-1, 0)) +
    theme_void() 
}

squarePlotsInPanel <- function(plotPanelWidth, plotPanelHeight, nCols, nRows){
    if(plotPanelWidth == plotPanelHeight){ return(NULL) }

    sqDim <- min(plotPanelWidth/nCols, plotPanelHeight/nRows) - 0.2
    hMarg <- (plotPanelWidth - (sqDim * nCols))/nCols
    vMarg <- (plotPanelHeight - (sqDim * nRows))/nRows

    theme(plot.margin = margin(t = vMarg * 0.5,
                               r = hMarg * 0.25,
                               b = vMarg * 0.5,
                               l = hMarg * 0.25,
                               unit = "in"))
}


combinePlots <- function(plotList, maxPlots = 9, xTitle = "", yTitle = "", pageTitle = "", name = NULL, square = T, 
                          pageWidth = 11, pageHeight = 8.5){

    maxPlots <- min(maxPlots, length(plotList))
    cols <- ceiling(sqrt(maxPlots))
    rows <- ceiling(maxPlots/cols)

    pt <- sharedPageTitle(pageTitle, 7) + theme(plot.margin = margin(t = 0.5, r = 0, b = 0, l = 0.5))
    xt <- sharedXtitle(xTitle, 5)
    yt <- sharedYtitle(yTitle, 5)

    rel_heights <- c(rows*0.1, rows, rows*0.1)
    rel_widths <- c(cols*0.05, cols, cols*0.05)

    plotPanelWidth <- (pageWidth/sum(rel_widths)) * rel_widths[2] 
    plotPanelHeight <- (pageHeight/sum(rel_heights)) * rel_heights[2]    

    if(square){
        plotList <- lapply(plotList, function(x) { 
            x + squarePlotsInPanel(plotPanelWidth, plotPanelHeight, cols, rows)
        })
    }

    pgList <- list()
    numPgs <- ceiling(length(plotList)/maxPlots)
    for(i in 1:numPgs){
        idxs <- seq(maxPlots * (i - 1) + 1, min(maxPlots * i, length(plotList))) 
        pg <- do.call("plot_grid", c(plotList[idxs], nrow = rows, ncol = cols, align = 'hv', axis = 'lrtb'))
        blank <- ggplot(tibble(1)) + theme_void()
        if(i > 1){
            pt <- sharedPageTitle(paste(pageTitle, "(continued)"), 7)
        }
        pgList[[paste0(name, "_", i)]] <- 
               plot_grid(blank, pt, blank, 
                         yt, pg, blank,
                         blank, xt, blank,
                         rel_widths = rel_widths, rel_heights = rel_heights,
                         axis = 'lr', align = 'vh', ncol = 3, nrow = 3)
    }    

    pgList
}


get.deseq.plots <- function(compSetRes, request, svnrev, pre = "TEMP", idMap = NULL, maxPlots = 9, maxLabels = 25, 
                             key = NULL, annColors = NULL, clusterDir = NULL, deFigDir = NULL, qcDir = NULL,
                             cellCompDir = NULL, height = 8.5, width = 11){

    cds <- compSetRes$cds
    nc  <- compSetRes$norm.counts
    labelPlots <- ifelse(ncol(nc) - 2 <= maxLabels, TRUE, FALSE)

    if(is.null(idMap)){
        idMap <- nc %>% select(ID, GeneSymbol)
    }

    ## get plot margin for single plot pages
    sqHW <- 7
    if(!is.null(key) && length(unique(key$Group)) > 5){
        sqHW <- max(height - 0.5, sqHW + 0.1 * length(unique(key$Group))/5)
    }
    sqMargin <- squareFigureMargin(width, height, 7, unit = 'in')
    xtheme <- theme(plot.margin = sqMargin)
    multiplotTheme <- theme(plot.title = element_text(size = 16),
                            axis.title.x = element_blank(),
                            axis.title.y = element_blank())

    report <- list()

    log_debug("creating report cover")
    report$cover <- deseq.cover.page(request, svnrev)

    ## hclustering     
    log_debug("hierarchical clustering")
    report$`Hierarchical clustering` <- bic.hclust(nc, clrs = annColors) 

    ## MDS            
    log_debug("MDS without labels")
    report$`MDS` <- bic.mds.clust.samples(nc, clrs = annColors, plotLog2 = FALSE, labels = FALSE) +  
                    xtheme 
    if(labelPlots){
        log_debug("MDS with labels")
        report$`MDS (samples labeled)` <- bic.mds.clust.samples(nc, clrs = annColors, plotLog2 = FALSE, labels = TRUE) +  
                                          xtheme 
    }

    ## PCA & loadings    
    log_debug("PCA")
    report$`PCA`       <- bic.deseq.plot.pca(cds, labels = FALSE, clrs = annColors) +
                          xtheme
    if(labelPlots){
        report$`PCA (samples labeled)` <- bic.deseq.plot.pca(cds, labels = TRUE, clrs = annColors) +
                                          xtheme
    }

    log_debug("PC loadings")
    report$`PC loadings` <- bic.plot.pc.loading(cds, idMap = idMap)

    ## top 50 heatmap      
    log_debug("top 50 heatmap")
    report$`Expression of top 50 genes` <- bic.high.expression.heatmap(cds, transform = TRUE, 
                                                                       num.gns = 50, idMap = idMap, 
                                                                       key = key, clrs = annColors) 

    ## sample to sample   
    log_debug("sample to sample distances")
    report$`Sample to sample distances` <- bic.sample.to.sample.distances(cds, key = key, clrs = annColors) 

    ## cibersort        
    if(is.null(compSetRes$cibersort.res)){
        log_warn("No CIBERSORT results to plot.")
    } else {
        log_debug("CIBERSORT")
        report$`Cell composition` <- bic.plot.cibersort(compSetRes$cibersort.res) 
    }

    ## per-comparison plots
    comps <- names(compSetRes)[grepl("_vs_", names(compSetRes))]
    for(comp in comps){
        log_debug(comp)

        dat <- tryCatch({ 
                   compSetRes[[comp]]$DE$all.res 
                 }, error = function(x){ 
                     print(x) 
                     log_error("No DE results for comparison [", gsub("_vs_", " vs ", comp), "].") 
                     NULL
               })
        if(is.null(dat)){ next }


        report[[comp]] <- list()

        # volcano plots
        log_debug("  volcano P")
        fcCol <- names(dat)[grepl("log2", names(dat))]
        title <- gsub("_vs_", " vs ", comp)
        pCut <- 10e-11
        minFC <- log2(compSetRes[[comp]]$DE$min.abs.fc)
        max.p <- compSetRes[[comp]]$DE$max.p

        p1 <- prep.for.volcano(dat, "pval", fcCol, pCut, minFC) %>%
              bic.volcano.plot(fcCol, pvalCol = "-Log10P",
                               maxSig = log10(pCut)*-1,
                               fcCut = minFC,
                               title = title,
                               labelGenes = F,
                               yLabel = bquote('-log'[10]~italic('P')),
                               xLabel = bquote('log'[2](.(title))),
                               pointLabelCol = "plotLabel") +
              theme(legend.position = c(-0.05, 1.15),
                    legend.text = element_text(margin = unit(c(0,0.15, 0, 0), "in")),
                    legend.spacing.x = unit(0, "in"),
                    plot.title = element_text(margin = unit(c(0, 0, 0.75, 0), "in")),
                    plot.margin = margin(b = 1.5, t = 1, l = 0.5, r = 0.25, unit = "in"))

        log_debug("  volcano Padj")
        p2 <- prep.for.volcano(dat, "P.adj", fcCol, max.p, minFC) %>%
              bic.volcano.plot(fcCol,
                               labelGenes = F,
                               showCutoffs = FALSE,
                               pvalCol = "-Log10P",
                               maxSig = log10(max.p)*-1,
                               fcCut = minFC,
                               title = title,
                               yLabel = expression('-log'[10]('adjusted'~italic('P'))),
                               xLabel = bquote('log'[2](.(title))),
                               pointLabelCol = "plotLabel",
                               maxUpLabels = 25,
                               maxDnLabels = 25) + 
              labs(title = "") +
              theme(plot.margin = margin(b = 1.5, t = 1, l = 0, r = 0.5, unit = "in"))

        report[[comp]]$`Volcano plots` <- plot_grid(p1, p2, ncol = 2, align = 'h', axis = 'bt') 

        ## pval hist             
        log_debug("  pval histograms")
        report[[comp]]$`P-value histogram` <- bic.pval.histogram(dat, title = title) 

        ## MA plots 
        if(is.null(compSetRes[[comp]]$DE$DESeq)){
            log_warn(paste0("No DESeq results for comparison ", gsub("_vs_", " vs ", comp)))
        } else {
            log_debug("  MA plots")
            report[[comp]]$`MA plot` <- bic.plot.ma(compSetRes[[comp]]$DE$DESeq, title = title)
        }

        ## DE heatmaps      
        if(is.null(compSetRes[[comp]]$DE$filtered) || nrow(compSetRes[[comp]]$DE$filtered) < 3 ){
            log_warn(paste0("Less than three filtered DESeq results for comparison ", gsub("_vs_", " vs ", comp), 
                            ". Can not make heatmap."))
        } else {
            log_debug("  DE heatmap")
            conds <- unlist(strsplit(title, " vs "))
            compKey <- key %>% filter(Group %in% conds)
            genes <- compSetRes[[comp]]$DE$filtered %>% pull(GeneSymbol)
            genes <- genes[1:min(50, length(genes))]
            report[[comp]]$`DE heatmap` <- bic.expression.heatmap(compSetRes$norm.counts, genes, 
                                                                  key = compKey, 
                                                                  annClrs = annColors[names(annColors) %in% compKey$Group],
                                                                  title = title)  
        }
    }

    ## dispersion estimates 
    log_debug("dispersion estimates")  ## TO DO: move this to complete.de.analysis
    report <- c(report,
                   lapply(bic.plot.dispersion.estimates(cds),
                          function(x){
                              x + multiplotTheme
                          }
                   ) %>%
                   combinePlots(maxPlots = maxPlots, xTitle = "mean of normalized counts",
                                yTitle = "dispersion", pageTitle = "Dispersion Estimates",
                                name = "Dispersion Estimates"))

    ## gather pval histograms from all comparisons
    report <- c(report,
                  lapply(comps, function(x){ 
                          if(!is.null(report[[x]]$`P-value histogram`)){
                              report[[x]]$`P-value histogram` + multiplotTheme
                          } else {
                              log_warn("No p-value histogram found for comparison [", 
                                              gsub("_vs_", " vs ", x), "].")
                              NULL
                          } 
                  }) %>%
                  combinePlots(maxPlots = maxPlots, xTitle = "P-value", pageTitle = "P-value Histograms",
                               yTitle = "Frequency", name = "P-value histograms"))


    ### gather MA plots from all comparisons
    report <- c(report,
                lapply(comps, function(x){
                        if(!is.null(report[[x]]$`MA plot`)){ 
                            report[[x]]$`MA plot` + multiplotTheme + theme(legend.position = "none") 
                        } else {
                            log_warn(paste0("No MA plot found for comparison ", gsub("_vs_", " vs ", x)))
                            NULL
                        } 
                }) %>%
                combinePlots(maxPlots = maxPlots, xTitle = "mean of normalized counts", pageTitle = "MA plots",
                             yTitle = expression(log[2] ~ fold ~ change), name = "MA plots")) 

    report
}

addPageNumber <- function(plt, pageNum, fontSize = 10){

    pg <- textGrob(as.character(pageNum), gp=gpar(fontfamily = font, fontsize = fontSize, hjust = 1))
    plot_grid(plt, pg, ncol = 1, align = 'v', rel_heights = c(8.00, 0.5))

}

getReportContents <- function(report){

    report <- report[!names(report) == "cover"]
print(names(report))
    ##
    ## print overall study plots first
    ##
    #contents <- c('Hierarchical clustering', 'MDS', 'MDS (samples labeled)',
    #              'PCA', 'PCA (samples labeled)', 'PC loadings',
    #              'Expression of top 50 genes', 'Sample to sample distances', 
    #              'Cell composition') 

    contents <- REPORT_CONTENTS$study

    tbl <- tibble()
    nxtPg <- 1
    for(pt in contents){
        #ptrn <- gsub("\\(", "\\\\(", gsub("\\)", "\\\\)", pt))
        #idxs <- grep(ptrn, names(report))
        # these should be exact matches
        idxs <- names(report)[names(report) == pt]
        if(!is.null(idxs) && length(idxs) > 0){
            tbl <- tbl %>%
                   bind_rows(tibble(Name = pt, Page = nxtPg, listName = names(report)[idxs[1]]))
            nxtPg <- nxtPg + length(idxs)
        }
    } 
    
    ##
    ## comparison level plots
    ##
    #contents <- c('Volcano plots', 'DE heatmap')
    contents <- REPORT_CONTENTS$comparisons

    comps <- names(report)[grep("_vs_", names(report))]
    for(comp in comps){
        plts <- contents[contents %in% names(report[[comp]])]
        nms <- paste(comp, plts)
        pgs <- seq(plts) + (nxtPg - 1) #tbl$Page[nrow(tbl)]
        tbl <- tbl %>% bind_rows(tibble(Name = nms, Page = pgs, listName = plts))
        nxtPg <- nxtPg + length(plts)
    }

    ##
    ## QC plots
    ##
    #contents <- c('Dispersion Estimates', 'P-value histograms', 'MA plots')
    contents <- REPORT_CONTENTS$comparisons_qc
    for(pt in contents){
        ## these report names will NOT be an exact match like the overall study plots
        ## because there are multiple plots with names that each start with pt
        ## e.g., 'Dispersion Estimates_1', 'Dispersion_estimates_2'
        idxs <- grep(paste0("^", pt, "_"), names(report)) 
        if(!is.null(idxs) && length(idxs) > 0){
            tbl <- tbl %>%
                   bind_rows(tibble(Name = pt, Page = nxtPg, listName = names(report)[idxs[1]]))
            nxtPg <- nxtPg + length(idxs)
        }
    } 

    tbl
}

reportTableOfContents <- function(toc, maxItems = 28){

    tbl <- toc %>%
           select(-listName) %>%
           mutate(tocPg = ceiling((1:nrow(toc))/maxItems)) %>%
           group_by(tocPg) %>%
           mutate(y = maxItems - row_number() + 1) %>%
           gather(c('Name', 'Page'), key = 'col', value = 'value') %>%
           mutate(col = ifelse(col == "Name", 1, 2), 
                  hjust = ifelse(col == 1, 0, 1))
  
    tocList <- list()

    for(pg in unique(tbl$tocPg)){
        tocList[[pg]] <- 
            ggplot(tbl %>% filter(tocPg == pg), 
                   aes(x = col, y = y, label = value, hjust = hjust)) +
              geom_text(size = 4) +
              scale_x_continuous(limits = c(1, 2.1), expand = c(0.1,0.1)) +
              scale_y_continuous(limits = c(-1 * (7.5 - (0.2 * maxItems)), maxItems + 1)) +
              labs(title = ifelse(pg == 1, "Report Contents", "Report Contents (continued)")) +
              theme_void() +
              theme(plot.title = element_text(hjust = 0.5, size = 20), 
                    plot.margin = margin(1, 1, 1, 1, unit = "in"))
    }
    tocList 
}

bic.deseq.report <- function(compSetRes, request, svnrev, pdfFile, pre = "TEMP", maxPlots = 9, idMap = NULL, 
                               key = NULL, annColors = NULL, clusterDir = NULL, deFigDir = NULL,
                               qcDir = NULL, cellCompDir = NULL){ 

    if(is.null(annColors) && !is.null(key)){
        grps <- unique(key$Group)
        annColors <- bic.get.colors(length(grps))
        names(annColors) <- grps
    }
    report <- tryCatch({
                  get.deseq.plots(compSetRes, request, svnrev, pre = pre, idMap = idMap, maxPlots = maxPlots, 
                                   key = key, annColors = annColors, 
                                   clusterDir = clusterDir, deFigDir = deFigDir, qcDir = qcDir, cellCompDir = cellCompDir)
                }, error = function(){
                  stop("REPORT NOT CREATED.")
              })

    log_debug("Compiling report contents")
    toc <- getReportContents(report)
    report$toc <- reportTableOfContents(toc)

    pdf(pdfFile, height = 8.5, width = 11, onefile = T)
    showtext_begin()
    grid.draw(report$cover)
    for(i in 1:length(report$toc)){
        grid.draw(report$toc[[i]])
    }
    pgNum <- 1
    for(i in 1:nrow(toc)){
        pltInfo <- toc[i,]
        ptrn <- paste0("^", gsub("\\(", "\\\\(", gsub("\\)", "\\\\)", pltInfo$Name)), "(_\\d+|$)")
        plts <- report[grep(ptrn, names(report))]

        if(length(plts) == 0 && grepl("_vs_", pltInfo$Name)){
            comp <- gsub(" .*", "", pltInfo$Name)
            ptrn <- gsub(".* ", "", pltInfo$Name)
            plts <- report[[comp]][grep(ptrn, names(report[[comp]]))]
        }
        for(plt in plts){
            tryCatch({
                grid.draw(plt %>% addPageNumber(pgNum))
                pgNum <- pgNum + 1
            }, error = function(e){
               log_warn("Could not print ", pltInfo$listName)
            })
        }
    }
    showtext_end()
    dev.off()

}
