#' Run CIBERSORT
#'
#' Run CIBERSORT to get cell type proportions in each sample
#'
#' @param counts  table of counts for each sample
#' @param outFile path to output file
#' @param javeExe java executable
#' @param jar     CIBERSORT jar file
#' @param sigFile signatures file (downloaded from cibersort site)
#' @export
bic.run.cibersort <- function(counts, outFile, cibersortR, sigFile = NULL, plot = TRUE){
    # format table
    if(!"GeneSymbol" %in% names(counts)){
        log_error("No 'GeneSymbol' column found in normalized counts file. Can not run CIBERSORT.")
        stop("No 'GeneSymbol' column found in normalized counts file. Can not run CIBERSORT.")
    } 
    if("ID" %in% names(counts)){
        counts <- counts %>% select(-ID)
    }
    dat <- counts %>% select(GeneSymbol, everything())

    # save table to tmp file
    of <- tempfile(pattern = "tmp", fileext = ".txt")
    bic.write.dat(dat, file = of)

    source(cibersortR)
    results <- CIBERSORT(sigFile, of)
    smps <- rownames(results)
    results <- results %>%
               as_tibble() %>% 
               mutate(Sample = smps)
    rm(of)
    bic.write.dat(results, file = outFile)

    if(plot){
        bic.plot.cibersort(results, file = gsub(".txt", ".pdf", outFile), annotate = FALSE)
        #bic.plot.cibersort(results, file = gsub(".txt", "_annotated.pdf", outFile), 
        #                   annotate = TRUE, sortByGroup = TRUE)
    }
    return(results)
}

#' Run and save results of a single GSA comparison
#'
#' Run and save results of a single GSA comparison
#'
#' @param de.res     DESeq results
#' @param comp.name  comparison name (e.g., "GroupB_vs_GroupA")
#' @param species    character, either 'mouse' or 'human'
#' @param gmt.dir    path to GMT files
#' @param gsa.dir    output directory
#'
#' @return two-element list of GSA results, one for 'up' results and one for 'down' 
#'         results
#' @export
bic.complete.gsa.analysis <- function(de.res, comp.name, species, gmt.dir, gsa.dir){

    gsa.res <- bic.run.gsa(species, de.res, gmt.dir)

    bic.write.dat(gsa.res$dn,
                  file = file.path(gsa.dir,
                                   paste0("GeneSet_Dn_", comp.name, ".xlsx")))

    bic.write.dat(gsa.res$up,
                  file = file.path(gsa.dir,
                                   paste0("GeneSet_Up_", comp.name, ".xlsx")))

    gsa.res

}

#' Run DESeq comparison, write results to XLSX and save visualizations in PDF files
#'
#' For a single comparison, run DESeq, save all and filtered results in XLSX files, 
#' and save PCA, MA and expression heatmap figures in PDF files.
#'
#' @param cds               countDataSet object
#' @param counts            tibble of normalized counts
#' @param conds             vector of sample conditions, in same order as columns in counts tibble
#' @param all.gene.dir      directory where table of unfiltered DESeq results should be saved
#' @param diff.exp.dir      directory where table of filtered DESeq results should be saved
#' @param diff.exp.fig.dir  directory where DESeq figures should be saved
#' @param condA             condition A (denominator of comparison)
#' @param condB             condition B (numerator of comparison)
#' @param geneSymbols       optional; two-column tibble with columns ID and GeneSymbol
#' @param max.p             maximum p.value to pass filter
#' @param min.abs.fc        minumum absolute fold change to pass filter
#' @param min.count         minimum mean gene counts
#' @param zeroaddQ          keep genes with zero count in one condition
#' @param orderPvalQ        sort filtered results by adjusted p value
#' @param heatmaps          logical; when TRUE (default), generate heatmap of up to 100 most
#'                          differentially expressed significant conditions
#'
#' @return list of all DESeq results
#' @export
bic.complete.de.analysis <- function(cds, counts, conds, all.gene.dir, diff.exp.dir,
                                     diff.exp.fig.dir = NULL, qc.dir = NULL,
                                     condA = NULL, condB = NULL, geneSymbols = NULL, key = NULL, 
                                     max.p = 0.05, min.abs.fc = 2, min.count = 10,
                                     zeroaddQ = TRUE, orderPvalQ = TRUE, heatmaps = TRUE,
                                     maxPlotLabels = 25, annClrs = NULL){
    ## run DESeq comparison
    de.res <- bic.run.deseq.comparison(cds, conds, condA, condB,
                                       max.p = max.p,
                                       min.abs.fc = min.abs.fc,
                                       min.count = min.count,
                                       zeroaddQ = zeroaddQ,
                                       genes = geneSymbols)

    ##
    ## write All DE results to file
    ##
    log_debug("    Writing results for ALL genes to file...")
    file.name <- file.path(all.gene.dir, paste0("ALLResDESeq_",condB,"_vs_",condA,".xlsx"))
    bic.write.deseq.results(de.res$all.res, file.name, condA, condB, orderPvalQ)
    log_debug("Done.\n")

    ##
    ## format and write DE results to file
    ##
    log_debug("    Writing results for DE genes to file...")
    file.name <- file.path(diff.exp.dir, paste("ResDESeq_",condB,"_vs_",condA,".xlsx",sep=""))
    bic.write.deseq.results(de.res$filtered, file.name, condA, condB, orderPvalQ)
    log_debug("Done.\n")

    ##
    ## DESeq visualization
    ##
    fcCol <- names(de.res$all.res)[grepl("log2", names(de.res$all.res))]
    comp  <- gsub("log2\\[", "", gsub("\\]", "", fcCol))
    title <- gsub("/", " vs ", comp)

    if(!is.null(qc.dir) | !is.null(diff.exp.fig.dir)){
        log_debug("    Drawing DESeq plots...")
    }

    if(!is.null(qc.dir)){
log_debug("Writing QC plots to ", qc.dir)
        bic.plot.ma(de.res$DESeq,
                    file=file.path(qc.dir,
                                   paste0(pre, "_MAplot_", gsub(" ", "_", title), ".pdf")
                         ))
        bic.pval.histogram(de.res$DESeq,
                           title = title,
                           file = file.path(qc.dir,
                                            paste0(pre, "_pval_histogram_", gsub(" ", "_", title), ".pdf")
                           ))
        bic.plot.dispersion.estimates(cds, out.dir = "", 
                                      file.prefix = file.path(qc.dir, paste0(pre, "_dispersion_estimates")))
    }

    if(!is.null(diff.exp.fig.dir)){

        pCut  <- 10e-11
        file.name <- file.path(diff.exp.fig.dir,
                               paste0(pre, "_volcano_raw_p_", gsub("/", "_vs_", comp), ".pdf"))
        prep.for.volcano(de.res$all.res, "pval", fcCol, pCut, log2(min.abs.fc)) %>%
        bic.volcano.plot(fcCol,
                         pvalCol = "-Log10P",
                         maxSig = log10(pCut)*-1,
                         fcCut = log2(min.abs.fc),
                         title = title,
                         yLabel = bquote('-log'[10]~italic('P')),
                         xLabel = bquote('log'[2](.(comp))),
                         pointLabelCol = "plotLabel",
                         labelGenes = F,
                         maxUpLabels = maxPlotLabels,
                         maxDnLabels = maxPlotLabels,
                         file = file.name)

        file.name <- file.path(diff.exp.fig.dir, 
                               paste0(pre, "_volcano_adj_p_", gsub("/", "_vs_", comp), ".pdf"))
        prep.for.volcano(de.res$all.res, "P.adj", fcCol, max.p, log2(min.abs.fc)) %>%
        bic.volcano.plot(fcCol,
                         showCutoffs = FALSE,
                         pvalCol = "-Log10P",
                         maxSig = log10(max.p)*-1,
                         fcCut = log2(min.abs.fc),
                         title = title,
                         yLabel = expression('-log'[10]('adjusted'~italic('P'))),
                         xLabel = bquote('log'[2](.(comp))),
                         labelGenes = FALSE,
                         pointLabelCol = "plotLabel",
                         maxUpLabels = maxPlotLabels,
                         maxDnLabels = maxPlotLabels,
                         file = file.name)

        ##
        ## heatmap of top DE genes
        ##
        if(heatmaps & !is.null(de.res$DEgenes) & length(de.res$DEgenes) > 0){ 
            file.name <- file.path(diff.exp.fig.dir,
                                   paste0(pre,"_DE_heatmap_",condB,"_vs_",condA,".pdf"))
            tryCatch({
                genes <- de.res$filtered %>% pull(GeneSymbol)
                genes <- genes[1:min(50, length(genes))]
                htmpDat <- counts %>% 
                           select(GeneSymbol, 
                                  grep(paste(paste0("__", c(condA, condB), "$"), collapse = "|"), names(.)))
                bic.expression.heatmap(htmpDat,
                                       genes = genes, 
                                       key = key %>% filter(Group %in% conds),
                                       title = paste0(condB, "_vs_", condA),
                                       file = file.name, 
                                       annClrs = annClrs)
              }, error = function(e){
                  traceback()
                  warning(paste0("Could not generate heatmap for ",condB," vs ",condA))
            })
        }

    }

    return(de.res) 
}


#' Read in raw HTSeq count file and format for analysis
#'
#' Reformatting includes removing samples that are not in 
#' sample key (if key is provided), rounding counts to
#' nearest integer, assigning IDs to rownames of matrix,
#' and removing ID column from matrix if it also contains
#' a GeneSymbol column.
#'
#' @param htseq.files tab-delimited file containing raw HTseq
#'                    counts; columns are samples and rows are
#'                    genes; must contain ID column but may 
#'                    also contain a GeneSymbol column
#' @param key         two-column matrix where column one contains
#'                    all samples to be included in analysis and
#'                    column two contains condition groups to which
#'                    the samples belong
#' @return a list of two items: 1) a matrix of raw counts formatted 
#'         for BIC differential expression analysis and 2) a vector
#'         of gene symbols named by the ID in htseq file (this
#'         vector will be NULL if there is no GeneSymbol column in
#'         the htseq file.
#' @export
bic.format.htseq.counts <- function(htseq.file,key=NULL){

  log_debug("Reformatting raw counts...")

  HTSeq.dat <- read.csv(htseq.file, header = T, sep = "\t", check.names = F) %>%
               as_tibble()

  names(HTSeq.dat)[1] <- "ID" 

  HTSeq.dat <- HTSeq.dat %>%
               filter(!grepl("alignment_not_unique|ambiguous|no_feature|not_aligned|too_low_aQual", ID))

  ## subset/sort data according to key
  if(!is.null(key)){
    smps <- key$Sample
    if(!all(smps %in% names(HTSeq.dat))){
        print(paste0("WARNING: These samples are NOT in count data... ", 
                     paste(smps[!smps %in% names(HTSeq.dat)], collapse=", ")))
        smps <- smps[smps %in% names(HTSeq.dat)] 
    }
    HTSeq.dat = HTSeq.dat %>% select_at(c("ID", "GeneSymbol", smps))
  } 

  formatted.counts <- list()
  formatted.counts$raw <- HTSeq.dat %>% mutate_if(is.numeric, round)
  formatted.counts$ids <- HTSeq.dat %>% select(ID, GeneSymbol)

  return(formatted.counts)
}

#' Non-deseq way to normalize counts(?)
#'
#' TO DO: INSERT DESCRIPTION HERE
#' 
#' @param raw.counts  Raw count matrix, where a column is one sample and a row is one gene
#' @param percentile  A string representing percentile; e.g., "100\%". Default: "75\%"
#' @return TO DO
#' @export
bic.scale.factor.quantile <- function(raw.counts,percentile="75%"){
  qu <- quantile(raw.counts, probs = seq(0,1,0.05))
  qu <- qu[grep(percentile, names(qu))]
  return(sum(raw.counts[which(raw.counts <= q)]))
}

#' Run a differential expression comparison between two
#' conditions.
#'
#' Compare expression in condition B vs expression in condition B. Results include
#' log2FC, adjusted p-val, and mean expression in each condition. Two result sets
#' are included in returned list: one including ALL genes, and one including only 
#' genes differentially expressed as determined by \code{max.p} and \code{min.abs.fc}. 
#' User has the option to include in the results those genes with zero counts in one
#' condition and count greater than \code{min.count/mean(sizeFactors)} in the other
#' condition, even if the p-value is greater than \code{max.p}, as those genes may still
#' be of interest.
#'
#' @param countDataSet  A countDataSet object created by running DESeq's newCountDataSet()
#' @param conds         Vector of sample conditons, in the same order as samples 
#'                      in raw count matrix
#' @param condA         The first condition in comparison (must be in conditions vector)
#' @param condB         The second condition in comparison (must be in conditions vector)
#' @param max.p         A number, the max p-value cutoff to determine which genes are 
#'                      considered differentially expressed; Default: 0.05
#' @param min.abs.fc    A number, the minimum absolute fold change cutoff to determine 
#'                      which genes are considered differentially expressed; Default: 2
#' @param min.count     A number, minimum average number of reads in one condition; 
#'                      Default: 10
#' @param zeroaddQ      Logical, indicates whether to include genes that have zero 
#'                      counts in one condition and insignificant p-value, but at 
#'                      least min.count/mean(sizeFactors) in the other condition; 
#'                      Default: FALSE
#' @param genes         A named vector of gene symbols, where the names are equivalent 
#'                      to the IDs in the original raw count matrix (e.g., Ensembl IDs); 
#'                      Default: NULL
#' @return List containing a matrix of ALL results, a matrix containing results
#'         for differentially expressed genes (as determined by cutoffs), a vector
#'         of DE genes, the cutoffs used in analysis, and norm factors (TO DO: reword this) 
#' @export
bic.run.deseq.comparison <- function(countDataSet, conds, condA, condB,
                                     max.p=0.05, min.abs.fc=2, min.count=10,
                                     zeroaddQ=F, genes=NULL){

  log_debug(paste0("Comparing [", condB, "] vs [", condA, "] ... "))

  cds <- countDataSet
  ng  <- c()  ## significant genes 
  
  rawRes <- tryCatch({ 
                nbinomTest(cds, condA, condB) %>% as_tibble() %>% rename(ID = id) 
              }, error = function(e){
                print(e)
                return(NULL)
            }) 
         
  if(is.null(rawRes) || nrow(rawRes) == 0){ return(NULL) }

  res <- rawRes
  if(!is.null(genes)){
    res <- res %>% 
           left_join(genes, by = intersect(names(.), names(genes))) %>%
           select(ID, dplyr::matches("Gene"), everything())
  }

  norm.factors <- sizeFactors(cds)

  ## get all significant genes
  ng <- res %>%
        filter(padj <= max.p, 
               abs(log2FoldChange) >= log2(min.abs.fc),
               (baseMeanA >= min.count/mean(norm.factors) | 
                  baseMeanB >= min.count/mean(norm.factors))) %>%
        pull(ID)

  if(zeroaddQ){
    jj <- res %>%
          filter((baseMeanA == 0 & baseMeanB >= min.count/mean(norm.factors)) |
                 (baseMeanB == 0 & baseMeanA >= min.count/mean(norm.factors))) %>%
          pull(ID)
    ng <- unique(c(ng, jj))
  }

  allRes <- res %>%
            mutate(log2FoldChange = ifelse(log2FoldChange %in% c(Inf,-Inf), 
                                             log2((baseMeanB + 1)/(baseMeanA + 1)), 
                                             log2FoldChange)) %>%
            select(dplyr::matches("ID|Gene"),
                   pval, 
                  `P.adj` = padj,
                  !!as.name(paste0("log2[", condB, "/", condA, "]")) := log2FoldChange,
                  !!as.name(paste0("Mean_at_cond_", condA)) := baseMeanA,
                  !!as.name(paste0("Mean_at_cond_", condB)) := baseMeanB) %>%
            arrange(desc(!!as.name(paste0("log2[", condB, "/", condA, "]")))) %>%
            filter(!!as.name(paste0("Mean_at_cond_", condA)) >= min.count | 
                    !!as.name(paste0("Mean_at_cond_", condB)) >= min.count)

  ## prepare filtered output
  resSig <- allRes %>% filter(ID %in% ng) 

  if(nrow(resSig) == 0){
    log_debug("\n===================================================\n")
    log_debug(paste0("No genes pass the significant cutoff of ", max.p,
              "\n  AND have sufficient mean number of reads across samples\n AND at least ", 
              min.count, " reads in one condition\n",
              "===================================================\n"))
    log_debug(paste0("Largest log2FC = ", max(abs(res$log2FoldChange)), "\n"))
    log_debug(paste0("Best corrected p.value = ", min(res$padj), "\n"))
    log_debug("\n\n")
  }

  list(DESeq = rawRes,
       filtered = resSig,
       DEgenes = resSig$ID,
       all.res = allRes,
       max.p = max.p,
       min.abs.fc = min.abs.fc,
       min.count = min.count,
       norm.factors = norm.factors)
  
} 



#' Normalize raw count matrix using DESeq scaling method
#' 
#' @param raw.counts  Raw count matrix, where each column is a sample and each
#'                    row is a gene 
#' @param conds       vector of sample conditions, in the same order as
#'                    column names in \code{raw.counts}
#' @param min.count   minimum total number of reads for a gene across samples; Default: 10
#' @param libsizeQ    logical to indicate we should use library size to normalize counts 
#'                    rather than DESeq's method [TO DO: RECHECK THIS]
#' @param percentile  string indicating percentile to use when normalizing by library size
#' @param fitType     See DESeq documentation. Default: "parameteric"
#' @param method      See DESeq documentation. Default: "per-condition"
#' @param sharingMode See DESeq documentation. Default: "maximum"
#' @return DESeq's countDataSet object
#' @export
bic.get.deseq.cds <- function(raw.counts, conds,
                              min.count=10, libsizeQ=F, percentile="100%",
                              fitType="parametric", method="per-condition",
                              sharingMode="maximum"){

  log_debug("Getting DESeq countDataSet...")

  counts.tmp <- raw.counts %>%
                select(ID, dplyr::matches("^s_")) %>%
                mutate_at(vars(matches("^s_")), funs(as.integer(.))) %>%
                mutate(Total = rowSums(.[,2:ncol(.)])) %>%
                filter(Total >= min.count) %>%
                select(-Total)

  ids <- counts.tmp$ID
  mat <- as.matrix(counts.tmp[,-1])
  rownames(mat) <- ids

  ## create CountDataSet, DESeq's main data structure
  cds <- newCountDataSet(mat, conds)

  ## by default, normalize counts using DESeq's method
  ## if otherwise specified, normalize by library size
  if(!libsizeQ){
    cds = estimateSizeFactors(cds)
  }

  if(libsizeQ | length(which(is.na(sizeFactors(cds))))>0){
    log_debug("\nEstimating sizeFactors using library size\n")
    if(percentile == "100%"){
      libsizes <- apply(mat, 2, sum)
    } else {
      libsizes <- apply(mat, 2, util.scale.factor.quantile, percentile=percentile)
    }
    libsizes <- libsizes / median(libsizes)
    sizeFactors(cds) <- libsizes
  }
  
  cds <- tryCatch({
           estimateDispersions(cds, fit=fitType, method=method, sharingMode=sharingMode)
         }, error = function(err){
           if(!method == "pooled"){
             warning(paste("method='",method,"' did not work. Trying method='pooled'",sep=""))
             cds <- tryCatch({
                  estimateDispersions(cds, fit=fitType, method='pooled', sharingMode=sharingMode)
               }, error = function(x){
                  if(!method == "blind"){
                    warning(paste("method='",method,"' did not work. Trying method='blind'",sep=""))
                    cds <- tryCatch({
                         estimateDispersions(cds, fit=fitType, method='blind', sharingMode='fit-only')
                    }, error = function(x){
                         warning(paste("All methods '",method,"', 'pooled' and 'blind' failed.",sep=""))
                         return(NULL)
                    })
                  }
               })
           } else {
             warning("method 'pooled' did not work. Please try another method.")
             return(NULL)
           }   
         })
  log_debug("Done.\n")
  return(cds) 
}


bic.standard.clustering <- function(cds, scaled.counts, out.dir, idMap = idMap, pre = "TEMP", annColors = NULL){

    file.name <- file.path(out.dir, paste0(pre,"_counts_scaled_hclust.pdf"))
    tryCatch({
        tmp <- bic.hclust(scaled.counts, annColors, file = file.name)
      }, error = function(){
          log_error("Error running hclust")
      })

    file.name <- file.path(out.dir, paste0(pre,"_MDS.pdf"))
    file.name2 <- gsub("MDS", "MDS_labeled", file.name)
    tryCatch({
        tmp <- bic.mds.clust.samples(scaled.counts, file = file.name, labels = FALSE, clrs = annColors)
        tmp <- bic.mds.clust.samples(scaled.counts, file = file.name2, labels = TRUE, clrs = annColors)
      }, error = function(){
          log_error("Error running MDS clustering") 
      })

    file.name <- file.path(out.dir, paste0(pre, "_PCA.pdf"))
    file.name2 <- gsub("PCA", "PCA_labeled", file.name)
    file.name3 <- gsub("PCA", "PC_loadings", file.name)
    tryCatch({
        tmp <- bic.deseq.plot.pca(cds, file = file.name, labels = FALSE, clrs = annColors)
        tmp <- bic.deseq.plot.pca(cds, file = file.name2, labels = TRUE, clrs = annColors)
        tmp <- bic.plot.pc.loading(cds, idMap = idMap, file = file.name3) 
    }, error = function(err){
        log_error("Error running PCA")
    })

}


#' Write normalized counts matrix to file
#'
#' @param  norm.counts.mat     matrix containing normalized counts
#' @param  file.name           File to which results should be written
#' @export
bic.write.normalized.counts.file <- function(norm.counts.mat,file.name){
    if(grepl("\\.txt$", file.name)){
        write.table(norm.counts.mat,
                    file      = file.name, 
                    quote     = F,
                    col.names = T,
                    row.names = F,
                    sep       = "\t")
    } else if(grepl("\\.xlsx$", file.name)){
        openxlsx::write.xlsx(as.data.frame(norm.counts.mat), file=file.name)
    } else {
        stop(paste0("Unrecognized file type: .", gsub(".*\\.", "", file.name)))
    }
}

#' Write results for ALL genes and for DE genes
#'
#' Two files are written: One containing log2FC, adjusted p-val, mean(condA),
#' mean(condB) for ALL genes, and one containing all of that same information
#' for ONLY DE genes. Results are sorted by abs(log2FC) by default but may be
#' sorted by adjusted p-val. 
#'
#' @param res         Matrix to write to file, either all.res or filtered, 
#'                    from \code{bic.get.norm.deseq.counts}
#' @param file.name   File to which results should be written
#' @param orderPvalQ  Logical indicating results should be sorted by pValue instead of FC 
#' @return Nothing.   Files are written.
#' @export
bic.write.deseq.results <- function(res, file.name, condA, condB, orderPvalQ=FALSE){

  if(is.null(res) || nrow(res) == 0){
      bic.write.empty.results(file.name, condA, condB)
      return()
  }

  ## Sort results  
  if(orderPvalQ){
    res <- res[order(res[,"P.adj"],decreasing=F),]
  } else {
    res <- res[order(abs(res[,grep("log",colnames(res))]),decreasing=T),]
  }
  bic.write.dat(res, file=file.name) 
}

#' Use DESeq to normalize raw HTSeq counts 
#' 
#' Use DESeq scaling method to normalize raw HTSeq counts and return a list containing both 
#' raw and normalized counts. Takes in either formatted.counts matrix from bic.format.htseq.counts()
#' OR an HTSeq counts file. Optionally takes in previously generated countDataSet. If not given,
#' also optionally takes in parameters for generating one. Sample key will be used in generation of
#' countDataSet but if not given, it will be assumed that all samples belong to one condition group.
#'
#' @param formatted.counts  matrix containing raw counts from HTSeq. One column per sample,
#'                          one row per gene; may be used INSTEAD of htseq.file, but NOT in
#'                          conjunction with
#' @param htseq.file        file containing raw HTSeq counts; One column per sample, one row
#'                          per gene; may be used INSTEAD of formatted counts, but NOT in 
#'                          conjunction with
#' @param cds               DESeq countDataSet object; Default: NULL, if null will compute
#'                          on the fly
#' @param key               sample key where column one contains sample ID that matches
#'                          that in counts file, and column two contains the condition 
#'                          group to which it belongs. Only those samples that are in the 
#'                          key matrix will be included in the returned matrices
#' @param min.count         minimum total number of reads for a gene across samples; Default: 10
#' @param libsizeQ          logical to indicate we should use library size to normalize counts 
#'                          rather than DESeq's method [TO DO: RECHECK THIS]
#' @param percentile        string indicating percentile to use when normalizing by library size
#' @param fitType           See DESeq documentation. Default: "parameteric"
#' @param method            See DESeq documentation. Default: "per-condition"
#' @param sharingMode       See DESeq documentation. Default: "maximum"#' 
#' @return list containing two matrices: one of raw counts and one of
#'         DESeq-scaled, log2 transformed counts. Column names will have condition appended
#' @export
bic.deseq.normalize.htseq.counts <- function(formatted.counts=NULL,htseq.file=NULL,cds=NULL,key=NULL,
                                             min.count=10, libsizeQ=F, percentile="100%",
                                             fitType="parametric", method="per-condition",
                                             sharingMode="maximum"){

  if(!is.null(formatted.counts) & !is.null(htseq.file)){ 
    warning("Both formatted.counts and htseq.file are not NULL. Using htseq file for raw counts")
    formatted.counts=NULL
  }  

  if(is.null(formatted.counts)){
    formatted.counts <- bic.format.htseq.counts(htseq.file, key)
  }

  idsAndGns <- formatted.counts$ids
  
  conds <- NULL
  if(is.null(key)){
    nsamp <- ncol(formatted.counts$raw) - 2
    smps <- names(formatted.counts$raw)
    smps <- smps[!smps %in% c("ID", "GeneSymbol")]
    key = tibble(Sample = smps, Group = "s")
    warning("No sample key given. Creating countDataSet using one condition for all samples.")
  }

  if(is.null(cds)){
    cds <- bic.get.deseq.cds(formatted.counts$raw, conds, min.count=min.count, libsizeQ=libsizeQ, 
                             percentile=percentile, fitType=fitType, method=method,
                             sharingMode=sharingMode)
  }

  log_debug("Getting DESeq scaled counts...")
  counts.scaled <- counts(cds, norm = T) %>%
                   as_tibble(rownames = "ID") %>%
                   left_join(formatted.counts$ids, by = "ID") %>%
                   select(ID, dplyr::matches("Gene"), key$Sample)
  log_debug("Done.\n")

  log_debug("Formatting scaled counts...")
  if(!is.null(key)){
    counts.scaled <- counts.scaled %>%
                     gather(names(.)[!grepl("^ID|^Gene", names(.))], key = "Sample", val = "Count") %>%
                     left_join(key, by = "Sample") %>%
                     unite("Sample_Group", "Sample", "Group", sep = "__") %>%
                     spread(Sample_Group, Count)  
  } 

  counts <- list()
  counts$scaled <- counts.scaled
  counts$raw <- formatted.counts$raw
  counts$ids <- idsAndGns
  return(counts) 
}


#' Use quantile normalization to normalize raw HTSeq counts 
#' 
#' Run \code{normalizeBetweenArrays()} function to normalize
#' raw HTSeq counts and return a list containing both raw
#' and normalized counts
#'
#' @param htseq.counts  matrix containing raw counts from HTSeq. One column per 
#'                      sample, one row per gene
#' @param key           matrix where column one contains sample IDs that matches 
#'                      those in counts file, and column two contains the condition 
#'                      groups to which they belong. Only those samples that are in 
#'                      the key matrix will be included in the normalization
#' @return list containing two matrices (one of raw counts and one of quantile 
#'         normalized counts) and a vector of gene symbols if one was included in 
#'         the raw counts matrix
#' @export 
bic.quantile.normalize.htseq.counts <- function(htseq.counts,key=NULL){
  ## read and reformat raw counts file
  ## use unique identifiers as rownames (e.g., ensembl IDs)
  ## but keep gene symbols (if given) for later use)
  log_debug("Reformatting raw counts...")
  HTSeq.dat <- bic.file2matrix(htseq.counts,header=T, sep="\t")
  HTSeq.dat <- bic.make.rownames(HTSeq.dat)

  gns <- NULL
  idsAndGns <- NULL

  if ("GeneSymbol" %in% colnames(HTSeq.dat)){
    gns <- HTSeq.dat[,1]
    idsAndGns <- as.matrix(gns)
    rownames(idsAndGns) <- rownames(HTSeq.dat)
    HTSeq.dat <- HTSeq.dat[,-1]
  }

  ## only keep columns for samples in key
  if(!is.null(key)){
    HTSeq.dat <- HTSeq.dat[,key[,1]]
  }

  samps <- colnames(HTSeq.dat)

  ## for kallisto counts, round floats to integers
  counts.dat <- round(bic.matrix2numeric(HTSeq.dat))

  log_debug("Done.\n")

  log_debug("Getting quantile-normalized counts...")
  counts.log.dat=log2(counts.dat+1)
  counts.log.norm.dat=normalizeBetweenArrays(counts.log.dat,method='quantile')
  dat=2^counts.log.norm.dat
  ## append gene symbols to matrix if there are any
  if (exists("gns")){
    dat=as.matrix(cbind(rownames(dat),gns,dat))
    colnames(dat)[1:2]=c("ID","GeneSymbol")
  } else {
    dat=as.matrix(cbind(rownames(dat),dat))
    colnames(dat)[1]="GeneSymbol"
  }
  log_debug("Done.\n")

  counts <- list()
  counts$normalized <- dat
  counts$raw <- counts.dat
  counts$gns <- gns
  return(counts)   
}


#' Get the gene sets that meet criteria set by cut offs
#'
#' @param gsaRes      results from \code{piano::runGSA()}
#' @param gsa.tab     BIC-ified GSA results
#' @param fc2keep     minimum absolute value fold change; Default: log2(1.5)
#' @param frac2keep   at least 1/frac2keep genes in set must have FC >= fc2keep
#' @param dir         direction ["Up"|"Dn"] of the results in gsaRes
#'                     
#' @return filtered gsa.tab
bic.get.gene.sets <- function(gsaRes,gsa.tab,fc2keep=log2(1.5),frac2keep=2,dir="Up"){

  tmp <- gsa.tab %>% mutate(GenesInGeneSet = "", GenesAndFC = "")

  geneSets2remove <- c()
  for(nn in 1:nrow(gsa.tab)){
    fc <- geneSetSummary(gsaRes, gsa.tab$Name[nn])$geneLevelStats

    o <- order(abs(fc),decreasing=T)
    fc <- fc[o]
    gn.names <- names(fc) 
    gn.names.str <- paste0(gn.names, collapse = ",") 
    gn.names.fc.str <- paste0(paste(gn.names,"=",fc), collapse = ";")
    tmp$GenesInGeneSet[nn] <- gn.names.str
    tmp$GenesAndFC[nn] <- gn.names.fc.str

    ## Add a condition that at least half of the genes have to have 
    ## FC >= log2(1.5) ! 

    if(dir == "Dn"){ jj <- length(which(fc <= -fc2keep)) }
    if(dir == "Up"){ jj <- length(which(fc >= fc2keep)) }
    if(jj <= length(fc)/frac2keep){ geneSets2remove=c(geneSets2remove,nn) }
  }
  if(length(geneSets2remove) > 0){ 
      tmp <- tmp %>% filter(!row_number() %in% geneSets2remove) 
  } 

  return(tmp)
}

#' Format and filter results from piano::runGSA()
#' 
#' @param gsaRes    result object from piano::runGSA()
#' @param max.p     maximum p-value to use as cutoff; Default: 0.1
#' @param fc2keep   minimum fold change
#' @param frac2keep at least 1/frac2keep genes in set must have FC >= fc2keep
#' @param fcQ       remove gene sets that do not meet cutoffs set 
#' 
#' @return list of two matrices: one with up-regulated gene sets and one with down-regulated
#'         gene sets
bic.process.gsa.res <- function(gsaRes,max.p=0.1,fc2keep=log2(1.5),frac2keep=2,fcQ=T){

  gsa.tab.dn <- gsa.tab.up <- NULL
  up.cols <- c("Name","Genes (tot)", "Stat (dist.dir)",
             "p adj (dist.dir.up)","Genes (up)")
  dn.cols <- c("Name","Genes (tot)", "Stat (dist.dir)",
             "p adj (dist.dir.dn)","Genes (down)")

  tmp <- GSAsummaryTable(gsaRes) %>% 
         as_tibble() %>%
         arrange(desc(abs(`Stat (dist.dir)`)))

  if(!is.null(tmp) && nrow(tmp) > 0){
    dn <- tmp %>%
          select_at(dn.cols) %>%
          filter(`p adj (dist.dir.dn)` <= max.p)

    if(nrow(dn) > 0){
       gsa.tab.dn <- dn
       if(fcQ){
           gsa.tab.dn <- bic.get.gene.sets(gsaRes, 
                                           gsa.tab = dn, 
                                           fc2keep = fc2keep, 
                                           frac2keep = frac2keep, 
                                           dir = "Dn")
       }
    }

    up <- tmp %>% 
          select_at(up.cols) %>%
          filter(`p adj (dist.dir.up)` <= max.p) 

    if(nrow(up) > 0){
       gsa.tab.up <- up
       if(fcQ){
           gsa.tab.up <- bic.get.gene.sets(gsaRes, 
                                           gsa.tab = up, 
                                           fc2keep = fc2keep, 
                                           frac2keep = frac2keep,
                                           dir = "Up")
       }
    }
  }

  tab.res <- list()
  tab.res$gsa.tab.dn <- gsa.tab.dn
  tab.res$gsa.tab.up <- gsa.tab.up
  return(tab.res)
}

#' Run gene set analysis using the package 'piano'
#'
#' Run gene set analysis for one pairwise comparison of sample groups. 
#' GSA is run using the R package 'piano'. For each gene set, means 
#' are calculated for the genes with positive fold changes and for
#' genes with negative fold changes. P-values and adjusted p-values 
#' are then calculated and gene sets with a significant number of
#' of up- or down-regulated genes are reported in the output. Gene Sets
#' included are from the MSigDB Collections:
#'
#' H:  hallmark gnee sets 
#' C2: curated gene sets
#' C3: motif gene sets
#' C5: GO gene sets
#' C6: oncogenic signatures
#' C7: immunologic signatures
#' 
#' Note: Mouse gene sets were downloaded from 
#'       http://bioinf.wehi.edu.au/software/MSigDB/
#'
#' @param species   currently only human and mouse are supported
#' @param deseq.res the matrix of ALL DESeq results (unfiltered)
#'                  in BIC format (ID,GeneSymbol,pval,P.adj,
#'                  log2[condB/condA],Mean_at_cond_condA,
#'                  Mean_at_cond_condB
#' @param min.gns   minimum number of genes a set must have in 
#'                  order to be included in analysis; Default: 5
#' @param max.gns   maximum number of genes a set must have in
#'                  order to be included in analysis; Default: 1000
#' @param max.p     maximum p-value to use as cut off; Default: 0.1
#' @param nPerm     the number of times the genes are randomized;
#'                  Default: 1e4
#' @param fc2keep   fold change cutoff (see frac2keep); Default: 
#'                  log2(1.5) 
#' @param frac2keep at least \code{1/frac2keep} of genes in a set
#'                  must have \code{fold change >= abs(fc2keep)};
#'                  Default: 4
#' @param fcQ       filter gene sets based on fc2keep and frac2keep       
#'                  Default: TRUE
#' @return a list of two matrices: one with gene sets found to be
#'         significantly UP-regulated and one with gene sets found
#'         to be significantly DOWN-regulated
#' @export
bic.run.gsa <- function(species, deseq.res, gmt.dir, min.gns=5, max.gns=1000,
                        max.p=0.1, nPerm=1e4, fcQ=T,
                        fc2keep=log2(1.5), frac2keep=4){

  spMap <- list(human = "human", hybrid = "human", 
                b37 = "human", hg19 = "human", hg18 = "human", 
                mouse = "mouse", mm9 = "mouse", mm10 = "mouse")

  gs.names <- list(human = c("h.all.v7.1.symbols.gmt", "c2.all.v7.1.symbols.gmt",
                             "c3.all.v7.1.symbols.gmt", "c5.all.v7.1.symbols.gmt",
                             "c6.all.v7.1.symbols.gmt", "c7.all.v7.1.symbols.gmt"),
                   mouse = c("Mm.h.all.v7.1.gmt", "Mm.c2.all.v7.1.gmt",
                             "Mm.c3.all.v7.1.gmt", "Mm.c5.all.v7.1.gmt",
                             "Mm.c6.all.v7.1.gmt", "Mm.c7.all.v7.1.gmt"))

  species <- spMap[[species]]
  gs.names <- gs.names[[species]]

  if(is.null(species) || length(species) == 0){
    log_error("    Only human and mouse data supported.") 
    log_error("        Can not run Gene Set Analysis\n")
    return(NULL)
  }

  all.up.res <- all.dn.res <- tibble()

  for (gs in gs.names){
    log_debug(paste0(gs, "\n"))
    gs.name <- gsub("\\.gmt", "", gs)
    gs.cat <- toupper(gsub("Mm\\.", "", gsub("\\.all.*", "", gs.name)))

    gsc <- loadGSC(file.path(gmt.dir, gs)) 

    fcCol <- names(deseq.res)[grep("log2", names(deseq.res))]
    fc <- deseq.res %>% 
          select(GeneSymbol, !!as.name(fcCol)) %>% 
          group_by(GeneSymbol) %>%
          summarize(FC = mean(!!as.name(fcCol))) 

    vec <- fc$FC
    names(vec) <- fc$GeneSymbol
   
    gsa.res <- runGSA(vec,
                      geneSetStat = "mean",
                      gsc = gsc,
                      gsSizeLim = c(min.gns,max.gns),
                      nPerm = nPerm)

    tab.res <- bic.process.gsa.res(gsa.res,
                                   max.p = max.p,
                                   fc2keep = fc2keep,
                                   frac2keep = frac2keep,
                                   fcQ = fcQ)

    if(!is.null(tab.res$gsa.tab.up)){
        tab.res$gsa.tab.up <- tab.res$gsa.tab.up %>% 
                              mutate(`Gene Set Category` = toupper(gs.name)) %>%
                              select(`Gene Set Category`, everything())
        all.up.res <- all.up.res %>% bind_rows(tab.res$gsa.tab.up)
    } 
    if(!is.null(tab.res$gsa.tab.dn)){
        tab.res$gsa.tab.dn <- tab.res$gsa.tab.dn %>% 
                              mutate(`Gene Set Category` = toupper(gs.name)) %>%
                              select(`Gene Set Category`, everything())
        all.dn.res <- all.dn.res %>% bind_rows(tab.res$gsa.tab.dn)
    }
  }

  noRes <- tibble(`Gene Set Category` = "",
                  Name = "",
                  `Genes (tot)` = "",
                  `Stat (dist.dir)` = "",
                  `p adj (dist.dir.up)` = "",
                  `Genes (up)` = "",
                  `GenesInGeneSet` = "",
                  GenesAndFC = "")

  if(nrow(all.up.res) == 0){
    all.up.res <- noRes
  }
  if(nrow(all.dn.res) == 0){
    all.dn.res <- noRes %>%
                  rename(`p adj (dist.dir.dn)` = `p adj (dist.dir.up)`,
                         `Genes (down)` = `Genes (up)`)               
  }
  list(up = all.up.res, dn = all.dn.res)

}
