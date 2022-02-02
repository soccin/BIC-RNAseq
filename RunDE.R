#!/juno/opt/common/bic/R/R-4.0.2/bin/R

usage <- function(){

  cat("\nUsage: Rscript RunDE.R

  [REQUIRED]

    --counts.file     absolute path to htseq counts file


  [CONDITIONALLY REQUIRED]

    --comp.file         absolute path to file of sample group assignments and comparisons; when not 
                        provided, script will look for individual [pre]_sample_key*.txt and 
                        [pre]_sample_comparisons*.txt files 
    --sample.key        path to a single key file; to be used when $PWD contains multiple key files
                        matching pattern [pre]_sample_key{\\d}.txt but only one is to be included OR
                        when the key file name does not match that pattern 
    --comps.only        path to a single comparisons file; to be used under the same circumstances as
                        --sample.key (to specify a single file or one that does not match pattern)
    --species           required when running GSA; [hg19|human|mm9|mm10|mouse|hybrid] Note: only human 
                        or mouse currently supported; specific build does not matter, as long as 
                        it is clearly human or mouse (hybrid will be run as human)
    --gmt.dir           required when running GSA; directory containing GMT files for correct species
    --gsea.sh           default = NULL; path to GSEA shell script
    --cibersortR        required when running CIBERSORT; path to cibersort R script
    --cibersortSigFile  required when running CIBERSORT; txt file containing gene signature data required to run CIBERSORT
    --request           required when creating PDF report; project request file for generating PDF report
    --svnrev            required when creating PDF report; SVN revision number (to be included in PDF report)


  [OPTIONAL]

    --forceRerun     default = FALSE; when TRUE, script will NOT attempt to read previously run results
                     to save time and will run ALL analyses for ALL comparisons no matter what
    --Rlibs          default = NULL: path to R libraries 
    --pre            default = 'TEMP'; prefix for all output files
    --diff.exp       default = T; run differential expression analysis
    --GSA            default = T; run gene set analysis; if running GSA but not DESeq, BE SURE TO SET 
                     {diff.exp.dir} (see below) to point to existing DESeq results directory.
    --cibersort      default = T; run CIBERSORT; when TRUE, must also specify path to CIBERSORT Rscript file
                     AND path to CIBERSORT signatures file
    --heatmaps       default = T; generate a heatmap for each comparison showing relative expression of the 
                     top ~100 DE genes
    --counts.dir     default = ./counts_gene
    --clustering.dir default = ./clustering
    --cibersort.dir  default = ./cell_composition
    --diff.exp.dir   default = ./differentialExpression_gene; if running DESeq, this is where 
                     output will go; if NOT running DESeq, this is wherever DESeq results already 
                     exist. Must be correct in order for GSA to run properly.
    --gsa.dir        default = ./GSA
    --gsea.dir       default = ./GSEA
    --qc.dir         default = ./QC
    --no.replicates  default = F: T|F, automatically sets fitType, method, and sharingMode to 
                     accomodate comparison of single samples
    --min.abs.fc     default = 2: Fold change cutoff
    --min.count      default = 15: minimum total reads required for a gene
                     to be included in analysis
    --fitType        default = 'parametric': ['parametric'|'local', See DESeq docs]
    --orderPvalQ     default = T: T|F, sort final results by adjusted Pvalue
    --method         default = 'per-condition': ['pooled'|'per-condition'|'blind']; See DESeq docs.
    --sharingMode    default = 'maximum': ['maximum'|'fit-only'|'gene-est-only']; See DESeq docs.
    --zeroaddQ       default = T: [T|F], include genes with zero counts in one condition even if 
                     their pValues are insignificant
    --libsizeQ       default = F: [T|F], normalize data using library size
    \n\n")
}


################################################################################
####   get user input
################################################################################

## unset LD_LIBRARY_PATH in order for GSEA to work
Sys.setenv("LD_LIBRARY_PATH" = "")

## need to read args with base commandArgs first, just to get Rlibs dir, then
## need to load R.utils in order to be able to read command args with names and
## to load package 'logger' to start logging immediately

tmpargs <- commandArgs(trailingOnly = TRUE)
if(!any(tmpargs == '--bin')){ stop("No 'bin' directory provided.") }
bin <- tmpargs[which(tmpargs == '--bin') + 1]
Rlibs <- tmpargs[grep("Rlibs$", tmpargs) + 1]

if(!is.null(Rlibs)){
    print(paste0("Loading libraries from ", Rlibs)) 
    .libPaths(c(Rlibs))
}

source(file.path(bin, 'bicrnaseqR/defaults.R'))
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(logger))
log_threshold(DEBUG)

log_info(c("\n++++++++++++++++ BIC RNA-Seq Counts Analysis ++++++++++++++++\n\n"))

args <- commandArgs(trailingOnly = TRUE, asValue = TRUE, adhoc = TRUE)
tmp <- lapply(1:length(args), function(x){
         assign(names(args)[x], args[[x]], inherits = T) 
       })

log_info("Loading bicrnaseq library...")
source(file.path(bin, "bicrnaseqR/source_bicrnaseq.R"))
log_info("Done.\n")

report    <- bic.report.runnable(report, request, svnrev)
GSA       <- bic.gsa.runnable(GSA, species)
GSEA      <- bic.gsea.runnable(GSEA, gsea.sh, gmt.dir)
cibersort <- bic.cibersort.runnable(cibersort, cibersortR, cibersortSigFile)

if(!diff.exp){ diff.exp.dir <- NULL }
if(!GSA){ gsa.dir <- NULL }
if(!GSEA){ gsea.dir <- NULL }
if(!cibersort){ cibersort.dir <- NULL }

bic.log.params(pre, forceRerun, cibersort, diff.exp, GSA, GSEA, report, counts.file, cibersortR,
             cibersortSigFile, gmt.dir, species, gsea.sh, request, svnrev,
             counts.dir, clustering.dir, cibersort.dir, 
             all.gene.dir, diff.exp.dir, gsa.dir, gsea.dir, qc.dir, 
             max.p, lfc, min.abs.fc, min.count, percentile, fitType,
             orderPvalQ, method, sharingMode, zeroaddQ, libsizeQ, no.replicates, heatmaps,
             comps.only, comp.file, sample.key)

##
## set up sample key, comparisons and output directories
##

## For now, need to read one or more pairs of key/comp files at a time. 
## In the future, possibly provide a single key file that will contain a table where each
## column contains exactly two unique values that represent groups to be compared
compSetNums <- 1
keysAndComps <- bic.get.keys.and.comparisons(pre, path = ".", comp.file = comp.file, 
                                             comps.only = comps.only, sample.key = sample.key) 

if(!is.null(keysAndComps)){
    compSetNums  <- lapply(keysAndComps$keys, function(x) !is.null(x)) %>% unlist() %>% which()
    multi.comps  <- length(keysAndComps$keys) > 1 || !1 %in% compSetNums 
}

out.dirs <- bic.setup.directories(counts.dir, clustering.dir,
                                  all.gene.dir  = all.gene.dir,
                                  cibersort.dir = cibersort.dir, 
                                  diff.exp.dir  = diff.exp.dir, 
                                  gsa.dir       = gsa.dir,
                                  gsea.dir      = gsea.dir,
                                  qc.dir        = qc.dir,
                                  multi.comps   = multi.comps, 
                                  comp.set.nums = compSetNums) 

##
## Run all comparison sets 
##  
## A set of comparisons is defined as the group of pair-wise sample group comparisons
## where each sample is involved in ONE comparison (i.e, all comparisons that can be made
## using one sample key and one comparison file as input)
##

nc <- min(detectCores()/2, 6, length(compSetNums))
cl <- makeCluster(nc, type = "FORK", outfile = "")
clusterExport(cl, 
              c('pre', 'counts.file', 'keysAndComps', 'out.dirs', 
                'diff.exp', 'GSA', 'GSEA', 'cibersort', 'report',
                'min.count', 'min.abs.fc', 'libsizeQ', 'zeroaddQ', 'percentile',
                'fitType', 'method', 'sharingMode', 'orderPvalQ',
                'heatmaps', 'species', 'request', 'svnrev', 
                'gsea.sh', 'gmt.dir', 'cibersortR', 'cibersortSigFile'), 
              envir = environment())

### each process will return exit code either 0 or 1
res <-
parLapply(cl, compSetNums, function(compSet){

    ec <<- 0
    if(diff.exp){
        log_info(paste0("Analyzing comparison set [", compSet, "]\n"))
        log <- file.path(dirname(out.dirs[[compSet]]$DEdir), paste0("comparisons_", compSet, ".log"))
    } else {
        log <- file.path(dirname(out.dirs[[compSet]]$clusterDir), paste0("clustering.log"))
    }

    log_debug("Output directories for comparison set [", compSet, "]")
    lapply(names(out.dirs[[compSet]]), function(x){
        log_debug("  [", x, "] ", out.dirs[[compSet]][[x]])
    })

    key <- conds <- comps <- NULL
    if(!is.null(keysAndComps)){
        key   <- keysAndComps$keys[[compSet]]
        comps <- keysAndComps$comparisons[[compSet]]
        conds <- key$Group
        sColors <- bic.get.colors(length(unique(key$Group)))
        names(sColors) <- unique(key$Group)
    }

    log_info(paste0("Output to be captured in log file [", log, "]\n"))
    
    if(!interactive()){
        lf <- file(log)
        sink(lf, type = "output")
        sink(lf, type = "message") 
    }

    results_updated <- FALSE
    resFile <- file.path(dirname(out.dirs[[compSet]]$countsDir), "all_results.rda")
    all_results <- list()
    idMap <- NULL
    if(!forceRerun && file.exists(resFile)){
        log_info("Loading results from previously-run DE analysis...")
        all_results <- readRDS(resFile)
        idMap <- all_results$norm.counts %>% select(!dplyr::matches("^s_"))
        cds <- all_results$cds
        norm.counts <- all_results$norm.counts
    }


    if(is.null(all_results$cds)){    
        ## format counts
        formatted.counts <- bic.format.htseq.counts(counts.file, key=key) ## this function can handle null key
        idMap <- formatted.counts$ids

        ## if no conditions given, create a vector of one mock condition, 
        ## needed for normalization
        if(ncol(key) < 2 || is.null(conds)){
            conds <- rep("s", ncol(formatted.counts$raw) - 2)
        }

        ## get DESeq countDataSet
        cds <- bic.get.deseq.cds(formatted.counts$raw,
                                 conds,
                                 min.count = min.count,
                                 libsizeQ = libsizeQ,
                                 percentile = percentile,
                                 fitType = fitType,
                                 method = method,
                                 sharingMode = sharingMode)

        all_results$cds <- cds
        if(is.null(cds)){
            log_info("\nZERO RESULTS. Could not cluster data, writing empty DE results files.\n")
            return 
         }
    
        log_info("Writing normalized counts to file...")
        norm.counts <- bic.deseq.normalize.htseq.counts(formatted.counts = formatted.counts, 
                                                        cds = cds, 
                                                        key = key) 
        all_results$norm.counts <- norm.counts$scaled
        idMap <- all_results$norm.counts %>% select(!dplyr::matches("^s_"))

        file.name <- file.path(out.dirs[[compSet]]$countsDir, paste0(pre, "_counts_scaled_DESeq.xlsx"))
        bic.write.normalized.counts.file(norm.counts$scaled, file.name)
        log_info("Done.\n")

        if(is.null(norm.counts$scaled)){
            log_error("Can not find normalized counts matrix. Can not cluster samples or run DESeq.\n")
            return
        }
        log_info("Clustering all samples...")
        bic.standard.clustering(cds, norm.counts$scaled, 
                                out.dirs[[compSet]]$clusterDir,
                                idMap = idMap, 
                                pre = pre,
                                annColors = sColors)
        log_info("Done.\n")

        log_info("Plotting high expression heatmap...")
        file.name <- file.path(out.dirs[[compSet]]$clusterDir, paste0(pre, "_top_50_genes.pdf"))
        tmp <- bic.high.expression.heatmap(cds, transform = TRUE,
                                           num.gns = 50, idMap = idMap,
                                           key = key, clrs = sColors, 
                                           file = file.name) 

        log_debug("Plotting sample to sample distances...")
        file.name <- file.path(out.dirs[[compSet]]$clusterDir, paste0(pre, "_sample_to_sample.pdf"))
        tmp <- bic.sample.to.sample.distances(cds, 
                                              key = key, 
                                              clrs = sColors, file = file.name) 

        cbsDat <- NULL
        if(cibersort){
            log_info("Running CIBERSORT...")
            file.name <- file.path(out.dirs[[compSet]]$cellCompDir, "cibersort.txt")
            tryCatch({
                cbsDat <- bic.run.cibersort(norm.counts$scaled,
                                            file.name,
                                            cibersortR,
                                            sigFile = cibersortSigFile)
                all_results$cibersort.res <- cbsDat
              }, error = function(e){
                log_error(e)
                log_error("CIBERSORT failed.")
                ec <<- 1
                return(NULL)
            })
            log_info("Done.\n")
        }

        saveRDS(all_results, resFile)
        results_updated <- TRUE
    }    

    if(diff.exp){

        for (comp in comps){
            compName <- paste(unlist(strsplit(comp, " - "))[2:1], collapse = "_vs_")
            log_info("Working on comparison [", compName, "].")

            if(!compName %in% names(all_results)){
                all_results[[compName]] <- list()
            }

            ## DESeq
            if(is.null(all_results[[compName]]$DE) || forceRerun){
                all_results[[compName]]$DE <- 
                     tryCatch({
                         bic.complete.de.analysis(cds,
                                                  norm.counts, 
                                                  conds, 
                                                  out.dirs[[compSet]]$allGeneDir,
                                                  out.dirs[[compSet]]$DEdir,
                                                  diff.exp.fig.dir = out.dirs[[compSet]]$DEfigDir, 
                                                  qc.dir      = out.dirs[[compSet]]$QCdir,
                                                  condA       = gsub(".*_vs_", "", compName), 
                                                  condB       = gsub("_vs_.*", "", compName),
                                                  max.p       = max.p,
                                                  min.abs.fc  = min.abs.fc,
                                                  min.count   = min.count,
                                                  zeroaddQ    = zeroaddQ,
                                                  orderPvalQ  = orderPvalQ,
                                                  geneSymbols = idMap, 
                                                  heatmaps    = heatmaps,
                                                  annClrs     = sColors,
                                                  key = key %>% filter(Group %in% c(condA, condB)))
                       }, error = function(e){
                           print(str(e))
                           log_error("\nERROR running DESeq\n")
                           ec <<- 1
                           return(NULL) ## need this or else ec returned 
                      })
                if(is.null(all_results[[compName]]$DE)){ return }
                results_updated <- TRUE
            } else {
                log_warn("Found existing DE results for ", compName, ". Will NOT rerun DESeq for this comparison.")
            }

            ## Gene Set Enrichment Analysis
            if(GSEA){
                if(is.null(all_results[[compName]]$GSEA) || forceRerun){
                    log_info("Running GSEA...\n")
                    all_results[[compName]]$GSEA <- 
                        tryCatch({
                            fcCol <- paste0("log2[", gsub("_vs_", "/", compName), "]")
                            bic.run.gsea.preranked(all_results[[compName]]$DE$all.res,
                                                   fcCol,
                                                   gsea.sh,
                                                   compName,
                                                   out.dirs[[compSet]]$GSEAdir,
                                                   gmt.dir = gmt.dir)
                            ## GSEA doesn't return a value, so just save a message pointing to the results
                            paste0("GSEA output: ", out.dirs[[compSet]]$GSEAdir)
                        }, error = function(e){
                            print(str(e))
                            log_error("\nERROR running GSEA\n")
                            ec <<- 1 
                            return(NULL)
                        })
                    results_updated <- TRUE
                } else {
                    log_warn("Found existing GSEA results for ", compName, ". Will NOT rerun GSEA for this comparison.")
                }
            } else {
                log_info("GSEA is turned OFF.")
            }
  
            ## Gene Set Analysis (R package 'piano')   TO BE REMOVED?
            if(GSA){
                if(is.null(all_results[[compName]]$GSA) || forceRerun){
                    log_info("Running GSA...\n")
                    fcCol <- paste0("log2[", gsub("_vs_", "/", compName), "]")
                    allRes <- all_results[[compName]]$DE$all.res %>% filter(!is.na(!!as.name(fcCol)))
                    all_results[[compName]]$GSA <- 
                        tryCatch({
                            bic.complete.gsa.analysis(allRes,
                                                      compName, 
                                                      bic.gsa.species(species),                          
                                                      gmt.dir, 
                                                      out.dirs[[compSet]]$GSAdir)
                          }, error = function(e){
                            print(str(e))
                            log_error("\nError running GSA")
                            ec <<- 1
                            return(NULL)
                       })
                    log_info("Done.\n")
                    results_updated <- TRUE
                } else {
                    log_warn("Found existing GSA results for ", compName, ". Will NOT rerun.")
                }
            } else {
                log_info("GSA is turned OFF.")
            }

            if(results_updated){
                log_info("Saving updated results...")
                ## save after every comparison is complete
                saveRDS(all_results, resFile) 
            } else {
                log_info("No new results to save for [", compName, "].")
            }
        }

        if(report){
            log_info("Making DESeq PDF report...")
            tryCatch({
                fn <- ifelse(length(compSetNums) > 1, 
                             paste0(pre, "_DE_Report_comparisons",compSet,".pdf"), 
                             paste0(pre, "_DE_Report.pdf"))
                fileName <- file.path(out.dirs[[compSet]]$DEdir, fn) 

                rpt <- bic.deseq.report(all_results, request, svnrev, fileName, pre = pre, 
                                        maxPlots = 9, idMap = idMap, key = key,
                                        clusterDir = out.dirs[[compSet]]$clusterDir, 
                                        deFigDir = out.dirs[[compSet]]$DEfigDir,
                                        qcDir = out.dirs[[compSet]]$QCdir,
                                        annColors = sColors)
                log_info("Done.\n")
              }, error = function(e){
                print(str(e))
                log_error("ERROR making DESeq PDF report.")
                ec <<- 1 
                return(NULL)
            })
        }

    } else {
        log_info("Differential expression is OFF. No DE and no GSA analysis to run.")
    }

    sink(type = "output")
    sink(type = "message")

    if(ec == 1){
        log_error(paste0("One or more errors were thrown during analysis of comparison set [", compSet, 
                             "]. See comparisons_", compSet, ".log for details."))
        return(1)
    } else {
        log_info(paste0("All analyses in comparison set [",compSet,"] completed successfully."))
        return(0)
    }

}) ## end for compSet in allcomp sets

stopCluster(cl)

if(any(unlist(res) == 1)){
    q(save = "no", status = 1, runLast = FALSE)
}
