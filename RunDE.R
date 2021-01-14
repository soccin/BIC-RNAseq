#!/opt/common/CentOS_6/R/R-4.0.0/bin/R

usage <- function(){

  cat("\nUsage: Rscript RunDE.R

  [REQUIRED]

    --counts.file     absolute path to htseq counts file


  [CONDITIONALLY REQUIRED]

    --comp.file       absolute path to file of sample group assignments and comparisons; when not 
                      provided, script will look for individual [pre]_sample_key*.txt and 
                      [pre]_sample_comparisons*.txt files 
    --species         required when running GSA; [hg19|human|mm9|mm10|mouse|hybrid] Note: only human 
                      or mouse currently supported; specific build does not matter, as long as 
                      it is clearly human or mouse (hybrid will be run as human)
    --gmt.dir         required when running GSA; directory containing GMT files for correct species


  [OPTIONAL]

    --Rlibs          default = NULL: path to R libraries 
    --pre            default = 'TEMP'; prefix for all output files
    --diff.exp       default = T; run differential expression analysis
    --GSA            default = T; run gene set analysis; if running GSA but not DESeq, BE SURE TO SET 
                     {diff.exp.dir} (see below) to point to existing DESeq results directory.
    --heatmaps       default = T; generate a heatmap for each comparison showing relative expression of the 
                     top ~100 DE genes
    --counts.dir     default = ./counts_gene
    --clustering.dir default = ./clustering
    --diff.exp.dir   default = ./differentialExpression_gene; if running DESeq, this is where 
                     output will go; if NOT running DESeq, this is wherever DESeq results already 
                     exist. Must be correct in order for GSA to run properly.
    --gsa.dir        default = ./GSA
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
    --java           default = NULL: path to java for generating PDF report
    --javacp         default = NULL: java class path for generating PDF report
    --pdfjar         default = NULL: jar file for generating PDF report
    --request        default = NULL: project request file for generating PDF report
    \n\n")
}


################################################################################
####   Set defaults
################################################################################
pd                <- getwd()
counts.dir        <- file.path(pd,"counts_gene")
clustering.dir    <- file.path(pd,"clustering")
diff.exp.dir      <- file.path(pd,"differentialExpression_gene")
all.gene.dir      <- file.path(pd,"all_gene")
gsa.dir           <- file.path(pd,"GSA")
gmt.dir           <- NULL
max.p             <- 0.05
lfc               <- 1   
min.abs.fc        <- 2   
min.count         <- 15
percentile        <- "100%"
fitType           <- "parametric"
orderPvalQ        <- T
method            <- "per-condition"
sharingMode       <- "maximum"
zeroaddQ          <- T
libsizeQ          <- F
no.replicates     <- F
key               <- NULL
conds             <- NULL
GSA               <- T
diff.exp          <- T
heatmaps          <- T
pre               <- "TEMP"

key.file          <- NULL
comp.file         <- NULL
counts.file       <- NULL

norm.counts       <- NULL
cds               <- NULL
mds.labels        <- FALSE
Rlibs             <- NULL

java              <- NULL
javacp            <- NULL
pdfjar            <- NULL
request           <- NULL
svnrev            <- NULL 
report            <- TRUE
################################################################################
####   get user input
################################################################################

tmpargs <- commandArgs(trailingOnly = TRUE)

if(any(grepl("Rlibs", tmpargs))){
    Rlibs = tmpargs[grep("Rlibs$", tmpargs) + 1]
    print(paste0("Loading libraries from ", Rlibs)) 
    .libPaths(Rlibs)
    suppressPackageStartupMessages(library(R.utils))
    suppressPackageStartupMessages(library(logger))
    log_threshold(DEBUG)
}

log_info(c("\n++++++++++++++++ BIC RNA-Seq Counts Analysis ++++++++++++++++\n\n"))

## after loading R.utils, get command args again in list format
args <- commandArgs(trailingOnly = TRUE, asValue = TRUE, adhoc = TRUE)
for(i in 1:length(args)){
    assign(names(args)[i], args[[i]]) 
}

log_debug(paste0("LD_LIBRARY_PATH  ", Sys.getenv("LD_LIBRARY_PATH")))
log_debug(paste0("PATH             ", Sys.getenv("PATH")))

log_debug(paste0("pre: ", pre))
log_debug(paste0("counts.file: ", counts.file))
log_debug(paste0("counts.dir: ", counts.dir))
log_debug(paste0("clustering.dir: ", clustering.dir))
log_debug(paste0("all.gene.dir: ", all.gene.dir))
log_debug(paste0("diff.exp.dir: ", diff.exp.dir))
log_debug(paste0("gsa.dir: ", gsa.dir))
log_debug(paste0("gmt.dir: ", gmt.dir))
log_debug(paste0("max.p: ", max.p))
log_debug(paste0("lfc: ", lfc))
log_debug(paste0("min.abs.fc: ", min.abs.fc))
log_debug(paste0("min.count: ", min.count))
log_debug(paste0("percentile: ", percentile))
log_debug(paste0("fitType: ", fitType))
log_debug(paste0("orderPvalQ: ", orderPvalQ))
log_debug(paste0("method: ", method))
log_debug(paste0("sharingMode: ", sharingMode))
log_debug(paste0("zeroaddQ: ", zeroaddQ))
log_debug(paste0("libsizeQ: ", libsizeQ))
log_debug(paste0("no.replicates: ", no.replicates))
log_debug(paste0("key: ", key))
log_debug(paste0("conds: ", conds))
log_debug(paste0("GSA: ", GSA))
log_debug(paste0("diff.exp: ", diff.exp))
log_debug(paste0("heatmaps: ", heatmaps))
log_debug(paste0("key.file: ", key.file))
log_debug(paste0("comp.file: ", comp.file))
log_debug(paste0("java: ", java))
log_debug(paste0("javacp: ", javacp))
log_debug(paste0("pdfjar: ", pdfjar))
log_debug(paste0("request: ", request))
log_debug(paste0("svnrev: ", svnrev))
log_debug(paste0("report: ", report))

if(any(is.null(java), is.null(javacp), is.null(pdfjar), is.null(request), is.null(svnrev))){
    log_warn("Missing parameter(s) needed to make PDF report is here. Turning 'report' OFF.\n")
    log_warn("Parameters required for creating report:\n")
    log_warn("    java, javacp, pdfjar, request, svnrev\n")
    report <- FALSE
}

log_info("Loading bicrnaseq library...")
source(file.path(bin, "bicrnaseqR/source_bicrnaseq.R"))
log_info("Done.\n")

## do some validation
if((GSA && (!is.character(species) || is.null(species) || 
             (is.character(species) && 
                !species %in% c("human", "mouse", "hg19", "b37", "mm10", "hybrid")))) || GSA == FALSE){
    log_warn("Invalid or missing species (GSA only supported for human and mouse). Turning off GSA.\n")
    GSA <- FALSE
    gsa.dir <- NULL
}
if(!diff.exp){ diff.exp.dir = NULL }

## For now, automatically search working directory for 'sample_key*.txt' and 'comparisons*.txt';
## In the future, possibly provide a single key file that will contain a table where each
## column contains exactly two unique values that represent groups to be compared
keysAndComps <- bic.get.keys.and.comparisons(pre, path = ".") 
compSetNums  <- lapply(keysAndComps$keys, function(x) !is.null(x)) %>% unlist() %>% which()
multi.comps  <- length(keysAndComps$keys) > 1 || !1 %in% compSetNums 

out.dirs <- bic.setup.directories(counts.dir, clustering.dir,
                                  all.gene.dir  = all.gene.dir, 
                                  diff.exp.dir  = diff.exp.dir, 
                                  gsa.dir       = gsa.dir, 
                                  multi.comps   = multi.comps, 
                                  comp.set.nums = compSetNums) 

nc <- min(detectCores()/2, 6, length(keysAndComps$keys))
cl <- makeCluster(nc, type = "FORK", outfile = "")
clusterExport(cl, 
              c('pre', 'keysAndComps', 'out.dirs', 'diff.exp', 'GSA', 'counts.file', 
                'min.count', 'min.abs.fc', 'libsizeQ', 'zeroaddQ', 'percentile',
                'fitType', 'mds.labels', 'method', 'sharingMode', 'orderPvalQ',
                'heatmaps', 'report', 
                'java', 'javacp', 'pdfjar', 'request', 'svnrev'), 
              envir = environment())

errorsThrown <<- FALSE

parLapply(cl, compSetNums, function(compSet){

    log_info(paste0("Analyzing comparison set [", compSet, "]\n"))
    log <- file.path(dirname(out.dirs[[compSet]]$DEdir), paste0("comparisons_", compSet, ".log"))
    log_info(paste0("Output to be captured in log file [", log, "]\n"))
    lf <- file(log)
    sink(lf, type = "output")
    sink(lf, type = "message") 

    key   <- keysAndComps$keys[[compSet]]
    comps <- keysAndComps$comparisons[[compSet]]
    conds <- key$Group

    ## format counts
    formatted.counts <- bic.format.htseq.counts(counts.file, key=key)

    ## if no conditions given, create a vector of one mock condition, 
    ## needed for normalization
    if(ncol(key) < 2 || is.null(conds)){
        conds <- rep("s", ncol(formatted.counts$raw))
        mds.labels <- TRUE
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

    if(is.null(cds)){
        log_info("\nZERO RESULTS\n")
        bic.write.all.empty.results(comps, out.dirs[[compSet]]$DEdir)
        next
    }

    log_info("Running clustering and QC...")
    bic.deseq.qc(cds, out.dirs[[compSet]]$clusterDir)
    log_info("Done.\n")

    log_info("Writing normalized counts to file...")
    norm.counts <- bic.deseq.normalize.htseq.counts(formatted.counts = formatted.counts, 
                                                    cds = cds, 
                                                    key = key) 
    file.name <- file.path(out.dirs[[compSet]]$countsDir, paste0(pre, "_counts_scaled_DESeq.xlsx"))
    bic.write.normalized.counts.file(norm.counts$scaled, file.name)
    log_info("Done.\n")

    if(!is.null(norm.counts$scaled)){
        log_info("Clustering all samples...")
        bic.standard.clustering(norm.counts$scaled, 
                                out.dirs[[compSet]]$clusterDir, 
                                conds, 
                                mds.labels = mds.labels, 
                                pre = pre)
        log_info("Done.\n")
    } else {
       log_info("ERROR: Can not find normalized counts matrix. Can not cluster samples or run DESeq.\n")
        next
    }

    
    if(diff.exp == T){

        log_info("Running differential expression analysis...\n")
        all_results <- list(norm.counts = norm.counts$scaled,
                            cds = cds)

        for (comp in comps){
            compName <- paste(unlist(strsplit(comp, " - "))[2:1], collapse = "_vs_")

            all_results[[compName]] <- list()
            if(file.exists(file.path(dirname(out.dirs[[compSet]]$countsDir), "all_results.rda"))){
                all_results <- readRDS(file.path(dirname(out.dirs[[compSet]]$countsDir), "all_results.rda"))
            }

            if(is.null(all_results[[compName]]$DE)){
                all_results[[compName]]$DE <- 
                     tryCatch({
                         bic.complete.de.analysis(cds, 
                                                  norm.counts$scaled, 
                                                  conds, 
                                                  out.dirs[[compSet]]$allGeneDir,
                                                  out.dirs[[compSet]]$DEdir,
                                                  out.dirs[[compSet]]$DEfigDir,
                                                  condA       = gsub(".*_vs_", "", compName), 
                                                  condB       = gsub("_vs_.*", "", compName),
                                                  max.p       = max.p,
                                                  min.abs.fc  = min.abs.fc,
                                                  min.count   = min.count,
                                                  zeroaddQ    = zeroaddQ,
                                                  orderPvalQ  = orderPvalQ, 
                                                  geneSymbols = norm.counts$ids,
                                                  heatmaps = heatmaps)
                   }, error = function(e){
                       print(e)
                       log_error("\nERROR running DESeq\n")
                       errorsThrown <- TRUE
                  })
            }

            if(is.null(all_results[[compName]]$DE)){ next }

            if(GSA){
                if(is.null(all_results[[compName]]$GSA)){
                    log_info("Running GSA...\n")
                    fcCol <- paste0("log2[", gsub("_vs_", "/", compName), "]")
                    allRes <- all_results[[compName]]$DE$all.res %>% filter(!is.na(!!as.name(fcCol)))
                    all_results[[compName]]$GSA <- 
                        tryCatch({
                            bic.complete.gsa.analysis(allRes,
                                                      compName, 
                                                      species,                          
                                                      gmt.dir, 
                                                      out.dirs[[compSet]]$GSAdir)
                          }, error = function(e){
                            print(e)
                            log_error("\nError running GSA")
                            errorsThrown <- TRUE
                       })
                    log_info("Done.\n")
                } else {
                    log_info("GSA is turned OFF.")
                }
            }

            ## save after every comparison is complete
            saveRDS(all_results, file.path(dirname(out.dirs[[compSet]]$countsDir), "all_results.rda"))

        }

        if(report){
            log_info("Making DESeq PDF report...")
            tryCatch({
                system(paste(java, "-cp", javacp, "-jar", pdfjar,
                             "-rf", request,
                             "-v", svnrev,
                             "-cd", out.dirs[[compSet]]$clusterDir,
                             "-df", out.dirs[[compSet]]$DEfigDir,
                             "-o", out.dirs[[compSet]]$DEdir))
                log_info("Done.\n")
              }, error = function(e){
                print(e)
                log_error("ERROR making DESeq PDF report.")
                errorsThrown <- TRUE
            })
        }

    } else {
        log_info("Differential expression is OFF. No DE and no GSA analysis to run.")
    }
    sink(type = "output")
    sink(type = "message")
}) ## end for compSet in allcomp sets

stopCluster(cl)

if(errorsThrown){
    stop("One or more errors were thrown during analysis. See comparisons_*.log for details")
} else {
    log_info("All analyses complete.")
}
