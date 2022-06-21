bic.de.to.rnk <- function(dat, statCol, file = NULL){

    idCol <- ifelse("GeneSymbol" %in% names(dat), "GeneSymbol", "ID")

    rnk   <- dat %>%
             select_at(c(idCol, statCol)) %>%
             filter(complete.cases(.)) %>%
             group_by_at(idCol) %>%
             summarize(stat = mean(!!as.name(statCol)))

    if(!is.null(file)){
        write.table(rnk, file = file, quote = F, sep="\t", col.names = F, row.names = F, na = '')    
    }
}

bic.write.gct.file <- function(counts, key, file.name){
    idCol <- ifelse("GeneSymbol" %in% names(counts), "GeneSymbol", "ID")
    clsKey <- key %>% 
              mutate(Sample = paste0(Sample, "__", Group))

    header <- rbind(c('#1.2', rep(NA, nrow(key) + 1)), c(nrow(counts), nrow(key), rep(NA, nrow(key))))

    dat <- counts %>%
           mutate(Description = "NA") %>%
           select(NAME = idCol, Description, clsKey$Sample) %>%
           mutate_if(is.numeric, list(~format(., scientific = FALSE, digits = 4))) %>%
           as.matrix
    gct <- rbind(header, colnames(dat), dat)
    write.table(gct, file = file.name, row.names = F, sep="\t", col.names = F, quote = F, na = '')
}


bic.write.cls.file <- function(key, file.name){

    smpNum <- nrow(key)
    grpNum <- length(unique(key$Group))
    cls <- rbind(c(smpNum, grpNum, 1, rep(NA, smpNum - 3)),
                 c("#", unique(key$Group), rep(NA, smpNum - (grpNum + 1))),
                 key %>% select(Group) %>% t)
    write.table(cls, file = file.name, row.names=F, sep="\t", col.names=F, quote=F, na = '')

}

bic.write.gmx.file <- function(gmt.files = NULL, gmt.dir = NULL){
    if((is.null(gmt.files) || length(gmt.files) == 0) && is.null(gmt.dir)){
        log_error("Must specify either a vector of *.gmt files or directory containing *.gmt files")
        stop()
    }
    if(!is.null(gmt.dir)){
        gmt.files <- dir(gmt.dir, full.names = T, pattern = '.gmt')
    }
    tf <- tempfile(pattern = "tmp", fileext = ".txt")
    write.table(gmt.files, tf, row.names = F, quote = F, col.names = F, na = '')
    return(tf)
}

bic.get.newest.file <- function(files){
    file.info(files) %>%
    filter(mtime == max(mtime)) %>%
    rownames()
}

bic.link.gsea.report <- function(resDir, reportLabel, linkLoc){

    gseaDir <- dir(resDir, full.names = T, pattern = paste0(reportLabel, "\\.Gsea\\.*"))
    if(length(gseaDir) > 1){
        log_warn("Found multiple results directories. Linking to latest.. ")
        gseaDir <- bic.get.newest.file(gseaDir)
    }
    relDir <- getRelativePath(gseaDir, linkLoc)
    htmlLink <- file.path(linkLoc, paste0(reportLabel, ".html"))
    if(file.exists(htmlLink)){
        log_warn("Removing existing link: ", htmlLink)
        unlink(htmlLink)
    }
    file.symlink(file.path(relDir, "index.html"), file.path(linkLoc, paste0(reportLabel, ".html")))   
}


bic.run.gsea <- function(ncounts, key, comp, gsea.sh, reportLabel, outDir,
                         gmt.dir = NULL, gmt.files = NULL,
                         gs.min = 15, gs.max = 500){

    allResDir <- file.path(outDir, "all_comparisons")
    dir.create(allResDir, recursive = T, mode = "0755", showWarnings = F)

    ## prep input
    clsFile <- file.path(outDir, paste0(reportLabel, ".cls"))
    bic.write.cls.file(key, clsFile)

    gctFile <- file.path(outDir, paste0(reportLabel, ".gct"))
    bic.write.gct.file(ncounts, key, gctFile)
print(4)
    tf <- bic.write.gmx.file(gmt.files = gmt.files, gmt.dir = gmt.dir)
    gs.db <- paste("-gmx_list", tf)
print(5)
    cmd = paste("sh", gsea.sh, "GSEA",
                "-res", gctFile, 
                "-cls", paste0(clsFile, "#", gsub("_vs_", "_versus_", comp)),
                gs.db,
                "-rpt_label", reportLabel,
                "-out", allResDir)
print(6)
    log_debug(cmd)

    tryCatch({
        system(cmd)
        bic.link.gsea.report(allResDir, reportLabel, outDir) ## returns true or false, or NULL if multiple directories found
     }, error = function(e){
        str(e)
        FALSE
     }, finally = function(){
        if(!is.null(tf)){ rm(tf) }
    })

}

bic.gsea.per.collection <- function(collections, ncounts, key, comp, gsea.sh, reportLabel, outDir, 
                                    gmt.files, gs.min = 15, gs.max = 500){
    errFound <- FALSE
    dn <-
     lapply(tolower(collections), function(coll){

        reportLabel <- paste0(reportLabel, "__", coll)
        gmt_file <- gmt.files[grep(paste0("(Mm\\.|^)(", coll, ")\\."), basename(gmt.files))]

        tryCatch({
            log_info("Running GSEA on comparison [", reportLabel, "]...")
            bic.run.gsea(ncounts, key, comp, gsea.sh, reportLabel, outDir,
                         gmt.files = gmt_file, gs.min = gs.min, gs.max = gs.max)
            log_info("Done!")
         }, error = function(e){
            log_error("An error occurred while running GSEA for comparison [", reportLabel, "].")
            errFound <<- TRUE
        })
    })
    if(errFound){
        stop("One or more errors thrown while running GSEA per collection.")
    }
}


bic.gsea.preranked.per.collection <- function(collections, dat, gmt.files, statCol, gsea.sh, reportLabel, outDir,
                                gs.min = 15, gs.max = 500){
    errFound <- FALSE
    dn <- 
     lapply(tolower(collections), function(coll){

        reportLabel <- paste0(reportLabel, "__", coll)
        gmt_file <- gmt.files[grep(paste0("(Mm\\.|^)(", coll, ")\\."), basename(gmt.files))]

        tryCatch({
            log_info("Running GSEAPreranked on comparison [", reportLabel, "]...")
            bic.run.gsea.preranked(dat, statCol, gsea.sh, reportLabel, outDir,
                                   gmt.files = gmt_file, gs.min = gs.min,
                                   gs.max = gs.max)
            log_info("Done!")
          }, error = function(e){
            log_error("An error occurred while running GSEAPreranked for comparison [", reportLabel, "].")
            errFound <<- TRUE
        })
    })
    if(errFound){
        stop("One or more errors thrown while running GSEAPreranked per collection.")
    }
}

bic.run.gsea.preranked <- function(de.dat, statCol, gsea.sh, reportLabel, outDir, 
                                   gmt.dir = NULL, gmt.files = NULL,
                                   gs.min = 15, gs.max = 500){

    allResDir <- file.path(outDir, "all_comparisons")
    dir.create(allResDir, recursive = T, mode = "0755")

    tf <- bic.write.gmx.file(gmt.files = gmt.files, gmt.dir = gmt.dir)
    gs.db <- paste("-gmx_list", tf)

    rnkFile <- file.path(outDir, paste0(reportLabel, ".rnk"))
    bic.de.to.rnk(de.dat, statCol, rnkFile)

    cmd = paste("sh", gsea.sh, "GSEAPreranked", 
                "-rnk", rnkFile,
                gs.db,   
                "-rpt_label", reportLabel,
                "-out", allResDir)

    tryCatch({
        system(cmd)
        bic.link.gsea.report(allResDir, reportLabel, outDir)
        unlink(rnkFile)
     }, error = function(e){
        str(e) 
     }, finally = function(){
        if(!is.null(tf)){ rm(tf) }
    })

}

