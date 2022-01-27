bic.counts.to.gct <- function(norm.counts, idCol = "GeneSymbol"){
    print(1)
}

bic.de.to.rnk <- function(dat, statCol, file = NULL){

    idCol <- ifelse("GeneSymbol" %in% names(dat), "GeneSymbol", "ID")

    rnk   <- dat %>%
             select_at(c(idCol, statCol)) %>%
             filter(complete.cases(.)) %>%
             group_by_at(idCol) %>%
             summarize(stat = mean(!!as.name(statCol)))

    if(!is.null(file)){
        write.table(rnk, file = file, quote = F, sep="\t", col.names = F, row.names = F)    
    }
}

bic.write.cls.file <- function(key, file.name){

    smpNum <- nrow(key)
    grpNum <- length(unique(key$Group))
    cls <- rbind(c(smpNum, grpNum, 1, rep(NA, smpNum - 3)),
                 key %>% select(Group, Sample) %>% t)
    write.table(cls, file = file.name, row.names=F, sep="\t", col.names=F, quote=F)

}


bic.run.gsea.preranked <- function(de.dat, statCol, gsea.sh, reportLabel, outDir, 
                                   gs = NULL, gmt.dir = NULL,
                                   gs.min = 15, gs.max = 500){

    if(is.null(gs) && is.null(gmt.dir)){
        log_error("Must specify either single *.gmt file or directory of *.gmt files")
        stop()
    }

    allResDir <- file.path(outDir, "all_comparisons")
    dir.create(allResDir, recursive = T, mode = "0755")

    rnkFile <- file.path(outDir, paste0(reportLabel, ".rnk"))
    bic.de.to.rnk(de.dat, statCol, rnkFile)

    tf <- NULL
    if(!is.null(gmt.dir)){
        tf <- tempfile(pattern = "tmp", fileext = ".txt")
        write.table(dir(gmt.dir, full.names=T, pattern=".gmt"), tf, row.names=F, quote=F, col.names=F)
        gs.db <- paste("-gmx_list", tf)
    } else {
        gs.db <- paste("-gmx", gs)
    }

    cmd = paste("sh", gsea.sh, "GSEAPreranked", 
                "-rnk", rnkFile,
                gs.db,   
                "-rpt_label", reportLabel,
                "-out", allResDir)

    tryCatch({
        system(cmd)
        resDir <- dir(allResDir, full.names = T, pattern = paste0(reportLabel, "\\.GseaPreranked\\.*"))
        relDir <- getRelativePath(resDir, outDir)
        file.symlink(file.path(relDir, "index.html"), file.path(outDir, paste0(reportLabel, ".html")))
        unlink(rnkFile)
     }, error = function(e){
        str(e) 
     }, finally = function(){
        if(!is.null(tf)){ rm(tf) }
    })

}

