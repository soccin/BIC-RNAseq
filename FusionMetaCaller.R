#!/usr/bin/env Rscript
library(data.table)
## library(FusionMetaCaller)

### ./FusionMetaCaller.R merged.txt merged.annotated.txt
###
### This script applies FusionMetaCaller ranking to fusions,
### after ignoring STAR results and
### grouping both by breakpoint and by gene pair.
### It adds columns and writes an annotated file.

FusionMetaCaller <- function (countMatrix, vote = 2, plot = F, trueFusion = NA){
  filterInd <- apply(countMatrix > 0, 1, sum) >= vote
  if(sum(filterInd) == 1){
    filterMatrix <- t(as.matrix(countMatrix[filterInd, ]))
    rownames(filterMatrix) <- rownames(countMatrix)[filterInd]
  } else {
    filterMatrix <- countMatrix[filterInd, ]
  }
  orderMatrix <- matrix(0, nrow(filterMatrix), ncol(filterMatrix))
  for (j in 1:ncol(filterMatrix)) {
    ind <- filterMatrix[, j] != 0
    orderMatrix[ind, j] <- order(filterMatrix[ind, j])
  }
  rankSum <- apply(orderMatrix, 1, sum)
  orderOfRankSum <- order(-rankSum)
  if(sum(filterInd) == 1){
    sortMatrix <- t(as.matrix(filterMatrix[orderOfRankSum, ]))
    rownames(sortMatrix) <- rownames(filterMatrix)[orderOfRankSum]
  } else {
    sortMatrix <- filterMatrix[orderOfRankSum, ]
  }
  if (plot == T) {
    co <- 2
    colLabel <- c("red", "blue", "green", "yellow", "purple",
                  "chocolate1", "hotpink", "orange", "darkblue", "darkgreen")
    ltyLabel <- c(rep(c(1, 5), each = 10))
    mark <- rep(0, times = nrow(countMatrix))
    mark[rownames(countMatrix) %in% trueFusion] <- 1
    pdf("precision_recall_plot.pdf")
    par(mar = c(5.1, 4.1, 4.1, 8.1), xpd = TRUE)
    plot(x = NULL, y = NULL, xlim = c(0, 1), ylim = c(0,
                                                      1), xlab = "Recall", ylab = "Precision", main = "Precison-Recall Plot",
         cex.lab = 1.5, cex.axis = 1.3)
    for (j in 1:ncol(countMatrix)) {
      sortCount <- sort(countMatrix[, j], decreasing = T)
      numFusion <- sum(sortCount > 0)
      sortMark <- mark[order(countMatrix[, j], decreasing = T)]
      sortMark <- sortMark[1:numFusion]
      CUTOFF <- seq(from = 1, to = length(sortMark), by = co)
      if (CUTOFF[length(CUTOFF)] != length(sortMark)) {
        CUTOFF <- c(CUTOFF, length(sortMark))
      }
      precision <- rep(0, times = length(CUTOFF))
      recall <- rep(0, times = length(CUTOFF))
      for (i in 1:length(CUTOFF)) {
        cutoff <- CUTOFF[i]
        TP <- sum(sortMark[1:cutoff])
        precision[i] <- TP/cutoff
        recall[i] <- TP/length(trueFusion)
      }
      collabel <- colLabel[j%%10]
      ltylabel <- ltyLabel[j%%20]
      lines(recall, precision, col = collabel, lwd = 2,
            lty = ltylabel)
    }
    sortMark <- rownames(sortMatrix) %in% trueFusion
    CUTOFF <- seq(from = 1, to = length(sortMark), by = co)
    if (CUTOFF[length(CUTOFF)] != length(sortMark)) {
      CUTOFF <- c(CUTOFF, length(sortMark))
    }
    precision <- rep(0, times = length(CUTOFF))
    recall <- rep(0, times = length(CUTOFF))
    for (i in 1:length(CUTOFF)) {
      cutoff <- CUTOFF[i]
      TP <- sum(sortMark[1:cutoff])
      precision[i] <- TP/cutoff
      recall[i] <- TP/length(trueFusion)
    }
    collabel <- "black"
    ltylabel <- 2
    lines(recall, precision, col = collabel, lwd = 2, lty = ltylabel)
    legend("topright", inset = c(-0.35, 0), legend = c(colnames(countMatrix),
                                                       "meta-caller"), lty = c(rep(ltyLabel, length = ncol(countMatrix)),
                                                                               2), lwd = 3, col = c(rep(colLabel, length = ncol(countMatrix)),
                                                                                                    "black"), bty = "n")
    dev.off()
  }
  return(list(sortMatrix = sortMatrix))
}

annotate_merged_file <- function(merged,
                                 TAG_name = "TAG_breakpoints",
                                 FusionMetaCaller_rank_name = "FusionMetaCaller_rank_breakpoints"){

  ### uses FusionMetaCaller to rank fusion calls
  ### adds a column to a data.table (merged) with name FusionMetaCaller_rank_name
  ### requires columns with names "TOTAL_COVERAGE" and TAG_name

  merged_dc <- dcast.data.table(merged[!is.na(get(TAG_name))],
                                paste(TAG_name, "~", "CALLER_ID"),
                                value.var = "TOTAL_COVERAGE",
                                fill = 0)
  merged_m <- as.matrix(merged_dc[, -1, with = F]);
  rownames(merged_m) <- merged_dc[[1]]
  rm(merged_dc)

  fusions <- rownames(FusionMetaCaller(merged_m)[["sortMatrix"]])

  merged[, FusionMetaCaller_rank_name := match(get(TAG_name), fusions), with = F]

  merged
}

write.maf <- function (...){
  write.table(..., quote = F, col.names = T, row.names = F,
              sep = "\t")
}

if(!interactive()){

  args <- commandArgs(TRUE)
  filename <- args[[1]]; args <- args[-1]
  ### /ifs/e63data/schultzlab/wangq/rna-seq-normal/fusion/bladder/fusion/Proj_GTEX_Bladder_merged_fusions_s_SRR1071717.txt
  output_filename <- args[[1]]; args <- args[-1]

  ### load file
  merged <- suppressWarnings(fread(filename))
  merged <- merged[CALLER_ID != "STAR"]
  merged <- merged[order(-TOTAL_COVERAGE)]

  merged[!duplicated(merged,
                     by = c("CALLER_ID",
                            "GENE1",
                            "GENE2")),
         TAG_genes := paste(GENE1,
                            GENE2,
                            sep = ":")]

  merged <- annotate_merged_file(merged,
                                 "TAG_genes",
                                 "FusionMetaCaller_rank_genes")


  merged[!duplicated(merged,
                     by = c("CALLER_ID",
                            "GENE1",
                            "GENE2",
                            "CHR1",
                            "BREAK_POINT1",
                            "CHR2",
                            "BREAK_POINT2")),
         TAG_breakpoints := paste(GENE1,
                                  GENE2,
                                  CHR1,
                                  BREAK_POINT1,
                                  CHR2,
                                  BREAK_POINT2,
                                  sep = ":")]

  merged <- annotate_merged_file(merged,
                                 "TAG_breakpoints",
                                 "FusionMetaCaller_rank_breakpoints")

  write.maf(merged, output_filename)
}
