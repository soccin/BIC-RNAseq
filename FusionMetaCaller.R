#!/usr/bin/env Rscript
library(data.table)
library(FusionMetaCaller)

### ./FusionMetaCaller.R merged.txt merged.annotated.txt
###
### This script applies FusionMetaCaller ranking to fusions,
### after ignoring STAR results and
### grouping both by breakpoint and by gene pair.
### It adds columns and writes an annotated file.

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

  fusions <- rownames(FusionMetaCaller::FusionMetaCaller(merged_m)[["sortMatrix"]])

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
