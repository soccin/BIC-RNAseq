#! /juno/opt/common/bic/R/R-4.0.2 

suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(rlang))
suppressPackageStartupMessages(library(lattice))
suppressPackageStartupMessages(library(locfit))
suppressPackageStartupMessages(library(DESeq))
suppressPackageStartupMessages(library(dendextend))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(piano))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(reshape))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(ggdendro))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(GSA))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(Rserve))
suppressPackageStartupMessages(library(RSclient))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(logger))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(showtext))
suppressPackageStartupMessages(library(matrixStats))

log_threshold(DEBUG)
font <- "Roboto"
font_add_google(font, font)
showtext_opts(dpi = 300)
showtext_auto(enable = TRUE)
#font <- "Roboto Condensed"
#font <- "Heebo"
#font <- "Source Sans Pro"
#font <- "Manrope"
#font <- "Arial OS"

source(file.path(bin,"bicrnaseqR/bic_analyze_counts.R"))
source(file.path(bin,"bicrnaseqR/bic_plots.R"))
source(file.path(bin,"bicrnaseqR/bic_util.R"))
source(file.path(bin,"bicrnaseqR/bic_deseq_report.R"))
source(file.path(bin,"bicrnaseqR/bic_gsea.R"))
