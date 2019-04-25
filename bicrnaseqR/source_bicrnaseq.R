#! /opt/common/CentOS_6/R/R-3.2.0/bin/R

suppressPackageStartupMessages(library(rlang))
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
suppressPackageStartupMessages(library(scales))


source(file.path(bin,"bicrnaseqR/bic_analyze_counts.R"))
source(file.path(bin,"bicrnaseqR/bic_plots.R"))
source(file.path(bin,"bicrnaseqR/bic_util.R"))
