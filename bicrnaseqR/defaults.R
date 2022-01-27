pd                <- getwd()
forceRerun        <- FALSE
counts.dir        <- file.path(pd, "counts_gene")
clustering.dir    <- file.path(pd, "clustering")
diff.exp.dir      <- file.path(pd, "differentialExpression_gene")
all.gene.dir      <- file.path(pd, "all_gene")
cibersort.dir     <- file.path(pd, "cell_composition")
gsa.dir           <- file.path(pd, "GSA")
gsea.dir          <- file.path(pd, "GSEA")
qc.dir            <- file.path(pd, "QC")
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
GSEA              <- T
gsea.sh           <- NULL
diff.exp          <- T
cibersort         <- T
heatmaps          <- T
pre               <- "TEMP"

sample.key        <- NULL
comps.only        <- NULL
comp.file         <- NULL
counts.file       <- NULL
multi.comps       <- FALSE

norm.counts       <- NULL
cds               <- NULL
Rlibs             <- NULL

cibersortR        <- NULL
cibersortSigFile  <- NULL
request           <- NULL
svnrev            <- NULL
report            <- TRUE

