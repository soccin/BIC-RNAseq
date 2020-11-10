## -----------------------------------------------------------------------------
library(piano)
data("gsa_input")

## -----------------------------------------------------------------------------
head(gsa_input$gsc,10)
head(gsa_input$pvals, 10)
head(gsa_input$directions, 10)

## -----------------------------------------------------------------------------
geneSets <- loadGSC(gsa_input$gsc)
geneSets


## ----results='hide', message=F------------------------------------------------
gsares <- runGSA(gsa_input$pvals,
                 gsa_input$directions,
                 gsc = geneSets,
                 nPerm = 500) # set to 500 for fast run

## ----eval=F-------------------------------------------------------------------
#  exploreGSAres(gsares)

## ----echo=FALSE---------------------------------------------------------------
sessionInfo()

