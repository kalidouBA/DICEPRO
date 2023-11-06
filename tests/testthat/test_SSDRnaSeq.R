test_that("SSDRnaSeq", {
  skip_on_cran()
  library(clusterGeneration)
  library(NMF)

  nCellsType <- 5
  nGenes <- 50
  nSample <- 10

  # GENERATE A RANDOM CORRELATION MATRIX BASED ON RANDOM PARTIAL CORRELATIONS
  corr <- rcorrmatrix(nCellsType)

  simulation <- simulation(loi = "gauss", scenario = " ", bias = TRUE, corr = corr, nSample = nSample, prop = NULL,
                           nGenes = nGenes, nCellsType = nCellsType)

  cellTypeOut <- sample(1:ncol(simulation$reference), 1)
  refDataIncomplet <- simulation$reference[,-cellTypeOut]

  res <- SSDRnaSeq(reference = refDataIncomplet, bulk = simulation$bulk, k_folds = 2, nIteration = 2,
                   methodDeconv = "DCQ")
})
