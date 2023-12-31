test_that("DICEPRO", {
  skip_on_cran()
  library(clusterGeneration)
  library(NMF)

  nCellsType <- 5
  nGenes <- 50
  nSample <- 10

  simulation <- simulation(loi = "gauss", scenario = " ", bias = TRUE, nSample = nSample, prop = NULL,
                           nGenes = nGenes, nCellsType = nCellsType)

  cellTypeOut <- sample(1:ncol(simulation$reference), 1)
  refDataIncomplet <- simulation$reference[,-cellTypeOut]

  # res <- DICEPRO(reference = refDataIncomplet, bulk = simulation$bulk, nIteration = 2,
  #                  methodDeconv = "DCQ")
})
