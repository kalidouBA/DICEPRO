#' This function is used for generating synthetic data representing the proportions of different cell
#' types in biological samples for various scenarios, which can be useful in various data analysis and
#' modeling tasks in the field of biology and genetics.
#'
#' @param n_cell_types The number of different cell types.
#' @param nSample The number of samples to generate.
#' @param nCell The total number of cells in each sample.
#' @param scenario A scenario that determines how the proportions are generated.
#' It can be one of three values: "even," "uniform," or any other value.
#'
#' @export
#'
#' @importFrom rBeta2009 rdirichlet
#' @importFrom stats rnorm runif aggregate
#'
#' @return a matrix representing proportions of different cell types in samples
#'
#' @examples
#' if(interactive()){
#' set.seed(2101)
#' prop <- generateProp(nSample = 20, n_cell_types = 22, scenario = "even")
#' }

generateProp <- function(n_cell_types, nSample, nCell,scenario=NULL){
  if(scenario == "even")
    m <- round(matrix(abs(rnorm(n_cell_types, mean = 1 / n_cell_types,
                                       sd = 0.01)),  ncol = n_cell_types), 3)
  else if(scenario == "uniform"){
    min.percentage = 1
    max.percentage = 99

    CT <- paste0("CellType_",seq(n_cell_types))
    P_all <- data.frame()
    for (i in 1:nSample) {
      num.CT = sample(x = 2:length(CT),1)
      selected.CT = sample(CT, num.CT, replace = FALSE)

      P = runif(num.CT, min.percentage, max.percentage)
      P = round(P/sum(P), digits = log10(nCell))
      P = data.frame(CT = selected.CT, expected = P, stringsAsFactors = FALSE)

      missing.CT = CT[!CT %in% selected.CT]
      missing.CT = data.frame(CT = missing.CT,
                              expected = rep(0, length(missing.CT)),
                              stringsAsFactors = FALSE)

      P = rbind.data.frame(P, missing.CT)

      P_all = rbind.data.frame(P_all, P)
    }

    m <- matrix(aggregate(P_all$expected, list(P_all$CT), FUN=mean)$x,
                ncol = n_cell_types)
  }
  else
    m <- rdirichlet(n=nSample, shape = 1:n_cell_types)

  dimnames(m) <- list(paste0("Sample_",1:nrow(m)),paste0("CellType_",1:ncol(m)))
  return(m)
}


#' Generate a Reference Matrix
#'
#' This function generates a reference matrix for use in simulating single-cell RNA-seq data.
#'
#' @param loi The law of integration to use: "rpois" for Poisson
#' @param tpm Logical indicating whether to convert the counts to transcripts per million (TPM).
#' @param bloc Logical indicating whether to introduce block sparsity.
#' @param nGenesByCellType Number of genes per cell type (used for block sparsity).
#' @param nCell Total number of cells.
#' @param nCellsType Number of cell types.
#' @param nGenes Number of genes.
#' @param lam The lambda parameter for the Poisson distribution (used when loi = "rpois").
#' @param pois_eps The epsilon parameter for the Poisson distribution (used when loi = "rpois").
#' @param corr The correlation matrix for generating gene expressions (used when loi = "rpois" is FALSE).
#' @param method The method for generating data: "Polynomial" or other.
#' @param sparse Logical indicating whether to introduce sparsity.
#' @param prob_sparse Probability of sparsity (used when sparse = TRUE).
#'
#'
#' @export
#'
#' @importFrom matrixStats rowVars
#' @importFrom dplyr group_by summarise across everything
#' @importFrom SimMultiCorrData rcorrvar2
#' @importFrom clusterGeneration rcorrmatrix
#' @importFrom SimMultiCorrData calc_theory
#' @importFrom MASS mvrnorm
#' @importFrom stats runif rbinom
#' @import dplyr
#'
#' @return A reference matrix with cell types as columns and genes as rows.
#'
#' @seealso \code{\link{rcorrvar2}}, \code{\link{mvrnorm}}
#'
#' @examples
#' if(interactive()){
#' set.seed(2101)
#' ref <- generate_ref_matrix(loi = "gauss", nGenes = 50, nCellsType = 5)
#' }

generate_ref_matrix <- function(loi="rpois", tpm = FALSE, bloc = FALSE, nGenesByCellType=50,
                                nCell = 500, nCellsType = 10, nGenes = 500, lam = NULL,
                                pois_eps = NULL,corr = NULL, method = "Polynomial",
                                sparse = FALSE, prob_sparse=NULL){
 if(is.null(corr))
    corr <- rcorrmatrix(nCellsType)

  # rpois
 if(loi == "rpois"){
    Dist <- c("Logistic", "Weibull")
    Params <- list(c(0, 1), c(3, 5))
    Stcum1 <- calc_theory(Dist[1], Params[[1]])
    Stcum2 <- calc_theory(Dist[2], Params[[2]])
    Stcum <- rbind(Stcum1, Stcum2)
    rownames(Stcum) <- Dist
    colnames(Stcum) <- c("mean", "sd", "skew", "skurtosis", "fifth", "sixth")
    Six <- list(seq(1.7, 1.8, 0.01), seq(0.10, 0.25, 0.01))

    size <- 2
    prob <- 0.75
    nb_eps <- 0.0001

    # Simulate variables with error loop
    counts <- 40*rcorrvar2(n = nGenes, k_pois = nCellsType,
                           method = method, means = Stcum[, 1],
                           vars = Stcum[, 2]^2, skews = Stcum[, 3],
                           skurts = Stcum[, 4], fifths = Stcum[, 5],
                           sixths = Stcum[, 6], Six = Six,
                           lam = lam, pois_eps = pois_eps, size = size,
                           prob = prob, nb_eps = nb_eps, rho = corr,
                           seed = 1234, errorloop = TRUE)$Poisson_variables

    colnames(counts) <- paste0("CellType_",rep(1:nCellsType))
    rownames(counts) <- paste0("Gene_",rep(1:nGenes))
  }
  else{
    counts <- matrix(NA,nGenes,nCellsType)
    mu <- runif(nCellsType,-1,.5)
    counts <- exp(mvrnorm(n=nGenes, mu=mu, Sigma=corr,empirical = TRUE)+6)

    colnames(counts) <- paste0("CellType_", 1:nCellsType)
    rownames(counts) <- paste0("Gene_", 1:nGenes)
  }

  if(sparse)
    counts <- apply(counts, 2, FUN = function(x) {
      x*rbinom(n=nGenes, size=1, prob=prob_sparse)
    })

  if(bloc){
      listNGE <- seq(1,nGenes+1,nGenesByCellType)
      for (ind in seq(nCellsType)){
        counts[-(listNGE[ind]:(listNGE[ind+1]-1)),ind] <- 0
      }
  }
  if(tpm)
    counts <- t(1e6*t(counts)/colSums(counts))
  return(as.matrix(counts))
}


#' Simulate bulk RNA-Seq Data
#'
#' This function simulates bulk RNA-seq data using a reference matrix and cell type proportions.
#
#' @param W A reference matrix, or NULL if you want to generate one.
#' @param prop Cell type proportions, or NULL if you want to generate them.
#' @param nSample Number of samples to generate.
#' @param nCell Total number of cells per sample.
#' @param nCellsType Number of cell types.
#' @param nGenes Number of genes.
#' @param lam The lambda parameter for the Poisson distribution (used when loi = "rpois").
#' @param pois_eps The epsilon parameter for the Poisson distribution (used when loi = "rpois").
#' @param corr The correlation matrix for generating gene expressions (used when loi = "rpois" is FALSE).
#' @param method The method for generating data: "Polynomial" or other.
#' @param scenario The scenario for generating cell type proportions: "even," "uniform," or other.
#' @param loi The law of integration to use: "rpois" for Poisson.
#' @param tpm Logical indicating whether to convert the counts to transcripts per million (TPM).
#' @param bloc Logical indicating whether to introduce block sparsity.
#' @param nGenesByCellType Number of genes per cell type (used for block sparsity).
#' @param missCellTypes Vector of cell types to exclude in the simulation.
#' @param sparse Logical indicating whether to introduce sparsity.
#' @param prob_sparse Probability of sparsity (used when sparse = TRUE).
#' @param bias Logical indicating whether to add bias to the generated data.
#' @param mu_bruit Mean value for bias (used when bias = TRUE).
#' @param sigma_bruit Standard deviation for bias (used when bias = TRUE).
#'
#' @export
#'
#' @importFrom scater mockSCE
#'
#' @return A list data simulated.
#' @details The function calculates and returns the simulations:
#' \itemize{
#'   \item \code{prop}: A matrix representing proportions of different cell types in samples.
#'   \item \code{reference}: A reference matrix with cell types as columns and genes as rows.
#'   \item \code{bulk}: A bulk RnaSeq simulated
#'}
#'
#' @examples
#' if(interactive()){
#' set.seed(2101)
#' simulation <- simulation(loi = "gauss", scenario = " ", bias = TRUE, nSample = 10, prop = NULL,
#'                          nGenes = 50, nCellsType = 5)
#'}

simulation <- function(W = NULL, prop = NULL, nSample = 50, nCell = 500, nCellsType = 50,
                        nGenes = 500, lam = NULL, pois_eps = NULL, corr = NULL,
                        method = "Polynomial", scenario = NULL,
                        loi = ' ', tpm = FALSE, bloc = FALSE,
                        nGenesByCellType = 50, missCellTypes = NA,
                        sparse=FALSE, prob_sparse=0.5,
                        bias = FALSE, mu_bruit = 0, sigma_bruit = .025){

  # Generate ref
  if(is.null(W)){
    W <- generate_ref_matrix(loi = loi, nCellsType = nCellsType,nGenes = nGenes,
                             lam = lam,pois_eps = pois_eps,corr = corr)
  }else{
    nGenes <- nrow(W)
    cellTypeName <- colnames(W)
    nCellsType <- ncol(W)
    W[,cellTypeName] <- sapply(W[,cellTypeName],as.numeric)
  }
  dimname_ <- dimnames(W)
  W <- apply(W, 2, as.numeric)
  dimnames(W) <- dimname_

  # Generate prop
  if(is.null(prop)){
    prop = matrix(0,nrow = nSample, ncol = nCellsType)
    if(scenario %in% c("even","uniform"))
      for (i in 1:nSample)
        prop[i,] <- generateProp(nCellsType,scenario = scenario)

    else
      prop <- generateProp(nSample = nSample,n_cell_types = nCellsType,scenario = scenario)

    dimnames_ <- dimnames(prop)
    prop <- apply(prop, 2, as.numeric)
    rownames(prop) <- dimnames_[[1]]
  }
  else
    prop <- apply(prop, 2, as.numeric)
  prop <- as.data.frame(prop)

  # Generate bulk matrix
  GEO <- as.matrix(W)%*%t(prop)

  # Add bias
  if(bias){
    biasAdd <- matrix(NA,nrow = nGenes,ncol = nSample)
    for (nS in 1:nSample) {
      if(loi == "rpois")
        biasAdd[,nS] <- rpois(nGenes,runif(1,0,.3))
      else
        biasAdd[,nS] <- rnorm(nGenes,0,runif(1,1,2))
    }
    biasAdd <- apply(biasAdd, 2, as.integer)
    GEO <- GEO + biasAdd
  }
  GEO <- as.data.frame(apply(GEO, 2, as.numeric))
  rownames(GEO) <- paste0("Gene_", 1:nGenes)
  return(list("prop" = prop, "reference" = as.data.frame(W), "bulk" = GEO))
}
