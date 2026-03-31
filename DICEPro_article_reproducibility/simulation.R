#' Generate Cell-Type Proportion Matrix
#'
#' Simulates cell-type proportion matrices for bulk RNA-seq deconvolution
#' benchmarking. Three scenarios are supported: equal proportions across
#' cell types, uniform random sampling, or hierarchical Dirichlet sampling
#' reflecting realistic tissue compartment organization.
#'
#' @param n_cell_types Integer. Number of distinct cell types.
#' @param nSample      Integer. Number of samples to generate.
#' @param nCell        Integer. Total number of cells per sample (used for
#'                     rounding precision in the uniform scenario).
#' @param scenario     Character. Proportion generation strategy. One of:
#'                     \code{"even"} for near-equal proportions,
#'                     \code{"uniform"} for random uniform sampling, or
#'                     \code{"hierarchical"} for two-level Dirichlet sampling
#'                     reflecting tissue compartment structure.
#'                     Any other value defaults to a flat Dirichlet draw.
#'
#' @return A numeric matrix of dimensions \code{nSample x n_cell_types},
#'         where each row sums to 1. Row names are \code{Sample_1}, ...,
#'         \code{Sample_nSample}; column names are \code{CellType_1}, ...,
#'         \code{CellType_n_cell_types}.
#'
#' @details
#' The \code{"hierarchical"} scenario implements the two-level Dirichlet model
#' described in the DicePro simulation framework. At the first level, samples
#' are drawn from five major tissue compartments (Immune, Stromal, Endothelial,
#' Epithelial, Muscle) using compartment-specific concentration parameters
#' reflecting expected relative abundances. At the second level, each
#' compartment is further decomposed into its constituent cell sub-populations.
#' This hierarchical structure produces realistic covariance patterns between
#' cell types that are absent from flat Dirichlet draws.
#'
#' @export
#'
#' @importFrom MCMCpack rdirichlet
#' @importFrom stats rnorm runif aggregate
#'
#' @examples
#' if (interactive()) {
#'   set.seed(2101)
#'   # Equal proportions
#'   prop_even <- generateProp(nSample = 20, n_cell_types = 10,
#'     scenario = "even")
#'
#'   # Hierarchical Dirichlet (recommended for realistic benchmarking)
#'   prop_hier <- generateProp(nSample = 20, n_cell_types = 34,
#'     scenario = "hierarchical")
#'   all(abs(rowSums(prop_hier) - 1) < 1e-8)  # TRUE
#' }
generateProp <- function(n_cell_types,
                         nSample,
                         nCell    = 500,
                         scenario = NULL) {

  if (scenario == "even") {
    # ── Near-equal proportions with small Gaussian perturbation ──────────
    m <- round(
      matrix(
        abs(rnorm(n_cell_types,
          mean = 1 / n_cell_types,
          sd   = 0.01)),
        ncol = n_cell_types
      ), 3
    )

  } else if (scenario == "uniform") {
    # ── Uniform random sampling with variable number of active cell types ─
    min.percentage <- 1
    max.percentage <- 99
    CT             <- paste0("CellType_", seq(n_cell_types))
    P_all          <- data.frame()

    for (i in seq_len(nSample)) {
      num.CT      <- sample(x = 2:length(CT), 1)
      selected.CT <- sample(CT, num.CT, replace = FALSE)

      P <- runif(num.CT, min.percentage, max.percentage)
      P <- round(P / sum(P), digits = log10(nCell))
      P <- data.frame(CT       = selected.CT,
        expected = P,
        stringsAsFactors = FALSE)

      missing.CT <- data.frame(
        CT       = CT[!CT %in% selected.CT],
        expected = rep(0, length(CT) - num.CT),
        stringsAsFactors = FALSE
      )

      P_all <- rbind.data.frame(P_all, rbind(P, missing.CT))
    }

    m <- matrix(
      aggregate(P_all$expected,
        list(P_all$CT),
        FUN = mean)$x,
      ncol = n_cell_types
    )

  } else if (scenario == "hierarchical") {
    # ── Two-level hierarchical Dirichlet model ────────────────────────────
    # Level 1: major tissue compartments
    alpha_groups <- c(
      Immune      = 6.0,
      Stromal     = 2.0,
      Endothelial = 1.5,
      Epithelial  = 1.5,
      Muscle      = 1.0
    )
    groups <- MCMCpack::rdirichlet(nSample, alpha_groups)
    colnames(groups) <- names(alpha_groups)

    # Level 2: immune sub-populations
    alpha_immune <- c(
      B_cell_naive      = 1.8,
      B_cell_memory     = 1.2,
      T_CD4_naive       = 2.0,
      T_CD4_memory      = 1.5,
      T_CD8_memory      = 1.5,
      NK_cell           = 1.0,
      Monocyte_CD14     = 1.0,
      Macrophage        = 1.0,
      Mature_neutrophil = 1.0
    )
    immune <- MCMCpack::rdirichlet(nSample, alpha_immune) *
      groups[, "Immune"]

    # Level 2: stromal sub-populations
    stromal <- MCMCpack::rdirichlet(nSample, rep(2, 8)) *
      groups[, "Stromal"]

    # Level 2: endothelial sub-populations
    endothelial <- MCMCpack::rdirichlet(nSample, rep(2, 3)) *
      groups[, "Endothelial"]

    # Level 2: epithelial sub-populations
    epithelial <- MCMCpack::rdirichlet(nSample, c(2, 2, 2, 1.5, 1.5)) *
      groups[, "Epithelial"]

    # Level 2: muscle sub-populations
    muscle <- MCMCpack::rdirichlet(nSample, rep(1, 9)) *
      groups[, "Muscle"]

    m <- cbind(immune, stromal, endothelial, epithelial, muscle)

    # adapter au nombre demandé
    if (ncol(m) > n_cell_types) {
      m <- m[, seq_len(n_cell_types), drop = FALSE]
    }

    if (ncol(m) < n_cell_types) {
      extra <- MCMCpack::rdirichlet(nSample, rep(1, n_cell_types - ncol(m)))
      m <- cbind(m, extra)
    }
    # Override column naming to match generic CellType_k convention
    n_cell_types <- ncol(m)

  } else {
    # ── Flat Dirichlet draw (default fallback) ────────────────────────────
    m <- MCMCpack::rdirichlet(n = nSample, alpha = rep(1, n_cell_types))

  }

  dimnames(m) <- list(
    paste0("Sample_",   seq_len(nrow(m))),
    paste0("CellType_", seq_len(ncol(m)))
  )

  return(m)
}


#' Generate a Reference Signature Matrix
#'
#' Constructs a synthetic gene-by-cell-type reference matrix using either
#' Poisson-based or log-normal gene expression models, with optional block
#' sparsity and TPM normalization.
#'
#' @param loi             Character. Expression model: \code{"rpois"} for
#'                        Poisson-based simulation or any other value for
#'                        log-normal simulation. Default: \code{"rpois"}.
#' @param tpm             Logical. If \code{TRUE}, columns are TPM-normalized.
#'                        Default: \code{FALSE}.
#' @param bloc            Logical. If \code{TRUE}, introduces block sparsity
#'                        so that each cell type expresses only a subset of
#'                        genes. Default: \code{FALSE}.
#' @param nGenesByCellType Integer. Number of genes per cell-type block
#'                        (used when \code{bloc = TRUE}).
#' @param nCell           Integer. Total number of cells (unused directly;
#'                        retained for API compatibility).
#' @param nCellsType      Integer. Number of cell types. Default: \code{10}.
#' @param nGenes          Integer. Number of genes. Default: \code{500}.
#' @param lam             Numeric vector. Lambda parameters for Poisson
#'                        simulation (passed to
#'                        \code{SimMultiCorrData::rcorrvar2}).
#' @param pois_eps        Numeric vector. Epsilon parameters for Poisson
#'                        simulation.
#' @param corr            Numeric matrix. Inter-cell-type correlation matrix.
#'                        If \code{NULL}, generated via
#'                        \code{clusterGeneration::rcorrmatrix}.
#' @param method          Character. Simulation method passed to
#'                        \code{rcorrvar2}. Default: \code{"Polynomial"}.
#' @param sparse          Logical. If \code{TRUE}, introduces random sparsity
#'                        by zeroing entries with probability
#'                        \code{1 - prob_sparse}. Default: \code{FALSE}.
#' @param prob_sparse     Numeric. Probability of a non-zero entry when
#'                        \code{sparse = TRUE}. Default: \code{NULL}.
#'
#' @return A numeric matrix of dimensions \code{nGenes x nCellsType}.
#'         Row names are \code{Gene_1}, ..., \code{Gene_nGenes};
#'         column names are \code{CellType_1}, ..., \code{CellType_nCellsType}.
#'
#' @export
#'
#' @importFrom dplyr group_by summarise across everything
#' @importFrom SimMultiCorrData rcorrvar2 calc_theory
#' @importFrom clusterGeneration rcorrmatrix
#' @importFrom MASS mvrnorm
#' @importFrom stats runif rbinom
#'
#' @examples
#' if (interactive()) {
#'   set.seed(2101)
#'   ref <- generate_ref_matrix(loi = "gauss", nGenes = 50, nCellsType = 5)
#'   dim(ref)  # 50 x 5
#' }
generate_ref_matrix <- function(loi            = "rpois",
                                tpm            = FALSE,
                                bloc           = FALSE,
                                nGenesByCellType = 50,
                                nCell          = 500,
                                nCellsType     = 10,
                                nGenes         = 500,
                                lam            = NULL,
                                pois_eps       = NULL,
                                corr           = NULL,
                                method         = "Polynomial",
                                sparse         = FALSE,
                                prob_sparse    = NULL) {

  if (is.null(corr))
    corr <- clusterGeneration::rcorrmatrix(nCellsType)

  if (loi == "rpois") {
    # ── Poisson-based simulation via correlated variable generation ───────
    Dist   <- c("Logistic", "Weibull")
    Params <- list(c(0, 1), c(3, 5))
    Stcum  <- do.call(rbind, mapply(
      SimMultiCorrData::calc_theory,
      Dist, Params, SIMPLIFY = FALSE
    ))
    rownames(Stcum) <- Dist
    colnames(Stcum) <- c("mean", "sd", "skew", "skurtosis", "fifth", "sixth")

    Six  <- list(seq(1.7, 1.8, 0.01), seq(0.10, 0.25, 0.01))
    size <- 2
    prob <- 0.75
    nb_eps <- 0.0001

    counts <- 40 * SimMultiCorrData::rcorrvar2(
      n        = nGenes,
      k_pois   = nCellsType,
      method   = method,
      means    = Stcum[, 1],
      vars     = Stcum[, 2]^2,
      skews    = Stcum[, 3],
      skurts   = Stcum[, 4],
      fifths   = Stcum[, 5],
      sixths   = Stcum[, 6],
      Six      = Six,
      lam      = lam,
      pois_eps = pois_eps,
      size     = size,
      prob     = prob,
      nb_eps   = nb_eps,
      rho      = corr,
      seed     = 1234,
      errorloop = TRUE
    )$Poisson_variables

  } else {
    # ── Log-normal simulation via multivariate normal ─────────────────────
    mu     <- runif(nCellsType, -1, 0.5)
    counts <- exp(
      MASS::mvrnorm(n = nGenes, mu = mu, Sigma = corr, empirical = TRUE) + 6
    )
  }

  colnames(counts) <- paste0("CellType_", seq_len(nCellsType))
  rownames(counts) <- paste0("Gene_",     seq_len(nGenes))

  if (sparse)
    counts <- apply(counts, 2, function(x)
      x * rbinom(n = nGenes, size = 1, prob = prob_sparse))

  if (bloc) {
    listNGE <- seq(1, nGenes + 1, nGenesByCellType)
    for (ind in seq_len(nCellsType))
      counts[-(listNGE[ind]:(listNGE[ind + 1] - 1)), ind] <- 0
  }

  if (tpm)
    counts <- t(1e6 * t(counts) / colSums(counts))

  return(as.matrix(counts))
}


#' Simulate Bulk RNA-seq Data with Biological and Technical Noise
#'
#' Generates synthetic bulk RNA-seq expression matrices by linearly mixing
#' reference cell-type profiles according to simulated proportions, with
#' optional biological and technical noise. The noise model is consistent
#' with the Gaussian likelihood framework underlying DicePro's objective
#' function: biological noise is applied multiplicatively and technical
#' noise additively, so that normalization to a Gaussian scale prior to
#' deconvolution is a formal requirement.
#'
#' @param W               Numeric matrix (genes x cell types) or \code{NULL}.
#'                        If \code{NULL}, a reference matrix is generated
#'                        internally via \code{\link{generate_ref_matrix}}.
#' @param prop            Numeric matrix (samples x cell types) or \code{NULL}.
#'                        If \code{NULL}, proportions are generated via
#'                        \code{\link{generateProp}}.
#' @param nSample         Integer. Number of samples. Default: \code{50}.
#' @param nCell           Integer. Total cells per sample. Default: \code{500}.
#' @param nCellsType      Integer. Number of cell types. Default: \code{50}.
#' @param nGenes          Integer. Number of genes. Default: \code{500}.
#' @param lam             Numeric. Lambda for Poisson model (see
#'                        \code{\link{generate_ref_matrix}}).
#' @param pois_eps        Numeric. Epsilon for Poisson model.
#' @param corr            Numeric matrix. Inter-cell-type correlation.
#' @param method          Character. Simulation method. Default:
#'                        \code{"Polynomial"}.
#' @param scenario        Character. Proportion scenario passed to
#'                        \code{\link{generateProp}}. Use
#'                        \code{"hierarchical"} for realistic tissue
#'                        composition (recommended).
#' @param loi             Character. Expression law for reference generation.
#'                        Default: \code{"rpois"}.
#' @param tpm             Logical. TPM normalization of reference.
#' @param bloc            Logical. Block sparsity in reference.
#' @param nGenesByCellType Integer. Genes per block (when \code{bloc = TRUE}).
#' @param sparse          Logical. Random sparsity in reference.
#' @param prob_sparse     Numeric. Sparsity probability.
#' @param sigma_bio       Numeric. Standard deviation of multiplicative
#'                        biological noise. Default: \code{0.07}.
#' @param sigma_tech      Numeric. Standard deviation of additive technical
#'                        noise (proportional to expression level).
#'                        Default: \code{0.07}.
#' @param seed            Integer. Random seed for noise reproducibility.
#'                        Default: \code{1234}.
#'
#' @return A named list with three elements:
#' \describe{
#'   \item{\code{prop}}{data.frame. True cell-type proportions
#'         (samples x cell types).}
#'   \item{\code{reference}}{data.frame. Reference signature matrix
#'         (genes x cell types).}
#'   \item{\code{bulk_noise}}{data.frame. Noisy bulk expression matrix
#'         (genes x samples), with biological and technical noise applied.}
#' }
#'
#' @details
#' The bulk expression matrix is generated as:
#' \deqn{\mathbf{B} = \mathbf{W} \cdot \mathbf{P}^\top}
#' where \eqn{\mathbf{W}} is the reference matrix and \eqn{\mathbf{P}}
#' the proportion matrix. Noise is then applied as:
#' \deqn{B^*_{gn} = B_{gn}(1 + \varepsilon^{\text{bio}}_{gn}) +
#'       \varepsilon^{\text{tech}}_{gn} \cdot B_{gn}}
#' where \eqn{\varepsilon^{\text{bio}} \sim \mathcal{N}(0, \sigma^2_{\text{bio}})}
#' and \eqn{\varepsilon^{\text{tech}} \sim \mathcal{N}(0, \sigma^2_{\text{tech}})}.
#' Negative values are clipped to zero. This noise structure motivates
#' Gaussian normalization of both bulk and reference inputs prior to
#' deconvolution with DicePro.
#'
#' @export
#'
#' @importFrom stats rpois rnorm
#'
#' @examples
#' if (interactive()) {
#'   set.seed(2101)
#'   sim <- simulation(
#'     scenario   = "hierarchical",
#'     nSample    = 20,
#'     nGenes     = 100,
#'     nCellsType = 10,
#'     sigma_bio  = 0.07,
#'     sigma_tech = 0.07,
#'     seed       = 2101
#'   )
#'   dim(sim$prop)       # 20 x 10
#'   dim(sim$bulk_noise) # 100 x 20
#' }
simulation <- function(W              = NULL,
                       prop           = NULL,
                       nSample        = 50,
                       nCell          = 500,
                       nCellsType     = 50,
                       nGenes         = 500,
                       lam            = NULL,
                       pois_eps       = NULL,
                       corr           = NULL,
                       method         = "Polynomial",
                       scenario       = NULL,
                       loi            = "rpois",
                       tpm            = FALSE,
                       bloc           = FALSE,
                       nGenesByCellType = 50,
                       sparse         = FALSE,
                       prob_sparse    = 0.5,
                       sigma_bio      = 0.07,
                       sigma_tech     = 0.07,
                       seed           = 1234) {
  if (is.null(W)) {
    W <- generate_ref_matrix(
      loi              = loi,
      nCellsType       = nCellsType,
      nGenes           = nGenes,
      lam              = lam,
      pois_eps         = pois_eps,
      corr             = corr,
      tpm              = tpm,
      bloc             = bloc,
      nGenesByCellType = nGenesByCellType,
      sparse           = sparse,
      prob_sparse      = prob_sparse
    )
  } else {
    nGenes     <- nrow(W)
    nCellsType <- ncol(W)
    dimnames_ <- dimnames(W)
    W          <- apply(W, 2, as.numeric)
    dimnames(W) <- dimnames_
  }

  if (is.null(prop)) {
    prop <- generateProp(
      n_cell_types = nCellsType,
      nSample      = nSample,
      nCell        = nCell,
      scenario     = scenario
    )
  }
  prop <- as.matrix(prop)

  bulk <- as.matrix(W) %*% t(prop)
  rownames(bulk) <- rownames(W)

  set.seed(seed)

  bio_noise  <- matrix(
    rnorm(length(bulk), mean = 0, sd = sigma_bio),
    nrow = nrow(bulk), ncol = ncol(bulk)
  )
  bulk_noisy <- bulk * (1 + bio_noise)

  tech_noise <- matrix(
    rnorm(length(bulk), mean = 0, sd = sigma_tech),
    nrow = nrow(bulk), ncol = ncol(bulk)
  )
  bulk_noisy <- bulk_noisy + tech_noise * bulk
  bulk_noisy[bulk_noisy < 0] <- 0

  rownames(bulk_noisy) <- rownames(bulk)
  colnames(bulk_noisy) <- colnames(bulk)
  colnames(prop) <- colnames(W)
  rownames(prop) <- colnames(bulk_noisy) <- paste("sample_", 1:nSample)
  return(list(
    p       = prop,
    W  = as.data.frame(W),
    B = as.data.frame(bulk_noisy)
  ))
}
