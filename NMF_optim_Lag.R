library(pracma)  # Pour la fonction `expm`
library(MASS)    # Pour la fonction `mvrnorm`
library(ggplot2)
library(DICEPRO)
library(reshape2)




obj_fun <- function(theta, V_obs, W, nUnknownSup=1) {

  ngenes <- nrow(W)
  npop <- ncol(W)
  nsamples <- ncol(V_obs)

  Wunknown <- matrix(theta[1:(ngenes * nUnknownSup)], nrow = ngenes, ncol = nUnknownSup)
  H <- matrix(theta[(ngenes * nUnknownSup) + 1:((npop + nUnknownSup)*nsamples)], ncol = nsamples, nrow = npop + nUnknownSup)
  lambda <- theta[ngenes*nUnknownSup + nsamples*(npop + nUnknownSup) + 1]
  gamma <- theta[ngenes*nUnknownSup + nsamples* (npop + nUnknownSup) + 2]
  sigma2 <- theta[ngenes*nUnknownSup + nsamples* (npop + nUnknownSup) + 3]

  V_optim <- cbind(W, Wunknown) %*% H
  constraints <- colSums(H) - 1
  penalty <- lambda * sum(constraints) + (gamma / 2) * sum(constraints^2)

  frobenius_norm <- sum((V_optim - V_obs)^2)
  mlog_likelihood_value <- 1/(2 * sigma2) * frobenius_norm + prod(dim(V_optim)) * log(2 * pi * sigma2) + penalty

  return(mlog_likelihood_value)
}

grad_obj_fun <- function(theta, V_obs, W, nUnknownSup=1) {

  ngenes <- nrow(W)
  npop <- ncol(W)
  nsamples <- ncol(V_obs)

  Wunknown <- matrix(theta[1:(ngenes * nUnknownSup)], nrow = ngenes, ncol = nUnknownSup)
  H <- matrix(theta[(ngenes * nUnknownSup) + 1:((npop + nUnknownSup)*nsamples)], ncol = nsamples, nrow = npop + nUnknownSup)
  lambda <- theta[ngenes*nUnknownSup + nsamples*(npop + nUnknownSup) + 1]
  gamma <- theta[ngenes*nUnknownSup + nsamples* (npop + nUnknownSup) + 2]
  sigma2 <- theta[ngenes*nUnknownSup + nsamples* (npop + nUnknownSup) + 3]

  Wall <- cbind(W, Wunknown)
  V_optim <- Wall %*% H

  grad_Wunknown <- (1 / sigma2) * (V_obs - V_optim) %*% t(H[npop + 1:nUnknownSup, , drop=FALSE])

  sum_diff <- rowSums(H) - 1
  grad_lambda <- sum(sum_diff)

  grad_gamma <- 0.5 * sum(sum_diff^2)

  error <- sum((V_obs - V_optim)^2)
  grad_sigma <- error / (2*sigma2^2) - prod(dim(V_obs)) / (2*sigma2)

  constraints <- colSums(H) - 1
  penalty_grad <- lambda + gamma * constraints

  grad_H <- (1 / sigma2) * t(Wall) %*% (V_obs - V_optim) + penalty_grad
  # for (s in 1:length(penalty_grad)) {
  #   grad_H[, s] <- grad_H[, s] - penalty_grad[s]
  # }

  grads <- c(as.vector(grad_Wunknown),
            as.vector(grad_H),
            grad_lambda,
            grad_gamma,
            grad_sigma
  )
  return(grads)
}


nmf_conjugate_gradient <- function(V, Winit, Hinit = NULL, nUnknownSup = 1,  ncores = detectCores() - 1,
                                   con = list(fnscale = 1, maxit = 100, tmax = 10)) {
  k <- ncol(Winit)
  if (is.null(Winit)) {
    stop("Winit must be provided and cannot be NULL.")
  }

  if (is.null(Hinit)) {
    Hinit <- matrix(runif(k * ncol(V)), nrow = k, ncol = ncol(V))
  }

  Winit <- as.matrix(Winit)
  Hinit <- as.matrix(Hinit)

  k_CT <- ncol(Hinit) + nUnknownSup
  dimnames_Comp <- list(NULL, paste0("Unknown_", 1:nUnknownSup))
  Winit_C <- matrix(abs(rnorm(n=nrow(Winit)*nUnknownSup, m=mean(Winit), sd=sd(Winit))),
                    nrow = nrow(Winit), ncol = nUnknownSup, dimnames = dimnames_Comp)

  parNMF <- function(i) {
    lambda <- gamma <- 1

    Hinit_C <- matrix(i, nrow = nrow(Hinit), ncol = nUnknownSup, dimnames = dimnames_Comp)
    H <- cbind(Hinit, Hinit_C)
    H <- H/rowSums(H)

    dimnames_H <- dimnames(t(H))
    sigma2 <-  var(as.vector(Winit %*% t(Hinit) - V))
    theta <- c(as.vector(Winit_C), as.vector(t(H)),
               lambda, gamma, sigma2)

    browser()
    result <- optim(par = theta, fn = obj_fun, #gr = grad_obj_fun,
                    W = Winit, V_obs = V, nUnknownSup = nUnknownSup,
                    lower = c(rep(0, length(theta)-3), 0, 0, 0),
                    control =  list(fnscale = 1, maxit = 100, trace=3),
                    method = "L-BFGS-B")

      theta_opt <- result$par

      ngenes <- nrow(Winit)
      npop <- ncol(Winit)
      nsamples <- ncol(V)
      Wunknown <- matrix(theta_opt[1:(ngenes * nUnknownSup)], nrow = ngenes, ncol = nUnknownSup)
      H <- matrix(theta_opt[(ngenes * nUnknownSup) + 1:((npop + nUnknownSup)*nsamples)], ncol = nsamples, nrow = npop + nUnknownSup)
      lambda <- theta_opt[ngenes*nUnknownSup + nsamples*(npop + nUnknownSup) + 1]
      gamma <- theta_opt[ngenes*nUnknownSup + nsamples*(npop + nUnknownSup) + 2]
      sigma2 <- theta_opt[ngenes*nUnknownSup + nsamples*(npop + nUnknownSup) + 3]
      normF <- result$value
      constraints <- colSums(H) - 1
      opt <- list("Wunknown"  = Wunknown,
                  "H" = H,
                  "lambda" = lambda,
                  "gamma" = gamma,
                  "sigma2" = sigma2,
                  "normF" = normF,
                  "constraints" = constraints)
      #sqrt(sum((simulation$pred_CSx[,1:47] -simulation$prop[,1:47])^2))
      #sqrt(sum((t(H)[,1:47] -simulation$prop[,1:47])^2))
      #cor(as.vector(t(H)[,1:47]), as.vector(as.matrix(simulation$prop[,1:47])))
      #cor(as.vector(as.matrix(simulation$pred_CSx[,1:47])), as.vector(as.matrix(simulation$prop[,1:47])))
      return(opt)
    }

  result_list <- pbapply::pblapply(X = seq(0.05, 0.95, 0.05), FUN = parNMF,
                                   cl = ncores)
  return(result_list)
}



set.seed(2101)
simulation <- DICEPRO::simulation(loi = "gauss", scenario = " ", bias = TRUE,nSample = nSample, prop = NULL,
                                  nGenes = nGene, nCellsType = nCellsType)
prop <- simulation$prop
reference <- simulation$reference
bulk <- simulation$B

Winit <- simulation$Winit

pred_CSx <- simulation$pred_CSx
intersectCT <- intersect(colnames(pred_CSx), colnames(prop))
truthProp <- prop[,intersectCT]/rowSums(prop[,intersectCT])

perfCSx <- t(as.matrix(DICEPRO::perfFunction(melt(truthProp, id.vars = NULL)$value,
                                             melt(pred_CSx, id.vars = NULL)$value)))

perfCSx

resultsDICEPRO <- nmf_conjugate_gradient(V = as.matrix(bulk), Winit = Winit, Hinit = as.matrix(pred_CSx),
                                         ncores = NULL)


data_all <- list(B = bulk, Winit = Winit, pred_CSx = pred_CSx, prop = prop)
saveRDS(data_all , "dataSave.rds")
PERF <- data.frame()
for (ind in 1:length(resultsDICEPRO)) {
  pred_DICEPRO <- as.data.frame(resultsDICEPRO[[ind]]$H[,-(ncol(Winit) + 1)])
  pred_DICEPRO <- pred_DICEPRO/rowSums(pred_DICEPRO)
  RES_Perf <- DICEPRO:::perfFunction(melt(truthProp, id.vars = NULL)$value,
                                    melt(pred_DICEPRO, id.vars = NULL)$value)
  PERF <- rbind(PERF, RES_Perf)
}
ind <- which.max(PERF[,1])
perfDICEPRO <- PERF[ind,]
perfDICEPRO

losses <- resultsDICEPRO[[ind]]$normF
iterations <- seq_along(losses)

df <- data.frame(iteration = iterations, loss = losses)

ggplot(df, aes(x = iteration, y = losses)) +
  geom_line() +
  labs(title = "Loss function variation",
       x = "Iteration",
       y = "Loss") +
  theme_minimal()

resultsDICEPRO[[ind]]$H
resultsDICEPRO[[ind]]$lambda
resultsDICEPRO[[ind]]$gamma
resultsDICEPRO[[ind]]$sigma
