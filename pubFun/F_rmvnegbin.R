#' Generates correlated NB data, given a covariance matrix
#'
#'@param n number of observations
#'@param mu: n-by-p matrix with means of NB distribution
#'@param Sigma a p-by-p positive definite covariance matrix
#'@param ks vector of length p with overdispersion parameters (size), in the parametrization of the qnbinom() function in base R
#'@param ... additional arguments passed on tot the qnbinom function
#'
#'mu is expected to be a matrix with taxa in the columns and samples in the rows
#'
#'@return A n-by-p matrix of correlated negative binomial dta in the columns
rmvnegbin = function (n, mu, Sigma, ks, ...)
{
  Cor <- cov2cor(Sigma)
  SDs <- sqrt(diag(Sigma))
  if (missing(mu))
    stop("mu is required")
  if (dim(mu)[2] != dim(Sigma)[2])
    stop("Sigma and mu dimensions don't match")
  if (missing(ks)) {
    ks <- unlist(lapply(1:length(SDs), function(i) .negbin_getK(mu[i],
                                                                SDs[i])))
  }
  d <- dim(mu)[2]
  normd <- rmvnorm(n, rep(0, d), Sigma = Cor) #'The normal-to-anything framework
  unif <- pnorm(normd)
  data <- t(qnbinom(t(unif), mu = t(mu), size = ks, ...))
  data <- .fixInf(data)
  return(data)
}