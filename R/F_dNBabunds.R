#'A score function for the column components of the independence model
#' (mean relative abundances)
#'
#'@param beta a vector of length p with current abundance estimates
#'@param X a n-by-p count matrix
#'@param reg a vector of length n with library sizes estimates
#'@param thetas a n-by-p matrix with overdispersion estimates in the rows
#' @param allowMissingness A boolean, are missing values present
#' @param naId The numeric index of the missing values in X
#'
#'@return a vector of length p with evaluations of the score function
dNBabundsOld = function(beta, X, reg, thetas, allowMissingness, naId) {
    mu = exp(outer(reg, beta, "+"))
    X = correctXMissingness(X, mu, allowMissingness, naId)
    score = colSums((X - mu)/(1 + (mu/thetas)))
}
dNBabunds = function(beta, X, reg, thetas, allowMissingness, naId) {
  mu = exp(reg + beta)
  X = correctXMissingness(X, mu, allowMissingness, naId)
  sum((X - mu)/(1 + (mu/thetas)))
}
