#' The score function of the response function for 1 taxon at the time
#' @param betas a vector of v parameters of the
#'  response function of a single taxon
#' @param X the count vector of length n
#' @param reg a n-by-v matrix of regressors
#' @param theta The dispersion parameter of this taxon
#' @param muMarg offset of length n
#' @param psi a scalar, the importance parameter
#' @param allowMissingness A boolean, are missing values present
#' @param naId The numeric index of the missing values in X
#'
#'  Even though this approach does not imply normalization over the parameters
#'  of all taxa, it is very fast and they can be normalized afterwards
#'
#' @return A vector of length v with the evaluation of the score functions
dNBllcol_constr = function(betas, X, reg,
    theta, muMarg, psi, allowMissingness, naId) {
    mu = exp(c(reg %*% betas) * psi) * muMarg
    X = correctXMissingness(X, mu, allowMissingness, naId)
    crossprod((X - mu)/(1 + mu/theta), reg) *
        psi
}
