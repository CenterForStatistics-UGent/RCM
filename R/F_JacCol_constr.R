#' Jacobian of the constrained analysis with linear response function.
#'
#' @param betas a vector of v parameters of the response function
#'  of a single taxon
#' @param X the count vector of length n
#' @param reg a n-by-v matrix of regressors
#' @param theta The dispersion parameter of this taxon
#' @param muMarg offset of length n
#' @param psi a scalar, the importance parameter
#' @param allowMissingness A boolean, are missing values present
#' @param naId The numeric index of the missing values in X
#'
#' Even though this approach does not imply normalization over
#' the parameters of all taxa, it is very fast
#'  and they can be normalized afterwards
#'
#' @return The jacobian, a square symmetric matrix of dimension v
JacCol_constr = function(betas, X, reg, theta,
    muMarg, psi, allowMissingness, naId) {
    mu = exp(c(reg %*% betas) * psi) * muMarg
    X = correctXMissingness(X, mu, allowMissingness, naId)
    tmp = (1 + X/theta) * mu/(1 + mu/theta)^2
    -crossprod(tmp * reg, reg) * psi^2  #Don't forget to square psi!

}
