#' Derivative of the Lagrangian of the parametric response function
#'
#' @param betas a vector of length (deg+1)*(p+1) with regression parameters with
#'  deg the degree of the response function and the lagrangian multipliers
#' @param X the nxp data matrix
#' @param reg a matrix of regressors with the dimension nx(deg+1)
#' @param thetaMat The n-by-p matrix with dispersion parameters
#' @param muMarg offset matrix of size nxp
#' @param psi a scalar, the importance parameter
#' @param v an integer, one plus the degree of the response function
#' @param p an integer, the number of taxa
#' @param allowMissingness A boolean, are missing values present
#' @param naId The numeric index of the missing values in X
#' @param ... further arguments passed on to the jacobian
#'
#' The parameters are restricted to be normalized, i.e. all squared intercepts,
#'  first order and second order parameters sum to 1
#'
#' @return The evaluation of the score functions, a vector of length (p+1)*
#' (deg+1)
#'
respFunScoreMat = function(betas, X, reg,
    thetaMat, muMarg, psi, p, v, allowMissingness, naId,...) {
    NBparams = matrix(betas[seq_len(p * v)],
        ncol = p)
    mu = exp((reg %*% NBparams) * psi) *
        muMarg
    X = correctXMissingness(X, mu, allowMissingness, naId)
    score = crossprod(reg, (X - mu)/(1 +
        mu/thetaMat)) * psi + 2 * betas[seq_len(v) +
        p * v] * NBparams
    norm = rowSums(NBparams^2) - 1
    return(c(score, norm))  #Taxon per taxon
}
