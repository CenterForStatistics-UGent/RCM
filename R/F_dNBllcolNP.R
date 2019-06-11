#'Estimation of the parameters of a third degree GLM
#'
#' @param beta A vector of any length
#' @param X the data vector of length n
#' @param reg a nxlength(beta) regressor matrix
#' @param theta a scalar, the overdispersion
#' @param muMarg the offset of length n
#' @param allowMissingness A boolean, are missing values present
#' @param naId The numeric index of the missing values in X
#' @param ... further arguments passed on to the jacobian

#' @return A vector of the same length as beta with evaluations
#'  of the score function
dNBllcolNP = function(beta, X, reg, theta,
    muMarg, allowMissingness, naId, ...) {
    mu = exp(reg %*% beta) * muMarg
    X = correctXMissingness(X, mu, allowMissingness, naId)
    crossprod(reg, ((X - mu)/(1 + mu/theta)))
}
