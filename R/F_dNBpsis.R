#' A score function for the psi of a given dimension
#'
#' @param beta a scalar, the initial estimate
#' @param X the n-by-p count matrix
#' @param muMarg the nxp offset matrix
#' @param reg the regressor matrix, the outer product of current row
#'  and column scores
#' @param theta a n-by-p matrix with the dispersion parameters
#' @param ... other arguments passed on to the jacobian

#' @return The evaluation of the score function at beta, a scalar
dNBpsis = function(beta, X, reg, theta, muMarg, allowMissingness,
    ...) {
    mu = muMarg * exp(reg * beta)
    X = correctXMissingness(X, mu, allowMissingness)
    sum(reg * (X - mu)/(1 + (mu/theta)))
}
