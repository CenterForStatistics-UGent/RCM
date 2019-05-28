#' Jacobian function for the estimation of a third degree GLM
#'
#' @param beta vector of any length
#' @param X the data vector of length n
#' @param reg a nxlength(beta) regressor matrix
#' @param theta a scalar, the overdispersion
#' @param muMarg the offset of length n
#'
#' @return A matrix of dimension 8-by-8
NBjacobianColNP = function(beta, X, reg,
    theta, muMarg, allowMissingness) {
    # Calculate the mean
    mu = exp(reg %*% beta) * muMarg
    X = correctXMissingness(X, mu, allowMissingness)
    # Return the Jacobian
    -crossprod(reg * c((1 + X/theta) * mu/(1 +
        mu/theta)^2), reg)
}
