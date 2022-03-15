#'Jacobian for the column components of the independence model
#'
#'@param beta a vector of length p with current abundance estimates
#'@param X a n-by-p count matrix
#'@param reg a vector of length n with library sizes estimates
#'@param thetas a n-by-p matrix with overdispersion estimates in the rows
#' @param allowMissingness A boolean, are missing values present
#' @param naId The numeric index of the missing values in X
#'
#'@return a diagonal matrix of dimension p with evaluations
#'of the jacobian function
NBjacobianAbundsOld = function(beta, X, reg,
    thetas, allowMissingness, naId) {
    mu = exp(outer(reg, beta, "+"))
    X = correctXMissingness(X, mu, allowMissingness, naId)
    -diag(colSums((1 + (X/thetas)) * mu/(1 +
        (mu/thetas))^2))
}
NBjacobianAbunds = function(beta, X, reg,
                            thetas, allowMissingness, naId) {
    mu = exp(reg+beta)
    X = correctXMissingness(X, mu, allowMissingness, naId)
    -as.matrix(sum((1 + (X/thetas)) * mu/(1 +
                                             (mu/thetas))^2))
}
