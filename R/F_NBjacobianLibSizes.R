#' Jacobian for the raw components of the independence model
#'
#'@param beta a vector of length n with current library size estimates
#'@param X a n-by-p count matrix
#'@param reg a vector of length p with relative abundance estimates
#'@param thetas a n-by-p matrix with overdispersion estimates in the rows
#' @param allowMissingness A boolean, are missing values present
#' @param naId The numeric index of the missing values in X
#'
#'@return a diagonal matrix of dimension n: the Fisher information matrix
NBjacobianLibSizes = function(beta, X, reg,
    thetas, allowMissingness, naId) {
    mu = exp(outer(beta, reg, "+"))
    X = correctXMissingness(X, mu, allowMissingness, naId)
    diag(-rowSums((1 + (X/thetas)) * mu/(1 +
        (mu/thetas))^2))
}
