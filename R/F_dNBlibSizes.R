#'A score function for the row components of the independence model
#'(library sizes)
#'
#'@param beta a vector of length n with current library size estimates
#'@param X a n-by-p count matrix
#'@param reg a vector of length p with relative abundance estimates
#'@param thetas a n-by-p matrix with overdispersion estimates in the rows
#' @param allowMissingness A boolean, are missing values present
#' @param naId The numeric index of the missing values in X
#'
#'@return a vector of length n with evaluations of the score function
dNBlibSizes = function(beta, X, reg, thetas, allowMissingness, naId) {
    mu = exp(outer(beta, reg, "+"))
    X = correctXMissingness(X, mu, allowMissingness, naId)
    rowSums((X - mu)/(1 + (mu/thetas)))
}
