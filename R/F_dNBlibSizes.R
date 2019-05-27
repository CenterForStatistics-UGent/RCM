#'A score function for the row components of the independence model
#'(library sizes)
#'
#'@param beta a vector of length n with current library size estimates
#'@param X a n-by-p count matrix
#'@param reg a vector of length p with relative abundance estimates
#'@param thetas a n-by-p matrix with overdispersion estimates in the rows
#'
#'@return a vector of length n with evaluations of the score function
dNBlibSizes = function(beta, X, reg, thetas, allowMissingness) {
    mu = exp(outer(beta, reg, "+"))
    X = correctXMissingness(X, mu, allowMissingness)
    rowSums((X - mu)/(1 + (mu/thetas)))
}
