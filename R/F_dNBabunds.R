#'A score function for the column components of the independence model
#' (mean relative abundances)
#'
#'@param beta a vector of length p with current abundance estimates
#'@param X a n-by-p count matrix
#'@param reg a vector of length n with library sizes estimates
#'@param thetas a n-by-p matrix with overdispersion estimates in the rows
#'
#'@return a vector of length p with evaluations of the score function
dNBabunds = function(beta, X, reg, thetas) {
    mu = exp(outer(reg, beta, "+"))
    score = colSums((X - mu)/(1 + (mu/thetas)))
}
