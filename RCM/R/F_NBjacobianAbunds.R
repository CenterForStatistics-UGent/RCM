#'A jacobian function for the column components of the independence model (mean relative abundances)
#'
#'@param beta a vector of length p with current abundance estimates
#'@param X a n-by-p count matrix
#'@param reg a vector of length n with library sizes estimates
#'@param thetas a n-by-p matrix with overdispersion estimates in the rows
#'
#'@return a diagonal matrix of dimension p with evaluations of the jacobian function
NBjacobianAbunds = function(beta, X, reg, thetas){
mu = exp(outer(reg,beta, "+"))
-diag(colSums((1+(X/thetas))*mu/(1+(mu/thetas))^2))
}