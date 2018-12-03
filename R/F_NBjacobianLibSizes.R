#'Evaluates the jacobian for the row components of the independence model
#' (library sizes)
#'
#'@param beta a vector of length n with current library size estimates
#'@param X a n-by-p count matrix
#'@param reg a vector of length p with relative abundance estimates
#'@param thetas a n-by-p matrix with overdispersion estimates in the rows
#'
#'@return a diagonal matrix of dimension n: the Fisher information matrix
NBjacobianLibSizes = function(beta, X, reg, thetas){
  mu = exp(outer(beta,reg, "+"))
  diag(-rowSums((1+(X/thetas))*mu/(1+(mu/thetas))^2))
}
