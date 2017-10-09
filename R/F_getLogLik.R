#' A function to extract the logged likelihood of every cell of the data matrix from an RCM object
#'
#' @param rcm an RCM object
#' @param Dim A vector of integers indicating which dimensions to take along, or Inf for the saturated model, or 0 for the independence model
#'
#' @return A matrix with logged likelihood of the size of the data matrix
getLogLik = function(rcm, Dim){
  E = extractE(rcm, Dim)
  thetaMat = extractDisp(rcm, E, Dim)
  dnbinom(x = rcm$X, mu = E, size = thetaMat, log = TRUE)
}