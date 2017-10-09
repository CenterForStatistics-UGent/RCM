#' A function that returns dispersions belonging to a RCM object
#'
#' @param rcm an RCM object
#' @param E a matrix of expectations
#' @param k an integer, the required dimensions
#'
#' @return a matrix of dispersions
extractDisp = function(rcm, E, k = rcm$k){
 if(k==rcm$k) {
    matrix(rcm$thetas, byrow = TRUE, nrow = nrow(rcm$X), ncol = ncol(rcm$X))
  } else {
    matrix(estDisp(X = rcm$X, muMarg = E), byrow = TRUE, nrow = nrow(rcm$X), ncol = ncol(rcm$X))
  }
}