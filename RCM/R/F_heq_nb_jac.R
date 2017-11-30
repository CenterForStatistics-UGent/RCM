#' The jacobian of the linear equality constraints
#'
#' @param Alpha the current estimate of the environmental gradient
#' @param alphaK a matrix with the environmental gradients of the lower dimensions
#' @param d an integer, the number of environmental variables, including dummies
#' @param k an integer, the current dimension
#' @param centMat a centering matrix
#'
#' @return The jacobian matrix
heq_nb_jac = function(Alpha, alphaK, d, k, centMat, ...){
  if(k==1) {return(rbind( centMat, 2*Alpha))
  } else {rbind(centMat, 2*Alpha, t(alphaK))}
}