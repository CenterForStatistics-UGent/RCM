#' A small auxiliary function for the length of the lambdas
#'
#' @param y an integer, the current dimension
#' @param nLambda1s the number of centering restrictions
#'
#' @return a vector containing the ranks of the current lagrangian multipliers
seq_k = function(y, nLambda1s = 1){
  (y-1)*(1+nLambda1s+(y-2)/2) + seq_len(y+nLambda1s)
}