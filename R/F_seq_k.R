#' A small auxiliary function for the length of the lambdas
#'
#' @param y an integer, the current dimension
#'
#' @return a vector that goes incrementally from 1 to the number of lagrangian multipliers
seq_k = function(y){
  (y-1)*(2+(y-2)/2) + seq_len(y+1)
}