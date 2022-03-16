#'Gram-Schmidt orthogonalization of vectors
#'
#' @param x The vector that is to be orthogonalized
#' @param otherVecs a matrix; x is orthogonalized with respect to its rows
#' @param weights The weights used in the orthogonalization
#' @return The orthogonalized vector
GramSchmidt = function(x, otherVecs, weights = rep(1, length(x))){
  for(i in seq_len(nrow(otherVecs))){
    x = x- sum(x*otherVecs[i,]*weights)/sum(otherVecs[i,]^2*weights)*otherVecs[i,]
    x = x/sqrt(sum(x^2*weights))
  }
  return(x)
}
