#' A function to calculate correlations of coordinates with marginal sums
#'
#' @param coord a coordinate matrix
#' @param margins the marginal sums to calculate the correlation with
#' @param Dims the dimensions to use
#'
#' @return A vector of correlations with the same length as Dims
libCor = function(coord, margins, Dim = 1:3){
  if(length(margins)!=nrow(coord)) {return(rep(NA, length(Dim)))}
  apply(coord[, Dim, drop=FALSE], 2, cor, margins)
}