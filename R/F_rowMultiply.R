#' A function to efficiently row multiply a a-by-b matrix
#'  by a vector of length b.
#'  More memory intensive but that does not matter with given matrix sizes
#'
#' @param matrix a numeric matrix of dimension a-by-b
#' @param vector a numeric vector of length b
#'
#' t(t(matrix)*vector) but then faster
#'
#' @return a matrix, row multplied by the vector
rowMultiply = function(matrix, vector){
  matrix * rep(vector, rep(nrow(matrix),length(vector)))
}
