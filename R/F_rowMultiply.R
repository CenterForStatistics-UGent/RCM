#' A function to efficiently row multiply a matrix and a vector
#'
#' @details Memory intensive but that does not matter with given matrix sizes
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
