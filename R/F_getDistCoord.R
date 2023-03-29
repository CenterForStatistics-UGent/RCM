#' Get coordinates of a distance object of n observations for the provided indices
#'
#' @param indices The row indices for which distance indices are wanted
#' @param n The total number of objects in the distance matrix
#'
#' @return a vector of coordinates
getDistCoord = function(indices, n){
    indexI = n*(indices-1) - indices*(indices+1)/2# See details of ?dist
    indexMat = outer(indexI, indices, FUN = "+")
    indexMat[upper.tri(indexMat)]
}