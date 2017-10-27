#' A function to extract the influence of all observations on a given row score
#'
#' @param score the vector of row scores
#' @param InvJac The inverse jacobian
#' @param sample the row score or  sample index
#'
#' @return A matrix with all observations' influence on the row score
getInflRow = function(score, InvJac, sample){
  score* InvJac[, sample]
}