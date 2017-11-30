#' A function to extract the influence of all observations on a given row score
#'
#' @param score the score function evaluated for every observation
#' @param InvJac The inverse jacobian
#' @param sample the row score or  sample index
#'
#' @return A matrix with all observations' influence on the row score
getInflRow = function(score, InvJac, sample){
  score* InvJac[, sample]
}