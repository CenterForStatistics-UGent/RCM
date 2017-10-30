#' A function to extract the influence of all observations on a given component of alpha
#'
#' @param score the score function evaluated for every observation
#' @param InvJac The inverse jacobian
#' @param coef a character string or index referring to the component of alpha
#'
#' @return A matrix with all observations' influence on the row score
getInflRow = function(score, InvJac, coef){
  score* InvJac[, coef]
}