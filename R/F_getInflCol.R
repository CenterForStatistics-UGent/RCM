#' A function to extract the influence for a given parameter index
#'
#' @param score a score matrix
#' @param InvJac The inverted jacobian
#' @param taxon The taxon name or index
#'
#' @return
getInflCol = function(score, InvJac, taxon){
  rowMultiply(score* InvJac[, taxon])
}