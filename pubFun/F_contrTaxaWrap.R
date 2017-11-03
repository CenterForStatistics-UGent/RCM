#' A function to wrap the calculation of the contribution of the taxa to the separation of the clusters over a whole pair of lists
#'
#' @param resListRows the list of row scores
#' @param resListCols the list of column scores
#' @param IDlist a list with id's of taxa that were selected to be DA in that particular instance
#' @param groupFactor the grouping factor that defines the clusters of the samples
#' @param upDown a boolean, were some taxa upregulated and other ones downregulated? Almost always true for non-parametric simulation
#' @param centerFun a function that determines which summary statistic the contribution of the two groups undergoes before taking a ratio
#'
#' @return a dataframe of contributions by method
contrTaxaWrap = function(resListRows, resListCols, IDlist, groupFactor, upDown = FALSE, centerFun = median){
  dat = mapply(resListRows, resListCols, IDlist, SIMPLIFY = TRUE, FUN = function(row, col, id){
    mapply(row, col, SIMPLIFY = TRUE, MoreArgs = list(idTax = id, groupFactor = groupFactor, upDown = upDown, centerFun = centerFun), FUN = contrTaxa)
  })
  as.data.frame(t(dat))
}