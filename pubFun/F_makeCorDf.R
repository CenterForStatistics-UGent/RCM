#' A function that makes the correlation dataframe for plotting
#'
#' @param scores the list of scores
#' @param datList the list of data matrices, or a list of lists containing the data matrix
#' @param Dims which dimensions to consider
#' @param scoreDim a character vector, "rows" or "columns": which margins to use to calculate the correlations
#' @param dataMat a boolean, is datList a list of data matrices? Otherwise it is a list of lists
#' @param groupMeth a factor defining groups of the methods provided, based on their ordination paradigms
#' @param dispersion a boolean, have overdispersion parameters been provided?
#'
#'@return a list of dataframes with "Correl" and "Method" components
makeCorDf = function(scores, datList, Dims = 1:3, scoreDim = "rows", dataMat = TRUE,
                     groupMeth = droplevels(factor(c(as.character(groupsMeth[names(groupsMeth) %in% names(scores[[1]])]), "Control"), levels = c(levels(groupsMeth), "Control"))), dispersion = FALSE){
  if(!dataMat){
    datList = lapply(datList, function(x){x$dataMat})
  }
  cor0 = lapply(Dims, function(Dim){
      mapply(scores, datList, FUN = function(x,z){
        margins = switch(scoreDim, "rows" = rowSums(z), "columns" = colSums(z))
        x$Control = matrix(rnorm(length(margins)), length(margins), length(Dims))
        sapply(x , function(y){libCor(y, margins = margins, Dim = Dim)}) #Include a control
      }, SIMPLIFY = TRUE)
      })

  cor0df = lapply(cor0, melt, value.name = "Correl", varnames = c("Method"))
  cor0df
}