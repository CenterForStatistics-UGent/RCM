#' A function that makes the correlation dataframe for plotting
#'
#' @param scores the list of scores
#' @param datList a list of lists containing the data matrix
#' @param Dims which dimensions to consider
#' @param scoreDim a character vector, "rows", "columns" or "dispersions": which margins to use to calculate the correlations, or the dispersions
#' @param groupMeth a factor defining groups of the methods provided, based on their ordination paradigms

#'
#'@return a list of dataframes with "Correl" and "Method" components
makeCorDf = function(scores, datList, Dims = 1:3, scoreDim = "rows", groupMeth = droplevels(factor(c(as.character(groupsMeth[names(groupsMeth) %in% names(scores[[1]])]), "Control"), levels = c(levels(groupsMeth), "Control")))){

  cor0 = lapply(Dims, function(Dim){
      mapply(scores, datList, FUN = function(x,z){
        margins = switch(scoreDim, "rows" = rowSums(z$dataMat), "columns" = colSums(z$dataMat), "dispersions"=z$thetasSampled[colnames(z$dataMat)])
        x$Control = matrix(rnorm(length(margins)), length(margins), length(Dims))
        sapply(x , function(y){libCor(y, margins = margins, Dim = Dim)}) #Include a control
      }, SIMPLIFY = TRUE)
      })

  cor0df = lapply(cor0, function(x){melt(x, value.name = "Correl", varnames = c("Method"))[,c("Correl","Method")]})
  cor0df
}