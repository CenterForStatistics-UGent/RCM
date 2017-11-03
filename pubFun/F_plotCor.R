#' A function to plot the correlations of the library sizes with the row scores, or of the abundances with the column scores, as boxplots
#'
#' @param scores the list of scores
#' @param datList the list of data matrices, or a list of lists containing the data matrix
#' @param Dims which dimensions to consider
#' @param scoreDim a character vector, "rows" or "columns": which margins to use to calculate the correlations
#' @param dataMat a boolean, is datList a list of data matrices? Otherwise it is a list of lists
#' @param groupMeth a factor defining groups of the methods provided, based on their ordination paradigms
#'
plotCor = function(scores, datList, Dims = 1:3, scoreDim = "rows", dataMat = TRUE, groupMeth = factor(c(as.character(groupsMeth[names(groupsMeth) %in% names(scores[[1]])]), "Control"), levels = c(levels(groupsMeth), "Control"))){
  if(!dataMat){
    datList = lapply(datList, function(x){x$dataMat})
  }
  margins = switch(scoreDim, "rows" = rowSums(datList[[1]]), "columns" = colSums(datList[[1]]))
  cor0 = lapply(Dims, function(Dim){
    rbind(
      mapply(scores, datList, FUN = function(x,z){
        margins = switch(scoreDim, "rows" = rowSums(z), "columns" = colSums(z))
        sapply(x, function(y){libCor(y, margins = margins, Dim = Dim)})
      }, SIMPLIFY = TRUE),
      "Control" = replicate(length(datList),cor(rnorm(length(margins)), margins)))}) #Include a control

  cor0df = lapply(cor0, melt, value.name = "Correl", varnames = c("Method"))
  par(mfrow = c(1,max(Dims)))
  for (i in Dims){
    if(i==max(Dims)){
      parTmp = par(no.readonly = TRUE)
      par(mar = c(4,4,4,5))
    }
    data = orderDF(cor0df[[i]])
    boxplot(Correl ~ Method, main = paste(
      switch(scoreDim, "rows" = "Correlations of library sizes \n with row scores of dimension",
             "columns" = "Correlations of species abundances \n with column scores of dimension"
      ),i),
      data = data, ylim = c(-1,1), ylab = "Pearson correlation", col = groupMeth, las = 2)
    meansCor = tapply(data$Correl, data$Method, mean)
    points(meansCor, col="red", pch=18, cex = 1.6)
    abline(h=0, lty = "dashed")
    if(i==max(Dims)){
      par(parTmp)
    }
  }
  addLegend(groupMeth = groupMeth, x = length(groupMeth), cex = 0.6)
  par(mfrow = c(1,1))

}