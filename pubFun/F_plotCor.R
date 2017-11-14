#' A function to plot the correlations of the library sizes with the row scores, or of the abundances with the column scores, as boxplots
#'
#' @param corDF a list of dataframes as resulting from the makeCorDf function
#' @param groupMeth a factor defining groups of the methods provided, based on their ordination paradigms
#' @param bordercol The colour of the borders of the boxplot
#' @param cex an expansion factor for the legend
#' @param diamCol The color of the diamond for the mean
plotCor = function(corDF, scoreDim = "rows",
                   groupMeth = droplevels(factor(c(as.character(groupsMeth),"Control")[levelsMeth %in% levels(corDF[[1]]$Method)], levels = c(levels(groupsMeth), "Control"))),
                   bordercol = borderCol, cex = 1, diamCol = "orange"){
  parTmp = par(no.readonly = TRUE)
  Dims = seq_along(corDF)
  par(mfrow = c(1,max(Dims)), oma = c(2,2,2,3), mar = c(7,1,3,7.5), pty = "m")
  for (i in Dims){
    data = orderDF(corDF[[i]], groupsMeth = groupMeth)
    boxplot(Correl ~ Method, main = paste(
      switch(scoreDim, "rows" = "Correlations of library sizes \n with row scores of dimension",
             "columns" = "Correlations of species abundances \n with column scores of dimension"
      ),i),
      data = data, ylim = c(-1,1), ylab = "Pearson correlation", col = groupMeth, las = 2, cex.main = 0.75)
    meansCor = tapply(data$Correl, data$Method, mean)
    points(meansCor, col = diamCol, pch=18, cex = 1.6)
    abline(h=0, lty = "dashed")
    if(i==max(Dims)){
      addLegend(groupMeth = groupMeth, x = length(unique(data$Method))+1, cex = cex)
    }
  }
  par(parTmp)

}