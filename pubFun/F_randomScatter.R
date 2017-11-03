#'A function to make a scatterplot of a randomly selected MC instance scores, either sample or taxon scores
#'
#'@param list the list with scores of all possible ordination methods
#'@param Dim the dimensions to plot
#'@param col the colour of the plot
randomScatter = function(list, Dim = c(1,2), col = "black"){
  parTmp = par(no.readonly = TRUE)
  if(length(list)<=9) {
    par(mfrow = c(3,3), mar = c(3,3,3,3))
  } else if(length(list)<=12){
    par(mfrow = c(3,4), mar = c(3,3,3,3))
  } else {
    par(mfrow = c(4,4), mar = c(3,3,3,3))
  }
  ID = sample(size=1, seq_along(list))
  obj = list[[ID]]
  lapply(names(obj),function(x){plot(obj[[x]][,Dim], main=x, col = col, xlab = paste("Dim", Dim[1]),ylab = paste("Dim", Dim[2]), pty = "s")})
  par(parTmp)
}