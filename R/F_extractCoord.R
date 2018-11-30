#' A function to extract plotting coordinates, either for plot.RCM or to export to other plotting software
#'
#' @param RCM an RCm object
#' @param Dim an integer vector of required dimensions
#'
#' The parameters for the ellipses of the quadratic response function come from the parametrization f(x) = a*x^2 + b*x + c
#' For an unconstrained object the row and column coordinates are returned in separate matrices. The row names will correspond to the labels. For a constrained analysis also the variable points are returned. All variables still need to be scaled to optimally fill hte available space
#'
#' @return A list with components
#' \item{samples}{A dataframe of sample scores}
#' \item{species}{A dataframe of column scores, with origin, slope, end and ellipse coordinates as needed}
#' \item{variables}{A dataframe of variable scores, loadings of the environmental gradient}
#' @export
#' @examples
#' data(Zeller)
#' require(phyloseq)
#' tmpPhy = prune_taxa(taxa_names(Zeller)[1:100],
#' prune_samples(sample_names(Zeller)[1:50], Zeller))
#' zellerRCM = RCM(tmpPhy, k = 2, round = TRUE)
#' coordsZeller = extractCoord(zellerRCM)
extractCoord = function(RCM, Dim = c(1,2)){
  # Samples
  constrained = !is.null(RCM$covariates)
  dataSam <- if(constrained){
     data.frame(RCM$covariates %*% RCM$alpha[,Dim] %*% diag(RCM$psis[Dim]))
  } else {
     data.frame(RCM$rMat[, Dim] %*% diag(RCM$psis[Dim]))
  }
  names(dataSam)=paste0("Dim", Dim)

  # Species
  if(constrained) {
    if(RCM$responseFun=="linear"){
      dataTax = data.frame(
        origin1 = -RCM$NB_params[1,,Dim[1]]/RCM$NB_params[2,,Dim[1]] ,
        origin2 = -RCM$NB_params[1,,Dim[2]]/RCM$NB_params[2,,Dim[2]] ,
        slope1 = RCM$NB_params[2,,Dim[1]] , #The gradient
        slope2 = RCM$NB_params[2,,Dim[2]]
      )
       dataTax$end1 = dataTax$origin1 + dataTax$slope1
       dataTax$end2 = dataTax$origin2 + dataTax$slope2
       rownames(dataTax) = colnames(RCM$X)
    } else if (RCM$responseFun == "quadratic"){
      dataTax = data.frame(
        apply(RCM$NB_params[c(2,3),,Dim],c(2,3), function(x){
          a = x[2]; b=x[1]
          -b/(2*a)
        })) #The location of the extrema
      names(dataTax) = c("end1","end2")

      peakHeights = apply(RCM$NB_params[,,Dim],2, function(x){
        A = x[3,]; B = x[2,]; C = x[1,];
        vapply(FUN.VALUE = numeric(1), exp(B^2 -4*A*C)/(4*A), function(y){max(y,1/y)}) #select largest relative departure
      })
      rownames(peakHeights) = c("peak1","peak2")
      dataTax = cbind(dataTax, t(peakHeights))

      #Get ellipse parameters
      dataEllipse = t(Reduce(x = lapply(Dim, function(x) {RCM$NB_params[,,x]}), f = rbind))
      colnames(dataEllipse) = c(vapply(FUN.VALUE = character(length(Dim)), Dim, function(x){paste0(c("c","b","a"),x)}))

      #Rescale peak heights for plotting
      dataTax[, c("peak1","peak2")] = rowMultiply(dataTax[, c("peak1","peak2")], apply(dataTax[, c("end1","end2")],2,function(x){max(abs(x))})/apply(dataTax[, c("peak1","peak2")], 2,max))
      # dataTax[, c("peak1","peak2")] = apply(dataTax[, c("peak1","peak2")], c(1,2), max, 0.0075) # Make sure a line always appears
      dataTax = cbind(dataTax, dataEllipse)
      rownames(dataTax) = colnames(RCM$X)

    } else if(RCM$responseFun  == "nonparametric"){ #For non-parametric response function we cannot plot the taxa
      dataTax = NULL
    } else {
      stop("No valid response function present in this RCM object!")
    }
  } else { #If not constrained
    dataTax = data.frame(cbind(t(RCM$cMat[Dim,]),0,0))
    names(dataTax) = c("end1","end2", "origin1","origin2")
    rownames(dataTax) = colnames(RCM$X)
  }

  #Variables
  if(!constrained) {dataVar =NULL
  } else {
    dataVar = data.frame(RCM$alpha)[,Dim]
  }
  list(samples = dataSam, species = dataTax, variables = dataVar)
}
