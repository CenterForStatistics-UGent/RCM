#' A function to build a centering matrix based on a dataframe or an rcm object. It also drops factirs with one level from the dataframe
#'
#' @param object an rcm object or dataframe
#'
#' @return a centering matrix consisting of ones and zeroes, or a list with components
#' \item{centMat}{a centering matrix consisting of ones and zeroes}
#' \item{datFrame}{The dataframe with factors with one level removed}
buildCentMat = function(object){
if(is.data.frame(object)){
  nFactorLevels = sapply(object, function(x){if(is.factor(x)) nlevels(x) else 1}) #Number of levels per factor
  oneLevelID = sapply(object, function(x){length(unique(x))==1})
  object[,oneLevelID] = NULL #Drop factors with one level
  if(any(!oneLevelID)){
    warning("The following variables were not included in the analyses because they have only one value: \n", paste(covariates[oneLevelID], sep = " \n"),immediate. = TRUE)
  }
} else if(class(object) == "RCM"){
  nFactorLevels = sapply(attr(object$covariates, "assign"), function(x){sum(attr(object$covariates, "assign") == x)}) #Number of levels per factor

} else {stop("Invalid object supplied! \n")}
  #Already prepare the matrix that defines the equations for centering the coefficients of the dummy variables
  centMat  = t(sapply(seq_along(nFactorLevels), function(i){
    c(rep.int(0, sum(nFactorLevels[seq(0,i-1)])), #Zeroes before
      rep.int(if(nFactorLevels[i]==1) 0 else 1, nFactorLevels[i]), #Ones within the categorical variable
      rep.int(0, sum(nFactorLevels[-seq(1,i)]))) #Zeroes after
  }))
  centMat = if(all(rowSums(centMat)==0)) {matrix(0, 1,sum(nFactorLevels))} else {centMat[rowSums(centMat)>0,, drop = FALSE]} #Remove zero rows, corresponding to the non-factors
  if(is.data.frame(object)) return(list(centMat = centMat, datFrame  = object)) else return(centMat)
}