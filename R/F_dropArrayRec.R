#' A function to drop the recorded parameter estimates throughout the iterations to save space
#'
#' @param RCM an object of the RCM class
#'
#' @return The same object, but with all records from every iteration dropped
dropArrayRec = function(RCM){
  class(RCM)= "list"
  tmp =  within(RCM, {
    rowRec = colRec = psiRec = thetaRec = alphaRec = NULL
  })
  class(tmp)= "RCM"
  tmp
}