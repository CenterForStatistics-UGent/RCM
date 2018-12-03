#' A function that returns the coordinates of an ellipse
#'
#' @param a,b,c parameters of the quadratic function a^2x+bx+c
#' @param quadDrop A scalar, fraction of peak height at which to draw
#'  the ellipse
#' @param nPoints an integer, number of points to use to draw the ellipse
#'
#' @return a matrix with x and y coordinates of the ellipse
ellipseCoord = function(a,b,c, quadDrop = 0.95, nPoints = 100){
  center = -b/(2*a)
  roots = mapply(a,b,c, FUN = function(A,B,C){
    Re(polyroot(c(C - log(ifelse(A>0, 1/quadDrop,quadDrop)) +
                    (B^2-4*A*C)/(4*A), B, A))[1]) #Log to go to count scale.
    #Since function is symmetric, does not matter which root we pick,
    #only distance to vertex matters so we just take the first
  })
  anglesEval = seq(0,2*pi, length.out = nPoints)
  ab = center-roots
  radii = prod(ab)/sqrt(sin(anglesEval)^2*ab[1]^2 + ab[2]^2* cos(anglesEval)^2)
  cbind(x = radii*cos(anglesEval) + center[1],
        y = radii * sin(anglesEval) + center[2])
}
