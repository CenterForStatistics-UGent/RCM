#' A function to estimate the taxon-wise NB-params
#'
#' @param design an n-by-v design matrix
#' @param thetas a vector of dispersion parameters of length p
#' @param muMarg an offset matrix
#' @param psi a scalar, the importance parameter
#' @param X the data matrix
#' @param nleqslv.control a list of control elements, passed on to nleqslv()
#' @param ncols an integer, the number of columns of X
#' @param initParam a v-by-p matrix of initial parameter estimates
#' @param v an integer, the number of parameters per taxon
#' @param dynamic a boolean, should response function be determined dynamically?
#'  See details
#' @param envRange a vector of length 2, giving the range of observed
#' environmental scores
#'
#' If dynamic is TRUE, quadratic response functions are fitted for every taxon.
#' If the optimum falls outside of the observed range of environmental scores,
#' a linear response function is fitted instead
#'
#' @return a v-by-p matrix of parameters of the response function
estNBparams = function(design, thetas, muMarg, psi, X, nleqslv.control, ncols,
                       initParam,v, dynamic = FALSE, envRange){
  vapply(seq_len(ncols), FUN.VALUE = vector("numeric",v), function(i){
    nleq = nleqslv(initParam[,i] , reg = design,  fn = dNBllcol_constr,
                   theta = thetas[i], muMarg = muMarg[,i], psi = psi, X = X[,i],
                   control = nleqslv.control, jac = JacCol_constr)$x
    if(dynamic && ((-nleq[2]/(2*nleq[3]) < envRange[1]) ||
                   (-nleq[2]/(2*nleq[3]) > envRange[2]))){
      #If out of observed range, fit a linear model
    nleq = c(nleqslv(initParam[-3,i] , reg = design[,-3], fn = dNBllcol_constr,
                     theta = thetas[i], muMarg = muMarg[,i], psi = psi,
                     X = X[,i], control = nleqslv.control,
                     jac = JacCol_constr)$x,0)
    }
    return(nleq)
  })
}
