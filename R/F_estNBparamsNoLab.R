#' A function to estimate the NB-params ignoring the taxon labels
#'
#' @param design an n-by-v design matrix
#' @param thetasMat A matrix of dispersion parameters
#' @param muMarg an offset matrix
#' @param psi a scalar, the importance parameter
#' @param X the data matrix
#' @param nleqslv.control a list of control elements, passed on to nleqslv()
#' @param initParam a vector of length v of initial parameter estimates
#' @param n an integer, the number of samples
#' @param v an integer, the number of parameters per taxon
#' @param dynamic a boolean, should response function be determined dynamically?
#'  See details
#' @param envRange a vector of length 2,
#'  giving the range of observed environmental scores
#' @param preFabMat a pre-fabricated auxiliary matrix
#' @param allowMissingness A boolean, are missing values present
#' @param naId The numeric index of the missing values in X
#'
#' If dynamic is TRUE, quadratic response functions are fitted for every taxon.
#'  If the optimum falls outside of the observed range of environmental scores,
#'  a linear response function is fitted instead
#'
#' @return a v-by-p matrix of parameters of the response function
estNBparamsNoLab = function(design, thetasMat,
    muMarg, psi, X, nleqslv.control, initParam,
    n, v, dynamic, envRange, preFabMat, allowMissingness, naId) {
    # Without taxon Labels
    nleq = nleqslv(x = initParam, reg = design,
        fn = dNBllcol_constr_noLab, thetasMat = thetasMat,
        muMarg = muMarg, psi = psi, X = X,
        control = nleqslv.control, jac = JacCol_constr_noLab,
        n = n, v = v, preFabMat = preFabMat,
        allowMissingness = allowMissingness, naId = naId)$x
    if (dynamic && ((-nleq[2]/(2 * nleq[3]) <
        envRange[1]) || (-nleq[2]/(2 * nleq[3]) >
        envRange[2]))) {
        # If out of observed range, fit a linear
        # model
        nleq = c(nleqslv(initParam[-3], reg = design[,
            -3], fn = dNBllcol_constr_noLab,
            thetasMat = thetasMat, muMarg = muMarg,
            psi = psi, X = X, control = nleqslv.control,
            jac = JacCol_constr_noLab, preFabMat = preFabMat,
            n = n, v = v - 1, allowMissingness = allowMissingness,
            naId = naId)$x, 0)
    }
    return(nleq)
}
