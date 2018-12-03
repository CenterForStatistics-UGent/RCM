#' Fit the RC(M) model with the negative binomial distribution.
#'
#' @details Includes fitting of the independence model, filtering out the
#' effect of confounders and fitting the RC(M) components in a constrained
#'  or an unconstrained way for any dimension k.
#'
#' @param X a nxp data matrix
#' @param k an scalar, number of dimensions in the RC(M) model
#' @param rowWeights a character string, either "uniform" or "marginal"
#'  row weights.
#' @param colWeights a character string, either "uniform" or "marginal"
#'  column weights.
#' @param tol a scalar, the relative convergende tolerance for the row scores
#'  and column scores parameters.
#' @param maxItOut an integer, the maximum number
#'  of iterations in the outer loop.
#' @param Psitol a scalar, the relative convergence tolerance
#'  for the psi parameters.
#' @param verbose a boolean, should information on iterations be printed?
#' @param global global strategy for solving non-linear systems, see ?nleqslv
#' @param nleqslv.control a list with control options, see nleqslv
#' @param jacMethod Method for solving non-linear equations, ?see nleqslv.
#' Defaults to Broyden. The difference with the newton method is that
#'  the Jacobian is not recalculated at every iteration,
#'  thereby speeding up the algorithm
#' @param dispFreq an integer, how many iterations the algorithm should wait
#'  before reestimationg the dispersions.
#' @param convNorm a scalar, the norm to use to determine convergence
#' @param prior.df an integer, see estDisp()
#' @param marginEst a character string, either "MLE" or "marginSums",
#' indicating how the independence model should be estimated
#' @param confounders a list with
#' -confounders an nxg matrix with confounders
#' -confoundersFilt an nxh matrix with confounders for filtering,
#' with all levels and without intercept
#' @param covariates an nxd matrix with covariates.
#' If set to null an unconstrained analysis is carried out,
#' otherwise a constrained one.
#' Factors must have been converted to dummy variables already
#' @param centMat a fxd matrix containing the contrasts to center
#'  the categorical variables. f equals the number of continuous variables +
#'  the total number of levels of the categorical variables.
#' @param prevCutOff a scalar the minimum prevalence needed to retain a taxon
#'  before the the confounder filtering
#' @param minFraction a scalar, total taxon abundance should equal minFraction*n
#'  if it wants to be retained before the confounder filtering
#' @param responseFun a characters string indicating the shape
#'  of the response function
#' @param record A boolean, should intermediate parameter estimates be stored?
#' @param control.outer a list of control options
#'  for the outer loop constrOptim.nl function
#' @param control.optim a list of control options for the optim() function
#' @param envGradEst a character string, indicating how the
#' environmental gradient should be fitted. "LR" using the likelihood-ratio
#'  criterion, or "ML" a full maximum likelihood solution
#' @param dfSpline a scalar, the number of degrees of freedom for the splines
#'  of the non-parametric response function, see VGAM::s()
#' @param vgamMaxit an integer,
#'  the maximum number of iteration in the vgam() function
#' @param degree an integer,
#'  the degree of the polynomial fit if the spline fit fails
#'
#' #'@seealso \code{\link{RCM}}
#'
#' Not intended to be called directly but only through the RCM() function
#'
#' @return A list with elements
#' \item{converged}{a vector of booleans of length k indicating if the algorithm
#'  converged for every dimension}
#' \item{rMat}{(if not constrained a nxk matrix with estimated row scores}
#' \item{cMat}{ a kxp matrix with estimated column scores}
#' \item{psis}{a vector of length k
#'  with estimates for the importance parameters psi}
#' \item{thetas}{a vector of length p with estimates for the overdispersion}
#' \item{rowRec}{(if not constrained) a n x k x maxItOut array with a record
#'  of all rMat estimates through the iterations}
#' \item{colRec}{a k x p x maxItOut array with a record of all cMat
#'  estimates through the iterations}
#' \item{psiRec}{a k x maxItOut array with a record of all psi estimates
#'  through the iterations}
#' \item{thetaRec}{ a matrix of dimension pxmaxItOut with estimates for
#'  the overdispersion along the way}
#' \item{iter}{ number of iterations}
#' \item{Xorig}{ (if confounders provided) the original fitting matrix}
#' \item{X}{ the trimmed matrix if confounders provided,
#'  otherwise the original one}
#' \item{fit}{ type of fit, either "RCM_NB" or "RCM_NB_constr"}
#' \item{lambdaRow}{(if not constrained)
#'  vector of Lagrange multipliers for the rows}
#' \item{lambdaCol}{ vector of Lagrange multipliers for the columns}
#' \item{rowWeights}{(if not constrained) the row weights used}
#' \item{colWeights}{ the column weights used}
#' \item{alpha}{(if constrained) the kxd matrix of environmental gradients}
#' \item{alphaRec}{(if constrained) the kxdxmaxItOut array of alpha estimates
#'  along the iterations}
#' \item{covariates}{(if constrained) the matrix of covariates}
#' \item{libSizes}{ a vector of length n with estimated library sizes}
#' \item{abunds}{ a vector of length p with estimated mean relative abundances}
#' \item{confounders}{(if provided) the confounder matrix}
#' \item{confParams}{ the parameters used to filter out the confounders}
#' \item{nonParamRespFun}{A list of the non parametric response functions}
#' \item{degree}{The degree of the alternative parametric fit}
#' \item{devFilt}{The deviance after filtering confounders}
#' \item{llFilt}{The likelihood of the model after filtering on confounders}
#' @export
#' @note Plotting is not supported for quadratic response functions
#' @examples
#' data(Zeller)
#' require(phyloseq)
#' tmpPhy = prune_taxa(taxa_names(Zeller)[seq_len(100)],
#' prune_samples(sample_names(Zeller)[seq_len(50)], Zeller))
#' mat = matrix(otu_table(tmpPhy), nsamples(tmpPhy), ntaxa(tmpPhy))
#' mat = mat[rowSums(mat)>0, colSums(mat)>0]
#' zellerRCM = RCM_NB(mat, k = 2)
#' #Needs to be called directly onto a matrix
RCM_NB = function(X, k, rowWeights = "uniform", colWeights = "marginal",
                  tol = 1e-3, maxItOut = 1000L, Psitol = 1e-3, verbose = FALSE,
                  global = "dbldog",
                  nleqslv.control = list(maxit = 500L, cndtol = 1-16),
                  jacMethod = "Broyden", dispFreq = 10L, convNorm = 2,
                  prior.df=10, marginEst = "MLE", confounders = NULL,
                  prevCutOff, minFraction = 0.1, covariates = NULL,
                  centMat = NULL,
                  responseFun = c("linear", "quadratic","dynamic",
                                  "nonparametric"),
                  record = FALSE, control.outer = list(trace=FALSE),
                  control.optim = list(), envGradEst = "LR", dfSpline = 3,
                  vgamMaxit = 100L,
                  degree = switch(responseFun[1], "nonparametric" = 3, NULL)){
  Xorig = NULL #An original matrix, not returned if no trimming occurs
  responseFun = responseFun[1]
  if(!responseFun %in% c("linear", "quadratic","dynamic","nonparametric")){
    stop("Unknown response function provided! See ?RCM_NB for details.")
  }

  if(!is.null(confounders$confounders)){
    #First and foremost: filter on confounders
      Xorig = X
      X = trimOnConfounders(X, confounders = confounders$confoundersTrim,
                            prevCutOff = prevCutOff, n=nrow(Xorig),
                            minFraction = minFraction)
    }

    n=NROW(X)
    p=NCOL(X)

    thetas = matrix(0,p, k+1+(!is.null(confounders$confounders)),
                    dimnames = list(colnames(X),
                                    c("Independence",
                                      if(!is.null(confounders$confounders))
                                        "Filtered" else NULL,
                                      paste0("Dim",seq_len(k)))))
    #Track the overdispersions,also the one associated to the independence model

    #Initialize margin parameters
    abunds = colSums(X)/sum(X)
    libSizes = rowSums(X)

    #Get the trended-dispersion estimates. Very insensitive to the offset
    #so only need to be calculated once per dimension
    trended.dispersion <- edgeR::estimateGLMTrendedDisp(y = t(X),
                                                        design = NULL,
                                                        method = "bin.loess",
                                                        offset = t(log(outer(
                                                          libSizes, abunds))),
                                                        weights = NULL)

    if(marginEst == "MLE"){
      logLibSizesMLE = log(libSizes)
      logAbundsMLE = log(abunds)
      initIter = 1
      if(verbose) cat("\nEstimating the independence model \n\n")

      while((initIter ==1) || ((initIter <= maxItOut) && (!convergenceInit))){
        logLibsOld = logLibSizesMLE
        logAbsOld = logAbundsMLE

        thetas[,"Independence"] = estDisp(X = X, rMat = as.matrix(rep(0,n)),
                                          cMat = t(as.matrix(rep(0,p))),
                                          muMarg=exp(outer(logLibSizesMLE,
                                                           logAbundsMLE, "+")),
                                          psis = 0, prior.df = prior.df,
                                        trended.dispersion = trended.dispersion)
        thetasMat = matrix(thetas[,"Independence"], n, p, byrow=TRUE)

        logLibSizesMLE = nleqslv(fn = dNBlibSizes, x = logLibSizesMLE,
                                 theta = thetasMat, X = X, reg=logAbundsMLE,
                                 global=global, control = nleqslv.control,
                                 jac=NBjacobianLibSizes, method=jacMethod)$x
        logAbundsMLE = nleqslv(fn = dNBabunds, x = logAbundsMLE,
                               theta = thetasMat, X = X, reg=logLibSizesMLE,
                               global=global, control = nleqslv.control,
                               jac=NBjacobianAbunds, method=jacMethod)$x
        initIter = initIter + 1

        convergenceInit = ((initIter <= maxItOut) &&
                             ((sum(abs(1-logLibSizesMLE/logLibsOld)^
                                     convNorm)/n)^(1/convNorm) < tol) &&
                             ((sum(abs(1-logAbundsMLE/logAbsOld)^convNorm)/p)^
                                (1/convNorm) < tol) )
      }
      #Converges very fast, even when dispersions are re-estimated.
      #For the library sizes there is a big difference,
      #for the abundances less so
      muMarg = exp(outer(logLibSizesMLE, logAbundsMLE, "+"))
      #The marginals to be used as expectation. These are augmented with the
      #previously estimated dimensions every time

    } else if(marginEst=="marginSums"){
      muMarg = outer(libSizes, abunds)
    } else{
      stop("No valid margin estimation paradigm provided! \n")
    }

    rowWeights = switch(paste(marginEst, rowWeights, sep = "_"),
                        "marginSums_marginal" = libSizes/sum(libSizes),
                        "MLE_marginal" = exp(logLibSizesMLE)/
                          sum(exp(logLibSizesMLE)),
                        rep.int(1/n,n) #For uniform weights
    )
    colWeights = switch(paste(marginEst, colWeights, sep = "_"),
                        "marginSums_marginal" = abunds,
                        "MLE_marginal" = exp(logAbundsMLE),
                        rep.int(1/p,p) #For uniform weights
    )

    nLambda = 2*k+k*(k-1)/2
    if(record){
      # Pre-allocate arrays to track iterations
      rowRec = array(0,dim=c(n,k, maxItOut))
      colRec = thetaRec = array(0,dim=c(k,p, maxItOut))
      psiRec = matrix(0, nrow=k,ncol=maxItOut)
    } else {
      rowRec = colRec = thetaRec = psiRec = NULL
    }
    convergence = rep(FALSE, k)
    iterOut = rep(1,k)
    if(!is.null(confounders$confounders)){
      ## Filter out the confounders by adding them to the intercept,
      #also adapt overdispersions
      filtObj = filterConfounders(muMarg = muMarg,
                                  confMat = confounders$confounders, p=p, X=X,
                                  thetas = thetas[,1],
                                  nleqslv.control = nleqslv.control, n=n,
                                  trended.dispersion = trended.dispersion)
      muMarg = muMarg * exp(confounders$confounders %*% filtObj$NB_params)
      thetas[,"Filtered"] = filtObj$thetas
      confParams = filtObj$NB_params
    } else {
      confParams = NULL
    }
    ## 1) Initialization
    svdX = svd(diag(1/libSizes) %*% (X-muMarg) %*% diag(1/colSums(X)))
    rMat = svdX$u[,seq_len(k), drop=FALSE]
    cMat = t(svdX$v[,seq_len(k), drop=FALSE])
    psis = svdX$d[seq_len(k)]

    lambdaRow =  rep.int(0,nLambda)
    lambdaCol =  rep.int(0,nLambda )
    #Center
    cMat = t(apply(cMat, 1, function(colS){
      colS-sum(colS*colWeights)/sum(colWeights)
    }))
    rMat = apply(rMat, 2, function(rowS){
      rowS-sum(rowS*rowWeights)/sum(rowWeights)
    })

    #Redistribute some weight to fit the constraints
    psis = c(psis *t(apply(cMat, 1, function(colS){
      sqrt(sum(colWeights * colS^2))
    })) * apply(rMat, 2, function(rowS){
      sqrt(sum(rowWeights * rowS^2))
    }))

    #Normalize
    cMat = t(apply(cMat, 1, function(colS){
      colS/sqrt(sum(colWeights * colS^2))
    }))

    rMat = apply(rMat, 2, function(rowS){
      rowS/sqrt(sum(rowWeights * rowS^2))
    })

  if(is.null(covariates)){
    #If no covariates provided, perform an unconstrained analysis
    for (KK in seq_len(k)){

      if(verbose) cat("Dimension" ,KK, "is being esimated \n")

      #Modify offset if needed
      if(KK>1){muMarg = muMarg * exp(rMat[,(KK-1), drop=FALSE] %*% (cMat[(KK-1),
                        , drop=FALSE]*psis[(KK-1)]))}
      idK = seq_k(KK) #prepare an index

      #Re-estimate the trended dispersions, once per dimensions
      trended.dispersion <- edgeR::estimateGLMTrendedDisp(y = t(X),
                                                          design = NULL,
                                                          method = "bin.loess",
                                                          offset = t(log(
                                                            muMarg)),
                                                          weights = NULL)

      JacR = matrix(0, nrow = n+KK+1, ncol = n+KK+1) #Prepare sparse Jacobians,
      #and prefill what we can
      JacR[seq_len(n), n+1] = rowWeights
      if(KK>1){
        JacR[seq_len(n),(n+3):(n+KK+1)] = rMat[,seq_len(KK-1), drop=FALSE]*
          rowWeights
      }
      #Symmetrize
      JacR = JacR + t(JacR)

      JacC = matrix(0, nrow = p+KK+1, ncol = p+KK+1)
      JacC[seq_len(p), p+1] = colWeights
      if(KK>1){
        JacC[seq_len(p),(p+3):(p+KK+1)] = t(cMat[seq_len(KK-1),, drop=FALSE])*
          colWeights
      }
      #Symmetrize
      JacC = JacC + t(JacC)

      while((iterOut[KK] ==1) || ((iterOut[KK] <= maxItOut) &&
                                  (!convergence[KK])))
      {

        if(verbose && iterOut[KK]%%1 == 0){
          cat("\n","Outer Iteration", iterOut[KK], "\n","\n")
          if(iterOut[KK]!=1){
            cat("Old psi-estimate: ", psisOld, "\n")
            cat("New psi-estimate: ", psis[KK], "\n")
          }
        }
        ## 2)a. Store old parameters
        psisOld = psis[KK]
        rMatOld = rMat[,KK]
        cMatOld = cMat[KK,]

        #Overdispersions (not at every iterations to speed things up,
        #the estimates do not change a lot anyway)
        if((iterOut[KK] %% dispFreq) ==0 || (iterOut[KK]==1)){
          if (verbose) cat(" Estimating overdispersions \n")
          thetas[,paste0("Dim",KK)] = estDisp(X = X,
                                              rMat = rMat[,KK,drop=FALSE],
                                              cMat = cMat[KK,,drop=FALSE],
                                              muMarg=muMarg, psis = psis[KK],
                                              prior.df = prior.df,
                                              trended.dispersion =
                                                trended.dispersion)
          thetasMat = matrix(thetas[,paste0("Dim",KK)], n, p, byrow=TRUE)
          #Make a matrix for numerical reasons,
          #it avoids excessive use of the t() function
          preFabMat = 1+X/thetasMat # Another matrix that can be pre-calculated
        }
        #Psis
        if (verbose) cat("\n Estimating psis \n")
        regPsis = outer(rMat[,KK] ,cMat[KK,])

        psiTry = try(abs(nleqslv(fn = dNBpsis, x = psis[KK], theta = thetasMat,
                                 X = X, reg=regPsis, muMarg=muMarg,
                                 global=global, control = nleqslv.control,
                                 jac=NBjacobianPsi, method=jacMethod,
                                 preFabMat = preFabMat)$x))
        if(inherits(psiTry,"try-error")){
          stop("Fit failed, likely due to numeric reasons.
               Consider more stringent filtering by increasing
               the prevCutOff parameter.\n")
        } else {psis[KK] = psiTry}

        #Column scores
        if (verbose) cat("\n Estimating column scores \n")
        regCol = rMat[,KK, drop=FALSE]*psis[KK]
        tmpCol = nleqslv(fn = dNBllcol, x = c(cMat[KK,], lambdaCol[idK]),
                         thetas = thetasMat, X = X, reg = regCol,
                         muMarg = muMarg, k = KK,  global = global,
                         control = nleqslv.control, n=n, p=p,
                         jac = NBjacobianCol, method = jacMethod,
                         colWeights = colWeights, nLambda = (KK+1),
                         cMatK = cMat[seq(1,(KK-1)),,drop=FALSE],
                         preFabMat = preFabMat, Jac = JacC)

        if(verbose) cat(ifelse(tmpCol$termcd==1, "Column scores converged \n",
                               "Column scores DID NOT converge \n"))
        cMat[KK,] = tmpCol$x[seq_len(p)]
        lambdaCol[idK] = tmpCol$x[p + seq_along(idK)]

        #Normalize (speeds up algorithm if previous step had not converged)
        cMat[KK,] = cMat[KK,] - sum(cMat[KK,] * colWeights)/sum(colWeights)
        cMat[KK,] = cMat[KK,]/sqrt(sum(colWeights * cMat[KK,]^2))
        #Row scores
        if (verbose) cat("\n Estimating row scores \n")
        regRow = cMat[KK,,drop=FALSE]*psis[KK]
        tmpRow = nleqslv(fn = dNBllrow, x = c(rMat[,KK], lambdaRow[idK]),
                         thetas=thetasMat, X = X, reg = regRow, muMarg = muMarg,
                         k = KK,  global = global, control = nleqslv.control,
                         n = n, p = p, jac = NBjacobianRow, method = jacMethod,
                         rowWeights = rowWeights, nLambda = (KK+1),
                         rMatK = rMat[,seq(1,(KK-1)), drop=FALSE],
                         preFabMat = preFabMat, Jac = JacR)

        if(verbose) cat(ifelse(tmpRow$termcd==1, "Row scores converged \n",
                               "Row scores DID NOT converge \n"))
        rMat[,KK] = tmpRow$x[seq_len(n)]
        lambdaRow[idK] = tmpRow$x[n + seq_along(idK)]

        #Normalize (speeds up algorithm if previous step had not converged)
        rMat[,KK] = rMat[,KK] - sum(rMat[,KK] * rowWeights)/sum(rowWeights)
        rMat[,KK] = rMat[,KK]/sqrt(sum(rowWeights * rMat[,KK]^2))

        if(record){
          #Store intermediate estimates
          rowRec[,KK, iterOut[KK]] = rMat[,KK]
          colRec[KK,, iterOut[KK]] = cMat[KK,]
          thetaRec [KK,, iterOut[KK]] = thetas[,paste0("Dim",KK)]
          psiRec[KK, iterOut[KK]] = psis[KK]
        }

        ## Change iterator
        iterOut[KK] = iterOut[KK] + 1

        ##Check convergence  (any numbered norm for row and column scores)
        convergence[KK] = ((iterOut[KK] <= maxItOut) &&
                             (abs(1-psis[KK]/psisOld) < Psitol) &&
                             #Infinity norm for the psis
                             ((mean(abs(1-rMat[,KK]/rMatOld)^convNorm))^
                                (1/convNorm) < tol) &&
                             ((mean(abs(1-cMat[KK,]/cMatOld)^convNorm))^
                                (1/convNorm) < tol) )
      } # END while-loop until convergence

    }# END for-loop over dimensions

    ## 3) Termination

    rownames(rMat) = rownames(X)
    colnames(cMat) = colnames(X)
    rownames(cMat) = colnames(rMat) = paste0("Dim",seq_len(k))

    returnList = list(rMat = rMat, cMat=cMat, rowRec = rowRec, colRec = colRec,
                      psiRec = psiRec, thetaRec = thetaRec, fit = "RCM_NB",
                      lambdaRow = lambdaRow, lambdaCol = lambdaCol)

  } else { #If covariates provided, do a constrained analysis
    d = ncol(covariates)
    CCA = vegan::cca(X = X, Y = covariates)$CCA
    #Constrained correspondence analysis for starting values
    if(sum(!colnames(covariates) %in% CCA$alias)<k) {
      k = sum(!colnames(covariates) %in% CCA$alias)
      warning(immediate. = TRUE, paste("Can only fit an ordination with",
                                       k,"dimensions with so few covariates!"))
    }
    alpha = matrix(0,d,k)
    alpha[!colnames(covariates) %in% CCA$alias,] = CCA$biplot[,seq_len(k)]
    #Leave the sum constraints for the factors alone for now,
    #may or may not speed up the algorithm
    alpha = t(t(alpha)-colMeans(alpha))
    alpha = t(t(alpha)/sqrt(colSums(alpha^2)))
    psis = CCA$eig[seq_len(k)]
    alphaRec = if(record){array(0, dim=c(d, k, maxItOut))} else {NULL}
    v = switch(responseFun, linear = 2, quadratic = 3, dynamic = 3, 1)
    #Number of parameters per taxon
    NB_params = array(0.1,dim=c(v,p,k))
    #Initiate parameters of the response function,
    #taxon-wise. No zeroes or trivial fit!
    #Improved starting values may be possible.
    NB_params = if(responseFun != "nonparametric") {
      vapply(seq_len(k),FUN.VALUE = matrix(0,v,p),
             function(x){x = NB_params[,,x, drop=FALSE]
             x/sqrt(rowSums(x^2))})}else NULL
    NB_params_noLab = if(responseFun != "nonparametric" && envGradEst == "LR") {
      matrix(0.1,v,k)} else NULL
    #Initiate parameters of the response function, ignoring taxon-labels
    if(responseFun == "nonparametric") {
      nonParamRespFun = lapply(seq_len(k),
                               function(x){list(taxonWise = lapply(integer(p),
      function(d){list(fit =list(coef=NULL))}), overall = NULL)})
    } else {nonParamRespFun =NULL}
    rowMat = NULL

    #Number of lambda parameters for centering
    nLambda1s = NROW(centMat)

    lambdasAlpha = rep(0, (k*(1+nLambda1s+(k-1)/2)))
    for (KK in seq_len(k)){
      if(verbose) cat("Dimension" ,KK, "is being esimated \n")

      #Modify offset if needed
      if(KK>1){
        muMarg = if(responseFun %in% c("linear","quadratic", "dynamic")){
          exp(getRowMat(responseFun = responseFun, sampleScore = covariates %*%
                          alpha[,KK-1, drop = FALSE],
                        NB_params = NB_params[,,KK-1])*psis[KK-1])*muMarg
        } else {
          exp(nonParamRespFun[[KK-1]]$rowMat) * muMarg
        }
      }

      idK = seq_k(KK)
      ## 2) Propagation
      while((iterOut[KK] ==1) || ((iterOut[KK] <= maxItOut) &&
                                  (!convergence[KK])))
      {

        if(verbose && iterOut[KK]%%1 == 0){
          cat("\n","Outer Iteration", iterOut[KK], "\n","\n")
          if(iterOut[KK]!=1 && responseFun!="nonparametric"){
            cat("Old psi-estimate: ", psisOld, "\n")
            cat("New psi-estimate: ", psis[KK], "\n")
          }
        }
        ## 2)a. Store old parameters to check for convergence
        psisOld = psis[KK]
        alphaOld = alpha[,KK]
        NBparamsOld = NB_params[,,KK]

        sampleScore = covariates %*% alpha[,KK]
        envRange = range(sampleScore)
        if(responseFun %in% c("linear","quadratic", "dynamic")){
          design = buildDesign(sampleScore, responseFun)
          rowMat = design %*% NB_params[,,KK]
        } else{
          rowMat = nonParamRespFun[[KK]]$rowMat
        }
        #Overdispersions (not at every iterations to speed things up,
        #doesn't change a lot anyway)
        if((iterOut[KK] %% dispFreq) == 0 || iterOut[KK] == 1){
          if (verbose) cat(" Estimating overdispersions \n")
          thetas[,paste0("Dim", KK)] = estDisp(X = X, muMarg = muMarg,
                                               psis = psis[KK],
                                               prior.df = prior.df,
                                               trended.dispersion =
                                               trended.dispersion,
                                               rowMat = rowMat)
          thetasMat = matrix(thetas[,paste0("Dim",KK)], n, p, byrow=TRUE)
          preFabMat = 1+X/thetasMat
        }

        if(responseFun %in% c("linear","quadratic", "dynamic")){
          #Psis
          if (verbose) cat("\n Estimating psis (k = ", KK, ") \n", sep="")
          psis[KK]  = abs(nleqslv(fn = dNBpsis, x = psis[KK], theta = thetasMat,
                                  X = X, reg = rowMat, muMarg = muMarg,
                                  global = global, control = nleqslv.control,
                                  jac = NBjacobianPsi, method = jacMethod,
                                  preFabMat = preFabMat)$x)

          if (verbose) cat("\n Estimating response function \n")
          NB_params[,,KK] = estNBparams(design = design,
                                        thetas = thetas[,paste0("Dim",KK)],
                                        muMarg = muMarg, psi = psis[KK], X = X,
                                        nleqslv.control = nleqslv.control,
                                        ncols = p, initParam = NB_params[,,KK],
                                        v = v, dynamic = responseFun=="dynamic",
                                        envRange = envRange)
          NB_params[,,KK] = NB_params[,,KK]/sqrt(rowSums(NB_params[,,KK]^2))

          if(envGradEst == "LR") {
            NB_params_noLab[, KK] = estNBparamsNoLab(design = design,
                                                     thetasMat = thetasMat,
                                                     muMarg = muMarg,
                                                     psi = psis[KK], X = X,
                                                     nleqslv.control =
                                                       nleqslv.control,
                                                     initParam =
                                                       NB_params_noLab[,KK],
                                                     v = v,
                                                     dynamic = responseFun ==
                                                       "dynamic",
                                                     envRange = envRange,
                                                     preFabMat = preFabMat,
                                                     n=n)}

          if (verbose) cat("\n Estimating environmental gradient \n")
          AlphaTmp = nleqslv(x = c(alpha[,KK],
                                   lambdasAlpha[seq_k(KK, nLambda1s)]),
                             fn = dLR_nb, jac = LR_nb_Jac,
                             X = X, CC = covariates, responseFun = responseFun,
                             cMat = cMat, psi = psis[KK],
                             NB_params = NB_params[,,KK],
                             NB_params_noLab = NB_params_noLab[, KK],
                             alphaK = alpha[, seq_len(KK-1), drop=FALSE],
                             k = KK, d = d, centMat = centMat,
                             nLambda = nLambda1s+KK, nLambda1s = nLambda1s,
                             thetaMat = thetasMat, muMarg = muMarg,
                             control = nleqslv.control, n=n, v=v, ncols = p,
                             preFabMat = preFabMat, envGradEst = envGradEst)$x
          alpha[,KK] = AlphaTmp[seq_len(d)]
          lambdasAlpha[seq_k(KK, nLambda1s)] = AlphaTmp[-seq_len(d)]

        } else {
          if (verbose) cat("\n Estimating response functions \n")
          nonParamRespFun[[KK]] = estNPresp(sampleScore = sampleScore,
                                            muMarg = muMarg, X = X, ncols = p,
                                            thetas = thetas[,paste0("Dim",KK)],
                                            n=n, coefInit =
                                              nonParamRespFun[[KK]]$taxonCoef,
                                            coefInitOverall =
                                            nonParamRespFun[[KK]]$overall$coef,
                                            vgamMaxit = vgamMaxit,
                                            dfSpline = dfSpline,
                                            verbose = verbose,
                                            degree = degree)

          if (verbose) cat("\n Estimating environmental gradient \n")
          AlphaTmp = constrOptim.nl(par = alpha[,KK], fn = LR_nb, gr = NULL,
                                    heq = heq_nb, heq.jac = heq_nb_jac,
                                    alphaK = alpha[, seq_len(KK-1), drop=FALSE],
                                    X=X, CC=covariates,
                                    responseFun = responseFun,
                                    muMarg = muMarg, d = d, ncols=p,
                                    control.outer = control.outer,
                                    control.optim = control.optim,
                                    k = KK, centMat = centMat, n=n,
                                    nonParamRespFun = nonParamRespFun[[KK]],
                                    thetaMat = thetasMat,
                                    envGradEst = envGradEst)
          alpha[,KK] = AlphaTmp$par
          lambdasAlpha[seq_k(KK, nLambda1s)] = AlphaTmp$lambda
        }

        #Store intermediate estimates
        if(record){
          alphaRec[,KK, iterOut[KK]] = alpha[,KK]
          thetaRec [KK,, iterOut[KK]] = thetas[,paste0("Dim",KK)]
          psiRec[KK, iterOut[KK]] = psis[KK]
        }
        ## Change iterator
        iterOut[KK] = iterOut[KK] + 1

        ##Check convergence  (any numbered norm for row and column scores)
        convergence[KK] = ((iterOut[KK] <= maxItOut) &&
                             (abs(1-psis[KK]/psisOld) < Psitol) &&
                             #Infinity norm for the psis
                             ((mean(abs(1-alpha[,KK]/alphaOld)^convNorm))^
                                (1/convNorm) < tol) &&
                             #Env gradient
                             if(responseFun=="nonparametric") TRUE else
                               (mean(abs(1-(NB_params[,,KK]/NBparamsOld)[
                                 NBparamsOld*NB_params[,,KK] != 0])^convNorm)^
                                  (1/convNorm) < tol)
                           #Parameters of the response function
        )
      } # END while-loop until convergence

    }# END for-loop over dimensions
    ## 3) Termination

    if(responseFun=="nonparametric"){#Post hoc calculate integrals and psis
      nonParamRespFun = lapply(seq_len(k), function(KK){
        samScore = covariates %*% alpha[,KK]
        nonPar = nonParamRespFun[[KK]]
        nonPar$intList = vapply(FUN.VALUE = numeric(1),
                                colnames(X), function(tax){
          getInt(coef = nonParamRespFun[[KK]]$taxonCoef[[tax]],
                 spline = nonParamRespFun[[KK]]$splineList[[tax]],
                 sampleScore = samScore)
        })
        return(nonPar)
      })
      psis = vapply(nonParamRespFun, function(x){sqrt(mean(x$rowMat^2))},
                    numeric(1))
      names(nonParamRespFun) = names(psis) = paste0("Dim",seq_len(k))
    }

    rownames(alpha) = colnames(covariates)
    colnames(alpha) = paste0("Dim",seq_len(k))

    returnList = list( fit = "RCM_NB_constr", lambdaCol = lambdaCol,
                       alpha = alpha, alphaRec = alphaRec,
                       covariates = covariates, NB_params = NB_params,
                       NB_params_noLab = NB_params_noLab,
                       responseFun = responseFun,
                       nonParamRespFun = nonParamRespFun,
                       envGradEst = if(is.null(covariates)) NULL else
                         envGradEst, lambdasAlpha = lambdasAlpha,
                       degree = degree)
  }
  if(!all(convergence)){
    warning(paste0("Algorithm did not converge for dimensions ",
                   paste(which(!convergence), collapse = ","),
                   "! Check for errors or consider changing tolerances
                   or number of iterations"))
  }
  return(
    c(returnList, list(converged = convergence, psis = psis, thetas = thetas,
                       psiRec = psiRec, thetaRec = thetaRec, iter = iterOut-1,
                       X = X, Xorig = Xorig, rowWeights = rowWeights,
                       colWeights = colWeights,
                       libSizes = switch(marginEst,
                                         "MLE" = exp(logLibSizesMLE),
                                         "marginSums" = libSizes),
                       abunds = switch(marginEst,
                                       "MLE" = exp(logAbundsMLE),
                                       "marginSums" = abunds),
                       confounders = confounders, confParams = confParams))
  )
}
