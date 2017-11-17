#' A function to fit the RC(M) model with the negative binomial distribution. Includes fitting of the independence model, filtering out the effect of confounders and fitting the RC(M) components in a constrained or an unconstrained way for any dimension k.
#'
#' @param X a nxp data matrix
#' @param k an scalar, number of dimensions in the RC(M) model
#' @param rowWeights a character string, either "uniform" or "marginal" row weights. Defaults to "uniform"
#' @param colWeights a character string, either "uniform" or "marginal" column weights. Defaults to "uniform "marginal"
#' @param tol a scalar, the relative convergende tolerance for the row scores and column scores parameters, defaults to 1e-3
#' @param maxItOut an integer, the maximum number of iteration in the outer loop, defaults to 500
#' @param Psitol a scalar, the relative convergence tolerance for the psi parameters, defaults to 1e-4
#' @param verbose a boolean, should information on iterations be printed? Defaults to TRUE
#' @param NBRCM a previously fitted NBRCM object, from which the lower dimensions can be extracted. Only useful if NBRCM$xk < k
#' @param global global strategy for solving non-linear systems, see ?nleqslv
#' @param nleqslv.control a list with control options, see nleqslv
#' @param method Method for solving non-linear equations, ?see nleqslv. Defaults to Broyden. The difference with the newton method is that the Jacobian is not recalculated at every iteration, thereby speeding up the algorithm
#' @param dispFreq an integer, how many iterations the algorithm should wait before reestimationg the dispersions. Defaults to 20
#' @param convNorm a scalar, the norm to use to determine convergence
#' @param prior.df an integer, see estDisp()
#' @param marginEst a character string, either "MLE" or "marginSums", indicating how the independence model should be estimated
#' @param confounders a list with
#' -confounders an nxg matrix with confounders
#' -confoundersFilt an nxh matrix with confounders for filtering, with all levels and without intercept
#' @param covariates an nxd matrix with covariates. If set to null an unconstrained analysis is carried out, otherwise a constrained one. Factors must have been converted to dummy variables already
#' @param centMat a fxd matrix containing the contrasts to center the categorical variables. f equals the number of continuous variables + the total number of levels of the categorical variables.
#' @param responseFun a character string, either "linear", "quadratic", "dynamic" or "non-parametric"
#' @param prevCutOff a scalar the minimum prevalence needed to retain a taxon before the the confounder filtering
#' @param minFraction a scalar, total taxon abundance should equal minFraction*n if it wants to be retained before the confounder filtering
#' @param responseFun a characters string indicating the shape of the response function
#' @param record A boolean, should intermediate parameter estimates be stored?
#' @param control.outer a list of control options for the outer loop constrOptim.nl function
#' @param control.optim a list of control options for the optim() function
#' @param envGradEst a character string, indicating how the environmental gradient should be fitted. "LR" using the likelihood-ratio criterion, or "ML" a full maximum likelihood solution
#' @param dfSpline a scalar, the number of degrees of freedom for the splines of the non-parametric response function, see ?VGAM::s()
#' @param vgamMaxit an integer, the maximum number of iteration in the vgam() function
#'
#' Not intended to be called directly but only through the RCM() function
#'
#' @return A list with elements
#' \item{converged}{a vector of booleans of length k indicating if the algorithm converged for every dimension}
#' \item{rMat}{(if not constrained a nxk matrix with estimated row scores}
#' \item{cMat}{ a kxp matrix with estimated column scores}
#' \item{psis}{ a vector of length k with estimates for the importance parameters psi}
#' \item{thetas}{ a vector of length p with estimates for the overdispersion}
#' \item{rowRec}{(if not constrained) a n x k x maxItOut array with a record of all rMat estimates through the iterations}
#' \item{colRec}{ a k x p x maxItOut array with a record of all cMat estimates through the iterations}
#' \item{psiRec}{. a k x maxItOut array with a record of all psi estimates through the iterations}
#' \item{thetaRec}{ a matrix of dimension pxmaxItOut with estimates for the overdispersion along the way}
#' \item{iter}{ number of iterations}
#' \item{Xorig}{ (if confounders provided) the original fitting matrix}
#' \item{X}{ the trimmed matrix if confounders provided, otherwise the original one}
#' \item{fit}{ type of fit, either "RCM_NB" or "RCM_NB_constr"}
#' \item{lambdaRow}{(if not constrained) vector of Lagrange multipliers for the rows}
#' \item{lambdaCol}{ vector of Lagrange multipliers for the columns}
#' \item{rowWeights}{(if not constrained) the row weights used}
#' \item{colWeights}{ the column weights used}
#' \item{alpha}{(if constrained) the kxd matrix of environmental gradients}
#' \item{alphaRec}{(if constrained) the kxdxmaxItOut array of alpha estimates along the iterations}
#' \item{covariates}{(if constrained) the matrix of covariates}
#' \item{libSizes}{ a vector of length n with estimated library sizes}
#' \item{abunds}{ a vector of length p with estimated mean relative abundances}
#' \item{confounders}{(if provided) the confounder matrix}
#' \item{confParams}{ the parameters used to filter out the confounders}
#' \item{nonParamRespFun}{A list of the non parametric response functions}
RCM_NB = function(X, k, rowWeights = "uniform", colWeights = "marginal", tol = 1e-3, maxItOut = 2000L, Psitol = 1e-3, verbose = FALSE, NBRCM = NULL, global = "dbldog", nleqslv.control=list(maxit = 500L, cndtol = 1-16), jacMethod = "Broyden", dispFrec = 20L, convNorm = 2, prior.df=10, marginEst = "MLE", confounders = NULL, prevCutOff = 2.5e-2, minFraction = 0.1, covariates = NULL, centMat = NULL, responseFun = c("linear", "quadratic","dynamic","nonparametric"), record = FALSE, control.outer = list(trace=FALSE), control.optim = list(), envGradEst = "LR", dfSpline = 4, vgamMaxit = 100L){

  Xorig = NULL #An original matrix, not returned if no trimming occurs
  responseFun = responseFun[1]

  #If previous fit provided with higher or equal dimension, stop here
  if((!is.null(NBRCM)) ){
    if((is.null(covariates) & NBRCM$fit != "RCM_NB") | (!is.null(covariates) & NBRCM$fit != "RCM_NB_constr")){
      stop("Fit provided is not of same type as the one requested! \n Make sure you used exactly the same trimming parameters, otherwise start a new fit. \n")
    }else if((k <= NBRCM$k)) {
      stop("Fit provided is already of the required dimension or higher! \n")
    } else{
      #Read in all starting values
      rMat = NBRCM$rMat
      cMat = NBRCM$cMat
      psis = NBRCM$psis
      Kprev = NBRCM$k
      if(!all(X == NBRCM$X)){ stop("Different count matrix provided from original fit! \n")
      } else {}
      n=NROW(X)
      p=NCOL(X)
      thetas = NBRCM$thetas
      thetasMat = matrix(thetas[,Kprev+1], byrow = TRUE, n, p)
      rowWeights = NBRCM$rowWeights
      colWeights = NBRCM$colWeights
      libSizesMLE = NBRCM$libSizes
      logLibSizesMLE = log(libSizesMLE)
      abundsMLE = NBRCM$abunds
      lambdaRow = NBRCM$lambdaRow
      lambdaCol = NBRCM$lambdaCol
      muMarg = outer(libSizesMLE, abundsMLE)
      logAbundsMLE = log(abundsMLE)
      if(record){
        #Extend the record matrices, force same number of outer iterations as previous fit
        if(dim(NBRCM$colRec)[3]!=maxItOut){
          maxItOut = dim(NBRCM$colRec)[3]
          warning("Using same number of outer iterations as old fit! \n")
        }
      }
      newK = (Kprev+1):k
      if(is.null(covariates)){
        if(record){
          rowRec = array(0, dim = c(n,k,maxItOut))
          rowRec[,1:Kprev,] = NBRCM$rowRec
        } else {rowRec =  NULL}
        lambdaRow[1:(Kprev*(2+(Kprev-1)/2))] = NBRCM$lambdaRow
        svdX = svd(diag(1/rowSums(X)) %*% (X-muMarg) %*% diag(1/colSums(X)))
        muMarg = muMarg * exp(rMat %*% (cMat*psis))
        rMat = cbind(rMat, svdX$u[,newK, drop=FALSE])
        cMat = rbind(cMat, t(svdX$v[,newK, drop=FALSE]))
        psis = c(psis, svdX$d[newK])
        thetas = cbind(thetas, matrix(0,p,k-Kprev))
        #Center
        if(length(newK)==1){
          cMat[newK,] = c(apply(cMat[newK,, drop=FALSE], 1, function(colS){
            colS-sum(colS*colWeights)/sum(colWeights)
          }))
          rMat[,newK] = c(apply(rMat[,newK, drop=FALSE], 2, function(rowS){
            rowS-sum(rowS*rowWeights)/sum(rowWeights)
          }))

          #Redistribute some weight to fit the constraints
          psis[newK] = c(psis[newK] *t(apply(cMat[newK,, drop=FALSE], 1, function(colS){
            sqrt(sum(colWeights * colS^2))
          })) * apply(rMat[,newK, drop=FALSE], 2, function(rowS){
            sqrt(sum(rowWeights * rowS^2))
          }))

          #Normalize
          cMat[newK,] = c(apply(cMat[newK,, drop=FALSE], 1, function(colS){
            colS/sqrt(sum(colWeights * colS^2))
          }))

          rMat[,newK] = c(apply(rMat[,newK, drop=FALSE], 2, function(rowS){
            rowS/sqrt(sum(rowWeights * rowS^2))
          }))
        } else {
        cMat[newK,, drop=FALSE] = t(apply(cMat[newK,, drop=FALSE], 1, function(colS){
          colS-sum(colS*colWeights)/sum(colWeights)
        }))
        rMat[,newK, drop=FALSE] = apply(rMat[,newK, drop=FALSE], 2, function(rowS){
          rowS-sum(rowS*rowWeights)/sum(rowWeights)
        })

        #Redistribute some weight to fit the constraints
        psis[newK] = c(psis[newK] *t(apply(cMat[newK,, drop=FALSE], 1, function(colS){
          sqrt(sum(colWeights * colS^2))
        })) * apply(rMat[,newK, drop=FALSE], 2, function(rowS){
          sqrt(sum(rowWeights * rowS^2))
        }))

        #Normalize
        cMat[newK,, drop=FALSE] = t(apply(cMat[newK,, drop=FALSE], 1, function(colS){
          colS/sqrt(sum(colWeights * colS^2))
        }))

        rMat[,newK, drop=FALSE] = apply(rMat[,newK, drop=FALSE], 2, function(rowS){
          rowS/sqrt(sum(rowWeights * rowS^2))
        })
        }
      }
      if(record){
        colRec = array(0, dim = c(k,p,maxItOut))
        colRec[1:Kprev,,] = NBRCM$colRec
        thetaRec = array(0, dim = c(k,p,maxItOut))
        thetaRec[1:Kprev,,] = NBRCM$thetaRec
        psiRec = array(0, dim = c(k,maxItOut))
        psiRec[1:Kprev,] = NBRCM$psiRec
      } else {colRec = psiRec = thetaRec = NULL}
      lambdaCol = lambdaRow = rep(0, k*(2+(k-1)/2))
      lambdaCol[1:(Kprev*(2+(Kprev-1)/2))] = NBRCM$lambdaCol
      convergence = c(NBRCM$converged, rep(FALSE, k-Kprev))
      iterOut = c(NBRCM$iter, rep(1,k-Kprev))
      trended.dispersion <- estimateGLMTrendedDisp(y = t(X), design = NULL, method = "bin.loess", offset = t(outer(logLibSizesMLE, logAbundsMLE, FUN = "+")), weights = NULL)
      confParams = NBRCM$confParams
    }
    #Otherwise start the fit from scratch
  } else {
    if(!is.null(confounders[[1]])){ #First and foremost: filter on confounders
      Xorig = X
      X = trimOnConfounders(X, confounders = confounders$confoundersTrim, prevCutOff = prevCutOff, n=nrow(Xorig), minFraction = minFraction)
    }

    n=NROW(X)
    p=NCOL(X)

    thetas = matrix(0,p,k+1) #Track the overdispersions, also the one associated to the independence model

    #Initialize some parameters
    abunds = colSums(X)/sum(X)
    libSizes = rowSums(X)

    #Get the trended-dispersion estimates. Very insensitive to the offset so only need to be calculated once
    trended.dispersion <- estimateGLMTrendedDisp(y = t(X), design = NULL, method = "bin.loess", offset = t(log(outer(libSizes, abunds))), weights = NULL)

    if(marginEst == "MLE"){

      logLibSizesMLE = log(libSizes)
      logAbundsMLE = log(abunds)

      initIter = 1

      if(verbose) cat("Estimating the independence model \n")

      while((initIter ==1) || ((initIter <= maxItOut) && (!convergenceInit))){
        logLibsOld = logLibSizesMLE
        logAbsOld = logAbundsMLE

        thetas[,1] = estDisp(X = X, rMat = as.matrix(rep(0,n)), cMat = t(as.matrix(rep(0,p))),  muMarg=exp(outer(logLibSizesMLE, logAbundsMLE, "+")), psis = 0, prior.df = prior.df, trended.dispersion = trended.dispersion)
        thetasMat = matrix(thetas[,1], n, p, byrow=TRUE)

        logLibSizesMLE = nleqslv(fn = dNBlibSizes, x = logLibSizesMLE, theta = thetasMat, X = X, reg=logAbundsMLE, global=global, control = nleqslv.control, jac=NBjacobianLibSizes, method=jacMethod)$x
        logAbundsMLE = nleqslv(fn = dNBabunds, x = logAbundsMLE, theta = thetasMat, X = X, reg=logLibSizesMLE, global=global, control = nleqslv.control, jac=NBjacobianAbunds, method=jacMethod)$x
        initIter = initIter + 1

        convergenceInit = ((initIter <= maxItOut) &&
                             ((sum(abs(1-logLibSizesMLE/logLibsOld)^convNorm)/n)^(1/convNorm) < tol) &&
                             ((sum(abs(1-logAbundsMLE/logAbsOld)^convNorm)/p)^(1/convNorm) < tol) )
      }
      #Converges very fast, even when dispersions are re-estimated. For the library sizes there is a big difference, for the abundances less so
      muMarg = exp(outer(logLibSizesMLE, logAbundsMLE, "+")) #The marginals to be used as expectation. These are augmented with the previously estimated dimensions every time

    } else if(marginEst=="marginSums"){
      muMarg = outer(libSizes, abunds)
    } else{
      stop("No valid margin estimation paradigm provided! \n")
    }

    rowWeights = switch(paste(marginEst, rowWeights, sep = "_"),
                        "marginSums_marginal" = libSizes/sum(libSizes),
                        "MLE_marginal" = exp(logLibSizesMLE)/sum(exp(logLibSizesMLE)),
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
    if(!is.null(confounders[[1]])){
      ## Filter out the confounders by adding them to the intercept, also adapt overdispersions
      filtObj = filterConfounders(muMarg = muMarg, confMat = confounders$confounders, p=p, X=X, thetas = thetas[,1], nleqslv.control = nleqslv.control, n=n, trended.dispersion = trended.dispersion)
      muMarg = filtObj$muMarg
      thetas[,1] = filtObj$thetas
      confParams = filtObj$NB_params
    } else {
      confParams=NULL
    }
    ## 1) Initialization
    svdX = svd(diag(1/libSizes) %*% (X-muMarg) %*% diag(1/colSums(X)))
    rMat = svdX$u[,1:k, drop=FALSE]
    cMat = t(svdX$v[,1:k, drop=FALSE])
    psis = svdX$d[1:k]

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
  } # END if-else: no previous fit provided

  if(is.null(covariates)){ #If no covariates provided, perform an unconstrained analysis

    minK = ifelse(is.null(NBRCM),1,Kprev+1)
    for (KK in minK:k){

      if(verbose) cat("Dimension" ,KK, "is being esimated \n")

      #Modify offset if needed
      if(KK>1 && (if(!is.null(NBRCM)) {Kprev != (KK-1)} else {TRUE})){muMarg = muMarg * exp(rMat[,(KK-1), drop=FALSE] %*% (cMat[(KK-1),, drop=FALSE]*psis[(KK-1)]))}
      idK = seq_k(KK) #prepare an index

      JacR = matrix(0, nrow = n+KK+1, ncol = n+KK+1) #Prepare sparse Jacobians, and prefill what we can
      JacR[1:n, n+1] = rowWeights
      if(KK>1){
        JacR[1:n,(n+3):(n+KK+1)] = rMat[,1:(KK-1), drop=FALSE]*rowWeights
      }
      #Symmetrize
      JacR = JacR + t(JacR)

      JacC = matrix(0, nrow = p+KK+1, ncol = p+KK+1)
      JacC[1:p, p+1] = colWeights
      if(KK>1){
        JacC[1:p,(p+3):(p+KK+1)] = t(cMat[1:(KK-1),, drop=FALSE])*colWeights
      }
      #Symmetrize
      JacC = JacC + t(JacC)

      while((iterOut[KK] ==1) || ((iterOut[KK] <= maxItOut) && (!convergence[KK])))
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

        #Overdispersions (not at every iterations to speed things up, the estimates do not change a lot anyway)
        if((iterOut[KK] %% dispFrec) ==0 || (iterOut[KK]==1)){
          if (verbose) cat(" Estimating overdispersions \n")
          thetas[,KK+1] = estDisp(X = X, rMat = rMat[,KK,drop=FALSE], cMat = cMat[KK,,drop=FALSE], muMarg=muMarg, psis = psis[KK], prior.df = prior.df, trended.dispersion = trended.dispersion)
          thetasMat = matrix(thetas[,KK+1], n, p, byrow=TRUE) #Make a matrix for numerical reasons, it avoids excessive use of the t() function
          preFabMat = 1+X/thetasMat # Another matrix that can be pre-calculated
        }
        #Psis
        if (verbose) cat("\n Estimating psis \n")
        regPsis = outer(rMat[,KK] ,cMat[KK,])

        psis[KK]  = abs(nleqslv(fn = dNBpsis, x = psis[KK], theta = thetasMat, X = X, reg=regPsis, muMarg=muMarg, global=global, control = nleqslv.control, jac=NBjacobianPsi, method=jacMethod, preFabMat = preFabMat)$x)
        #Column scores
        if (verbose) cat("\n Estimating column scores \n")
        regCol = rMat[,KK, drop=FALSE]*psis[KK]
        tmpCol = nleqslv(fn = dNBllcol, x = c(cMat[KK,], lambdaCol[idK]), thetas = thetasMat, X = X, reg = regCol, muMarg = muMarg, k = KK,  global = global, control = nleqslv.control, n=n, p=p, jac = NBjacobianCol, method = jacMethod, colWeights = colWeights, nLambda = (KK+1), cMatK = cMat[1:(KK-1),,drop=FALSE], preFabMat = preFabMat, Jac = JacC)

        if(verbose) cat(ifelse(tmpCol$termcd==1, "Column scores converged \n", "Column scores DID NOT converge \n"))
        cMat[KK,] = tmpCol$x[1:p]
        lambdaCol[idK] = tmpCol$x[p + seq_along(idK)]

        #Normalize (speeds up algorithm if previous step had not converged)
        cMat[KK,] = cMat[KK,] - sum(cMat[KK,] * colWeights)/sum(colWeights)
        cMat[KK,] = cMat[KK,]/sqrt(sum(colWeights * cMat[KK,]^2))
        #Row scores
        if (verbose) cat("\n Estimating row scores \n")
        regRow = cMat[KK,,drop=FALSE]*psis[KK]
        tmpRow = nleqslv(fn = dNBllrow, x = c(rMat[,KK], lambdaRow[idK]), thetas=thetasMat, X = X, reg = regRow, muMarg = muMarg, k = KK,  global = global, control = nleqslv.control, n = n, p = p, jac = NBjacobianRow, method = jacMethod, rowWeights = rowWeights, nLambda = (KK+1), rMatK = rMat[,1:(KK-1), drop=FALSE], preFabMat = preFabMat, Jac = JacR)

        if(verbose) cat(ifelse(tmpRow$termcd==1, "Row scores converged \n", "Row scores DID NOT converge \n"))
        rMat[,KK] = tmpRow$x[1:n]
        lambdaRow[idK] = tmpRow$x[n + seq_along(idK)]

        #Normalize (speeds up algorithm if previous step had not converged)
        rMat[,KK] = rMat[,KK] - sum(rMat[,KK] * rowWeights)/sum(rowWeights)
        rMat[,KK] = rMat[,KK]/sqrt(sum(rowWeights * rMat[,KK]^2))

        if(record){
          #Store intermediate estimates
          rowRec[,KK, iterOut[KK]] = rMat[,KK]
          colRec[KK,, iterOut[KK]] = cMat[KK,]
          thetaRec [KK,, iterOut[KK]] = thetas[,KK+1]
          psiRec[KK, iterOut[KK]] = psis[KK]
        }

        ## Change iterator
        iterOut[KK] = iterOut[KK] + 1

        ##Check convergence  (any numbered norm for row and column scores)
        convergence[KK] = ((iterOut[KK] <= maxItOut) &&
                             (all(abs(1-psis[KK]/psisOld) < Psitol)) && #Infinity norm for the psis
                             ((sum(abs(1-rMat[,KK]/rMatOld)^convNorm)/n)^(1/convNorm) < tol) &&
                             ((sum(abs(1-cMat[KK,]/cMatOld)^convNorm)/p)^(1/convNorm) < tol) )
      } # END while-loop until convergence

    }# END for-loop over dimensions

    ## 3) Termination

    rownames(rMat) = rownames(X)
    colnames(cMat) = colnames(X)
    rownames(cMat) = colnames(rMat) = paste0("Dim",1:k)

    returnList = list(converged = convergence, rMat = rMat, cMat=cMat, psis = psis, thetas = thetas,  rowRec = rowRec, colRec = colRec, psiRec = psiRec, thetaRec = thetaRec, iter = iterOut-1, X=X, Xorig = Xorig, fit = "RCM_NB", lambdaRow = lambdaRow, lambdaCol = lambdaCol, rowWeights = rowWeights, colWeights = colWeights,
                      libSizes = switch(marginEst, "MLE" = exp(logLibSizesMLE), "marginSums" = libSizes), abunds = switch(marginEst, "MLE" = exp(logAbundsMLE), "marginSums" = abunds),
                      confounders = confounders, confParams = confParams)

  } else { #If covariates provided, do a constrained analysis

    d = ncol(covariates)
    require(vegan)
    CCA = cca(X = X, Y = covariates)$CCA #Constrained correspondence analysis for starting values
    alpha = matrix(0,d,k)
    alpha[!colnames(covariates) %in% CCA$alias,] = CCA$biplot[,1:k] #Leave the sum constraints for the factors alone for now, may or may not speed up the algorithm
    alpha = t(t(alpha)-colMeans(alpha))
    alpha = t(t(alpha)/sqrt(colSums(alpha^2)))
    if(!is.null(NBRCM)){
      newK = (Kprev+1):k
      psis = c(psis, CCA$eig[newK])
    } else {
      psis = CCA$eig[1:k]
    }
    alphaRec = if(record){array(0, dim=c(d, k, maxItOut))} else {NULL}
    v = switch(responseFun, linear = 2, quadratic = 3, dynamic = 3, 1) #Number of parameters per taxon
    NB_params = array(0.1,dim=c(v,p,k)) #Initiate parameters of the response function, taxon-wise. No zeroes or trivial fit! Improved starting values may be possible.
    NB_params = if(responseFun != "nonparametric") vapply(seq_len(k),FUN.VALUE = matrix(0,v,p), function(x){x = NB_params[,,x, drop=FALSE];x/sqrt(rowSums(x^2))}) else NULL
    NB_params_noLab = if(responseFun != "nonparametric" && envGradEst == "LR") matrix(0.1,v,k) else NULL #Initiate parameters of the response function, ignoring taxon-labels
    if(responseFun == "nonparametric") {nonParamRespFun =lapply(1:k,function(x){list(taxonWiseFitted = muMarg, taxonCoef = NULL, overallFitted = matrix(0,n,p), overallCoef = NULL)}); names(nonParamRespFun) = paste0("Dim",1:k)} else {nonParamRespFun =NULL}
    rowMat = NULL

    if(!is.null(NBRCM)){ #If fit provided, replace lower dimension starting values
      alpha[,1:Kprev] = NBRCM$alpha
      if(record) alphaRec[,1:Kprev,] = NBRCM$alphaRec
      NB_params[,,1:Kprev] = NBRCM$NB_params
      if(envGradEst == "LR") {NB_params_noLab[,1:Kprev] = NBRCM$NB_params_noLab}
    }

    #Number of lambda parameters for centering
    nLambda1s = NROW(centMat)

    minK = ifelse(is.null(NBRCM),1,Kprev+1) #Next dimension to fit
    lambdasAlpha = c(NBRCM$lambdasAlpha, rep(0, (k*(1+nLambda1s+(k-1)/2) - (minK-1)*(1+nLambda1s+(minK-2)/2))))
    for (KK in minK:k){
      if(verbose) cat("Dimension" ,KK, "is being esimated \n")

      #Modify offset if needed
      if(KK>1){
        muMarg = if(responseFun %in% c("linear","quadratic", "dynamic")){
          exp(getRowMat(responseFun = responseFun, sampleScore = covariates %*% alpha[,KK-1, drop = FALSE], NB_params = NB_params[,,KK-1])*psis[KK-1])*muMarg
        } else {
          nonParamRespFun[[KK-1]]$taxonWiseFitted
        }
      }

      idK = seq_k(KK)
      ## 2) Propagation
      while((iterOut[KK] ==1) || ((iterOut[KK] <= maxItOut) && (!convergence[KK])))
      {

        if(verbose && iterOut[KK]%%1 == 0){
          cat("\n","Outer Iteration", iterOut[KK], "\n","\n")
          if(iterOut[KK]!=1){
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
          rowMat = NULL
        }
        #Overdispersions (not at every iterations to speed things up, doesn't change a lot anyway)
        if((iterOut[KK] %% dispFrec) == 0 || iterOut[KK] == 1){
          if (verbose) cat(" Estimating overdispersions \n")
          thetas[,KK+1] = estDisp(X = X, muMarg = if(responseFun == "nonparametric") nonParamRespFun[[KK]]$taxonWiseFitted else muMarg, psis = psis[KK], prior.df = prior.df, trended.dispersion = trended.dispersion, rowMat = rowMat)
          thetasMat = matrix(thetas[,KK+1], n, p, byrow=TRUE)
          preFabMat = 1+X/thetasMat
        }

        if(responseFun %in% c("linear","quadratic", "dynamic")){
        #Psis
        if (verbose) cat("\n Estimating psis (k = ", KK, ") \n", sep="")
        psis[KK]  = abs(nleqslv(fn = dNBpsis, x = psis[KK], theta = thetasMat , X = X, reg = rowMat, muMarg = muMarg, global = global, control = nleqslv.control, jac = NBjacobianPsi, method = jacMethod, preFabMat = preFabMat)$x)

        if (verbose) cat("\n Estimating response function \n")
        NB_params[,,KK] = estNBparams(design = design, thetas = thetas[,KK+1], muMarg = muMarg, psi = psis[KK], X = X, nleqslv.control = nleqslv.control, ncols = p, initParam = NB_params[,,KK], v = v, dynamic = responseFun=="dynamic", envRange = envRange)
        #NB_params[,,KK] = NB_params[,,KK]/sqrt(rowSums(rowMultiply(NB_params[,,KK]^2,colWeights))) #The post-hoc normalization is much more efficient, since the equations are easier to solve. Crucially, we do not need orthogonality with other dimensions, which makes this approach feasible
        NB_params[,,KK] = NB_params[,,KK]/sqrt(rowSums(NB_params[,,KK]^2))

        if(envGradEst == "LR") {NB_params_noLab[, KK] = estNBparamsNoLab(design = design, thetasMat = thetasMat, muMarg = muMarg, psi = psis[KK], X = X, nleqslv.control = nleqslv.control, initParam = NB_params_noLab[,KK], v = v, dynamic = responseFun == "dynamic", envRange = envRange, preFabMat = preFabMat, n=n)}

        if (verbose) cat("\n Estimating environmental gradient \n")
      AlphaTmp = nleqslv(x = c(alpha[,KK],lambdasAlpha[seq_k(KK, nLambda1s)]), fn = dLR_nb, jac = LR_nb_Jac, X = X, CC = covariates, responseFun = responseFun, cMat = cMat, psi = psis[KK], NB_params = NB_params[,,KK], NB_params_noLab = NB_params_noLab[, KK], alphaK = alpha[, seq_len(KK-1), drop=FALSE], k = KK, d = d, centMat = centMat, nLambda = nLambda1s+KK, nLambda1s = nLambda1s, thetaMat = thetasMat, muMarg = muMarg, control = nleqslv.control, n=n, v=v, ncols = p, preFabMat = preFabMat, envGradEst = envGradEst)
      if(AlphaTmp$termcd!=1) warning("Environmental gradient could not be fitted!")
      alpha[,KK] = AlphaTmp$x[seq_len(d)]
      lambdasAlpha[seq_k(KK, nLambda1s)] = AlphaTmp$x[-seq_len(d)]

        } else {
          if (verbose) cat("\n Estimating response functions \n")
          nonParamRespFun[[KK]] = estNPresp(sampleScore = sampleScore, muMarg = muMarg, X = X, ncols = p, thetas = thetas[,KK+1], n=n, coefInit = nonParamRespFun[[KK]]$taxonCoef, coefInitOverall = nonParamRespFun[[KK]]$overallCoef, vgamMaxit = vgamMaxit, degreeSpline = dfSpline, colWeights = colWeights, verbose = verbose)
          if (verbose) cat("\n Estimating environmental gradient \n")
          AlphaTmp = constrOptim.nl(par = alpha[,KK], fn = LR_nb, gr = NULL, heq = heq_nb, heq.jac = heq_nb_jac, alphaK = alpha[, seq_len(KK-1), drop=FALSE], X=X, CC=covariates, responseFun = responseFun, muMarg = muMarg, d = d, ncols=p, control.outer = control.outer, control.optim = control.optim, nleqslv.control = nleqslv.control, k = KK, centMat = centMat, n=n, nonParamRespFun = nonParamRespFun[[KK]], thetaMat = thetasMat, envGradEst = envGradEst)
          alpha[,KK] = AlphaTmp$par
          lambdasAlpha[seq_k(KK, nLambda1s)] = AlphaTmp$lambda
          psis[KK] = nonParamRespFun[[KK]]$psi #Get the psis based on the integrals
        }

        #Store intermediate estimates
        if(record){
          alphaRec[,KK, iterOut[KK]] = alpha[,KK]
          thetaRec [KK,, iterOut[KK]] = thetas[,KK+1]
          psiRec[KK, iterOut[KK]] = psis[KK]
        }
        ## Change iterator
        iterOut[KK] = iterOut[KK] + 1

        ##Check convergence  (any numbered norm for row and column scores)
        convergence[KK] = ((iterOut[KK] <= maxItOut) &&
                             (all(abs(1-psis[KK]/psisOld) < Psitol)) && #Infinity norm for the psis
                             ((sum(abs(1-alpha[,KK]/alphaOld)^convNorm)/p)^(1/convNorm) < tol) && #Env gradient
                             if(responseFun=="nonparametric") TRUE else (mean(abs(1-(NB_params[,,KK]/NBparamsOld)[NBparamsOld*NB_params[,,KK] != 0])^convNorm)^(1/convNorm) < tol) #Parameters of the response function
        )
      } # END while-loop until convergence

    }# END for-loop over dimensions
    ## 3) Termination

    rownames(alpha) = colnames(covariates)
    colnames(cMat) = colnames(X)
    rownames(cMat) = colnames(alpha) = paste0("Dim",1:k)

    returnList = list(converged = convergence, psis = psis, thetas = thetas, psiRec = psiRec, thetaRec = thetaRec, iter = iterOut-1, X = X, Xorig = Xorig, fit = "RCM_NB_constr", lambdaCol = lambdaCol, rowWeights = rowWeights, colWeights = colWeights, alpha = alpha, alphaRec = alphaRec, covariates = covariates, NB_params = NB_params, NB_params_noLab = NB_params_noLab, libSizes = switch(marginEst, "MLE" = exp(logLibSizesMLE), "marginSums" = libSizes), abunds = switch(marginEst, "MLE" = exp(logAbundsMLE), "marginSums" = abunds), confounders = confounders, confParams = confParams, responseFun = responseFun, nonParamRespFun = nonParamRespFun, envGradEst = if(is.null(covariates)) NULL else envGradEst, lambdasAlpha = lambdasAlpha)
  }
  if(!all(convergence)){
    warning("Algorithm did not converge for all dimensions! Check for errors or consider changing tolerances or number of iterations")
  }
  return(returnList)
}