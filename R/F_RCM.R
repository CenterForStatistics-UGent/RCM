#'A simple wrapper function for the RC(M) object for phyloseq objects and matrices, for all posible distributions it passes all argument onto the correct function.
#'
#' @param dat: an nxp count matrix or a phyloseq object with an otu_table slot
#' @param k: an integer, the number of dimensions of the RCM solution
#' @param round : a boolean, whether to round to nearest integer. Defaults to FALSE.
#' @param distribution : a character string, the error distribution of the RC(M) model. Defaults to NB
#' @param prevCutoff : a scalar, the prevalance cutoff for the trimming. Defaults to 2.5e-2
#' @param minFraction : a scalar, each taxon's total abundance should equal at leat the number of samples n times minFraction, otherwise it is trimmed. Defaults to 10%
#' @param rowWeights, colWeights : character strings, the weighting procedures for the normalization of row and column scores. Defaults to "uniform" and "marginal" respectively
#' @param covariates : In case "dat" is a phyloseq object, the names of the sample variables to be used as covariates in the constrained analysis, or "all" to indicate all variables to be used. In case dat is a matrix, a nxf matrix or dataframe of covariates. Character variables will be converted to factors, with a warning. Defaults to NULL, in which case an unconstrained analysis is carried out.
#' @param confounders : In case "dat" is a phyloseq object, the names of the sample variables to be used as confounders to be filtered out. In case dat is a matrix, a nxf matrix or dataframe of confounders Character variables will be converted to factors, with a warning. Defaults to NULL, in which case no filtering occurs.
#' @param prevFit An object with a previous fit, normally from a lower dimension, that should be extended.
#'
#'Trim on prevalence and total abundance to avoid instability of the algorithm. We cannot conlcude much anyway on lowly abundant taxa
#'
#'@export
RCM = function(dat, k, round=FALSE, distribution= "NB", prevCutOff = 0.025, minFraction = 0.1, rowWeights = "uniform", colWeights = "marginal", covariates = NULL, confounders = NULL, prevFit = NULL, ...){

  classDat = class(dat) #The class provided

  ##The count data##
  if (classDat=="matrix"){
    X = dat
  }else  if(classDat=="phyloseq"){
    X = if (taxa_are_rows(dat)) t(otu_table(dat)@.Data) else otu_table(dat)@.Data
  } else {stop("Please provide a matrix or a phyloseq object! \n")}
  p=ncol(X); n=nrow(X)

  if(is.null(colnames(X))){colnames(X)=1:ncol(X)}
  if(is.null(rownames(X))){rownames(X)=1:nrow(X)}
  colNames = colnames(X); rowNames =rownames(X)
  # names(rowWeightsNum) =rowNames; names(colWeightsNum) = colNames
  if(round) {X=round(X, 0) }#Round to integer

  #Check X type
  if(!all(sapply(X, function(x){(x%%1)==0}))){stop("Please provide integer count matrix! \n")
  } else{X=matrix(as.integer(X), ncol=ncol(X), nrow=nrow(X))}

  colnames(X)=colNames; rownames(X)=rowNames
  X=X[rowSums(X)>0, ] #Remove empty samples
  X=X[, (colMeans(X==0)<(1-prevCutOff)) & (colSums(X)>(n*minFraction))] #Remove txa with low prevalence and those who do not meet the minFraction requirement
  if (distribution %in% c("ZIP","ZINB")) X = X[rowSums(X==0)>0, colSums(X==0)>0] #For a zero-inflated model, make sure every row and column has zeroes

  ##Build confounder matrix if applicable. For the preliminary trimming, we do not include an intercept, but we do include all the levels of the factors using contrasts=FALSE: we want to do the trimming in every subgroup, so no hidden reference levels. For the filtering we just use a model with an intercept and treatment coding, here the interest is only in adjusting the offset##
  if(!is.null(confounders)){
    if(class(confounders) %in% c("vector","matrix")){
      if(NROW(X)!=NROW(confounders)){ #Check dimensions
        stop("Data and confounder matrix do not have the same number of samples! \n")
      }
      if(is.vector(confounders)){
        confounders = as.matrix(confounders) #Convert to matrix if only 1 variable
      }
      if(is.null(colnames(confounders))){ #assign names if needed
        colnames(confounders) = paste("var",seq_len(NCOL(confounders)))
      }
      if(class(confounders[,1]) %in% c("character", "factor")){ #If character or factor, build model matrix
        if(class(confounders[,1])=="character"){
          warning("Converting character vector to factor! \n")
        }
        confModelMatTrim = model.matrix(object = as.formula(paste("~" , paste(colnames(confounders), collapse="+"),"-1")), contrasts.arg = apply(colnames(confounders),2,contrasts, contrasts=FALSE)) #No intercept for preliminary trimming
        confModelMat = model.matrix(object = as.formula(paste("~", paste(colnames(confounders), collapse="+"))), contrasts.arg = apply(colnames(confounders),2,contrasts, contrasts=TRUE)) #With intercept for filtering
      } else {
        confModelMat = confModelMatTrim = confounders
      }
    } else if(class(confounders) == "data.frame"){
      if(NROW(X)!=NROW(confounders)){ #Check dimensions
        stop("Data and confounder matrix do not have the same number of samples! \n")
      }
      confModelMatTrim = model.matrix( #No intercept for preliminary trimming
        object = as.formula(paste("~", paste(names(confounders), collapse="+"),"-1")),
        data = datFrame,
        contrasts.arg = lapply(confounders[sapply(confounders, is.factor)],
                               contrasts, contrasts=FALSE))
      confModelMat = model.matrix( #With intercept for filtering
        object = as.formula(paste("~", paste(names(confounders), collapse="+"))),
        data = datFrame,
        contrasts.arg = lapply(confounders[sapply(confounders, is.factor)],
                               contrasts, contrasts = TRUE))
    } else if(class(confounders)=="character"){
      if(classDat != "phyloseq"){
        stop("Providing confounders through variable names is only allowed if phyloseq object is provided! \n")
      }
      datFrame = data.frame(sample_data(dat)) # The dataframe with the confounders
      if(anyNA(datFrame)){stop("Confounders contin missing values!\n")}
      confModelMatTrim = model.matrix(
        object = formula(paste("~", paste(confounders, collapse="+"),"-1")),
        data = datFrame,
        contrasts.arg = lapply(datFrame[sapply(datFrame, is.factor), drop=FALSE][, confounders,drop=FALSE],
                               contrasts, contrasts=FALSE))
      confModelMat = model.matrix(
        object = formula(paste("~", paste(confounders, collapse="+"))),
        data = datFrame,
        contrasts.arg = lapply(datFrame[sapply(datFrame, is.factor), drop=FALSE][, confounders,drop=FALSE],
                               contrasts, contrasts = TRUE))
    } else{
      stop("Please provide the confounders either as matrix, dataframe, or character string! \n")
    }
  } else{
    confModelMat = confModelMatTrim = NULL
  }

  ##Build covariate matrix if applicable##
  # In this case we will 1) Include dummy's for every level of the categorical variable, and force them to sum to zero. This is needed fro plotting and required for normalization. 2) Exclude an intercept. The density function f() will provide this already. See introduction.
  if(!is.null(covariates)){
    if(class(covariates) == "data.frame"){
      if(NROW(X)!=NROW(covariates)){ #Check dimensions
        stop("Data and covariate matrix do not have the same number of samples! \n")
      }
      datFrame = covariates
      covariatesNames = names(covariates)
    } else if(class(covariates)=="character"){
      if(classDat != "phyloseq"){
        stop("Providing covariates through variable names is only allowed if phyloseq object is provided! \n")
      }
      if(covariates[[1]]=="all"){covariates = sample_variables(dat)} #Enable the "all" option if phyloseq object is provided
      datFrame = data.frame(sample_data(dat))[,covariates, drop=FALSE] # The dataframe with the covariates
      covariatesNames = covariates
    } else{
      stop("Please provide the covariates either as matrix, dataframe, or character string! \n")
    }

    if(any(sapply(datFrame, is.logical))){
      datFrame[,sapply(datFrame, is.logical)] = lapply(datFrame[sapply(datFrame, is.logical)], as.factor) #Convert logicals to factors
      #warning("Logicals converted to factors! \n", immediate. = TRUE). No warning needed
    }
    if(any(sapply(datFrame, is.integer))){
      datFrame[,sapply(datFrame, is.integer)] = lapply(datFrame[sapply(datFrame, is.integer)], as.numeric) #Convert integers to numeric
      warning("Integer values treated as numeric! \n", immediate. = TRUE)
    }
    if(any(sapply(datFrame, is.character))){
      datFrame[,sapply(datFrame, is.character)] = lapply(datFrame[sapply(datFrame, is.character)], factor) #Convert characters to factor
      warning("Character vectors treated as factors! \n", immediate. = TRUE)
    }
    nFactorLevels = sapply(datFrame, function(x){if(is.factor(x)) nlevels(x) else 1}) #Number of levels per factor
    datFrame[,sapply(datFrame, is.factor) & (nFactorLevels < 2)] = NULL #Drop factors with one level
    if(any(sapply(datFrame, is.factor) & (nFactorLevels < 2))){
      warning("The following variables were not included in the analyses because they are factors with only one level: \n", paste(covariates[sapply(datFrame, is.factor) & (nFactorLevels < 2)], sep = " \n"),immediate. = TRUE)
    }
    if(anyNA(datFrame)){
      NArows = apply(datFrame, 1, anyNA)
      if(all(NArows)){ stop("All samples have missing covariates")}
      X = X[!NArows,]
      datFrame = datFrame[!NArows,]
      if(!is.null(confModelMat)){
        confModelMat = confModelMat[!NArows,]
        confModelMatTrim = confModelMatTrim[!NArows,]
      }
      warning(paste("Some covariates contain missing values. We removed samples \n", paste(which(NArows), collapse = ", "), "\n prior to analysis." ),immediate. = TRUE)
    }

    # Center and scale the continuous covariates
    datFrame[sapply(datFrame, is.numeric)] = scale(datFrame[sapply(datFrame, is.numeric)])

    covModelMat = model.matrix(
      object = formula(paste("~", paste(covariatesNames, collapse="+"),"-1")),
      data = datFrame,
      contrasts.arg = lapply(datFrame[sapply(datFrame, is.factor)],
                             contrasts, contrasts=FALSE))

    # centMat  = t(sapply(seq_along(nFactorLevels), function(i){
    #   c(rep.int(0, sum(nFactorLevels[seq(0,i-1)])), #Zeroes before
    #     rep.int(if(nFactorLevels[i]==1) 0 else 1, nFactorLevels[i]), #Ones within the categorical variable
    #     rep.int(0, sum(nFactorLevels[-seq(1,i)]))) #Zeroes after
    # }))
    #
    # centMat = if(all(rowSums(centMat)==0)) {matrix(0, 1,sum(nFactorLevels))} else {centMat[rowSums(centMat)>0,, drop = FALSE]} #Remove zero rows, corresponding to the non-factors

    #Already prepare the matrix that defines the equations for centering the coefficients of the dummy variables
    tmp = buildCentMat(datFrame)
    centMat = tmp$centMat
    datFrame = tmp$datFrame

  } else {
    covModelMat = centMat = NULL
  }

  tic = proc.time() #Time the calculation
  tmp = switch(distribution,
               NB=RCM_NB(X, rowWeights = rowWeights, colWeights = colWeights, k = k,
                         confounders  = list(confounders = confModelMat, confoundersTrim = confModelMatTrim),
                         covariates = covModelMat, prevCutOff = prevCutOff, minFraction = minFraction, centMat = centMat, NBRCM = prevFit,...),
               ZIP=RCM_ZIP(X, rowWeights = rowWeights, colWeights = colWeights, k = k,
                           confounders  = list(confounders = confModelMat, confoundersTrim = confModelMatTrim),
                           covariates = covModelMat, prevCutOff = prevCutOff, minFraction = minFraction, ZIPRCM = prevFit,...),
               ZINB=RCM_ZINB(X,  rowWeights = rowWeights, colWeights = colWeights, k = k,
                             confounders  = list(confounders = confModelMat, confoundersTrim = confModelMatTrim),
                             covariates = covModelMat, prevCutOff = prevCutOff, minFraction = minFraction, ZINBRCM = prevFit, ...))
  if(classDat=="phyloseq"){tmp$physeq = prune_samples(rownames(X),prune_taxa(colnames(X),dat)) }
  tmp = within(tmp, {
    runtimeInMins = (proc.time()-tic)[1]/60 + if(is.null(prevFit)) {0} else {prevFit$runtimeInMins} # Sum the runtimes
    k = k #Store number of dimensions
  })
  class(tmp) = "RCM"
  tmp
}