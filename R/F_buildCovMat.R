#' A function to build the covariate matrix of the constraints
#'
#' @param covariates the covariates, either as dataframe or as character string
#' @param n the number of samples
#' @param dat the phyloseq object
#'
#' In this case we will 1) Include dummy's for every level of the categorical variable, and force them to sum to zero. This is needed for plotting and required for reference level indepenent normalization. 2) Exclude an intercept. The density function f() will provide this already. See introduction.
#'
#' @return a list with components
#' \item{covModelMat}{The model matrix}
#' \item{datFrame}{The dataframe used to construct the model matrix}
buildCovMat = function(covariates, n,  dat){
  if(is.data.frame(covariates)){
    datFrame = covariates
    covariatesNames = names(covariates)
  } else if(is.character(covariates)){
    if(!is(dat,"phyloseq")){
      stop("Providing covariates through variable names is only allowed if phyloseq object is provided! \n")
    }
    if(covariates[[1]]=="all"){covariates = sample_variables(dat)} #Enable the "all" option if phyloseq object is provided
    datFrame = data.frame(sample_data(dat))[,covariates, drop=FALSE] # The dataframe with the covariates
    covariatesNames = covariates
  } else {
    stop("Please provide the covariates either as dataframe or as character string! \n")
  }
  if(n!=NROW(datFrame)){ #Check dimensions
    stop("Data and covariate matrix do not have the same number of samples! \n")
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
  covariatesNames = covariatesNames[!(sapply(datFrame, is.factor) & (nFactorLevels < 2))] #Drop factors with one level
  nFactorLevels = nFactorLevels[covariatesNames]
  datFrame = datFrame[,covariatesNames, drop=FALSE]
  if(any(sapply(datFrame, is.factor) & (nFactorLevels < 2))){
    warning("The following variables were not included in the analyses because they are factors with only one level: \n", paste(covariates[sapply(datFrame, is.factor) & (nFactorLevels < 2)], sep = " \n"),immediate. = TRUE, call. = FALSE)
      }
  # Center and scale the continuous covariates
  datFrame[sapply(datFrame, is.numeric)] = scale(datFrame[sapply(datFrame, is.numeric)])

  covModelMat = model.matrix(
    object = formula(paste("~", paste(covariatesNames, collapse="+"),"-1")),
    data = datFrame,
    contrasts.arg = lapply(datFrame[sapply(datFrame, is.factor)],
                           contrasts, contrasts=FALSE))
  if(NCOL(covModelMat)==1) stop("A constrained ordination with only one variable is meaningless.\nPlease provide more covariates or perform an unconstrained analysis.", call. = FALSE)

  list(covModelMat = covModelMat, datFrame = datFrame)
}
