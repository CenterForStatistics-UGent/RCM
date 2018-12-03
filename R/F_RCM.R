#'A simple wrapper function for the RC(M) object for phyloseq objects and matrices.
#'
#' @param dat an nxp count matrix or a phyloseq object with an otu_table slot
#' @param k an integer, the number of dimensions of the RCM solution
#' @param round a boolean, whether to round to nearest integer. Defaults to FALSE.
#' @param prevCutOff a scalar, the prevalance cutoff for the trimming. Defaults to 2.5e-2
#' @param minFraction a scalar, each taxon's total abundance should equal at leat the number of samples n times minFraction, otherwise it is trimmed. Defaults to 10\%
#' @param rowWeights,colWeights character strings, the weighting procedures for the normalization of row and column scores. Defaults to "uniform" and "marginal" respectively
#' @param covariates In case "dat" is a phyloseq object, the names of the sample variables to be used as covariates in the constrained analysis, or "all" to indicate all variables to be used. In case "dat" is a matrix, a nxf matrix or dataframe of covariates. Character variables will be converted to factors, with a warning. Defaults to NULL, in which case an unconstrained analysis is carried out.
#' @param confounders In case "dat" is a phyloseq object, the names of the sample variables to be used as confounders to be filtered out. In case "dat" is a matrix, a nxf matrix or dataframe of confounders. Character variables will be converted to factors, with a warning. Defaults to NULL, in which case no filtering occurs.
#' @param ... Further arguments passed on to the RCM.NB() function
#'
#'@description This is a wrapper function, which currently only fits the negative binomial distribution, but which could easily be extended to other ones.
#'
#'@details This function should be called on a raw count matrix, without rarefying or normalization to proportions. This functions trims on prevalence and total abundance to avoid instability of the algorithm. Covariate and confounder matrices are constructed, so that everything is passed on to the workhorse function RCM.NB() as matrices.
#'@seealso \code{\link{RCM_NB}},\code{\link{plot.RCM}},\code{\link{residualPlot}},\code{\link{plotRespFun}}
#'
#' @return see \code{\link{RCM_NB}}
#' @importFrom nleqslv nleqslv
#' @importFrom tensor tensor
#' @importFrom alabama constrOptim.nl
#' @import VGAM
#' @import phyloseq
#' @importFrom stats as.formula contrasts dnbinom dpois glm integrate quantile
#' @export
#' @examples
#' data(Zeller)
#' require(phyloseq)
#' tmpPhy = prune_taxa(taxa_names(Zeller)[1:100],
#' prune_samples(sample_names(Zeller)[1:50], Zeller))
#' zellerRCM = RCM(tmpPhy, round = TRUE)
RCM = function(dat, k = 2, round=FALSE, prevCutOff = 0.05, minFraction = 0.1,
               rowWeights = "uniform", colWeights = "marginal",
               covariates = NULL, confounders = NULL, ...){
  classDat = class(dat) #The class provided

  ##The count data##
  if (classDat=="matrix"){
    X = dat
  }else  if(classDat=="phyloseq"){
    X = if (taxa_are_rows(dat)) {t(otu_table(dat)@.Data)
      }else {otu_table(dat)@.Data}
  } else {stop("Please provide a matrix or a phyloseq object! \n")}
  p=ncol(X); n=nrow(X)

  if(is.null(colnames(X))){colnames(X)=seq_len(ncol(X))}
  if(is.null(rownames(X))){rownames(X)=seq_len(nrow(X))}
  colNames = colnames(X); rowNames =rownames(X)
  if(round) {X=round(X, 0) }#Round to integer

  #Check X type
  if(!all(floor(X)==X)){
    stop("Please provide integer count matrix (not a matrix of proportions!),
         or set 'round' to TRUE! \n")
  } else{X=matrix(as.integer(X), ncol=p, nrow=n)}

  colnames(X)=colNames; rownames(X)=rowNames
  X=X[, (colMeans(X==0)<(1-prevCutOff)) & (colSums(X)>(n*minFraction))]
  #Remove txa with low prevalence and
  #those who do not meet the minFraction requirement
  if(any(rowSums(X)==0)){
    warning(immediate. = TRUE,
            paste0("Samples \n",paste(collapse = ", ",
                                      rownames(X)[rowSums(X)==0]),
                   "\n contained no more reads after trimming taxa
                   with low prevalence, and were excluded from the fit"))}
  rowIDkeep = rowSums(X)>0
  X=X[rowIDkeep, ] #Remove empty samples
  n=nrow(X)
  if(classDat=="phyloseq"){
    dat = prune_samples(x = dat, samples = rowIDkeep)
  } else {
  confounders = confounders[rowIDkeep, ]
  covariates = covariates[rowIDkeep, ]
  }

  ##Build confounder matrix if applicable. ##
  if(!is.null(confounders)){
tmp = buildConfMat(confounders, n, dat)
confModelMatTrim = tmp$confModelMatTrim
confModelMat = tmp$confModelMat
rm(tmp)
  } else{
    confModelMat = confModelMatTrim = NULL
  }

  ##Build covariate matrix if applicable##
  if(!is.null(covariates)){
    tmp = buildCovMat(covariates, n, dat)
    covModelMat = tmp$covModelMat
    covariates = tmp$covariates
    #Already prepare the matrix that defines the equations for
    #centering the coefficients of the dummy variables
    tmp = buildCentMat(tmp$datFrame)
    centMat = tmp$centMat
    datFrame = tmp$datFrame
    rm(tmp)

  #Remove rows with NA's, we might want to find something better or
    #leave it to the end user
  if(anyNA(datFrame)){
    NArows = apply(datFrame, 1, anyNA)
    if(all(NArows)){ stop("All samples have missing covariates")}
    X = X[!NArows,]
    datFrame = datFrame[!NArows,]
    if(!is.null(confModelMat)){
      confModelMat = confModelMat[!NArows,]
      confModelMatTrim = confModelMatTrim[!NArows,]
    }
    warning(paste("Some covariates contain missing values.
                  We removed samples \n",
                  paste(which(NArows), collapse = ", "),
                  "\n prior to analysis." ),immediate. = TRUE)
    }
  } else {
    covModelMat = centMat = NULL
  }
if(anyNA(X)) {stop("NA values present in count matrix,
                   please filter these out first!\n")}
  tic = proc.time() #Time the calculation
  tmp = RCM_NB(X, rowWeights = rowWeights, colWeights = colWeights, k = k,
  confounders  = list(confounders = confModelMat,
                      confoundersTrim = confModelMatTrim),
  covariates = covModelMat, prevCutOff = prevCutOff, minFraction = minFraction,
  centMat = centMat,...)
  if(classDat=="phyloseq"){
    tmp$physeq = prune_samples(rownames(tmp$X),prune_taxa(colnames(tmp$X),dat))}
  tmp = within(tmp, {
    runtimeInMins = (proc.time()-tic)[1]/60 # Sum the runtimes
    k = k #Store number of dimensions
  })
  tmp$call = match.call()
  class(tmp) = "RCM"
  tmp
}
