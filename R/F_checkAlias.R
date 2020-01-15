#' Check for alias structures in a dataframe, and throw an error when one is found
#' @param datFrame the data frame to be checked for alias structure
#' @param covariatesNames The names of the variables to be considered
#'
#' @return Throws an error when an alias structure is detected,
#' returns invisible otherwise
#'
#' @importFrom stats alias
#' @export
#' @examples
#' #Make a dataframe with aliased variables
#' df = data.frame(foo = rnorm(10), baa = rep(c(TRUE, FALSE), each = 5),
#' foo2 = factor(rep(c("male", "female"), each = 5)))
#' checkAlias(df, c("foo", "baa"))
#' #Check test files for the error being thrown
checkAlias = function(datFrame, covariatesNames){
    mockDf = cbind("Mock" = 1, datFrame)
    #Fake dataframe for syntax purposes
    Alias = alias(object = formula(paste("Mock~",
                                paste(covariatesNames, collapse = "+"),
                                "-1")), mockDf)
    if(!is.null(Alias$Complete)){
        stop("Sample variables\n'", paste(rownames(Alias$Complete),
                                          collapse ="' and '"),
            "'\nare aliased with other variables.
            Drop some sample-variables and try again.")
    } else {
        return(invisible())
    }
}
