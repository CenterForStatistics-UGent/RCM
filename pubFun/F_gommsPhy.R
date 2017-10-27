#' A phyloseq wrapper for the gomms() function
#'
#' @param physeq a phyloseq object
#' @param k an integer, the dimension of the solution
#' @param ... passed on to the gomms function
#'
#' @return See the gomms() function
gommsPhy = function(physeq, k = 3,...){
  require(gomms)
  otuTab = otu_table(physeq)@.Data
  if(taxa_are_rows(physeq)) gomms(t(otuTab),...) else gomms(otuTab,...)
  }