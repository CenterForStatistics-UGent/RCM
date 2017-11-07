#' A phyloseq wrapper for the gomms() function
#'
#' @param physeq a phyloseq object
#' @param k an integer, the dimension of the solution
#' @param ... passed on to the gomms function
#'
#' @return See the gomms() function
gommsPhy = function(physeq, k = 3,...){
  require(gomms)
  otuTab = if(class(physeq) == "phyloseq") {
    if(taxa_are_rows(physeq)) t(otu_table(physeq)@.Data) else otu_table(physeq)@.Data
    } else {physeq}
  gomms(otuTab, n.factors = k,...)
  }