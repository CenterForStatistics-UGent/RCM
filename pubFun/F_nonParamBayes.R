#' A function to implement a bayesian approach to ordination, this is quite slow
nonParamBayes = function(dat, hyper = list( nv = 3, a.er = 1, b.er = 0.3, a1 = 3, a2 = 4, m = 10, alpha = 10, beta = 0 )){
  require(DirFactor)
  classDat = class(dat) #The class provided
  ##The count data##
  if (classDat=="matrix"){
    X = dat
  }else  if(classDat=="phyloseq"){
    X = if (taxa_are_rows(dat)) t(otu_table(dat)@.Data) else otu_table(dat)@.Data
  } else {stop("Please provide a matrix or a phyloseq object! \n")}
  DirFactor(t(X), hyper, step = 100)
}