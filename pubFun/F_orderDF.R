#' A function to reorder the methods
orderDF = function(df, groupsMeth){
  df$Method = factor(df$Method, labels = labelsMeth, levels = levelsMeth, ordered = TRUE)
  df$groupMethod = groupsMeth[match(df$Method, unique(groupsMeth))]
  droplevels(df)
}