GetTopGenesByPeriod <- function(dat, jperiod){
  # dat from dat.complex.all_T
  dat.sub <- subset(dat, period == jperiod)
  return(dat.sub[order(Mod(dat.sub$exprs.transformed), decreasing = TRUE), ])
}

OrderDecreasing <- function(dat, jfactor, jval){
  # Reorder factors by decreasing value of jval
  dat[[jfactor]] <- factor(dat[[jfactor]], levels = dat[[jfactor]][order(dat[[jval]], decreasing = TRUE)])
  return(dat)
}

GetIndex <- function(x){
  # Get ith estimate given vector
  # make a hash table
  x.o <- order(x, decreasing = TRUE)
  sapply(x, function(y) which(x.o == which(x == y)))
}