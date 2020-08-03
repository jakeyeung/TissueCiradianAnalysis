# Jake Yeung
# Dec 4 2014

GetTissueSpecificMatrix <- function(n.tissues, n.timepoints){
  # set my tissue specific factors, matches my response y to a 
  # specific tissue. 12 tissues, 24 time points
  # n.tissues <- 12  # 12 conditions
  # n.timepoints <- 24  # 24 time points
  
  # initialize first tissue vector, then cbind for next tissues
  tissue.factors <- rep(0, n.tissues * n.timepoints)
  tissue.factors[1:n.timepoints] <- rep(1, n.timepoints)
  # print(length(tissue.factors))
  
  # do same for all other tissues, cbind to tissue.factors
  for (c in 2:n.tissues){  # start at 2 because we did 1 already
    start.i <- (c - 1) * n.timepoints + 1  # starts at 25, if c = 2
    end.i <- start.i + n.timepoints - 1  # ends at 48, if c = 2 
    t.fac <- rep(0, n.tissues * n.timepoints)
    t.fac[start.i:end.i] <- rep(1, n.timepoints)
    # print(length(t.fac))
    tissue.factors <- cbind(tissue.factors, t.fac)
  }
  return(tissue.factors)
}