# RemoveProblemGenes.R
# Jake Yeung
# Dec 5 2014
# After adjusting, there may be one or two problem genes.
# Define problem genes as expression that goes negative. That should 
# not happen if we used nls with proper constraints.
# These problem genes arise from defaulting to "lm" fit.

RemoveProblemGenes <- function(dat, verbose=TRUE){
  negs <- apply(dat, 1, function(x){
    if (min(x) < 0){
      return(1)
    } else {
      return(0)
    }
  })
  problem.genes <- names(negs[which(negs == 1)])
  rm(negs)
  
  genes <- rownames(dat)
  filtered.genes <- genes[which(!genes %in% problem.genes)]
  dat <- dat[filtered.genes, ]
  
  # Tell user how many problem genes existed.
  if (verbose == TRUE){
    print(paste0(length(problem.genes), " genes because they had negative expression."))
    print(paste("Problem genes:", problem.genes))
  }
  return(dat)
}
