WriteListToFile <- function(lst, outf){
  sink(file = outf)
  for (g in lst){
    cat(g)
    cat("\n")
  }
  sink()
  return(NA)
}