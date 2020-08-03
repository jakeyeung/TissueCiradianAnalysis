ShannonEntropy <- function(p, normalize = FALSE){
  if (normalize){
    p <- p / sum(p)
  }
  s <- p * log2(1 / p)
  return(sum(s))
}