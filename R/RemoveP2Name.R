RemoveP2Name <- function(n){
  # remove MOTIF.2.p2 -> MOTIF.2
  n.split <- strsplit(n, split = "\\.")[[1]]
  n.split <- n.split[1:(length(n.split) - 1)]
  n.fixed <- paste(n.split, collapse = ".")
  return(n.fixed)
}
