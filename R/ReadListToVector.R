# ReadListToVector.R

ReadListToVector <- function(fname, HEADER=FALSE){
  # read file name, one column with no colnames as a vector.
  vec <- read.table(fname, header = HEADER)
  vec <- as.character(unlist(vec))
  return(vec)
}