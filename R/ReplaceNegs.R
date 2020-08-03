# replace negative values with 0 in a vector.
# use "apply" for dataframes.
# ReplaceNegs.R
# jakeyeung
# December 4 2014

ReplaceNegs <- function(vector, replace.with=0){
  x <- vector
  x[which(x < 0)] <- replace.with
  return(x)
}