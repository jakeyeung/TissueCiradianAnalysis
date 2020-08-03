# RemoveExtension.R
# February 26 2014
# Jake Yeung

RemoveExtension <- function(fname){
  # removes any .txt or .pdf, returns the fname without the extension
  fname.split <- strsplit(fname, '\\.')
  # remove last element from split string
  fname.split <- fname.split[[1]][-length(fname.split[[1]])]
  
  # reassemble any potential dots in filename from before
  fname.noext <- paste(fname.split, collapse = ".")
  return(fname.noext)
}