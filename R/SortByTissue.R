SortByTissue <- function(dat, by.var = "tpm"){
  # Sort signal by tissue by renaming factors
  dat.ordered <- dat[order(dat[[by.var]], decreasing = TRUE), ]
  dat.ordered$tissue <- factor(dat.ordered$tissue, levels = dat.ordered$tissue)
  return(dat.ordered)
}