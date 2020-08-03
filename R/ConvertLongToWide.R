ConvertLongToWide <- function(long.df, measurement.var = "exprs.transformed"){
  library(plyr)
  wide.df <- dcast(long.df, gene ~ tissue, value.var = measurement.var)
  # first row is gene name, let's make them rowname and remove first column.
  
  rownames(wide.df) <- wide.df$gene
  
  wide.df <- subset(wide.df, select = -c(gene))
  
  return(wide.df)
}