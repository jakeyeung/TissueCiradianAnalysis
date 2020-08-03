# Calculations on gene and tissues

GetMeanVarByTissues <- function(exprs, tissue.names){
  # Calculate mean and variance for each gene per tissue
  N <- nrow(exprs) * length(tissue.names)  # one measurement for each gene for all tissues.
  mean.var <- list(mean=matrix(NA, nrow=nrow(exprs), ncol=length(tissue.names),
                               dimnames = list(rownames(exprs), 
                                               tissue.names)), 
                   var=matrix(NA, nrow=nrow(exprs), ncol=length(tissue.names),
                              dimnames = list(rownames(exprs), tissue.names)))
  for (j in 1:length(tissue.names)){
    tissue <- tissue.names[j]
    gene.tissue.exprs <- exprs[, grepl(tissue, colnames(exprs))]
    # calculate mean and var, by row
    exprs.mean <- apply(gene.tissue.exprs, 1, mean)
    exprs.var <- apply(gene.tissue.exprs, 1, var)
    # append to matrix 
    mean.var$mean[, j] <- exprs.mean
    mean.var$var[, j] <- exprs.var
  }
  # Make dataframe for ggplot2 and other analyses
  mean.var.df <- data.frame(mean=as.vector(mean.var$mean), var=as.vector(mean.var$var))
  return(mean.var.df)
}
