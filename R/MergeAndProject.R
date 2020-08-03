MergeAndProject <- function(array.normalized, rna.seq.exprs, T = 24){
  # given array.normalized and rna.seq.exprs, project onto frequency domain.
  # use T = 24 for circadian.
  
  # source functions
  source(file.path("scripts", "functions", "GetTissueTimes.R"))
  
  merged.dat <- cbind(array.normalized, rna.seq.exprs)  
  
  times.array <- as.integer(GetTimes(samp.names = colnames(array.normalized)))
  times.rnaseq <- as.integer(GetTimes(samp.names = colnames(rna.seq.exprs)))
  tissues <- GetTissues(colnames(array.normalized))
  
  omega <- 2 * pi / T
  n.tissues <- 12
  n.timepts.array <- 24
  interval.array <- 2
  n.timepts.rnaseq <- 8
  interval.rnaseq <- 6
  
  # init output matrix
  Y.gcs <- matrix(nrow = length(filtered.genes), ncol = length(tissues), dimnames = list(filtered.genes, tissues))
  
  # Begin projection: iterate for each tissue.
  for (tissue in tissues){
    # Get subset of merged.dat containing tissue
    merged.tissue <- as.matrix(merged.dat[, grepl(tissue, colnames(merged.dat))])
    # multiply by times.array and rna.seq array (unnormalizd fourier transform-ish)
    transformed <- merged.tissue %*% exp(-1i * omega * c(times.array, times.rnaseq))
    # normalize by number of data points
    transformed <- 1 / length(c(times.array, times.rnaseq)) * transformed
    Y.gcs[, tissue] <- transformed
  }
  return(Y.gcs)
}


