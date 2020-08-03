OrderPhaseMatrix <- function(complex.mat, order.by = "Liver", order.tissues = FALSE){
  # Given complex.mat, order rows such that the average angles are
  # increasing.
  # 
  # Args:
  # complex.mat -> complex matrix
  # 
  # Returns:
  # complex.mat.ordered <- complex matrix, rows ordered by 
  # average angle across columns
#   avg.phases <- apply(complex.mat, 1, function(x){
#     # re.avg <- sum(cos(Re(x))) / length(x)
#     # im.avg <- sum(sin(Re(x))) / length(x)
#     # vec.avg <- complex(real = re.avg, imaginary = im.avg)
#     return(Arg(vec.avg))
#   })
#   print(avg.phases)
  tissue.phases <- Arg(complex.mat[, order.by])
  # order by genes
  complex.mat <- complex.mat[order(tissue.phases, decreasing = FALSE), , drop = FALSE]
  # order by tissues: decreasing magnitude
  if (order.tissues){
    ordered.tissue.loadings <- names(sort(apply(complex.mat, 2,
                                           function(x){
                                            sum(Mod(x))
                                           }), decreasing = FALSE))
    complex.mat <- complex.mat[, ordered.tissue.loadings, drop = FALSE]
  }
  return(complex.mat)
}
