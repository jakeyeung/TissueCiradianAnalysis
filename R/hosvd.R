# November 5 2014
# taken from James Li, author of rTensor package
# contact: jamesyili@gmail.com
# 
# hosvd.R
# calculate hosvd. Differs from rTensor v 1.1 because
# this function will calculate all eigenvectors.
# 
# It's changed so it will calculate the ranks[m] each time.

hosvd <- function (tnsr, ranks = NULL){
  num_modes <- tnsr@num_modes
  if (is.null(ranks)) {
    cat("!ranks not provided so left singular matrices will not be truncated.\n")
    ranks <- tnsr@modes
  }
  #  original code from v 1.1 is commented
  # 
  #   pb <- txtProgressBar(min = 0, max = num_modes, style = 3)
  #   U_list <- vector("list", num_modes)
  #   for (m in 1:num_modes) {
  #     U_list[[m]] <- (svd(rs_unfold(tnsr, m = m)@data)$u)[, 
  #                                                         1:ranks[m]]
  #     setTxtProgressBar(pb, m)
  #   }
  
  # new fix from James Li
  pb <- txtProgressBar(min = 0, max = num_modes, style = 3)
  U_list <- vector("list", num_modes)
  for (m in 1:num_modes) {
    temp_mat <- rs_unfold(tnsr, m = m)@data
    U_list[[m]] <- svd(temp_mat,nu=ranks[m])$u
    setTxtProgressBar(pb, m)
  }
  close(pb)
  Z <- ttl(tnsr, lapply(U_list, t), ms = 1:num_modes)
  est <- ttl(Z, U_list, ms = 1:num_modes)
  resid <- fnorm(est - tnsr)
  list(Z = Z, U = U_list, est = est, fnorm_resid = resid)
}


