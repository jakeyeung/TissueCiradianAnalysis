OuterComplex <- function(u, v){
  # Given two complex vectors, get the outer product.
  # 
  # given u = x + i*y and v = x2 + i*y2, we want to get
  # u * Conj(v) = (x + i*y) * (x2 - i*y2) = x*x2 - y*y2 - 1i (x*y2 - x2*y)
  # 
  # Note in R: using outer(u, v) returns (x + i*y) * (x2 + i*y2) which
  # is not equivalent.
  # 
  # Args:
  # u -> complex number/vector
  # v -> complex number/vector
  # 
  # Returns
  # o -> outer product of u and v (u * Conj(v))
  # 
  o <- u %*% Conj(v)
  return(o)
}