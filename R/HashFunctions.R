Vectorize(AssignHash <- function(x, jhash, null.fill = NA){
  # assign hash key to hash value, handle NULLs
  # null.fill = "original", returns original value x into jhash
  x.mapped <- jhash[[as.character(x)]]
  if (is.null(x.mapped)){
    if (is.na(null.fill)){
      x.mapped <- null.fill
    } else if (as.character(null.fill) == "original"){
      x.mapped <- x
    } else {
      x.mapped <- null.fill
    }
  } 
  return(x.mapped)
}, vectorize.args = "x")
