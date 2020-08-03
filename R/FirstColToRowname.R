# FirstColToRowname.R

FirstColToRowname <- function(dat){
  # if dat first column contains rownames you want, move it to
  # rownames, then return dat
  rownames(dat) <- dat[, 1]
  dat <- dat[, -1]
  return(dat)
}