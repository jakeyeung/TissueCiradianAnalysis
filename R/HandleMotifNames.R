source("/home/yeung/projects/tissue-specificity/scripts/functions/RemoveP2Name.R")

RemoveCommasBraces <- function(m){
  return(gsub("\\{|\\}|\\,", replacement = "\\.", m))
} 

RemoveDashes <- function(m){
  return(gsub("-", replacement = "_", m))
} 

RemoveP2AndCommaBracesDashes <- function(m){
  m <- RemoveP2Name(m)
  m <- RemoveCommasBraces(m)
  m <- RemoveDashes(m)
  return(m)
}

PcToMotif <- function(v, pc, top.n = 3){
  # for PMD -> MARA analysis 
  x <- v[, pc]
  x <- x[which(x != 0)]
  motif.names <- names(sort(x, decreasing = TRUE))
  return(paste(motif.names[1:top.n], collapse = "-"))
}
