# Jake Yeung
# DataHandlingFunctions.R
# November 10 2014
# Simple functions for easier manipulating of data

CopyZT0ToZT24 <- function(enrichment, twindow = 6, jorder = TRUE, convert.cname=FALSE){
  # complete the circle
  enrichment$tmid <- enrichment$tstart + twindow / 2
  enrichment$tmid <- sapply(enrichment$tmid, function(tmid) ifelse(tmid >= 24, tmid - 24, tmid))
  
  enrichment.24 <- subset(enrichment, tmid == 0)
  enrichment.24$tmid <- 24
  enrichment <- bind_rows(enrichment, enrichment.24)
  enrichment <- enrichment %>% arrange(Term, tmid)
  if (convert.cname){
    enrichment <- dplyr::rename_(enrichment, phase = "tmid")
    enrichment <- dplyr::rename_(enrichment, amp = "minuslogpval")
  }
  return(enrichment)
}


Vectorize(IsBtwnTimes <- function(phase, tstart, tend){
  # check if phase (between 0 to 24, is between tstart and tend, considering the modulo)
  if (tend > tstart){
    # easy case
    is.btwn <- phase >= tstart & phase <= tend
  } else {
    # harder case, must consider the modulo
    is.btwn <- phase >= tstart | phase <= tend
  }
  # replace NAs with FALSE
  is.btwn[which(is.na(is.btwn))] <- FALSE
  return(is.btwn)
}, vectorize.args="phase")

Peek <- function(dat, N=5){
  # Peek at first N rows and N columns of data
  # 
  # ARGS
  #   dat: dataframe or matrix
  # 
  # OUTPUT
  #   print to user the first N rows and N columns of data
  if (N > nrow(dat)){
    N.row <- nrow(dat)
  } else {
    N.row <- N
  }
  if (N > ncol(dat)){
    N.col <- ncol(dat)
  } else {
    N.col <- N
  }
  print(dat[1:N.row, 1:N.col])
}

GetRowIndx <- function(dat, cname){
  indx <- which(colnames(dat) == cname)
  if (length(indx) == 0) warning("Warning, cname not found in column names of dat")
  return(indx)
}

MatchGroup <- function(x, group, no.match.str = "Flat"){
  # group is list of comma separated strings
  # return group[i] which contains x
  group.lst <- sapply(group, function(g) strsplit(g, ","), USE.NAMES = FALSE)
  match <- group.lst[which(sapply(group.lst, function(r) x %in% r))]
  if (length(match) == 0){
    return(no.match.str)
  } else {
    return(paste0(match[[1]], collapse=","))
  }
}