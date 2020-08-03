# Functions to get tissue and times: new after column fixed to array.
# Jake Yeung
# Dec 4 2014
# GetTissueTimes.R

GetTissues <- function(samp.names, get_unique=TRUE){
  # Samp names of form: WFAT48 (as vector)
  # return WFAT (as vector, unique)
  tissues <- unlist(lapply(samp.names, function(samp.name){
    substr(samp.name, 1, nchar(samp.name) - 2)
  }))
  if (get_unique){
    return(unique(tissues))
  } else {
    return(tissues)
  }
}

GetTimes <- function(samp.names, get_unique=TRUE){
  # Samp names of form: WFAT48 (as vector)
  # return 48 (as vector, unique)
  times <- unlist(lapply(samp.names, function(samp.name){
    substr(samp.name, nchar(samp.name) - 1, nchar(samp.name))
  }))
  if (get_unique){
    return(as.numeric(unique(times)))
  } else {
    return(as.numeric(times))
  }
}

GetTissues.merged <- function(samp.names){
  # Sampnames of form: Kidney60.array
  # return NON-unique
  tissues <- unlist(lapply(samp.names, function(samp.name.full){
    samp.name <- strsplit(samp.name.full, "[.]")[[1]][[1]]
    substr(samp.name, 1, nchar(samp.name) - 2)
  }))
  return(tissues)
}

GetTimes.merged <- function(samp.names){
  # Sampnames of form: Kidney60.array
  # return NON-unique
  times <- unlist(lapply(samp.names, function(samp.name.full){
    samp.name <- strsplit(samp.name.full, "[.]")[[1]][[1]]
    substr(samp.name, nchar(samp.name) - 1, nchar(samp.name))
  }))
  return(as.numeric(times))}

GetExperiments.merged <- function(samp.names){
  experiments <- unlist(lapply(samp.names, function(samp.name.full){
    strsplit(samp.name.full, "[.]")[[1]][[2]]
  }))
}