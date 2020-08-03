# Jake Yeung
# Nov 5 2015
# Functions to handle sample names

ShortenSampNames <- function(long.names, show="tissue") {
  # Function to shorten sample names of a certain form. 
  # User options to determine how much to shorten
  #
  # INPUT: 
  # [GSM1321058_BS58_MoGene1.0ST.CEL, GSM1321058_BS59_MoGene1.0ST.CEL ..., ]
  # 
  # OUTPUT:
  # [BS58, BS59, ..., ]  if show == "tissue.time"
  # or
  # [BS, BS, ..., ]  if show == "tissue"
  # or
  # [58, 59, ..., ] if show == "time"
  # 
  # Whether you want the number afterwards is user defined.
  
  long.names <- strsplit(long.names, '_')
  short.names <- sapply(long.names, function(x) {
    sample.id <- x[2]    # sample name with time point
    
    # THREE WAYS OF SHOWING SAMPLES NAMES:
    
    if (show == "tissue.time") {
      # 1.
      # show sample name with timepoint e.g. Liver34
      samp.name <- sample.id 
    }
    else if (show == "tissue") {
      # 2.
      # show opnly tissue component in sample name e.g. "Adr", "Liver"
      samp.name <- substring(sample.id, 1, nchar(sample.id) - 2)      
    }
    else if (show == "time") {
      # 3.
      # show only time component in sample name e.g. 34
      samp.name <- substring(sample.id, nchar(sample.id) - 1, nchar(sample.id)) 
    }
  })
  return(short.names)
}

GetTimes <- function(dat.colnames, n.digits=2){
  # Parse colnames and extract times. e.g. Adr30 -> 30
  # Expects colnames with times at last part of string. Default is last two characters.
  
  char.from.end <- n.digits - 1  # if n.digits=2, extract second last digit (n.digit - 1) to last digit
  
  times <- sapply(dat.colnames, function(colname, char.from.end){
    time <- substring(colname, nchar(colname) - char.from.end, nchar(colname))
  }, char.from.end)
  return(unique(times))
}

GetTissueNames <- function(dat.colnames, dat.type){
  # Parse colnames and extract tissue names. 
  # e.g. Aorta62 -> Aorta or Adr_CT40 -> Adr
  # Invariant is that there are constant number of digits from
  # end of tissue name to end of string.
  # 
  # dat.type -> "rna.seq" or "array". 
  # rna.seq has colnames like Adr_CT40
  # array has colnames like Aorta62
  
  if (missing(dat.type)){
    warning('dat.type unspecified, defaulting to colnames of form: "Cere28"...')
    dat.type <- "array"
  }
  tissues <- sapply(dat.colnames, function(colname, dat.type){
    if (dat.type == "rna.seq"){
      tissue <- substring(colname, 1, nchar(colname) - 5)
    } else if (dat.type == "array"){
      tissue <- substring(colname, 1, nchar(colname) - 2)
    }
    return(tissue)
  }, dat.type)
  return(unique(tissues))
}


