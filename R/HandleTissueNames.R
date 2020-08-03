GetTissueNames <- function(tissue.names.all){
  # Convert gene names Adr_CT22 to Adr for all tissue.names.all
  # 
  # Args:
  # tissue.names.all List of names like Adr_CT22 ...
  # 
  # Returns:
  # tissue.names List of names like Adr ...
  
  tissues.names <- strsplit(tissue.names.all, '_')
  tissues.names <- unique(unlist(lapply(tissues.names, '[[', 1)))
  return(tissues.names)
} 