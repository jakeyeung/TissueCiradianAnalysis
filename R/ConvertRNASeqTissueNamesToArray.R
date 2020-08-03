TissueNames <- function(){
  array.names <- c("Adr", "Aorta", "BFAT", "BS", "Cere", "Heart", "Hypo", "Kidney", "Liver", "Lung", "Mus", "WFAT")
  return(array.names)
}

TimesNames <- function(){
  # Begin CT22, end CT64, take every 6 hours
  return(seq(22, 64, 6))
}

RnaseqNames <- function(){
  rnaseq.names <- c('Adr','Aor','BFat','Bstm','Cer','Hrt','Hyp','Kid','Liv','Lun','Mus','WFat')
  return(rnaseq.names)
}

ConvertRNASeqTissueNamesToArray <- function(jnames){
  array.names <- TissueNames()
  rnaseq.names <- RnaseqNames()
  names.dic <- setNames(object = array.names, nm = rnaseq.names)
  return(names.dic[jnames])
}

