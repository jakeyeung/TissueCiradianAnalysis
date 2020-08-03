ModelStrToModel <- function(jmod){
  # Liver_SV129,Kidney_SV129.Liver_BmalKO,Kidney_BmalKO-Liver_SV129,Kidney_SV129 - >
  # c(Liver_SV129,Kidney_SV129, Liver_SV129,Kidney_SV129;Liver_BmalKO,Kidney_BmalKO)  # need some reordering magic
  jmod.long <- gsub("\\.", ";", jmod)
  jmod.long <- strsplit(jmod.long, "-")[[1]]
  # rearrange each mod so that SV129 goes before BmalKO
  jmod.long.sorted <- rep(NA, length(jmod.long))
  i <- 1
  for (j in jmod.long){
    j.sorted <- paste(sort(strsplit(j, ";")[[1]], decreasing = TRUE), collapse = ";")
    jmod.long.sorted[i] <- j.sorted
    i <- i + 1
  }
  return(jmod.long.sorted)
}

ModelToTissue <- function(jmod, jsplit = "_", jsplit.diffparam = ";", jsplit.sameparam = ","){
  # Try to convert model name to tissue name (Liver, Kidney)
  # mod names: Liver_SV129, Kidney_BmalKO etc (so split by "_" by default)
  tissue <- gsub(pattern = jsplit.diffparam, replacement = jsplit.sameparam, x = jmod)
  tissue <- strsplit(tissue, jsplit.sameparam)[[1]]
  tissue <- sapply(tissue, function(tiss) strsplit(tiss, jsplit)[[1]][[1]], USE.NAMES = FALSE)
  
  # tissue <- unique(sapply(jmod, function(m) strsplit(m, jsplit)[[1]][[1]], USE.NAMES = FALSE))
  if (length(tissue) == 1){
    return(tissue)
  } else {
    return(unique(tissue))
  }
}