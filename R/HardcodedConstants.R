GetSuffix <- function(jweight, use.sql, jmodstr, jcutoffstr, include.promoter=FALSE, mindist=0, do.pairs=FALSE){
  if (!include.promoter){
    suffix <- paste0(".weight.", jweight, ".sql.", use.sql, ".mod.", jmodstr, ".dhscutoff.", jcutoffstr, ".final")
  } else {
    suffix <- paste0(".mindist.", mindist, ".inclprom.", include.promoter, ".weight.", jweight, ".sql.", use.sql, ".mod.", jmodstr, ".dhscutoff.", jcutoffstr, ".final")
  }
  if (do.pairs){
    suffix <- paste0(suffix, ".dopairsInDHSFinal.", do.pairs)
  }
  return(suffix)
}

GetESubDir <- function(do.center, jmodstr, jweight){
  E.subdir <- paste0("centered.", do.center, ".mod.", jmodstr, ".weightcutoff.", jweight)
}