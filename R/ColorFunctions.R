#library(PhaseHSV)

FactorToHex <- function(i, colhash){
  return(colhash[[as.character(i)]])
}

HexColorToRgb <- function(hex){
  # HEX input: #FFC8DD
  R <- as.integer(paste0("0x", substr(hex, 2, 3)))
  G <- as.integer(paste0("0x", substr(hex, 4, 5)))
  B <- as.integer(paste0("0x", substr(hex, 6, 7)))
  return(paste(R, G, B, sep = ","))
}

AddAlphaToHexColor <- function(hex, alpha){
  # add alpha to hex, return hex
  R <- as.integer(paste0("0x", substr(hex, 2, 3)))
  G <- as.integer(paste0("0x", substr(hex, 4, 5)))
  B <- as.integer(paste0("0x", substr(hex, 6, 7)))
  HSV <- rgb2hsv(R, G, B, maxColorValue = 255)
  # add alpha and return
  return(hsv(HSV[1], HSV[2], HSV[3], alpha = alpha))
}

SaturationCurve <- function(x, Vmax, k, x0){
  y <- (Vmax * 2 / (1 + exp(-k * (x - x0)))) - 1
}

RotatePhase <- Vectorize(function(phase, rotate.hr=0){
  # rotate phase, in hours. If negative add 24. If > 24 minus 24.
  if (is.na(phase)) return(phase)
  phase <- phase + rotate.hr
  if (phase < 0) phase <- phase + 24
  if (phase > 24) phase <- phase - 24
  return(phase)
}, "phase")

PhaseAmpPvalToColor <- function(phase, amp, pval, rotate.hr=0, pval.k = 0.25){
  # rotate.phase: rotate phase by some hours
  phase <- RotatePhase(phase, rotate.hr)
  phase.col <- PhaseToHsv(phase, min.phase = 0, max.phase = 24)
  amp.col <- SaturationCurve(amp, Vmax = 1, k = 2, x0 = 0)
  pval.col <- SaturationCurve(-log10(pval), Vmax = 1, k = pval.k, x0 = 0)
  # convert bad phase amp pvals to 0,0,0
  bad.i <- which(is.na(phase.col) | is.na(amp.col) | is.na(pval.col))
  # make it 0 0 0
  if (length(bad.i) > 0){
    warning(paste0("Converting to 0,0,0 for:", paste0(bad.i, collapse = ",")))
    phase.col[bad.i] <- 0; amp.col[bad.i] <- 0; pval.col[bad.i] <- 0
  }
  hsv.col <- hsv(phase.col, amp.col, pval.col)
  return(hsv.col)
}