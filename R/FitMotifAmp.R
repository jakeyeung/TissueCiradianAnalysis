FitMotifAmp <- function(dat){
  fit <- lm(motevo.value ~ relamp, dat)
  int <- coef(fit)[["(Intercept)"]]
  relamp <- coef(fit)[["relamp"]]
  if (is.na(relamp)){
    pval <- NA
  } else {
    pval <- summary(fit)$coefficients[["relamp", "Pr(>|t|)"]]
  }
  return(data.frame(int = int, relamp = relamp, pval = pval))
  #   return(data.frame(int = coef(fit)[["(Intercept)"]], 
  #                     relamp = coef(fit)[["motevo.value.norm"]],
  #                     pval = summary(fit)$coefficients[["motevo.value.norm", "Pr(>|t|)"]]))
}
