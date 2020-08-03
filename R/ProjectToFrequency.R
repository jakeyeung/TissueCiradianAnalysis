ProjectToFrequency <- function(df, my.omega, normalize = TRUE, rhythmic.only = FALSE, method = "ANOVA", pval.cutoff = 5e-3){
  # Perform fourier transform and normalize across all frequencies (squared and square root)
  # 
  # Input:
  # df: long dataframe containing expression and time columns for one condition. Expect 'exprs' and 'time' columns.
  # to be transformed to frequency domain.
  # my.omega: the omega of interest.
  # normalize: converts transform into a sort of z-score
  # rhytmic.only: transforms only genes that are rhythmic (by BIC method)
  # 
  # omega in which we are interested.  
  if (rhythmic.only){
    if (IsRhythmic(df, my.omega, pval.cutoff = pval.cutoff, method = "ANOVA")){
      df.transformed <- Transform(df, my.omega, normalize)
    } else {
      df.transformed <- data.frame(exprs.transformed = 0)
    }
  } else {
    df.transformed <- Transform(df, my.omega, normalize)
  }
  return(df.transformed)
}
