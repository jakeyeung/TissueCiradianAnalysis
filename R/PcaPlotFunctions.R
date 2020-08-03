# Jake Yeung
# PCA plotting functions
# Periodogram plotting functions
# Nov 5 2014
# 
PlotLoadings <- function(Loadings, title="Plot title", plot.colors, cex = 1) {
  # Given vector from PCA, plot vector and color by tissue.
  # Maybe give fancy legend
  if (missing(plot.colors)){
    plot.colors <- rep(1:12, each=24)
  }
  plot(Loadings, main=title, col=plot.colors, type='o',
       cex.axis = cex,
       cex.main = cex,
       cex.lab = cex)
}
