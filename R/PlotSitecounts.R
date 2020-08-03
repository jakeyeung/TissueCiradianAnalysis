PlotSitecounts <- function(vec, filter.top = 40, title = 'MyTitle', jcex = 0.5){
  # Given vector of sitecounts for a gene and labels for its motifs,
  # plot a labelled bar graph
  #
  # Remove zeros
  vec <- vec[which(vec > 0), drop = FALSE]
  i.ordered <- order(vec, decreasing = TRUE)
  # Filter top 20
  i.ordered <- head(i.ordered, n = filter.top)
  vec.ord <- vec[i.ordered]
  heights <- rev(vec.ord)  # first bar at top
  barplot(height = heights, names.arg = names(heights), main = title, horiz = TRUE, las = 1, cex.names=jcex)
}
