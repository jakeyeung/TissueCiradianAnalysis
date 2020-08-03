GetTopMotifs <- function(type="rhythmic"){
  # type: rhythmic | tissue
  if (type == "rhythmic"){
    inf <- "/home/yeung/projects/tissue-specificity/data/gene_lists/motif_lists/top.rhyth.txt"
  } else if (type == "tissue"){
    inf <- "/home/yeung/projects/tissue-specificity/data/gene_lists/motif_lists/top.tiss.txt"
  }
  motifs <- unlist(read.table(inf, col.names = "motif", stringsAsFactors = FALSE), use.names = FALSE)
  return(motifs)
}