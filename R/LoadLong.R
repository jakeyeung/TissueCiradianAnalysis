LoadLong <- function(array.path, rna.seq.path, scale.factor = 100, pseudocount = 1){
  if (missing(array.path) & missing(rna.seq.path)){
    array.path <- "/home/yeung/projects/tissue-specificity/data/exprs_matrices/array_adj_to_kallisto.slope07.txt"
    rna.seq.path <- "/home/yeung/projects/tissue-specificity/data/kallisto/abundance.genecounts.matrix.txt"
  }
  source("/home/yeung/projects/tissue-specificity/scripts/functions/LoadKallistoGene.R")
  kallisto.wide <- LoadKallistoGene(rna.seq.path, form = "wide")  # adjusts colnames to match array
  array.wide <- read.table(array.path)
  # remove rows with negatives (should be 7 of them)
  problem.genes <- rownames(array.wide[which(apply(array.wide, 1, min) < 0), ])
  print("Problem genes:")
  print(problem.genes)
  good.genes <- setdiff(rownames(array.wide), problem.genes)
  array.wide <- array.wide[good.genes, ]
  common.genes <- intersect(rownames(kallisto.wide), rownames(array.wide))
  kallisto.sub <- kallisto.wide[common.genes, ]
  array.sub <- array.wide[common.genes, ]
  
  tissues.rnaseq <- GetTissues(colnames(kallisto.sub), get_unique = FALSE)
  times.rnaseq <- GetTimes(colnames(kallisto.sub), get_unique = FALSE)
  tissues.array <- GetTissues(colnames(array.sub), get_unique = FALSE)
  times.array <- GetTimes(colnames(array.sub), get_unique = FALSE)
  
  ka.long <- data.frame(gene = c(rep(rownames(kallisto.sub), ncol(kallisto.sub)), rep(rownames(array.sub), ncol(array.sub))),
                        tissue = c(rep(tissues.rnaseq, each = nrow(kallisto.sub)), rep(tissues.array, each = nrow(array.sub))),
                        time = as.numeric(c(rep(times.rnaseq, each = nrow(kallisto.sub)), rep(times.array, each = nrow(array.sub)))),
                        exprs = c(unlist(kallisto.sub), unlist(array.sub)),
                        experiment = c(rep("rnaseq", nrow(kallisto.sub) * ncol(kallisto.sub)), rep("array", nrow(array.sub) * ncol(array.sub))))
  ka.long$exprs <- log2(scale.factor * ka.long$exprs + pseudocount)
  return(ka.long)
}